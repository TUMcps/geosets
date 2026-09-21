import numpy as np
import pytest
from definitions import GEOSET_CLASSES

import geosets_py as gsp

# Sampler name prefixes; the concrete function is per set type:
#   getattr(gsp.sampling, f"{prefix}_sampling_{typestr}")
SAMPLERS = ["gaussian_rejection", "gaussian_rdhr", "gaussian_rejection_rdhr_fallback"]
BATCHED_SAMPLERS = ["batched_" + s for s in SAMPLERS]


def get_sampler(prefix, SetType):
    return getattr(gsp.sampling, f"{prefix}_sampling_{SetType.__name__.lower()}")


def run_or_skip(fn):
    """Mirror of the C++ skip_if_unimplemented: RDHR/combined need line_intersection,
    which some set types (VPolytope, Zonotope) do not implement -> skip those."""
    try:
        return fn()
    except gsp.NotImplemented:
        pytest.skip("line intersection not implemented for this set type")


def assert_single_shape_and_containment(samples, set_object, samples_per_set):
    assert isinstance(samples, np.ndarray)
    assert samples.shape == (samples_per_set, 2)
    assert samples.flags.c_contiguous
    for sample in samples:
        # loosened tol: MCMC samples can land essentially on the boundary
        assert set_object.contains_point(sample, 1e-6)


def assert_batched_shape_and_containment(samples, sets, samples_per_set):
    # Batched layout (pack_samples_to_numpy): (samples_per_set, num_sets, dim), C-contiguous.
    assert isinstance(samples, np.ndarray)
    assert samples.shape == (samples_per_set, len(sets), 2)
    assert samples.flags.c_contiguous
    for sample_idx in range(samples_per_set):
        for set_idx, set_object in enumerate(sets):
            assert set_object.contains_point(samples[sample_idx, set_idx], 1e-6)


# Every returned sample must lie inside the set (a real check for the MCMC paths, not just
# shape), and n == 0 must yield an empty (0, dim) array.
@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
@pytest.mark.parametrize("prefix", SAMPLERS)
def test_gaussian_sampling_shape_and_containment(SetType, prefix):
    set_object = SetType.from_unit_box(2)
    sample = get_sampler(prefix, SetType)
    loc = np.zeros(2)
    scale = np.full(2, 0.5)

    gsp.random_seed(17)
    samples = run_or_skip(lambda: sample(set_object, loc, scale, 16))
    assert_single_shape_and_containment(samples, set_object, 16)

    empty = sample(set_object, loc, np.ones(2), 0)
    assert empty.shape == (0, 2)


# The one distribution-level test. The unit box and a loc=0 Gaussian are both symmetric under
# x -> -x, and RDHR starts from the (symmetric) center, so every sample's marginal is mean-zero
# for ANY walk length -- only sampling noise remains. This sidesteps RDHR's weak mixing
# (walk_length = dim**3 = 8 at dim 2) that would make an off-center mean check flaky.
@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
@pytest.mark.parametrize("prefix", SAMPLERS)
def test_gaussian_sampling_symmetric_set_is_centered(SetType, prefix):
    set_object = SetType.from_unit_box(2)
    sample = get_sampler(prefix, SetType)
    loc = np.zeros(2)
    scale = np.full(2, 0.6)

    gsp.random_seed(101)
    samples = run_or_skip(lambda: sample(set_object, loc, scale, 256))
    mean = samples.mean(axis=0)
    assert np.all(np.abs(mean) < 0.12)  # ~4 * standard error


# The only source of nondeterminism is the global RNG: same seed => identical output.
@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
@pytest.mark.parametrize("prefix", SAMPLERS)
def test_gaussian_sampling_is_reproducible(SetType, prefix):
    set_object = SetType.from_unit_box(2)
    sample = get_sampler(prefix, SetType)
    loc = np.zeros(2)
    scale = np.full(2, 0.5)

    gsp.random_seed(7)
    a = run_or_skip(lambda: sample(set_object, loc, scale, 8))
    gsp.random_seed(7)
    b = sample(set_object, loc, scale, 8)
    np.testing.assert_array_equal(a, b)


# When acceptance is easy no fallback fires, so the combined sampler consumes the RNG stream
# exactly like pure rejection and produces an identical result.
@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_combined_equals_rejection_without_fallback(SetType):
    set_object = SetType.from_unit_box(2)
    rejection = get_sampler("gaussian_rejection", SetType)
    combined = get_sampler("gaussian_rejection_rdhr_fallback", SetType)
    loc = np.zeros(2)
    scale = np.full(2, 0.5)

    gsp.random_seed(31)
    r = rejection(set_object, loc, scale, 8, 10000)
    gsp.random_seed(31)
    c = run_or_skip(lambda: combined(set_object, loc, scale, 8, 10000))
    np.testing.assert_allclose(c, r, atol=1e-12)


# Differential behaviour justifying three functions: with the Gaussian mass parked far outside
# the set and a tiny budget, exact rejection gives up (RuntimeError) while RDHR and the combined
# fallback still deliver contained samples via the interior-point (center) fallback -- except on
# set types whose line intersection is unimplemented, where RDHR/combined skip.
@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_rejection_exhausts_budget_where_rdhr_and_combined_survive(SetType):
    set_object = SetType.from_unit_box(2)
    rejection = get_sampler("gaussian_rejection", SetType)
    rdhr = get_sampler("gaussian_rdhr", SetType)
    combined = get_sampler("gaussian_rejection_rdhr_fallback", SetType)
    loc = np.full(2, 5.0)  # outside [-1, 1]^2
    scale = np.full(2, 0.1)

    gsp.random_seed(41)
    with pytest.raises(RuntimeError, match="sample index"):
        rejection(set_object, loc, scale, 8, 1)

    rdhr_samples = run_or_skip(lambda: rdhr(set_object, loc, scale, 8))
    assert_single_shape_and_containment(rdhr_samples, set_object, 8)

    combined_samples = run_or_skip(lambda: combined(set_object, loc, scale, 8, 1))
    assert_single_shape_and_containment(combined_samples, set_object, 8)


# Shared input validation (all three route through validate_sample_inputs). Non-finite loc,
# non-positive scale, and loc/scale dimension mismatch all raise ValueError before any sampling.
# (Empty/infinite sets are covered in C++; they are not constructible from Python.)
@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
@pytest.mark.parametrize("prefix", SAMPLERS)
def test_gaussian_sampling_rejects_invalid_inputs(SetType, prefix):
    set_object = SetType.from_unit_box(2)
    sample = get_sampler(prefix, SetType)
    loc = np.zeros(2)
    scale = np.ones(2)

    bad_loc = loc.copy()
    bad_loc[0] = np.inf
    with pytest.raises(ValueError):
        sample(set_object, bad_loc, scale, 4)

    bad_scale = scale.copy()
    bad_scale[0] = 0.0
    with pytest.raises(ValueError):
        sample(set_object, loc, bad_scale, 4)

    with pytest.raises(ValueError):
        sample(set_object, np.zeros(3), scale, 4)
    with pytest.raises(ValueError):
        sample(set_object, loc, np.ones(3), 4)


# Batched output has shape (n_samples_per_set, num_sets, dim), is C-contiguous, and every
# sample lies in its set. locs/scales carry one row per set.
@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
@pytest.mark.parametrize("prefix", BATCHED_SAMPLERS)
def test_batched_gaussian_sampling_shape_and_containment(SetType, prefix):
    sets = [SetType.from_unit_box(2), SetType.from_unit_box(2)]
    sampler = get_sampler(prefix, SetType)
    locs = np.zeros((2, 2))
    scales = np.full((2, 2), 0.25)
    n = 4

    gsp.random_seed(19)
    samples = run_or_skip(lambda: sampler(sets, locs, scales, n))
    assert_batched_shape_and_containment(samples, sets, n)


# Attribution: three disjoint boxes (10 units apart) each with its Gaussian centered inside.
# samples[:, i] must land in sets[i] and in no other set -- this fails if the packing layout
# (output_row = sample_idx * num_sets + set_idx) attributes a column to the wrong set.
@pytest.mark.parametrize("prefix", BATCHED_SAMPLERS)
def test_batched_gaussian_sampling_attributes_samples_to_correct_set(prefix):
    centers = [np.array([0.0, 0.0]), np.array([10.0, 0.0]), np.array([0.0, 10.0])]
    sets = [gsp.HPolytope.from_unit_box(2).translate(c) for c in centers]
    sampler = getattr(gsp.sampling, f"{prefix}_sampling_hpolytope")
    locs = np.array(centers)  # Gaussian centered in each box
    scales = np.full((3, 2), 0.25)
    n = 8

    gsp.random_seed(23)
    samples = run_or_skip(lambda: sampler(sets, locs, scales, n))
    assert samples.shape == (n, 3, 2)

    for set_idx, set_object in enumerate(sets):
        for sample_idx in range(n):
            point = samples[sample_idx, set_idx]
            assert set_object.contains_point(point, 1e-6)
            for other_idx, other in enumerate(sets):
                if other_idx != set_idx:
                    assert not other.contains_point(point, 1e-6)
