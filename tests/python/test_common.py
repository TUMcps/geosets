import numpy as np
import pytest
from definitions import GEOSET_CLASSES, get_dim_range

from geosets_py import NotImplemented


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_dim_common(SetType):
    for dim in get_dim_range(SetType):
        s = SetType.from_unit_box(dim)
        assert s.dim() == dim


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_center_common(SetType):
    for dim in get_dim_range(SetType):
        s = SetType.from_unit_box(dim)
        assert np.allclose(s.center(), np.zeros(dim))


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_volume_common(SetType):
    for dim in get_dim_range(SetType):
        s = SetType.from_unit_box(dim)
        try:
            volume = s.volume()
            assert np.isclose(volume, 2**dim)
        except NotImplemented:
            pytest.skip()


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_empty_common(SetType):
    for dim in get_dim_range(SetType):
        s = SetType.from_unit_box(dim)
        assert not s.is_empty()


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_bounded_common(SetType):
    for dim in get_dim_range(SetType):
        s = SetType.from_unit_box(dim)
        assert s.bounded()


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_translate_center_common(SetType):
    for dim in get_dim_range(SetType):
        s = SetType.from_unit_box(dim)
        assert np.allclose(s.center(), np.zeros(dim))

        for _ in range(10):
            vec = np.random.rand(dim)
            s_translated = s.translate(vec)

            assert np.allclose(s_translated.center(), vec)


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_degenerate_common(SetType):
    for dim in get_dim_range(SetType):
        s = SetType.from_unit_box(dim)
        try:
            is_degenerate = s.degenerate()
            assert not is_degenerate
        except NotImplemented:
            pytest.skip()


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_linear_transform_common(SetType):
    for dim in get_dim_range(SetType):
        s = SetType.from_unit_box(dim)
        try:
            matrix = np.eye(dim) * 2.0
            s.linear_transform_(matrix)
            center = s.center()
            assert np.allclose(center, np.zeros(dim))
        except NotImplemented:
            pytest.skip()


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_support_function_common(SetType):
    s = SetType.from_unit_box(2)
    try:
        support = s.support_function(np.array([1.0, 0.0]))
        assert support.vector[0] == 1.0
        assert -1.0 <= support.vector[1] <= 1.0
        assert support.value == 1.0

        support = s.support_function(np.array([0.0, 1.0]))
        assert support.vector[1] == 1.0
        assert -1.0 <= support.vector[0] <= 1.0
        assert support.value == 1.0

        support = s.support_function(np.array([1.0, 1.0]))
        assert np.allclose(support.vector, np.array([1.0, 1.0]))
        assert support.value == 2.0

    except NotImplemented:
        pytest.skip()


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_tangent_hyperplane_axis_aligned(SetType):
    s = SetType.from_unit_box(2)
    try:
        # Direction (1, 0): plane x = 1
        hp = s.tangent_hyperplane(np.array([1.0, 0.0]))
        assert np.allclose(hp.a, [1.0, 0.0])
        assert np.isclose(hp.b, 1.0)

        # Direction (-1, 0): plane -x = 1
        hp = s.tangent_hyperplane(np.array([-1.0, 0.0]))
        assert np.allclose(hp.a, [-1.0, 0.0])
        assert np.isclose(hp.b, 1.0)

        # Direction (0, 1): plane y = 1
        hp = s.tangent_hyperplane(np.array([0.0, 1.0]))
        assert np.allclose(hp.a, [0.0, 1.0])
        assert np.isclose(hp.b, 1.0)
    except NotImplemented:
        pytest.skip()


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_tangent_hyperplane_support_on_plane(SetType):
    for dim in get_dim_range(SetType):
        s = SetType.from_unit_box(dim)
        try:
            for _ in range(20):
                d = np.random.randn(dim)
                if np.linalg.norm(d) < 1e-9:
                    continue
                hp = s.tangent_hyperplane(d)
                sp = s.support_function(d)
                # Support point lies on tangent hyperplane
                assert np.isclose(hp.a @ sp.vector, hp.b, atol=1e-6)
                # Normal points outward
                assert hp.a @ d > -1e-9
        except NotImplemented:
            pytest.skip()


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_contains_point_common(SetType):
    for dim in get_dim_range(SetType):
        s = SetType.from_unit_box(dim)
        try:
            samples = np.random.rand(10, dim) * 2 - 1
            for i in range(samples.shape[0]):
                point = samples[i, :]
                assert s.contains_point(point)
        except NotImplemented:
            pytest.skip()


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_boundary_point_common(SetType):
    for dim in get_dim_range(SetType):
        s = SetType.from_unit_box(dim)
        try:
            # Axis-aligned direction
            d = np.zeros(dim)
            d[0] = 1.0
            bp = s.boundary_point(d)
            expected = np.zeros(dim)
            expected[0] = 1.0
            assert np.allclose(bp, expected, atol=1e-6)

            # Random directions: boundary of unit box has inf-norm 1
            for _ in range(20):
                d = np.random.randn(dim)
                if np.linalg.norm(d) < 1e-9:
                    continue
                bp = s.boundary_point(d)
                assert np.isclose(np.max(np.abs(bp)), 1.0, atol=1e-6)
                assert s.contains_point(bp)
        except NotImplemented:
            pytest.skip()


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_line_intersection_common(SetType):
    box = SetType.from_unit_box(2)

    try:
        lower, upper = box.line_intersection(np.zeros(2), np.array([1.0, 0.0]))
    except NotImplemented:
        pytest.skip()

    assert np.allclose(lower, [-1.0, 0.0])
    assert np.allclose(upper, [1.0, 0.0])
