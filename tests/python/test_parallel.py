import numpy as np
import pytest
from definitions import GEOSET_CLASSES, get_dim_range

import geosets_py as gsp
from geosets_py import NotImplemented

BATCH = 32


def _typestr(SetType):
    return SetType.__name__.lower()


def _make_sets(SetType, dim, n=BATCH, seed=0):
    rng = np.random.default_rng(seed)
    sets = []
    for _ in range(n):
        s = SetType.from_unit_box(dim)
        s.translate_(rng.uniform(-0.5, 0.5, dim))
        sets.append(s)
    return sets


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_batched_volume(SetType):
    for dim in get_dim_range(SetType):
        sets = _make_sets(SetType, dim)
        fn = getattr(gsp, f"batched_volume_{_typestr(SetType)}")
        try:
            expected = [s.volume() for s in sets]
        except NotImplemented:
            pytest.skip()
        got = fn(sets)
        assert np.allclose(np.asarray(got), np.asarray(expected))


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_batched_is_empty(SetType):
    for dim in get_dim_range(SetType):
        sets = _make_sets(SetType, dim)
        fn = getattr(gsp, f"batched_is_empty_{_typestr(SetType)}")
        expected = [s.is_empty() for s in sets]
        got = list(fn(sets))
        assert got == expected


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_batched_bounded(SetType):
    for dim in get_dim_range(SetType):
        sets = _make_sets(SetType, dim)
        fn = getattr(gsp, f"batched_bounded_{_typestr(SetType)}")
        expected = [s.bounded() for s in sets]
        got = list(fn(sets))
        assert got == expected


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_batched_degenerate(SetType):
    for dim in get_dim_range(SetType):
        sets = _make_sets(SetType, dim)
        fn = getattr(gsp, f"batched_degenerate_{_typestr(SetType)}")
        try:
            expected = [s.degenerate() for s in sets]
        except NotImplemented:
            pytest.skip()
        got = list(fn(sets))
        assert got == expected


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_batched_to_vertices(SetType):
    for dim in get_dim_range(SetType):
        sets = _make_sets(SetType, dim)
        fn = getattr(gsp, f"batched_to_vertices_{_typestr(SetType)}")
        try:
            expected = [s.to_vertices() for s in sets]
        except NotImplemented:
            pytest.skip()
        got = fn(sets)
        for g, e in zip(got, expected):
            assert np.allclose(np.asarray(g), np.asarray(e))


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_batched_center(SetType):
    for dim in get_dim_range(SetType):
        sets = _make_sets(SetType, dim)
        fn = getattr(gsp, f"batched_center_{_typestr(SetType)}")
        expected = np.stack([s.center(validate=False) for s in sets], axis=0)
        got = fn(sets, validate=False)
        assert np.allclose(np.asarray(got), expected)


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_batched_support_function(SetType):
    for dim in get_dim_range(SetType):
        sets = _make_sets(SetType, dim)
        rng = np.random.default_rng(1)
        directions = rng.normal(size=(BATCH, dim))
        fn = getattr(gsp, f"batched_support_function_{_typestr(SetType)}")
        expected = [s.support_function(directions[i], validate=False) for i, s in enumerate(sets)]
        got = fn(sets, directions, validate=False)
        assert len(got) == len(expected)
        for g, e in zip(got, expected):
            assert np.isclose(g.value, e.value)
            assert np.allclose(np.asarray(g.vector), np.asarray(e.vector))


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_batched_tangent_hyperplane(SetType):
    for dim in get_dim_range(SetType):
        sets = _make_sets(SetType, dim)
        rng = np.random.default_rng(2)
        directions = rng.normal(size=(BATCH, dim))
        fn = getattr(gsp, f"batched_tangent_hyperplane_{_typestr(SetType)}")
        try:
            expected = [
                s.tangent_hyperplane(directions[i], validate=False) for i, s in enumerate(sets)
            ]
        except NotImplemented:
            pytest.skip()
        got = fn(sets, directions, validate=False)
        assert len(got) == len(expected)
        for g, e in zip(got, expected):
            assert np.isclose(g.b, e.b)
            assert np.allclose(np.asarray(g.a), np.asarray(e.a))


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_batched_boundary_point(SetType):
    for dim in get_dim_range(SetType):
        sets = _make_sets(SetType, dim)
        rng = np.random.default_rng(3)
        directions = rng.normal(size=(BATCH, dim))
        fn = getattr(gsp, f"batched_boundary_point_{_typestr(SetType)}")
        try:
            expected = np.stack(
                [s.boundary_point(directions[i], validate=False) for i, s in enumerate(sets)],
                axis=0,
            )
        except NotImplemented:
            pytest.skip()
        got = fn(sets, directions, validate=False)
        assert np.allclose(np.asarray(got), expected)


@pytest.mark.parametrize("SetType", GEOSET_CLASSES)
def test_batched_contains_point(SetType):
    for dim in get_dim_range(SetType):
        sets = _make_sets(SetType, dim)
        rng = np.random.default_rng(4)
        points = rng.uniform(-1.5, 1.5, size=(BATCH, dim))
        fn = getattr(gsp, f"batched_contains_point_{_typestr(SetType)}")
        expected = [s.contains_point(points[i], validate=False) for i, s in enumerate(sets)]
        got = list(fn(sets, points, validate=False))
        assert got == expected
