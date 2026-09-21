import numpy as np
import pytest
from definitions import TOL, get_dim_range, vertices_to_set

from geosets_py import (
    DimensionMismatchException,
    GeometricSet,
    HPolytope,
    Interval,
    NotImplemented,
    Polygon,
    VPolytope,
    Zonotope,
    do_intersect,
    intersection,
    minkowski_sum,
)

# Test types that support binary operations
BINARY_OP_GEOSET_CLASSES: list[type[GeometricSet]] = [
    Interval,
    VPolytope,
    HPolytope,
    Polygon,
    Zonotope,
]


@pytest.mark.parametrize("GeoSet", BINARY_OP_GEOSET_CLASSES)
def test_dimension_compatibility(GeoSet):
    """Test that binary operations check dimension compatibility."""
    if GeoSet is Polygon:
        return
    try:
        a2d = GeoSet.from_unit_box(2)
        b3d = GeoSet.from_unit_box(3)

        with pytest.raises(DimensionMismatchException):
            intersection(a2d, b3d)

        with pytest.raises(DimensionMismatchException):
            do_intersect(a2d, b3d)

        with pytest.raises(DimensionMismatchException):
            minkowski_sum(a2d, b3d)
    except NotImplemented:
        pytest.skip()


@pytest.mark.parametrize("GeoSet", BINARY_OP_GEOSET_CLASSES)
def test_intersection_identical_sets(GeoSet):
    """Intersection of identical sets returns the same set."""
    for dim in get_dim_range(GeoSet):
        a = GeoSet.from_unit_box(dim)
        try:
            result = intersection(a, a)

            assert result.dim() == a.dim()
            assert np.allclose(result.center(), a.center())
            assert np.isclose(result.volume(), a.volume())
            assert not result.is_empty()
        except NotImplemented:
            pytest.skip()


@pytest.mark.parametrize("GeoSet", BINARY_OP_GEOSET_CLASSES)
def test_intersect_predicate(GeoSet):
    for dim in get_dim_range(GeoSet):
        a = GeoSet.from_unit_box(dim)
        b = GeoSet.from_unit_box(dim).translate(np.zeros(dim) + 0.5)
        disjoint = GeoSet.from_unit_box(dim).translate(np.zeros(dim) + 3.0)
        touching = GeoSet.from_unit_box(dim).translate(np.eye(1, dim, 0).reshape(-1) * 2.0)

        try:
            assert do_intersect(a, b)
            assert do_intersect(b, a)
            assert not do_intersect(a, disjoint)
            assert not do_intersect(disjoint, a)
            assert do_intersect(a, touching)
        except NotImplemented:
            pytest.skip()


@pytest.mark.parametrize("GeoSet", BINARY_OP_GEOSET_CLASSES)
def test_intersection_contained_in_operands(GeoSet):
    """Intersection result is contained in both operands."""
    for dim in get_dim_range(GeoSet):
        a = GeoSet.from_unit_box(dim)
        b = GeoSet.from_unit_box(dim).translate(np.zeros(dim) + 0.5)

        try:
            result = intersection(a, b)
            assert a.contains_point(result.center())
            assert b.contains_point(result.center())

            # Should be commutative
            ab = intersection(a, b)
            ba = intersection(b, a)
            assert np.allclose(ab.center(), ba.center())
            assert np.isclose(ab.volume(), ba.volume())
        except NotImplemented:
            pytest.skip()


@pytest.mark.parametrize("GeoSet", BINARY_OP_GEOSET_CLASSES)
def test_minkowski_sum_2d_unit_boxes(GeoSet):
    """Test Minkowski sum of 2D unit boxes."""
    try:
        a = GeoSet.from_unit_box(2)
        b = GeoSet.from_unit_box(2)

        ab = minkowski_sum(a, b)
        ba = minkowski_sum(b, a)

        expected = np.array([[2.0, -2.0], [2.0, 2.0], [-2.0, 2.0], [-2.0, -2.0]])

        assert vertices_to_set(ab.to_vertices()) == vertices_to_set(expected)
        assert vertices_to_set(ba.to_vertices()) == vertices_to_set(expected)

        assert np.allclose(ab.center(), np.zeros(2))
        assert np.allclose(ba.center(), np.zeros(2))
        assert abs(ab.volume() - 16.0) < TOL
        assert abs(ba.volume() - 16.0) < TOL
    except NotImplemented:
        pytest.skip()


@pytest.mark.parametrize("GeoSet", BINARY_OP_GEOSET_CLASSES)
def test_minkowski_sum_3d_unit_boxes(GeoSet):
    """Test Minkowski sum of 3D unit boxes."""
    if GeoSet is Polygon:
        return
    try:
        a = GeoSet.from_unit_box(3)
        b = GeoSet.from_unit_box(3)

        ab = minkowski_sum(a, b)
        ba = minkowski_sum(b, a)

        expected = np.array(
            [
                [2.0, -2.0, -2.0],
                [2.0, -2.0, 2.0],
                [2.0, 2.0, 2.0],
                [2.0, 2.0, -2.0],
                [-2.0, 2.0, 2.0],
                [-2.0, 2.0, -2.0],
                [-2.0, -2.0, 2.0],
                [-2.0, -2.0, -2.0],
            ]
        )

        assert vertices_to_set(ab.to_vertices()) == vertices_to_set(expected)
        assert vertices_to_set(ba.to_vertices()) == vertices_to_set(expected)

        assert np.allclose(ab.center(), np.zeros(3))
        assert np.allclose(ba.center(), np.zeros(3))
        assert abs(ab.volume() - 64.0) < TOL
        assert abs(ba.volume() - 64.0) < TOL
    except NotImplemented:
        pytest.skip()
