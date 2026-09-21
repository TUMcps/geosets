import numpy as np

from geosets_py import GeometricSet, HPolytope, Interval, Polygon, VPolytope, Zonotope

GEOSET_CLASSES: list[type[GeometricSet]] = [HPolytope, VPolytope, Polygon, Interval, Zonotope]

TOL = 1e-9


def get_dim_range(T: type[GeometricSet]) -> range:
    if T is Polygon:
        return range(2, 3)
    return range(2, 5)


def vertices_to_set(vertices: np.ndarray) -> set:
    """Convert vertices matrix to a set for order-independent comparison."""
    result = set()
    for i in range(vertices.shape[0]):
        # Convert row to tuple of floats for hashing
        row_tuple = tuple(vertices[i, :].tolist())
        result.add(row_tuple)
    return result
