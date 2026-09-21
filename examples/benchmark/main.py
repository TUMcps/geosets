import time
from collections.abc import Callable

import numpy as np
import shapely

import geosets_py

from .python_hpolytope import PythonHPolytope


def creation(n: int, set_creation_fn: Callable) -> float:
    start = time.time()
    for _ in range(n):
        _ = set_creation_fn()
    end = time.time()
    return end - start


def center(n: int, set_creation_fn: Callable) -> float:
    sets = [set_creation_fn() for _ in range(n)]
    start = time.time()
    for i in range(n):
        _ = sets[i].center()
    end = time.time()
    return end - start


def center_batched(n: int, set_creation_fn: Callable) -> float:
    geosets_py.set_omp_num_threads(12)
    sets = [set_creation_fn() for _ in range(n)]
    start = time.time()
    _ = geosets_py.batched_center_hpolytope(sets)
    end = time.time()
    return end - start


def translation(n: int, set_creation_fn: Callable) -> float:
    start = time.time()
    set_o = set_creation_fn()

    for _ in range(n):
        vector = np.random.uniform(-10, 10, size=(set_o.dim(),))
        set_o.translate_(vector)
    end = time.time()

    return end - start


def data_copy(n: int, set_creation_fn: Callable) -> float:
    sets = [set_creation_fn() for _ in range(n)]
    start = time.time()
    for i in range(n):
        vertices = sets[i].to_vertices()
        _ = shapely.geometry.Polygon(vertices)
    end = time.time()
    return end - start


if __name__ == "__main__":

    def hpoly_creation():
        return geosets_py.HPolytope.from_random(2, 10)

    def polygon_creation():
        return geosets_py.Polygon.from_random(2, 10)

    def python_hpoly_creation():
        return PythonHPolytope.from_random(2, 10)

    n_sets = 1_000

    print("-- Ours --")
    print(f"Creation {creation(n_sets, set_creation_fn=hpoly_creation)}")
    print(f"Translation {translation(n_sets, set_creation_fn=hpoly_creation)}")
    print(f"Center {center(n_sets, set_creation_fn=hpoly_creation)}")
    print(f"Center batched {center_batched(n_sets, set_creation_fn=hpoly_creation)}")
    print()

    print("-- Python --")
    print(f"Creation {creation(n_sets, set_creation_fn=python_hpoly_creation)}")
    print(f"Translation {translation(n_sets, set_creation_fn=python_hpoly_creation)}")
    print(f"Center {center(n_sets, set_creation_fn=python_hpoly_creation)}")
    print(f"Copying data to shapely {data_copy(n_sets, set_creation_fn=polygon_creation)}")
