import time

import numpy as np

import geosets_py as gsp


def main():
    batch_size = 128
    dim = 3
    # gsp.set_omp_num_threads(4)

    sets = [gsp.HPolytope.from_random(dim, dim * 3) for _ in range(batch_size)]
    # sets = [gsp.Zonotope.from_random(dim, dim * 3) for _ in range(batch_size)]
    points = np.random.uniform(-1, 1, (batch_size, dim))
    print(points.shape)

    # Python loop
    start = time.time()
    distances = []
    for i in range(batch_size):
        _ = sets[i].center(validate=False)
        boundary_point = sets[i].boundary_point(points[i], validate=False)
        distances.append(np.linalg.norm(points[i] - boundary_point))

    end = time.time()
    print(f"Mean distance to boundary: {np.mean(distances):.4f}")
    print(f"Python computation time: {end - start:.4f} seconds")

    # C++ loop
    # sets = [s.copy() for s in sets]
    sets = [gsp.HPolytope(s.A, s.b) for s in sets]
    start = time.time()
    _ = gsp.batched_center_hpolytope(sets)
    # centers = gsp.batched_center_zonotope(sets)
    boundary_points = gsp.batched_boundary_point_hpolytope(sets, points)
    # boundary_points = gsp.batched_boundary_point_zonotope(sets, points)
    print(type(boundary_points))
    distances = np.linalg.norm(points - boundary_points, axis=1)
    end = time.time()
    print(f"Mean distance to boundary: {np.mean(distances):.4f}")
    print(f"C++ computation time: {end - start:.4f} seconds")


if __name__ == "__main__":
    main()
