import numpy as np

import geosets_py as gsp


def main():
    a = gsp.Zonotope.from_random(2, 4)
    locs = np.zeros((1, 2))
    scales = np.ones((1, 2))

    gsp.batched_is_empty_zonotope([a])

    samples = gsp.sampling.batched_gaussian_rejection_sampling_zonotope([a], locs, scales)
    print(samples.shape)


if __name__ == "__main__":
    main()
