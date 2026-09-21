import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, truncnorm

import geosets_py as gsp


def main():
    offset = 5
    print(offset)
    poly = gsp.HPolytope.from_unit_box(2).translate(np.full(2, offset))
    loc = np.zeros((1, 2))
    scale = np.ones((1, 2))

    trunc_samples = gsp.sampling.batched_gaussian_rdhr_sampling_hpolytope([poly], loc, scale, 1000)
    print(trunc_samples.shape)

    samples = np.random.normal(loc=loc, scale=scale, size=(1000, 1, 2))

    poly.plot()

    # plot samples (10, 1, 1)
    plt.scatter(
        trunc_samples[:, 0, 0], trunc_samples[:, 0, 1], color="red", s=5, label="truncated samples"
    )
    plt.scatter(samples[:, 0, 0], samples[:, 0, 1], color="green", s=5, label="standard samples")

    # mean of truncated (1 dimensions)
    trunc_mean_mc = np.mean(trunc_samples, axis=0)
    trunc_variance_mc = np.var(trunc_samples, axis=0)

    # analytic mean of truncated
    mu, sigma = 0.0, 1.0
    a, b = offset - 1, offset + 1
    alpha = (a - mu) / sigma
    beta = (b - mu) / sigma
    dist = truncnorm(alpha, beta, loc=mu, scale=sigma)
    trunc_mean = dist.mean()
    trunc_var = dist.var()
    print(norm.pdf(alpha) - norm.pdf(beta))

    print("truncated mean (MC):", trunc_mean_mc)
    print("truncated mean (analytic):", trunc_mean)

    print("truncated variance (MC):", trunc_variance_mc)
    print("truncated variance (analytic):", trunc_var)

    # plot means
    plt.scatter(
        trunc_mean_mc[0, 0], trunc_mean_mc[0, 1], color="blue", s=50, label="truncated mean (MC)"
    )
    plt.scatter(trunc_mean, trunc_mean, color="orange", s=50, label="truncated mean (analytic)")

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
