"""CORA-COMP instance runner: python run.py <params-json> <results-file>.

Implements the operations of https://github.com/CORA-COMP/benchmarks (cpu only).
"""

import json
import sys

import numpy as np

import geosets_py as gs


def rand_set(rng, p):
    """Random set following CORA's generateRandom."""
    n = p["dim"]
    if p["set"] == "interval":
        c = rng.uniform(-2, 2, n)
        r = rng.uniform(0, 10, n) * rng.uniform(0, 1, n) / 2
        return gs.Interval(c - r, c + r)
    G = rng.standard_normal((n, p["generators"]))
    G *= rng.uniform(0, 1, G.shape[1]) / np.linalg.norm(G, axis=0)
    return gs.Zonotope(10 * rng.standard_normal(n), G)


def rand_points(rng, S, k):
    """k points of CORA's randPoint 'standard', shape (k, dim)."""
    if isinstance(S, gs.Interval):
        return rng.uniform(S.lb, S.ub, (k, S.dim()))
    return S.c + rng.uniform(-1, 1, (k, S.n_generators())) @ S.G.T


def run(p):
    rng = np.random.default_rng(0)
    op, n, reps = p["operation"], p["dim"], p["repetition"]
    bs = p.get("batch_size", 1)
    support = getattr(gs, f"batched_support_function_{p['set']}")
    contains = getattr(gs, f"batched_contains_point_{p['set']}")

    if op == "startup":
        rand_set(rng, p)
        return
    if op == "generateRandom":
        for _ in range(reps):
            [rand_set(rng, p) for _ in range(bs)]
        return

    sets = [rand_set(rng, p) for _ in range(bs)]
    if op == "randPoint":
        for _ in range(reps):
            [rand_points(rng, S, p["points"]) for S in sets]
    elif op == "supportFunc":
        D = rng.standard_normal((bs, n))
        D /= np.linalg.norm(D, axis=1, keepdims=True)
        for _ in range(reps):
            [s.value for s in support(sets, D)]
    elif op == "matMul":
        M = rng.standard_normal((n, n))
        for _ in range(reps):
            [S.linear_transform(M) for S in sets]
    elif op == "minkSum":
        others = [rand_set(rng, p) for _ in range(bs)]
        for _ in range(reps):
            [gs.minkowski_sum(a, b) for a, b in zip(sets, others)]
    elif op == "contains":
        # (points, bs, dim): row j holds the j-th point of every set
        P = np.stack([rand_points(rng, S, p["points"]) for S in sets], axis=1)
        for _ in range(reps):
            if not all(all(contains(sets, pts)) for pts in P):
                raise RuntimeError("sampled point reported outside its set")
    else:
        return "unsupported"


def main():
    params, results_file = json.loads(sys.argv[1]), sys.argv[2]
    try:
        verdict = run(params) or "finished"
    except Exception as e:  # noqa: BLE001
        print(f"error: {e!r}", file=sys.stderr)
        verdict = "error"
    with open(results_file, "w") as f:
        f.write(f"result\n{verdict}\n")


if __name__ == "__main__":
    main()
