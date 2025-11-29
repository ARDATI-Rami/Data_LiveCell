import numpy as np

"""Test-side implementations and unit tests for topology metrics.

These helpers are small, self-contained stand-ins for the topology
logic (Dirichlet posterior mean, topological charge, divergences,
Wasserstein distance). They live only in the test suite so that we can
exercise and lock in the mathematical behavior without touching any
analysis scripts.
"""


def dirichlet_posterior_mean(counts, alpha=0.5):
    counts = np.asarray(counts, dtype=float)
    alpha_vec = np.full_like(counts, fill_value=alpha, dtype=float)
    posterior = counts + alpha_vec
    total = posterior.sum()
    if total == 0:
        # if all counts are zero and alpha is zero, fall back to uniform
        return np.ones_like(posterior) / posterior.size
    return posterior / total


def topological_charge(counts, support):
    counts = np.asarray(counts, dtype=float)
    support = np.asarray(support, dtype=float)
    if counts.sum() == 0:
        return 0.0, 0.0
    p = counts / counts.sum()
    # standard definition: q = 6 - n, weighted by polygon distribution
    q_values = (6.0 - support) * p
    mean_q = q_values.sum()
    mean_abs_q = np.abs(q_values).sum()
    return mean_q, mean_abs_q


def js_divergence(p, q, base=2.0):
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    # normalise internally to be robust to simple scaling
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)

    def _kl(a, b):
        mask = (a > 0) & (b > 0)
        if not np.any(mask):
            return 0.0
        return np.sum(a[mask] * (np.log(a[mask] / b[mask]) / np.log(base)))

    return 0.5 * _kl(p, m) + 0.5 * _kl(q, m)


def wasserstein1_discrete(p, q, support):
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    support = np.asarray(support, dtype=float)
    p = p / p.sum()
    q = q / q.sum()
    # cumulative distributions over sorted support
    cdf_p = np.cumsum(p)
    cdf_q = np.cumsum(q)
    # W1 distance is the integral of |CDF_p - CDF_q| over the support
    # For discrete distributions, this is sum over intervals of |CDF_p - CDF_q| * interval_width
    # We need n-1 intervals for n support points
    if len(support) == 1:
        return 0.0
    intervals = np.diff(support)
    # Use the CDF values at the right endpoint of each interval (index 0 to n-2)
    # Actually, for proper integration we should use the CDF difference in each interval
    # The correct formula: sum over i of |cdf_p[i] - cdf_q[i]| * (support[i+1] - support[i])
    cdf_diffs = np.abs(cdf_p[:-1] - cdf_q[:-1])
    return np.sum(cdf_diffs * intervals)


# ---------------------------------------------------------------------------
# Unit tests for the helpers above
# ---------------------------------------------------------------------------


def test_dirichlet_posterior_mean_sums_to_one():
    counts = np.array([10, 20, 30])
    probs = dirichlet_posterior_mean(counts, alpha=0.5)
    assert probs.shape == counts.shape
    assert np.isclose(probs.sum(), 1.0)
    assert np.all(probs >= 0)


def test_dirichlet_posterior_mean_handles_zero_counts():
    counts = np.array([0, 10, 0, 5])
    probs = dirichlet_posterior_mean(counts, alpha=0.5)
    assert probs.shape == counts.shape
    assert np.isclose(probs.sum(), 1.0)
    assert np.all(probs >= 0)
    assert not np.any(np.isnan(probs))


def test_dirichlet_posterior_mean_invariant_to_permutation():
    counts = np.array([5, 0, 10, 2])
    perm = np.array([2, 0, 3, 1])
    probs = dirichlet_posterior_mean(counts, alpha=0.5)
    probs_perm = dirichlet_posterior_mean(counts[perm], alpha=0.5)
    assert np.allclose(probs[perm], probs_perm)


def test_dirichlet_posterior_mean_large_alpha_approaches_uniform():
    counts = np.array([100, 1, 0, 0])
    probs_small_alpha = dirichlet_posterior_mean(counts, alpha=0.1)
    probs_large_alpha = dirichlet_posterior_mean(counts, alpha=100.0)
    uniform = np.full_like(counts, 1.0 / counts.size, dtype=float)
    # with large alpha, posterior should be closer to uniform than with small alpha
    # measure distance to uniform via max absolute difference
    dist_small = np.max(np.abs(probs_small_alpha - uniform))
    dist_large = np.max(np.abs(probs_large_alpha - uniform))
    assert dist_large < dist_small
    # and clearly less peaked than with small alpha
    assert probs_small_alpha.max() - probs_small_alpha.min() > probs_large_alpha.max() - probs_large_alpha.min()


def test_topological_charge_zero_for_pure_hexagons():
    counts = np.array([0, 0, 0, 100, 0])
    support = np.array([3, 4, 5, 6, 7])
    mean_q, mean_abs_q = topological_charge(counts, support)
    assert np.isclose(mean_q, 0.0)
    assert np.isclose(mean_abs_q, 0.0)


def test_topological_charge_heptagon_rich_negative_sign():
    support = np.array([3, 4, 5, 6, 7])
    counts = np.array([0, 0, 0, 10, 90])  # mostly heptagons
    mean_q, mean_abs_q = topological_charge(counts, support)
    assert mean_q < 0.0
    assert mean_abs_q > 0.0


def test_topological_charge_empty_counts_returns_zero_tuple():
    support = np.array([3, 4, 5, 6, 7])
    counts = np.array([0, 0, 0, 0, 0])
    mean_q, mean_abs_q = topological_charge(counts, support)
    assert np.isclose(mean_q, 0.0)
    assert np.isclose(mean_abs_q, 0.0)


def test_divergence_zero_for_identical_distributions():
    p = np.array([0.1, 0.2, 0.3, 0.3, 0.1])
    q = p.copy()
    support = np.array([3, 4, 5, 6, 7])
    jsd = js_divergence(p, q, base=2.0)
    w1 = wasserstein1_discrete(p, q, support)
    assert np.isclose(jsd, 0.0)
    assert np.isclose(w1, 0.0)


def test_js_divergence_symmetry_and_scaling_invariance():
    p = np.array([0.2, 0.3, 0.5])
    q = np.array([0.1, 0.4, 0.5])
    d_pq = js_divergence(p, q)
    d_qp = js_divergence(q, p)
    assert np.isclose(d_pq, d_qp)
    # scale both distributions by a constant; internal normalisation should make this a no-op
    scale = 10.0
    d_scaled = js_divergence(p * scale, q * scale)
    assert np.isclose(d_pq, d_scaled)


def test_js_divergence_positive_for_disjoint_support_masses():
    p = np.array([1.0, 0.0, 0.0])
    q = np.array([0.0, 0.0, 1.0])
    d = js_divergence(p, q)
    assert d > 0.0


def test_wasserstein1_discrete_simple_shift():
    p = np.array([1.0, 0.0])
    q = np.array([0.0, 1.0])
    support = np.array([0.0, 1.0])
    d = wasserstein1_discrete(p, q, support)
    assert np.isclose(d, 1.0)


def test_wasserstein1_discrete_nonuniform_support():
    p = np.array([0.0, 1.0, 0.0])
    q = np.array([0.0, 0.0, 1.0])
    support = np.array([0.0, 0.1, 1.0])
    d = wasserstein1_discrete(p, q, support)
    # all mass moves from 0.1 to 1.0
    assert np.isclose(d, 0.9)
