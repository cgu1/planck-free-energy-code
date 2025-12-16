#!/usr/bin/env python3
"""
cascade_model.py

Minimal Poisson-cluster realization for Appendix D (cavity mode).

Model:
  - Parent (kick) events: Poisson process of rate r_kick.
  - Each parent injects an integer burst N with geometric law:
        P(N=n) = (1-p) p^n,  n=0,1,2,... where p = exp(-gamma).
    Then E[N] = p/(1-p) = 1/(exp(gamma)-1).
  - Each of the N photons is emitted after an independent exponential delay with rate kappa.
  - Detection: independent thinning with efficiency eta.

This produces photon time tags suitable for g2 / Fano estimation.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class CascadeParams:
    gamma: float   # ħω/(k_B T)
    r_kick: float  # parent (kick) rate
    kappa: float   # cavity decay rate (sets correlation time ~ 1/kappa)
    eta: float = 1.0  # detection efficiency


def mean_burst_size(gamma: float) -> float:
    """E[N] for the geometric law on {0,1,2,...} with p=exp(-gamma)."""
    p = float(np.exp(-gamma))
    if p >= 1.0:
        return float("inf")
    return p / (1.0 - p)


def sample_burst_sizes(rng: np.random.Generator, gamma: float, size: int) -> np.ndarray:
    """
    Sample N~Geom0 with P(N=n)=(1-p)p^n, n>=0.
    If U~Unif(0,1), N = floor(log(U)/log(p)).
    """
    p = float(np.exp(-gamma))
    if p <= 0.0:
        return np.zeros(size, dtype=int)
    u = rng.random(size=size)
    # avoid log(0)
    u = np.clip(u, 1e-15, 1.0 - 1e-15)
    n = np.floor(np.log(u) / np.log(p)).astype(int)
    n[n < 0] = 0
    return n


def simulate_poisson_parents(T_total: float, rate: float, rng: np.random.Generator) -> np.ndarray:
    """Simulate parent event times for a Poisson process on [0,T_total)."""
    if rate <= 0 or T_total <= 0:
        return np.array([], dtype=float)
    # simulate count then uniform times (equivalent for homogeneous Poisson)
    K = rng.poisson(rate * T_total)
    if K <= 0:
        return np.array([], dtype=float)
    times = rng.random(K) * T_total
    times.sort()
    return times


def simulate_cascade_process(T_total: float, params: CascadeParams, seed: int = 0) -> np.ndarray:
    """
    Return detected photon times (sorted) on [0,T_total).
    """
    rng = np.random.default_rng(seed)

    parent_times = simulate_poisson_parents(T_total, params.r_kick, rng)
    if parent_times.size == 0:
        return np.array([], dtype=float)

    Ns = sample_burst_sizes(rng, params.gamma, parent_times.size)
    total_photons = int(Ns.sum())
    if total_photons == 0:
        return np.array([], dtype=float)

    # For each parent i, generate Ns[i] exponential delays and add to parent time
    # Vectorize by expanding parent times
    parents_rep = np.repeat(parent_times, Ns)
    delays = rng.exponential(scale=1.0 / params.kappa, size=total_photons)
    photon_times = parents_rep + delays

    # Keep within window
    photon_times = photon_times[(photon_times >= 0.0) & (photon_times < T_total)]
    if photon_times.size == 0:
        return photon_times

    # Detection thinning
    if params.eta < 1.0:
        keep = rng.random(photon_times.size) < params.eta
        photon_times = photon_times[keep]

    photon_times.sort()
    return photon_times


def simulate_poisson_photons(T_total: float, rate: float, seed: int = 0) -> np.ndarray:
    """Direct Poisson photon stream on [0,T_total)."""
    rng = np.random.default_rng(seed)
    return simulate_poisson_parents(T_total, rate, rng)


def binned_counts(times: np.ndarray, T_total: float, dt: float) -> np.ndarray:
    """Histogram counts into bins of width dt over [0,T_total)."""
    if dt <= 0 or T_total <= 0:
        return np.array([], dtype=int)
    nbins = int(np.ceil(T_total / dt))
    edges = np.linspace(0.0, nbins * dt, nbins + 1)
    counts, _ = np.histogram(times, bins=edges)
    return counts.astype(int)


def g2_from_counts(counts: np.ndarray, max_lag_bins: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Estimate binned g2(k) for k=0..max_lag_bins.
    Uses factorial moment for k=0: <n(n-1)>/<n>^2.
    For k>=1: <n_j n_{j+k}>/<n>^2.
    """
    counts = np.asarray(counts, dtype=float)
    m = counts.mean()
    if m <= 0 or counts.size == 0:
        lags = np.arange(max_lag_bins + 1, dtype=int)
        return lags, np.full_like(lags, np.nan, dtype=float)

    lags = np.arange(max_lag_bins + 1, dtype=int)
    g2 = np.empty_like(lags, dtype=float)

    # k=0
    g2[0] = np.mean(counts * (counts - 1.0)) / (m * m)

    # k>=1
    for k in range(1, max_lag_bins + 1):
        a = counts[:-k]
        b = counts[k:]
        g2[k] = np.mean(a * b) / (m * m)

    return lags, g2


def g2_0_from_times(times: np.ndarray, T_total: float, dt: float) -> float:
    c = binned_counts(times, T_total, dt).astype(float)
    m = c.mean()
    if m <= 0 or c.size == 0:
        return float("nan")
    return float(np.mean(c * (c - 1.0)) / (m * m))


def fano_from_times(times: np.ndarray, T_total: float, dt: float) -> float:
    c = binned_counts(times, T_total, dt).astype(float)
    m = c.mean()
    if m <= 0 or c.size == 0:
        return float("nan")
    return float(c.var(ddof=0) / m)
