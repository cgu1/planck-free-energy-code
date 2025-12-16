#!/usr/bin/env python3
"""
make_appendixD_figures.py

Generates Appendix D figures into figs/.

This version uses the equilibrium calibration r_kick = kappa so that the
continuous-time g2 peak satisfies g2(0)=2 in the high-time-resolution limit.
"""

from __future__ import annotations

import os, csv
import numpy as np
import matplotlib.pyplot as plt

from cascade_model import (
    CascadeParams,
    simulate_cascade_process,
    mean_burst_size,
    simulate_poisson_photons,
    binned_counts,
    g2_from_counts,
    g2_0_from_times,
    fano_from_times,
)

OUTPUT_DIR = "figs"

# --- Simulation length (increase for smooth curves) ---
T_total = 300000.0

# --- Cavity parameters (dimensionless time units) ---
gamma_demo = 1.5
kappa = 2.0            # cavity linewidth
r_kick = kappa         # equilibrium calibration: r_kick = kappa  -> g2(0)->2 as dt->0
eta = 1.0

# --- g2(tau) estimation settings ---
dt_g2 = 0.25
max_tau = 10.0

# --- dt sweep ---
dt_values = np.logspace(np.log10(0.05), np.log10(20.0), 28)

# --- gamma sweep ---
gammas = np.linspace(0.4, 3.0, 14)
dt_sweep = 1.0 / kappa  # order of correlation time

SEED = 1


def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)


def outpath(name: str) -> str:
    return os.path.join(OUTPUT_DIR, name)


# -------- Analytic expectations for this cluster model --------
def analytic_g2_tau(tau: np.ndarray, r: float, kappa: float) -> np.ndarray:
    return 1.0 + (kappa / r) * np.exp(-kappa * np.abs(tau))


def analytic_g2_binned0(dt: np.ndarray, r: float, kappa: float) -> np.ndarray:
    dt = np.asarray(dt, dtype=float)
    num = dt - (1.0 - np.exp(-kappa * dt)) / kappa
    return 1.0 + (2.0 / r) * (num / (dt ** 2))


def analytic_g2_binned_lags(max_lag_bins: int, dt: float, r: float, kappa: float) -> np.ndarray:
    g = np.empty(max_lag_bins + 1, dtype=float)
    g[0] = analytic_g2_binned0(np.array([dt]), r=r, kappa=kappa)[0]
    if max_lag_bins >= 1:
        k = np.arange(1, max_lag_bins + 1, dtype=float)
        Ik = np.exp(-kappa * (k - 1.0) * dt) * (1.0 - np.exp(-kappa * dt)) ** 2 / (kappa ** 2)
        g[1:] = 1.0 + (kappa / r) * (Ik / (dt ** 2))
    return g


def analytic_fano(dt: np.ndarray, lam: float, g2bin0: np.ndarray) -> np.ndarray:
    dt = np.asarray(dt, dtype=float)
    mu = lam * dt
    return 1.0 + mu * (g2bin0 - 1.0)


def main() -> None:
    ensure_dir(OUTPUT_DIR)

    params = CascadeParams(gamma=gamma_demo, r_kick=r_kick, kappa=kappa, eta=eta)
    times = simulate_cascade_process(T_total, params, seed=SEED)

    mean_rate_emp = len(times) / T_total if T_total > 0 else 0.0
    times_p = simulate_poisson_photons(T_total, mean_rate_emp, seed=SEED + 1)

    # Theoretical mean photon rate for our model
    lam_theory = r_kick * mean_burst_size(gamma_demo) * eta

    # --- Fig 1: g2(tau) ---
    max_lag = int(max_tau / dt_g2)
    c = binned_counts(times, T_total, dt_g2)
    cp = binned_counts(times_p, T_total, dt_g2)

    lags, g2c = g2_from_counts(c, max_lag_bins=max_lag)
    _, g2p = g2_from_counts(cp, max_lag_bins=max_lag)

    tau = lags * dt_g2
    g2c_ana = analytic_g2_binned_lags(max_lag, dt_g2, r=r_kick, kappa=kappa)
    g2p_ana = np.ones_like(g2c_ana)

    plt.figure()
    plt.plot(tau, g2c, marker="o", label="Cascade simulation")
    plt.plot(tau, g2c_ana, linestyle="--", label="Cascade analytic (binned)")
    plt.plot(tau, g2p, marker="o", label="Poisson simulation (matched mean)")
    plt.plot(tau, g2p_ana, linestyle="--", label="Poisson analytic")
    plt.xlabel(r"Delay $\tau$")
    plt.ylabel(r"$g^{(2)}(\tau)$ (binned)")
    plt.xlim(0, max_tau)
    #ymax = max(2.5, np.nanmax(g2c_ana) * 1.10, np.nanmax(g2c) * 1.10)
    plt.ylim(0.75, 2.1)
    plt.legend()
    plt.savefig(outpath("appendixD_g2_tau.pdf"), bbox_inches="tight")

    # --- Fig 2: g2(0) vs dt ---
    g2c0 = np.array([g2_0_from_times(times, T_total, float(dtv)) for dtv in dt_values])
    g2p0 = np.array([g2_0_from_times(times_p, T_total, float(dtv)) for dtv in dt_values])
    g2c0_ana = analytic_g2_binned0(dt_values, r=r_kick, kappa=kappa)
    g2p0_ana = np.ones_like(dt_values)

    plt.figure()
    plt.plot(dt_values, g2c0, marker="o", label="Cascade simulation")
    plt.plot(dt_values, g2c0_ana, linestyle="--", label="Cascade analytic")
    plt.plot(dt_values, g2p0, marker="o", label="Poisson simulation (matched mean)")
    plt.plot(dt_values, g2p0_ana, linestyle="--", label="Poisson analytic")
    plt.xscale("log")
    plt.xlabel(r"Bin width $\Delta t$")
    plt.ylabel(r"$g^{(2)}(0)$ (binned)")
    plt.ylim(0.8, max(2.25, np.nanmax(g2c0_ana) * 1.10))
    plt.legend()
    plt.savefig(outpath("appendixD_g2_vs_dt.pdf"), bbox_inches="tight")

    # --- Fig 3: Fano vs dt ---
    fano_c = np.array([fano_from_times(times, T_total, float(dtv)) for dtv in dt_values])
    fano_p = np.array([fano_from_times(times_p, T_total, float(dtv)) for dtv in dt_values])
    fano_c_ana = analytic_fano(dt_values, lam=lam_theory, g2bin0=g2c0_ana)
    fano_p_ana = np.ones_like(dt_values)

    plt.figure()
    plt.plot(dt_values, fano_c, marker="o", label="Cascade simulation")
    plt.plot(dt_values, fano_c_ana, linestyle="--", label="Cascade analytic")
    plt.plot(dt_values, fano_p, marker="o", label="Poisson simulation (matched mean)")
    plt.plot(dt_values, fano_p_ana, linestyle="--", label="Poisson analytic")
    plt.axhline(1.0, linestyle=":", label="F=1")
    plt.xscale("log")
    plt.xlabel(r"Bin width $\Delta t$")
    plt.ylabel("Fano factor Var(N)/E(N)")
    plt.ylim(0.8, max(1.5, np.nanmax(fano_c_ana) * 1.10))
    plt.legend()
    plt.savefig(outpath("appendixD_fano_vs_dt.pdf"), bbox_inches="tight")

    # --- Gamma sweep ---
    g2_sweep = []
    fano_sweep = []
    lam_sweep = []
    for i, g in enumerate(gammas):
        p = CascadeParams(gamma=float(g), r_kick=r_kick, kappa=kappa, eta=eta)
        t = simulate_cascade_process(T_total, p, seed=1000 + i)
        g2_sweep.append(g2_0_from_times(t, T_total, dt_sweep))
        fano_sweep.append(fano_from_times(t, T_total, dt_sweep))
        lam_sweep.append(r_kick * mean_burst_size(float(g)) * eta)

    g2_sweep = np.asarray(g2_sweep)
    fano_sweep = np.asarray(fano_sweep)
    lam_sweep = np.asarray(lam_sweep)

    # Analytic: g2 depends only on r,kappa,dt (not on gamma) in this minimal model; Fano depends via mean rate
    g2_sweep_ana = np.full_like(gammas, analytic_g2_binned0(np.array([dt_sweep]), r=r_kick, kappa=kappa)[0])
    fano_sweep_ana = analytic_fano(np.full_like(gammas, dt_sweep), lam=lam_sweep, g2bin0=g2_sweep_ana)

    plt.figure()
    plt.plot(gammas, g2_sweep, marker="o", label="Cascade simulation")
    plt.plot(gammas, g2_sweep_ana, linestyle="--", label="Cascade analytic")
    plt.axhline(1.0, linestyle="--", label="Poisson analytic")
    plt.xlabel(r"$\gamma=\hbar\omega/(k_B T)$")
    plt.ylabel(r"$g^{(2)}(0)$ at fixed $\Delta t$")
    plt.ylim(0.8, max(3.0, np.nanmax(g2_sweep_ana) * 1.10, np.nanmax(g2_sweep) * 1.10))
    plt.legend()
    plt.savefig(outpath("appendixD_g2_vs_gamma.pdf"), bbox_inches="tight")

    plt.figure()
    plt.plot(gammas, fano_sweep, marker="o", label="Cascade simulation")
    plt.plot(gammas, fano_sweep_ana, linestyle="--", label="Cascade analytic")
    plt.axhline(1.0, linestyle="--", label="Poisson analytic")
    plt.xlabel(r"$\gamma=\hbar\omega/(k_B T)$")
    plt.ylabel("Fano factor at fixed $\Delta t$")
    plt.ylim(0.8, max(1.5, np.nanmax(fano_sweep_ana) * 1.10, np.nanmax(fano_sweep) * 1.10))
    plt.legend()
    plt.savefig(outpath("appendixD_fano_vs_gamma.pdf"), bbox_inches="tight")

    with open(outpath("appendixD_gamma_sweep.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gamma", "g2_0_sim", "g2_0_ana", "fano_sim", "fano_ana", "mean_rate_theory", "dt_sweep"])
        for g, g2s, g2a, fs, fa, lam in zip(gammas, g2_sweep, g2_sweep_ana, fano_sweep, fano_sweep_ana, lam_sweep):
            w.writerow([float(g), float(g2s), float(g2a), float(fs), float(fa), float(lam), float(dt_sweep)])

    print("Generated figures in", OUTPUT_DIR)
    print("Thermal calibration: r_kick =", r_kick, "kappa =", kappa, "=> g2(0)~2 as dt->0")


if __name__ == "__main__":
    main()
