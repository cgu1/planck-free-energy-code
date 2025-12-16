Code for “Planck’s Law from a Classical Free Energy Extremum Involving Fisher Information”

Carlos Gomez-Uribe — accepted at Quantum Studies: Mathematics and Foundations
Preprint: arXiv:2506.00586

Overview

This repository contains code to reproduce the time-domain simulation and photon-counting diagnostics described in Appendix D of the paper. Appendix D provides an explicit stochastic (time-tag) realization of the threshold-activated emission cascade mechanism introduced in Appendix A.

Core message (paper-level):
The Planck factor (exp(γ) − 1)⁻¹ can emerge without invoking quantized oscillator energy levels or Bose counting, using (1) a classical variational free-energy extremum over continuous probability densities, and (2) a complementary kinetic cascade picture. Here γ = ℏω / (kBT) is the key dimensionless ratio controlling the crossover between thermal and quantum regimes.

What this repository implements (Appendix D)

The code implements a minimal time-tag emission model and computes standard photon-counting diagnostics:

• g²(τ): second-order coherence (Hanbury Brown–Twiss bunching signature)
• g²₀(Δt): zero-delay binned estimate versus bin width Δt
• F(Δt): Fano factor versus bin width Δt

All figures are produced by a single self-contained script:

make_appendixD_figures.py

By default, it writes three PDF files to figs/:

appendixD_g2_tau.pdf
appendixD_g2_vs_dt.pdf
appendixD_fano_vs_dt.pdf

Model definition (Poisson cluster / cascade time-tag model)

The Appendix D model is a Poisson cluster process.

Parent process (thermal kicks)

Kick times Sᵢ form a Poisson process with rate r.

Burst size per kick

Each kick i produces an integer number of photons Nᵢ with a geometric distribution on {0,1,2,…}:

P(N = n) = (1 − q) qⁿ,  n = 0,1,2,…
q = exp(−γ)

This implies:

mean(N) = q / (1 − q) = 1 / (exp(γ) − 1)

Conditioning on N ≥ 1 yields the burst-size distribution on {1,2,…} used in Appendix A.

Emission / leakage kernel (time tags)

Conditional on Nᵢ, each of the Nᵢ photons receives an independent emission delay:

Dᵢⱼ ~ Exponential(rate = κ),  j = 1,…,Nᵢ

The photon time tags are:

tᵢⱼ = Sᵢ + Dᵢⱼ

Cascades can overlap in time: a new kick may occur before photons from earlier kicks have emitted.

Mean photon rate

The mean photon flux (intensity) is:

ν = r · mean(N) = r / (exp(γ) − 1)

Baseline comparison

A Poisson photon process on the same observation window [0,T] with rate ν.

Diagnostics from time tags (binning definitions)

Given photon time tags {tₖ}, choose a bin width Δt and define binned counts:

nⱼ = number of events tₖ in [jΔt, (j+1)Δt)

Binned g² at nonzero lag

For τ = kΔt with k ≥ 1:

g²(kΔt) = mean(nⱼ · nⱼ₊ₖ) / mean(nⱼ)²

Zero-delay binned g² (factorial moment form)

g²₀(0) = mean(nⱼ (nⱼ − 1)) / mean(nⱼ)²

Fano factor

F(Δt) = Var(nⱼ) / mean(nⱼ)

Poisson baseline

For a Poisson process, g²(τ) = 1 and F(Δt) = 1 for all Δt.

Analytic expectations used for the Appendix D plots

For the Poisson cluster process above with exponential delays (rate κ), the normalized second-order coherence has an exponential bunching form:

g²(τ) = 1 + (κ / r) · exp(−κ |τ|)

A commonly used thermal calibration in this minimal model is r = κ, which yields:

g²(0) = 2
g²(τ) = 1 + exp(−κ |τ|)

This is the standard single-mode thermal HBT form.

Binned zero-delay g² (bin width Δt)

g²₀(0) = 1 + (2 / r) · [ Δt − (1 − exp(−κΔt)) / κ ] / Δt²

Fano factor (closed form)

F(Δt) = 1 + (2ν / r) · [ 1 − (1 − exp(−κΔt)) / (κΔt) ]

where ν = r / (exp(γ) − 1)

These analytic expressions are exactly the dashed curves overlaid on the simulation results in the figures.

Installation

Python

Tested with Python 3.10+ (should work with 3.9+).

Dependencies:

numpy
matplotlib

Example setup (optional virtual environment):

python -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib

Reproducing the Appendix D figures

Run (default parameters; writes PDFs to figs/):

python make_appendixD_figures.py –outdir figs

The script exposes command-line parameters for r, κ, γ, and simulation size (target number of cascades), as well as plotting controls for Δt sweeps and g²(τ).

Citation

If you use this code, please cite:

Carlos Gomez-Uribe,
“Planck’s Law from a Classical Free Energy Extremum Involving Fisher Information,”
Quantum Studies: Mathematics and Foundations (accepted),
arXiv:2506.00586
