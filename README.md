Code for "Planck's Law from a Classical Free Energy Extremum Involving Fisher Information"

Carlos Gomez-Uribe — accepted at Quantum Studies: Mathematics and Foundations
Preprint: arXiv:2506.00586

This repository contains code to reproduce the time-domain simulation and photon-counting diagnostics
described in Appendix D of the paper. Appendix D provides an explicit stochastic (time-tag) realization
of the threshold-activated emission-cascade mechanism introduced in Appendix A.

Core message (paper-level): the Planck factor (exp(gamma) - 1)^(-1) can emerge without invoking
quantized oscillator energy levels or Bose counting, using:
(1) a classical variational free-energy extremum over continuous densities, and
(2) a complementary kinetic cascade picture.
Here gamma = (hbar * omega) / (kB * T) is the key dimensionless ratio.

What this repository implements (Appendix D)

The code implements a minimal "time-tag" emission model and computes standard photon-counting diagnostics:

g2(tau): second-order coherence (HBT bunching signature)

g2_dt0(0): the zero-delay binned estimate versus bin width dt

F(dt): the Fano factor versus bin width dt

All figures are produced by a single self-contained script:

make_appendixD_figures.py

It writes three PDF files (by default into figs/):

appendixD_g2_tau.pdf

appendixD_g2_vs_dt.pdf

appendixD_fano_vs_dt.pdf

Model definition (Poisson cluster / cascade time-tag model)

This Appendix D model is a Poisson cluster process:

Parent process (kicks):
Kick times S_i form a Poisson process with rate r.

Burst size per kick:
Each kick i produces an integer number of photons N_i with a geometric law on {0,1,2,...}:

P(N = n) = (1 - q) * q^n, n = 0,1,2,...
q = exp(-gamma)

This implies:
mean(N) = q/(1-q) = 1/(exp(gamma) - 1)

(Conditioning on N >= 1 gives the burst size on {1,2,...} used in Appendix A.)

Emission / leakage kernel (time tags):
Conditional on N_i, each of the N_i photons receives an independent emission delay:

D_ij ~ Exponential(rate = kappa), j = 1,...,N_i

The photon time tags are:
t_ij = S_i + D_ij

Cascades can overlap in time: a new kick can occur before photons from earlier kicks have emitted.

Mean photon rate:
The mean photon flux (intensity) is:
nu = r * mean(N) = r/(exp(gamma) - 1)

Baseline comparison:

A Poisson photon process on the same window [0,T] with rate nu.

Diagnostics from time tags (binning definitions)

Given photon time tags {t_k}, choose a bin width dt and define binned counts:

n_j = number of events t_k in [j*dt, (j+1)*dt)

Binned g2 at nonzero lag:
For tau = kdt with k >= 1,
g2(kdt) = mean( n_j * n_{j+k} ) / mean(n_j)^2

Zero-delay binned g2 (factorial moment form):
g2_dt0(0) = mean( n_j * (n_j - 1) ) / mean(n_j)^2

Fano factor:
F(dt) = Var(n_j) / mean(n_j)

Poisson baseline:
For a Poisson process, g2(tau) = 1 and F(dt) = 1 for all dt.

Analytic expectations used for the Appendix D plots

For the Poisson cluster process above with exponential delays (rate kappa), the normalized
second-order coherence has an exponential bunching form:

g2(tau) = 1 + (kappa/r) * exp(-kappa * |tau|)

A commonly used "thermal calibration" in this minimal model is r = kappa, which gives:
g2(0) = 2
g2(tau) = 1 + exp(-kappa * |tau|)
(the standard single-mode thermal HBT form).

Binned zero-delay g2 (bin width dt):
g2_dt0(0) = 1 + (2/r) * [ dt - (1 - exp(-kappa*dt))/kappa ] / dt^2

Fano factor in closed form:
F(dt) = 1 + (2nu/r) * [ 1 - (1 - exp(-kappadt))/(kappa*dt) ]
where nu = r/(exp(gamma) - 1)

(These are exactly the analytic curves overlaid as dashed lines in the generated figures.)

Installation

Python:

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

python make_appendixD_figures.py --outdir figs

The script exposes command-line parameters for r, kappa, gamma, and simulation size (target number
of cascades), plus plotting controls for dt sweeps and g2(tau).

Citation

If you use this code, please cite:

Carlos Gomez-Uribe,
"Planck's Law from a Classical Free Energy Extremum Involving Fisher Information,"
Quantum Studies: Mathematics and Foundations (accepted),
arXiv:2506.00586
