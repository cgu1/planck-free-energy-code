Code for "Planck’s Law from a Classical Free Energy Extremum Involving Fisher Information"

Carlos Gomez-Uribe — accepted at Quantum Studies: Mathematics and Foundations
Preprint: arXiv:2506.00586

This repository contains a single Python script to reproduce the time-domain simulation and photon-counting diagnostics described in Appendix D of the paper. Appendix D provides an explicit stochastic time-tag realization of the threshold-activated emission-cascade mechanism introduced in Appendix A.

The paper’s core message is that the Planck factor 1/(exp(gamma)-1) can emerge without invoking quantized oscillator energy levels or Bose counting, using:

a classical variational free-energy functional over continuous densities, and

a complementary kinetic cascade picture.

Here

gamma = hbar*omega/(k_B*T)

is the key dimensionless ratio.

What the paper proves (high level)
Main variational result (paper body + Appendices B–C)

A generalized free-energy functional over configuration-space densities, containing Shannon entropy and Fisher information with gamma-dependent weights, has an extremum (under a Gaussian ansatz for a harmonic oscillator) whose mean energy reproduces Planck’s law. The only explicitly quantum input is a threshold assumption: emission at frequency omega requires a minimum energy transfer hbar*omega.

For a single oscillator mode of frequency omega, the radiative mean energy per mode is

<E_rad> = (hbar*omega)/(exp(gamma) - 1).

Complementary kinetic result (Appendix A)

Thermal “kicks” arrive with energies E drawn from an exponential (Boltzmann) distribution with mean k_B*T. The probability a kick exceeds the emission threshold hbar*omega is

p_th = P(E > hbar*omega) = exp(-gamma).

A supra-threshold kick initiates an emission cascade (burst): emissions continue while subsequent interaction opportunities remain supra-threshold, ending at the first sub-threshold event. This yields the Planck factor via classical steady-state balance.

What this repository implements

This repo implements the Appendix D time-resolved / time-tag model and computes standard photon-counting diagnostics from simulated time tags:

g2(tau) (second-order coherence / HBT bunching signature)

binned g2_dt(0) vs bin width dt

Fano factor F(dt)

It also simulates a matched-mean Poisson photon stream as a baseline comparison.

Quick start

Tested with Python 3.10+ (should work on 3.9+).

Create a virtual environment and install dependencies:

python -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib


Generate the Appendix D figures:

python make_appendixD_figures.py --outdir figs


This writes three PDFs to figs/ by default:

appendixD_g2_tau.pdf

appendixD_g2_vs_dt.pdf

appendixD_fano_vs_dt.pdf

Model implemented in Appendix D (Poisson cluster / cascade time tags)

This is a Poisson cluster process:

Parent process (“kicks”)
Kick times S_i form a homogeneous Poisson process with rate r.

Burst size per kick
Each kick produces an integer number of emitted quanta N_i with a geometric law determined by gamma.

Emission/leakage delays
Conditional on N_i, the emission delays are i.i.d. exponential:

D_ij ~ Exp(rate = kappa) for j = 1,...,N_i

Photon time tags
Photon time tags are

t_ij = S_i + D_ij.

Cascades overlap naturally: a new kick can occur before earlier emissions have all occurred, since the observed photon stream is the superposition of all clusters.

This script assumes unit detection efficiency (perfect detection). Independent thinning can be added later if desired; normalized quantities like g2(tau) are unchanged by uniform thinning, while raw count rates rescale.

Burst-size law from thresholded kick energy

Appendix D uses the mapping

N = floor( E / (hbar*omega) )

E ~ Exponential(mean = k_B*T).

Let gamma = hbar*omega/(k_B*T) and q = exp(-gamma). Then N has a geometric distribution on {0,1,2,...}:

P(N >= n) = exp(-n*gamma)

P(N = n) = (1-q) q^n = (1-exp(-gamma)) exp(-n*gamma), for n = 0,1,2,...

and

<N> = 1/(exp(gamma)-1).

The mean photon flux (intensity) of the cluster process is

nu = <I(t)> = r * <N> = r/(exp(gamma)-1).

The Poisson baseline is simulated with the same mean photon rate nu.

Diagnostics computed from time tags

Let {t_i} be photon times in an observation window [0,T]. Choose a bin width dt and define binned counts

n_j = #{ i : t_i in [j*dt, (j+1)*dt) }.

Binned g2(tau)

For lag tau = k*dt with k >= 1:

g2(k*dt) = < n_j * n_{j+k} > / <n_j>^2.

Zero-delay g2_dt(0) (factorial moment form)

g2_dt(0) = < n_j (n_j - 1) > / <n_j>^2

equivalently g2_dt(0) = 1 + (Var(n_j) - <n_j>)/<n_j>^2.

Fano factor

F(dt) = Var(n_j) / <n_j>.

For a matched-mean Poisson stream, g2(tau) = 1 and F(dt) = 1 for all dt.

Analytic expectations used in Appendix D

For this Poisson cluster model (Poisson parents of rate r, i.i.d. exponential delays with rate kappa, and geometric burst sizes as above), the normalized second-order coherence is

g2(tau) = 1 + (kappa/r) * exp(-kappa*|tau|).

A common “thermal calibration” in this minimal model is r = kappa, which gives

g2(0) = 2 and g2(tau) = 1 + exp(-kappa*|tau|),

i.e. the standard single-mode thermal (HBT) bunching form.

For finite bin width dt, the analytic bin-averaged zero-delay value is

g2_dt(0) = 1 + (2/r) * [ dt - (1 - exp(-kappa*dt))/kappa ] / dt^2.

Using the identity F(dt) = 1 + <n_j> * (g2_dt(0) - 1) with <n_j> = nu*dt, we get the closed form

F(dt) = 1 + (2*nu/r) * [ 1 - (1 - exp(-kappa*dt))/(kappa*dt) ],

nu = r/(exp(gamma)-1).

These are the analytic dashed curves plotted by the script.

Script outputs and options

By default:

the script chooses an observation window length T so that E[C] = r*T equals --target_cascades (default: 1,000,000),

it simulates cluster photons using a padded kick window [-pad_mult/kappa, T] (default pad_mult = 10) to better approximate stationarity near t = 0,

it generates the three PDFs listed above.

Run:

python make_appendixD_figures.py --outdir figs

If you use this code, please cite

Carlos Gomez-Uribe,
“Planck’s Law from a Classical Free Energy Extremum Involving Fisher Information,”
Quantum Studies: Mathematics and Foundations,
arXiv:2506.00586
