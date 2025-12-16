# Code for “Planck’s Law from a Classical Free Energy Extremum Involving Fisher Information”
Carlos Gomez-Uribe — accepted at *Quantum Studies: Mathematics and Foundations*  
Preprint: arXiv:2506.00586

This repository contains code to reproduce the **time-domain simulation and photon-counting diagnostics**
described in **Appendix D** of the paper, which provides an explicit stochastic realization of the
**threshold-activated emission cascade** mechanism introduced in **Appendix A**.

The core message of the paper is that the **Planck factor** \((e^\gamma-1)^{-1}\) can emerge without
invoking quantized oscillator energy levels or Bose counting, using (i) a classical **variational** free
energy functional over continuous densities, and (ii) a complementary **kinetic** cascade picture.
Here \(\gamma \equiv \hbar\omega/(k_B T)\) is the key dimensionless ratio.

---

## What the paper proves (high level)

### Main variational result (paper body + Appendices B–C)
A generalized free-energy functional over configuration-space densities,
containing Shannon entropy and Fisher information with \(\gamma\)-dependent weights, has an extremum (under a
Gaussian ansatz for a harmonic oscillator) whose mean energy reproduces Planck’s law. The *only* explicitly
quantum input is the **threshold assumption**: emission at frequency \(\omega\) requires a minimum energy
transfer \(\hbar\omega\).

For a single oscillator mode of frequency \(\omega\),
the radiative mean energy per mode is recovered as
\[
\langle E_{\mathrm{rad}}\rangle = \frac{\hbar\omega}{e^\gamma - 1}.
\]

### Complementary kinetic result (Appendix A)
Thermal “kicks” arrive with energies \(E\) drawn from an exponential (Boltzmann) distribution with mean
\(k_B T\). The probability a kick exceeds the emission threshold \(\hbar\omega\) is
\[
p_{\mathrm{th}} = \Pr(E>\hbar\omega) = e^{-\gamma}.
\]
A supra-threshold kick initiates an **emission cascade** (burst): emissions continue while subsequent
interaction opportunities remain supra-threshold, ending at the first sub-threshold event.
This yields the Planck factor via classical steady-state balance.

---

## What this repository implements

This code implements the **Appendix D** time-resolved model (a minimal single-mode cavity realization) and
computes standard photon-counting diagnostics from simulated time tags:

- \(g^{(2)}(\tau)\) (second-order coherence, HBT bunching signature)
- \(g^{(2)}(0)\) vs bin width \(\Delta t\)
- Fano factor \(F(\Delta t)\)

The simulation is a **Poisson cluster process**:

1. **Parent process (kicks):** kick times form a Poisson process of rate \(r\).
2. **Burst size per kick:** each kick injects an integer number of quanta \(N\) with a geometric law
   determined by \(\gamma\).
3. **Cavity leakage:** each injected photon exits after an independent exponential delay with rate \(\kappa\).
4. **Detection:** thinning with overall efficiency \(\eta\) (optional; statistics below are for detected tags).

For comparison, the code can also simulate a **matched-mean Poisson** photon stream.

---

## Model details and formulas (Appendix D)

### Geometric burst-size law from thresholded kick energy
Appendix D uses the mapping
\[
N = \left\lfloor \frac{E}{\hbar\omega} \right\rfloor,\qquad
E \sim \text{Exponential}(\text{mean }k_B T).
\]
Then
\[
\Pr(N\ge n) = \Pr(E\ge n\hbar\omega)=e^{-n\gamma},\qquad
\Pr(N=n) = (1-e^{-\gamma})e^{-n\gamma},\quad n=0,1,2,\dots
\]
and the mean burst size is
\[
\langle N\rangle = \frac{1}{e^\gamma-1}.
\]
Conditioning on \(N\ge1\) gives the “cascade size” distribution on \(\{1,2,\dots\}\) with mean
\(\langle N\mid N\ge1\rangle = 1/(1-e^{-\gamma})\).

### Time-domain emission / leakage model
Given a kick at time \(t_0\), each of the \(N\) photons exits after an independent delay
\(\Delta \sim \text{Exponential}(\kappa)\), so the photon time tags are \(t_0+\Delta_i\).

A “thermal calibration” often used in the Appendix D plots is \(r=\kappa\), yielding the standard
single-mode thermal bunching baseline.

---

## Diagnostics computed from time tags

Let \(\{t_i\}\) be detected photon times. Choose a bin width \(\Delta t\), and define binned counts
\[
n_j = \#\{i: t_i\in [j\Delta t,(j+1)\Delta t)\}.
\]

### Binned \(g^{(2)}(\tau)\)
For lag \(\tau = k\Delta t\) with \(k\ge1\),
\[
g^{(2)}(k\Delta t) = \frac{\langle n_j\,n_{j+k}\rangle}{\langle n_j\rangle^2}.
\]

### Zero-delay \(g^{(2)}(0)\) (factorial moment form)
\[
g^{(2)}(0)=\frac{\langle n_j(n_j-1)\rangle}{\langle n_j\rangle^2}
=1+\frac{\mathrm{Var}(n_j)-\langle n_j\rangle}{\langle n_j\rangle^2}.
\]

### Fano factor
\[
F(\Delta t) = \frac{\mathrm{Var}(n_j)}{\langle n_j\rangle}.
\]
For a matched-mean Poisson stream, \(g^{(2)}(\tau)=1\) and \(F=1\) for all \(\Delta t\).

---

## Analytic expectations for the cluster realization (Appendix D)

For the Poisson cluster process described above (Poisson parents of rate \(r\), exponential leakage with
rate \(\kappa\)), the second-order coherence has an exponential bunching form:
\[
g^{(2)}(\tau)=1+\frac{\kappa}{r}\,e^{-\kappa|\tau|}.
\]
With the thermal calibration \(r=\kappa\), this reduces to
\[
g^{(2)}(\tau)=1+e^{-\kappa|\tau|},
\]
the standard single-mode thermal (HBT) form.

**Note:** Different micro-realizations can share the same geometric burst law (hence the same Planck mean)
while differing in short-\(\tau\) structure; Appendix D chooses the minimal cavity-consistent model where
\(\kappa\) sets the only correlation time.

---

## Repository layout (edit filenames to match this repo)

- `src/simulate_cluster.py` — simulate the Poisson cluster process (kicks → burst sizes → leakage → tags)
- `src/simulate_poisson.py` — matched-mean Poisson baseline
- `src/diagnostics.py` — compute binned \(g^{(2)}(\tau)\), \(g^{(2)}(0)\), Fano
- `scripts/run_appendixD.py` — generate the Appendix D plots
- `figs/` — generated figures (optional; can be excluded from the repo)
- `requirements.txt` — Python dependencies

---

## Installation

### Python
Tested with Python 3.10+.

Install dependencies:
bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

## Reproducing Appendix D figures

Run:
bash
python scripts/run_appendixD.py

If you use this code, please cite:

Carlos Gomez-Uribe,
“Planck’s Law from a Classical Free Energy Extremum Involving Fisher Information,”
Quantum Studies: Mathematics and Foundations (2026),
arXiv:2506.00586
