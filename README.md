# stochastic-dark-energy-ou

Jesús Morales Souhail  
[ORCID 0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818) · jmskjym@gmail.com  
[github.com/jesus-morales-souhail](https://github.com/jesus-morales-souhail)

July 2026 · independent research

I use this repository for the **data** side: public DESI DR2 BAO, residual models (OU / QNM), what is excluded, and how large the soft-amplification gap is if you start from a pure Sorkin seed.

Theory lives in [measurable-stochastic-vacuum](https://github.com/jesus-morales-souhail/measurable-stochastic-vacuum).  
Method / optics notes only: [stochastic-de-exploratory-notes](https://github.com/jesus-morales-souhail/stochastic-de-exploratory-notes).

---

## Read these first

1. [`manuscript/PREPRINT.md`](manuscript/PREPRINT.md) — main write-up  
2. [`manuscript/CLAIMS.md`](manuscript/CLAIMS.md) — what I claim and what I do not  
3. [`papers/amplification-gap.md`](papers/amplification-gap.md) — why soft \(10^{56}\) is not free  
4. [`papers/HONEST_HEADLINES.md`](papers/HONEST_HEADLINES.md) — wording I allow myself  
5. [`papers/EXPLORATORY_BOUNDARY.md`](papers/EXPLORATORY_BOUNDARY.md) — what is not part of the paper  
6. [`papers/scope-and-mixups.md`](papers/scope-and-mixups.md) — residual \(\sigma_X\) is not \(\omega_b\), not EHT  

---

## Results in one place

Under the OU/QNM residual kernel on public DESI DR2 BAO, a stochastic residual is not required. Profile 95% CL (free \(\theta\), 7 isotropic \(\alpha\) bins, full measurement covariance):

\[
\sigma_X < 2.5\times 10^{-2}.
\]

Best fit still goes to \(\sigma_X\to 0\). A coherent tachyonic growing mode is excluded when it is active (rank-1 covariance). Soft amplification of a pure Sorkin seed \(\sim 10^{-61}\) to Euclid-scale residuals needs \(\sim 10^{56}\) (to the DESI ceiling, \(\sim 10^{59}\)) and is not free — see the amplification note.

Main long notes: `papers/stochastic-dark-energy-desi-dr2.md`, `papers/quantum-fluid-instabilities-desi-dr2.md`.  
Geometry / unimodular discussion: `papers/unimodular-gravity-vacuum-smoothness.md` and the short principle `papers/principle-of-vacuum-smoothness.md`.  
Numbers snapshot: `papers/resume.txt`. Bound history: `papers/BOUND_CORRECTION_SIGMA_X.md`.

---

## Layout

```
manuscript/   preprint and claims
papers/       technical notes
scripts/      analysis code
data/         local DESI products (public provenance)
results/      numerical outputs
notes/        working notes
```

```bash
python scripts/profile_sigma_x_desi.py
python scripts/tachyonic_rank1_mle.py
python scripts/amplifier_audit.py
```
