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
2. [`manuscript/CLAIMS.md`](manuscript/CLAIMS.md) — what I claim (C1–C9) and what I do not (N1–N13)  
3. [`papers/amplification-gap.md`](papers/amplification-gap.md) — why soft \(10^{56}\) is not free  
4. [`papers/HONEST_HEADLINES.md`](papers/HONEST_HEADLINES.md) — wording I allow myself  
5. [`theory/EXPLORATORY_BOUNDARY.md`](theory/EXPLORATORY_BOUNDARY.md) — what is not part of the paper  
6. [`papers/scope-and-mixups.md`](papers/scope-and-mixups.md) — residual \(\sigma_X\) is not \(\omega_b\), not EHT  
7. [`analisis_cmb/RESULTADO.md`](analisis_cmb/RESULTADO.md) — CMB consistency check (Planck 2018)

---

## Results in one place

Under the OU/QNM residual kernel on public DESI DR2 BAO, a stochastic residual is not required. Profile 95% CL (free \(\theta\), 7 isotropic \(\alpha\) bins, full measurement covariance):

\[
\sigma_X < 2.6\times 10^{-2}.
\]

Best fit still goes to \(\sigma_X\to 0\). A coherent tachyonic growing mode is excluded when it is active (rank-1 covariance). Soft amplification of a pure horizon-area Sorkin seed \(\sim 10^{-61}\) (not Everpresent Λ, which is a different, order-one model) to Euclid-scale residuals needs \(\sim 10^{56}\) (to the DESI ceiling, \(\sim 10^{59}\)) and is not free — see the amplification note.

Planck 2018 (TT ℓ<30 + lensing) does not distinguish this residual from ΛCDM at the DESI ceiling (p95 Δχ² ≤ 0.067) and adds no bound; its H0 imprint, 0.48 km/s/Mpc, does not explain the H0 tension. A percent-level harmonic residual is equally allowed (58/60) and CMB-invisible. See `analisis_cmb/RESULTADO.md`.

Main long notes: [`papers/stochastic-dark-energy-desi-dr2.md`](papers/stochastic-dark-energy-desi-dr2.md), [`papers/quantum-fluid-instabilities-desi-dr2.md`](papers/quantum-fluid-instabilities-desi-dr2.md).  
Geometry / unimodular discussion (theory, outside the claim set): [`theory/unimodular-gravity-vacuum-smoothness.md`](theory/unimodular-gravity-vacuum-smoothness.md) and the short principle [`theory/principle-of-vacuum-smoothness.md`](theory/principle-of-vacuum-smoothness.md).  
Numbers snapshot: [`papers/resume.txt`](papers/resume.txt). Bound history: [`papers/BOUND_CORRECTION_SIGMA_X.md`](papers/BOUND_CORRECTION_SIGMA_X.md).

---

## Layout

What supports the preprint (claims C1–C7) lives in `manuscript/`, `papers/`, `scripts/`, `results/` and `figures/`. Everything else is labelled as theory, exploratory or archived.

```
manuscript/    preprint, claim checklist, cover letter
papers/        technical notes that support C1–C7
scripts/       analysis code for the claim set (DESI DR2 BAO, amplification, slip bridge)
results/       numerical outputs of scripts/
figures/       the two figures copied there by their scripts (profile_sigma_x, joint_w0wa_sigma)
theory/        theory-programme notes (unimodular / SDiff, vacuum smoothness, ν–DE, Euclid protocol); not C1–C7
exploratory/   desqueezing programme: scripts/, results/, notes/; not cited as DESI results
analisis_cmb/  OU at the DESI ceiling vs Planck 2018 (CMB), everpresent Λ, audit and correction logs
archive/       retired material (synthetic Euclid MCMC, stubs, DR1 script, orphan figures), kept for history
data/          local DESI products (public provenance; not versioned)
```

Each of `theory/`, `exploratory/`, `archive/` and `analisis_cmb/` has its own README.

---

## Reproduce

```bash
pip install -r requirements.txt
python scripts/profile_sigma_x_desi.py
python scripts/tachyonic_rank1_mle.py
python scripts/amplifier_audit.py
```

OU / QNM residual likelihood used in the preprint: `python scripts/ou_bao_likelihood.py`.  
Consistency check across the BAO pipelines: `python scripts/verify_bao_pipelines.py`.
