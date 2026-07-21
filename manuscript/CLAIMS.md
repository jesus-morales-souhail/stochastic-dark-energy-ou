# Referee checklist — claims vs non-claims

| | |
|:--|:--|
| **Manuscript** | [`PREPRINT.md`](PREPRINT.md) |
| **Author** | Jesús Morales Souhail · ORCID [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818) |
| **Code** | [stochastic-dark-energy-ou](https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou) |

---

## Claims (supported by public pipeline)

| ID | Claim | Evidence |
|:---|:------|:---------|
| **C1** | No preference for stationary OU/QNM residual over smooth BAO residuals | MLE floor; $\Delta\ln\mathcal{L}=0$ vs $\Lambda$CDM; AIC not preferred |
| **C2** | $\sigma_X < 1.5\times 10^{-4}$ (95% CL) under OU kernel + public BAO summary stats | Profile / residual pipeline in repo |
| **C3** | Nested free background + stochastic extension not preferred ($\Delta\mathrm{AIC}\approx +4$) | `eos_efectiva.py` / joint fits |
| **C4** | Coherent tachyonic growth (rank-1 cov) excluded: $\Delta\ln\mathcal{L}\approx -11.35$ | Quantum-fluid MLE |
| **C5** | Gap $\sigma_0\sim 10^{-61}\to\sim 10^{-5}$ is $\sim 10^{56}$; linear amplifiers fail audit | `amplifier_audit.py` + PREPRINT §5 |
| **C6** | Late freeze-out: gain $=1$; soft avalanche gain $\sim 2$; only $N_{\mathrm{eff}}$ redefinition hits amplitude | [`VERDICT.md`](../results/amplification_routes/VERDICT.md) |
| **C7** | Even $\varepsilon=1$ anisotropic fraction of BAO-bounded residual gives $\|\gamma-1\|\sim 10^{-4}$ | `scripts/slip_bridge.py` |

---

## Explicit non-claims

| ID | Non-claim |
|:---|:----------|
| **N1** | No detection of Planck-scale vacuum noise |
| **N2** | No claim that Euclid “will see” unamplified Sorkin seeds |
| **N3** | No unique proof of unimodular / SDiff gravity |
| **N4** | No homemade Boltzmann code results for $\eta(a,k)$ |
| **N5** | Cross-correlation DESI×SN is **preliminary**, not a primary result |
| **N6** | Optical tesseract / pupil / wavefront lab metaphors are **out of claim set** ([exploratory repo](https://github.com/jesus-morales-souhail/stochastic-de-exploratory-notes)) |
| **N7** | Not peer reviewed |

---

## Dataset boundary

- **In:** public DESI DR2 BAO summary statistics (baseline: diagonal measurement errors).
- **Out of primary claim:** full multi-probe DESI+CMB+SN dynamical DE preference (cited as context only).

---

## Reproducibility commands

```bash
python scripts/ou_bao_likelihood.py
python scripts/profile_sigma_x_desi.py
python scripts/amplifier_audit.py
python scripts/amplification/run_all_routes.py
python scripts/slip_bridge.py
```
