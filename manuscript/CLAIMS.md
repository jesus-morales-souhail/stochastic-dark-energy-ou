# Referee checklist — claims vs non-claims

| | |
|:--|:--|
| **Manuscript** | [`PREPRINT.md`](PREPRINT.md) |
| **Author** | Jesús Morales Souhail · ORCID [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818) |
| **Code** | [stochastic-dark-energy-ou](https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou) |

This is what I stand behind, and what I do not.
---

## Claims (supported by the public pipeline)

| ID | Claim | Evidence |
|:---|:------|:---------|
| **C1** | No preference for a stationary OU/QNM residual over smooth BAO residuals | MLE floor; $\Delta\ln\mathcal{L}=0$ vs $\Lambda$CDM; AIC not preferred |
| **C2** | $\sigma_X < 1.5\times 10^{-4}$ (95% CL) under the OU kernel and public BAO summary statistics | Profile / residual pipeline in this repo |
| **C3** | Nested free background plus stochastic extension is not preferred ($\Delta\mathrm{AIC}\approx +4$) | `eos_efectiva.py` / joint fits |
| **C4** | Coherent tachyonic growth (rank-1 cov) is excluded: $\Delta\ln\mathcal{L}\approx -12.7$ | Quantum-fluid MLE |
| **C5** | Gap $\sigma_0\to 10^{-5}$ is $\sim 10^{56}$ (**Euclid target**); gap $\to 1.5\times 10^{-4}$ is $\sim 10^{57}$ (**DESI ceiling**); linear amplifiers fail both | `gap_two_targets.py` + `amplifier_audit.py` + PREPRINT §5 |
| **C6** | Late freeze-out: gain $=1$; soft avalanche gain $\sim 2$; only $N_{\mathrm{eff}}$ redefinition hits amplitude | [`VERDICT.md`](../results/amplification_routes/VERDICT.md) |
| **C7** | Even $\varepsilon=1$ anisotropic fraction of the BAO-bounded residual gives $\|\gamma-1\|\sim 10^{-4}$ | `scripts/slip_bridge.py` |

---

## Explicit non-claims

| ID | Non-claim |
|:---|:----------|
| **N1** | I do not claim detection of Planck-scale vacuum noise |
| **N2** | I do not claim that Euclid “will see” unamplified Sorkin seeds |
| **N3** | I do not claim a unique proof of unimodular / SDiff gravity |
| **N4** | I do not report homemade Boltzmann-code results for $\eta(a,k)$ |
| **N5** | The DESI×SN cross-correlation is **preliminary**, not a primary result |
| **N6** | Optical tesseract / pupil / wavefront lab metaphors are **outside the claim set** ([exploratory repo](https://github.com/jesus-morales-souhail/stochastic-de-exploratory-notes)) |
| **N7** | This is independent work |
| **N8** | I do not claim a bound on $\sum m_\nu$; neutrino mass is not in the likelihood. See [`papers/neutrino-de-degeneracy-and-sigma-x.md`](../papers/neutrino-de-degeneracy-and-sigma-x.md) |

---

## Dataset boundary

- **In:** public DESI DR2 BAO summary statistics (baseline: diagonal measurement errors).
- **Out of primary claim:** full multi-probe DESI+CMB+SN dynamical DE preference (I cite it as context only).

---

## Reproducibility commands

```bash
python scripts/ou_bao_likelihood.py
python scripts/profile_sigma_x_desi.py
python scripts/amplifier_audit.py
python scripts/amplification/run_all_routes.py
python scripts/slip_bridge.py
```

