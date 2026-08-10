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
| **C2** | $\sigma_X < 2.5\times 10^{-2}$ (95% CL, profile over $\theta$, 7 isotropic $\alpha$ bins, full meas. cov). Best fit $\sigma_X\to 0$ (null preferred). Former $1.5\times 10^{-4}$ was a working target **not** equal to this profile — corrected. | `profile_sigma_x_desi.py` / `results/profile_sigma_x/` |
| **C3** | Nested free background plus stochastic extension is not preferred ($\Delta\mathrm{AIC}\approx +4$) | `eos_efectiva.py` / joint fits |
| **C4** | Coherent tachyonic growth (rank-1 cov) is excluded: $\Delta\ln\mathcal{L}\approx -12.7$ | Quantum-fluid MLE |
| **C5** | Gap $\sigma_0\to 10^{-5}$ is $\sim 10^{56}$ (**Euclid target**); gap $\to 2.5\times 10^{-2}$ is $\sim 10^{59}$ (**DESI ceiling**); linear amplifiers fail both | `gap_two_targets.py` + `amplifier_audit.py` + PREPRINT §5 |
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
| **N9** | Synthetic Euclid MCMC folders (`results/euclid_mcmc/`, `results/euclid_joint_mcmc/`) are **retired**; not data claims |
| **N10** | I do not identify the residual ceiling \(\sigma_X<2.5\times 10^{-2}\) with the baryon density \(\omega_b\approx 0.022\); different quantities |
| **N11** | VLBI / EHT / photon sphere (\(1.5\,R_s\)) are **outside** the residual claim set; no map to \(\ell_\ast\) or \(\sigma_X\). Hygiene: [`papers/operator-scale-separation.md`](../papers/operator-scale-separation.md) |

---

## Dataset boundary

- **In:** public DESI DR2 BAO summary statistics + official Gaussian BAO 13×13 cov projected to α
  (`papers/covariance-desi-dr2-full.md`). Cross-bin α correlations are zero in that public product.
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

