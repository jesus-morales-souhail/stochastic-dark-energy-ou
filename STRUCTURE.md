# Repository structure

```
manuscript/ Unified preprint + referee claim checklist (start here for publication)
papers/ Technical notes and data packs (English) — narrative arc in README.md
scripts/ Analysis code (BAO, amplification routes, desqueezing, GPE)
 desqueezing/ Open-system desqueezing and timescale mapping
 gpe/ Gross–Pitaevskii / Bogoliubov numerics
 amplification/ Routes 1–3 (seed / freeze / avalanche)
figures/ Figures
results/ Numerical outputs (BAO fits, forecasts, desqueezing, routes; euclid_mcmc* RETIRED synthetic)
notes/ Technical notes (desqueezing synthesis, mapping tables)
local_archive/ Offline drafts (not part of the public scientific claim)
```

**Reading guide:** `manuscript/PREPRINT.md` · **Non-claims fence:** `papers/EXPLORATORY_BOUNDARY.md`
**Exploratory repo (optics / method only):** https://github.com/jesus-morales-souhail/stochastic-de-exploratory-notes

This layout is how I keep the claim set readable.
---

## Narrative index (read with README.md)

### I — Empirical BAO / models

| Path | Role |
|------|------|
| `papers/resume.txt` | Compact numerical summary |
| `papers/covariance-desi-dr2-full.md` | Official 13×13 BAO cov → α 7×7 |
| `papers/stochastic-dark-energy-desi-dr2.md` | Main OU/QNM vs DESI DR2 |
| `papers/quantum-fluid-instabilities-desi-dr2.md` | Tachyonic fluid exclusion |
| `papers/sensitivity_kernel_table.md` | BAO sensitivity kernel $S(z)$ |
| `scripts/ou_bao_likelihood.py` | OU / QNM residual kernel |
| `scripts/eos_efectiva.py` | CPL background + nested $(\sigma,\theta)$ |
| `scripts/ou_bao_stochastic_test.py` | Alternate OU/QNM pipeline |
| `scripts/joint_w0wa_sigma_desi.py` | Joint $\{w_0,w_a,\sigma_X\}$ |
| `scripts/profile_sigma_x_desi.py` | Profile likelihood for $\sigma_X$ |
| `scripts/cross_correlation_DESI.py` | Galaxy × residual cross-check |
| `scripts/gpe/gpe_sim.py` | GPE illustration |
| `scripts/tachyonic_rank1_mle.py` | Rank-1 coherent growth MLE on real DESI $\alpha$ |
| `papers/quantum-fluid-instabilities-desi-dr2.md` | Tachyonic / GPE exclusion note |
| `results/eos_cpl_desi_dr2/` | CPL fit outputs |

### II — Geometry of vacuum smoothness

| Path | Role |
|------|------|
| `papers/principle-of-vacuum-smoothness.md` | Late-time vacuum homogeneity |
| `papers/smoothness-of-the-vacuum-unimodular.md` | Unimodular + smoothness |
| `papers/unimodular-gravity-vacuum-smoothness.md` | Extended geometric synthesis |
| `papers/sdiff-fundamental-vs-emergent.md` | Fundamental vs emergent SDiff |

### III — Desqueezing / amplification gap / Euclid

| Path | Role |
|------|------|
| `papers/amplification-gap.md` | **amplification closure:** $G_{\mathrm{Euclid}}\sim 10^{56}$, $G_{\mathrm{DESI}}\sim 10^{57}$; amplifier audit |
| `papers/amplification-no-free-lunch.md` | Redirect stub → `amplification-gap.md` (old bookmarks) |
| `papers/HONEST_HEADLINES.md` | Title discipline: do not claim vacuum energy was calculated |
| `scripts/gap_two_targets.py` | Exact two-gap arithmetic |
| `scripts/amplifier_audit.py` | Runnable gain table (Sorkin → DESI/Euclid) |
| `scripts/amplification/route1_local_causal_set_seed.py` | Route 1: N_eff / local seed scan |
| `scripts/amplification/route2_late_horizon_exit.py` | Route 2: θ(x) freeze-out Monte Carlo |
| `scripts/amplification/route3_nonlinear_avalanche.py` | Route 3: double-well avalanche scan |
| `scripts/amplification/run_all_routes.py` | Batch runner (`--heavy` for full CPU) |
| `results/amplification_routes/` | CSV outputs of routes 1–3 |
| `results/amplification_routes/VERDICT.md` | Machine-checked amplification table (N_eff, freeze gain, avalanche) |
| `notes/desqueezing-relaxation-vacuum-fluctuations-note.md` | Open-system half-life |
| `papers/fundamental-vs-emergent-vacuum-relaxation.md` | Path-integrated residuals; regions F/E0–E3 |
| `papers/euclid-protocol-vacuum-relaxation.md` | Euclid BAO protocol |
| `scripts/desqueezing/desqueezing_relax_time.py` | QuTiP Lindblad scan |
| `scripts/desqueezing/cosmological_mapping.py` | $\gamma \leftrightarrow \theta H(z)$ |
| `scripts/desi_dr2_real_bao_test.py` | DESI DR2 real alpha residual test |
| `scripts/joint_w0wa_sigma_desi.py` | Joint w0,wa,sigma_X on DESI BAO |

### IV — Option 0 (slip / anisotropic gap)

| Path | Role |
|------|------|
| `papers/anisotropic-slip-option0.md` | $\eta$ / shear vs $\sigma_X$; amplification inheritance |
| `papers/data-pack-option0-internet.md` | arXiv numbers (Maus, Sakr, DESI MG/BAO) |
| `scripts/slip_bridge.py` | Runnable $\sigma_X\to\|\gamma-1\|$ map; amplitude-starved verdict |

### V — Boundary (no fantasy copies in this repo)

| Path | Role |
|------|------|
| `papers/EXPLORATORY_BOUNDARY.md` | Fence: what is claim vs digression |
| `papers/HONEST_HEADLINES.md` | Allowed wording |
| `papers/THEORY_PROGRAMME_POINTER.md` | Points at theory repo |
| [stochastic-de-exploratory-notes](https://github.com/jesus-morales-souhail/stochastic-de-exploratory-notes) | **Canonical** home for optics / wrong-scale demos |

Optics, car–drone pupil, tesseract \(B_4\), wavefront analogies, Maxwell device, self-shielding essays, topological edge analogy: **only** in the exploratory repo. They are not under `papers/` here.

---

## Cohesion rule

Every claim in this repository should be classifiable as one of:

1. **BAO residual / model constraint**,
2. **Geometric interpretation** of smoothness,
3. **Amplification / Euclid protocol**,
4. **Slip Option 0 with amplitude honesty**, or
5. **Explicit boundary** (this section + CLAIMS.md).

Lab optics metaphors belong in **stochastic-de-exploratory-notes** — not in the DESI likelihood narrative.

See **`README.md`** for the full story and reading order.

