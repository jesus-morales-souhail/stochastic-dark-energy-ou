# Repository structure

```
manuscript/ Unified preprint + referee claim checklist (START HERE for publication)
papers/ Technical notes and data packs (English) — narrative arc in README.md
scripts/ Analysis code and small demos
 desqueezing/ Open-system desqueezing and timescale mapping
 gpe/ Gross–Pitaevskii / Bogoliubov numerics
 amplification/ Routes 1–3 (seed / freeze / avalanche)
figures/ Figures
results/ Numerical outputs (BAO fits, forecasts, desqueezing, routes)
notes/ Technical notes (desqueezing synthesis, mapping tables)
local_archive/ Offline drafts (not part of the public scientific claim)
```

**Reading guide:** `manuscript/PREPRINT.md` · **Non-claims fence:** `papers/EXPLORATORY_BOUNDARY.md` 
**Exploratory repo (separate):** `../stochastic-de-exploratory-notes/`

---

## Narrative index (read with README.md)

### I — Empirical BAO / models

| Path | Role |
|------|------|
| `papers/resume.txt` | Compact numerical summary |
| `papers/stochastic-dark-energy-desi-dr2.md` | Main OU/QNM vs DESI DR2 |
| `papers/quantum-fluid-instabilities-desi-dr2.md` | Tachyonic fluid exclusion |
| `papers/sensitivity_kernel_table.md` | BAO sensitivity kernel \(S(z)\) |
| `scripts/ou_bao_likelihood.py` | OU / QNM residual kernel |
| `scripts/eos_efectiva.py` | CPL background + nested \((\sigma,\theta)\) |
| `scripts/ou_bao_stochastic_test.py` | Alternate OU/QNM pipeline |
| `scripts/joint_w0wa_sigma_desi.py` | Joint \(\{w_0,w_a,\sigma_X\}\) |
| `scripts/profile_sigma_x_desi.py` | Profile likelihood for \(\sigma_X\) |
| `scripts/cross_correlation_DESI.py` | Galaxy × residual cross-check |
| `scripts/gpe/gpe_sim.py` | GPE illustration |
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
| `papers/amplification-gap.md` | **Act III closure:** gap \(10^{56}\), amplifier audit, only (d) worth new theory |
| `papers/amplification-no-free-lunch.md` | Redirect stub → `amplification-gap.md` (old bookmarks) |
| `scripts/amplifier_audit.py` | Runnable gain table (Sorkin → DESI/Euclid) |
| `scripts/amplification/route1_local_causal_set_seed.py` | Route 1: N_eff / local seed scan |
| `scripts/amplification/route2_late_horizon_exit.py` | Route 2: θ(x) freeze-out Monte Carlo |
| `scripts/amplification/route3_nonlinear_avalanche.py` | Route 3: double-well avalanche scan |
| `scripts/amplification/run_all_routes.py` | Batch runner (`--heavy` for full CPU) |
| `results/amplification_routes/` | CSV outputs of routes 1–3 |
| `results/amplification_routes/VERDICT.md` | Machine-checked Act III table (N_eff, freeze gain, avalanche) |
| `notes/desqueezing-relaxation-vacuum-fluctuations-note.md` | Open-system half-life |
| `papers/fundamental-vs-emergent-vacuum-relaxation.md` | Path-integrated residuals; regions F/E0–E3 |
| `papers/euclid-protocol-vacuum-relaxation.md` | Euclid BAO protocol |
| `scripts/desqueezing/desqueezing_relax_time.py` | QuTiP Lindblad scan |
| `scripts/desqueezing/cosmological_mapping.py` | \(\gamma \leftrightarrow \theta H(z)\) |
| `scripts/euclid_mock_mcmc.py` | Euclid mock MCMC |
| `scripts/euclid_joint_bao_sne_mcmc.py` | Joint BAO+SN mock MCMC |

### IV — Option 0 (slip / anisotropic gap)

| Path | Role |
|------|------|
| `papers/anisotropic-slip-option0.md` | \(\eta\) / shear vs \(\sigma_X\); amplification inheritance |
| `papers/data-pack-option0-internet.md` | arXiv numbers (Maus, Sakr, DESI MG/BAO) |
| `scripts/slip_bridge.py` | Runnable \(\sigma_X\to\|\gamma-1\|\) map; amplitude-starved verdict |

### V — Method hygiene and closed wrong paths

| Path | Role |
|------|------|
| `papers/scale-operator-experiment-map.md` | Master map: closed vs next (scale + operator) |
| `papers/self-shielding-vs-untestability.md` | Shielding vs amplitude honesty |
| `papers/self-shielding-triple-barrier.md` | Born / band-limit / OU: 0% physical abuse, 100% math |
| `papers/no-go-superoscillation-tesseract.md` | No-go sub-Airy + B4/Sp(4,R) seal |
| `papers/maxwell-device-to-M-derivation.md` | Maxwell+device→M; only Sp(4,R) is legitimate 4D |
| `papers/optics-ou-analogies-and-limits.md` | Thin lens OK; tesseract ABCD no; \(f\leftrightarrow 1/\theta\) analogy only |
| `papers/wavefront-shaping-vs-ou-vacuum.md` | Lab \(T\) ≠ BAO OU |
| `papers/pattern-undeclared-physical-power.md` | \(\ln 4\), \(B_4\) pattern filter |
| `papers/car-drone-pupil-newton-einstein.md` | Highway pupil: **CLOSED** wrong scale/operator |
| `scripts/superoscillation_energy_cost_demo.py` | 1D energy-tax demo |
| `scripts/b4_symplectic_count.py` | Seal: 32/384 B4 optical-symplectic; commutant dim 1 |
| `scripts/derive_M_maxwell_device.py` | Maxwell→Fresnel→ABCD + B4 restriction check |
| `scripts/car_drone_pupil_newton_einstein.py` | Runnable Newton vs SR numbers |

### Separated (not the peer-review claim set)

| Path | Content |
|------|---------|
| `../stochastic-de-exploratory-notes/` | Optics / scale-operator pedagogy and no-gos |
| `~/Proyectos/04_Optica_medios_complejos/` | Lab transmission-matrix wavefront shaping experiment |
| `papers/EXPLORATORY_BOUNDARY.md` | Fence between claims and digressions |

---

## Cohesion rule

Every claim in this repository should be classifiable as one of:

1. **BAO residual / model constraint** (Act I), 
2. **Geometric interpretation** of smoothness (Act II), 
3. **Amplification / Euclid protocol** (Act III), 
4. **Slip Option 0 with amplitude honesty** (Act IV), or 
5. **Explicit no-go / boundary** (Act V).

Anything that is only a lab optics metaphor belongs in Act V or in the **separated** optics project — not in the DESI likelihood narrative.

See **`README.md`** for the full story and reading order.
