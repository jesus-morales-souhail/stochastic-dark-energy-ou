# analisis_cmb — OU at the DESI maximum variance vs Planck 2018

Question (current version): with σ_X at the 95 % limit of the DESI DR2 BAO profile
(`../results/profile_sigma_x/profile_sigma_x.json`; 2.633e-2 since 2026-10-17, previously the grid
node 2.506e-2), does any CMB observable (TT at l<30, C_L^φφ, θ_*) separate the 68 % band of 200 OU
realisations from the w0-wa curve by more than Planck 2018 allows? What bound on σ_X would Planck
impose? Current method: the H0 of each curve is solved so that θ_* matches ΛCDM; observables TT l<30,
C_L^φφ and the H0 distribution.

It does not touch the repository's BAO fitting scripts. All decisions, provenance of values and open
problems: `NOTAS_ANALISIS_CMB.md`. Result: `RESULTADO.md`.

**Scope.** This analysis is not part of the preprint's claim set (C1–C7 in
[`../manuscript/CLAIMS.md`](../manuscript/CLAIMS.md)); the preprint does not cite it. The working
documents in this folder (`NOTAS_ANALISIS_CMB.md`, `RESULTADO.md`, the audit and the block reports)
are in Spanish.

```
venv/                 CAMB 2.0.4 (ignored by git)
datos_planck/         files from the Planck Legacy Archive (ignored by git)
ou_cmb.py             X(ln a) → ρ_DE → w(a) → CAMB (PPF)
planck_datos.py       loads the Planck errors (TT, lensing bandpowers, θ_*)
realizacion_unica.py  one realisation at θ=1 + convergence test with 500/2000 points
diagnosticos/         path-grid / table convergence tests, band per Δ
realizaciones.py      200 realisations per (θ, Δ, σ_X), H0 solved → resultados/<halofit>/sigma_<σ>/
estadisticos.py       Δχ², percentiles, H0 → resumen.json, tabla.md
resultados/           outputs (resultados/sigma_2.506e-02/: old method, fixed H0 — obsolete)
```

The old fixed-H0 scripts (formerly `diagnosticos/obsoleto/`) are in
[`../archive/analisis_cmb/diagnosticos/obsoleto/`](../archive/README.md).

```bash
./venv/bin/python realizacion_unica.py
./venv/bin/python realizaciones.py --halofit takahashi            # σ_X = 2.506e-2
./venv/bin/python estadisticos.py resultados/takahashi/sigma_2.506e-02
```

Status (2026-10-03): analysis finished; result in `RESULTADO.md` (re-run at σ_X = 2.633e-2 on
2026-10-17, `resultados/takahashi/sigma_2.633e-02/`).

## Other contents

| Path | Content |
|:--|:--|
| `everpresent/` | Everpresent Λ (Ahmed–Dodelson–Greene–Sorkin and follow-ups) against DESI DR2 BAO; `EVERPRESENT_vs_DESI.md`. The CMB enters through a recombination proxy and through θ* computed for each realization, with the expansion rate at recombination (§3f, `everpresent/cmb_real/theta_star.py`); the full Planck likelihood with a sign-changing ρ_DE was not run (`everpresent/cmb_real/EXPLORACION.md`) |
| `everpresent/cmb_real/valvula.py`, `valvula.json`, `VALVULA.md` | Role of the tanh map of Model 2 (Zwane et al. eq. 2.10): fraction of realizations where it acts, survivors and criterion (c) without it (our modification, not a published model) |
| `reloj/` | σ_X profile with four stationarity clocks (ln a, cosmic time, conformal time, z): `perfil_reloj.py`, `perfil_reloj.{json,txt}`, `RELOJ.md`; before/after of the related text in `COMPARACION_RELOJ_VALVULA.md` |
| `fondo/` | The DESI σ_X profile with the background free (`FONDO.md`, `perfil_fondo.py`): the α reference is the DESI DR2 best fit, the θ → 0 limit is the Ω_m direction (ln det term), θ ≥ 1 limits by background treatment; σ_X needed to resolve θ with DESI and a Euclid-like survey (`RESOLVER_THETA.md`, `resolver_theta.py`) |
| `fondo/TRECE.md` | The same profile with the 13 DESI numbers (D_M, D_H separately, D_V for BGS) instead of the 7 isotropic α (`perfil_13.py`): limits tighter by a factor 1.11–1.25, same meaning |
| `sn/` | σ_X profile with Pantheon+ (`SN.md`, `perfil_sn.py`, data in `sn/datos/` with sha256); effect of the limit on D_L for standard sirens |
| `cmb/` | High-ℓ CMB check (`PICOS.md`, `picos.py`): TT, TE, EE distance to ΛCDM over ℓ ∈ [30, 2500], lensing in the troughs, ∝ σ_X² |
| `oscilador/` | Deterministic oscillating X(ln a) against DESI DR2 BAO and Planck 2018; `oscilador/README.md` |
| `AUDITORIA.md` | Audit of the numbers across the repository (2026-10-03) |
| `bloque_A/`, `bloque_B/` | Correction blocks A, B1, B1b, B2: before/after comparisons, snapshots (`antes/`) and logs |
| `verificacion_numeros.md`, `verificacion_numeros.py` | Re-computation of the flagged numbers |
| `REORGANIZACION.md`, `REORGANIZACION_CHECK.md` | Repository reorganisation of 2026-10-17 and its checks |

The audit and the block reports quote paths as they were at commit `d71c4d2`, before the
reorganisation.

## Downloading the data again (`datos_planck/` is in .gitignore)

```bash
mkdir -p datos_planck
curl -L -o datos_planck/COM_PowerSpect_CMB-TT-full_R3.01.txt \
  "http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_PowerSpect_CMB-TT-full_R3.01.txt"
# sha256: ccf3113604020536f6f13ccf51680a7316ad0f32da558eee7f625e613bdd5522
./venv/bin/pip install camb cobaya numpy scipy matplotlib
./venv/bin/cobaya-install planck_2018_lensing.native -p "$PWD/datos_planck/cobaya_packages"
```
