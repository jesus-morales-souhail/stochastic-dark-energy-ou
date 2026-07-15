"""
Cosmological mapping FROM FIRST PRINCIPLES using ONLY:
  - canonical numbers from stochastic-dark-energy-ou (GitHub + local paper)
  - standard LCDM integrals with the paper fiducial (H0, Om)
  - QuTiP law measured in this project: t_1/2(|<a^2>|) = ln2 / gamma

NO illustrative A0=0.01. Seed amplitude from repo Axiom A2 / quantum-fluid paper:
  sigma_0 ~ 1/sqrt(N) ~ 10^{-61}  (N ~ 10^{122} Bekenstein-Hawking)
"""

from __future__ import annotations

import os
import csv
import numpy as np
from scipy.integrate import quad

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results", "cosmological_mapping_from_repo")
DOCS = os.path.join(ROOT, "docs")
os.makedirs(RESULTS, exist_ok=True)

# =============================================================================
# CANONICAL INPUTS — stochastic-dark-energy-ou (resume.txt, papers, scripts)
# =============================================================================

# Fiducial cosmology (sensitivity_kernel_table.md, ou_bao_stochastic_test.py)
H0_KMS = 67.4          # km/s/Mpc
OM = 0.315             # Omega_m
OL = 1.0 - OM          # flat LCDM

# DESI BAO path (sensitivity_kernel_table.md)
Z_DESI_MIN, Z_DESI_MAX = 0.295, 2.330
DX_DESI = 0.94          # ln(1+zmax)-ln(1+zmin) as in repo table
Z_EFF = np.array([0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.330])
S_Z = np.array([-0.284, -0.462, -0.595, -0.719, -0.870, -0.917, -1.070])

# DESI DR2 MLE (resume.txt, stochastic-dark-energy-desi-dr2.md)
SIGMA_X_LIMIT_95 = 1.5e-4
SIGMA_X_MLE_FLOOR = 5.0e-5
THETA_MLE_FLOOR = 1.0e-3
THETA_MIN_OPT = 1.0e-3

# Superseded calibration still present in script (for comparison ONLY)
THETA_BASE_OLD = 1.2
SIGMA_OU_OLD = 2.31e-2

# Euclid targets (sdiff-fundamental-vs-emergent.md, main paper)
SIGMA_X_EUCLID_TARGET = 1.0e-5
SIGMA_X_EUCLID_DETECT_WINDOW = (0.0, 1.0e-4)  # 0 < sigma_X ≲ 10^{-4}

# Axiom A1/A2 (stochastic-dark-energy-desi-dr2.md, quantum-fluid paper)
N_BH = 1.0e122                    # Bekenstein-Hawking DOF scale
SIGMA_0_SORKIN = 1.0 / np.sqrt(N_BH)  # ~ 1e-61 ; also sigma_0 ~ 10^{-61} in GPE paper
# quantum-fluid-instabilities-desi-dr2.md: "sigma_0 ~ 10^{-61}"
SIGMA_0_REPO = 1.0e-61

# Recombination (standard cosmology; z_rec ~ 1090)
Z_RECOMB = 1090.0

# QuTiP result (this project, desqueezing_relax_time.py)
# t_1/2 = ln2 / gamma  (gamma in inverse lab time)
LN2 = np.log(2.0)

# Conversion: H0^{-1} in Gyr
# Standard: t_H [Gyr] = 977.8 / (H0 in km/s/Mpc)
T_H0_GYR = 977.8 / H0_KMS
H0_PER_GYR = 1.0 / T_H0_GYR


def E_z(z):
    return np.sqrt(OM * (1 + z) ** 3 + OL)


def H_z_per_Gyr(z):
    return H0_PER_GYR * E_z(z)


def age_to_z(z):
    """Cosmic time from BB to redshift z (Gyr), flat LCDM."""
    # t = int_z^inf dz' / [(1+z') H(z')]
    def integrand(zp):
        return 1.0 / ((1.0 + zp) * H_z_per_Gyr(zp))

    val, _ = quad(integrand, z, 1.0e5, epsabs=1e-8)
    return val


def lookback_to_z(z):
    """Lookback time from z=0 to redshift z (Gyr)."""
    def integrand(zp):
        return 1.0 / ((1.0 + zp) * H_z_per_Gyr(zp))

    val, _ = quad(integrand, 0.0, z, epsabs=1e-8)
    return val


def delta_x(z1, z2):
    """Delta x = ln a2 - ln a1 = ln((1+z1)/(1+z2)) with z1>z2 (past to present)."""
    return np.log((1.0 + z1) / (1.0 + z2))


def main():
    lines = []

    def p(s=""):
        lines.append(s)
        print(s)

    p("=" * 70)
    p("FIRST-PRINCIPLES COSMOLOGICAL MAPPING (repo + LCDM fiducial)")
    p("=" * 70)
    p("Sources: stochastic-dark-energy-ou (GitHub main + local paper_OU)")
    p("         QuTiP law from this project: t_1/2 = ln2/gamma")
    p("")

    # ---- LCDM ages from paper fiducial ----
    t0 = age_to_z(0.0)
    t_rec = age_to_z(Z_RECOMB)
    t_since_rec = t0 - t_rec
    dx_recomb = delta_x(Z_RECOMB, 0.0)
    dx_desi_check = delta_x(Z_DESI_MAX, Z_DESI_MIN)

    p("--- Fiducial cosmology (repo) ---")
    p(f"H0 = {H0_KMS} km/s/Mpc, Omega_m = {OM}, Omega_L = {OL}")
    p(f"t_H0 = H0^{{-1}} = {T_H0_GYR:.4f} Gyr")
    p(f"Age of Universe t0 = {t0:.4f} Gyr  (integral LCDM)")
    p(f"Age at recombination t(z={Z_RECOMB}) = {t_rec*1e3:.3f} Myr")
    p(f"Time since recombination = {t_since_rec:.4f} Gyr")
    p(f"Delta x since recomb = ln(1+z_rec) = {dx_recomb:.4f}")
    p(f"Delta x DESI (repo table) = {DX_DESI:.2f}")
    p(f"Delta x DESI (from z_eff ends) = {dx_desi_check:.4f}")
    p("")

    # ---- Repo observational anchors ----
    p("--- Repo observational anchors ---")
    p(f"sigma_X < {SIGMA_X_LIMIT_95:.1e} (95% CL, DESI DR2)")
    p(f"MLE floor: theta ~ {THETA_MLE_FLOOR}, sigma_X ~ {SIGMA_X_MLE_FLOOR:.1e}")
    p(f"Old (superseded) calibration: theta={THETA_BASE_OLD}, sigma_X={SIGMA_OU_OLD}")
    p(f"N_BH ~ {N_BH:.0e}, sigma_0 = 1/sqrt(N) = {SIGMA_0_SORKIN:.3e}")
    p(f"sigma_0 (repo GPE text) = {SIGMA_0_REPO:.1e}")
    p(f"Euclid target ~ {SIGMA_X_EUCLID_TARGET:.0e}; detect window 0 < sigma_X ≲ 1e-4")
    p("")

    # ---- Mapping dictionary FROM REPO EQUATIONS ----
    # Paper: dX/dt + ... = noise; in x=ln a: dX = -theta X dx + sigma dW
    # => physical damping rate Gamma_phys = theta * H(t)
    # Identify lab gamma with Gamma_phys at z=0: gamma = theta * H0
    p("--- Mapping (derived from repo OU, not ad hoc) ---")
    p("From paper: dx = H dt  =>  dX/dt = -theta H X + ...")
    p("=> Gamma_phys(t) = theta * H(t)")
    p("Mapping A (z=0): gamma = theta * H0")
    p("Mapping B (OU-native): |<X(x)X(0)>| ~ exp(-theta |x|)  =>  Dx_1/2 = ln2/theta")
    p("QuTiP: t_1/2(|a^2|) = ln2/gamma  =>  consistent if gamma <-> Gamma_phys")
    p("")

    # ---- Table 1: theta values from REPO and their physical half-lives ----
    p("--- Table 1: theta FROM REPO -> t_1/2 and gamma ---")
    theta_cases = [
        ("MLE floor (resume/paper)", THETA_MLE_FLOOR),
        ("Old base THETA_BASE=1.2 (superseded)", THETA_BASE_OLD),
        ("Half-decay over DESI path Dx=0.94", LN2 / DX_DESI),
        ("Half-decay since recombination", LN2 / dx_recomb),
        ("Hubble-scale: t_1/2 = t_H0 (Map A)", LN2),  # theta = ln2 when t_1/2=t_H0
        ("Age-scale: t_1/2 = t0 (Map A)", LN2 * T_H0_GYR / t0),
    ]

    table1 = []
    for name, th in theta_cases:
        gamma = th * H0_PER_GYR  # Mapping A
        t_half_A = LN2 / gamma if gamma > 0 else np.nan  # Gyr
        dx_half = LN2 / th if th > 0 else np.nan
        # residual of Sorkin seed over DESI and recomb paths
        res_desi = SIGMA_0_REPO * np.exp(-th * DX_DESI)
        res_rec = SIGMA_0_REPO * np.exp(-th * dx_recomb)
        # stationary OU variance if sigma chosen so Var=sigma^2/(2theta) equals seed^2 scale
        # (order of magnitude only)
        row = {
            "case": name,
            "theta": th,
            "gamma_per_Gyr": gamma,
            "t_half_MapA_Gyr": t_half_A,
            "Delta_x_half": dx_half,
            "sigma_res_Sorkin_DESI": res_desi,
            "sigma_res_Sorkin_recomb": res_rec,
        }
        table1.append(row)
        p(
            f"{name}\n"
            f"  theta={th:.6e}  gamma={gamma:.6e}/Gyr  t_1/2(A)={t_half_A:.4e} Gyr\n"
            f"  Dx_1/2={dx_half:.4e}  "
            f"sigma_res(Sorkin,DESI)={res_desi:.3e}  "
            f"sigma_res(Sorkin,recomb)={res_rec:.3e}"
        )

    # ---- Table 2: What A0 is needed to hit DESI/Euclid after path damping ----
    p("")
    p("--- Table 2: Required kick A0 so residual = target after path ---")
    p("sigma_res = A0 * exp(-theta * Dx)  =>  A0 = target * exp(+theta * Dx)")
    p("Compare to Sorkin A0 = 1e-61")
    table2 = []
    for th_name, th in [
        ("theta_MLE=0.001", THETA_MLE_FLOOR),
        ("theta=1.2 (old)", THETA_BASE_OLD),
        ("theta=ln2~0.693 (Hubble half-life)", LN2),
        ("theta=0.099 (half over recomb path)", LN2 / dx_recomb),
    ]:
        for path_name, dx in [("DESI_0.94", DX_DESI), ("recomb", dx_recomb)]:
            for target_name, target in [
                ("DESI_limit_1.5e-4", SIGMA_X_LIMIT_95),
                ("Euclid_1e-5", SIGMA_X_EUCLID_TARGET),
            ]:
                A0_need = target * np.exp(th * dx)
                ratio = A0_need / SIGMA_0_REPO
                row = {
                    "theta_case": th_name,
                    "theta": th,
                    "path": path_name,
                    "Delta_x": dx,
                    "target": target_name,
                    "target_value": target,
                    "A0_required": A0_need,
                    "A0_over_Sorkin": ratio,
                }
                table2.append(row)
                p(
                    f"  {th_name:35s} {path_name:10s} -> {target_name:20s}: "
                    f"A0_req={A0_need:.3e}  (= {ratio:.3e} x Sorkin)"
                )

    # ---- Table 3: Page / scrambling with S from repo ----
    p("")
    p("--- Table 3: Horizon scales with S=N_BH=1e122 (repo Axiom A1) ---")
    t_scr = T_H0_GYR * np.log(N_BH)  # dS scrambling ~ t_H ln S
    t_page_sqrt = T_H0_GYR * np.sqrt(N_BH)
    t_page_S = T_H0_GYR * N_BH
    table3 = []
    for name, T in [
        ("t_H0", T_H0_GYR),
        ("t0_LCDM", t0),
        ("t_since_recomb", t_since_rec),
        ("t_scr = t_H ln S", t_scr),
        ("t_Page ~ t_H sqrt(S)", t_page_sqrt),
        ("t_Page ~ t_H S", t_page_S),
    ]:
        gamma = LN2 / T
        theta_A = gamma / H0_PER_GYR
        table3.append(
            {
                "scale": name,
                "T_Gyr": T,
                "gamma_per_Gyr": gamma,
                "theta_MapA": theta_A,
                "t_half_equals_T": True,
            }
        )
        p(f"  {name}: T={T:.6e} Gyr  => gamma={gamma:.6e}/Gyr  theta_A={theta_A:.6e}")

    # ---- Table 4: DESI BAO bins from repo (for documentation) ----
    table4 = []
    for z, s in zip(Z_EFF, S_Z):
        table4.append({"z_eff": z, "S_z": s, "x_ln_a": -np.log(1 + z)})

    # ---- Critical conclusions ----
    p("")
    p("--- CRITICAL FIRST-PRINCIPLES CONCLUSIONS ---")
    p("1) Pure Sorkin seed A0~1e-61: after ANY cosmological path,")
    p("   residual sigma_X ~ 1e-61 << DESI (1.5e-4) and Euclid (1e-5).")
    p("   => Euclid CANNOT see pure Sorkin Poisson noise without amplification.")
    p("2) DESI null (sigma_X < 1.5e-4) does NOT constrain Sorkin seed;")
    p("   it constrains effective late-time amplitude after any growth/amplification.")
    p("3) For t_1/2 ~ age of Universe: theta ~ O(1), gamma ~ H0 (Map A).")
    p("4) MLE theta~0.001: essentially NO damping on DESI path (Dx=0.94);")
    p("   damping factor exp(-0.001*0.94)~0.999.")
    p("5) Old theta=1.2: strong damping; half-life Dx=ln2/1.2~0.58 e-folds.")
    p("6) Emergent SDiff with finite tau is viable only if there is a mechanism")
    p("   that amplifies seed to A0 ≳ 1e-5..1e-4 before (or instead of) pure Sorkin.")
    p("7) Fundamental SDiff: sigma_X=0 exactly — residual identically zero.")
    p("")

    # Stationary OU: if sigma diffusion from Sorkin
    # Var_st = sigma^2/(2theta). If we identify rms with sigma_0:
    p("--- Stationary OU variance (repo formula Var=sigma^2/(2theta)) ---")
    for th in [THETA_MLE_FLOOR, THETA_BASE_OLD, LN2]:
        # If diffusion sigma_diff set so that sqrt(Var)=sigma_0 when balanced:
        # need sigma_diff = sigma_0 * sqrt(2*theta)
        sig_diff = SIGMA_0_REPO * np.sqrt(2 * th)
        var = sig_diff**2 / (2 * th)
        p(
            f"  theta={th:.4f}: if rms_st=sigma_0={SIGMA_0_REPO:.1e}, "
            f"then sigma_diff={sig_diff:.3e}, Var={var:.3e}"
        )

    # Write outputs
    def write_csv(path, rows):
        if not rows:
            return
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    write_csv(os.path.join(RESULTS, "table1_theta_from_repo.csv"), table1)
    write_csv(os.path.join(RESULTS, "table2_A0_required_vs_Sorkin.csv"), table2)
    write_csv(os.path.join(RESULTS, "table3_page_scrambling_repo_S.csv"), table3)
    write_csv(os.path.join(RESULTS, "table4_DESI_S_z_repo.csv"), table4)

    summary_path = os.path.join(RESULTS, "FROM_REPO_mapping_resumen.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    # Update mapping_tables.md with FIRST PRINCIPLES version
    md_path = os.path.join(DOCS, "mapping_tables_from_repo.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("# First-Principles Mapping Tables (repo anchors only)\n\n")
        f.write("**No illustrative \(A_0=0.01\).** Seed: \(\\sigma_0\\sim 10^{-61}\) (repo Axiom A2 / GPE paper).\n\n")
        f.write(f"Fiducial: \(H_0={H0_KMS}\), \(\\Omega_m={OM}\), \(t_0={t0:.3f}\\,\\mathrm{{Gyr}}\).\n\n")
        f.write("## Table 1 — \(\\theta\) from repo → \(\\gamma\), \(t_{1/2}\), Sorkin residual\n\n")
        f.write("| Case | \(\\theta\) | \(\\gamma\) [1/Gyr] | \(t_{1/2}\) [Gyr] | \(\\sigma_{\\rm res}\) (DESI path) |\n")
        f.write("|------|----------:|-------------------:|------------------:|----------------------------------:|\n")
        for r in table1:
            f.write(
                f"| {r['case']} | {r['theta']:.4e} | {r['gamma_per_Gyr']:.4e} | "
                f"{r['t_half_MapA_Gyr']:.4e} | {r['sigma_res_Sorkin_DESI']:.3e} |\n"
            )
        f.write("\n## Table 2 — \(A_0\) required to hit DESI/Euclid (vs Sorkin \(10^{-61}\))\n\n")
        f.write("| \(\\theta\) case | path | target | \(A_0\) required | \(A_0/\\sigma_0\) |\n")
        f.write("|----------------|------|--------|-----------------:|-----------------:|\n")
        for r in table2:
            f.write(
                f"| {r['theta_case']} | {r['path']} | {r['target']} | "
                f"{r['A0_required']:.3e} | {r['A0_over_Sorkin']:.3e} |\n"
            )
        f.write("\n## Table 3 — Page/scrambling with \(S=10^{122}\)\n\n")
        f.write("| Scale | T [Gyr] | \(\\gamma\) | \(\\theta_A\) |\n")
        f.write("|-------|--------:|-----------:|-------------:|\n")
        for r in table3:
            f.write(
                f"| {r['scale']} | {r['T_Gyr']:.6e} | {r['gamma_per_Gyr']:.6e} | {r['theta_MapA']:.6e} |\n"
            )

    p(f"\nWrote {RESULTS}")
    p(f"Wrote {md_path}")
    p(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
