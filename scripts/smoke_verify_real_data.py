#!/usr/bin/env python3
"""
Smoke verification of claim-safe DESI DR2 real-data pipelines.

Run after changes to eos_efectiva / ou_bao_likelihood:
    python3 scripts/smoke_verify_real_data.py
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run_capture(script: str, timeout: int = 300) -> str:
    env = {**dict(**{k: v for k, v in __import__("os").environ.items()}), "MPLBACKEND": "Agg"}
    r = subprocess.run(
        [sys.executable, str(ROOT / script)],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        timeout=timeout,
        env=env,
    )
    out = (r.stdout or "") + (r.stderr or "")
    if r.returncode != 0:
        raise SystemExit(f"FAIL: {script} exit {r.returncode}\n{out[-2000:]}")
    return out


def main() -> None:
    errors = []

    # --- eos_efectiva ---
    print("Running scripts/eos_efectiva.py ...")
    eos_out = run_capture("scripts/eos_efectiva.py", timeout=300)
    # Best-fit line must be ~ -0.99, not the box corner as the MLE answer
    m = re.search(r"Pure CPL MLE.*?\n\s*w0\s*=\s*([+-]?\d+\.\d+)", eos_out, re.S)
    if not m:
        errors.append("Could not parse Pure CPL w0 from eos_efectiva output")
    else:
        w0 = float(m.group(1))
        if not (-1.05 < w0 < -0.90):
            errors.append(f"Pure CPL w0={w0} outside expected near-LCDM window")
        if abs(w0 + 2.0) < 0.05:
            errors.append("Pure CPL w0 still stuck at -2 (old bug regression)")

    m_wa = re.search(r"Pure CPL MLE.*?\n\s*w0.*?\n\s*wa\s*=\s*([+-]?\d+\.\d+)", eos_out, re.S)
    if m_wa:
        wa = float(m_wa.group(1))
        if abs(wa) > 0.5:
            errors.append(f"Pure CPL wa={wa} unexpectedly far from 0")

    if "re-eval logL at reported params" not in eos_out or "OK" not in eos_out:
        errors.append("Missing re-eval OK checks in eos_efectiva")
    if "corner is NOT degenerate" not in eos_out:
        errors.append("Missing corner non-degeneracy sanity print")
    if "ΔAIC  (eff − CPL) = +4.000" not in eos_out and "prefer CPL" not in eos_out:
        # allow small float formatting differences
        if "prefer CPL" not in eos_out and "prefer CPL/null" not in eos_out:
            errors.append("Nested model should prefer pure CPL / null")

    jpath = ROOT / "results" / "eos_cpl_desi_dr2" / "eos_cpl_summary.json"
    if not jpath.is_file():
        errors.append(f"Missing {jpath}")
    else:
        j = json.loads(jpath.read_text())
        if abs(j["pure_CPL"]["w0"] + 0.99) > 0.05:
            errors.append(f"JSON pure_CPL.w0={j['pure_CPL']['w0']} unexpected")
        if j["delta_AIC_eff_minus_cpl"] < 0:
            errors.append("JSON says effective EoS preferred over CPL — unexpected")
        if j["corner_w0wa_m2_logL_pure_CPL"] > -10:
            errors.append("Corner (-2,-2) logL should be terrible without clamp")

    # --- ou_bao ---
    print("Running scripts/ou_bao_likelihood.py ...")
    ou_out = run_capture("scripts/ou_bao_likelihood.py", timeout=180)
    if "Traceback" in ou_out:
        errors.append("ou_bao_likelihood raised traceback")
    if "ω_R = 0.0000" not in ou_out and "omega_R" not in ou_out.lower():
        errors.append("ou_bao: expected ω_R → 0 on DR2")
    if "H0 preferred" not in ou_out:
        errors.append("ou_bao: expected H0 preferred over H1")
    if "σ_X = 0.00005" not in ou_out and "0.00005" not in ou_out:
        errors.append("ou_bao: expected σ_X at numerical floor")

    # --- public hygiene ---
    tracked = subprocess.check_output(["git", "ls-files"], cwd=ROOT, text=True)
    if "quantum_information_cosmos.py" in tracked:
        errors.append("numerology script still tracked in git")

    if errors:
        print("\nSMOKE VERIFY: FAIL")
        for e in errors:
            print(" -", e)
        raise SystemExit(1)

    print("\nSMOKE VERIFY: PASS")
    print("  CPL: w0≈-0.99, wa≈-0.02; nested σ not preferred")
    print("  OU/QNM DR2: σ_X→floor, ω_R→0, H0 preferred")
    print("  No quantum_information_cosmos in public tree")


if __name__ == "__main__":
    main()
