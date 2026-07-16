#!/usr/bin/env python3
"""
Regression checks for the DESI DR2 BAO MLE scripts.

  python3 scripts/verify_bao_pipelines.py

Confirms that ``eos_efectiva.py`` returns a near-ΛCDM CPL background and that
``ou_bao_likelihood.py`` drives the residual stochastic amplitude to the floor.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run_script(rel: str, timeout: int = 300) -> str:
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
    env["MPLBACKEND"] = "Agg"
    r = subprocess.run(
        [sys.executable, str(ROOT / rel)],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        timeout=timeout,
        env=env,
    )
    out = (r.stdout or "") + (r.stderr or "")
    if r.returncode != 0:
        raise SystemExit(f"{rel} failed (exit {r.returncode})\n{out[-2000:]}")
    return out


def main() -> None:
    errors: list[str] = []

    print("eos_efectiva.py ...")
    eos = run_script("scripts/eos_efectiva.py", timeout=300)
    m = re.search(r"Pure CPL MLE.*?\n\s*w0\s*=\s*([+-]?\d+\.\d+)", eos, re.S)
    if not m:
        errors.append("Could not parse pure CPL w0")
    else:
        w0 = float(m.group(1))
        if not (-1.05 < w0 < -0.90):
            errors.append(f"Pure CPL w0 = {w0} (expected near −1)")

    jpath = ROOT / "results" / "eos_cpl_desi_dr2" / "eos_cpl_summary.json"
    if not jpath.is_file():
        errors.append(f"Missing {jpath}")
    else:
        j = json.loads(jpath.read_text())
        w0j = j.get("pure_CPL", {}).get("w0")
        if w0j is None or abs(w0j + 0.99) > 0.05:
            errors.append(f"JSON pure_CPL.w0 unexpected: {w0j}")
        if j.get("delta_AIC_extension_minus_cpl", -1) < 0:
            errors.append("Nested extension should not beat pure CPL on AIC")

    print("ou_bao_likelihood.py ...")
    ou = run_script("scripts/ou_bao_likelihood.py", timeout=180)
    if "Traceback" in ou:
        errors.append("ou_bao_likelihood raised an exception")
    if "ω_R = 0.0000" not in ou and "ω_R = 0" not in ou:
        errors.append("Expected ω_R → 0 for free QNM on DR2")
    if "H0 preferred" not in ou:
        errors.append("Expected H0 preferred over H1")

    if errors:
        print("VERIFY: FAIL")
        for e in errors:
            print(" -", e)
        raise SystemExit(1)

    print("VERIFY: OK")
    print("  CPL background near ΛCDM; nested σ not preferred")
    print("  OU/QNM residual amplitude at floor; ω_R → 0")


if __name__ == "__main__":
    main()
