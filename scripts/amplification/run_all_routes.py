#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run amplifier audit + routes 1–3 (heavy optional)."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run(cmd: list[str]) -> None:
    print("\n>>>", " ".join(cmd), flush=True)
    r = subprocess.run(cmd, cwd=str(ROOT))
    if r.returncode != 0:
        sys.exit(r.returncode)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--heavy", action="store_true")
    args = ap.parse_args()
    py = sys.executable

    run([py, "scripts/amplifier_audit.py"])
    run([py, "scripts/amplification/route1_local_causal_set_seed.py"])
    cmd2 = [py, "scripts/amplification/route2_late_horizon_exit.py"]
    cmd3 = [py, "scripts/amplification/route3_nonlinear_avalanche.py"]
    if args.heavy:
        cmd2.append("--heavy")
        cmd3.append("--heavy")
    run(cmd2)
    run(cmd3)
    print("\nAll amplification route scans finished.")
    print("Results under results/amplification_routes/")


if __name__ == "__main__":
    main()
