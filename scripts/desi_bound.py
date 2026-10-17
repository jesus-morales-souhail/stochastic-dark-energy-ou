#!/usr/bin/env python3
"""
DESI DR2 BAO 95% CL bound on sigma_X, read from the profile likelihood output.

Single source for every script that uses the bound as a constant:
results/profile_sigma_x/profile_sigma_x.json -> sigma_X_95CL_upper_profile
(root of Delta lnL = -1.92 of the profile over theta; written by scripts/profile_sigma_x_desi.py).
"""

from __future__ import annotations

import json
from pathlib import Path

PROFILE_JSON = Path(__file__).resolve().parents[1] / "results" / "profile_sigma_x" / "profile_sigma_x.json"


def desi_sigma_x_bound() -> float:
    """sigma_X_95CL_upper_profile (2.633e-2 since 2026-10-17)."""
    return float(json.loads(PROFILE_JSON.read_text(encoding="utf-8"))["sigma_X_95CL_upper_profile"])
