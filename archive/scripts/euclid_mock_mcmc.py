#!/usr/bin/env python3
"""
DEPRECATED stub (was synthetic Euclid BAO).

Real DESI DR2 BAO residual test (public alpha on disk):

  python scripts/desi_dr2_real_bao_test.py

Also: joint_w0wa_sigma_desi.py, profile_sigma_x_desi.py, ou_bao_likelihood.py
on the published DESI DR2 isotropic alpha table.
"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

REAL = Path(__file__).resolve().parent / "desi_dr2_real_bao_test.py"


def main() -> int:
    print("euclid_mock_mcmc.py is retired (synthetic). Running DESI DR2 real BAO test:")
    print(" ", REAL)
    if not REAL.is_file():
        print("ERROR: missing", REAL)
        return 2
    sys.argv = [str(REAL)] + sys.argv[1:]
    runpy.run_path(str(REAL), run_name="__main__")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
