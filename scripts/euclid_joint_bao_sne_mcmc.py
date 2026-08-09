#!/usr/bin/env python3
"""
DEPRECATED stub (was synthetic Euclid BAO+SN).

Real DESI DR2 BAO residual tests:

  python scripts/desi_dr2_real_bao_test.py
  python scripts/joint_w0wa_sigma_desi.py
  python scripts/profile_sigma_x_desi.py
"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

REAL = Path(__file__).resolve().parent / "desi_dr2_real_bao_test.py"


def main() -> int:
    print("euclid_joint_bao_sne_mcmc.py is retired (synthetic). Running DESI DR2 real BAO test:")
    print(" ", REAL)
    if not REAL.is_file():
        print("ERROR: missing", REAL)
        return 2
    sys.argv = [str(REAL)] + sys.argv[1:]
    runpy.run_path(str(REAL), run_name="__main__")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
