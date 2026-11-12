#!/usr/bin/env python3
"""Percentiles of the real Pantheon+ statistics within the pure-LCDM universes of INYECCION (scenario A), read from
inyeccion.json (per_universe.A and real_check); nothing is recalculated. Writes inyeccion_percentiles.json."""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
d = json.loads((HERE / "inyeccion.json").read_text(encoding="utf-8"))
A = {k: np.asarray(v, float) for k, v in d["per_universe"]["A"].items() if k in ("b5", "marg", "dlnL_max", "tension")}
rc = d["real_check"]
real = {"b5": rc["5pct_linear"], "marg": rc["marg_linear"], "dlnL_max": rc["dlnL_max_linear"], "tension": rc["tension_linear"]}
out = {"source": "inyeccion.json (per_universe.A, real_check)", "n_universes": int(len(A["b5"])), "real_linear": real,
       "percentile_of_real_in_A": {k: float(np.mean(A[k] <= real[k])) for k in real},
       "A_p95": {k: float(np.percentile(A[k], 95)) for k in real},
       "percentile_of_exact_5pct_0.169_in_A": float(np.mean(A["b5"] <= 0.169)), "percentile_of_exact_marg_0.287_in_A": float(np.mean(A["marg"] <= 0.287))}
(HERE / "inyeccion_percentiles.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(json.dumps(out, indent=1))
