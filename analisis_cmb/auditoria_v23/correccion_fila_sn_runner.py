#!/usr/bin/env python3
"""
Runner for CORRECCION_FILA_SN: re-runs the A-dependent linear parts of positividad_lin, familias, muestras_sn and inyeccion with the OLD
build_A (reproduction check) or with the corrected respuesta_v3A (mode 'new'), writing into an output directory. The target module is
imported FIRST so that its own thread settings apply as in the original runs; nothing in the old scripts is modified (monkeypatching of
the A builders and of the output directory only).

    ../venv/bin/python -I correccion_fila_sn_runner.py <positividad_lin|familias|muestras_sn|inyeccion> <old|new> <outdir>
"""
import sys, importlib, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
target, mode, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
mod = importlib.import_module(target)  # first import: the module's os.environ.setdefault for BLAS threads applies as originally
import json, os, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import marginal_guiada as MG
import muestras_sn as MS
import respuesta_v3A as R
out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
if mode == "new":
    MG.build_A = R.build_A_v3
    MS.make_linear = R.make_linear_v3
t0 = time.time()
if target == "positividad_lin":
    mod.HERE = out
    mod.main()
elif target == "familias":
    F = mod
    F.init_linear()
    pairs = [(fam, l) for fam in F.FAMILIES for l in F.ELLS]
    with Pool(min(16, len(pairs))) as pool:
        lin = pool.map(F._linear, pairs, chunksize=1)
    (out / "familias_linear.json").write_text(json.dumps({f"{r['fam']}|l_{r['l']:g}": r for r in lin}, indent=1, default=F._json_default), encoding="utf-8")
elif target == "muestras_sn":
    samples = MS.load_samples()
    for samp, sn in samples.items():
        MS.LIN[samp] = MS.make_linear(sn)   # patched to make_linear_v3 in mode 'new'
    ljobs = [(samp, fam, l) for samp in MS.SAMPLES for fam, l in MS.CELLS]
    with Pool(16) as pool:
        lin = pool.map(MS._linear, ljobs, chunksize=1)
    (out / "muestras_sn_linear.json").write_text(json.dumps({f"{r['sample']}|{r['cell']}": r for r in lin}, indent=1, default=float), encoding="utf-8")
elif target == "inyeccion":
    mod.HERE = out
    mod.main(False)
print(f"{target} {mode}: {time.time() - t0:.0f} s, threads OMP={os.environ.get('OMP_NUM_THREADS')}", flush=True)
