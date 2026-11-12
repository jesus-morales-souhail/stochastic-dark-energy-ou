#!/usr/bin/env python3
"""
Downloads the supernova data files that are not redistributed in this repository (no explicit licence in their official
source: Pantheon+, DES-Dovekie, and the DES-SN5YR README) and verifies the sha256 of every file listed in datos/FUENTES.json,
including the ones kept in the repository (DES-SN5YR data, CC-BY-4.0 via Zenodo; Union3, MIT).

    python3 -I descargar_datos.py            # download what is missing or wrong, then verify everything
    python3 -I descargar_datos.py --verify   # only verify, download nothing

Standard library only. Downloaded files are untrusted input: they are only ever read as numbers by the analysis scripts.
"""

from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
FUENTES = HERE / "datos" / "FUENTES.json"
FOLDER = {"des_sn5yr": HERE / "datos" / "des_sn5yr", "dovekie": HERE / "datos" / "dovekie", "union3": HERE / "datos" / "union3",
          "pantheon_plus": HERE.parent / "sn" / "datos"}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "descargar_datos.py (stochastic-dark-energy-ou)"})
    with urllib.request.urlopen(req, timeout=600) as r, tmp.open("wb") as f:
        for chunk in iter(lambda: r.read(1 << 20), b""):
            f.write(chunk)
    tmp.replace(dest)


def main(verify_only: bool) -> int:
    fuentes = json.loads(FUENTES.read_text(encoding="utf-8"))
    bad = 0
    for samp, folder in FOLDER.items():
        for name, info in fuentes[samp]["ficheros"].items():
            dest = folder / name
            want = info["sha256"]
            status = "ok" if dest.exists() and sha256(dest) == want else ("missing" if not dest.exists() else "sha256 mismatch")
            if status != "ok" and not verify_only:
                print(f"{samp}/{name}: {status}, downloading {info['url']}", flush=True)
                folder.mkdir(parents=True, exist_ok=True)
                fetch(info["url"], dest)
                got = sha256(dest)
                if got != want:
                    dest.unlink()
                    print(f"{samp}/{name}: sha256 after download {got} != expected {want}; file removed", flush=True)
                    bad += 1
                    continue
                status = "ok (downloaded)"
            if status.startswith("ok"):
                print(f"{samp}/{name}: {status}")
            else:
                print(f"{samp}/{name}: {status}")
                bad += 1
    print("all files verified" if bad == 0 else f"{bad} file(s) missing or wrong")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main("--verify" in sys.argv[1:]))
