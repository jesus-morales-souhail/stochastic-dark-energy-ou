#!/usr/bin/env python3
"""Check that every relative Markdown link [..](path) in the repo points to an existing file.

Usage: python analisis_cmb/check_links.py  -> prints a Markdown report (used for REORGANIZACION_CHECK.md).
External links (http, https, mailto) and pure anchors (#...) are skipped; fenced code blocks and inline code are ignored.
Broken links inside analisis_cmb/bloque_*/antes/ (historical snapshots, not edited) are listed separately.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INLINE_CODE = re.compile(r"`[^`]*`")
LINK = re.compile(r"\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")


def md_files():
    out = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard", "*.md"],
                         cwd=ROOT, capture_output=True, text=True, check=True).stdout.split()
    return sorted(set(out))


def main() -> int:
    files = md_files()
    n_links, broken = 0, []
    per_file = {}
    for f in files:
        in_code = False
        count = 0
        for i, line in enumerate((ROOT / f).read_text(encoding="utf-8").splitlines(), 1):
            if line.lstrip().startswith("```"):
                in_code = not in_code
                continue
            if in_code:
                continue
            for target in LINK.findall(INLINE_CODE.sub("", line)):
                if re.match(r"^(https?:|mailto:|#)", target):
                    continue
                count += 1
                path = target.split("#", 1)[0]
                if not (ROOT / f).parent.joinpath(path).resolve().exists():
                    broken.append((f, i, target))
        per_file[f] = count
        n_links += count
    snap = [b for b in broken if "/antes/" in b[0]]
    broken = [b for b in broken if "/antes/" not in b[0]]
    print(f"Ficheros .md revisados: {len(files)}; enlaces relativos: {n_links}; "
          f"rotos: {len(broken)} (más {len(snap)} en snapshots históricos)\n")
    for title, rows in (("Rotos", broken), ("Rotos en snapshots históricos (no se editan)", snap)):
        if rows:
            print(f"{title}:\n\n| Fichero | Línea | Enlace |\n|:--|--:|:--|")
            for f, i, t in rows:
                print(f"| `{f}` | {i} | `{t}` |")
            print()
    print("| Fichero | Enlaces relativos |\n|:--|--:|")
    for f, c in per_file.items():
        if c:
            print(f"| `{f}` | {c} |")
    return 1 if broken else 0


if __name__ == "__main__":
    sys.exit(main())
