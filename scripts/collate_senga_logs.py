#!/usr/bin/env python3
"""Extract OPS/SENGA timings from Archer2 sweep logs into CSV + markdown.

Usage:
  python3 scripts/collate_senga_logs.py [log_dir]

Default log_dir is ./results next to this script, or $PWD/results.
Copy the resulting summary.csv / summary.md (or the tarball) back for analysis.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path


def grab(text: str, pat: str, default=""):
    m = re.search(pat, text)
    return m.group(1).strip() if m else default


def parse_log(path: Path) -> dict:
    t = path.read_text(errors="replace")
    d = {
        "file": path.name,
        "config": path.stem,
        "ok": "total_ttime" in t,
        "nxglbl": grab(t, r"nxglbl=\s*([0-9]+)"),
        "exec1": "",
        "exec2": "",
        "dump": grab(t, r"SENGA_DUMP_LAST=\s+(\S+)"),
        "total": grab(t, r"total_ttime\s*=\s*([0-9.]+)"),
        "rhscal": grab(t, r"rhscal_ttime\s*=\s*([0-9.]+)"),
        "rhsvel": grab(t, r"rhsvel_ttime\s*=\s*([0-9.]+)"),
        "parfer": grab(t, r"parfer_ttime\s*=\s*([0-9.]+)"),
        "indata": grab(t, r"indata_ttime\s*=\s*([0-9.]+)"),
        "kernel": grab(t, r"Total kernel time:\s*([0-9.eE+-]+)"),
        "tiled_halo": grab(t, r"Total tiled halo exchange time:\s*([0-9.eE+-]+)"),
        "user_halo": grab(t, r"Total user halo exchange time:\s*([0-9.eE+-]+)"),
        "msgs": grab(t, r"Average number of MPI messages per process:\s*([0-9.]+)"),
        "mb": grab(t, r"Average volume of MPI communications per process:\s*([0-9.]+)"),
        "energy": grab(t, r"Total CPU energy consumed \(RAPL\):\s*([0-9.eE+-]+)"),
        "ispec": "",
        "srun_sec": grab(t, r"elapsed[^\n]*?([0-9]+:[0-9]+:[0-9]+|[0-9]+m[0-9.]+s)"),
    }
    for line in t.splitlines():
        if line.strip().startswith("SENGA_OPS_EXEC"):
            if not d["exec1"]:
                d["exec1"] = " ".join(line.split())
            elif not d["exec2"]:
                d["exec2"] = " ".join(line.split())
            if "ispec=" in line:
                im = re.search(r"ispec=\s*(\S+)", line)
                if im:
                    d["ispec"] = im.group(1)
    try:
        k = float(d["kernel"]) if d["kernel"] else 0.0
        h = float(d["tiled_halo"]) if d["tiled_halo"] else 0.0
        d["kernel_plus_halo"] = f"{k + h:.4f}" if (d["kernel"] or d["tiled_halo"]) else ""
    except ValueError:
        d["kernel_plus_halo"] = ""
    return d


COLS = [
    "config",
    "ok",
    "nxglbl",
    "ispec",
    "total",
    "kernel",
    "tiled_halo",
    "kernel_plus_halo",
    "user_halo",
    "msgs",
    "mb",
    "rhscal",
    "rhsvel",
    "indata",
    "file",
]


def write_reports(rows, outdir: Path):
    csv_path = outdir / "summary.csv"
    md_path = outdir / "summary.md"
    with csv_path.open("w") as f:
        f.write(",".join(COLS) + "\n")
        for r in rows:
            f.write(",".join(str(r.get(c, "")) for c in COLS) + "\n")

    lines = [
        "# SENGA2 Archer2 sweep summary",
        "",
        f"Logs: `{outdir}` ({len(rows)} files)",
        "",
        "| " + " | ".join(COLS) + " |",
        "|" + "|".join("---" for _ in COLS) + "|",
    ]
    for r in rows:
        lines.append("| " + " | ".join(str(r.get(c, "")) for c in COLS) + " |")
    lines += [
        "",
        "## Flag echo",
        "",
    ]
    for r in rows:
        if r.get("exec1"):
            lines.append(f"- `{r['config']}`: {r['exec1']}")
            if r.get("exec2"):
                lines.append(f"  {r['exec2']}")
    lines.append("")
    md_path.write_text("\n".join(lines))
    return csv_path, md_path


def main():
    if len(sys.argv) > 1:
        logdir = Path(sys.argv[1])
    else:
        here = Path(__file__).resolve().parent
        cand = Path.cwd() / "results"
        logdir = cand if cand.is_dir() else here.parent / "results"
    logs = sorted(logdir.glob("*.out")) + sorted(logdir.glob("*.txt"))
    logs = [p for p in logs if p.name not in ("summary.md",) and "slurm" not in p.name]
    if not logs:
        sys.exit(f"no log files in {logdir}")
    rows = [parse_log(p) for p in logs]
    csv_path, md_path = write_reports(rows, logdir)
    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")
    print(f"{'config':28} {'ok':5} {'nx':>5} {'total':>10} {'kernel':>10} {'halo':>10} {'msgs':>8} {'MB':>10}")
    for r in rows:
        print(
            f"{r['config'][:28]:28} {str(r['ok']):5} {r['nxglbl']:>5} "
            f"{r['total']:>10} {r['kernel']:>10} {r['tiled_halo']:>10} "
            f"{r['msgs']:>8} {r['mb']:>10}"
        )


if __name__ == "__main__":
    main()
