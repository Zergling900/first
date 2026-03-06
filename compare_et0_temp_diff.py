#!/usr/bin/env python3
"""
Compare temperature series between:
- mdmodular_resume Et0 files
- w-w2 Et files

Compute DeltaT = T_mdmodular - T_ww2 and plot DeltaT vs time.

Defaults:
- md files: mdmodular_resume/Data/BCC.Et0.*.md3
- ww files: w-w2/Data/BCC.Et.*.md3
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot temperature difference between mdmodular_resume Et0 and w-w2 Et."
    )
    parser.add_argument(
        "--md-pattern",
        default="mdmodular_resume/Data/BCC.Et0.*.md3",
        help="Glob for mdmodular_resume Et0 files.",
    )
    parser.add_argument(
        "--ww-pattern",
        default="w-w2/Data/BCC.Et.*.md3",
        help="Glob for w-w2 Et files.",
    )
    parser.add_argument(
        "--method",
        choices=("interp", "exact"),
        default="interp",
        help="Time alignment method: interp (recommended) or exact.",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-6,
        help="Tolerance (fs) for exact time matching.",
    )
    parser.add_argument(
        "--out-plot",
        default="temp_diff_md_et0_minus_ww.png",
        help="Output figure path.",
    )
    parser.add_argument(
        "--out-csv",
        default="temp_diff_md_et0_minus_ww.csv",
        help="Output CSV path.",
    )
    return parser.parse_args()


def file_sort_key(path: Path) -> Tuple[int, str]:
    m = re.search(r"\.(\d+)\.md3$", path.name)
    if m:
        return (int(m.group(1)), path.name)
    return (10**9, path.name)


def collect_files(pattern: str) -> List[Path]:
    files = sorted(Path(".").glob(pattern), key=file_sort_key)
    return [p for p in files if p.is_file()]


def read_series(
    files: Sequence[Path],
    time_col: int,
    temp_col: int,
) -> Tuple[np.ndarray, np.ndarray]:
    # Keep the latest value for duplicated time points (later files overwrite earlier ones).
    temp_by_time = {}

    for p in files:
        with p.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                if len(parts) <= max(time_col, temp_col):
                    continue
                try:
                    t = float(parts[time_col])
                    temp = float(parts[temp_col])
                except ValueError:
                    continue
                temp_by_time[t] = temp

    if not temp_by_time:
        return np.array([]), np.array([])

    times = np.array(sorted(temp_by_time.keys()), dtype=float)
    temps = np.array([temp_by_time[t] for t in times], dtype=float)
    return times, temps


def align_interp(
    md_t: np.ndarray,
    md_temp: np.ndarray,
    ww_t: np.ndarray,
    ww_temp: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    t_min = max(md_t[0], ww_t[0])
    t_max = min(md_t[-1], ww_t[-1])
    mask = (md_t >= t_min) & (md_t <= t_max)
    t = md_t[mask]
    md = md_temp[mask]
    ww = np.interp(t, ww_t, ww_temp)
    return t, md, ww


def align_exact(
    md_t: np.ndarray,
    md_temp: np.ndarray,
    ww_t: np.ndarray,
    ww_temp: np.ndarray,
    tol: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    t_out: List[float] = []
    md_out: List[float] = []
    ww_out: List[float] = []

    for t, mt in zip(md_t, md_temp):
        idx = np.searchsorted(ww_t, t)
        cand = []
        if idx < ww_t.size:
            cand.append(idx)
        if idx > 0:
            cand.append(idx - 1)
        if not cand:
            continue
        best = min(cand, key=lambda i: abs(ww_t[i] - t))
        if abs(ww_t[best] - t) <= tol:
            t_out.append(t)
            md_out.append(mt)
            ww_out.append(ww_temp[best])

    return np.array(t_out), np.array(md_out), np.array(ww_out)


def save_csv(
    out_csv: Path,
    t: np.ndarray,
    md: np.ndarray,
    ww: np.ndarray,
    dtemp: np.ndarray,
) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8") as f:
        f.write("time_fs,T_md_K,T_ww_K,DeltaT_K\n")
        for a, b, c, d in zip(t, md, ww, dtemp):
            f.write(f"{a:.10f},{b:.10f},{c:.10f},{d:.10f}\n")


def plot_delta_t(out_plot: Path, t: np.ndarray, dtemp: np.ndarray, method: str) -> None:
    out_plot.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 5.5))
    plt.plot(t, dtemp, lw=1.0)
    plt.axhline(0.0, color="black", lw=0.8, ls="--")
    plt.xlabel("Time (fs)")
    plt.ylabel("Temperature Difference (K)\nDeltaT = T_mdmodular - T_w-w2")
    plt.title(f"Temperature Difference vs Time ({method})")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_plot, dpi=150)
    plt.close()


def main() -> None:
    args = parse_args()

    md_files = collect_files(args.md_pattern)
    ww_files = collect_files(args.ww_pattern)

    if not md_files:
        raise SystemExit(f"No md files found: {args.md_pattern}")
    if not ww_files:
        raise SystemExit(f"No w-w2 files found: {args.ww_pattern}")

    # md Et0 format: time_fs(0), ... T_inst_K(5)
    md_t, md_temp = read_series(md_files, time_col=0, temp_col=5)
    # w-w2 Et format: time(fs)(0), ... T(4)
    ww_t, ww_temp = read_series(ww_files, time_col=0, temp_col=4)

    if md_t.size == 0:
        raise SystemExit("No valid rows parsed from md files.")
    if ww_t.size == 0:
        raise SystemExit("No valid rows parsed from w-w2 files.")

    if args.method == "interp":
        t, md, ww = align_interp(md_t, md_temp, ww_t, ww_temp)
    else:
        t, md, ww = align_exact(md_t, md_temp, ww_t, ww_temp, args.tol)

    if t.size == 0:
        raise SystemExit("No overlapping/matching time points after alignment.")

    dtemp = md - ww

    out_plot = Path(args.out_plot)
    out_csv = Path(args.out_csv)
    save_csv(out_csv, t, md, ww, dtemp)
    plot_delta_t(out_plot, t, dtemp, args.method)

    print(f"md files: {len(md_files)}")
    print(f"w-w2 files: {len(ww_files)}")
    print(f"aligned points: {t.size}")
    print(f"time range: [{t.min():.6f}, {t.max():.6f}] fs")
    print(f"DeltaT stats (K): min={dtemp.min():.6f}, max={dtemp.max():.6f}, mean={dtemp.mean():.6f}")
    print(f"saved plot: {out_plot}")
    print(f"saved csv : {out_csv}")


if __name__ == "__main__":
    main()
