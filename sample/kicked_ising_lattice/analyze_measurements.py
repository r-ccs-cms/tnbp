#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parse measurement logs from output.dat and produce:
- data_meas_x.txt  (sorted by time step, then site)
- data_meas_z.txt  (sorted by time step, then site)
- data_meas_j.txt  (sorted by time step, then i, then j)
Optionally, generate simple heatmaps/line plots into ./figs.
"""

import re
import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


# ---------- Regex patterns ----------
# X: complex tuple "(real, imag)"; we store real & imag separately
PAT_X  = re.compile(
    r"measurement of X\s+at site\s+(\d+)\s+time step\s+(\d+)\s*=\s*\(\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?),\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*\)"
)
# Z: single real
PAT_Z  = re.compile(
    r"measurement of Z\s+at site\s+(\d+)\s+time step\s+(\d+)\s*=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)"
)
# ZZ: sites (i,j)
PAT_ZZ = re.compile(
    r"measurement of ZZ\s+at sites\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)\s+time step\s+(\d+)\s*=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)"
)


def parse_file(path: Path):
    xs, zs, zzs = [], [], []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            # X
            m = PAT_X.search(line)
            if m:
                site = int(m.group(1))
                step = int(m.group(2))
                real = float(m.group(3))
                imag = float(m.group(4))
                xs.append((step, site, real, imag, line.rstrip("\n")))
                continue

            # Z
            m = PAT_Z.search(line)
            if m:
                site = int(m.group(1))
                step = int(m.group(2))
                val  = float(m.group(3))
                zs.append((step, site, val, line.rstrip("\n")))
                continue

            # ZZ
            m = PAT_ZZ.search(line)
            if m:
                i    = int(m.group(1))
                j    = int(m.group(2))
                step = int(m.group(3))
                val  = float(m.group(4))
                zzs.append((step, i, j, val, line.rstrip("\n")))
                continue

    return xs, zs, zzs


def save_tables(xs, zs, zzs, outdir: Path):
    # ----- X -----
    if xs:
        df_x = pd.DataFrame(xs, columns=["time_step", "site", "real", "imag", "raw"])
        df_x.sort_values(["time_step", "site"], inplace=True)
        # Keep the original raw lines (spec matches your current logging)
        (outdir / "data_meas_x.txt").write_text("\n".join(df_x["raw"]), encoding="utf-8")

        # Also provide a clean TSV (optional, handy for analysis)
        df_x[["time_step", "site", "real", "imag"]].to_csv(outdir / "data_meas_x.tsv",
                                                           sep="\t", index=False)
    # ----- Z -----
    if zs:
        df_z = pd.DataFrame(zs, columns=["time_step", "site", "value", "raw"])
        df_z.sort_values(["time_step", "site"], inplace=True)
        (outdir / "data_meas_z.txt").write_text("\n".join(df_z["raw"]), encoding="utf-8")
        df_z[["time_step", "site", "value"]].to_csv(outdir / "data_meas_z.tsv",
                                                    sep="\t", index=False)
    # ----- ZZ -----
    if zzs:
        df_zz = pd.DataFrame(zzs, columns=["time_step", "i", "j", "value", "raw"])
        df_zz.sort_values(["time_step", "i", "j"], inplace=True)
        (outdir / "data_meas_j.txt").write_text("\n".join(df_zz["raw"]), encoding="utf-8")
        df_zz[["time_step", "i", "j", "value"]].to_csv(outdir / "data_meas_j.tsv",
                                                       sep="\t", index=False)


def make_plots(xs, zs, zzs, figdir: Path):
    figdir.mkdir(parents=True, exist_ok=True)

    # Heatmap helpers
    def heatmap(df, value_col, row="time_step", col="site", fname="heatmap.png", title=None):
        if df.empty:
            return
        pivot = df.pivot(index=row, columns=col, values=value_col).sort_index().sort_index(axis=1)
        plt.figure()
        plt.imshow(pivot.values, aspect="auto", origin="lower")
        plt.colorbar(label=value_col)
        plt.xlabel(col)
        plt.ylabel(row)
        if title:
            plt.title(title)
        plt.tight_layout()
        plt.savefig(figdir / fname, dpi=160)
        plt.close()

    # X heatmaps: real and imag
    if xs:
        dfx = pd.DataFrame(xs, columns=["time_step", "site", "real", "imag", "raw"])
        heatmap(dfx, "real", fname="X_real_heatmap.png", title="X (real)")
        heatmap(dfx, "imag", fname="X_imag_heatmap.png", title="X (imag)")

    # Z heatmap
    if zs:
        dfz = pd.DataFrame(zs, columns=["time_step", "site", "value", "raw"])
        heatmap(dfz, "value", fname="Z_heatmap.png", title="Z")

    # ZZ: if data looks like nearest-neighbor pairs (j == i+1), we can plot a heatmap over bonds
    if zzs:
        dfzz = pd.DataFrame(zzs, columns=["time_step", "i", "j", "value", "raw"])
        # Try to map (i,j) → bond index i when j==i+1
        df_bonds = dfzz[dfzz["j"] == dfzz["i"] + 1].copy()
        if not df_bonds.empty:
            df_bonds = df_bonds.rename(columns={"i": "bond"})
            heatmap(df_bonds, "value", col="bond", fname="ZZ_nn_heatmap.png", title="ZZ (nearest-neighbor)")
        else:
            # Fallback: simple line plot per time step over sorted pairs
            # (X-axis: pair index after sorting; just to have a quick visual)
            dfzz_sorted = dfzz.sort_values(["time_step", "i", "j"]).copy()
            for step, grp in dfzz_sorted.groupby("time_step"):
                plt.figure()
                xs_plot = range(len(grp))
                plt.plot(xs_plot, grp["value"])
                plt.xlabel("pair index (sorted by i,j)")
                plt.ylabel("ZZ value")
                plt.title(f"ZZ values at time step {step}")
                plt.tight_layout()
                plt.savefig(figdir / f"ZZ_pairs_t{step}.png", dpi=160)
                plt.close()


def main():
    ap = argparse.ArgumentParser(description="Parse output.dat and export sorted tables & plots.")
    ap.add_argument("input", nargs="?", default="output.dat", help="input log file (default: output.dat)")
    ap.add_argument("--no-plots", action="store_true", help="skip figure generation")
    ap.add_argument("--outdir", default=".", help="directory to write data_*.txt/tsv")
    ap.add_argument("--figdir", default="figs", help="directory to write figures")
    args = ap.parse_args()

    inpath = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    xs, zs, zzs = parse_file(inpath)
    save_tables(xs, zs, zzs, outdir)

    if not args.no_plots:
        make_plots(xs, zs, zzs, Path(args.figdir))

    print("Done.")
    print(f"- data_meas_x.txt / data_meas_x.tsv written to: {outdir.resolve()}")
    print(f"- data_meas_z.txt / data_meas_z.tsv written to: {outdir.resolve()}")
    print(f"- data_meas_j.txt / data_meas_j.tsv written to: {outdir.resolve()}")
    if not args.no_plots:
        print(f"- Figures written to: {Path(args.figdir).resolve()}")


if __name__ == "__main__":
    main()
