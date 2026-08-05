#!/usr/bin/env python3
"""
Metadynamics post-processing pipeline.
Produces three plots matching the conventions of the existing
plot_fes_new.py, plot_colvar.py, and plot_convergence.py scripts.
"""
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects


# ====================================================================
# Helpers  (extracted verbatim from the existing scripts)
# ====================================================================

def _read_fes_header(filename):
    """Parse PLUMED ``#! SET`` lines for grid metadata (plot_convergence style)."""
    meta = {}
    with open(filename) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#! SET"):
                parts = line.split()
                for k in ("min_phi", "max_phi", "nbins_phi",
                          "min_psi", "max_psi", "nbins_psi"):
                    if k in parts:
                        v = parts[parts.index(k) + 1].lower()
                        meta[k] = (-np.pi if v == "-pi" else
                                   np.pi if v == "pi" else float(v))
            elif line and not line.startswith("#"):
                break
    return meta


def _read_fes_grid(filename):
    """Return (phi_vals, psi_vals, F_2d) from a 3-column PLUMED FES file."""
    data = np.loadtxt(filename, comments="#")
    phi_vals = np.unique(data[:, 0])
    psi_vals = np.unique(data[:, 1])
    n_phi, n_psi = len(phi_vals), len(psi_vals)
    F = data[:, 2].reshape(n_psi, n_phi)          # psi rows, phi cols
    F -= np.nanmin(F)
    return phi_vals, psi_vals, F


def _project_1d(F, axis, kT=0.6):
    """Boltzmann-weighted 1D free-energy projection (plot_convergence style)."""
    beta = 1.0 / kT
    s = np.exp(-beta * F)
    s = np.sum(s, axis=axis)
    f = -np.log(s) / beta
    return f - np.nanmin(f)


def _pi_ticks(ax, x=True, y=False):
    """Set axis ticks to -pi, -pi/2, 0, pi/2, pi with LaTeX labels."""
    pi_labels = [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"]
    tick_pos = [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi]
    if x:
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(pi_labels)
    if y:
        ax.set_yticks(tick_pos)
        ax.set_yticklabels(pi_labels)


# ====================================================================
# Plot 1 — 2D FES with global-minimum marker  (plot_fes_new style)
# ====================================================================

def plot_fes_labeled(fes_file, outpath):
    """
    Reads a 3-column FES (phi psi energy), finds global minimum,
    draws a filled contour with 40 levels / viridis, marks the minimum
    with a red dot + LaTeX label, pi ticks.  Saved at 300 dpi.
    """
    print(f"Loading FES data from {fes_file} ...")
    data = np.loadtxt(fes_file, comments="#")

    phi_raw = data[:, 0]
    psi_raw = data[:, 1]
    fes_raw = data[:, 2]

    min_idx = np.argmin(fes_raw)
    min_phi = phi_raw[min_idx]
    min_psi = psi_raw[min_idx]
    raw_min_energy = fes_raw[min_idx]
    print(f"  Global minimum: phi={min_phi:.4f} rad, psi={min_psi:.4f} rad")

    fes_raw = fes_raw - raw_min_energy          # shift so min = 0

    n_phi = len(np.unique(phi_raw))
    n_psi = len(np.unique(psi_raw))
    print(f"  Grid: {n_phi}x{n_psi}")

    PHI = phi_raw.reshape(n_psi, n_phi)
    PSI = psi_raw.reshape(n_psi, n_phi)
    FES = fes_raw.reshape(n_psi, n_phi)

    plt.figure(figsize=(8, 6), dpi=100)
    contour = plt.contourf(PHI, PSI, FES, levels=40, cmap="viridis")

    plt.scatter([min_phi], [min_psi], color="red", marker="o", s=100,
                edgecolor="black", linewidth=1.5, zorder=10,
                label="Global minimum")

    label_text = f"$\\phi$: {min_phi:.2f}\n$\\psi$: {min_psi:.2f}"
    offset_x, offset_y = 0.15, 0.15
    text_x = min_phi + offset_x if min_phi < np.pi / 2 else min_phi - offset_x * 4
    text_y = min_psi + offset_y if min_psi < np.pi / 2 else min_psi - offset_y * 3
    ha = "left" if min_phi < np.pi / 2 else "right"

    txt = plt.text(text_x, text_y, label_text, color="white", fontsize=11,
                   fontweight="bold", zorder=11, ha=ha, va="bottom")
    txt.set_path_effects([path_effects.withStroke(linewidth=3,
                                                   foreground="black")])

    plt.legend(loc="upper right", framealpha=0.9, fontsize=10)

    cbar = plt.colorbar(contour)
    cbar.set_label("Free Energy (kJ/mol)", fontsize=12, fontweight="bold")

    plt.title("Free Energy Surface (Ramachandran Plot)", fontsize=14,
              fontweight="bold")
    plt.xlabel(r"$\phi$ (rad)", fontsize=12)
    plt.ylabel(r"$\psi$ (rad)", fontsize=12)

    _pi_ticks(plt.gca(), x=True, y=True)
    plt.xlim(-np.pi, np.pi)
    plt.ylim(-np.pi, np.pi)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {outpath}")


# ====================================================================
# Plot 2 — COLVAR visualisation (plot_colvar style)
# ====================================================================

def plot_colvar_panel(colvar_file, outpath):
    """
    Reads a PLUMED COLVAR, expects columns (time phi psi bias).
    Draws a 2×2 panel: phi vs t, psi vs t, 2D hexbin density, bias vs t.
    """
    print(f"Loading COLVAR from {colvar_file} ...")
    data = np.loadtxt(colvar_file, comments=["#", "@"])
    t = data[:, 0]
    phi = data[:, 1]
    psi = data[:, 2]
    bias = data[:, 3]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    axes[0, 0].plot(t, phi, "b-", lw=0.5)
    axes[0, 0].set_ylabel(r"$\phi$ (rad)")
    axes[0, 0].set_xlabel("Time (ps)")
    axes[0, 0].set_title(r"$\phi$ vs Time")

    axes[1, 0].plot(t, psi, "r-", lw=0.5)
    axes[1, 0].set_ylabel(r"$\psi$ (rad)")
    axes[1, 0].set_xlabel("Time (ps)")
    axes[1, 0].set_title(r"$\psi$ vs Time")

    hb = axes[0, 1].hexbin(phi, psi, gridsize=100, cmap="viridis")
    axes[0, 1].set_xlabel(r"$\phi$ (rad)")
    axes[0, 1].set_ylabel(r"$\psi$ (rad)")
    axes[0, 1].set_title("2D Density")
    plt.colorbar(hb, ax=axes[0, 1])

    axes[1, 1].plot(t, bias, "g-", lw=0.5)
    axes[1, 1].set_ylabel("Metadynamics Bias")
    axes[1, 1].set_xlabel("Time (ps)")
    axes[1, 1].set_title("Bias vs Time")

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {outpath}")


# ====================================================================
# Plot 3 — FES convergence  (plot_convergence style)
# ====================================================================

def _stride_sort_key(path):
    """Extract stride integer from filename:  fes_0  or  fes_1000.dat  etc."""
    base = os.path.basename(path)
    try:
        num_part = base.split("_", 1)[1]
        num_part = num_part.split(".")[0]
        return int(num_part)
    except (IndexError, ValueError):
        return 0


def plot_convergence_panel(pattern, outpath):
    """
    Finds FES stride files matching *pattern*, reads them with
    _read_fes_grid + _read_fes_header, produces a 2×3 panel:

      [phi 1D evolution]  [phi RMSD]  [final 2D FES]
      [psi 1D evolution]  [psi RMSD]  [initial 2D FES]
    """
    files = sorted(glob.glob(pattern), key=_stride_sort_key)
    if not files:
        print(f"  No convergence files matching '{pattern}'. Skipping.")
        return
    print(f"  Found {len(files)} stride file(s)")

    phi, psi, F0 = _read_fes_grid(files[0])
    n = len(files)

    stride_idx = []
    F_phi_all, F_psi_all = [], []
    rmsd_phi, rmsd_psi = [0.0], [0.0]

    for i, f in enumerate(files):
        sid = _stride_sort_key(f)
        stride_idx.append(sid)
        _pv, _psiv, F = _read_fes_grid(f)
        F_phi = _project_1d(F, axis=0, kT=0.6)
        F_psi = _project_1d(F, axis=1, kT=0.6)
        F_phi_all.append(F_phi)
        F_psi_all.append(F_psi)
        if i > 0:
            rmsd_phi.append(
                np.sqrt(np.nanmean((F_phi - F_phi_all[i - 1]) ** 2)))
            rmsd_psi.append(
                np.sqrt(np.nanmean((F_psi - F_psi_all[i - 1]) ** 2)))

    F_phi_all = np.array(F_phi_all)
    F_psi_all = np.array(F_psi_all)

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # --- phi 1D evolution ---
    ax = axes[0, 0]
    step = max(1, n // 20)
    for i in range(0, n, step):
        lbl = (f"step {stride_idx[i]}"
               if i == 0 or i == n - 1 or i == n // 2 else "")
        ax.plot(phi, F_phi_all[i], lw=0.8, alpha=0.6, label=lbl)
    ax.plot(phi, F_phi_all[-1], "k-", lw=2.5,
            label=f"final ({stride_idx[-1]})")
    ax.set_xlabel(r"$\phi$ (rad)")
    ax.set_ylabel("FE (kJ/mol)")
    ax.set_title(r"$\phi$ 1D Profile Convergence")

    # --- psi 1D evolution ---
    ax = axes[1, 0]
    for i in range(0, n, step):
        ax.plot(psi, F_psi_all[i], lw=0.8, alpha=0.6)
    ax.plot(psi, F_psi_all[-1], "k-", lw=2.5,
            label=f"final ({stride_idx[-1]})")
    ax.set_xlabel(r"$\psi$ (rad)")
    ax.set_ylabel("FE (kJ/mol)")
    ax.set_title(r"$\psi$ 1D Profile Convergence")

    # --- RMSD convergence ---
    ax = axes[0, 1]
    ax.plot(stride_idx, rmsd_phi, "b-o", ms=3)
    ax.set_xlabel("Stride step")
    ax.set_ylabel("RMSD (kJ/mol)")
    ax.set_title(r"$\phi$ RMSD between consecutive blocks")
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    ax.plot(stride_idx, rmsd_psi, "r-o", ms=3)
    ax.set_xlabel("Stride step")
    ax.set_ylabel("RMSD (kJ/mol)")
    ax.set_title(r"$\psi$ RMSD between consecutive blocks")
    ax.grid(alpha=0.3)

    # --- final / initial 2D FES ---
    Phi, Psi = np.meshgrid(phi, psi)

    _, _, F_final = _read_fes_grid(files[-1])
    for row, (title, F_mat) in enumerate([
        (f"Final FES (step {stride_idx[-1]})", F_final),
        (f"Initial FES (step {stride_idx[0]})", F0),
    ]):
        ax = axes[row, 2]
        cf = ax.contourf(Phi, Psi, F_mat, levels=50, cmap="viridis")
        ax.set_xlabel(r"$\phi$ (rad)")
        ax.set_ylabel(r"$\psi$ (rad)")
        ax.set_title(title)
        plt.colorbar(cf, ax=ax, label="kJ/mol")

    # --- pi tick labels on relevant axes ---
    for ax_row in axes:
        for ax in ax_row:
            xlab = ax.get_xlabel() or ""
            ylab = ax.get_ylabel() or ""
            if xlab and "phi" in xlab.lower():
                _pi_ticks(ax, x=True, y=False)
            if ylab and ("psi" in ylab.lower() or "phi" in ylab.lower()):
                _pi_ticks(ax, x=False, y=True)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {outpath}")


# ====================================================================
# Main
# ====================================================================

def main():
    p = argparse.ArgumentParser(description="Metadynamics post-processing")

    p.add_argument("--colvar", help="PLUMED COLVAR file "
                   "(skips plot_colvar_panel if omitted)")
    p.add_argument("--fes2d", help="3-column FES (from sum_hills) "
                   "→ plot_fes_labeled")
    p.add_argument("--fes-pattern", default="tmp/fes_[0-9]*",
                   help="Glob for convergence stride files "
                   "(default: tmp/fes_[0-9]*)")
    p.add_argument("--outdir", default="results",
                   help="Output directory for PNG files")

    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out = lambda stem: os.path.join(args.outdir, stem)

    if args.colvar and os.path.isfile(args.colvar):
        plot_colvar_panel(args.colvar, out("meta_colvar.png"))

    if args.fes2d and os.path.isfile(args.fes2d):
        plot_fes_labeled(args.fes2d, out("meta_fes2d.png"))

    plot_convergence_panel(args.fes_pattern, out("meta_convergence.png"))


if __name__ == "__main__":
    main()
