#!/usr/bin/env python3
"""Make the README figures from the output of scripts/run_paper.sh.

    python scripts/make_figures.py [rundir] [figdir]

defaults: rundir = data/runs, figdir = docs/figures
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# palette (light surface): diverging blue <-> red with a neutral ink midpoint,
# sequential blue ramp for the ordered size ratio
INK, INK2, GRID, AXIS = "#191919", "#666666", "#e7e7e7", "#c2c2c2"
BLUE, RED, ORANGE = "#257de0", "#ec4647", "#f47230"
RAMP = ["#84b9f3", "#338ced", "#2370cd", "#1757ac", "#0e3d7e"]

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK,
    "axes.titlesize": 11,
    "axes.titlecolor": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "legend.frameon": False,
    "legend.fontsize": 9.5,
    "lines.linewidth": 2.0,
    "mathtext.fontset": "dejavusans",
})


@dataclass
class Run:
    q: float
    eps2: float
    gamma: float
    rho2l: float
    rho2v: float
    profile: np.ndarray = field(repr=False)
    veff: np.ndarray = field(repr=False)

    @property
    def symmetric(self) -> bool:
        return abs(self.rho2v / self.rho2l - 1.0) < 1e-6


def _header(path: Path) -> dict[str, float]:
    head = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        for key, val in re.findall(r"([\w*^/ ]+?)\s*=\s*([-+0-9.eE]+)", line[1:]):
            head[key.strip()] = float(val)
    return head


def load_runs(rundir: Path) -> list[Run]:
    runs = []
    for prof in sorted(rundir.glob("profile_q*_eps*.dat")):
        veff = prof.with_name(prof.name.replace("profile_", "veff_"))
        h = _header(prof)
        runs.append(Run(q=h["sigma2/sigma1"], eps2=h["beta*eps2"], gamma=h["beta*gamma*s1^2"],
                        rho2l=h["rho2^l s1^3"], rho2v=h["rho2^v s1^3"],
                        profile=np.loadtxt(prof), veff=np.loadtxt(veff)))
    if not runs:
        sys.exit(f"no runs found in {rundir}; run scripts/run_paper.sh first")
    return runs


def triplet(runs: list[Run], q: float) -> list[Run]:
    """The runs below, at and above the symmetric eps2 for one size ratio."""
    sel = sorted((r for r in runs if abs(r.q - q) < 1e-9), key=lambda r: r.eps2)
    sym = next(i for i, r in enumerate(sel) if r.symmetric)
    return [sel[sym - 1], sel[sym], sel[sym + 1]]


def style_triplet(i: int) -> dict:
    return [dict(color=BLUE, ls="--"), dict(color=INK, ls="-"), dict(color=RED, ls=(0, (1, 1.5)))][i]


def fig_interface(runs: list[Run], out: Path) -> None:
    r = next(r for r in runs if abs(r.q - 1) < 1e-9 and abs(r.eps2 - 0.6) < 1e-9)
    z, rho1, rho2 = r.profile[:, 0], r.profile[:, 1], r.profile[:, 2]

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 5.0), sharex=True, gridspec_kw=dict(hspace=0.38))
    a1.plot(z, rho1, color=INK)
    a1.set_ylabel(r"$\rho_1(z)\,\sigma_1^3$")
    a1.set_ylim(0, 0.68)
    a1.text(-6.6, 0.47, "vapor", color=INK2)
    a1.text(4.8, 0.47, "liquid", color=INK2)
    a1.set_title(r"solvent ($\beta\epsilon_1=1,\ \lambda_1=1.5$)", loc="left")

    a2.plot(z, rho2 / r.rho2l, color=BLUE)
    a2.axhline(1, color=AXIS, lw=0.8, zorder=0)
    a2.set_ylabel(r"$\rho_2(z)/\rho_{2,0}$")
    a2.set_xlabel(r"$z/\sigma_1$")
    a2.set_title(rf"solute ($\sigma_2/\sigma_1=1,\ \beta\epsilon_2={r.eps2:g}$), "
                 r"$\rho_{2,0}\sigma_2^3=6\cdot10^{-6}/\pi$", loc="left")
    a2.set_xlim(-7, 7)
    for a in (a1, a2):
        a.axvline(0, color=AXIS, lw=0.8, zorder=0)
    fig.savefig(out / "interface.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_adsorption(runs: list[Run], out: Path) -> None:
    qs = [1.0, 1.5, 2.0, 2.5, 3.0]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.0), gridspec_kw=dict(hspace=0.42, wspace=0.28))
    for ax, q in zip(axes.flat, qs):
        trip = triplet(runs, q)
        for i, r in enumerate(trip):
            ax.plot(r.profile[:, 0], r.profile[:, 2] / r.rho2l, **style_triplet(i))
        ax.set_yscale("log")
        lo, hi = ax.get_ylim()
        if hi / lo < 300:  # narrow range: label 1-2-5 steps in plain numbers
            ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        else:
            ax.yaxis.set_major_locator(LogLocator(base=10, numticks=12))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_xlim(-7.5, 7.5)
        eps = " / ".join(f"{r.eps2:.4g}" for r in trip)
        ax.set_title(rf"$\sigma_2/\sigma_1={q:g}$,   $\beta\epsilon_2$ = {eps}", loc="left")
    for ax in axes[1]:
        ax.set_xlabel(r"$z/\sigma_1$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\rho_2(z)/\rho_{2,0}$")
    handles = [Line2D([], [], **style_triplet(i)) for i in range(3)]
    fig.legend(handles, [r"$\epsilon_2$ below symmetric (solvophobic, prefers vapor)",
                         r"symmetric $\epsilon_2$ ($\rho_2^l=\rho_2^v$)",
                         r"$\epsilon_2$ above symmetric (solvophilic, prefers liquid)"],
               loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.02), fontsize=10.5)

    # last panel: peak enrichment at the interface for the symmetric solute
    ax = axes[1, 2]
    sym = sorted((r for r in runs if r.symmetric), key=lambda r: r.q)
    peak = [np.max(r.profile[:, 2]) / r.rho2l for r in sym]
    ax.plot([r.q for r in sym], peak, color=INK, marker="o", ms=7, mfc=INK, mec="white", mew=1.5)
    for r, p in zip(sym, peak):
        ax.annotate(f"{p:,.0f}" if p >= 10 else f"{p:.2f}", (r.q, p), textcoords="offset points",
                    xytext=(-8, 4), ha="right", va="bottom", color=INK2, fontsize=9)
    ax.set_yscale("log")
    ax.set_xlabel(r"$\sigma_2/\sigma_1$")
    ax.set_ylabel(r"$\max_z\,\rho_2(z)/\rho_{2,0}$")
    ax.set_title("peak enrichment, symmetric case", loc="left")
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=12))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(1, 1e5)
    ax.set_xlim(0.6, 3.2)
    fig.savefig(out / "adsorption.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_veff(runs: list[Run], out: Path) -> None:
    sel = sorted((r for r in runs if abs(r.q - 2) < 1e-9 and (r.symmetric or abs(r.eps2 - 1.5) < 1e-9)),
                 key=lambda r: r.eps2)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), sharey=True, gridspec_kw=dict(wspace=0.08))
    for ax, r, tag in zip(axes, sel, ["(a)", "(b)"]):
        z, vmix, vwid = r.veff[:, 0], r.veff[:, 1], r.veff[:, 2]
        ax.plot(z, vwid, color=INK, label=r"Widom insertion, $\rho_{2,0}\to0$", zorder=2)
        ax.plot(z[::24], vmix[::24], ls="none", marker="o", ms=6.5, mfc=ORANGE, mec="white", mew=1.2,
                label=r"binary mixture, $-\ln[\rho_2(z)/\rho_{2,0}]$", zorder=3)
        ax.set_xlim(-7, 7)
        ax.set_xlabel(r"$z/\sigma_1$")
        ax.set_title(rf"{tag}  $\sigma_2/\sigma_1=2,\ \beta\epsilon_2={r.eps2:.4g}$"
                     + (" (symmetric)" if r.symmetric else ""), loc="left")
        ax.axhline(0, color=AXIS, lw=0.8, zorder=0)
    axes[0].set_ylabel(r"$\beta V_{\mathrm{eff}}(z)$")
    axes[1].legend(loc="upper right")
    fig.savefig(out / "effective_potential.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_veff_sizes(runs: list[Run], out: Path) -> None:
    sym = sorted((r for r in runs if r.symmetric), key=lambda r: r.q)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for r, c in zip(sym, RAMP):
        ax.plot(r.veff[:, 0], r.veff[:, 2], color=c,
                label=rf"$\sigma_2/\sigma_1={r.q:g}$, $\beta\epsilon_2={r.eps2:.4f}$")
    ax.axhline(0, color=AXIS, lw=0.8, zorder=0)
    ax.set_xlim(-7, 7)
    ax.set_xlabel(r"$z/\sigma_1$")
    ax.set_ylabel(r"$\beta V_{\mathrm{eff}}(z)$")
    ax.set_title("Widom route, symmetric solutes", loc="left")
    ax.legend(loc="lower left")
    fig.savefig(out / "effective_potential_sizes.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    rundir = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "data" / "runs"
    figdir = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "docs" / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    runs = load_runs(rundir)

    fig_interface(runs, figdir)
    fig_adsorption(runs, figdir)
    fig_veff(runs, figdir)
    fig_veff_sizes(runs, figdir)

    print(f"{'q':>5} {'beta eps2':>11} {'sym':>4} {'beta gamma s1^2':>16} {'max rho2/rho20':>15}")
    for r in sorted(runs, key=lambda r: (r.q, r.eps2)):
        print(f"{r.q:5.2f} {r.eps2:11.6f} {'*' if r.symmetric else '':>4} {r.gamma:16.6f} "
              f"{np.max(r.profile[:, 2]) / r.rho2l:15.4g}")


if __name__ == "__main__":
    main()
