import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from style import COPPER, DATA, GOLD, INK, MUTED, NAVY, TEAL, TEAL_SOFT, save, summary


def main():
    s = summary()
    direct = np.genfromtxt(os.path.join(DATA, "direct_milestone.csv"), delimiter=",", names=True)
    se = np.genfromtxt(os.path.join(DATA, "schrodinger_ho.csv"), delimiter=",", names=True)
    mask = (direct["x"] >= -0.02) & (direct["x"] <= 2.02)
    d = direct[mask]
    sm = (se["x"] >= -0.05) & (se["x"] <= 2.05)
    sch = se[sm]

    fig = plt.figure(figsize=(8.4, 7.1))
    grid = GridSpec(2, 1, height_ratios=[2.55, 1.0], hspace=0.08)
    ax = fig.add_subplot(grid[0])
    axr = fig.add_subplot(grid[1], sharex=ax)

    ax.fill_between(d["x"], d["psi2"] - d["err"], d["psi2"] + d["err"], color=TEAL_SOFT, linewidth=0, label="Monte Carlo uncertainty")
    ax.plot(d["x"], d["exact"], color=NAVY, lw=1.7, label="exact sum on this lattice")
    ax.plot(sch["x"], sch["psi2"], color=INK, lw=2.0, label="Schrödinger")
    ax.plot(sch["x"], sch["analytic"], color=COPPER, lw=1.35, ls=(0, (5, 2)), label="analytic ground state")
    ax.plot(d["x"], d["psi2"], color=TEAL, lw=1.15, label="direct path integral")
    ax.axvline(1.0, color=MUTED, lw=0.7, ls=(0, (1, 2)))
    ax.plot([1.0], [s["direct_psi2_1"]], "o", color=TEAL, ms=6.5, zorder=6, markeredgecolor=INK, markeredgewidth=0.4)

    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, 0.64)
    ax.set_ylabel(r"$|\phi(x)|^2$")
    ax.tick_params(labelbottom=False)
    ax.legend(loc="upper right", fontsize=9)
    ax.set_title("Harmonic oscillator ground state,  x from 0 to 2")

    note = (
        "value at  x = 1\n"
        f"analytic       {s['analytic_psi2_1']:.5f}\n"
        f"Schrödinger    {s['schrodinger_psi2_1']:.5f}\n"
        f"exact lattice  {s['lattice_psi2_1']:.5f}\n"
        f"direct paths   {s['direct_psi2_1']:.5f} ± {s['direct_psi2_1_err']:.5f}"
    )
    ax.text(
        0.03, 0.05, note, transform=ax.transAxes, ha="left", va="bottom", fontsize=8.4,
        family="Consolas", color=INK,
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#f7f1e4", edgecolor="#d9cbb6", linewidth=0.6),
    )

    axr.axhline(0.0, color=MUTED, lw=0.7)
    axr.fill_between(d["x"], d["psi2"] - d["exact"] - d["err"], d["psi2"] - d["exact"] + d["err"], color=TEAL_SOFT, linewidth=0)
    axr.plot(d["x"], d["psi2"] - d["exact"], color=TEAL, lw=1.2, label="Monte Carlo − exact lattice")
    axr.plot(d["x"], d["finiteT"] - d["analytic"], color=GOLD, lw=1.3, label="finite-T continuum − ground state")
    axr.set_xlim(0.0, 2.0)
    axr.set_xlabel("x")
    axr.set_ylabel("difference")
    axr.legend(loc="lower left", fontsize=8.5)
    axr.set_title("The Monte Carlo hits the lattice. The lattice is still a short time T = 4.", fontsize=11, loc="left", color=MUTED, pad=6)

    fig.text(0.125, 0.005, r"N = 7 time slices,   T = 4,   $\epsilon$ = 4/7,   m = $\hbar$ = 1,   V = $x^2/2$", color=MUTED, fontsize=9)
    save(fig, "milestone.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
