import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from style import COPPER, DATA, INK, MUTED, NAVY, TEAL, save, summary


def main():
    s = summary()
    hist = np.genfromtxt(os.path.join(DATA, "metro_ho_hist.csv"), delimiter=",", names=True)
    scan = np.genfromtxt(os.path.join(DATA, "hbar_scan.csv"), delimiter=",", names=True)
    pots = np.genfromtxt(os.path.join(DATA, "potentials.csv"), delimiter=",", names=True, dtype=None, encoding=None)

    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    x = np.linspace(-0.985, 0.985, 700)
    classical = 1.0 / (np.pi * np.sqrt(1.0 - x * x))
    ax.axvspan(-2.4, -1.0, color="#eadccb", zorder=0)
    ax.axvspan(1.0, 2.4, color="#eadccb", zorder=0)
    ax.plot(x, classical, color=COPPER, lw=1.7, label="classical, energy 1/2")
    ax.fill_between(hist["x"], hist["psi2"] - hist["err"], hist["psi2"] + hist["err"], color="#c5ddd9", linewidth=0, zorder=2)
    ax.plot(hist["x"], hist["psi2"], color=TEAL, lw=1.8, label="path integral", zorder=3)
    xs = np.linspace(-2.4, 2.4, 500)
    ax.plot(xs, np.exp(-xs * xs) / np.sqrt(np.pi), color=INK, lw=1.2, ls=(0, (5, 2)), label="analytic ground state", zorder=4)
    ax.set_xlim(-2.4, 2.4)
    ax.set_ylim(0.0, 1.15)
    ax.set_xlabel("x")
    ax.set_ylabel("probability density")
    ax.set_title("Classical particle piles up at the edges. The quantum one does not.")
    ax.legend(loc="upper right", fontsize=9)
    ax.text(-2.28, 0.08, "forbidden", color=MUTED, fontsize=8)
    ax.text(1.55, 0.08, "forbidden", color=MUTED, fontsize=8)
    ax.annotate(
        "turning\npoints",
        xy=(1.0, 0.9),
        xytext=(1.35, 0.72),
        fontsize=8,
        color=COPPER,
        arrowprops=dict(arrowstyle="-", color=COPPER, lw=0.6),
    )
    fig.text(
        0.125, -0.03,
        f"Chance of being past |x| = 1:  path integral {100 * s['forbidden_prob']:.1f}%,   "
        f"exact ground state {100 * s['forbidden_analytic']:.1f}%.  A classical particle with this energy never gets there.",
        color=MUTED, fontsize=9,
    )
    save(fig, "classical.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    hbars = np.unique(scan["hbar"])
    colors = LinearSegmentedColormap.from_list("hb", ["#6e2414", "#c24d2c", "#a97820", "#0c6e6e", "#243656", "#6a3a58"])
    offset = 0.95
    for i, hb in enumerate(hbars):
        sl = scan[np.isclose(scan["hbar"], hb)]
        base = i * offset
        color = colors(i / max(1, len(hbars) - 1))
        ax.fill_between(sl["x"], base, base + sl["psi2"], color=color, alpha=0.88, linewidth=0)
        ax.plot(sl["x"], base + sl["analytic"], color=INK, lw=0.8, alpha=0.75)
        ax.text(-6.7, base + 0.08, rf"$\hbar$ = {hb:g}", color=color, fontsize=10, va="bottom")
    ax.set_xlim(-6.9, 6.5)
    ax.set_yticks([])
    ax.set_xlabel("x")
    ax.set_title(r"Turn $\hbar$ down and the particle is squeezed toward the origin")
    fig.text(0.125, 0.01, r"Filled shape is the Metropolis density. The thin line is the exact ground state at that $\hbar$. Same potential, same mass.", color=MUTED, fontsize=9)
    for spine in ("left",):
        ax.spines[spine].set_visible(False)
    save(fig, "hbar.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.5))
    specs = [
        ("anharm", r"$V=\frac{1}{2}x^2+\frac{1}{4}x^4$", (-3.2, 3.2)),
        ("doublewell", r"$V=0.8\,(x^2-2.25)^2$", (-3.2, 3.2)),
        ("morse", r"$V=3\,(1-e^{-0.65x})^2$", (-1.5, 5.5)),
    ]
    for ax, (name, formula, xlim) in zip(axes, specs):
        sl = pots[pots["name"] == name]
        m = (sl["x"] >= xlim[0]) & (sl["x"] <= xlim[1])
        sl = sl[m]
        twin = ax.twinx()
        twin.fill_between(sl["x"], sl["V"], color="#eadccb", linewidth=0, alpha=0.95, zorder=0)
        energy = float(sl["energy"][0])
        twin.axhline(energy, color=COPPER, lw=1.0, ls=(0, (3, 2)))
        twin.set_ylim(0.0, max(float(sl["V"].max()) * 1.05, energy * 1.4))
        twin.set_ylabel("V(x)", color=MUTED, fontsize=9)
        twin.tick_params(colors=MUTED, labelsize=8)
        twin.spines["right"].set_color("#d9d0c2")
        ax.plot(sl["x"], sl["se_psi2"], color=INK, lw=1.6, label="Schrödinger", zorder=3)
        ax.fill_between(sl["x"], sl["psi2"] - sl["err"], sl["psi2"] + sl["err"], color="#c5ddd9", linewidth=0, zorder=2)
        ax.plot(sl["x"], sl["psi2"], color=TEAL, lw=1.3, label="path integral", zorder=4)
        ax.set_xlim(*xlim)
        ax.set_ylim(0.0, max(float(sl["se_psi2"].max()), float(sl["psi2"].max())) * 1.25)
        ax.set_xlabel("x")
        ax.set_title(formula, fontsize=11)
        twin.set_zorder(0)
        ax.set_zorder(1)
        ax.patch.set_alpha(0.0)
    axes[0].set_ylabel(r"$|\phi(x)|^2$")
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle("The same engine, three other potentials", fontsize=14, x=0.02, ha="left")
    fig.text(0.02, -0.02, "Shaded hill is the potential. Dashed copper line is the ground-state energy on that scale. Teal is Metropolis, black is the Schrödinger solution of the same V.", color=MUTED, fontsize=9)
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    save(fig, "potentials.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
