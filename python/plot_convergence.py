import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

from style import DATA, GOLD, INK, MUTED, NAVY, TEAL, save, summary


def gauss(x, alpha):
    return np.sqrt(alpha / np.pi) * np.exp(-alpha * x * x)


def main():
    s = summary()
    table = np.genfromtxt(os.path.join(DATA, "exact_alpha.csv"), delimiter=",", names=True)
    Ns = np.array(sorted(set(table["N"].astype(int))))
    Ts = np.array(sorted(set(np.round(table["T"], 5))))
    ground = s["analytic_psi2_1"]
    rel = np.zeros((len(Ns), len(Ts)))
    for i, N in enumerate(Ns):
        for j, T in enumerate(Ts):
            row = table[(table["N"] == N) & np.isclose(table["T"], T)]
            rel[i, j] = 100.0 * (row["psi2_1"][0] - ground) / ground

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.3), gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.28})
    ax = axes[0]
    cmap = LinearSegmentedColormap.from_list("err", ["#f7f3ea", "#f0c9a0", "#c24d2c", "#6e2414"])
    image = ax.imshow(rel, origin="lower", cmap=cmap, aspect="auto", vmin=0.0, vmax=max(12.0, rel.max()))
    ax.set_xticks(range(len(Ts)), [f"{T:g}" for T in Ts])
    ax.set_yticks(range(len(Ns)), [str(int(N)) for N in Ns])
    ax.set_xlabel("imaginary time  T")
    ax.set_ylabel("number of slices  N")
    ax.set_title(r"percent error in $|\phi(1)|^2$", loc="left")
    for i in range(len(Ns)):
        for j in range(len(Ts)):
            value = rel[i, j]
            color = "#fbf7f0" if value > 6 else INK
            ax.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=7.4, color=color)
    mark_i = list(Ns).index(7)
    mark_j = int(np.where(np.isclose(Ts, 4.0))[0][0])
    ax.add_patch(Rectangle((mark_j - 0.5, mark_i - 0.5), 1, 1, fill=False, edgecolor=NAVY, lw=1.6))
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#d9d0c2")
    bar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    bar.outline.set_visible(False)
    bar.set_label("% above the true ground state")

    ax2 = axes[1]
    xs = np.linspace(0.0, 2.0, 400)
    shown = [(7, 4.0, TEAL, "N = 7, T = 4   milestone"), (7, 16.0, NAVY, "N = 7, T = 16"), (160, 4.0, GOLD, "N = 160, T = 4")]
    for N, T, color, label in shown:
        row = table[(table["N"] == N) & np.isclose(table["T"], T)]
        ax2.plot(xs, gauss(xs, row["alpha"][0]), color=color, lw=1.8, label=label)
    ax2.plot(xs, gauss(xs, 1.0), color=INK, lw=1.3, ls=(0, (5, 2)), label="true ground state")
    ax2.set_xlim(0.0, 2.0)
    ax2.set_ylim(0.0, 0.62)
    ax2.set_xlabel("x")
    ax2.set_ylabel(r"$|\phi(x)|^2$")
    ax2.set_title("what actually moves the curve", loc="left")
    ax2.legend(fontsize=8.6, loc="upper right")
    fig.suptitle("Coarse N is not the main error. Short T is.", fontsize=14, x=0.02, ha="left")
    fig.text(0.02, -0.02, "The box on the left is the milestone lattice. Waiting longer (bigger T) pulls the answer onto the ground state. Adding slices at T = 4 only reaches the exact answer for that short time.", color=MUTED, fontsize=9)
    save(fig, "convergence.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
