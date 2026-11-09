import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.collections import LineCollection

from style import COPPER, DATA, INK, MUTED, PAPER, TEAL, save


def main():
    world = np.genfromtxt(os.path.join(DATA, "metro_ho_world.csv"), delimiter=",", names=True)
    paths = np.genfromtxt(os.path.join(DATA, "metro_ho_paths.csv"), delimiter=",", names=True)
    hist = np.genfromtxt(os.path.join(DATA, "metro_ho_hist.csv"), delimiter=",", names=True)
    bridges = np.genfromtxt(os.path.join(DATA, "bridges.csv"), delimiter=",", names=True)

    tau = np.round(world["tau"], 8)
    xx = np.round(world["x"], 8)
    taus = np.unique(tau)
    xs = np.unique(xx)
    grid = np.zeros((taus.size, xs.size))
    grid[np.searchsorted(taus, tau), np.searchsorted(xs, xx)] = world["density"]
    keep = np.abs(xs) <= 3.15
    xs = xs[keep]
    grid = grid[:, keep]

    dark = "#141820"
    ink = "#f3efe6"
    cmap = LinearSegmentedColormap.from_list("heat", ["#141820", "#1c3f4a", "#0c6e6e", "#d7a441", "#f4e2b8"])

    fig = plt.figure(figsize=(11.2, 8.4), facecolor=dark)
    outer = fig.add_gridspec(2, 1, height_ratios=[1.55, 1.0], left=0.07, right=0.98, top=0.86, bottom=0.07, hspace=0.38)
    top = outer[0].subgridspec(1, 2, width_ratios=[3.3, 1.05], wspace=0.06)
    bot = outer[1].subgridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.28)
    heat = fig.add_subplot(top[0, 0], facecolor=dark)
    side = fig.add_subplot(top[0, 1], facecolor=dark, sharey=heat)
    polar = fig.add_subplot(bot[0, 0], projection="polar", facecolor=dark)
    bridge_ax = fig.add_subplot(bot[0, 1], facecolor=dark)

    heat.pcolormesh(taus, xs, grid.T, shading="auto", cmap=cmap, rasterized=True)
    ids = np.unique(paths["path"].astype(int))
    for pid in ids[:12]:
        sl = paths[paths["path"].astype(int) == pid]
        heat.plot(sl["tau"], sl["x"], color="#f6e7c1", lw=0.7, alpha=0.55)
    heat.set_ylim(-3.1, 3.1)
    heat.set_xlim(taus.min(), taus.max())
    heat.set_xlabel("imaginary time  τ", color=ink)
    heat.set_ylabel("x", color=ink)
    heat.tick_params(colors=ink)
    for spine in heat.spines.values():
        spine.set_color("#3a4150")
    heat.set_title("Where the Metropolis paths actually sit", color=ink, loc="left")
    heat.axhline(1.0, color="#c24d2c", lw=0.6, ls=(0, (2, 2)), alpha=0.8)
    heat.axhline(-1.0, color="#c24d2c", lw=0.6, ls=(0, (2, 2)), alpha=0.8)

    side.fill_betweenx(hist["x"], 0.0, hist["psi2"], color="#1f8a84", alpha=0.85, linewidth=0)
    side.plot(hist["psi2"], hist["x"], color="#f4e2b8", lw=1.0)
    side.set_ylim(-3.1, 3.1)
    side.set_xlim(0.0, 0.72)
    side.set_xlabel(r"$|\phi|^2$", color=ink)
    side.tick_params(colors=ink)
    side.set_yticklabels([])
    for spine in side.spines.values():
        spine.set_color("#3a4150")
    side.set_title("marginal", color=ink, loc="left", fontsize=11)

    T = paths["tau"].max()
    for pid in ids[:8]:
        sl = paths[paths["path"].astype(int) == pid]
        theta = 2.0 * np.pi * sl["tau"] / T
        radius = 3.4 + sl["x"]
        points = np.column_stack([theta, radius]).reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        colors = plt.cm.copper(np.linspace(0.15, 0.9, len(sl) - 1))
        polar.add_collection(LineCollection(segments, colors=colors, linewidths=0.9, alpha=0.85))
    polar.set_ylim(0.0, 6.2)
    polar.set_yticklabels([])
    polar.set_xticklabels([])
    polar.spines["polar"].set_color("#3a4150")
    polar.grid(color="#2a3140", lw=0.5)
    polar.set_title("same paths, time as an angle", color=ink, fontsize=11, pad=14)

    bids = np.unique(bridges["path"].astype(int))
    weights = []
    for pid in bids:
        w = bridges[bridges["path"].astype(int) == pid]["weight"][0]
        weights.append(w if w > 0 else np.nan)
    weights = np.array(weights)
    logs = np.log10(weights)
    finite = np.isfinite(logs)
    norm = Normalize(vmin=np.nanpercentile(logs, 8), vmax=np.nanpercentile(logs, 92))
    cmap_b = LinearSegmentedColormap.from_list("w", ["#3d4d66", "#d7a441", "#fff6df"])
    for pid, logw in zip(bids, logs):
        sl = bridges[bridges["path"].astype(int) == pid]
        color = "#3d4d66" if not np.isfinite(logw) else cmap_b(norm(logw))
        bridge_ax.plot(sl["tau"], sl["x"], color=color, lw=0.9, alpha=0.9)
    bridge_ax.set_ylim(-1.5, 3.4)
    bridge_ax.set_xlabel("τ", color=ink)
    bridge_ax.set_ylabel("x", color=ink)
    bridge_ax.tick_params(colors=ink)
    for spine in bridge_ax.spines.values():
        spine.set_color("#3a4150")
    bridge_ax.set_title("milestone bridges pinned at x = 1, coloured by weight", color=ink, fontsize=11, loc="left")
    fig.suptitle("Imaginary-time paths for the harmonic oscillator", color=ink, fontsize=15, x=0.06, ha="left")
    fig.text(0.06, 0.915, "Top: Metropolis, N = 64, T = 12. Dashed lines are the classical turning points. Bottom right: the direct N = 7 paths. Pale lines carry almost all the weight.", color="#b7b1a6", fontsize=9)
    save(fig, "worldlines.png")
    plt.close(fig)

    fig2, ax = plt.subplots(figsize=(8.2, 4.6), facecolor=PAPER)
    order = np.argsort(weights)
    shown = order[:: max(1, len(order) // 14)]
    for rank, pid in enumerate(shown):
        sl = bridges[bridges["path"].astype(int) == int(bids[pid])]
        logw = logs[pid]
        color = cmap_b(norm(logw)) if np.isfinite(logw) else "#8d877c"
        ax.plot(sl["tau"], sl["x"], color=color, lw=1.35, alpha=0.95)
    ax.scatter([0, sl["tau"].max()], [1, 1], color=COPPER, s=28, zorder=5)
    ax.set_xlabel("imaginary time  τ")
    ax.set_ylabel("x")
    ax.set_title("A handful of the random paths used in the milestone sum")
    ax.set_ylim(-2.2, 4.2)
    fig2.text(0.125, -0.02, "Every path starts and ends at x = 1. Colour is log(weight). The pale ones dominate the average.", color=MUTED, fontsize=9)
    save(fig2, "bridges.png")
    plt.close(fig2)


if __name__ == "__main__":
    main()
