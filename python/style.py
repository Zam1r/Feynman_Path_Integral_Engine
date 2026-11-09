import os

import matplotlib as mpl

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA = os.path.join(ROOT, "data")
FIG = os.path.join(ROOT, "figures")

PAPER = "#f3efe6"
PANEL = "#fbf7f0"
INK = "#1c1915"
MUTED = "#5e584f"
TEAL = "#0c6e6e"
TEAL_SOFT = "#c5ddd9"
COPPER = "#c24d2c"
GOLD = "#a97820"
NAVY = "#243656"
GREEN = "#2c6b45"
PLUM = "#6a3a58"

mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Palatino Linotype", "Georgia", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "text.color": INK,
    "figure.facecolor": PAPER,
    "axes.facecolor": PANEL,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
    "savefig.dpi": 170,
    "pdf.fonttype": 42,
})


def ensure_fig():
    os.makedirs(FIG, exist_ok=True)


def save(fig, name):
    ensure_fig()
    path = os.path.join(FIG, name)
    fig.savefig(path, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(path)
    return path


def summary():
    out = {}
    with open(os.path.join(DATA, "summary.csv"), encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            key, value = line.strip().split(",")
            out[key] = float(value)
    return out
