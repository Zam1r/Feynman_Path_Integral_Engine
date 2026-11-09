import os

import matplotlib.pyplot as plt
import numpy as np

from style import COPPER, DATA, GOLD, INK, MUTED, NAVY, TEAL, TEAL_SOFT, save, summary


def acf(series, lags):
    x = np.asarray(series, dtype=float)
    x = x - x.mean()
    var = np.dot(x, x) / len(x)
    out = np.empty(lags)
    for lag in range(lags):
        if lag == 0:
            out[lag] = 1.0
        else:
            out[lag] = np.dot(x[:-lag], x[lag:]) / (len(x) - lag) / var
    return out


def main():
    s = summary()
    hist = np.genfromtxt(os.path.join(DATA, "metro_ho_hist.csv"), delimiter=",", names=True)
    trace = np.genfromtxt(os.path.join(DATA, "metro_ho_trace.csv"), delimiter=",", names=True)
    corr = np.genfromtxt(os.path.join(DATA, "metro_ho_corr.csv"), delimiter=",", names=True)
    paths = np.genfromtxt(os.path.join(DATA, "metro_ho_paths.csv"), delimiter=",", names=True)

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.4))
    fig.subplots_adjust(hspace=0.38, wspace=0.28)

    ax = axes[0, 0]
    meas = trace["measuring"] > 0.5
    ax.plot(trace["sweep"][~meas], trace["x2"][~meas], color=GOLD, lw=0.9, label="still settling")
    ax.plot(trace["sweep"][meas], trace["x2"][meas], color=TEAL, lw=0.9, label="used for the average")
    ax.axhline(0.5, color=INK, lw=1.0, ls=(0, (4, 2)), label=r"true $\langle x^2\rangle = 1/2$")
    ax.set_xlabel("sweep")
    ax.set_ylabel(r"$\langle x^2\rangle$ on the path")
    ax.set_title("Does the chain settle?", loc="left")
    ax.legend(fontsize=8)
    ax.set_ylim(0.0, max(1.2, np.percentile(trace["x2"], 99)))

    ax = axes[0, 1]
    live = trace["x2"][meas]
    lags = acf(live, 80)
    ax.plot(np.arange(len(lags)), lags, color=NAVY, lw=1.6)
    ax.axhline(0.0, color=MUTED, lw=0.6)
    ax.set_xlabel("lag, in saved sweeps")
    ax.set_ylabel("autocorrelation of  " + r"$\langle x^2\rangle$")
    ax.set_title("How independent are the samples?", loc="left")
    ax.set_ylim(-0.15, 1.05)

    ax = axes[1, 0]
    half = corr["tau"] <= 6.0
    ax.plot(corr["tau"][half], corr["theory"][half], color=COPPER, lw=1.8, label="exact continuum correlator")
    ax.fill_between(
        corr["tau"][half],
        corr["xx"][half] - corr["err"][half],
        corr["xx"][half] + corr["err"][half],
        color=TEAL_SOFT, linewidth=0,
    )
    ax.plot(corr["tau"][half], corr["xx"][half], color=TEAL, lw=1.3, label="measured on the lattice")
    ax.set_xlabel(r"separation  $\tau$")
    ax.set_ylabel(r"$\langle x(\tau)\,x(0)\rangle$")
    ax.set_title(f"Excited-state gap comes out near {s['corr_gap']:.2f}", loc="left")
    ax.legend(fontsize=8)
    ax.set_xlim(0.0, 6.0)

    ax = axes[1, 1]
    ids = np.unique(paths["path"].astype(int))
    left, right = [], []
    for pid in ids:
        sl = paths[paths["path"].astype(int) == pid]
        left.append(sl["x"][:-1])
        right.append(sl["x"][1:])
    left = np.concatenate(left)
    right = np.concatenate(right)
    ax.hist2d(left, right, bins=70, range=[[-2.6, 2.6], [-2.6, 2.6]], cmap="YlGnBu", cmin=1)
    ax.plot([-2.6, 2.6], [-2.6, 2.6], color=COPPER, lw=0.8, ls=(0, (3, 2)))
    ax.set_xlabel(r"$x_j$")
    ax.set_ylabel(r"$x_{j+1}$")
    ax.set_title("Neighbours stay close", loc="left")
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#d9d0c2")

    fig.suptitle("Health of the Metropolis run", fontsize=15, x=0.02, ha="left")
    fig.text(
        0.02, -0.01,
        f"Acceptance during the measurement was {100 * s['metro_accept']:.0f}%.  "
        f"Mean square position {s['metro_x2']:.3f} ± {s['metro_x2_err']:.3f}, against 0.500.  "
        f"The slope of the correlator gives a gap of {s['corr_gap']:.2f}; the oscillator gap is 1.",
        color=MUTED, fontsize=9,
    )
    save(fig, "diagnostics.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
