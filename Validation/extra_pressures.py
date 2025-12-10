import matplotlib
matplotlib.use("Agg")  # for cluster / non-GUI use

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ----------------- user settings -----------------
sim_file   = "../output.csv"  # path to solver output
T_cycle    = 1.0              # cardiac period [s]
n_cycles   = 5                # how many cycles to show
fold_cycles = True            # if True: plot each cycle on 0–1 s
out_png    = "P_1_2_3_subplots.png"
# -------------------------------------------------

# Font / style (Times-like, black)
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"]  = ["Times", "DejaVu Serif", "STIXGeneral"]
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["xtick.color"] = "black"
plt.rcParams["ytick.color"] = "black"
plt.rcParams["text.color"]  = "black"
plt.rcParams["axes.labelcolor"] = "black"

# 1. Load output.csv
sim = pd.read_csv(sim_file)

t   = sim["t"].to_numpy()
P1  = sim["P_1"].to_numpy()   # Thoracic Aorta
P2  = sim["P_2"].to_numpy()   # Iliac Artery
P3  = sim["P_3"].to_numpy()   # Common Carotid (mapped to node 3 in your chain)

pressures = [P1, P2, P3]
labels    = [
    "Thoracic Aorta (P$_1$)",
    "Iliac Artery (P$_2$)",
    "Common Carotid (P$_3$)",
]

# Restrict in time for n_cycles
t_max = n_cycles * T_cycle
mask  = (t >= 0.0) & (t <= t_max)
t     = t[mask]
pressures = [p[mask] for p in pressures]

# Global y-limits across P1, P2, P3
P_all = np.concatenate(pressures)
P_min = np.min(P_all)
P_max = np.max(P_all)
margin = 0.05 * (P_max - P_min + 1e-6)
ymin = P_min - margin
ymax = P_max + margin

# 2. Figure with 3 columns
fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.0), sharex=True, sharey=True)

for ax, P, lab in zip(axes, pressures, labels):

    if fold_cycles:
        # Plot each cycle mapped to x in [0, T_cycle]
        n_avail_cycles = int(np.floor(t[-1] / T_cycle))
        for k in range(n_avail_cycles):
            t_start = k * T_cycle
            t_end   = (k + 1) * T_cycle
            cyc_mask = (t >= t_start) & (t <= t_end)
            if not np.any(cyc_mask):
                continue
            t_cyc = t[cyc_mask] - t_start  # fold into [0, T_cycle]
            P_cyc = P[cyc_mask]
            ax.plot(
                t_cyc, P_cyc,
                color="black",
                linewidth=1.2,
                alpha=0.7
            )
        ax.set_xlim(0.0, T_cycle)
    else:
        # Plot continuous time series 0..n_cycles*T_cycle
        ax.plot(
            t, P,
            color="black",
            linewidth=1.2
        )
        ax.set_xlim(0.0, t_max)

    ax.set_ylim(ymin, ymax)
    ax.set_title(lab, fontsize=10)
    ax.grid(alpha=0.3, linestyle=":")

# Axis labels
axes[0].set_ylabel("Pressure [mmHg]")
for ax in axes:
    ax.set_xlabel("Time [s]")

plt.tight_layout()
plt.savefig(out_png, dpi=600)
