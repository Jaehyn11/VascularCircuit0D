import matplotlib
matplotlib.use("Agg")  # non-interactive backend (ok for saving animations)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, FFMpegWriter

# ----------------- user settings -----------------
sim_file = "../output.csv"  # path to solver output
T_cycle  = 1.0              # one cardiac cycle [s]
out_mp4  = "P_1_2_3_anim.mp4"
stride   = 10                # use every 5th time step
fps      = 100               # frames per second in final video
# -------------------------------------------------

# Font / style
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"]  = ["Times", "DejaVu Serif", "STIXGeneral"]
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["xtick.color"] = "black"
plt.rcParams["ytick.color"] = "black"
plt.rcParams["text.color"]  = "black"
plt.rcParams["axes.labelcolor"] = "black"

# 1. Load data
sim = pd.read_csv(sim_file)

t_all  = sim["t"].to_numpy()
P1_all = sim["P_1"].to_numpy()   # Thoracic Aorta
P2_all = sim["P_2"].to_numpy()   # Iliac Artery
P3_all = sim["P_3"].to_numpy()   # Common Carotid

# Keep only first cycle [0, T_cycle] and thin by stride
mask = (t_all >= 0.0) & (t_all <= T_cycle)
t  = t_all[mask][::stride]
P1 = P1_all[mask][::stride]
P2 = P2_all[mask][::stride]
P3 = P3_all[mask][::stride]

pressures = [P1, P2, P3]
titles = [
    "Thoracic Aorta (P$_1$)",
    "Iliac Artery (P$_2$)",
    "Common Carotid (P$_3$)"
]

# y-limits shared across panels
P_all = np.concatenate(pressures)
P_min = P_all.min()
P_max = P_all.max()
margin = 0.05 * (P_max - P_min + 1e-6)
ymin = P_min - margin
ymax = P_max + margin

# 2. Figure & empty lines
fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.5), sharex=True, sharey=True)
fig.patch.set_facecolor("white")
for ax in axes:
    ax.set_facecolor("white")

lines = []
markers = []
for ax, P, title in zip(axes, pressures, titles):
    line, = ax.plot([], [], color="black", linewidth=1.5)
    marker, = ax.plot([], [], "o", color="black", markersize=3)
    lines.append(line)
    markers.append(marker)

    ax.set_xlim(0.0, T_cycle)
    ax.set_ylim(ymin, ymax)
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3, linestyle=":")

axes[0].set_ylabel("Pressure [mmHg]")
for ax in axes:
    ax.set_xlabel("Time [s]")

# 3. Animation functions
n_frames = len(t)

def init():
    for line, marker in zip(lines, markers):
        line.set_data([], [])
        marker.set_data([], [])
    return lines + markers

def update(frame):
    t_slice = t[: frame + 1]
    for P, line, marker in zip(pressures, lines, markers):
        P_slice = P[: frame + 1]
        line.set_data(t_slice, P_slice)
        marker.set_data([t[frame]], [P[frame]])
    return lines + markers

interval_ms = 1000.0 / fps

anim = FuncAnimation(
    fig,
    update,
    init_func=init,
    frames=n_frames,
    interval=interval_ms,
    blit=True,
)

# 4. Save animation (ffmpeg)
writer = FFMpegWriter(
    fps=fps,
    bitrate=8000,  # higher than default -> cleaner video
)
anim.save(out_mp4, writer=writer, dpi=400)

print(f"Saved animation to {out_mp4}")
