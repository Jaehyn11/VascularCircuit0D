import matplotlib
matplotlib.use("Agg")  # non-interactive backend (for cluster use)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------- user settings ----------
sim_file   = "../output.csv"
val_file   = "Validation_Pressure_AAo.csv"
inflow_file = "AAo4.csv"
T_cycle    = 1.0       # length of one cardiac cycle [s]
out_png    = "P_AAo_with_flow_validation.png"
# -----------------------------------

# ---- Font & style setup (no Times New Roman string!) ----
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times", "DejaVu Serif", "STIXGeneral"]  # Times-like
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["xtick.color"] = "black"
plt.rcParams["ytick.color"] = "black"
plt.rcParams["text.color"] = "black"
plt.rcParams["axes.labelcolor"] = "black"
plt.rcParams["legend.edgecolor"] = "black"


# 1. Load model output (has header: t, P_0, ...)
sim = pd.read_csv(sim_file)
t_sim = sim["t"].to_numpy()
P_sim = sim["P_0"].to_numpy()

mask_sim = (t_sim >= 0.0) & (t_sim <= T_cycle)
t_sim = t_sim[mask_sim]
P_sim = P_sim[mask_sim]
# t_sim = t_sim[::15]
# P_sim = P_sim[::15]

# 2. Load validation pressure (no header: first col t, second col P_val)
val = pd.read_csv(val_file, header=None, names=["t_val", "P_val"])
t_val = val["t_val"].to_numpy()
P_val = val["P_val"].to_numpy()
mask_val = (t_val >= 0.0) & (t_val <= T_cycle)
t_val = t_val[mask_val]
P_val = P_val[mask_val]
t_val = t_val[::3]
P_val = P_val[::3]

# 3. Load inflow waveform (no header: first col t, second col Q_in [mL/s])
qin = pd.read_csv(inflow_file, header=None, names=["t_in", "Q_in"])
t_in = qin["t_in"].to_numpy()
Q_in = qin["Q_in"].to_numpy()
mask_in = (t_in >= 0.0) & (t_in <= T_cycle)
t_in = t_in[mask_in]
Q_in = Q_in[mask_in]

# 4. Plot: pressure on left axis, flow on right axis
fig, axP = plt.subplots(figsize=(5.0, 3.5))
axQ = axP.twinx()

# Validation pressure: dashed
axP.plot(
    t_val, P_val,
    linestyle="None",
    marker="s",
    markersize=2.5,
    linewidth=1.0,
    color="black",
    label="Stergiopulos et al.",
)

# Model pressure: solid with star markers
axP.plot(
    t_sim, P_sim,
    linestyle="-",
    marker="None",
    markersize=1.5,
    linewidth=2.0,
    color="black",
    label="VascularCircuit0D",
)

# Inflow: solid line on secondary axis
axQ.plot(
    t_in, Q_in,
    linestyle="-",
    linewidth=1.5,
    color="red",
    label="AAo Inflow",
)

# Axes labels / limits
axP.set_xlabel("Time [s]")
axP.set_ylabel("Pressure [mmHg]")
axQ.set_ylabel("Flow [mL/s]")
axP.set_xlim(0.0, T_cycle)

axP.grid(alpha=0.0, color="black", linestyle=":")

# Combined legend (both axes)
linesP, labelsP = axP.get_legend_handles_labels()
linesQ, labelsQ = axQ.get_legend_handles_labels()
axP.legend(
    linesP + linesQ,
    labelsP + labelsQ,
    frameon=False,
    loc="upper right",
    fontsize=9,
)

plt.tight_layout()
plt.savefig(out_png, dpi=600)

