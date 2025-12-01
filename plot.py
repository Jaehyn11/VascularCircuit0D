import matplotlib

matplotlib.use("Agg")  # non-interactive backend
import math
import os

import matplotlib.pyplot as plt
import numpy as np

# ---------------- User settings ----------------
output_csv = "output.csv"  # solver output: t, P_0.., Q_0..
# Build inflow path from $HOME
home = os.path.expanduser("~")
base = os.path.join(home, "NERS570", "VascularCircuit0D")  # setme
inflow_csv = os.path.join(base, "input", "AAo.csv")
inflow_csv = "./input/AAo.csv"  # override for local testing
# ------------------------------------------------


def latex_name(varname):
    """Converts column name like 'P_0' → r'$P_{\\mathrm{0}}$'"""
    if "_" in varname:
        base, sub = varname.split("_", 1)
        return rf"${base}_{{\mathrm{{{sub}}}}}$"
    return rf"${varname}$"


# 1. Load solver output
data = np.genfromtxt(output_csv, delimiter=",", names=True)
t = data["t"]

# Separate pressure and flow columns by prefix
all_names = data.dtype.names
pressure_names = [name for name in all_names if name.startswith("P_")]
flow_names = [name for name in all_names if name.startswith("Q_")]

# 2. Load inflow CSV: time [s], Q [mL/s]
inflow_data = np.genfromtxt(inflow_csv, delimiter=",", names=True)
t_in_raw = inflow_data[inflow_data.dtype.names[0]]
Q_in_raw = inflow_data[inflow_data.dtype.names[1]]

# 3. Build periodic inflow waveform on solver's time grid
t0 = t_in_raw[0]
t1 = t_in_raw[-1]
T_cycle = t1 - t0

t_mod = np.mod(t - t0, T_cycle) + t0
Qin = np.interp(t_mod, t_in_raw, Q_in_raw)

# ---------------- Figure 1: Inflow + Pressures ----------------
n_pressures = len(pressure_names)
n_total_panels_P = 1 + n_pressures  # 1 inflow + N pressures
ncols = 2
nrows_P = math.ceil(n_total_panels_P / ncols)

figP, axsP = plt.subplots(nrows_P, ncols, figsize=(12, 3.2 * nrows_P), sharex=True)
axsP = np.atleast_1d(axsP).flatten()

# Panel 0: Inflow
axsP[0].plot(t, Qin)
axsP[0].set_ylabel(r"$Q_{\mathrm{in}}$ [mL/s]")
axsP[0].set_title("Inflow and Pressures")
axsP[0].grid(True, alpha=0.3)

# Panels for pressures
for i, name in enumerate(pressure_names, start=1):
    ax = axsP[i]
    ax.plot(t, data[name])
    ax.set_ylabel(latex_name(name) + r" [mmHg]")
    ax.grid(True, alpha=0.3)

# X-labels
for ax in axsP:
    ax.set_xlabel(r"$t$ [s]")

# Hide unused axes if any
for j in range(n_total_panels_P, len(axsP)):
    figP.delaxes(axsP[j])

plt.tight_layout()
plt.savefig("pressures.png", dpi=600)

# ---------------- Figure 2: Inflow + Flows ----------------
n_flows = len(flow_names)
n_total_panels_Q = 1 + n_flows  # 1 inflow + N flows
nrows_Q = math.ceil(n_total_panels_Q / ncols)

figQ, axsQ = plt.subplots(nrows_Q, ncols, figsize=(12, 3.2 * nrows_Q), sharex=True)
axsQ = np.atleast_1d(axsQ).flatten()

# Panel 0: Inflow again for reference
axsQ[0].plot(t, Qin)
axsQ[0].set_ylabel(r"$Q_{\mathrm{in}}$ [mL/s]")
axsQ[0].set_title("Inflow and Segment Flows")
axsQ[0].grid(True, alpha=0.3)

# Panels for flows
for i, name in enumerate(flow_names, start=1):
    ax = axsQ[i]
    ax.plot(t, data[name])
    ax.set_ylabel(latex_name(name) + r" [mL/s]")
    ax.grid(True, alpha=0.3)

# X-labels
for ax in axsQ:
    ax.set_xlabel(r"$t$ [s]")

# Hide unused axes if any
for j in range(n_total_panels_Q, len(axsQ)):
    figQ.delaxes(axsQ[j])

plt.tight_layout()
plt.savefig("flows.png", dpi=600)
