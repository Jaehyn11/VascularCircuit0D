import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------- user settings ----------------
sim_file = "../output.csv"                  # path to your model output
val_file = "Validation_Pressure_AAo.csv"    # Stergiopulos data
T_cycle = 1.0                               # length of one cardiac cycle [s]
# ------------------------------------------------

# 1. Load model output (has header)
sim = pd.read_csv(sim_file)

# column names from your solver
SIM_T_COL = "t"
SIM_P_COL = "P_0"

t_sim = sim[SIM_T_COL].to_numpy()
P_sim = sim[SIM_P_COL].to_numpy()

# restrict to first cycle
mask_sim = (t_sim >= 0.0) & (t_sim <= T_cycle)
t_sim = t_sim[mask_sim]
P_sim = P_sim[mask_sim]

# 2. Load validation data (NO header; columns: t, P_val)
val = pd.read_csv(val_file, header=None, names=["t_val", "P_val"])
t_val = val["t_val"].to_numpy()
P_val = val["P_val"].to_numpy()

mask_val = (t_val >= 0.0) & (t_val <= T_cycle)
t_val = t_val[mask_val]
P_val = P_val[mask_val]

# 3. Make comparison plot (academic style)
plt.figure(figsize=(5.0, 3.5))

plt.plot(t_val, P_val, label="Stergiopulos (1994)", linewidth=2)
plt.plot(t_sim, P_sim, label="0D RLC model", linestyle="--", linewidth=2)

plt.xlabel(r"$t$ [s]")
plt.ylabel(r"$P_{\mathrm{AAo}}$ [mmHg]")
plt.xlim(0.0, T_cycle)
plt.grid(alpha=0.3)
plt.legend(frameon=False)
plt.tight_layout()

plt.savefig("P_AAo_model_vs_data.png", dpi=600)
