#!/usr/bin/env python3
"""
Compute RMSE, NRMSE, MAE for every (case, run, variable), average over
stations, and draw heat-maps.  Designed for the folder+filename scheme:

  csv_files/case_20231102_mp4_pbl1_station_17061.csv
  obs_csvs/202412_precip_mm_17061.csv
  …
"""

import os
import sys
from glob import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ────────────────────────────────────────────────────────────────────
# 0.  CONFIGURATION
# ────────────────────────────────────────────────────────────────────
CASES = ["20231102", "20231129", "20240120", "20241123"]

RUNS = [(4, 1), (4, 5), (4, 7),
        (6, 1), (6, 5), (6, 7),
        (38, 1), (38, 5), (38, 7),
        (10, 1), (10, 5), (10, 7),
        (24, 1), (24, 5), (24, 7)]
RUN_LABELS = [f"MP={m}, PBL={p}" for m, p in RUNS]

# model-column  → canonical variable id
VAR_MAP_MODEL = {
    "Yagis(mm)":              "precip_mm",
    "2_Metre_Sicaklik(C)":    "t2_degC",
    "10_Metre_Ruzgar(m/s)":   "ws10_mps",
}

# filename-token → canonical variable id
VAR_MAP_OBS = {
    "precip_mm":              "precip_mm",
    "t2_degC":                "t2_degC",
    "wdir_ws_deg_mps":        "ws10_mps",
    # extra tokens (unused in current task, but harmless):
    "mslp_hPa":               "mslp_hPa",
    "td2_degC":               "td2_degC",
}

METRICS = ["RMSE", "NRMSE", "MAE"]


def nrmse(pred, obs, eps=1e-12):
    """NRMSE with zero-range protection.

    • If obs.max() == obs.min() (constant series) return NaN.
    • Otherwise RMSE divided by (max-min).
    """
    rmse = np.sqrt(((pred - obs) ** 2).mean())
    rng  = obs.max() - obs.min()
    return np.nan if abs(rng) < eps else rmse / rng

# ────────────────────────────────────────────────────────────────────
# 1.  FILE DISCOVERY
# ────────────────────────────────────────────────────────────────────
model_files = glob(os.path.join("csv_files", "case_*_mp*_pbl*_station_*.csv"))
if not model_files:
    sys.exit("No model CSVs found under csv_files/ – aborting.")

obs_files = glob(os.path.join("obs_csvs", "*.csv"))
if not obs_files:
    sys.exit("No observation CSVs found under obs_csvs/ – aborting.")

# Build {(canonical_var, station_id, yyyymm): path}
obs_lookup = {}
for f in obs_files:
    base = os.path.basename(f)
    if len(base) < 11 or not base[:6].isdigit():
        continue
    yyyymm = base[:6]
    station_id = base.split("_")[-1].split(".")[0]

    canonical_var = None
    for token, canon in VAR_MAP_OBS.items():
        if token in base:
            canonical_var = canon
            break

    if canonical_var:
        obs_lookup[(canonical_var, station_id, yyyymm)] = f

if not obs_lookup:
    sys.exit("Observation lookup table is empty – check filenames.")

# ────────────────────────────────────────────────────────────────────
# 2.  HELPERS
# ────────────────────────────────────────────────────────────────────
def parse_model_filename(path):
    """Return case, mp, pbl, station_id"""
    base = os.path.basename(path)
    parts = base.split("_")
    case  = parts[1]          # 20231102
    mp    = int(parts[2][2:]) # mp4 → 4
    pbl   = int(parts[3][3:]) # pbl1 → 1
    stid  = parts[5].split(".")[0]
    return case, mp, pbl, stid


def metric_values(pred, obs):
    rmse = np.sqrt(((pred - obs) ** 2).mean())
    mae  = np.abs(pred - obs).mean()
    return rmse, nrmse(pred, obs), mae


# ────────────────────────────────────────────────────────────────────
# 3.  MAIN LOOP – COLLECT ERRORS
# ────────────────────────────────────────────────────────────────────
arrays = [[m for m in METRICS for _ in CASES],
          [c for _ in METRICS for c in CASES]]
results = pd.DataFrame(index=RUN_LABELS,
                       columns=pd.MultiIndex.from_arrays(arrays,
                                                         names=["Metric", "Case"]),
                       dtype=float)

missing_obs = 0
empty_merge = 0

for mp_run in model_files:
    case, mp, pbl, stid = parse_model_filename(mp_run)
    run_label = f"MP={mp}, PBL={pbl}"

    df_mod = (pd.read_csv(mp_run, parse_dates=["Tarih"])
                .rename(columns={"Tarih": "datetime"}))

    for model_col, canonical_var in VAR_MAP_MODEL.items():
        if model_col not in df_mod.columns:
            continue

        yyyymm  = case[:6]
        obs_path = obs_lookup.get((canonical_var, stid, yyyymm))
        if not obs_path:
            missing_obs += 1
            continue

        df_obs = pd.read_csv(obs_path, parse_dates=["datetime"])
        df_obs = df_obs.rename(columns={"value": canonical_var})

        df = pd.merge(df_mod[["datetime", model_col]],
                      df_obs[["datetime", canonical_var]],
                      on="datetime", how="inner").dropna()

        if df.empty:
            empty_merge += 1
            continue

        pred, obs = df[model_col], df[canonical_var]
        rmse, nrmse_val, mae = metric_values(pred, obs)

        for mname, v in zip(METRICS, (rmse, nrmse_val, mae)):
            if np.isfinite(v):                                 # v is a real number
                current = results.at[run_label, (mname, case)]
                results.at[run_label, (mname, case)] = (
                    v if np.isnan(current) else np.nanmean([current, v])
    )

print(f"[INFO] done → {missing_obs} obs misses, {empty_merge} empty merges\n")

# ────────────────────────────────────────────────────────────────────
# 4.  PLOT HEAT-MAPS  (Fix 1)
# ────────────────────────────────────────────────────────────────────
valid_metrics = []
panels        = {}

for metric in METRICS:
    panel = results.xs(metric, axis=1, level="Metric")
    if panel.isna().all().all():
        print(f"[WARN] {metric} heat-map skipped – all NaN")
        continue
    valid_metrics.append(metric)
    panels[metric] = panel

if not valid_metrics:
    sys.exit("Nothing to plot – every metric panel was empty.")

plt.figure(figsize=(5 * len(valid_metrics), 8))

for i, metric in enumerate(valid_metrics, 1):
    ax = plt.subplot(1, len(valid_metrics), i)
    #cmap = {"RMSE": "rocket_r", "NRMSE": "crest_r", "MAE": "mako_r"}[metric]
    cmap = {"RMSE": "rocket_r", "NRMSE": "rocket_r", "MAE": "rocket_r"}[metric]
    sns.heatmap(panels[metric],
                annot=True, fmt=".2f", linewidths=.3,
                cmap=cmap, vmin=0, cbar=False, ax=ax)
    ax.set_title(metric)
    ax.set_xlabel("Case")

plt.tight_layout()
plt.savefig("error_heatmaps.png", dpi=300)
print("Saved → error_heatmaps.png")
