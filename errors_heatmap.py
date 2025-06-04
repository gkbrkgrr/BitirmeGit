#!/usr/bin/env python3
# errors_heatmap_fixed.py
"""Create per-variable error heat-maps (RMSE, NRMSE, MAE x 4 cases) averaged
over 13 WMO stations for 15 WRF (MP,PBL) runs."""

import os, sys, itertools
from glob import glob
from collections import defaultdict

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# ─────────────── SETTINGS ─────────────── #
CASES = ["20231102", "20231129", "20240120", "20241123"]
RUNS  = [(4,1),(4,5),(4,7),(6,1),(6,5),(6,7),
         (28,1),(28,5),(28,7),(10,1),(10,5),(10,7),
         (24,1),(24,5),(24,7)]
RUN_LABELS  = [f"MP={m}, PBL={p}" for m,p in RUNS]
METRICS     = ["RMSE","NRMSE","MAE"]

VAR_MAP_MODEL = {"Yagis(mm)"            :"precip_mm",
                 "2_Metre_Sicaklik(C)"  :"t2_degC",
                 "10_Metre_Ruzgar(m/s)" :"ws10_mps"}
VAR_MAP_OBS   = {"precip_mm"            :"precip_mm",
                 "t2_degC"              :"t2_degC",
                 "wdir_ws_deg_mps"      :"ws10_mps"}
VAR_TITLE_VAR = {"precip_mm"            :"Hourly Precipitation",
                 "t2_degC"              :"Temperature",
                 "ws10_mps"      :"Wind Speed"}

# ─────────────── HELPERS ─────────────── #
def nrmse(pred, obs, eps=1e-12):
    rmse  = np.sqrt(((pred-obs)**2).mean())
    max_ = obs.max()
    return np.nan if abs(max_)<eps else rmse/max_

def metric_values(pred, obs):
    rmse = np.sqrt(((pred-obs)**2).mean())
    return {"RMSE":rmse, "NRMSE":nrmse(pred,obs), "MAE":np.abs(pred-obs).mean()}

def parse_model_filename(p):
    base = os.path.basename(p).split("_")
    return base[1], int(base[2][2:]), int(base[3][3:]), base[5].split(".")[0]  # case, mp, pbl, stid

# ─────────────── DISCOVER FILES ─────────────── #
model_files = glob("csv_files/case_*_mp*_pbl*_station_*.csv")
obs_files   = glob("obs_csvs/*.csv")
if not (model_files and obs_files):
    sys.exit("Model or observation CSVs missing.")

obs_lookup = {}
for f in obs_files:
    base = os.path.basename(f)
    if not (base[:6].isdigit() and "_" in base):
        continue
    yyyymm        = base[:6]
    station_id    = base.split("_")[-1].split(".")[0]
    canonical_var = next((v for t,v in VAR_MAP_OBS.items() if t in base), None)
    if canonical_var:
        obs_lookup[(canonical_var, station_id, yyyymm)] = f

# ─────────────── COLLECT RAW SCORES ─────────────── #
scores = defaultdict(  # var →
           lambda: defaultdict(          # run →
             lambda: defaultdict(        # metric →
               lambda: defaultdict(list) # case → list[float]
           )))

for mp_csv in model_files:
    case, mp, pbl, stid = parse_model_filename(mp_csv)
    run_lbl = f"MP={mp}, PBL={pbl}"

    mod = pd.read_csv(mp_csv, parse_dates=["Tarih"]).rename(columns={"Tarih":"datetime"})
    for col, var in VAR_MAP_MODEL.items():
        if col not in mod.columns: continue
        yyyymm   = case[:6]
        obs_path = obs_lookup.get((var, stid, yyyymm))
        if not obs_path: continue

        obs = (pd.read_csv(obs_path, parse_dates=["datetime"])
                 .rename(columns={"value":var}))

        paired = (mod[["datetime",col]]
                    .merge(obs[["datetime",var]], on="datetime", how="inner")
                    .dropna())
        if paired.empty: continue

        mets = metric_values(paired[col], paired[var])
        for m,v in mets.items():
            if np.isfinite(v):
                scores[var][run_lbl][m][case].append(v)

# ─────────────── BUILD & PLOT HEAT-MAPS ─────────────── #
for var, run_dict in scores.items():
    arrays = [[m for m in METRICS for _ in CASES],
              list(itertools.chain.from_iterable([CASES]*len(METRICS)))]
    cols   = pd.MultiIndex.from_arrays(arrays, names=["Metric", "Case"])
    results = pd.DataFrame(index=RUN_LABELS, columns=cols, dtype=float)

    for run_lbl, m_dict in run_dict.items():
        for m, case_dict in m_dict.items():
            for case, vals in case_dict.items():
                if vals:
                    results.loc[run_lbl, (m, case)] = np.mean(vals)

    # 1⃣  use *either* constrained-layout *or* manual spacing → choose manual
    fig, axes = plt.subplots(
        1, 3,
        figsize=(18, 10),
        sharey=False,
        gridspec_kw={"wspace": 0.12},
    )

    for i, m in enumerate(METRICS):
        metric_df = results.xs(m, level="Metric", axis=1)[CASES]
        cmap = {"RMSE": "rocket_r", "NRMSE": "crest_r", "MAE": "mako_r"}[m]
        sns.heatmap(
            metric_df,
            ax=axes[i],
            annot=True, fmt=".2f",
            linewidths=.4, linecolor="0.8",
            cmap=cmap, cbar=False,
            yticklabels=metric_df.index if i == 0 else False
        )

        axes[i].set_title(m, fontsize=16, pad=12)
        axes[i].set_xlabel("Case")

        if i == 0:
            axes[i].set_yticklabels(metric_df.index, rotation=0, fontsize=11)
        else:
            axes[i].set_yticklabels([])

        axes[i].set_ylabel("")



    fig.suptitle(f"{VAR_TITLE_VAR[var]}", fontsize=18, y=0.97)

    heatmaps_dir = os.path.join(os.getcwd(), "figures", "heatmaps")
    os.makedirs(heatmaps_dir, exist_ok=True)
    out = os.path.join(heatmaps_dir, f"{var}_error_heatmap.png")
    fig.savefig(out, dpi=300)
    
    print(f"saved → {out}")

