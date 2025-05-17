import numpy as np
import pandas as pd
import os
import re

sets_mapping = pd.DataFrame([
    {"set_name": "SET1", "physics": "mp4_pbl1"},
    {"set_name": "SET2", "physics": "mp4_pbl5"},
    {"set_name": "SET3", "physics": "mp4_pbl7"},
    {"set_name": "SET4", "physics": "mp6_pbl1"},
    {"set_name": "SET5", "physics": "mp6_pbl5"},
    {"set_name": "SET6", "physics": "mp6_pbl7"},
    {"set_name": "SET7", "physics": "mp38_pbl1"},
    {"set_name": "SET8", "physics": "mp38_pbl5"},
    {"set_name": "SET9", "physics": "mp38_pbl7"},
    {"set_name": "SET10", "physics": "mp10_pbl1"},
    {"set_name": "SET11", "physics": "mp10_pbl5"},
    {"set_name": "SET12", "physics": "mp10_pbl7"},
    {"set_name": "SET13", "physics": "mp24_pbl1"},
    {"set_name": "SET14", "physics": "mp24_pbl5"},
    {"set_name": "SET15", "physics": "mp24_pbl7"},
])


def errors_to_csv(case: str, variable: str, observation):
    """
    Obs verinin eklenmesi lazım. Her istasyon için ayrı obs gelecek. Belki yeni bir mapping uygulanabilir?  
    """

    csv_files_path = os.path.join(os.getcwd(), "csv_files")
    results = pd.DataFrame(columns=["Physics", "MAE", "RMSE", "Bias"])
    safe_variable = re.sub(r'[<>:"/\\|?*]', '_', variable)

    for physics in sets_mapping["physics"]: 
        matching_files = []
        df = {}

        if os.path.isdir(csv_files_path):
            for filename in os.listdir(csv_files_path):
                if f"case_{case}" in filename and f"{physics}" in filename and filename.endswith(".csv"):
                    matching_files.append(os.path.join(csv_files_path, filename))

        for file in matching_files:
            station_id = file.split("station_")[1].split(".")[0]
            df[station_id] = pd.read_csv(file)

        results_for_each_physics = pd.DataFrame(columns=["StationId", "MAE", "RMSE", "Bias"])
        for station_id in df.keys():
            results_for_each_physics.loc[len(results_for_each_physics)] = {
            "StationId": str(station_id),
            "MAE": compute_mae(df[station_id][variable], observation),
            "RMSE": compute_rmse(df[station_id][variable], observation),
            "Bias": compute_bias(df[station_id][variable], observation)
            }
        
        results.loc[len(results)] = {
            "Physics": physics,
            "MAE": round(np.nanmean(results_for_each_physics["MAE"]), 4),
            "RMSE": round(np.nanmean(results_for_each_physics["RMSE"]), 4),
            "Bias": round(np.nanmean(results_for_each_physics["Bias"]), 4),
            }
        
        error_csvs_path = os.path.join(os.getcwd(), "error_csvs")
        os.makedirs(error_csvs_path, exist_ok=True)
        
        results.to_csv(os.path.join(error_csvs_path, f"{case}_{safe_variable}.csv"), index=False, float_format="%.4f")

def compute_bias(model, obs):
    return np.nanmean(model - obs)

def compute_rmse(model, obs):
    return np.sqrt(np.nanmean((model - obs)**2))

def compute_mae(model, obs):
    return np.nanmean(np.abs(model - obs))