from functions.plots import plot_timeseries_15
import os
import pandas as pd

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

case = "20241123"
station_id = "17061"
csv_files_path = os.path.join(os.getcwd(), "csv_files")

matching_files = []
df = {}

if os.path.isdir(csv_files_path):
    for filename in os.listdir(csv_files_path):
        if filename.startswith(f"case_{case}") and filename.endswith(f"{station_id}.csv"):
            matching_files.append(os.path.join(csv_files_path, filename))

for file in matching_files:
    physics = file.split(f"case_{case}")[1].split("_station")[0].lstrip("_")
    set_name = sets_mapping.loc[sets_mapping["physics"] == physics, "set_name"].values[0]
    df[set_name] = pd.read_csv(file)

plot_timeseries_15(df, "t2m", os.path.join(os.getcwd(), "figures", "23112024"))