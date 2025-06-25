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
    {"set_name": "SET7", "physics": "mp28_pbl1"},
    {"set_name": "SET8", "physics": "mp28_pbl5"},
    {"set_name": "SET9", "physics": "mp28_pbl7"},
    {"set_name": "SET10", "physics": "mp10_pbl1"},
    {"set_name": "SET11", "physics": "mp10_pbl5"},
    {"set_name": "SET12", "physics": "mp10_pbl7"},
    {"set_name": "SET13", "physics": "mp24_pbl1"},
    {"set_name": "SET14", "physics": "mp24_pbl5"},
    {"set_name": "SET15", "physics": "mp24_pbl7"},
])

cases = ["20241123", "20231102", "20231129", "20240120"]
stations = ["18101", "18396", "18399", "17061", "17062", "17064", "17065", "17610", "17636", "18099", "18100", "18101", "18396", "18399"]

case = "20241123"
station_id = "17061"
csv_files_path = os.path.join(os.getcwd(), "csv_files")

for case in cases:
    for station in stations:
        try:
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

            plot_timeseries_15(case = case, station_id = station, variable = "tp", fig_path = os.path.join(os.getcwd(), "figures", case))
        except:
            print("Error")