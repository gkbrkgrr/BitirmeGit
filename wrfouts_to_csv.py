from functions.get_point_data import get_point_data
from functions.listwrfouts import listWrfouts
import pandas as pd
from netCDF4 import Dataset
import os
from wrf import getvar
from datetime import timedelta
import gc

base_dir = r"D:\istanbul_wrfouts\20012024"
sets = [f"SET{i}" for i in range(1, 16)]
df = {}

for set in sets:
    df[set] = {"d02": []}
    
    wrfouts_list = listWrfouts(
        os.path.join(base_dir, set, f"{set}_wrfouts.tar.gz"),
        os.path.join(base_dir, set, "temp_extract_dir")
    )
    
    for f in wrfouts_list:
        if "d02" in os.path.basename(f):
            df[set]["d02"].append(Dataset(f))
            gc.collect()

station_coords = pd.read_csv("station_coords.csv", delimiter=",", header=None,
                             names=["StationName", "Latitude", "Longitude", "StationID"])

csv_path = os.path.join(os.getcwd(), "csv_files")
os.makedirs(csv_path, exist_ok=True)

for set in sets:
    case = pd.to_datetime(getvar(df[set]["d02"][0], "times", timeidx=0).values)
    case = (case + timedelta(days=1)).strftime("%Y%m%d")

    for index, row in station_coords.iterrows():
        df_name = f'case_{case}_mp{df[set]["d02"][0].MP_PHYSICS}_pbl{df[set]["d02"][0].BL_PBL_PHYSICS}_station_{row["StationID"]}.csv'

        tp = get_point_data(df[set]["d02"], "tp", row["Latitude"], row["Longitude"])
        t2 = get_point_data(df[set]["d02"], "t2", row["Latitude"], row["Longitude"])
        psfc = get_point_data(df[set]["d02"], "psfc", row["Latitude"], row["Longitude"])
        ws10 = get_point_data(df[set]["d02"], "ws10", row["Latitude"], row["Longitude"])

        df_point = pd.concat([tp, t2, psfc, ws10], axis=1)
        df_point.columns = ["Yagis(mm)", "2_Metre_Sicaklik(C)", "Yuzey_Basinci(mb)", "10_Metre_Ruzgar(m/s)"]
        df_point.index.name = "Time"

        df_point["Istasyon_Adi"] = row["StationName"]
        df_point["Istasyon_Numarasi"] = row["StationID"]
        df_point["Enlem"] = row["Latitude"]
        df_point["Boylam"] = row["Longitude"]

        df_point = df_point.reset_index().rename(columns={"Time": "Tarih"})
        df_point["Tarih"] = pd.to_datetime(df_point["Tarih"])  

        column_order = [
            "Istasyon_Adi", "Istasyon_Numarasi", "Enlem", "Boylam", "Tarih",
            "Yagis(mm)", "2_Metre_Sicaklik(C)",
            "Yuzey_Basinci(mb)", "10_Metre_Ruzgar(m/s)"
        ]
        df_point = df_point[column_order]

        df_point.to_csv(os.path.join(csv_path, df_name), sep=",", index=False)
