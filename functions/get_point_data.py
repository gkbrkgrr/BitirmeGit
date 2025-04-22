import numpy as np
import pandas as pd
from wrf import getvar, ALL_TIMES, ll_to_xy, to_np

def get_point_data(wrfin, var, lat, lon):
    supported_vars = ["tp", "t2", "psfc", "ws10"]
    if var == "tp":
        x_y = ll_to_xy(wrfin, latitude=lat, longitude=lon)
        rainc = getvar(wrfin, "RAINC", timeidx=ALL_TIMES)
        rainnc = getvar(wrfin, "RAINNC", timeidx=ALL_TIMES)
        values = rainc + rainnc

        point_values = values[:, x_y[1], x_y[0]]
        deacc_values = np.diff(to_np(point_values).astype(int), axis=0)

        times = getvar(wrfin, "times", timeidx=ALL_TIMES).to_index()
        time_index = times[:-1]  # match deacc shape

        return pd.Series(deacc_values, index=time_index, name="tp")
    elif var in supported_vars and var != "tp":
        x_y = ll_to_xy(wrfin, latitude=lat, longitude=lon)
        if var == "ws10":
            values = getvar(wrfin, "wspd_wdir10", timeidx=ALL_TIMES)[0]
        else: 
            values = getvar(wrfin, var.upper(), timeidx=ALL_TIMES)

        time_index = getvar(wrfin, "times", timeidx=ALL_TIMES).to_index()[:-1]
        point_values = to_np(values[:, x_y[1], x_y[0]])

        if var == "t2": return pd.Series(np.round(point_values - 273.15, 2)[:-1], time_index, name=var.upper())  
        elif var == "psfc":  return pd.Series(np.round(point_values/100).astype(int)[:-1], time_index, name=var.upper())  
        elif var == "ws10":  return pd.Series(np.round(point_values).astype(int)[:-1], time_index, name=var.upper())  
    else:
        raise ValueError(f"Variable '{var}' is not supported. Supported variables are: {supported_vars}")