from netCDF4 import Dataset

def inject_missing_coordinates_from_geo_em(wrfout_path, geo_em_path):
    wrf_ds = Dataset(wrfout_path, "r+")
    geo_ds = Dataset(geo_em_path, "r")

    coord_vars = {
        "XLAT": "XLAT_M",
        "XLONG": "XLONG_M",
        "XLAT_U": "XLAT_U",
        "XLONG_U": "XLONG_U",
        "XLAT_V": "XLAT_V",
        "XLONG_V": "XLONG_V",
        "HGT": "HGT_M",
    }

    for wrf_var, geo_var in coord_vars.items():
        if wrf_var not in wrf_ds.variables:
            print(f"Injecting missing variable: {wrf_var}")
            geo_data = geo_ds.variables[geo_var][0, :, :]  # Remove Time dim

            if "U" in wrf_var:
                dims = ("Time", "south_north", "west_east_stag")
            elif "V" in wrf_var:
                dims = ("Time", "south_north_stag", "west_east")
            else:
                dims = ("Time", "south_north", "west_east")
            var = wrf_ds.createVariable(wrf_var, "f4", dims)
            var[0, :, :] = geo_data

    wrf_ds.close()
    geo_ds.close()
