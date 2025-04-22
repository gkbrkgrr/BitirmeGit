import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

geo_d01 = xr.open_dataset(r"D:\istanbul_wrfouts\geo_em.d01.nc")
geo_d02 = xr.open_dataset(r"D:\istanbul_wrfouts\geo_em.d02.nc")  # if nested

attrs = geo_d01.attrs
proj_id = attrs.get("MAP_PROJ", 1)

if proj_id == 1:
    wrf_proj = ccrs.LambertConformal(
        central_longitude=attrs['STAND_LON'],
        central_latitude=attrs['CEN_LAT'],
        standard_parallels=[attrs['TRUELAT1'], attrs['TRUELAT2']]
    )
elif proj_id == 2:
    wrf_proj = ccrs.PolarStereo(
        central_longitude=attrs['STAND_LON'],
        true_scale_latitude=attrs['TRUELAT1']
    )
elif proj_id == 3:
    wrf_proj = ccrs.Mercator(
        central_longitude=attrs['STAND_LON'],
        latitude_true_scale=attrs['TRUELAT1']
    )
else:
    raise ValueError("Unsupported WRF projection MAP_PROJ")

def get_extent(ds):
    lats = ds['XLAT_M'][0, :, :]
    lons = ds['XLONG_M'][0, :, :]
    return [
        float(lons[0, 0]),     
        float(lons[0, -1]),   
        float(lats[-1, 0]),    
        float(lats[0, 0])      
    ]

extent_d01 = get_extent(geo_d01)
extent_d02 = get_extent(geo_d02)

fig = plt.figure(figsize=(12, 10))
ax = plt.axes(projection=wrf_proj)

ax.set_extent(extent_d01, crs=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAKES, alpha=0.4)
ax.gridlines(draw_labels=True)

ax.plot(geo_d01['XLONG_M'][0, 0, :], geo_d01['XLAT_M'][0, 0, :], 'k-', transform=ccrs.PlateCarree())
ax.plot(geo_d01['XLONG_M'][0, -1, :], geo_d01['XLAT_M'][0, -1, :], 'k-', transform=ccrs.PlateCarree())
ax.plot(geo_d01['XLONG_M'][0, :, 0], geo_d01['XLAT_M'][0, :, 0], 'k-', transform=ccrs.PlateCarree())
ax.plot(geo_d01['XLONG_M'][0, :, -1], geo_d01['XLAT_M'][0, :, -1], 'k-', transform=ccrs.PlateCarree())

ax.plot(geo_d02['XLONG_M'][0, 0, :], geo_d02['XLAT_M'][0, 0, :], 'r-', transform=ccrs.PlateCarree())
ax.plot(geo_d02['XLONG_M'][0, -1, :], geo_d02['XLAT_M'][0, -1, :], 'r-', transform=ccrs.PlateCarree())
ax.plot(geo_d02['XLONG_M'][0, :, 0], geo_d02['XLAT_M'][0, :, 0], 'r-', transform=ccrs.PlateCarree())
ax.plot(geo_d02['XLONG_M'][0, :, -1], geo_d02['XLAT_M'][0, :, -1], 'r-', transform=ccrs.PlateCarree())

ax.text(extent_d01[0], extent_d01[2]+0.05, "D01", fontsize=16, color='black', transform=ccrs.PlateCarree())
ax.text(extent_d02[0], extent_d02[2]+0.05, "D02", fontsize=16, color='red', transform=ccrs.PlateCarree())

plt.title("Domains", fontsize=26)
plt.tight_layout()
plt.savefig("wrf_domain_layout.png", dpi=300)
