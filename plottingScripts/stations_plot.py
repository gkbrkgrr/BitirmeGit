import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.ticker import FixedLocator

selected_codes = [17610, 17636, 17813, 17814, 18099, 17061, 17062, 17064, 17065, 17454, 18100, 18101, 18397]

df = pd.read_csv("station_coords.csv", delimiter=",", header=None, names=["name", "latitude", "longitude", "Station Code"])  
df = df[df["Station Code"].isin(selected_codes)].reset_index(drop=True)
df['Index'] = range(1, len(df) + 1)  # Reassign Index after filtering

lon_ticks = [28.0, 28.5, 29.0, 29.5, 30.0]
lat_ticks = [40.6, 40.8, 41.0, 41.2, 41.4, 41.6, 41.8]

fig = plt.figure(figsize=(12, 6), layout="compressed")
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([27.999, 30.001, 40.599, 41.803])

ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.COASTLINE)

gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(),
                  linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
gl.xlocator = FixedLocator(lon_ticks)
gl.ylocator = FixedLocator(lat_ticks)
gl.xlabel_style = {'size': 10}
gl.ylabel_style = {'size': 10}

for i, row in df.iterrows():
    ax.plot(row['longitude'], row['latitude'], marker='o', color='red', markersize=5, transform=ccrs.PlateCarree())
    ax.text(row['longitude'] + 0.009, row['latitude'] + 0.009, str(row['Index']),
            fontsize=11, color='black', transform=ccrs.PlateCarree())

handles = [
    plt.Line2D([0], [0], marker='o', color='w', label=f"{i+1:02d}. {name}",
               markerfacecolor='red', markersize=6)
    for i, name in enumerate(df['name'])
]

ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.07, 0.5),
          fontsize='medium', borderaxespad=0.)

fig.savefig("stations.jpg", format="jpg", dpi=75)