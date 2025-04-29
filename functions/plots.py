from matplotlib.colors import BoundaryNorm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
from functions.get_physics import get_mp_and_pbl
from wrf import getvar, to_np, latlon_coords, get_cartopy, extract_times
import numpy as np
import pandas as pd
import gc
import cmaps
import matplotlib.dates as mdates

def plot_T2_15facet(df, domain, timestep, levels = np.arange(-5, 22.5, 2.5), cmap = cmaps.precip2_17lev, fig_path = "./"):
    norm = BoundaryNorm(levels, ncolors=cmap.N)


    fig = plt.figure(figsize=(20, 12)) 
    gs = gridspec.GridSpec(4, 5, height_ratios=[1, 1, 1, 0.05], hspace=0.3)

    mesh = None

    for i in range(15):
        set_name = f"SET{i+1}"
        ds = df[set_name][domain][timestep]

        t2_raw = getvar(ds, "T2", timeidx=timestep)
        t2 = t2_raw - 273.15
        lats, lons = latlon_coords(t2_raw)
        proj = get_cartopy(t2_raw)

        row, col = divmod(i, 5)
        ax = fig.add_subplot(gs[row, col], projection=proj)
        ax.coastlines()
        ax.set_title(set_name, fontsize=14)

        mesh = ax.pcolormesh(to_np(lons), to_np(lats), to_np(t2),
                            transform=ccrs.PlateCarree(), shading="auto", cmap=cmap, norm=norm)

    cbar_ax = fig.add_subplot(gs[3, :])
    cbar = fig.colorbar(mesh, cax=cbar_ax, orientation='horizontal')
    cbar.set_label("°C", fontsize=18)
    cbar.ax.tick_params(labelsize=16)

    fig.suptitle(f"2 Meters Temperatures\n{pd.Timestamp(extract_times(df[set_name][domain], timeidx=timestep)).strftime('%d/%m/%Y %H')}Z", fontsize=22, ha='center')

    plt.show()
    fig.savefig(f"{fig_path}T2_{domain}_15facet.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    gc.collect()

def plot_TP_15facet(df, domain, timestep, levels = np.arange(0, 19*5, 5), cmap = cmaps.precip2_17lev, fig_path = "./"):
    norm = BoundaryNorm(levels, ncolors=cmap.N)


    fig = plt.figure(figsize=(20, 12)) 
    gs = gridspec.GridSpec(4, 5, height_ratios=[1, 1, 1, 0.05], hspace=0.3)

    mesh = None

    row_labels = ["PBL=1", "PBL=5", "PBL=7"]  # Customize this!

    for plot_index, i in enumerate([1, 4, 7, 10, 13, 2, 5, 8, 11, 14, 3, 6, 9, 12, 15]):
        set_name = f"SET{i}"
        ds = df[set_name][domain][timestep]
        mp = ds.MP_PHYSICS

        tp = getvar(ds, "RAINC") + getvar(ds, "RAINNC")
        lats, lons = latlon_coords(getvar(ds, "RAINC"))
        proj = get_cartopy(getvar(ds, "RAINC"))

        row, col = divmod(plot_index, 5)
        ax = fig.add_subplot(gs[row, col], projection=proj)
        ax.coastlines()
        ax.set_title(f"MP={mp}", fontsize=18)

        mesh = ax.pcolormesh(to_np(lons), to_np(lats), to_np(tp),
                            transform=ccrs.PlateCarree(), shading="auto", cmap=cmap, norm=norm)

        if col == 0:
            ax.text(-0.15, 0.5, row_labels[row], fontsize=18,
                    va='center', ha='right', transform=ax.transAxes, rotation=90)

    cbar_ax = fig.add_subplot(gs[3, :])
    cbar = fig.colorbar(mesh, cax=cbar_ax, orientation='horizontal', ticks=levels)
    cbar.set_label("mm", fontsize=18)
    cbar.ax.tick_params(labelsize=16)

    fig.suptitle(f"Total Precipitation\n{pd.Timestamp(extract_times(df[set_name][domain], timeidx=0)).strftime('%d/%m/%Y %H')} - {pd.Timestamp(extract_times(df[set_name][domain], timeidx=timestep)).strftime('%d/%m/%Y %H')}Z", fontsize=22, ha='center')

    plt.show()
    fig.savefig(f"{fig_path}TP_{domain}_15facet.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    gc.collect()

def plot_PSFC_15facet(df, domain, timestep, cmap = cmaps.precip2_17lev, fig_path = "./"):
    fig = plt.figure(figsize=(20, 12)) 
    gs = gridspec.GridSpec(4, 5, height_ratios=[1, 1, 1, 0.05], hspace=0.3)

    mesh = None

    for i in range(15):
        set_name = f"SET{i+1}"
        ds = df[set_name][domain][timestep]

        psfc = getvar(ds, "PSFC") / 100
        lats, lons = latlon_coords(getvar(ds, "PSFC") )
        proj = get_cartopy(getvar(ds, "PSFC") )

        row, col = divmod(i, 5)
        ax = fig.add_subplot(gs[row, col], projection=proj)
        ax.coastlines()
        ax.set_title(set_name, fontsize=14)

        mesh = ax.pcolormesh(to_np(lons), to_np(lats), to_np(psfc),
                            transform=ccrs.PlateCarree(), shading="auto", cmap=cmap)

    cbar_ax = fig.add_subplot(gs[3, :])
    cbar = fig.colorbar(mesh, cax=cbar_ax, orientation='horizontal')
    cbar.set_label("hPa", fontsize=18)
    cbar.ax.tick_params(labelsize=16)

    fig.suptitle(f"Surface Pressure\n{pd.Timestamp(extract_times(df[set_name][domain], timeidx=timestep)).strftime('%d/%m/%Y %H')}Z", fontsize=22, ha='center')

    plt.show()
    fig.savefig(f"{fig_path}PSFC_{domain}_15facet.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    gc.collect()

def plot_WS10_15facet(df, domain, timestep, levels = np.arange(0, 16.5, 1.5), cmap = cmaps.precip2_17lev, fig_path = "./"):
    norm = BoundaryNorm(levels, ncolors=cmap.N)

    fig = plt.figure(figsize=(20, 12)) 
    gs = gridspec.GridSpec(4, 5, height_ratios=[1, 1, 1, 0.05], hspace=0.3)

    mesh = None

    for i in range(15):
        set_name = f"SET{i+1}"
        ds = df[set_name][domain][timestep]

        ws10 = getvar(ds, "wspd_wdir10")[0]
        lats, lons = latlon_coords(ws10)
        proj = get_cartopy(ws10)

        row, col = divmod(i, 5)
        ax = fig.add_subplot(gs[row, col], projection=proj)
        ax.coastlines()
        ax.set_title(set_name, fontsize=14)

        mesh = ax.pcolormesh(to_np(lons), to_np(lats), to_np(ws10),
                            transform=ccrs.PlateCarree(), shading="auto", cmap=cmap, norm=norm)

    cbar_ax = fig.add_subplot(gs[3, :])
    cbar = fig.colorbar(mesh, cax=cbar_ax, orientation='horizontal', ticks=levels)
    cbar.set_label("m/s", fontsize=18)
    cbar.ax.tick_params(labelsize=16)

    fig.suptitle(f"10 Meters Wind Speed\n{pd.Timestamp(extract_times(df[set_name][domain], timeidx=timestep)).strftime('%d/%m/%Y %H')}Z", fontsize=22, ha='center')

    plt.show()
    fig.savefig(f"{fig_path}WS10_{domain}_15facet.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    gc.collect()

def plot_RH2_15facet(df, domain, timestep, levels = np.arange(0, 110, 10), cmap = cmaps.precip2_17lev, fig_path = "./"):
    norm = BoundaryNorm(levels, ncolors=cmap.N)

    fig = plt.figure(figsize=(20, 12)) 
    gs = gridspec.GridSpec(4, 5, height_ratios=[1, 1, 1, 0.05], hspace=0.3)

    mesh = None

    for i in range(15):
        set_name = f"SET{i+1}"
        ds = df[set_name][domain][timestep]
        
        t2 = getvar(ds, "T2")                
        psfc = getvar(ds, "PSFC")            
        qv = getvar(ds, "QVAPOR")[0, :, :]   

        e = (qv * psfc) / (0.622 + 0.378 * qv)  

        t2_C = t2 - 273.15
        es = 6.112 * np.exp((17.67 * t2_C) / (t2_C + 243.5)) * 100  

        rh2 = 100 * (e / es)
        rh2 = np.clip(rh2, 0, 100)  #

        lats, lons = latlon_coords(getvar(ds, "T2"))
        proj = get_cartopy(getvar(ds, "T2"))

        row, col = divmod(i, 5)
        ax = fig.add_subplot(gs[row, col], projection=proj)
        ax.coastlines()
        ax.set_title(set_name, fontsize=14)

        mesh = ax.pcolormesh(to_np(lons), to_np(lats), to_np(rh2),
                            transform=ccrs.PlateCarree(), shading="auto", cmap=cmap, norm=norm)

    cbar_ax = fig.add_subplot(gs[3, :])
    cbar = fig.colorbar(mesh, cax=cbar_ax, orientation='horizontal', ticks=levels)
    cbar.set_label("%", fontsize=18)
    cbar.ax.tick_params(labelsize=16)

    fig.suptitle(f"Relative Humidity\n{pd.Timestamp(extract_times(df[set_name][domain], timeidx=timestep)).strftime('%d/%m/%Y %H')}Z", fontsize=22, ha='center')

    plt.show()
    fig.savefig(f"{fig_path}RH2_{domain}_15facet.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    gc.collect()

def plot_timeseries_15(df, variable):
    if variable == "t2m":
        variable_in_use = "2_Metre_Sicaklik(C)"
        title = "2 Meters Temperatures"
        y_label = "Temperature (°C)"
    elif variable == "tp":
        variable_in_use = "Yagis(mm)"
        title = "Hourly Precipitations"
        y_label = "Precipitation (mm)"
    elif variable == "ws10":
        variable_in_use = "10_Metre_Ruzgar(m/s)"
        title = "10 Meters Wind Speeds"
        y_label = "Wind Speed (m/s)"
    else:
        raise ValueError("Invalid variable input.")
    
    fig, ax = plt.subplots(figsize=(20, 12))

    for set_name, data_frame in df.items():
        if not data_frame.empty and 'Tarih' in data_frame.columns and variable_in_use in data_frame.columns:
            try:
                mp_val, pbl_val = get_mp_and_pbl(set_name)

                plot_dates = pd.to_datetime(data_frame['Tarih'])
                plot_temps = data_frame[variable_in_use]
                ax.plot(plot_dates, plot_temps, label=f"MP={mp_val} PBL={pbl_val}")
            except IndexError:
                print(f"Warning: Could not find physics mapping for {set_name}. Skipping label generation.")
            except Exception as e:
                print(f"Warning: Could not plot data for {set_name}. Error: {e}")
        else:
            print(f"Warning: Skipping {set_name} due to missing data or columns.")


    ax.set_xlabel("Date", fontsize=18)
    ax.set_ylabel(y_label, fontsize=18)
    ax.set_title(f'{title}\nStation: {df["SET1"]["Istasyon_Numarasi"].iloc[0]}', fontsize=24)
    ax.legend(fontsize=14)
    ax.grid(True)

    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    start_date = pd.to_datetime("2024-11-22 00:00:00")
    end_date = pd.to_datetime("2024-11-25 00:00:00")
    ax.set_xlim(start_date, end_date)


    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 3))) 
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))

    fig.autofmt_xdate()

    plt.tight_layout()
    plt.show()