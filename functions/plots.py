from matplotlib.colors import BoundaryNorm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
from functions.get_physics import get_mp_and_pbl
from wrf import getvar, to_np, latlon_coords, get_cartopy, extract_times, ALL_TIMES
import numpy as np
import pandas as pd
import gc
import cmaps
import matplotlib.dates as mdates
import os
from datetime import timedelta

def plot_T2_15facet(df, domain, timestep, levels = np.arange(-5, 22.5, 2.5), cmap = cmaps.precip2_17lev, fig_path = os.getcwd()):
    norm = BoundaryNorm(levels, ncolors=cmap.N)


    fig = plt.figure(figsize=(20, 12)) 
    gs = gridspec.GridSpec(4, 5, height_ratios=[1, 1, 1, 0.05], hspace=0.3)

    mesh = None

    for i in range(15):
        set_name = f"SET{i+1}"
        ds = df[set_name][domain]

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
    cbar = fig.colorbar(mesh, cax=cbar_ax, orientation='horizontal', ticks=levels)
    cbar.set_label("°C", fontsize=18)
    cbar.ax.tick_params(labelsize=16)

    fig.suptitle(f"2 Meters Temperatures\n{pd.Timestamp(extract_times(ds, timeidx=timestep)).strftime('%d/%m/%Y %H')}Z", fontsize=22, ha='center')

    plt.show()
    os.makedirs(fig_path, exist_ok=True)
    fig.savefig(os.path.join(fig_path, f"T2_{domain}_15facet.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
    gc.collect()

def plot_TP_15facet(df, domain, timestep, levels = np.arange(0, 19*5, 5), cmap = cmaps.precip2_17lev, fig_path = os.getcwd()):
    norm = BoundaryNorm(levels, ncolors=cmap.N)

    fig = plt.figure(figsize=(20, 12)) 
    gs = gridspec.GridSpec(4, 5, height_ratios=[1, 1, 1, 0.05], hspace=0.3)
    mesh = None
    row_labels = ["PBL=1", "PBL=5", "PBL=7"]  

    for plot_index, i in enumerate([1, 4, 7, 10, 13, 2, 5, 8, 11, 14, 3, 6, 9, 12, 15]):
        set_name = f"SET{i}"
        ds = df[set_name][domain]
        mp = ds.MP_PHYSICS

        tp = getvar(ds, "RAINC", timeidx=timestep) + getvar(ds, "RAINNC", timeidx=timestep)
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

    fig.suptitle(f"Total Precipitation\n{pd.Timestamp(extract_times(ds, timeidx=0)).strftime('%d/%m/%Y %H')} - {pd.Timestamp(extract_times(ds, timeidx=timestep)).strftime('%d/%m/%Y %H')}Z", fontsize=22, ha='center')

    plt.show()
    os.makedirs(fig_path, exist_ok=True)
    fig.savefig(os.path.join(fig_path, f"TP_{domain}_15facet.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
    gc.collect()

def plot_PSFC_15facet(df, domain, timestep, cmap = cmaps.precip2_17lev, fig_path = os.getcwd()):
    fig = plt.figure(figsize=(20, 12)) 
    gs = gridspec.GridSpec(4, 5, height_ratios=[1, 1, 1, 0.05], hspace=0.3)

    mesh = None

    for i in range(15):
        set_name = f"SET{i+1}"
        ds = df[set_name][domain]

        psfc = getvar(ds, "PSFC", timeidx=timestep) / 100
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

    fig.suptitle(f"Surface Pressure\n{pd.Timestamp(extract_times(ds, timeidx=timestep)).strftime('%d/%m/%Y %H')}Z", fontsize=22, ha='center')

    plt.show()
    os.makedirs(fig_path, exist_ok=True)
    fig.savefig(os.path.join(fig_path, f"PSFC_{domain}_15facet.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
    gc.collect()

def plot_WS10_15facet(df, domain, timestep, levels = np.arange(0, 16.5, 1.5), cmap = cmaps.precip2_17lev, fig_path = os.getcwd()):
    norm = BoundaryNorm(levels, ncolors=cmap.N)

    fig = plt.figure(figsize=(20, 12)) 
    gs = gridspec.GridSpec(4, 5, height_ratios=[1, 1, 1, 0.05], hspace=0.3)

    mesh = None

    for i in range(15):
        set_name = f"SET{i+1}"
        ds = df[set_name][domain]

        ws10 = getvar(ds, "wspd_wdir10", timeidx=timestep)[0]
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

    fig.suptitle(f"10 Meters Wind Speed\n{pd.Timestamp(extract_times(ds, timeidx=timestep)).strftime('%d/%m/%Y %H')}Z", fontsize=22, ha='center')

    plt.show()
    os.makedirs(fig_path, exist_ok=True)
    fig.savefig(os.path.join(fig_path, f"WS10_{domain}_15facet.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
    gc.collect()

def plot_RH2_15facet(df, domain, timestep, levels = np.arange(0, 110, 10), cmap = cmaps.precip2_17lev, fig_path = os.getcwd()):
    norm = BoundaryNorm(levels, ncolors=cmap.N)

    fig = plt.figure(figsize=(20, 12)) 
    gs = gridspec.GridSpec(4, 5, height_ratios=[1, 1, 1, 0.05], hspace=0.3)

    mesh = None

    for i in range(15):
        set_name = f"SET{i+1}"
        ds = df[set_name][domain]
        
        t2 = getvar(ds, "T2", timeidx=timestep)                
        psfc = getvar(ds, "PSFC", timeidx=timestep)            
        qv = getvar(ds, "QVAPOR", timeidx=timestep)   

        e = (qv * psfc) / (0.622 + 0.378 * qv)  

        t2_C = t2 - 273.15
        es = 6.112 * np.exp((17.67 * t2_C) / (t2_C + 243.5)) * 100  

        rh2 = 100 * (e / es)
        rh2 = np.clip(rh2, 0, 100)  

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

    fig.suptitle(f"Relative Humidity\n{pd.Timestamp(extract_times(ds, timeidx=timestep)).strftime('%d/%m/%Y %H')}Z", fontsize=22, ha='center')

    plt.show()
    os.makedirs(fig_path, exist_ok=True)
    fig.savefig(os.path.join(fig_path, f"RH2_{domain}_15facet.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
    gc.collect()

def plot_PRECIP_15facet(df, domain, timestep, levels = [0.1, 0.2, 0.5, 1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50], cmap = cmaps.precip2_17lev, fig_path = os.getcwd()):
    norm = BoundaryNorm(levels, ncolors=cmap.N)

    fig = plt.figure(figsize=(20, 12)) 
    gs = gridspec.GridSpec(4, 5, height_ratios=[1, 1, 1, 0.05], hspace=0.3)
    mesh = None
    row_labels = ["PBL=1", "PBL=5", "PBL=7"]  

    for plot_index, i in enumerate([1, 4, 7, 10, 13, 2, 5, 8, 11, 14, 3, 6, 9, 12, 15]):
        set_name = f"SET{i}"
        ds = df[set_name][domain]
        mp = ds.MP_PHYSICS

        tp = getvar(ds, "RAINC", timeidx=ALL_TIMES) + getvar(ds, "RAINNC", timeidx=ALL_TIMES)
        tp_deacc = tp.diff(dim='Time', label='upper')
        lats, lons = latlon_coords(getvar(ds, "RAINC"))
        proj = get_cartopy(getvar(ds, "RAINC"))

        row, col = divmod(plot_index, 5)
        ax = fig.add_subplot(gs[row, col], projection=proj)
        ax.coastlines()
        ax.set_title(f"MP={mp}", fontsize=18)

        mesh = ax.pcolormesh(to_np(lons), to_np(lats), to_np(tp_deacc.isel(Time=timestep)),
                            transform=ccrs.PlateCarree(), shading="auto", cmap=cmap, norm=norm)

        if col == 0:
            ax.text(-0.15, 0.5, row_labels[row], fontsize=18,
                    va='center', ha='right', transform=ax.transAxes, rotation=90)

    cbar_ax = fig.add_subplot(gs[3, :])
    cbar = fig.colorbar(mesh, cax=cbar_ax, orientation='horizontal', ticks=levels, extend="max")
    cbar.set_label("mm", fontsize=18)
    cbar.ax.tick_params(labelsize=16)

    fig.suptitle(f"Hourly Precipitation\n{pd.Timestamp(extract_times(ds, timeidx=timestep)).strftime('%d/%m/%Y %HZ')}", fontsize=22, ha='center')

    os.makedirs(fig_path, exist_ok=True)
    fig.savefig(os.path.join(fig_path, f"PRECIP_{domain}_ts_{timestep}_15facet.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
    gc.collect()
    plt.close()

def plot_timeseries_15(df, variable, fig_path = os.getcwd()):
    """
    df must contain all the dataframes which belongs to a specific station's specific case. 
    """

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

    start_date = plot_dates[0]
    end_date = plot_dates[len(plot_dates) - 1]
    station_id = df["SET1"]["Istasyon_Numarasi"].iloc[0]
    case = (start_date + timedelta(1)).strftime("%Y%m%d")

    ax.set_ylabel(y_label, fontsize=22)
    ax.set_title(f'{title}\nStation: {station_id}', fontsize=28)
    ax.legend(fontsize=18)
    ax.grid(True)
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.6)
    ax.xaxis.grid(False)

    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.set_xlim((start_date + timedelta(hours = -3)), (end_date + timedelta(hours = 3)))


    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 3))) 
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d %HZ'))

    legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
              ncol=5, fontsize=20, frameon=False)

    for legline in legend.get_lines():
        legline.set_linewidth(3)

    fig.autofmt_xdate()

    plt.tight_layout()

    os.makedirs(fig_path, exist_ok=True)
    fig.savefig(os.path.join(fig_path, f"{case}_{station_id}_{variable}.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)

    plt.show()