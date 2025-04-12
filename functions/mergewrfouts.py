from collections import defaultdict
import os
import xarray as xr

def wrfoutMerger(wrfout_files: list[str]):
    """
    Girdi: wrfout dosyalarının listesi.

    Çıktı: her domain için ayrı xarray.Dataset.
    
    Girdi olarak d01 ve d02 wrfoutları aynı listede verildiyse
    2 ayrı xarray.Dataset çıktısı verilir.
    """
    domain_files = defaultdict(list)
    for f in wrfout_files:
        if "d01" in os.path.basename(f):
            domain_files["d01"].append(f)
        elif "d02" in os.path.basename(f):
            domain_files["d02"].append(f)

    for dom in domain_files:
        domain_files[dom] = sorted(domain_files[dom])

    ds1 = xr.open_mfdataset(
    domain_files["d01"],
    combine = "nested",
    concat_dim="Time",         
    parallel=True,
    preprocess=lambda d: d.drop_vars([v for v in d.data_vars if 'MemoryOrder' in d[v].attrs.get('coordinates', '')]) if 'MemoryOrder' in str(d) else d
    )

    if len(domain_files["d02"]) != 0:
        ds2 = xr.open_mfdataset(
        domain_files["d02"],
        combine = "nested",
        concat_dim="Time",         
        parallel=True,
        preprocess=lambda d: d.drop_vars([v for v in d.data_vars if 'MemoryOrder' in d[v].attrs.get('coordinates', '')]) if 'MemoryOrder' in str(d) else d
        )
        return ds1, ds2
    return ds1