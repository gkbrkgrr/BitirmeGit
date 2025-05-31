import re
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd


def parse_mgm_hourlies(
    xlsx_path: str | Path,
    out_dir: str | Path = "out",
    *,
    sheet_name: int | str = 0,
    verbose: bool = False,
) -> list[Path]:
    """
    Parse Turkish MGM hourly-observation workbooks and write one CSV per
    (year, month, station, variable).  Wind tables keep both direction & speed.

    Returns
    -------
    list[Path]          Paths of the CSVs written or updated.
    """
    xlsx_path = Path(xlsx_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # 1. Load sheet (raw strings, no header detection)
    # ------------------------------------------------------------------ #
    df = pd.read_excel(
        xlsx_path, sheet_name=sheet_name, header=None, engine="openpyxl", dtype=str
    )

    # ------------------------------------------------------------------ #
    # 2. Look-ups
    # ------------------------------------------------------------------ #
    var_lookup = {
        "Deniz Seviyesine İndirgenmiş Basınç": ("mslp", "hPa"),
        "İşba Sıcaklığı": ("td2", "degC"),
        "Rüzgar Yönü": ("wdir_ws", "deg_mps"),
        "Sıcaklık": ("t2", "degC"),
        "Toplam Yağış": ("precip", "mm"),
    }

    meta_re = re.compile(
        r"Yıl/Ay:\s*(\d{4})/(\d{1,2})\s+İstasyon Adı/No:\s*[^/]+/(\d+)"
    )

    # ------------------------------------------------------------------ #
    # 3. State variables
    # ------------------------------------------------------------------ #
    current: dict[str, str | int] = {}
    buffer: list[dict[str, str]] = []
    files_written: list[Path] = []

    # ------------------------------------------------------------------ #
    # Helper: split wind cell → (direction, speed)
    # ------------------------------------------------------------------ #
    _num_re = re.compile(r"\d+(?:[.,]\d+)?")

    def _split_wind_cell(cell: str) -> Tuple[Optional[float], Optional[float]]:
        """Return (direction, speed) with None when missing."""
        nums = _num_re.findall(cell.replace(",", "."))
        if not nums:
            return None, None

        nums = [float(n) for n in nums]
        if len(nums) == 2:                    # “140 3.2”
            return nums[0], nums[1]

        num = nums[0]
        if 0 <= num <= 360:
            return num, None                 # direction only
        return None, num                     # speed only

    # ------------------------------------------------------------------ #
    # Helper: write or append a completed table
    # ------------------------------------------------------------------ #
    def flush_buffer() -> None:
        nonlocal buffer
        if not buffer or "var" not in current:
            buffer.clear()
            return

        # wide → long
        tmp = pd.DataFrame(buffer).melt(
            id_vars="day", var_name="hour", value_name="raw"
        )
        tmp["raw"] = (
            tmp["raw"].replace(",", ".", regex=False)
            .fillna("")
            .astype(str)
            .str.strip()
        )

        if current["var"] == "wdir_ws":
            tmp[["wdir", "value"]] = tmp["raw"].apply(_split_wind_cell).tolist()
            tmp = tmp.dropna(subset=["wdir", "value"], how="all")
        else:
            tmp["value"] = pd.to_numeric(
                tmp["raw"].replace("", np.nan), errors="coerce"
            )
            tmp = tmp.dropna(subset=["value"], how="all")

        if tmp.empty:
            buffer.clear()
            return

        # make timestamps
        y, m = int(current["year"]), int(current["month"])
        tmp["day"] = pd.to_numeric(tmp["day"], errors="coerce").astype("Int64")
        tmp["hour"] = pd.to_numeric(tmp["hour"], errors="coerce").astype("Int64")
        tmp["datetime"] = pd.to_datetime(
            dict(year=y, month=m, day=tmp["day"], hour=tmp["hour"]), errors="coerce"
        )
        tmp = tmp.dropna(subset=["datetime"])

        cols = ["datetime", "value"] + (["wdir"] if "wdir" in tmp.columns else [])
        tmp = tmp[cols].sort_values("datetime")

        # write / append
        fname = (
            f"{y:04d}{m:02d}_{current['var']}_{current['unit']}_{current['station']}.csv"
        )
        path = out_dir / fname
        if path.exists():
            old = pd.read_csv(path, parse_dates=["datetime"])
            tmp = (
                pd.concat([old, tmp], ignore_index=True)
                .drop_duplicates(subset=["datetime"])
                .sort_values("datetime")
            )
        tmp.to_csv(path, index=False, float_format="%.3f")
        files_written.append(path)
        if verbose:
            print(f"✓ {path.name:38s} {len(tmp):4d} rows")

        buffer.clear()

    # ------------------------------------------------------------------ #
    # 4. Walk the sheet
    # ------------------------------------------------------------------ #
    for _, row in df.iterrows():
        cell = next((str(c).strip() for c in row if pd.notna(c) and str(c).strip()), "")
        if not cell:
            continue

        # a) metadata line
        m = meta_re.match(cell)
        if m:
            flush_buffer()
            current = {"year": m[1], "month": m[2], "station": m[3]}
            continue

        # b) variable title
        matched = next(
            ((vcode, unit) for frag, (vcode, unit) in var_lookup.items() if frag in cell),
            None,
        )
        if matched:
            flush_buffer()
            vcode, unit = matched
            current |= {"var": vcode, "unit": unit}
            continue

        # c) header row
        if cell.startswith("Gün/Saat"):
            offset = row.first_valid_index()
            hours = (
                row.iloc[offset + 1 : offset + 1 + 24]
                .astype(str)
                .str.extract(r"(\d+)")
                .dropna()[0]
                .astype(int)
                .tolist()
            )
            current |= {"hours": hours, "offset": offset}
            continue

        # d) data row
        if "var" in current and "hours" in current:
            offset = current["offset"]
            day = pd.to_numeric(row.iloc[offset], errors="coerce")
            if not np.isnan(day):
                vals = row.iloc[offset + 1 : offset + 1 + len(current["hours"])].tolist()
                buffer.append({"day": int(day), **dict(zip(current["hours"], vals))})

    flush_buffer()
    return files_written

if __name__ == "__main__":
    for file in Path("./obs_data").glob("*.xlsx"): output = parse_mgm_hourlies(
    file,
    out_dir="obs_csvs",
    verbose=True,
)
