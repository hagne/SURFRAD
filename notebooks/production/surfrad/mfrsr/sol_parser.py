from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import xarray as xr


SECTION_RE = re.compile(r"^(SN|WE)\s*(\d+)\s*$")
NM_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)\s*nm$", re.IGNORECASE)


def read_sol_to_xarray(path: str | Path) -> xr.Dataset:
    """Read a SURFRAD/MFRSR .sol cosine file into an xarray Dataset.

    Returns a dataset with:
    - variables: `SN(channel, angle)`, `WE(channel, angle)`
    - coordinate: `channel` from the numeric suffix in section names (e.g. SN1 -> 1)
    - coordinate: `angle` from -90 to 90 degrees
    - variable: `wavelength_nm(channel)` parsed from the LAMBDAS header block
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    lambda_start = None
    lambda_end = None
    for i, line in enumerate(lines):
        tag = line.strip().upper()
        if tag == "LAMBDAS":
            lambda_start = i
        elif lambda_start is not None and tag == "END":
            lambda_end = i
            break

    if lambda_start is None or lambda_end is None:
        raise ValueError(f"Could not find complete LAMBDAS header block in {path}")

    header_rows: list[tuple[str, float]] = []
    for line in lines[lambda_start + 1 : lambda_end]:
        parts = line.split()
        if len(parts) < 2:
            continue
        label = parts[0]
        table_tag = parts[-1].upper()
        if not table_tag.startswith("TABLE"):
            continue
        m = NM_RE.search(label)
        nm = float(m.group(1)) if m else np.nan
        header_rows.append((label, nm))

    section_values: dict[str, dict[int, list[float]]] = {"SN": {}, "WE": {}}
    i = lambda_end + 1
    while i < len(lines):
        m = SECTION_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue

        var, channel_s = m.groups()
        channel = int(channel_s)
        i += 1
        values: list[float] = []

        while i < len(lines) and lines[i].strip().upper() != "END":
            for token in lines[i].split():
                values.append(float(token))
            i += 1

        section_values[var][channel] = values
        i += 1

    channels = sorted(set(section_values["SN"]) | set(section_values["WE"]))
    if not channels:
        raise ValueError(f"No SN/WE channel sections found in {path}")

    max_samples = max(
        [len(v) for v in section_values["SN"].values()] +
        [len(v) for v in section_values["WE"].values()]
    )
    angle = np.linspace(-90.0, 90.0, max_samples, dtype=float)

    def to_2d(var: str) -> np.ndarray:
        arr = np.full((len(channels), max_samples), np.nan, dtype=float)
        for row, ch in enumerate(channels):
            vals = section_values[var].get(ch, [])
            if vals:
                arr[row, : len(vals)] = vals
        return arr

    ds = xr.Dataset(
        data_vars={
            "SN": (("channel", "angle"), to_2d("SN")),
            "WE": (("channel", "angle"), to_2d("WE")),
        },
        coords={
            "channel": np.asarray(channels, dtype=int),
            "angle": angle,
        },
        attrs={"source_file": str(path)},
    )

    wavelength_nm = np.full(len(channels), np.nan, dtype=float)
    for row, ch in enumerate(channels):
        idx = ch - 1
        if 0 <= idx < len(header_rows):
            wavelength_nm[row] = header_rows[idx][1]
    ds["wavelength_nm"] = ("channel", wavelength_nm)

    return ds
