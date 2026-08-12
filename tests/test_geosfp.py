import sqlite3

import numpy as np
import pandas as pd
import xarray as xr

from surfradpy.products import geosfp


SITE_ROWS = (
    (1, "inl", "Idaho Falls", 1496.5, 43.594919, -112.928561),
    (2, "bnd", "Bondville", 230.0, 40.05192, -88.37309),
    (3, "dra", "Desert Rock", 1007.0, 36.62373, -116.01947),
    (4, "gwn", "Goodwin Creek", 98.0, 34.2547, -89.8729),
    (5, "psu", "Penn. State Univ.", 376.0, 40.72012, -77.93085),
    (6, "sxf", "Sioux Falls", 473.0, 43.73403, -96.62328),
    (7, "tbl", "Table Mountain", 1689.0, 40.12498, -105.2368),
    (8, "fpe", "Fort Peck", 634.0, 48.30783, -105.1017),
)


def _create_site_database(path):
    with sqlite3.connect(path) as database:
        database.execute(
            """
            CREATE TABLE sites (
                site_id INTEGER,
                abb TEXT,
                name TEXT,
                elevation REAL,
                latitude REAL,
                longitude REAL,
                network TEXT
            )
            """
        )
        database.executemany(
            "INSERT INTO sites VALUES (?, ?, ?, ?, ?, ?, 'surfrad')",
            SITE_ROWS,
        )


def _source_dataset():
    datetime = pd.date_range("2024-01-02T00:30:00", periods=24, freq="h")
    lat = np.array([34.25, 48.5])
    lon = np.array([-116.25, -77.75])
    values = (
        np.arange(datetime.size)[:, np.newaxis, np.newaxis]
        + lat[np.newaxis, :, np.newaxis]
        + lon[np.newaxis, np.newaxis, :]
    )
    return xr.Dataset(
        {
            "to3": (("datetime", "lat", "lon"), values.copy()),
            "extra": (("datetime", "lat", "lon"), values.copy()),
        },
        coords={"datetime": datetime, "lat": lat, "lon": lon},
        attrs={
            "geos_fp_collection": "tavg1_2d_slv_Nx",
            "geos_fp_stream": "assim",
        },
    )


def test_product_downloads_only_to3_interpolates_and_writes(tmp_path):
    database_path = tmp_path / "surfrad.db"
    _create_site_database(database_path)
    calls = []

    class FakeGeosFp:
        paths = [tmp_path / "geos_fp_to3_subset.nc4"]

        def download(self, start, end, **kwargs):
            calls.append((start, end, kwargs))
            return _source_dataset()

    product = geosfp.GeosFpSurfrad(
        path_out=tmp_path / "output" / "v{version}",
        start="2024-01-02",
        end="2024-01-02",
        path2surfrad_database=database_path,
        download_dir=tmp_path / "download",
        verbose=False,
    )
    product._geos_fp = FakeGeosFp()
    result = product.process_row(iloc=0, save=True)

    start, end, download_kwargs = calls[0]
    assert start.isoformat() == "2024-01-02"
    assert end.isoformat() == "2024-01-02"
    assert download_kwargs["variables"] == ["to3"]
    assert download_kwargs["bbox"] == (-116.25, 34.25, -77.75, 48.5)

    expected_path = (
        tmp_path
        / "output"
        / "v0.1"
        / "2024"
        / "surfrad_geosfp_tavg1_2d_slv_Nx_20240102.nc"
    )
    assert result["p2f_out"] == expected_path
    assert expected_path.is_file()
    assert result["dsout"].sizes == {"datetime": 24, "site": 8}
    assert list(result["dsout"].data_vars) == ["to3"]

    with xr.open_dataset(expected_path) as saved:
        assert list(saved.data_vars) == ["to3"]
        assert saved.to3.dims == ("datetime", "site")
        assert saved.attrs["product_version"] == "0.1"
        assert saved.attrs["interpolation_method"] == (
            "linear latitude-longitude"
        )
        assert saved.attrs["parent_files"].endswith(
            "geos_fp_to3_subset.nc4"
        )

    assert product.workplan.empty


def test_interpolate_to_sites_requires_geos_fp_coordinates():
    sites = xr.Dataset(coords={"site": ["test"], "lat": ("site", [40.0])})

    with np.testing.assert_raises_regex(
        ValueError,
        "GEOS-FP dataset is missing coordinates: lat, lon",
    ):
        geosfp.interpolate_to_sites(
            xr.Dataset(coords={"datetime": [pd.Timestamp("2024-01-01")]}),
            sites,
        )
