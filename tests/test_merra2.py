import sqlite3

import numpy as np
import pandas as pd
import xarray as xr

from surfradpy.products import merra2


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
    lat = np.array([34.0, 49.0])
    lon = np.array([-117.0, -77.5])
    values = (
        np.arange(datetime.size)[:, np.newaxis, np.newaxis]
        + lat[np.newaxis, :, np.newaxis]
        + lon[np.newaxis, np.newaxis, :]
    )
    return xr.Dataset(
        {
            variable: (("datetime", "lat", "lon"), values.copy())
            for variable in merra2.M2T1NXSLV_VARIABLES
        },
        coords={"datetime": datetime, "lat": lat, "lon": lon},
        attrs={
            "merra2_short_name": "M2T1NXSLV",
            "merra2_version": "5.12.4",
        },
    )


def test_m2t1nxslv_variable_inventory_is_complete():
    assert len(merra2.M2T1NXSLV_VARIABLES) == 47
    assert len(set(merra2.M2T1NXSLV_VARIABLES)) == 47
    assert {
        "CLDPRS",
        "PS",
        "T2M",
        "TO3",
        "TQV",
        "U10M",
        "V10M",
        "ZLCL",
    }.issubset(merra2.M2T1NXSLV_VARIABLES)


def test_sites_and_spatial_subset_come_from_database(tmp_path):
    database_path = tmp_path / "surfrad.db"
    _create_site_database(database_path)

    sites = merra2.load_surfrad_sites(database_path)

    assert sites.site.to_numpy().tolist() == list(merra2.SURFRAD_SITE_CODES)
    assert sites.sizes["site"] == 8
    assert sites.site_name.sel(site="tbl").item() == "Table Mountain"
    assert merra2.get_spatial_subset(sites) == (
        -116.25,
        34.0,
        -77.5,
        48.5,
    )


def test_interpolate_to_sites_uses_paired_coordinates():
    sites = xr.Dataset(
        coords={
            "site": ["west", "east"],
            "site_id": ("site", [1, 2]),
            "site_name": ("site", ["West", "East"]),
            "lat": ("site", [34.25, 48.25]),
            "lon": ("site", [-116.0, -78.0]),
            "elevation": ("site", [100.0, 200.0]),
        }
    )

    dsout = merra2.interpolate_to_sites(_source_dataset(), sites)

    assert dsout.T2M.dims == ("datetime", "site")
    assert "lat" not in dsout.dims
    assert "lon" not in dsout.dims
    np.testing.assert_allclose(
        dsout.T2M.isel(datetime=0),
        sites.lat + sites.lon,
    )
    assert dsout.site_name.to_numpy().tolist() == ["West", "East"]


def test_product_downloads_bbox_interpolates_and_writes(tmp_path):
    database_path = tmp_path / "surfrad.db"
    _create_site_database(database_path)
    calls = []

    class FakeMerra2:
        def __init__(self, **kwargs):
            self.init_kwargs = kwargs
            self.auth = None
            self.paths = [tmp_path / "M2T1NXSLV_subset.nc4"]

        def download(self, start, end, **kwargs):
            calls.append((start, end, kwargs))
            self.auth = object()
            return _source_dataset()

    product = merra2.Merra2Surfrad(
        path_out=tmp_path / "output" / "v{version}",
        start="2024-01-02",
        end="2024-01-02",
        path2surfrad_database=database_path,
        download_dir=tmp_path / "download",
        verbose=False,
    )
    product._merra2 = FakeMerra2()
    result = product.process_row(iloc=0, save=True)

    start, end, download_kwargs = calls[0]
    assert start.isoformat() == "2024-01-02"
    assert end.isoformat() == "2024-01-02"
    assert download_kwargs["variables"] == list(
        merra2.M2T1NXSLV_VARIABLES
    )
    assert download_kwargs["bbox"] == (-116.25, 34.0, -77.5, 48.5)

    expected_path = (
        tmp_path
        / "output"
        / "v0.1"
        / "2024"
        / "surfrad_merra2_M2T1NXSLV_20240102.nc"
    )
    assert result["p2f_out"] == expected_path
    assert expected_path.is_file()
    assert result["dsout"].sizes == {"datetime": 24, "site": 8}

    with xr.open_dataset(expected_path) as saved:
        assert set(saved.data_vars) == set(
            merra2.M2T1NXSLV_VARIABLES
        )
        assert saved.T2M.dims == ("datetime", "site")
        assert saved.attrs["product_version"] == "0.1"
        assert saved.attrs["interpolation_method"] == (
            "linear latitude-longitude"
        )
        assert saved.attrs["parent_files"].endswith(
            "M2T1NXSLV_subset.nc4"
        )

    assert product.workplan.empty
