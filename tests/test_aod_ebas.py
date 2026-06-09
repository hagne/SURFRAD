from io import BytesIO
from math import sqrt
from types import ModuleType, SimpleNamespace
from unittest.mock import patch

import numpy as np

from surfradpy.products import aod_ebas


def _sample_file(day_line=b"01-Jan-2024"):
    return BytesIO(
        b"SURFRAD AOD\n"
        + day_line
        + b"\n"
        + b"415 500 673\n"
        + b"daily average\n"
        + b"ozone\n"
        + b"field names\n"
        + b"0000 0 0.100 -999.9 0.300\n"
        + b"0003 1 0.200 0.400\n"
    )


def test_normalize_station_code_accepts_aod_aliases():
    assert aod_ebas.normalize_station_code("bon") == "bnd"
    assert aod_ebas.normalize_station_code("tbl") == "bos"
    assert aod_ebas.normalize_station_code("psu") == "psu"


def test_integrate_file_parses_daily_aod_data():
    data = aod_ebas.AodEbasData("bnd")
    data.integrate_file(_sample_file(), "bon_20240101.aod", tz_offset=6)
    data.sort_data(2024)

    ds = data.ds
    assert ds.sizes["datetime"] == 2
    assert ds.sizes["wavelength"] == 3
    assert ds.wavelength.to_numpy().tolist() == [415, 500, 673]
    np.testing.assert_allclose(ds.aod.isel(wavelength=0).to_numpy(), [0.1, 0.2])
    assert np.isnan(ds.aod.isel(datetime=0, wavelength=1).item())
    assert ds.cloud_screening_failure.to_numpy().tolist() == [False, True]
    assert data.epoch_ms[0] == np.datetime64("2024-01-01T06:00:00", "ms").astype(np.int64)


def test_integrate_file_falls_back_to_date_in_file_name():
    data = aod_ebas.AodEbasData("bnd")
    data.integrate_file(_sample_file(day_line=b"not a date"), "bon_20240101.aod", tz_offset=6)
    data.sort_data(2024)

    assert data.epoch_ms[0] == np.datetime64("2024-01-01T06:00:00", "ms").astype(np.int64)


def test_dataset_ignores_unlabelled_older_extra_channels():
    newer_file = BytesIO(
        b"SURFRAD AOD\n"
        b"02-Jan-2024\n"
        b"415 500\n"
        b"daily average\n"
        b"ozone\n"
        b"field names\n"
        b"0000 0 0.100 0.200\n"
    )
    older_file = BytesIO(
        b"SURFRAD AOD\n"
        b"01-Jan-2024\n"
        b"415 500 673\n"
        b"daily average\n"
        b"ozone\n"
        b"field names\n"
        b"0000 0 0.300 0.400 0.500\n"
    )
    data = aod_ebas.AodEbasData("bnd")
    data.integrate_file(newer_file, "bon_20240102.aod", tz_offset=6)
    data.integrate_file(older_file, "bon_20240101.aod", tz_offset=6)
    data.sort_data(2024)

    assert data.ds.sizes["wavelength"] == 2
    assert data.ds.wavelength.to_numpy().tolist() == [415, 500]
    np.testing.assert_allclose(data.ds.aod.to_numpy(), [[0.3, 0.4], [0.1, 0.2]])


def test_sort_data_filters_to_requested_utc_year():
    data = aod_ebas.AodEbasData("bnd")
    data.integrate_file(_sample_file(day_line=b"31-Dec-2023"), "bon_20231231.aod", tz_offset=6)
    data.sort_data(2024)

    assert data.epoch_ms.shape[0] == 0


def test_fixed_interval_averager_matches_legacy_hourly_statistics():
    times = np.array(
        [
            np.datetime64("2024-01-01T00:00:00", "ms").astype(np.int64),
            np.datetime64("2024-01-01T00:03:00", "ms").astype(np.int64),
            np.datetime64("2024-01-01T01:00:00", "ms").astype(np.int64),
        ],
        dtype=np.int64,
    )
    values = np.array([1.0, 3.0, 10.0], dtype=np.float64)
    flags = np.array([False, True, False], dtype=np.bool_)
    average = aod_ebas.FixedIntervalAverager(60 * 60 * 1000, times)

    np.testing.assert_allclose(average(values), [2.0, 10.0])
    np.testing.assert_allclose(average.quantiles(values, 0.5)[:, 0], [2.0, 10.0])
    np.testing.assert_allclose(average.stddev(values), [sqrt(2.0), np.nan], equal_nan=True)
    assert average.bitwise_or(flags).tolist() == [True, False]


def test_write_level_files_with_fake_ebas_modules():
    class FakeCharacteristics:
        def __init__(self):
            self.entries = []

        def add_parse(self, *args):
            self.entries.append(args)

    class FakeDataObject(SimpleNamespace):
        pass

    class FakeNasaAmes:
        instances = []

        def __init__(self):
            self.metadata = SimpleNamespace()
            self.variables = []
            self.sample_times = []
            self.write_calls = []
            self.instances.append(self)

        def write(self, **kwargs):
            self.write_calls.append(kwargs)

    ebas = ModuleType("ebas")
    ebas.__version__ = "test"
    ebas_io = ModuleType("ebas.io")
    ebas_io_file = ModuleType("ebas.io.file")
    ebas_io_file.nasa_ames = SimpleNamespace(EbasNasaAmes=FakeNasaAmes)
    ebas_metadata = ModuleType("ebas.io.ebasmetadata")
    ebas_metadata.DatasetCharacteristicList = FakeCharacteristics
    nilutility = ModuleType("nilutility")
    nilutility_datatypes = ModuleType("nilutility.datatypes")
    nilutility_datatypes.DataObject = FakeDataObject

    modules = {
        "ebas": ebas,
        "ebas.io": ebas_io,
        "ebas.io.file": ebas_io_file,
        "ebas.io.ebasmetadata": ebas_metadata,
        "nilutility": nilutility,
        "nilutility.datatypes": nilutility_datatypes,
    }

    data = aod_ebas.AodEbasData("bnd")
    data.integrate_file(_sample_file(), "bon_20240101.aod", tz_offset=6)
    data.sort_data(2024)

    with patch.dict("sys.modules", modules):
        data.write_level1("/tmp")
        data.write_level2(2024, "/tmp")

    level1, level2 = FakeNasaAmes.instances
    assert level1.metadata.datalevel == "1"
    assert len(level1.variables) == 3
    assert level1.variables[0].values_ == [0.1, 0.2]
    assert level1.variables[1].values_ == [None, 0.4]
    assert level1.variables[0].flags == [[], [799]]
    assert level1.write_calls == [{"createfiles": True, "destdir": "/tmp"}]

    assert level2.metadata.datalevel == "2"
    assert len(level2.sample_times) == 366 * 24
    assert len(level2.variables) == 9
    assert level2.variables[0].metadata.statistics == "arithmetic mean"
    assert level2.variables[1].metadata.statistics == "median"
    assert level2.variables[2].metadata.statistics == "stddev"
    assert level2.variables[0].flags[6] == [799]
    assert level2.variables[0].flags[0] == [999]
