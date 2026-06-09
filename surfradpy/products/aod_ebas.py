"""SURFRAD AOD EBAS NASA Ames generation."""

from __future__ import annotations

import datetime as _dt
import gzip
import re
import shutil
from dataclasses import dataclass
from importlib import metadata as _importlib_metadata
from math import ceil, floor, isfinite, nan
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import BinaryIO

import numpy as np
import xarray as xr

AOD_STATION_REMAP: dict[str, str] = {
    "bos": "tbl",
    "bnd": "bon",
    "fpe": "fpk",
    "gcr": "gwn",
}

AOD_TIMEZONE: dict[str, int] = {
    "bnd": 6,
    "bos": 7,
    "dra": 8,
    "fpe": 7,
    "gcr": 6,
    "psu": 5,
    "sxf": 6,
}

DEFAULT_STATIONS: tuple[str, ...] = ("bnd", "fpe", "gcr", "bos", "dra", "psu", "sxf")

DEFAULT_SOURCE = Path("/outgoing_ftp/data/radiation/surfrad/aod")
DEFAULT_ARCHIVE = Path("/aer/archive/ebas")

_MONTHS = {
    b"jan": 1,
    b"feb": 2,
    b"mar": 3,
    b"apr": 4,
    b"may": 5,
    b"jun": 6,
    b"jul": 7,
    b"aug": 8,
    b"sep": 9,
    b"oct": 10,
    b"nov": 11,
    b"dec": 12,
}


@dataclass(frozen=True)
class StationMetadata:
    station_code: str
    platform_code: str
    name: str
    latitude: float
    longitude: float
    altitude: float
    country_code: str
    subdivision: str
    land_use: str
    setting: str
    wmo_region: int
    gaw_type: str
    other_identifiers: str | None = None
    lab_code: str = "US06L"


SURFRAD_AOD_STATIONS: dict[str, StationMetadata] = {
    "bnd": StationMetadata(
        station_code="US0035R",
        platform_code="US0035S",
        name="Bondville, Illinois",
        latitude=40.0499992371,
        longitude=-88.3666687012,
        altitude=213.0,
        country_code="US",
        subdivision="IL",
        land_use="Agricultural",
        setting="Rural",
        wmo_region=4,
        gaw_type="R",
    ),
    "fpe": StationMetadata(
        station_code="US0905R",
        platform_code="US0905S",
        name="Fort Peck, Montana",
        latitude=48.310001,
        longitude=-105.099998,
        altitude=634.0,
        country_code="US",
        subdivision="MT",
        land_use="Agricultural",
        setting="Rural",
        wmo_region=4,
        gaw_type="R",
    ),
    "gcr": StationMetadata(
        station_code="US0900R",
        platform_code="US0900S",
        name="Goodwin Creek, Mississippi",
        latitude=34.25,
        longitude=-89.870003,
        altitude=98.0,
        country_code="US",
        subdivision="MS",
        land_use="Agricultural",
        setting="Rural",
        wmo_region=4,
        gaw_type="C",
    ),
    "bos": StationMetadata(
        station_code="US0901R",
        platform_code="US0901S",
        name="Boulder Table Mountain, Colorado",
        latitude=40.125,
        longitude=-105.2369995117,
        altitude=1689.0,
        country_code="US",
        subdivision="CO",
        land_use="Desert",
        setting="Mountain",
        wmo_region=4,
        gaw_type="C",
        other_identifiers="34 (BSRN), Boulder (AERONET)",
    ),
    "dra": StationMetadata(
        station_code="US0902R",
        platform_code="US0902S",
        name="Desert Rock, Nevada",
        latitude=36.619999,
        longitude=-116.017998,
        altitude=1007.0,
        country_code="US",
        subdivision="NV",
        land_use="Desert",
        setting="Rural",
        wmo_region=4,
        gaw_type="C",
    ),
    "psu": StationMetadata(
        station_code="US0903R",
        platform_code="US0903S",
        name="Penn State, Pennsylvania",
        latitude=40.783333,
        longitude=-77.933334,
        altitude=393.0,
        country_code="US",
        subdivision="PA",
        land_use="Forest",
        setting="Rural",
        wmo_region=4,
        gaw_type="R",
    ),
    "sxf": StationMetadata(
        station_code="US0904R",
        platform_code="US0904S",
        name="Sioux Falls, South Dakota",
        latitude=43.73,
        longitude=-96.620003,
        altitude=473.0,
        country_code="US",
        subdivision="SD",
        land_use="Remote Park",
        setting="Rural",
        wmo_region=4,
        gaw_type="R",
    ),
}


def normalize_station_code(station: str) -> str:
    """Normalize a SURFRAD AOD or GAW station code to the GAW code."""
    station = station.lower()
    for gaw_station, aod_station in AOD_STATION_REMAP.items():
        if station == aod_station:
            return gaw_station
    return station


def is_available() -> bool:
    """Return whether the EBAS NASA Ames writer dependencies are importable."""
    try:
        from ebas.io.file import nasa_ames  # noqa: F401
        from ebas.io.ebasmetadata import DatasetCharacteristicList  # noqa: F401
        from nilutility.datatypes import DataObject  # noqa: F401
    except ImportError:
        return False
    return True


def _require_ebas_dependencies():
    try:
        from ebas.io.ebasmetadata import DatasetCharacteristicList
        from ebas.io.file import nasa_ames
        from nilutility.datatypes import DataObject
    except ImportError as exc:
        raise RuntimeError(
            "EBAS file generation requires ebas-io and nilutility. "
            "Install EBAS-IO from https://git.nilu.no/ebas/ebas-io/-/tree/master."
        ) from exc
    return nasa_ames, DatasetCharacteristicList, DataObject


def _ebas_version() -> str:
    try:
        from ebas import __version__ as version
    except ImportError:
        return "unknown"
    return version


def _surfradpy_version() -> str:
    try:
        return _importlib_metadata.version("surfradpy")
    except _importlib_metadata.PackageNotFoundError:
        return "unknown"


def _organization(DataObject):
    return DataObject(
        OR_CODE="US06L",
        OR_NAME=(
            "National Oceanic and Atmospheric Administration/Earth System Research "
            "Laboratory/Global Monitoring Division"
        ),
        OR_ACRONYM="NOAA/ESRL/GMD",
        OR_UNIT=None,
        OR_ADDR_LINE1="325 Broadway",
        OR_ADDR_LINE2=None,
        OR_ADDR_ZIP="80305",
        OR_ADDR_CITY="Boulder, CO",
        OR_ADDR_COUNTRY="USA",
    )


def _aod_originator(DataObject):
    return DataObject(
        PS_LAST_NAME="Augustine",
        PS_FIRST_NAME="John",
        PS_EMAIL="John.A.Augustine@noaa.gov",
        PS_ORG_NAME=(
            "National Oceanic and Atmospheric Administration/Earth System Research "
            "Laboratory/Global Monitoring Division"
        ),
        PS_ORG_ACR="NOAA/ESRL/GMD",
        PS_ORG_UNIT=None,
        PS_ADDR_LINE1="325 Broadway",
        PS_ADDR_LINE2=None,
        PS_ADDR_ZIP="80305",
        PS_ADDR_CITY="Boulder, CO",
        PS_ADDR_COUNTRY="USA",
        PS_ORCID="0000-0002-6645-7404",
    )


def _projects(tags: set[str]) -> list[str]:
    if "nrt" in tags:
        return ["GAW-WDCA_NRT", "NOAA-ESRL_NRT"]
    return ["GAW-WDCA", "NOAA-ESRL"]


def _wdca_id(gaw_station: str, station: StationMetadata) -> str:
    return f"GAWA{station.country_code.upper()}{station.subdivision.upper()}{gaw_station.upper()}"


def set_ebas_station_metadata(nas, gaw_station: str, tags: set[str] | None = None) -> None:
    """Set SURFRAD station metadata on an EBAS NASA Ames object."""
    _, _, DataObject = _require_ebas_dependencies()
    tags = tags or set()
    try:
        station = SURFRAD_AOD_STATIONS[gaw_station]
    except KeyError as exc:
        raise ValueError(f"Unsupported SURFRAD AOD EBAS station: {gaw_station!r}") from exc

    nas.metadata.timezone = "UTC"
    nas.metadata.revdate = _dt.datetime.now(tz=_dt.timezone.utc)
    nas.metadata.revision = "1"
    nas.metadata.revdesc = (
        f"Version numbering not tracked, generated by surfradpy {_surfradpy_version()} "
        f"and EBAS-IO {_ebas_version()}"
    )

    originator = [_aod_originator(DataObject)]
    submitter = [_aod_originator(DataObject)]
    nas.metadata.org = _organization(DataObject)
    nas.metadata.projects = _projects(tags)
    nas.metadata.originator = originator
    nas.metadata.submitter = submitter
    nas.metadata.station_code = station.station_code
    nas.metadata.platform_code = station.platform_code
    nas.metadata.station_name = station.name
    nas.metadata.station_wdca_id = _wdca_id(gaw_station, station)
    nas.metadata.station_gaw_id = gaw_station.upper()
    nas.metadata.station_gaw_name = station.name
    nas.metadata.station_other_ids = station.other_identifiers
    nas.metadata.station_state_code = station.subdivision
    nas.metadata.station_landuse = station.land_use
    nas.metadata.station_setting = station.setting
    nas.metadata.station_gaw_type = station.gaw_type
    nas.metadata.station_wmo_region = station.wmo_region
    nas.metadata.station_latitude = station.latitude
    nas.metadata.station_longitude = station.longitude
    nas.metadata.station_altitude = station.altitude
    nas.metadata.mea_latitude = station.latitude
    nas.metadata.mea_longitude = station.longitude
    nas.metadata.mea_altitude = station.altitude
    nas.metadata.mea_height = None


def _bin_weighted_average(bin_start: np.ndarray, values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    values = np.asarray(values)
    weights = np.asarray(weights)
    assert (values.shape[0],) == (weights.shape[0],)

    weighted_values = (values.T * weights).T
    shaped_weights = np.full(values.T.shape, weights.T, dtype=weights.dtype).T
    invalid_values = np.invert(np.isfinite(weighted_values))
    weighted_values[invalid_values] = 0
    shaped_weights[invalid_values] = 0

    sum_values = np.add.reduceat(weighted_values, bin_start, dtype=np.float64)
    sum_weights = np.add.reduceat(shaped_weights, bin_start, dtype=np.float64)

    valid_averages = sum_weights != 0
    sum_values[valid_averages] /= sum_weights[valid_averages]
    sum_values[np.invert(valid_averages)] = nan

    empty_bins = np.where(bin_start[:-1] == bin_start[1:])[0]
    sum_values[empty_bins + 1] = nan

    return sum_values


def _bin_stddev(
    bin_start: np.ndarray,
    values: np.ndarray,
    unweighted_mean: np.ndarray,
    mask: np.ndarray | None = None,
) -> np.ndarray:
    sq = values**2
    sq[np.invert(np.isfinite(sq))] = 0
    if mask is not None:
        sq[mask] = 0
    sq = np.add.reduceat(sq, bin_start, dtype=np.float64)

    count = np.full(values.shape, 1, dtype=np.float64)
    count[np.invert(np.isfinite(values))] = 0
    if mask is not None:
        count[mask] = 0
    count = np.add.reduceat(count, bin_start, dtype=np.float64)

    result = np.full(unweighted_mean.shape, nan, dtype=np.float64)
    valid_values = count >= 2
    result[valid_values] = sq[valid_values] / count[valid_values] - unweighted_mean[valid_values] ** 2
    result[result < 0.0] = nan
    valid_values = np.isfinite(result)

    result[valid_values] = np.sqrt(result[valid_values] * (count[valid_values] / (count[valid_values] - 1)))
    result[np.invert(valid_values)] = nan
    return result


def _fixed_interval_bins(times: np.ndarray, interval: int | float) -> tuple[np.ndarray, np.ndarray]:
    bin_numbers = np.empty_like(times, dtype=np.int64)
    np.floor(times / interval, out=bin_numbers, casting="unsafe")
    return np.unique(bin_numbers, return_index=True)


def _fixed_interval_times(bin_numbers: np.ndarray, dtype, interval: int | float) -> np.ndarray:
    bin_times = np.empty_like(bin_numbers, dtype=dtype)
    return np.multiply(bin_numbers, interval, out=bin_times, casting="unsafe")


def _bin_quantiles(
    bin_start: np.ndarray,
    values: np.ndarray,
    quantiles: np.ndarray | float | list[float],
) -> np.ndarray:
    values = np.asarray(values)
    quantiles = np.asarray(quantiles).flatten()
    result = np.empty((bin_start.shape[0], *values.shape[1:], quantiles.shape[0]), dtype=values.dtype)

    for bin_number, start in enumerate(bin_start):
        stop = bin_start[bin_number + 1] if bin_number + 1 < bin_start.shape[0] else None
        with np.errstate(all="ignore"):
            result[bin_number] = np.nanquantile(values[start:stop], quantiles, axis=0).T
    return result


class FixedIntervalAverager:
    """Fixed interval averaging used by the legacy AOD EBAS converter."""

    def __init__(self, interval_ms: int, times_epoch_ms: np.ndarray):
        assert len(times_epoch_ms.shape) == 1
        assert times_epoch_ms.shape[0] > 0

        self._interval = interval_ms
        self._original_times = times_epoch_ms
        self._weights = np.full(times_epoch_ms.shape, 1.0, dtype=np.float64)
        self._bin_numbers, self._bin_start = _fixed_interval_bins(times_epoch_ms, interval_ms)
        assert len(self._bin_start) > 0

    @property
    def times(self) -> np.ndarray:
        return _fixed_interval_times(self._bin_numbers, np.int64, self._interval)

    def __call__(self, values: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
        weights = self._weights
        if mask is not None:
            weights = np.array(weights)
            weights[mask] = 0
        return _bin_weighted_average(self._bin_start, values, weights)

    def bitwise_or(self, values: np.ndarray) -> np.ndarray:
        return np.bitwise_or.reduceat(values, self._bin_start, dtype=values.dtype)

    def unweighted_mean(self, values: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
        weights = np.full((values.shape[0],), 1, dtype=values.dtype)
        if mask is not None:
            weights[mask] = 0
        return _bin_weighted_average(self._bin_start, values, weights)

    def stddev(
        self,
        values: np.ndarray,
        unweighted_mean: np.ndarray | None = None,
        mask: np.ndarray | None = None,
    ) -> np.ndarray:
        if unweighted_mean is None:
            unweighted_mean = self.unweighted_mean(values)
        return _bin_stddev(self._bin_start, values, unweighted_mean, mask=mask)

    def quantiles(self, values: np.ndarray, quantiles: np.ndarray | float | list[float]) -> np.ndarray:
        return _bin_quantiles(self._bin_start, values, quantiles)


class AodEbasData:
    """Parsed SURFRAD AOD data ready for EBAS NASA Ames export."""

    DATA_LINE = re.compile(rb"^\s*(\d{2})(\d{2})\s+([01])\s+(\d.+)")
    MVC_VALUE = re.compile(rb"-?9+(\.9+)?")

    def __init__(self, station: str):
        self.epoch_ms = np.empty((0,), dtype=np.int64)
        self.cloud_screening_failure = np.empty((0,), dtype=np.bool_)
        self.aod: list[np.ndarray] = []
        self.wavelengths: list[float] = []
        self._latest_wavelengths = 0
        self._ds: xr.Dataset | None = None
        self.station = station

    @property
    def tags(self) -> set[str]:
        return {"aod"}

    @property
    def ds(self) -> xr.Dataset:
        if self._ds is None:
            labelled_aod = self.aod[: len(self.wavelengths)]
            if labelled_aod:
                values = np.column_stack(labelled_aod)
            else:
                values = np.empty((self.epoch_ms.shape[0], 0), dtype=np.float64)
            datetime = self.epoch_ms.astype("datetime64[ms]")
            self._ds = xr.Dataset(
                data_vars={
                    "aod": (("datetime", "wavelength"), values),
                    "cloud_screening_failure": ("datetime", self.cloud_screening_failure),
                },
                coords={"datetime": datetime, "wavelength": self.wavelengths},
                attrs={"station": self.station, "source": "SURFRAD AOD"},
            )
            self._ds.aod.attrs["units"] = "no unit"
            self._ds.cloud_screening_failure.attrs["values"] = "False = passed; True = failed"
        return self._ds

    def _invalidate_dataset(self) -> None:
        self._ds = None

    @staticmethod
    def _file_start_from_header(line: bytes) -> int | None:
        match = re.match(rb"^\s*(\d{2})-([a-z]{3})-(\d{4})\s+", line.strip(), re.IGNORECASE)
        if not match:
            return None
        day = int(match.group(1))
        month = _MONTHS.get(match.group(2).lower())
        year = int(match.group(3))
        if 1970 <= year <= 2999 and month and 1 <= day <= 31:
            return int(_dt.datetime(year, month, day, tzinfo=_dt.timezone.utc).timestamp())
        return None

    @staticmethod
    def _file_start_from_name(file_name: str) -> int | None:
        match = re.fullmatch(r"[a-z]{3}_(\d{4})(\d{2})(\d{2})\.aod", file_name, re.IGNORECASE)
        if not match:
            return None
        year = int(match.group(1))
        month = int(match.group(2))
        day = int(match.group(3))
        if 1970 <= year <= 2999 and 1 <= month <= 12 and 1 <= day <= 31:
            return int(_dt.datetime(year, month, day, tzinfo=_dt.timezone.utc).timestamp())
        return None

    def integrate_file(self, source_file: BinaryIO, file_name: str, tz_offset: int) -> None:
        """Read one daily AOD text file into this object."""
        source_file.readline()

        file_start = self._file_start_from_header(source_file.readline())
        if file_start is None:
            file_start = self._file_start_from_name(file_name)
        if file_start is None:
            print(f"Unable to determine file start time for {file_name}, skipping")
            return

        wavelengths: list[int] = []
        for field in source_file.readline().strip().split():
            try:
                wavelength = float(field)
            except (TypeError, ValueError):
                break
            if wavelength <= 0.0:
                break
            wavelengths.append(int(round(wavelength)))
        if not wavelengths:
            print(f"No channel wavelengths in {file_name}, skipping")
            return

        if file_start > self._latest_wavelengths:
            while len(self.wavelengths) < len(wavelengths):
                self.wavelengths.append(0)
            for index, wavelength in enumerate(wavelengths):
                self.wavelengths[index] = wavelength
            self._latest_wavelengths = file_start

        source_file.readline()
        source_file.readline()
        source_file.readline()

        file_times: list[int] = []
        cloud_screening_failure: list[bool] = []
        wavelength_values: list[list[float]] = [[] for _ in range(len(wavelengths))]

        for line in source_file:
            match = self.DATA_LINE.match(line)
            if not match:
                continue
            hour = int(match.group(1))
            minute = int(match.group(2))
            bad = int(match.group(3)) == 1
            fields = match.group(4).strip().split()
            line_time = file_start + hour * 60 * 60 + minute * 60 + tz_offset * 60 * 60
            line_time = int(round(line_time * 1000))

            file_times.append(line_time)
            cloud_screening_failure.append(bad)
            for index in range(len(wavelengths)):
                if index < len(fields):
                    value = fields[index]
                    if self.MVC_VALUE.fullmatch(value):
                        parsed = nan
                    else:
                        try:
                            parsed = float(value)
                        except (TypeError, ValueError):
                            parsed = nan
                else:
                    parsed = nan
                wavelength_values[index].append(parsed)

        while len(wavelengths) > len(self.aod):
            self.aod.append(np.full(self.epoch_ms.shape, nan, dtype=np.float64))
        for index in range(len(wavelengths), len(self.aod)):
            self.aod[index] = np.append(
                self.aod[index],
                np.full((len(file_times),), nan, dtype=np.float64),
            )

        self.epoch_ms = np.append(self.epoch_ms, file_times)
        self.cloud_screening_failure = np.append(self.cloud_screening_failure, cloud_screening_failure)
        for index in range(len(wavelengths)):
            self.aod[index] = np.append(self.aod[index], wavelength_values[index])
        self._invalidate_dataset()

    def sort_data(self, year: int) -> None:
        """Keep data in the requested UTC year and sort by timestamp."""
        possible_start = int(
            floor(_dt.datetime(year, 1, 1, tzinfo=_dt.timezone.utc).timestamp() * 1000)
        )
        possible_end = int(
            ceil(_dt.datetime(year + 1, 1, 1, tzinfo=_dt.timezone.utc).timestamp() * 1000)
        )
        valid = np.all((self.epoch_ms >= possible_start, self.epoch_ms < possible_end), axis=0)
        self.epoch_ms = self.epoch_ms[valid]
        self.cloud_screening_failure = self.cloud_screening_failure[valid]
        for index in range(len(self.aod)):
            self.aod[index] = self.aod[index][valid]

        idx = np.argsort(self.epoch_ms)
        self.epoch_ms = self.epoch_ms[idx]
        self.cloud_screening_failure = self.cloud_screening_failure[idx]
        for index in range(len(self.aod)):
            self.aod[index] = self.aod[index][idx]
        self._invalidate_dataset()

    def _file_metadata(self, nas) -> None:
        set_ebas_station_metadata(nas, self.station, self.tags)
        station = SURFRAD_AOD_STATIONS[self.station]

        nas.metadata.instr_type = "rotating_shadowband_spectral_radiometer"
        nas.metadata.lab_code = station.lab_code
        nas.metadata.instr_name = f"YES_MFR-7_{self.station.upper()}"
        nas.metadata.instr_manufacturer = "YES"
        nas.metadata.instr_model = "MFR-7"
        nas.metadata.rescode_sample = "3mn"
        nas.metadata.method = "US06L_MFR_JA"
        nas.metadata.matrix = "aerosol"
        nas.metadata.unit = "no unit"
        nas.metadata.hum_temp_ctrl = "None"
        nas.metadata.hum_temp_ctrl_desc = "in situ measurement, ambient conditions"
        nas.metadata.detection_limit = (0.001, "no unit")
        nas.metadata.std_method = "Augustine2003"
        nas.metadata.acknowledgements = "Request acknowledgement details from data originator"
        nas.metadata.mea_height = 2

    @staticmethod
    def _set_file_times(nas, start_epoch_ms: np.ndarray, interval_ms: int) -> None:
        nas.sample_times = [
            (
                _dt.datetime.fromtimestamp(int(start_epoch_ms[index]) / 1000.0, tz=_dt.timezone.utc),
                _dt.datetime.fromtimestamp(
                    int(start_epoch_ms[index + 1]) / 1000.0
                    if index + 1 < start_epoch_ms.shape[0]
                    and start_epoch_ms[index + 1] - start_epoch_ms[index] <= interval_ms
                    else (int(start_epoch_ms[index]) + interval_ms) / 1000.0,
                    tz=_dt.timezone.utc,
                ),
            )
            for index in range(start_epoch_ms.shape[0])
        ]

        nas.metadata.period = "1y"
        nas.metadata.reference_date = _dt.datetime(nas.sample_times[0][1].year, 1, 1, tzinfo=_dt.timezone.utc)

    @staticmethod
    def _rounded_values(values: np.ndarray) -> list[float | None]:
        return [(float(round(float(value), 3)) if isfinite(value) else None) for value in values]

    @staticmethod
    def _variable_metadata(wavelength: float, statistics: str | None = None, prefix: str = "AOD"):
        _, DatasetCharacteristicList, DataObject = _require_ebas_dependencies()
        metadata = DataObject()
        metadata.comp_name = "aerosol_optical_depth"
        metadata.unit = "no unit"
        if statistics is not None:
            metadata.statistics = statistics
        metadata.title = f"{prefix}{int(round(wavelength))}"
        metadata.characteristics = DatasetCharacteristicList()
        metadata.characteristics.add_parse(
            "Wavelength",
            f"{int(round(wavelength))} nm",
            "rotating_shadowband_spectral_radiometer",
            "aerosol_optical_depth",
        )
        return metadata

    def write_level1(self, output_directory: str | Path) -> None:
        """Write an EBAS level 1 NASA Ames file."""
        nasa_ames, _, DataObject = _require_ebas_dependencies()
        ds = self.ds

        nas = nasa_ames.EbasNasaAmes()
        self._file_metadata(nas)
        self._set_file_times(nas, self.epoch_ms, 3 * 60 * 1000)

        nas.metadata.datalevel = "1"
        nas.metadata.statistics = "arithmetic mean"
        nas.metadata.duration = "3mn"
        nas.metadata.resolution = "3mn"

        flags = [[799] if bool(bad) else [] for bad in ds.cloud_screening_failure.to_numpy()]

        for wavelength_index, wavelength in enumerate(ds.wavelength.to_numpy()):
            nas.variables.append(
                DataObject(
                    values_=self._rounded_values(ds.aod.isel(wavelength=wavelength_index).to_numpy()),
                    flags=flags,
                    metadata=self._variable_metadata(wavelength),
                )
            )

        nas.write(createfiles=True, destdir=str(output_directory))

    def write_level2(self, year: int, output_directory: str | Path) -> None:
        """Write an EBAS level 2 NASA Ames file."""
        nasa_ames, _, DataObject = _require_ebas_dependencies()
        ds = self.ds
        average = FixedIntervalAverager(60 * 60 * 1000, self.epoch_ms)
        output_times = np.arange(
            int(floor(_dt.datetime(year, 1, 1, tzinfo=_dt.timezone.utc).timestamp() * 1000)),
            int(ceil(_dt.datetime(year + 1, 1, 1, tzinfo=_dt.timezone.utc).timestamp() * 1000)),
            60 * 60 * 1000,
            dtype=np.int64,
        )
        output_time_origin = output_times[0]
        target_indices = (average.times - output_time_origin) // (60 * 60 * 1000)

        nas = nasa_ames.EbasNasaAmes()
        self._file_metadata(nas)
        self._set_file_times(nas, output_times, 60 * 60 * 1000)

        nas.metadata.datalevel = "2"
        nas.metadata.duration = "1h"
        nas.metadata.resolution = "1h"
        nas.metadata.type = "TU"

        flags = np.full(output_times.shape, 0xFF, dtype=np.uint8)
        flags[target_indices] = average.bitwise_or(self.cloud_screening_failure)

        def convert_flag(flag: int) -> list[int]:
            if flag == 0xFF:
                return [999]
            if flag:
                return [799]
            return []

        ebas_flags = [convert_flag(int(flag)) for flag in flags]

        def filled(values: np.ndarray) -> np.ndarray:
            result = np.full(output_times.shape, nan, dtype=np.float64)
            result[target_indices] = values
            return result

        for wavelength_index, wavelength in enumerate(ds.wavelength.to_numpy()):
            values = ds.aod.isel(wavelength=wavelength_index).to_numpy()
            nas.variables.append(
                DataObject(
                    values_=self._rounded_values(filled(average(values))),
                    flags=ebas_flags,
                    metadata=self._variable_metadata(wavelength, statistics="arithmetic mean"),
                )
            )

            nas.variables.append(
                DataObject(
                    values_=self._rounded_values(filled(average.quantiles(values, 0.5)[:, 0])),
                    flags=ebas_flags,
                    metadata=self._variable_metadata(wavelength, statistics="median", prefix="AOD_med"),
                )
            )

            nas.variables.append(
                DataObject(
                    values_=self._rounded_values(filled(average.stddev(values))),
                    flags=ebas_flags,
                    metadata=self._variable_metadata(wavelength, statistics="stddev", prefix="AOD_std"),
                )
            )

        nas.write(createfiles=True, destdir=str(output_directory))


def load_station_year(source: str | Path, gaw_station: str, year: int) -> AodEbasData:
    """Load one SURFRAD AOD station-year from daily text files."""
    gaw_station = normalize_station_code(gaw_station)
    aod_station = AOD_STATION_REMAP.get(gaw_station, gaw_station)
    tz_offset = AOD_TIMEZONE.get(gaw_station, 0)
    aod_data = AodEbasData(gaw_station)
    file_root = Path(source)
    year_dir = file_root / aod_station / str(year)

    if year_dir.is_dir():
        for aod_file in year_dir.iterdir():
            if not aod_file.is_file():
                continue
            with aod_file.open("rb") as f:
                aod_data.integrate_file(f, aod_file.name, tz_offset)

    for aod_file in (
        file_root / aod_station / str(year - 1) / f"{aod_station}_{year - 1}1231.aod",
        file_root / aod_station / str(year + 1) / f"{aod_station}_{year + 1}0101.aod",
    ):
        try:
            with aod_file.open("rb") as f:
                aod_data.integrate_file(f, aod_file.name, tz_offset)
        except FileNotFoundError:
            continue

    aod_data.sort_data(year)
    return aod_data


def _gzip_file(path: Path) -> Path:
    output = path.with_name(path.name + ".gz")
    with path.open("rb") as source, gzip.open(output, "wb") as target:
        shutil.copyfileobj(source, target)
    path.unlink()
    return output


def upload_ebas(source_directory: str | Path, archive_output: str | Path) -> None:
    """Compress, upload, and archive generated EBAS files."""
    try:
        import paramiko
    except ImportError as exc:
        raise RuntimeError("EBAS upload requires paramiko.") from exc

    source_directory = Path(source_directory)
    archive_output = Path(archive_output)

    for input_file in source_directory.iterdir():
        if input_file.is_file() and input_file.suffix != ".gz":
            _gzip_file(input_file)

    print(f"Creating archive directory {archive_output} ...")
    archive_output.mkdir(parents=True, exist_ok=True)

    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        ssh.connect(hostname="ebas-submissions.nilu.no", username="ebasftp", timeout=120.0)
        sftp = ssh.open_sftp()
        try:
            for input_file in source_directory.iterdir():
                if not input_file.is_file():
                    continue
                print(f"Uploading {input_file.name} ...")
                sftp.put(str(input_file), input_file.name)
                shutil.move(str(input_file), archive_output / input_file.name)
                print("  ... Done")
        finally:
            sftp.close()
    finally:
        ssh.close()


def process_aod_ebas(
    year: int,
    stations: list[str] | tuple[str, ...] | None = None,
    output: str | Path | None = None,
    source: str | Path = DEFAULT_SOURCE,
    archive: str | Path = DEFAULT_ARCHIVE,
) -> None:
    """Generate EBAS AOD files for one year and optionally upload them."""
    if not 1971 <= year <= 2999:
        raise ValueError("year must be between 1971 and 2999")

    process_stations = [normalize_station_code(station) for station in (stations or DEFAULT_STATIONS)]
    unknown_stations = sorted(set(process_stations) - set(SURFRAD_AOD_STATIONS))
    if unknown_stations:
        raise ValueError(f"Unsupported station code(s): {', '.join(unknown_stations)}")

    for gaw_station in process_stations:
        aod_station = AOD_STATION_REMAP.get(gaw_station, gaw_station)
        aod_data = load_station_year(source, gaw_station, year)
        if aod_data.epoch_ms.shape[0] == 0:
            print(f"No data available for {gaw_station.upper()}/{aod_station.upper()} in {year}")
            continue
        print(
            f"Loaded {aod_data.epoch_ms.shape[0]} data points "
            f"for {gaw_station.upper()}/{aod_station.upper()} in {year}"
        )

        if output:
            output_dir = Path(output)
            output_dir.mkdir(parents=True, exist_ok=True)
            aod_data.write_level1(output_dir)
            aod_data.write_level2(year, output_dir)
        else:
            archive_dir = Path(archive)
            with TemporaryDirectory() as output_dir:
                aod_data.write_level1(output_dir)
                upload_ebas(Path(output_dir), archive_dir / gaw_station.lower() / "lev1" / str(year))
            with TemporaryDirectory() as output_dir:
                aod_data.write_level2(year, output_dir)
                upload_ebas(Path(output_dir), archive_dir / gaw_station.lower() / "lev2" / str(year))
