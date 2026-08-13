"""Daily GEOS-FP total-column ozone interpolated to SURFRAD sites."""

import math
import pathlib as pl
import socket

import pandas as pd
import productomator.worker as prowo
import xarray as xr

import surfradpy.database as sfp_db


GEOS_FP_COLLECTION = "tavg1_2d_slv_Nx"
GEOS_FP_VARIABLES = ("to3",)

SURFRAD_SITE_CODES = (
    "inl",
    "bnd",
    "dra",
    "gwn",
    "psu",
    "sxf",
    "tbl",
    "fpe",
)

GEOS_FP_LATITUDE_RESOLUTION = 0.25
GEOS_FP_LONGITUDE_RESOLUTION = 0.25
GEOS_FP_LATITUDE_ORIGIN = -90.0
GEOS_FP_LONGITUDE_ORIGIN = -180.0


def load_surfrad_sites(path2surfrad_database=None):
    """Return coordinates and metadata for all operational SURFRAD sites."""
    database = sfp_db.SurfradDatabase(path2surfrad_database)
    rows = [
        database.find_site_info(abb=site)
        for site in SURFRAD_SITE_CODES
    ]

    invalid_sites = [
        row.abb
        for row in rows
        if row.network.lower() != "surfrad"
    ]
    if invalid_sites:
        raise ValueError(
            "Sites are not members of the SURFRAD network: "
            f"{', '.join(invalid_sites)}"
        )

    sites = xr.Dataset(
        coords={
            "site": list(SURFRAD_SITE_CODES),
            "site_id": ("site", [row.site_id for row in rows]),
            "site_name": ("site", [row["name"] for row in rows]),
            "lat": ("site", [row.latitude for row in rows]),
            "lon": ("site", [row.longitude for row in rows]),
            "elevation": ("site", [row.elevation for row in rows]),
        },
        attrs={"network": "SURFRAD"},
    )
    sites.site.attrs.update(
        {
            "long_name": "SURFRAD site abbreviation",
            "cf_role": "timeseries_id",
        }
    )
    sites.site_id.attrs["long_name"] = "SURFRAD database site identifier"
    sites.site_name.attrs["long_name"] = "SURFRAD site name"
    sites.lat.attrs.update(
        {
            "long_name": "site latitude",
            "standard_name": "latitude",
            "units": "degrees_north",
        }
    )
    sites.lon.attrs.update(
        {
            "long_name": "site longitude",
            "standard_name": "longitude",
            "units": "degrees_east",
        }
    )
    sites.elevation.attrs.update(
        {
            "long_name": "site elevation above sea level",
            "standard_name": "height_above_mean_sea_level",
            "units": "m",
        }
    )
    return sites


def _grid_floor(value, *, origin, resolution):
    return origin + math.floor((value - origin) / resolution) * resolution


def _grid_ceil(value, *, origin, resolution):
    return origin + math.ceil((value - origin) / resolution) * resolution


def get_spatial_subset(sites):
    """Return the smallest GEOS-FP grid-aligned bbox bracketing the sites."""
    west = _grid_floor(
        sites.lon.min().item(),
        origin=GEOS_FP_LONGITUDE_ORIGIN,
        resolution=GEOS_FP_LONGITUDE_RESOLUTION,
    )
    east = _grid_ceil(
        sites.lon.max().item(),
        origin=GEOS_FP_LONGITUDE_ORIGIN,
        resolution=GEOS_FP_LONGITUDE_RESOLUTION,
    )
    south = _grid_floor(
        sites.lat.min().item(),
        origin=GEOS_FP_LATITUDE_ORIGIN,
        resolution=GEOS_FP_LATITUDE_RESOLUTION,
    )
    north = _grid_ceil(
        sites.lat.max().item(),
        origin=GEOS_FP_LATITUDE_ORIGIN,
        resolution=GEOS_FP_LATITUDE_RESOLUTION,
    )
    return west, south, east, north


def interpolate_to_sites(ds, sites):
    """Linearly interpolate a gridded GEOS-FP dataset to paired sites."""
    missing_coordinates = [
        coordinate
        for coordinate in ("datetime", "lat", "lon")
        if coordinate not in ds.coords
    ]
    if missing_coordinates:
        raise ValueError(
            "GEOS-FP dataset is missing coordinates: "
            f"{', '.join(missing_coordinates)}"
        )

    dsout = ds.interp(
        lat=sites.lat,
        lon=sites.lon,
        method="linear",
    )
    for coordinate in sites.coords:
        if coordinate != "site":
            dsout = dsout.assign_coords({coordinate: sites.coords[coordinate]})
        dsout.coords[coordinate].attrs = sites.coords[coordinate].attrs.copy()
    return dsout


class GeosFpSurfrad(prowo.Workplanner):
    """Create daily GEOS-FP TO3 files interpolated to SURFRAD sites."""

    version = "0.1"
    collection = GEOS_FP_COLLECTION

    def __init__(
        self,
        path_out,
        start,
        end,
        path2surfrad_database=None,
        download_dir=None,
        force_ipv4=False,
        reporter=None,
        verbose=True,
    ):
        self.start = pd.to_datetime(start).normalize()
        self.end = pd.to_datetime(end).normalize()
        if self.end < self.start:
            raise ValueError("end must be on or after start")

        super().__init__(
            p2fld_in=None,
            p2fld_out=str(path_out),
            output_file_format=(
                f"{{year}}/surfrad_geosfp_{self.collection}_{{date}}.nc"
            ),
            start=self.start,
            end=self.end,
            file_complete_check=False,
            reporter=reporter,
            verbose=verbose,
            version=self.version,
        )

        self.download_dir = download_dir
        self.force_ipv4 = force_ipv4
        self.variables = GEOS_FP_VARIABLES
        self.sites = load_surfrad_sites(path2surfrad_database)
        self.spatial_subset = get_spatial_subset(self.sites)

        self._geos_fp = None
        self.tp_row = None
        self.tp_dsin = None
        self.tp_dsout = None

    @property
    def geos_fp(self):
        if self._geos_fp is None:
            import atmPy.data_archives.nasa.geos_fp as atm_geos_fp

            self._geos_fp = atm_geos_fp.GeosFp(
                collection=self.collection,
                stream='seamless',
                download_dir=self.download_dir,
                chunks="auto",
                force_ipv4=self.force_ipv4,
            )
        return self._geos_fp

    def process_row(self, row=None, iloc=None, loc=None, save=True):
        """Download, interpolate, and optionally save one workplan day."""
        if iloc is not None:
            row = self.workplan.iloc[iloc]
        elif loc is not None:
            row = self.workplan.loc[pd.to_datetime(loc).normalize()]
        if row is None:
            raise ValueError("Provide row, iloc, or loc")
        self.tp_row = row

        day = row.name.normalize()
        dsin = self.geos_fp.download(
            day.date(),
            day.date(),
            variables=list(self.variables),
            bbox=self.spatial_subset,
        )
        self.tp_dsin = dsin

        try:
            day_end = day + pd.Timedelta(days=1)
            ds = dsin[list(self.variables)].sel(
                datetime=(dsin.datetime >= day) & (dsin.datetime < day_end)
            )
            # if ds.sizes.get("datetime") != 24:
            #     day_complete = "False"
            # else:
            #     day_complete = "True"
            #     # raise ValueError(
            #     #     f"Expected 24 hourly GEOS-FP samples for "
            #     #     f"{day:%Y-%m-%d}, "
            #     #     f"found {ds.sizes.get('datetime', 0)}"
            #     # )

            dsout = interpolate_to_sites(ds, self.sites)
            dsout.attrs.update(
                {
                    "title": (
                        "GEOS-FP total-column ozone linearly interpolated to "
                        "SURFRAD sites"
                    ),
                    "processing_date": pd.Timestamp.now(tz="UTC").isoformat(),
                    "processing_server": socket.gethostname(),
                    "processing_class": (
                        f"{self.__class__.__module__}."
                        f"{self.__class__.__qualname__}"
                    ),
                    "product_version": self.version,
                    "interpolation_method": "linear latitude-longitude",
                    "spatial_subset_bbox": ", ".join(
                        str(value) for value in self.spatial_subset
                    ),
                    "spatial_subset_bbox_order": "west, south, east, north",
                    "site_count": self.sites.sizes["site"],
                    # "day_complete": day_complete,
                }
            )
            if self.geos_fp.paths:
                dsout.attrs["parent_files"] = ", ".join(
                    pl.Path(path).as_posix()
                    for path in self.geos_fp.paths
                )

            # The projected daily product is small (24 x 8 x 1), so detach it
            # from the downloaded source before closing that file handle.
            dsout.load()
        finally:
            dsin.close()

        self.tp_dsout = dsout
        if save:
            row.p2f_out.parent.mkdir(parents=True, exist_ok=True)
            temporary_path = row.p2f_out.with_name(
                f".{row.p2f_out.name}.tmp"
            )
            dsout.to_netcdf(temporary_path)
            temporary_path.replace(row.p2f_out)
        return {"dsout": dsout, "p2f_out": row.p2f_out}

    # def process(self, raise_errors=False):
    #     """Process all missing daily files and return the last result."""
    #     last_processed = None
    #     for date, row in self.workplan.iterrows():
    #         try:
    #             last_processed = self.process_row(row=row)
    #             if self.reporter is not None:
    #                 self.reporter.clean_increment()
    #         except Exception as error:
    #             if self.reporter is not None:
    #                 self.reporter.errors_increment()
    #             if raise_errors:
    #                 raise
    #             print(f"Error occurred while processing {date:%Y-%m-%d}: {error}")
    #             continue

    #         if self.verbose:
    #             print(".", end="", flush=True)
    #     if self.verbose and last_processed is not None:
    #         print()
    #     return last_processed
