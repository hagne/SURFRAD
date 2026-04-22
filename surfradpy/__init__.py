"""surfradpy package."""

from importlib import import_module

__all__ = ["NCEI", "mfr_raw2netcdf"]


def __getattr__(name):
    if name == "NCEI":
        return import_module(".NCEI", __name__)
    if name == "mfr_raw2netcdf":
        return import_module(".products.mfr_raw2netcdf", __name__)
    raise AttributeError(f"module 'surfradpy' has no attribute {name!r}")
