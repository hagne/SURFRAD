#!/usr/bin/env python3
"""Download GEOS-FP TO3 data and interpolate it to SURFRAD sites."""

import argparse
import inspect

import surfradpy.config as sfp_config


class _RawDefaultsHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    pass


def run(
    prefix="/nfs",
    db_path=None,
    path_out=None,
    download_dir=None,
    force_ipv4=False,
    log_folder="/home/grad/htelg/.processlogs/",
    start=None,
    end=None,
    noofdays=60,
    test=0,
    raise_errors=False,
    verbose=True,
):
    """Generate daily GEOS-FP TO3 files for all SURFRAD sites.

    Parameters
    ----------
    prefix : str, optional
        Filesystem prefix used for the default output path.
    db_path : str, optional
        Path to the SURFRAD database.
    path_out : str, optional
        Output folder. May contain a ``{version}`` placeholder.
    download_dir : str, optional
        Cache folder used by the atmPy GEOS-FP downloader.
    force_ipv4 : bool, optional
        Route NCCS OPeNDAP requests over IPv4.
    log_folder : str, optional
        Folder where process logs are written.
    start : str or pandas.Timestamp, optional
        First UTC day. If omitted, computed as ``end - noofdays``.
    end : str or pandas.Timestamp, optional
        Last UTC day. Defaults to the current day.
    noofdays : int, optional
        Number of days subtracted from ``end`` when ``start`` is omitted.
    test : int or bool, optional
        If 1, only create the workplan. If 2, process its first row without
        saving. If 3, process and save its first row. Zero processes all rows.
    raise_errors : bool, optional
        If True, raise processing errors from the product.
    verbose : bool, optional
        Print progress information.
    """
    import pandas as pd
    import productomator.lab as prolab
    import surfradpy.products.geosfp as sfp_geosfp

    if db_path is None:
        db_path = sfp_config.get_db_path()
    if db_path is None:
        cfg_path = sfp_config.get_config_path()
        raise FileNotFoundError(
            "No database path provided and no default configured. "
            f"Set SURFRAD_DB_PATH or add [database] path to {cfg_path}."
        )

    if end is None:
        end = pd.Timestamp.now()
    else:
        end = pd.to_datetime(end)
    if start is None:
        start = end - pd.to_timedelta(noofdays, "D")
    else:
        start = pd.to_datetime(start)

    if path_out is None:
        path_out = (
            f"{prefix}/stu3data2/Model_data/geos_fp/"
            "tavg1_2d_slv_Nx/surfrad/v{version}/"
        )

    reporter = prolab.Reporter(
        "surfrad_geosfp_0.1",
        log_folder=log_folder,
        reporting_frequency=(6, "h"),
    )
    product = sfp_geosfp.GeosFpSurfrad(
        path_out=path_out,
        start=start,
        end=end,
        path2surfrad_database=db_path,
        download_dir=download_dir,
        force_ipv4=force_ipv4,
        reporter=reporter,
        verbose=verbose,
    )

    if verbose:
        print("start surfrad_geosfp")
        print(f"workplan.shape: {product.workplan.shape}")

    last_processed = None
    if test == 1:
        pass
    elif test == 2:
        last_processed = product.process_row(iloc=0, save=False)
    elif test == 3:
        last_processed = product.process_row(iloc=0, save=True)
    else:
        last_processed = product.process(raise_errors=raise_errors)

    reporter.wrapup()
    return {
        "product_instance": product,
        "last_processed": last_processed,
        "reporter": reporter,
    }


def _build_parser():
    parser = argparse.ArgumentParser(
        description=inspect.getdoc(run) or "",
        formatter_class=_RawDefaultsHelpFormatter,
    )
    parser.add_argument("--prefix", default="/nfs")
    parser.add_argument("--db-path", default=None)
    parser.add_argument("--path-out", default=None)
    parser.add_argument("--download-dir", default=None)
    parser.add_argument("--force-ipv4", action="store_true")
    parser.add_argument(
        "--log-folder",
        default="/home/grad/htelg/.processlogs/",
    )
    parser.add_argument("--start", default=None)
    parser.add_argument("--end", default=None)
    parser.add_argument("--noofdays", type=int, default=60)
    parser.add_argument(
        "--test",
        type=int,
        choices=(0, 1, 2, 3),
        default=0,
        help=(
            "Test mode: 0=full run, 1=workplan only, "
            "2=process first row without saving, 3=process first row and save."
        ),
    )
    parser.add_argument("--raise-errors", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose")
    parser.add_argument("--no-verbose", action="store_false", dest="verbose")
    parser.set_defaults(verbose=True)
    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)
    return run(
        prefix=args.prefix,
        db_path=args.db_path,
        path_out=args.path_out,
        download_dir=args.download_dir,
        force_ipv4=args.force_ipv4,
        log_folder=args.log_folder,
        start=args.start,
        end=args.end,
        noofdays=args.noofdays,
        test=args.test,
        raise_errors=args.raise_errors,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
