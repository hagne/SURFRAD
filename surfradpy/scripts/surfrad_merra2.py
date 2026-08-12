#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2026-07-24

@author: hagen telg

This script downloads daily MERRA-2 M2T1NXSLV data and interpolates it to all
SURFRAD sites.
"""

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
    settings_path=None,
    auth_strategy="auto",
    log_folder="/home/grad/htelg/.processlogs/",
    start=None,
    end=None,
    noofdays=60,
    test=0,
    raise_errors=False,
    verbose=True,
):
    """Generate daily MERRA-2 M2T1NXSLV files for all SURFRAD sites.

    Parameters
    ----------
    prefix : str, optional
        Filesystem prefix used for the default output path.
    db_path : str, optional
        Path to the SURFRAD database.
    path_out : str, optional
        Output folder. May contain a ``{version}`` placeholder.
    download_dir : str, optional
        Cache folder used by the atmPy MERRA-2 downloader.
    settings_path : str, optional
        Path to an atmPy NASA Earthdata settings file.
    auth_strategy : str, optional
        atmPy Earthdata authentication strategy.
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
    import surfradpy.products.merra2 as sfp_merra2

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
        path_out = f"{prefix}/stu3data2/Model_data/merra_2/M2T1NXSLV/surfrad/v{{version}}/"
        
    reporter = prolab.Reporter(
        "surfrad_merra2_0.1",
        log_folder=log_folder,
        reporting_frequency=(6, "h"),
    )
    product = sfp_merra2.Merra2Surfrad(
        path_out=path_out,
        start=start,
        end=end,
        path2surfrad_database=db_path,
        download_dir=download_dir,
        settings_path=settings_path,
        auth_strategy=auth_strategy,
        reporter=reporter,
        verbose=verbose,
    )

    if verbose:
        print("start surfrad_merra2")
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
    parser.add_argument("--settings-path", default=None)
    parser.add_argument(
        "--auth-strategy",
        choices=("auto", "settings", "netrc", "environment", "interactive"),
        default="auto",
    )
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
        settings_path=args.settings_path,
        auth_strategy=args.auth_strategy,
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
