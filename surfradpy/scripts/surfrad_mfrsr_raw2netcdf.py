#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2026-01-15

@author: hagen telg

This script converts SURFRAD MFR raw (.xmd) data files to netcdf format. In the process it merges and truncates the data to daily files in UTC time.
"""

import surfradpy.config as sfp_config
import argparse
import inspect


class _RawDefaultsHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    pass


def run(prefix = '/nfs', 
        db_path=None, 
        log_folder='/home/grad/htelg/.processlogs/',
        start = None,
        noofdays = 60,
        end = None,
        test = False,
        raise_errors = False,
        verbose = True,):
    """Run SURFRAD MFRSR raw to netCDF conversion for all operational sites.

    Parameters
    ----------
    prefix : str, optional
        The prefix for the input and output paths. Default is '/nfs'.
    db_path : str, optional
        Path to the SURFRAD database. 
    log_folder : str, optional
        Folder to save the process logs. Default is '/home/grad/htelg/.processlogs/'.
    start : str or pandas.Timestamp, optional
        Start date/time. If not provided, computed as `end - noofdays`.
    noofdays : int, optional
        Number of days to process when `start` is not given. Default is 60.
    end : str or pandas.Timestamp, optional
        End date/time. Defaults to current time.
    test : int or bool, optional
        If 1, stop after generating one site workplan.
        If 2, only the first row of the workplan is processed. 
        `False` or 0 processes full workplans.
    raise_errors : bool, optional
        If True, raise processing errors from the worker.
    verbose : bool, optional
        Print progress information.
    """
    import pandas as pd
    import productomator.lab as prolab
    import surfradpy.products.mfr_raw2netcdf as mfr_r2nc


    out = {}
    if db_path is None:
        db_path = sfp_config.get_db_path()
    if db_path is None:
        cfg_path = sfp_config.get_config_path()
        raise FileNotFoundError(
            "No database path provided and no default configured. "
            f"Set SURFRAD_DB_PATH or add [database] path to {cfg_path}."
        )

    reporter = prolab.Reporter('surfrad_mfr_raw2netcdf0.4', 
                                log_folder=log_folder,
                                reporting_frequency=(6, 'h'),
                                
                            )

    sites = ['inl',
             'bnd',
            'dra',
            'gwn',
            'psu', 
            'sxf',
            'tbl',
            'fpe',
    ]
    if end is None:
        end = pd.Timestamp.now()
    else:
        end = pd.to_datetime(end)
    if start is None: 
        start = end - pd.to_timedelta(noofdays, 'D')
    else:
        start = pd.to_datetime(start)

    for site in sites:
        print(site)
        ci = mfr_r2nc.MfrsrRawToNetcdf(path_in = f'{prefix}/grad/Inst/MFR/SURFRAD/{{site}}/mfrsr/raw/',
                                        path_out = f'{prefix}/grad/Inst/MFR/SURFRAD/{{site}}/mfrsr/raw.netcdf/v{{version}}',
                                        # date_from_name = lambda name: pd.to_datetime(' '.join(name.split('.')[0].split('_')[-2:])),
                                        # name_pattern_netcdf = '{year}/{site}_mfrsr_raw_{date}.nc',
                                        # glob_pattern_raw = "*.xmd",
                                        start = start,
                                        end = end,
                                        site = site,
                                        # version = '0.4', # this is actually 0.3.1 but i don't want to rerun everything.
                                        path2surfrad_database = db_path,
                                        reporter = reporter,
                                        verbose = verbose,
                                )
        print(f'{site} workplan.shape: {ci.workplan.shape}')
        if test == 1:
            return ci
        if test == 2:
            last_processed = ci.process_row(iloc = 1, save=False)
            break
        else:
            last_processed = ci.process(raise_errors = raise_errors)
    out['product_instance'] = ci
    out['last_processed'] = last_processed
    reporter.wrapup()
    out['reporter'] = reporter
    return out


def _build_parser():
    parser = argparse.ArgumentParser(
        description=inspect.getdoc(run) or "",
        formatter_class=_RawDefaultsHelpFormatter,
    )
    parser.add_argument('--prefix', default='/nfs')
    parser.add_argument('--db-path', default=None)
    parser.add_argument('--log-folder', default='/home/grad/htelg/.processlogs/')
    parser.add_argument('--start', default=None)
    parser.add_argument('--end', default=None)
    parser.add_argument('--noofdays', type=int, default=60)
    parser.add_argument(
        '--test',
        type=int,
        choices=(0, 1, 2),
        default=0,
        help='Test mode: 0=full run, 1=workplan only, 2=process first row only.',
    )
    parser.add_argument('--raise-errors', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose')
    parser.add_argument('--no-verbose', action='store_false', dest='verbose')
    parser.set_defaults(verbose=True)
    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    return run(
        prefix=args.prefix,
        db_path=args.db_path,
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
