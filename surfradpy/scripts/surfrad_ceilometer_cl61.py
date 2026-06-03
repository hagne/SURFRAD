#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2026-06-02

@author: hagen telg
"""
import argparse
import inspect

import surfradpy.config as sfp_config


class _RawDefaultsHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    pass


def run(prefix = '/nfs', 
        db_path=None, 
        log_folder='/home/grad/htelg/.processlogs/',
        noofdays = 60,
        start = None,
        end = None,
        test = False,
        raise_errors = False,
        verbose = True,
):
    """Run SURFRAD CL61 ceilometer cloud product level 1 processing.

    Parameters
    ----------
    prefix : str, optional
        The prefix for the input and output paths. Default is '/nfs'.
    db_path : str, optional
        Path to the SURFRAD database.
    log_folder : str, optional
        Folder to save the process logs. Default is '/home/grad/htelg/.processlogs/'.
    noofdays : int, optional
        Number of days to process when `start` is not given. Default is 60.
    start : str or pandas.Timestamp, optional
        Start date/time. If not provided, computed as `end - noofdays`.
    end : str or pandas.Timestamp, optional
        End date/time. Defaults to current time.
    test : int or bool, optional
        If 1, return the product instance without processing.
        If 2, process the first workplan row without saving.
        If 3, process the first workplan row and save output.
        `False` or 0 processes the full workplan.
    raise_errors : bool, optional
        If True, raise processing errors from the worker.
    verbose : bool, optional
        Print progress information.
    """
    import pandas as pd
    import productomator.lab as prolab
    import surfradpy.database as srfdb
    import surfradpy.products.cl61_cloud_prod_lev1 as srfcl61

    
    out = {}
    if db_path is None:
        db_path = sfp_config.get_db_path()
    if db_path is None:
        cfg_path = sfp_config.get_config_path()
        raise FileNotFoundError(
            "No database path provided and no default configured. "
            f"Set SURFRAD_DB_PATH or add [database] path to {cfg_path}."
        )
    
    db = srfdb.SurfradDatabase(db_path)


    reporter = prolab.Reporter('surfrad_cl61_cloud_prod_lev1_0.1', 
                                log_folder=log_folder,
                                reporting_frequency=(6, 'h'),
                                
                            )    
    
    site = 'inl'

    site_info = db.find_site_info(site)

    if end is None:
        end = pd.Timestamp.now()
    else:
        end = pd.to_datetime(end)
    if start is None: 
        start = end - pd.to_timedelta(noofdays, 'D')
    else:
        start = pd.to_datetime(start)


    ci = srfcl61.Cl61CloudLevel1_v0_1(
        p2fld_in =  f'{prefix}/grad/Inst/Ceil/SURFRAD/{site.upper()}/CL61_Daily/',
        p2fld_out = f'{prefix}/grad/surfrad/products_level1/cl61_cloud_prod_lev1/v{{version}}/{site}/',
        date_from_name = lambda name: pd.to_datetime(name.split('_')[-1].split('.')[0]),
        output_file_format = f'{{year}}/{site}.cl61.cloud_prod{{date}}.nc',
        site_code = site_info.abb,
        site_name = site_info.name,
        lat = site_info.latitude,
        lon = site_info.longitude,
        start = start,
        end = end,
        glob_pattern='*.nc',
        verbose = verbose,
        reporter = reporter,
        )
    
    print(f'{site} workplan.shape: {ci.workplan.shape}')
    if test == 1:
        print('Testing:  returning the product instance without running any processing.')
        return ci
    if test == 2:
        print('Testing:  running first row in workplan without saving the output file.')
        last_processed = ci.process_row(iloc = 1, save=False)
    elif test == 3:
        print('Testing:  running first row in workplan including saving the output file.')
        last_processed = ci.process_row(iloc = 1, save=True)
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
        choices=(0, 1, 2, 3),
        default=0,
        help=(
            'Test mode: 0=full run, 1=workplan only, '
            '2=process first row without saving, 3=process first row and save.'
        ),
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
        noofdays=args.noofdays,
        start=args.start,
        end=args.end,
        test=args.test,
        raise_errors=args.raise_errors,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
