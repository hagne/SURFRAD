#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2026-09-22

@author: hagen telg

This script uses the results from the Radflux clearsky analysis for all permanent Surfrad sites and creates the actual radflux product.
todo: upcate the below
Requirements:
- pvlib
- xarray
- pandas
- netcdf4  
- atmpy
# - scikit-learn not sure if this is needed
- productomator
- dask
"""

import argparse
import inspect


class _RawDefaultsHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    pass


def run(prefix = '/nfs',
        log_folder='/home/grad/htelg/.processlogs/',
        start = '2025-09-03',
        end = None,
        days = None,
        site = None,
        real_time = False,
        test = 0,
        raise_errors = False,
        verbose = True,):
    """Run SURFRAD MFRSR spectral and cosine calibration across all sites.

    Parameters
    ----------
    prefix : str, optional
        Filesystem prefix used to build input/output paths.
    log_folder : str, optional
        Folder where process logs are written.
    start : str or pandas.Timestamp, optional
        Start date/time. If not provided, computed as `end - days`.
    end : str or pandas.Timestamp, optional
        End date/time. Defaults to current time.
    site : str, optional
        Run only for this site. If not provided, runs for all sites.
    days : int, optional
        Number of days to process when `start` is not given.
    real_time : bool, optional
        If True, run in real-time mode.
    test : bool, optional
        If 1, creates workplan at first site and returns si
        If 2, process  one test of the first site row and stop.
    raise_errors : bool, optional
        If True, raise processing errors from the worker.
    verbose : bool, optional
        Print progress information.
    """
    import pandas as pd
    import productomator.lab as prolab
    import surfradpy.products.radflux as srfrf
    # import surfradpy.database as srfdb
    import atmPy.general.measurement_site as atmms


    if verbose:
        print("start surfrad_mfrsr_cosinecalibration")
    out = {}
    if real_time:
        logname = 'clamps3_radflux_real_time'
        subfld = 'near_real_time'
    else:
        logname = 'clamps3_radflux'                
        subfld = 'final'

    reporter = prolab.Reporter(logname, 
                                log_folder=log_folder,
                                reporting_frequency=(6, 'h'),
                            )

    p2fld_in = f'{prefix}/grad/gradobs/raw/short_term/clamps3/radsys/netcdf/radiation/v1.1/'
    p2fld_out = f'{prefix}/grad/gradobs/raw/short_term/clamps3/radsys/netcdf/radiation/{{version}}/{subfld}/'
    site = atmms.Station(lat=39.94713, lon = -105.19635, alt = 1742.04, name = 'clamps3 at Marshall', abbreviation='clamps3', state = 'CO')

    radflux_parameters_db = f'{prefix}/grad/gradobs/raw/short_term/clamps3/radsys/netcdf/radiation/radflux_params_1.0.db'
    path2raflux_setting = f'{prefix}/grad/gradobs/raw/short_term/clamps3/radsys/netcdf/radiation/radflux_settings_1.0.toml'

    ci = srfrf.RadfluxClamps3(
        p2fld_in=p2fld_in,
        p2fld_out=p2fld_out,
        # database=None,
        file_name_format='*{date:%Y%m%d}*',
        output_file_format='{site}_radflux_{date}.nc',
        start=start,
        end=end,
        days = days,
        real_time=real_time,
        input_directory_structure='yearly',
        # output_directory_structure=None,
        # file_complete_check=True,
        reporter=None,
        verbose=True,
        radflux_parameters_db = radflux_parameters_db,
        path2raflux_setting = path2raflux_setting,
        site = site
    )

    print(f'{site} workplan.shape: {ci.workplan.shape}')
    if test == 1:
        last_processed = None
    elif test == 2:
        last_processed = ci.process_row(iloc = 1, save=False)
    else:
        # try:

        last_processed = ci.process(raise_errors = raise_errors)

    out['product_instance'] = ci
    out['last_processed'] = last_processed
    reporter.wrapup(print_dagster_report=True)
    out['reporter'] = reporter
    if test > 0:
        return out
    return 


def _build_parser():
    parser = argparse.ArgumentParser(
        description=inspect.getdoc(run) or "",
        formatter_class=_RawDefaultsHelpFormatter,
    )
    parser.add_argument('--prefix', default='/nfs')
    parser.add_argument('--log-folder', default='/home/grad/htelg/.processlogs/')
    parser.add_argument('--start', default='2025-09-03')
    parser.add_argument('--end', default=None)
    parser.add_argument('--days', type=int, default=None)
    parser.add_argument('--real-time', action='store_true', dest='real_time')
    parser.add_argument('--site', default=None)
    parser.add_argument('--test', type=int, default=0)
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
        log_folder=args.log_folder,
        start=args.start,
        end=args.end,
        days=args.days,
        test=args.test,
        raise_errors=args.raise_errors,
        site=args.site,
        real_time=args.real_time,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
