#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2026-04-10

@author: hagen telg

This script performs spectral and cosine calibrations on SURFRAD MFRSR raw data (netcdf version).
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
        start = None,
        end = None,
        noofdays = 60,
        test = False,
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
        Start date/time. If not provided, computed as `end - noofdays`.
    end : str or pandas.Timestamp, optional
        End date/time. Defaults to current time.
    noofdays : int, optional
        Number of days to process when `start` is not given.
    test : bool, optional
        If True, process only one test row and stop.
    raise_errors : bool, optional
        If True, raise processing errors from the worker.
    verbose : bool, optional
        Print progress information.
    """
    import pandas as pd
    import productomator.lab as prolab
    import surfradpy.products.mfrsr_cosinecorrection as srfcc

    if verbose:
        print("start surfrad_mfrsr_cosinecalibration")
    out = {}
    reporter = prolab.Reporter('surfrad_mfrsr_cosinecalibraton', 
                                log_folder=log_folder,
                                reporting_frequency=(6, 'h'),
                            )

    if end is None:
         end = pd.Timestamp.now()
    if start is None:
         start = end - pd.to_timedelta(noofdays, 'D')

    sites = ['bnd',
            'dra',
            'gwn',
            'psu', 
            'sxf',
            'tbl',
            'fpe',
    ]
    version_in = '0.4'
    # version_out = '0.2'

    for site in sites:
        if verbose:
            print(site)
            print('-----')
        p2fld_in = f'{prefix}/grad/Inst/MFR/SURFRAD/{site}/mfrsr/raw.netcdf/v{version_in}/'
        p2fld_out = f'{prefix}/grad/Inst/MFR/SURFRAD/{site}/mfrsr/cosine_corrected/v{{version}}/'
        ci = srfcc.CalibrateMFRSR(p2fld_in  = p2fld_in,
                                        p2fld_out = p2fld_out,
                                        date_from_name = lambda name: pd.to_datetime(name.split('_')[-1].replace('.nc', '')),
                                        output_file_format = f'{{year}}/{site}_mfrsr_cosinecorrected_{{date}}.nc',
                                        # glob_pattern_raw = "*.xmd",
                                        start = start,
                                        end = end,
                                        # site = site,
                                        # version = '0.3', # this is actually 0.3.1 but i don't want to rerun everything.
                                        # path2surfrad_database = db_path,
                                        reporter = reporter,
                                        verbose = verbose,
                                )
        print(f'{site} workplan.shape: {ci.workplan.shape}')
        if test:
            last_processed = ci.process_row(iloc = 1, save=False)
            break
        else:
            # try:

            last_processed = ci.process(raise_errors = raise_errors)
            # except:
            #     return ci
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
    parser.add_argument('--prefix', default='/nfs/grad/')
    parser.add_argument('--log-folder', default='/home/grad/htelg/.processlogs/')
    parser.add_argument('--start', default=None)
    parser.add_argument('--end', default=None)
    parser.add_argument('--noofdays', type=int, default=60)
    parser.add_argument('--test', action='store_true')
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
        noofdays=args.noofdays,
        test=args.test,
        raise_errors=args.raise_errors,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
