#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2026-09-22

@author: hagen telg

This peforms the Radflux clearsky analysis for all permanent Surfrad sites. It populates the database but does not generate the radflux files, which will be generated in a separate step.
todo: upcate the below
Requirements:
- pvlib
- xarray
- pandas
- netcdf4  
- atmpy
- sklearn
- productomator
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
        days = 60,
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
    import surfradpy.database as srfdb


    if verbose:
        print("start surfrad_mfrsr_cosinecalibration")
    out = {}
    reporter = prolab.Reporter('surfrad_radflux_clearsky', 
                                log_folder=log_folder,
                                reporting_frequency=(6, 'h'),
                            )


    sites = ['inl',
             'bon',
            'dra',
            'gwn',
            'psu', 
            'sxf',
            'tbl',
            'fpk',
    ]
    version_in = '0.4'
    # version_out = '0.2'

    for site in sites:
        if verbose:
            print(site)
            print('-----')

        db = srfdb.SurfradDatabase(srfdb.get_default_db_path())
        site_info = db.find_site_info(abb = site)
        p2fld_in = f'/{prefix}/grad/surfrad/products_level1/radiation_netcdf/v1.1/{site}'
        radflux_parameters_db = f'/{prefix}/grad/surfrad/products_level2/radflux_rd/radflux_params_{site}_{{version}}.db'
        path2raflux_setting = f'/{prefix}/grad/surfrad/products_level2/radflux_rd/radflux_settings_{site}.toml'
        ci = srfrf.RadfluxClearskyParameterAnalysis(
            p2fld_in=p2fld_in,
            p2fld_out=None,
            database=None,
            file_name_format='*{date:%Y%m%d}*',
            output_file_format=None,
            start=start,
            end=end,
            days = days,
            input_directory_structure='yearly',
            output_directory_structure=None,
            file_complete_check=False,
            reporter=reporter,
            verbose=verbose,
            radflux_parameters_db = radflux_parameters_db,
            path2raflux_setting = path2raflux_setting,
            site = site_info
        )
        print(f'{site} workplan.shape: {ci.workplan.shape}')
        if test == 1:
            last_processed = None
            break
        elif test == 2:
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
    parser.add_argument('--prefix', default='/nfs')
    parser.add_argument('--log-folder', default='/home/grad/htelg/.processlogs/')
    parser.add_argument('--start', default=None)
    parser.add_argument('--end', default=None)
    parser.add_argument('--days', type=int, default=60)
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
