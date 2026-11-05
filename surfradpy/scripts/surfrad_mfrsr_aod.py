#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Produce SURFRAD MFRSR aerosol optical depth products for one or all sites.

Required packages are Python 3.11+, NumPy, pandas, xarray, netCDF4, SciPy,
pvlib, Ibis, ``atmPy``, ``productomator``, and ``surfradpy``. Create the Conda
environment with::

    conda create -n surfrad_mfrsr_aod -c conda-forge 'python>=3.13' numpy pandas xarray netcdf4 scipy pvlib ibis-framework pip

Install the local ``atmPy``, ``productomator``, and ``surfradpy`` projects into
that environment with pip.
"""

import argparse
import inspect
# import pandas as pd
import productomator.lab as prolab
import surfradpy.products.aod as srfaod


class _RawDefaultsHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    pass


def run(prefix = '/nfs',
        log_folder='/home/grad/htelg/.processlogs/',
        start = None,
        end = None,
        days = 90,
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

    #todo: remove when non-realime is implemented
    real_time = True
    if verbose:
        print("start surfrad_mfrsr_aod")
    out = {}
    logname = 'surfrad_mfrsr_aod'
    if real_time:
        logname += '_real_time'

    reporter = prolab.Reporter(logname, 
                                log_folder=log_folder,
                                reporting_frequency=(6, 'h'),
                            )

    if site is not None:
        sites = [site]
    else:
        sites = [
             'bnd',
            'dra',
            'gwn', 
            'sxf',
            'tbl',
            'fpe',
            # 'inl',
            'psu',
    ]
    def site_alt(site):
        if site == 'bnd':
            alt = 'bon'
        elif site =='fpe':
            alt = 'fpk'
        else:
            alt = site
        return alt

    for site in sites:
        if verbose:
            print(site)
            print('-----')
        #todo: this try/except should not be necessary, fix in worker.
        if 1:

            if real_time:
                subfld = 'near_real_time'
            else:
                subfld = 'final'

            p2fld_in = f'{prefix}/grad/Inst/MFR/SURFRAD/{site}/mfrsr/cosine_corrected/v0.2/'
            p2fld_out = f'{prefix}/grad/surfrad/products_level2/aod/v{{version}}/{subfld}/{site}'
            p2fld_model = f'{prefix}/stu3data2/Model_data/geos_fp/tavg1_2d_slv_Nx/surfrad/v0.1/'
            p2fdb = f'{prefix}/grad/surfrad/products_level2/aod/v0_prediction_{site}_v0.1_nearrealtime.db'
            sup_model = {'name': 'model', 
                                'path2folder': p2fld_model, 
                                'file_name_format': '*{date:%Y%m%d}*',
                                'input_directory_structure':'yearly'}
            sup_met = {'name': 'met',
                    'path2folder': f'/Volumes/grad/surfrad/scaled/{site_alt(site)}/',# 2026/gradobs-scaled.{site}.20261019.dat',
                    'file_name_format': '*{date:%Y%m%d}*',
                    'input_directory_structure':'yearly'}
            supplemental_data = [sup_model, sup_met]
            ci = srfaod.AOD(p2fld_in=p2fld_in,
                    supplemental_data=supplemental_data,
                    p2fld_out=p2fld_out,
                    path2v0database = p2fdb,
                        output_file_format='{site}_mfrsr_aod_{date}.nc',
                        start=start,
                        end=end,
                        days=days,
                    verbose = verbose,
                    site = site,
                    reporter = reporter)




            
            print(f'{site} workplan.shape: {ci.workplan.shape}')
            if test == 1:
                last_processed = None
                break
            elif test == 2:
                last_processed = ci.process_row(iloc = 0, save=False)
                break

            elif test == 3:
                last_processed = ci.process_row(iloc = 0, save=True)
                break
            else:
                # try:

                last_processed = ci.process(raise_errors = raise_errors)
                # except:
                #     return ci
        # except:
        #     continue
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
    parser.add_argument('--start', default=None)
    parser.add_argument('--end', default=None)
    parser.add_argument('--days', type=int, default=60)
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
