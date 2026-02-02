#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2026-01-15

@author: hagen telg

This script converts SURFRAD MFR raw (.xmd) data files to netcdf format. In the process it merges and truncates the data to daily files in UTC time.
"""

import surfradpy.mfr_raw2netcdf as mfr_r2nc
import surfradpy.database as sfp_db
import surfradpy.config as sfp_config
import productomator.lab as prolab
import pandas as pd


def run(prefix = '/nfs/grad/', 
        db_path=None, 
        log_folder='/home/grad/htelg/.processlogs/',
        test = False,
        verbose = True,):
    
    if db_path is None:
        db_path = sfp_config.get_db_path()
    if db_path is None:
        cfg_path = sfp_config.get_config_path()
        raise FileNotFoundError(
            "No database path provided and no default configured. "
            f"Set SURFRAD_DB_PATH or add [database] path to {cfg_path}."
        )

    reporter = prolab.Reporter('surfrad_mfr_raw2netcdf', 
                                log_folder=log_folder,
                                reporting_frequency=(6, 'h'),
                                
                            )

    sites = ['bnd',
            'dra',
            'gwn',
            'psu', 
            'sxf',
            'tbl',
            'fpe',
    ]
    end = pd.Timestamp.now()
    start = end - pd.to_timedelta(60, 'D')
    for site in sites:
        print(site)
        path_in = f'{prefix}/Inst/MFR/SURFRAD/{site}/mfrsr/raw/'
        ci = mfr_r2nc.MfrsrRawToNetcdf(path_in,
                                    f'{prefix}/Inst/MFR/SURFRAD/{site}/mfrsr/raw.netcdf/v{{version}}/',
                                    '{year}/{site}_mfrsr_raw_{year}{month}{day}.nc',
                                    glob_pattern_raw='*.xmd',
                                    start = start,
                                    end = end,
                                    # verbose = True,
                                    site = site,
                                    path2surfrad_database  = db_path,
                                    reporter = reporter,
                                    verbose = verbose,
                                )
        print(f'{site} workplan.shape: {ci.workplan.shape}')
        if test:
            out = ci.process(verbose=verbose, save=False, justone = False)
            break
        else:
            out = ci.process(verbose=verbose)

    reporter.wrapup()
    out['reporter'] = reporter
    return out


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', default='/nfs/grad/')
    parser.add_argument('--db-path', default=None)
    parser.add_argument('--log-folder', default='/home/grad/htelg/.processlogs/')
    parser.add_argument('--test', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose')
    parser.add_argument('--no-verbose', action='store_false', dest='verbose')
    parser.set_defaults(verbose=True)
    args = parser.parse_args(argv)

    return run(
        prefix=args.prefix,
        db_path=args.db_path,
        log_folder=args.log_folder,
        test=args.test,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
