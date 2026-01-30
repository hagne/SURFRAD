#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2026-01-15

@author: hagen telg

This script converts SURFRAD MFR raw (.xmd) data files to netcdf format. In the process it merges and truncates the data to daily files in UTC time.
"""

import surfradpy.mfr_raw2netcdf as mfr_r2nc
import surfradpy.database as sfp_db
import productomator.lab as prolab
import pandas as pd


def run(prefix = '/nfs/grad/', 
        db_path='/home/grad/htelg/prog/SURFRAD/notebooks/databases/surfrad_database.db', 
        log_folder='/home/grad/htelg/.processlogs/',
        test = False,):

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
                                )
        print(f'{site} workplan.shape: {ci.workplan.shape}')
        if test:
            out = ci.process(verbose=True, save=False, justone = False)
            break
        else:
            out = ci.process(verbose=True)

    reporter.wrapup()
    out['reporter'] = reporter
    return out