#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2024

@author: hagen

TODO
=====

"""
import productomator.lab as prolab
import surfradpy.radflux2netcdf as srfr2nc

import pandas as pd
import warnings

#### Functions that are executed by the scripts!    
def main():
    """
    This is the product that tries to converd radflux to netcdf format ... daily.

    Returns
    -------
    None.

    """
    reporter = prolab.Reporter('radflux2netcdf',
                             log_folder='/home/grad/htelg/.processlogs/',
                             reporting_frequency=(6, 'h'),
                            )

    warnings.filterwarnings('ignore')
    # warnings.filterwarnings('ignore',category=RuntimeWarning, message = 'invalid value encountered in log10')
    # warnings.filterwarnings('ignore',category=RuntimeWarning, message = 'invalid value encountered in divide')
    # warnings.filterwarnings('ignore',category=RuntimeWarning, message = 'divide by zero encountered in log10')
    
    #### FIXME: address warnings below!!! .... uncommend to see them, then fix them
    # warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
    # warnings.filterwarnings('ignore',category=FutureWarning)
    
    rfi = srfr2nc.Convert(path2fld_in='/nfs/iftp/aftp/g-rad/surfrad/RadFlux/',    
                          path2fld_out='/nfs/grad/surfrad/products_level4/radflux/v{version}/',    
                          sites=['tbl', 'dra', 'fpk', 'gwn', 'psu', 'sxf', 'bon'],    
                          start='180 days',    
                          overwrite=False,
                          reporter = reporter)
    # run = sfraodinv.AODInversion( version=1.0,
    #                                 channels=[415,500, 670, 870, 1625],
    #                                 sites = ['tbl'],
    #                                 start = '20160101',
    #                                 end = None, #'20230101',
    #                                 p2fldaod='/nfs/grad/surfrad/products_level2/aod_netcdf/v1.0/',
    #                                 p2fldout = '/nfs/grad/surfrad/products_level2/aodinversion/1.0/',
    #                                 ignore_in_cloud  = True,
    #                                 test = False,
    #                                 verbose = False,
    #                                 reporter = reporter)
    # run.workplan = run.workplan.sample(frac = 1)
    
    rfi.workplan = rfi.workplan[::-1]
    print(f'workplan.shape: {rfi.workplan.shape}')
    
    rfi.process(error_handling = 'skip')
    
    reporter.wrapup()
    return reporter