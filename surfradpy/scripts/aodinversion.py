#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 16:17:01 2023

@author: hagen

TODO
=====
- implement exit if instance is running
- run in cronjob on pulsar4... weekly?
- address some of the warnings
"""
import productomator.lab as prolab
import surfradpy.aodinversion as sfraodinv
import pandas as pd

#### Functions that are executed by the scripts!    
def produce_aodinversion1_01():
    """
    This is the product based on the AOD1 which is John's AOD.

    Returns
    -------
    None.

    """
    reporter = prolab.Reporter('aodinversion1',
                             log_folder='/home/grad/htelg/.processlogs/',
                             reporting_frequency=(6, 'h'),
                            )
    
    import warnings
    warnings.filterwarnings('ignore',category=RuntimeWarning, message = 'invalid value encountered in log10')
    warnings.filterwarnings('ignore',category=RuntimeWarning, message = 'invalid value encountered in divide')
    warnings.filterwarnings('ignore',category=RuntimeWarning, message = 'divide by zero encountered in log10')
    
    #### FIXME: address warnings below!!! .... uncommend to see them, then fix them
    warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
    warnings.filterwarnings('ignore',category=pd.errors.FutureWarning)
    

    run = sfraodinv.AODInversion( version=1.0,
                                    channels=[415,500, 670, 870, 1625],
                                    sites = ['tbl'],
                                    start = '20160101',
                                    end = None, #'20230101',
                                    p2fldaod='/nfs/grad/surfrad/products_level2/aod_netcdf/v1.0/',
                                    p2fldout = '/nfs/grad/surfrad/products_level2/aodinversion/1.0/',
                                    ignore_in_cloud  = True,
                                    test = False,
                                    verbose = False,
                                    reporter = reporter)
    # run.workplan = run.workplan.sample(frac = 1)
    
    run.workplan = run.workplan[::-1]
    print(f'workplan.shape: {run.workplan.shape}')
    
    run.run_product(max_processes=10)
    
    reporter.wrapup()
    return

def produce_aodinversion1_01_catchup():
    """
    This is the product based on the AOD1 which is John's AOD.

    Returns
    -------
    None.

    """
    reporter = prolab.Reporter('aodinversion1_catchup',
                             log_folder='/home/grad/htelg/.processlogs/',
                             reporting_frequency=(6, 'h'),
                            )
    
    import warnings
    warnings.filterwarnings('ignore',category=RuntimeWarning, message = 'invalid value encountered in log10')
    warnings.filterwarnings('ignore',category=RuntimeWarning, message = 'invalid value encountered in divide')
    warnings.filterwarnings('ignore',category=RuntimeWarning, message = 'divide by zero encountered in log10')
    
    #### FIXME: address warnings below!!! .... uncommend to see them, then fix them
    warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
    warnings.filterwarnings('ignore',category=pd.errors.FutureWarning)
    

    run = sfraodinv.AODInversion( version=1.0,
                                    channels=[415,500, 670, 870, 1625],
                                    #sites = ['tbl'],
                                    sites = ['bon', 'dra', 'fpk', 'gwn', 'sxf', 'psu'],
                                    start = '20160101',
                                    end = None, #'20230101',
                                    p2fldaod='/nfs/grad/surfrad/products_level2/aod_netcdf/v1.0/',
                                    p2fldout = '/nfs/grad/surfrad/products_level2/aodinversion/1.0/',
                                    ignore_in_cloud  = True,
                                    test = False,
                                    verbose = False,
                                    reporter = reporter)
    run.workplan = run.workplan[::-1]
    # run.workplan = run.workplan.sample(frac = 1)
    print(f'workplan.shape: {run.workplan.shape}')
    
    run.run_product(max_processes=25)
    
    reporter.wrapup()
    return

def produce_aodinversion3_01():
    """This is the product based on AOD3, which is the realtime aod ... I think"""
    reporter = prolab.Reporter('aodinversion',
                             log_folder='/home/grad/htelg/.processlogs/',
                             reporting_frequency=(6, 'h'),
                            )
    run = sfraodinv.AODInversion(reporter = reporter)
    print(f'workplan.shape: {run.workplan.shape}')
    run.run_product()
