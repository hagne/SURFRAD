#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 16:17:01 2023

@author: hagen
"""
import productomator.lab as prolab
import SurfRadPy.aodinversion as sfraodinv

#### Functions that are executed by the scripts!    
def produce_aodinversion1_01():
    """
    This is the product based on the AOD1 which is John's AOD.

    Returns
    -------
    None.

    """
    import warnings
    warnings.filterwarnings('ignore',category=RuntimeWarning, message = 'invalid value encountered in log10')
    warnings.filterwarnings('ignore',category=RuntimeWarning, message = 'invalid value encountered in divide')
    warnings.filterwarnings('ignore',category=RuntimeWarning, message = 'divide by zero encountered in log10')
    
    reporter = prolab.Reporter('aodinversion1',
                             log_folder='/home/grad/htelg/.processlogs/',
                             reporting_frequency=(6, 'h'),
                            )
    run = sfraodinv.AODInversion( version=1.0,
                                    channels=[415,500, 670, 870, 1625],
                                    sites = ['tbl'],
                                    start = '20180101',
                                    end = '20220101',
                                    p2fldaod='/home/grad/htelg/data/grad/surfrad/aod1/v1.0',
                                    p2fldaod='',
                                    p2fldout = '/home/grad/htelg/data/grad/surfrad/aod1/v1.0',
                                    ignore_in_cloud  = True,
                                    test = False,
                                    verbose = False,
                                    reporter = reporter)
    print(f'workplan.shape: {run.workplan.shape}')
    run.run_product()
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

def run():
    produce_aodinversion1_01()