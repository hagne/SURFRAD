# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 16:17:01 2023

@author: hagen

This module is the backend of the aod2netcdf script. 

Objectives:
    - convert John's AOD product to a unified netcdf product

"""
import productomator.lab as prolab
import surfradpy.aod2netcdf as aod2nc
import pandas as pd

def run():
    starttime = pd.Timestamp.now()
    print(f'start time: {starttime}')
    a2n = aod2nc.Aod2Netcdf(site = 'all', 
                            path2basefld_in = '/nfs/grad/surfrad/aod/',
                            path2basefld_out='/nfs/grad/surfrad/products_level2/aod_netcdf/v{version}/',
                            # overwrite=True
                           )
    a2n.workplan
    
    print(f'no 2 be processed: {a2n.workplan.shape[0]}')
    
    reporter = prolab.Reporter('aod2netcdf', 
                               # log_folder='/export/htelg/tmp/', 
                               # reporting_frequency=(1,'min'),
                              )
    
    a2n.process(reporter = reporter)
    
    print(f'clean: {reporter.clean}')
    print(f'warnings: {reporter.warnings}')
    print(f'errors: {reporter.errors}')
    
    endtime = pd.Timestamp.now()
    print(f'end time: {endtime}')
    print(f'total execution time: {endtime - starttime}')
        
    
if __name__ == "__main__":
    run()
