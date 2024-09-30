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
import warnings

def run():
    reporter = prolab.Reporter('aod2netcdf', 
                               # log_folder='/export/htelg/tmp/', 
                               # reporting_frequency=(1,'min'),
                              )
    #### FIXME: address warnings below!!! .... uncommend to see them, then fix them
    warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
    
    a2n = aod2nc.Aod2Netcdf(site = 'all', 
                            path2basefld_in = '/nfs/grad/surfrad/aod/',
                            path2basefld_out='/nfs/grad/surfrad/products_level2/aod_netcdf/v{version}/',
                            # overwrite=True
                           )
    
    max2process = 3000
    a2n.workplan = a2n.workplan.iloc[-max2process:]
    
    print(f'no 2 be processed: {a2n.workplan.shape[0]}')
    
    a2n.process(reporter=reporter, if_error='skip')
    
    reporter.wrapup()
        
    
if __name__ == "__main__":
    run()
