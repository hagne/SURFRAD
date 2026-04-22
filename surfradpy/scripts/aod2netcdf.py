# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 16:17:01 2023

@author: hagen

This module is the backend of the aod2netcdf script. 

Objectives:
    - convert John's AOD product to a unified netcdf product

"""
import warnings
import subprocess
import argparse
import inspect

def run(verbose=True):
    """Convert SURFRAD AOD product files to unified netCDF output."""
    import pandas as pd
    import productomator.lab as prolab
    import surfradpy.products.aod2netcdf as aod2nc

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
    if verbose:
        print(a2n.workplan)
    nooffiles = a2n.workplan.shape[0]
    if verbose:
        print(f'no 2 be processed: {nooffiles}')
    
    a2n.process(reporter=reporter, if_error='skip')
    
    #### rsync
    if verbose:
        print(f'errors: {reporter.errors}')
    # if False:
    if (nooffiles - reporter.errors)> 0:
        if verbose:
            print('starting rsync', end = '...')
        subprocess.run(
                        ["rsync", "-av", f"{a2n.path2basefld_out}/", 
                                          f"/nfs/iftp/aftp/g-rad/surfrad/aod_netcdf/v{a2n.version}/"],
                        check=True
                        )
        if verbose:
            print('done')
    else:
        if verbose:
            print('no files processed => resync skpped.')
        
        
        
    reporter.wrapup()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=inspect.getdoc(run) or "",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose')
    parser.add_argument('--no-verbose', action='store_false', dest='verbose')
    parser.set_defaults(verbose=True)
    args = parser.parse_args(argv)
    return run(verbose=args.verbose)
    
if __name__ == "__main__":
    main()
