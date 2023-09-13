#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 13:44:17 2023

@author: htelg
"""
import xarray as xr
import pathlib as pl
import atmPy.aerosols.physics.column_optical_properties as atmcop
import pandas as pd

class AODInversion(object):
    def __init__(self, reporter = None):
        self.version = 0.1
        self.channels_used = [500, 670, 870, 1625]
        self.p2fldaod = pl.Path('/export/htelg/data/grad/surfrad/aod_3/3.0/')
        self.p2fldout = pl.Path('/export/htelg/data/grad/surfrad/AODinversion/')
        self.site = 'tbl'
        self.name_format = 'srf_aodinv_{site}_{year:04d}{month:02d}{day:02d}.nc'
        self.start = '20190101'
        self.end = '20200101'
        
        self.reporter = reporter
        self._workplan = None
        
    @property
    def workplan(self):
        if isinstance(self._workplan, type(None)):
            
            files = self.p2fldaod.glob('**/*')

            #remove folders
            files = [f for f in files if f.is_file()]
            
            df = pd.DataFrame(files, columns=['p2in',])
            
            df.index = df.apply(lambda row: pd.to_datetime(row.p2in.name.split('.')[0].split('_')[-1]), axis = 1)
            
            df.sort_index(inplace=True)
            
            df = df.truncate(self.start, self.end)
            
            # pars output path
            df['p2out'] = df.apply(lambda row: self.p2fldout.joinpath(f'{self.version:0.1f}').joinpath(self.site).joinpath(f'{row.name.year:04d}').joinpath(self.name_format.format(site = self.site, year = row.name.year, month = row.name.month, day = row.name.day)), axis = 1)
            
            # remove if file exists
            df = df[~(df.apply(lambda row: row.p2out.is_file(), axis = 1))]
            
            self._workplan = df
        
        return self._workplan
    
    def run_product(self):
        for idx, row in self.workplan.iterrows():
            if not isinstance(self.reporter, type(None)):
                self.reporter.log()
            self.run_single_row(row)
    
    def run_single_row(self, row, verbose = False):
        ds = xr.open_dataset(row.p2in)
               
        #### aod instance 
        
        if verbose:
            print('AOD_AOT class - ansgstrom')
        aod = atmcop.AOD_AOT(ds.aod)
        # anglist = []
        # assert(False), 'aod is not defined anywhere!?!?' #something like: aodinst = atmcop.AOD_AOT(aod)
        # for angcomb in angcombos:
            # col1 = 415
            # col2 = 673
            # ang = aod.aod2angstrom_exponent(column_1=angcomb[0], column_2=angcomb[1],)
            # coordval = f'{angcomb[0]}_{angcomb[1]}'
            # anglist.append(ang.expand_dims({'ang_channels': [coordval]}))
        # ang = xr.concat(anglist, dim = 'ang_channels')
        # ds['angstrom_exp'] = ang.transpose('datetime', 'ang_channels')
        if verbose:
            print('aod inversion')
        
        
        dsn = aod.invertAOD(width_of_aerosol_mode=(0.2, 0.25),
                                channels=self.channels_used,
                                all_valid=True,
                                verbose=False,
                                cut_off = 1e-10,
                                test=False,
                    )
        # ds = ds.merge(aod_inv.rename({var:f'aodinv_{var}' for var in aod_inv.variables if var not in aod_inv.coords}))
        # ds = ds.merge(aod_inv.rename({var:f'aodinv_{var}' for var in aod_inv.variables if var not in aod_inv.coords}))
        
        
        dsn = dsn.merge(ds.sun_position)
        dsn = dsn.merge(ds.cloudmask_michalsky)
        
        dsn.attrs['cannels_used_in_inversion'] = self.channels_used
        
        dsn.attrs['path2input_file'] = row.p2in.as_posix()
        
        dsn.attrs['product_version'] = self.version
        dsn.attrs['production_date'] = f'{pd.Timestamp.now()}'
        
        #### save
        ### make sure all folders exist
        row.p2out.parent.parent.parent.parent.mkdir(exist_ok=True)
        row.p2out.parent.parent.parent.mkdir(exist_ok=True)
        row.p2out.parent.parent.mkdir(exist_ok=True)
        row.p2out.parent.mkdir(exist_ok=True)
        #### save
        dsn.to_netcdf(row.p2out)
        
        return dsn