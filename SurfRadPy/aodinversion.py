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
    def __init__(self, 
                 reporter = None, 
                 version = 0.1, 
                 channels = [500, 670, 870, 1625],
                 sites = ['tbl',],
                 p2fldaod = '/export/htelg/data/grad/surfrad/aod_3/3.0/',
                 ignore_in_cloud = False,
                 test = False,
                 verbose = False,
                 ):
        self.version = version
        self.channels_used = channels
        self.p2fldaod = pl.Path(p2fldaod)
        self.p2fldout = pl.Path('/export/htelg/data/grad/surfrad/AODinversion/')
        self.sites = sites
        self.name_format = 'srf_aodinv_{site}_{year:04d}{month:02d}{day:02d}.nc'
        self.start = '20190101'
        self.end = '20200101'
        self.ignore_in_cloud = ignore_in_cloud
        self.reporter = reporter
        self.test = test
        self.verbose = verbose
        self._workplan = None
        
    @property
    def workplan(self):
        if isinstance(self._workplan, type(None)):
            
            #files = self.p2fldaod.glob('**/*')

            #remove folders
            #files = [f for f in files if f.is_file()]
            
            files = []
            p2fldsites = list(self.p2fldaod.glob('*'))
            for p2fs in p2fldsites:
                site = p2fs.name
                if not site in self.sites:
                    # print(f'{site} not in sites')
                    continue
                
                filest = list(p2fs.glob('*'))
                #remove folders
                filest = [f for f in filest if f.is_file()]
                files += filest
            
            df = pd.DataFrame(files, columns=['p2in',])
            
            df.index = df.apply(lambda row: pd.to_datetime(row.p2in.name.split('.')[0].split('_')[-1]), axis = 1)
            
            df.sort_index(inplace=True)
            
            df = df.truncate(self.start, self.end)
            
            # add site column
            df['site'] = df.apply(lambda row: row.p2in.parent.name, axis = 1)
            
            # pars output path
            df['p2out'] = df.apply(lambda row: self.p2fldout.joinpath(f'{self.version:0.1f}').joinpath(row.site).joinpath(f'{row.name.year:04d}').joinpath(self.name_format.format(site = row.site, year = row.name.year, month = row.name.month, day = row.name.day)), axis = 1)
            
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
        # unify the nominal channel centers
        ds.channel.values[ds.channel == 673] = 670
        # return ds        
                       
        #### aod instance 
        
        if verbose:
            print('AOD_AOT class - ansgstrom')
            
        if self.ignore_in_cloud:
            aodarr = ds.aod.where(ds.cloudmask == 0)
        else:
            aodarr = ds.aod
        aod = atmcop.AOD_AOT(aodarr)
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
                                verbose=self.verbose,
                                cut_off = 1e-10,
                                test=self.test,
                    )
        # ds = ds.merge(aod_inv.rename({var:f'aodinv_{var}' for var in aod_inv.variables if var not in aod_inv.coords}))
        # ds = ds.merge(aod_inv.rename({var:f'aodinv_{var}' for var in aod_inv.variables if var not in aod_inv.coords}))
        
        if 'sun' in ds:
            dsn = dsn.merge(ds.sun)
            dsn = dsn.rename({'sun': 'sun_position'})
        else:
            dsn = dsn.merge(ds.sun_position)
        
        if 'cloudmask' in ds:
            dsn = dsn.merge(ds.cloudmask)
        else:
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