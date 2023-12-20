#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 10:31:27 2023

@author: hagen

Converts John's AOD (AOD1) to netcdf.
"""

import pathlib as pl
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.surfrad.surfrad as srf
import xarray as xr
import numpy as np
import pandas as pd


class Aod2Netcdf(object):
    def __init__(self, 
                overwrite = False,
                site = 'psu', #'all'
                            # ['tbl', 'bon',
                         # 'dra', 'fpk', 'gwn', 'psu', 'sxf', 'tbl'
                        # ]
                start_date = '2019-01-01',
                end_date = '2024-01-01',
                version = '1.0',
                path2basefld_in = '/nfs/grad/surfrad/aod/',
                path2basefld_out = '/export/htelg/data/grad/surfrad/aod/netcdf/v{version}/',
                ):
        self.overwrite = overwrite
        self.site = site #'all'
                    # ['tbl', 'bon',
                 # 'dra', 'fpk', 'gwn', 'psu', 'sxf', 'tbl'
                # ]
        self.start_date = start_date
        self.end_date = end_date
        self.version = version
        self.path2basefld_in = pl.Path(path2basefld_in)
        self.path2basefld_out = pl.Path('/export/htelg/data/grad/surfrad/aod/netcdf/v{version}/'.format(version = version))
        self.path2basefld_out.mkdir(exist_ok=True)

        path2current = self.path2basefld_out.parent.joinpath('current')
        if path2current.resolve() != self.path2basefld_out:
            path2current.unlink()
            path2current.symlink_to(self.path2basefld_out)
            print('current was linked to new folder')

        # private varibles
        self._workplan = None
    
    def process(self, if_error = 'raise',
                verbose = False,
                reporter = None):
        """
        

        Parameters
        ----------
        if_error : str, optional
            What to do if an error is encountered. The default is 'raise'.
            raise: raise the error
            skip: file is skipped, the error will be forwarded to the reporter.
        reporter: productomator.Reporter instance
        verbose : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """
        # print(df.shape)
        save = True
        # verbose = False
        for idx,row in self.workplan.iterrows():
            try:
                if verbose:
                    print(idx, end = ', ')
                else:
                    print('|', end = '')
                if row.path2file_out.is_file():
                    if not self.overwrite:
                        if not isinstance(reporter, type(None)):
                            reporter.warnings_increment()
                        continue
            
                out = self.process_single(row.path2file_in, row.site, self.version)
                ds = out['ds']
            
                if not row.path2file_out.parent.is_dir():
                    row.path2file_out.parent.mkdir()
            
                if save:
                    ds.to_netcdf(row.path2file_out)
                if not isinstance(reporter, type(None)):      
                    reporter.clean_increment()
            except:
                if if_error == 'raise':
                    raise
                elif if_error == 'skip':
                    if not isinstance(reporter, type(None)):      
                        reporter.errors_increment()
                    continue
            
            if not isinstance(reporter, type(None)):      
                reporter.log()
                
        if not isinstance(reporter, type(None)): 
            # aod is rarely process, in order to not have it kicked of the graphic 
            # i generate a warning if nothing was processed.
            if reporter.clean == 0:
                reporter.warnings_increment()
            reporter.log(reset_counters=False, overwrite_reporting_frequency=True)
        
        print('Done')  
    
    def process_single(self, path2file_in, site, version):
        aod = srf.open_path(path2file_in,
                            product = 'None',
                            cloud_sceened = False,
                            local2UTC=True,)
        if site == 'bon':
            site = 'bnd'
        site = srf.network.stations.find_site(site)
        
        ds = xr.Dataset()
        aod.AOD.data.columns.name = 'channel'
        aod.AOD.data.index.name = 'datetime'
        ds['aod'] = aod.AOD.data.tz_localize(None)#.astype(np.float32)
    
        aod.ang_exp.data.index.name = 'datetime'
    
        ds['ang'] = aod.ang_exp.data.tz_localize(None).iloc[:,0]
    
        ds['badmask'] = aod.aerosolmask.badmask_angstrom
        ds['badmaskpad'] = aod.aerosolmask.badmask_angstrom_padded
    
        ds['minutes2bad_angstrom'] = aod.aerosolmask.minutes2bad_angstrom
        aod.cloudmask.cloudmask_nativ.data.index.name = 'datetime'
        ds['cloudmask'] = aod.cloudmask.cloudmask_nativ.data.tz_localize(None).iloc[:,0]
    
        ds['aod'] = ds.aod.astype(np.float32)
    
        ds['ang'] = ds.ang.astype(np.float32)
        ds['badmask'] = ds.badmask.astype(int)
        ds['badmaskpad'] = ds.badmaskpad.astype(int)
        ds['minutes2bad_angstrom'] = ds.minutes2bad_angstrom.astype(int)
        ds['cloudmask'] = ds.cloudmask.astype(int)
        
        #### Sun info
        df_sun = site.get_sun_position(aod.AOD.data.tz_localize(None))
        df_sun.columns.name = 'sun_param'
        df_sun.drop('ampm', axis = 1, inplace=True)
        df_sun.rename({'altitude': 'elevation'}, axis=1, inplace=True)
        
        ds['sun'] = df_sun.astype(np.float32)
        ####
        
        ds.attrs['version'] = version
        ds.cloudmask.attrs['values'] = '0 = good; 1 = cloud detected'
        out = {}
        out['ds'] = ds
        # out['aod'] = aod
        # out['site'] = site
        # out['df_sun'] = df_sun
        return out

    @property
    def workplan(self):
        if isinstance(self._workplan, type(None)):
            files = []
            
            for path2sitefld_in in self.path2basefld_in.glob('*'):
            #     site = path2sitefld_in.name
                for path2yearfld_in in path2sitefld_in.glob('*'):
                    # year = path2yearfld_in.name
                    files += list(path2yearfld_in.glob('*'))
            
            files = [f for f in files if f.is_file()]
            
            df = pd.DataFrame({'path2file_in':files})
            
            df['site'] = df.apply(lambda row: row.path2file_in.parent.parent.name, axis = 1)
            
            def fn2datetime(row): 
                try:
                    dt =  pd.to_datetime(row.path2file_in.name.split('_')[1].split('.')[0])
                except:
                    print(row.path2file_in.name)
                    print(row.path2file_in)
                    raise
                return dt
            # df['datetime'] = df.apply(lambda row: pd.to_datetime(row.path2file_in.name.split('_')[1].split('.')[0]), axis = 1)
            df['datetime'] = df.apply(fn2datetime, axis = 1)
            
            
            
            df['path2file_out'] = df.apply(lambda row: self.path2basefld_out.joinpath(row.site).joinpath(row.datetime.date().__str__().replace('-', '')+'.nc'), axis = 1)
            df.index = df.datetime
            # dft = df
            
            # select site 
            if self.site != 'all':
                df = df[df.site == self.site].copy()
                
            # sleect time range
            df.sort_index(inplace=True)
            df = df.truncate(before=self.start_date, after = self.end_date)
            
            # remove if file exists
            if not self.overwrite:
                df = df[~df.apply(lambda row: row.path2file_out.is_file(), axis = 1)]
            
            self._workplan = df
            
        return self._workplan