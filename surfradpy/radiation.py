#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 13:03:05 2023

@author: hagen
"""
import pandas as pd
import numpy as np

import pathlib as pl
# import ipywidgets as widgets
# from IPython.display import display
import xarray as xr

def read_surfrad(p2f):
    # https://gml.noaa.gov/aftp/data/radiation/surfrad/dra/README
    # make the column names
    collab = """year			integer	year, i.e., 1995
    jday			integer	Julian day (1 through 365 [or 366])
    month			integer	number of the month (1-12)
    day			integer	day of the month(1-31)
    hour			integer	hour of the day (0-23)
    min			integer	minute of the hour (0-59)
    dt			real	decimal time (hour.decimalminutes, e.g., 23.5 = 2330)
    zen			real	solar zenith angle (degrees)
    dw_solar		real	downwelling global solar (Watts m^-2)
    uw_solar		real	upwelling global solar (Watts m^-2)
    direct_n		real	direct-normal solar (Watts m^-2)
    diffuse		real	downwelling diffuse solar (Watts m^-2)
    dw_ir			real	downwelling thermal infrared (Watts m^-2)
    dw_casetemp		real	downwelling IR case temp. (K)
    dw_dometemp		real	downwelling IR dome temp. (K)
    uw_ir			real	upwelling thermal infrared (Watts m^-2)
    uw_casetemp		real	upwelling IR case temp. (K)
    uw_dometemp		real	upwelling IR dome temp. (K)
    uvb			real	global UVB (milliWatts m^-2)
    par			real	photosynthetically active radiation (Watts m^-2)
    netsolar		real	net solar (dw_solar - uw_solar) (Watts m^-2)
    netir			real	net infrared (dw_ir - uw_ir) (Watts m^-2)
    totalnet		real	net radiation (netsolar+netir) (Watts m^-2)
    temp			real	10-meter air temperature (?C)
    rh			real	relative humidity (%)
    windspd		real	wind speed (ms^-1)
    winddir		real	wind direction (degrees, clockwise from north)
    pressure		real	station pressure (mb)"""
    
    collab = collab.split('\n')
    collab = [c.split()[0] for c in collab]

    # read the file
    df = pd.read_csv(p2f, skiprows=2, sep= r'\s+',
                    #  delim_whitespace = True, 
                     names = range(48))#, names = collab)
    
    # there are some redundant columns full of zeros, no idea what they man?!?
    df = df.iloc[:,list(range(8)) + list(range(8,48,2))].copy()

    # assign the column name
    df.columns = collab
    
    # generate a datetime index
    df.index = df.apply(lambda row: pd.to_datetime(f'{row.year:0.0f}{row.jday:03.0f}', format = '%Y%j')+pd.to_timedelta(row['dt'], unit = 'h'), axis = 1)

    # remove columns that are no longer needed
    df.drop(['year', 'jday', 'month', 'day', 'hour', 'min', 'dt'], axis = 1, inplace = True)
    
    # set invalid values to nan
    df[df == -9999.9] = np.nan
    # convert to xarray dataset
    df.index.name = 'datetime'
    ds = df.to_xarray()
    return ds


def generate_netcdfs(p2fld = '/nfs/iftp/aftp/data/radiation/surfrad/',
                     p2fldout = '/nfs/grad/surfrad/products_level1/radiation_netcdf/',
                     gui = False,
                     verbose = False,
                     test = False):
    """
    Convertes all the radiation data from ascii to netcdf. In principle I could 
    use the NCEI data ... I think there was additional data that is not stored 
    in the NCEI data

    Parameters
    ----------
    gui : bool, optional
        If eccecuted in jupyter this will give you a progressbar. The default is False.

    Returns
    -------
    None.

    """
    version = '1.1'
    """
    changelog:
    1.0: initial version
    1.1: save in annual folders instead of all in one folder
    """
    out = {}
    p2fld = pl.Path(p2fld)
    
    
    p2fldout = pl.Path(p2fldout) / f'v{version}'
    
    sites = ['dra','bon','fpk','gwn','psu','sxf','tbl']
    # sites = sites[1:]
    truncate = False#['20140101', '20210101']
    overwrite = False
    progressbars = {}
    
    if gui:
        import ipywidgets as widgets
        from IPython.display import display
        for site in sites:
            pg = widgets.IntProgress(min=0) # instantiate the bar
            lab = widgets.Label(value = '0')
            sitelab = widgets.Label(value = site)
            f = widgets.HBox([sitelab,pg, lab])
            display(f) # display the bar
            progressbars[site] = [pg,lab]
    totlnum = 0   
    for site in sites:
        start = pd.Timestamp.now()
        
        if gui:
            pg,lab = progressbars[site]
        p2fldouts = p2fldout.joinpath(site)
        if verbose:
            print(f'output folder:{p2fldouts}')
        assert(p2fldouts.parent.exists()), f'Output parent folder does not exist: {p2fldouts.parent}.'
        p2fldouts.mkdir(exist_ok=True, 
                        # parents=True
                        )
        p2flds = p2fld.joinpath(site)
        if verbose:
            print(f'p2flds: {p2flds}')
            
        years = [int(y.name) for y in p2flds.glob('*') if y.is_dir()]
        years.sort()
        if verbose:
            print(f'years: {years}')
            
        dfy = pd.DataFrame(index = years)
        # dfy = dfy.truncate(2014, 2021)
        yearcontent = []
        for year in dfy.index:
            # break
    
            p2fldsy = p2flds.joinpath(str(year))
            fns = pd.DataFrame(p2fldsy.glob('*'), columns = ['p2f'])
            fns.index = fns.apply(lambda row: pd.to_datetime(f'{year}{row.p2f.name.split(".")[0][-3:]}', format = '%Y%j'), axis = 1)
            yearcontent.append(fns)
    
        fns = pd.concat(yearcontent)
        fns.sort_index(inplace=True)
    
        if truncate:
            fns = fns.truncate(*truncate)
        
        fns['p2foutfull'] = fns.apply(lambda row: p2fldouts / f'{row.name.year:04d}/srf_rad_full_{site}_{row.name.year:04d}{row.name.month:02d}{row.name.day:02d}.nc', axis = 1)
        
        if not overwrite:
            fns = fns[~(fns.apply(lambda row: row.p2foutfull.is_file(), axis = 1))]
            
        totlnum += fns.shape[0]
        
        if gui:
            pgmax = fns.shape[0]
            pg.max = pgmax
        # break
        for e,(idx, row) in enumerate(fns.iterrows()): 
            if gui:
                if e!=0:
                    now = pd.Timestamp.now()
                    i = e+1
                    minremain = ((now - start)/pd.to_timedelta(1,'m')/i)*(pgmax - i)
                    pg.value += 1
                    lab.value = f'{100*e/pgmax:0.1f}% {e}({pgmax}) t remaining: {minremain:0.0f} min'
            ds = read_surfrad(row.p2f)
            # assert(False)
            out['pw'] = fns
            out['ds_sample'] = ds
            if test:
                return out
            row.p2foutfull.parent.mkdir(exist_ok=True)
            ds.to_netcdf(row.p2foutfull)
    out['numprocessed'] = totlnum
    return out
            
            
def concat2hourlymean(gui = False):
    p2fld = pl.Path('/mnt/telg/data/grad/surfrad/radiation/')
    p2fldout = pl.Path('/mnt/telg/data/grad/surfrad/radiation_hourly/')
    truncate = False# ['20140101', '20230101']
    sites = ['dra','tbl','bon','dra','fpk','gwn','psu','sxf']
    
    if gui:
        progressbars = {}
        for site in sites:
            pg = widgets.IntProgress(min=0) # instantiate the bar
            lab = widgets.Label(value = '0')
            sitelab = widgets.Label(value = site)
            f = widgets.HBox([sitelab,pg, lab])
            display(f) # display the bar
            progressbars[site] = [pg,lab]
    
    start = pd.Timestamp.now()
    for site in sites:
        if gui:
            pg,lab = progressbars[site]
        df = pd.DataFrame(p2fld.joinpath(site).glob('*.nc'), columns=['p2fin'])
        dt = df.apply(lambda row: pd.to_datetime(row.p2fin.name.split('_')[-1].split('.')[0]), axis=1)
    
        df['date'] = dt
        df.index = dt
    
        p2fldouts = p2fldout.joinpath(site)
        p2fldouts.mkdir(exist_ok=True)
    
        df['p2fout'] =df.apply(lambda row: p2fldouts.joinpath(f'srf_rad_hourly_{site}_{row.date.year:04d}{row.date.month:02d}.nc'), axis =1)
    
        df = df[~(df.apply(lambda row: row.p2fout.is_file(), axis = 1))]
    
        df.sort_index(inplace=True)
        
        if truncate:
            df = df.truncate(*truncate)
        
        if gui:
            pgmax = len(df.groupby(df.p2fout))#df.shape[0]
            pg.max = pgmax
        for e,(p2fnout, grp) in enumerate(df.groupby(df.p2fout)):    
            if e!=0:
                now = pd.Timestamp.now()
                minremain = ((now - start)/pd.to_timedelta(1,'m')/e)*(pgmax - e)
                if gui:
                    pg.value += 1
                    lab.value = f'{100*e/pgmax:0.1f}% {e}({pgmax}) t remaining: {minremain:0.0f} min'
            if p2fnout == df.iloc[-1].p2fout:
                print('last one skipped')
                continue
            ds = xr.open_mfdataset(grp.p2fin, lock=False)
            dsrs = ds.resample(datetime = '1h').mean()
            # break
            # time.sleep(1)
            dsrs.to_netcdf(p2fnout)