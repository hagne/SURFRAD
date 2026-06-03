"""
created on 2026-02-06

@author: hagen telg

This module contains the code to process the CL61 ceilometer cloud product from row to level 1. 
"""

from attr import attrs
import xarray as xr
import pandas as pd
import productomator.worker as prowo
import numpy as np

nan = ''
df_ncvars_lit = [{'variable': 'cloud_base_heights',
  'use with name': 'cloud_base_heights',
  'units': 'm',
  'long_name': 'heights (range) of the detected cloud bases',
  'averaging_time_in_seconds': nan},
 {'variable': 'beta_att',
  'use with name': 'backscatter_profile',
  'units': '1/(m*sr)',
  'long_name': 'attenuated volume backscatter coefficient',
  'averaging_time_in_seconds': 120.0},
 {'variable': 'linear_depol_ratio',
  'use with name': 'linear_depol_ratio',
  'units': nan,
  'long_name': 'linear depolarisation ratio of the backscatter volume',
  'averaging_time_in_seconds': 120.0}]
df_ncvars = pd.DataFrame(df_ncvars_lit)

class Cl61CloudLevel1_v0_1(prowo.WorkplannerDaily):
    def __init__(self, *args, **kwargs):
        kwargs['version'] = '0.1'
        super().__init__(*args, **kwargs)
        self._ncvars = None
        self.site_name = kwargs['site_name']
        self.site_code = kwargs['site_code']
        self.lat = kwargs['lat']
        self.lon = kwargs['lon']

    @property
    def ncvariables(self):
        if self._ncvars is None:
            # load the table that describes how variables are handled
            p2fattrs = 'CL61_netcdf_variables.csv'
            df = pd.read_csv(p2fattrs, 
                             # index_col= 0,
                            )
            df.columns = [c.replace(' ', '_') for c in df.columns]
            #remove the empty lines and the group names
            df = df[~df.variable.isin(['/status/', '/monitoring/', np.nan])]
            
            # remove all where 'use_with_name' isnan, if that value is not set we are not including it in the product
            df_ncvars = df[~df['use_with_name'].isna()]
            df_ncvars = df_ncvars.dropna(axis = 1, how='all')
            self._ncvars = df_ncvars
        return self._ncvars

    def open_single_file(self, p2f):
        df_ncvars = self.ncvariables
        dsg = xr.open_groups(p2f) #open_datatree does not work as time is not aligned
        
        # loop over variables and pic/rename/average ... and add to a new tree
        tree = {}
        for g in dsg:
            ds_orig = dsg[g]
            try:
                ds_orig = ds_orig.drop_vars(['longitude', 'latitude'])
            except:
                pass
            ds = xr.Dataset()
            
            for idx, row  in df_ncvars.iterrows():
                var = row.variable
                if var not in ds_orig:
                    continue
                varnew = row.use_with_name
                attrs = row.drop(['variable','use_with_name', 'remove_time_coordinate', 'report_as_attribute','CL51_product_equivalent'])
                attrs = attrs.dropna()
                attrs['original_variable_name'] = var
                da = ds_orig[var]
                if row.remove_time_coordinate or row.report_as_attribute:
                    da = np.float32(da.mean())
                if row.report_as_attribute:
                    ds.attrs[varnew] = da
                else:
                    ds[varnew] = da.astype(np.float32)
                    ds[varnew].attrs = attrs
                # break
            if len(ds.variables) > 0:
                tree[g] = ds
            # break
        
        # The time in the monitoring is by ~ a second different than the main varibles, xarray does not like that
        if '/monitoring' in tree:
            tree["/monitoring"] = tree["/monitoring"].assign_coords(time=tree["/"].time)
        
        tree = xr.DataTree.from_dict(tree)
        tree.attrs.update(dsg['/'].attrs)
        return tree
    
    def process_row(self, row = None, iloc = None, loc = None, save = True):
        """This is the method that does the particular work and will need to be overwritten in your subclass.
        Typical components:
        1. read the input file(s) (row.p2f_in)
        3. convert to xarray dataset (if needed)
        2. format the netcdf file
            2.1 add dataset attributes, creation datetime, creation software, server, site details, etc.
            2.2 add variable attributes, units, long_name, standard_name, etc.
        3. save the output file (row.p2f_out)
        
        Parameters
        ----------
        row : pandas.Series, optional
            A row from the workplan dataframe. This is how the process method callse this function.
        iloc : int, optional
            An integer index to select a row from the workplan dataframe.
        loc : index label, optional
            select a row by timestamp.
            """
        
        if iloc is not None:
            row = self.workplan.iloc[iloc]
        elif loc is not None:
            row = self.workplan.loc[loc]
        self.tp_row = row

        #######
        ## Open input files
        # dsall = xr.open_mfdataset(row.p2f_in)
        forest = [self.open_single_file(pf) for pf in row.p2f_in]
        dsall = xr.concat(forest, 'time')
        ## Do some processing here, e.g. add attributes, format the dataset, etc.
        # slice out the single day
        start = row.name
        end = start + pd.to_timedelta(1, 'D')
        ds = dsall.sel(time = slice(start, end))

        gattrs = {'serial_number' : ds.attrs['instrument_serial_number']}
        
        gattrs |= {'title':'Ceilometer cloud product',
                     'data_product_version': self.version,
                     'institution':'NOAA/GML/GRAD',
                     'author':'hagen.telg@noaa.gov',
                     'source':'Vaisala CL61 ceilometer',
                     # 'serial_number': ds_orig.instrument_serial_number,
                     'input_files': ', '.join([fn.as_posix() for fn in row.p2f_in]),#[fn.name for fn in poutg.path2raw]),
                     'Conventions':'CF-1.8',
                     # 'comments': ("The 'time' coordinate has been re-indexed to the nearest full minute, " 
                     #              "using the closest valid data value from within 1 minute prior to each timestamp. " 
                     #              "The data remains unprocessed, except for the default processing by Vaisala's software. "
                     #              "Quality control is limited to visual checks for implausible values or clear errors.")
                                                             }

        gattrs['campaign_network'] = 'NOAA/GML/GRAD/surfrad'
        gattrs['site_name'] = self.site_name
        gattrs['site_code'] = self.site_code
        gattrs['site_latitude'] = self.lat
        gattrs['site_longitude'] = self.lon
        gattrs['day_complete'] = row.day_complete.__str__()
        ds.attrs = gattrs
        ## Save the output file
        if save:
            # create parent directories 3 levels up if they do not exist
            row.p2f_out.parent.parent.parent.mkdir(exist_ok=True)
            row.p2f_out.parent.parent.mkdir(exist_ok=True)
            row.p2f_out.parent.mkdir(exist_ok=True)
            if self.verbose:
                print(f'Saving output file to {row.p2f_out}')
            ds.to_netcdf(row.p2f_out)
            
        dsall.close()
        ds.close()
        return ds