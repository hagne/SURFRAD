"""
created on 2026-02-06

@author: hagen telg

This module contains the code to process the CL61 ceilometer cloud product from row to level 1. 
"""


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
        def __init__(self, 
                    p2fld_in,
                    p2fld_out,
                    date_from_name,
                    output_file_format,
                    glob_pattern_in='*.nc',
                    start=None,
                    end=None,
                    verbose=False,
                    path2surfrad_database = None,
                    site = None,
                    reporter = None
                     ):
            super().__init__(
                    p2fld_in,
                    p2fld_out,
                    date_from_name,
                    output_file_format,
                    glob_pattern_in=glob_pattern_in,
                    start=start,
                    end=end,
                    verbose=verbose,
                    reporter=reporter)
            
            self.path2surfrad_database = path2surfrad_database
            
            
        def process_row(self, row = None, iloc = None, loc = None):
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
            dsall = xr.open_mfdataset(row.p2f_in)
        
            ## Do some processing here, e.g. add attributes, format the dataset, etc.
            # slice out the single day
            start = row.name
            end = start + pd.to_timedelta(1, 'D')
            ds = dsall.sel(time = slice(start, end))

            # limit to relevant variables (also rename) and apply variable attributes
            ds_orig = ds
            ds = xr.Dataset()
            for idx, vrow  in df_ncvars.iterrows():
                attrs = {a:rowd[a] for a in vrow.to_dict() if rowd[a] != '' and a != 'variable'}
                dst = ds_orig[vrow.variable]
                dst.attrs = attrs
                ds[vrow['use with name']] = dst

            # lat and lon where empty when i developed this --> remove
            ds = ds.drop_vars(['longitude', 'latitude'])

            
            # assign global attributes
            ds.attrs = {'title':'Ceilometer cloud product',
                         'data_product_version':version,
                         'institution':'NOAA/GML/GRAD',
                         'author':'hagen.telg@noaa.gov',
                         'source':'Vaisala CL61 ceilometer',
                         'serial_number': ds_orig.instrument_serial_number,
                         'input_files': ', '.join([pf.as_posix() for pf in row.p2f_in]),#[fn.name for fn in poutg.path2raw]),
                         'Conventions':'CF-1.8',
                         # 'comments': ("The 'time' coordinate has been re-indexed to the nearest full minute, " 
                         #              "using the closest valid data value from within 1 minute prior to each timestamp. " 
                         #              "The data remains unprocessed, except for the default processing by Vaisala's software. "
                         #              "Quality control is limited to visual checks for implausible values or clear errors.")
                                                                 }
            
            ds.attrs['campaign_network'] = 'NOAA/GML/GRAD/surfrad'
            ds.attrs['site_name'] = site_name
            ds.attrs['site_code'] = site_code
            ds.attrs['site_latitude'] = lat
            ds.attrs['site_longitude'] = lon
                    
            ## Save the output file
            # ds.to_netcdf(row.p2f_out)
            dsall.close()
            ds.close()
            return ds