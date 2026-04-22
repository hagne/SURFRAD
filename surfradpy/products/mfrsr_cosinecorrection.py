import productomator as pm
import xarray as xr
import pandas as pd
import pathlib as pl
import atmPy.radiation.instrumentation.spectral as atmspc
import socket


class CalibrateMFRSR(pm.worker.Workplanner):
    def __init__(self,*args, **kwargs):
        kwargs['version'] = '0.2'
        super().__init__(*args, **kwargs)
        
    def process_row(self, row = None, iloc = None, loc = None, save = True):
        """This is the method that does the particular work and will need to be overwritten in your subclass.
        Typical components:
        1. read the input file(s) (row.p2f_in)
        3. convert to xarray dataset (if needed)
        2. format the netcdf file
            2.1 add dataset attributes, creation datetime, creation software, server, site details, etc.
            2.2 add variable attributes, units, long_name, standard_name, etc.
        3. save the output file (row.p2f_out)
        
        Changelog
        ---------
        v0.2: Bugfix in raw reader resulst in a 20 second shift in the timestamps. 

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

        #####
        # get last processed instance - usefull if processing depends on the previous day
        ######
        # lastrow = self.get_last_row_before_workplan()
        # if isinstance(lastrow, type(None)):
        #     assert(False), 'set defaults?'
        # dslast = xr.open_dataset(lastrow.p2f_out)

        #######
        ## Open input files
        #######
        ds = xr.open_dataset(row.p2f_in)
        ds.attrs['site_latitude'] = ds.attrs['latitude']
        ds.attrs['site_longitude'] = ds.attrs['longitude']
        ds.attrs['site_elevation'] = ds.attrs['elevation']
        ds_raw = ds
        # si_raw = atmspc.atmsi.CombinedGlobalDiffuseDirect(ds)
        
        ## Do some processing here, e.g. add attributes, format the dataset, etc.
        # get last calibration files

        p2fld_cal = pl.Path(f'/Volumes/grad/Calibration_facilities/cucf/NOAA_GRAD_VisMFRSRs/V0{ds_raw.sn_mfrsr}/Reported/')
        assert(p2fld_cal.is_dir()), f'Folder containing calibration files, {p2fld_cal}, does not exist.'
        
        # get the .cos and .spr file with the closest date prior to the currently processed file
        def get_calfile(caltype = 'cos'):
            cosfiles = list(p2fld_cal.glob(f'*.{caltype}'))
            if len(cosfiles) == 0:
                return False
                # raise ValueError(f'No .{caltype} files found in {p2fld_cal}')
        
            cosdates = [pd.to_datetime(c.name.split('_')[-1].split('.')[0]) for c in cosfiles]
            df = pd.DataFrame({'cosfile': cosfiles}, index = cosdates) 
            df.sort_index(inplace = True)
            df = df.truncate(after = row.name)
            if df.shape[0] == 0:
                return False
                # raise ValueError(f'No .{caltype} files found in {p2fld_cal} with date prior to {row.name}\n{df}')
            p2fcos = df.cosfile.iloc[-1]
            return p2fcos
        
        p2fcos = get_calfile('cos')
        if not p2fcos:
            if self.verbose:
                print(f'No .cos files found in {p2fld_cal} with date prior to {row.name}. Trying .sol (factory calibration) file instead')
            p2fcos = get_calfile('sol')
            if not p2fcos:
                raise ValueError(f'No .cos or .sol files found in {p2fld_cal} with date prior to {row.name}')
            
        p2fspr = get_calfile('spr')
        if not p2fspr:
            if self.verbose:
                print(f'No .spr files found in {p2fld_cal} with date prior to {row.name}. Trying .sol (factory calibration) file instead')
            p2fspr = get_calfile('sol')
            if not p2fspr:
                raise ValueError(f'No .spr or .sol files found in {p2fld_cal} with date prior to {row.name}')
            
        if self.verbose:
            print(f'Using cosine calibration file: {p2fcos}')
            print(f'Using spectral calibration file: {p2fspr}')
        # create instrument instance
        instrument = atmspc.Mfrsr(
                                spectral_calibration = p2fspr,
                                logger_calibration = None, # note needed, we are not doing absolute calibration, meaning consant cal-facotors will be corrected in the langley calibration step
                                head_calibration = None , # note needed, we are not doing absolute calibration, meaning consant cal-facotors will be corrected in the langley calibration step
                                cosine_responds = p2fcos
                                )

        # pass raw data for calibration
        dscc = instrument.raw2calibrated(ds_raw)
        ds = dscc.dataset

        # Format the dataset for export
        # dropvars = ['alltime', 'zenith_geometric', 'elevation_geometric', 'elevation', 'equation_of_time', 'airmass_absolute','direct_horizontal','azimuth','zenith']
        # ds = ds.drop_vars(dropvars)
        reorg = ['global_horizontal',
         'direct_normal',
         'diffuse_horizontal',
         'channel_wavelength',
         'cosine_calibraion_direct',
         'cosine_calibration_diffuse',
         'solar_zenith_angle',
         'solar_azimuth_angle',
         # 'airmass',
         # 'sun_earth_distance',
        ]

        ds = ds[reorg]

        ds.solar_azimuth_angle.attrs['long_name'] = 'solar azimuth angle'
        ds.solar_azimuth_angle.attrs['units'] = 'degrees'
        ds.solar_azimuth_angle.attrs['standard_name'] = 'solar_azimuth_angle'

        ds.solar_zenith_angle.attrs['long_name'] = 'solar zenith angle'
        ds.solar_zenith_angle.attrs['units'] = 'degrees'
        ds.solar_zenith_angle.attrs['standard_name'] = 'solar_zenith_angle'
        ds.solar_zenith_angle.attrs['comment'] = 'Solar zenith angle including standard atmospheric refraction correction. Calculated using the pvlib python package'


        #########
        # Format the dataset attributes
        #########
        dropattrs = [
                    'day_complete','latitude','longitude','elevation',
                    'path2file',
                    'source',
                    ]
        for a in dropattrs:
            ds.attrs.pop(a)

        ds.attrs['processing_date'] = pd.Timestamp.now().isoformat()
        ds.attrs['processing_server'] = socket.gethostname()
        ds.attrs['info'] = 'Cosine corrected MFRSR data. See documentation for details: https://surfradpy.readthedocs.io/en/latest/products/mfrsr_cosinecorrection.html'
        ds.attrs['processing_class'] = f"{self.__class__.__module__}.{self.__class__.__qualname__}"
        ds.attrs['parent_files'] = row.p2f_in.as_posix()
        ds.attrs['calibration_file_cosine'] = p2fcos.as_posix()
        ds.attrs['calibration_file_spectral'] = p2fspr.as_posix()
        ds.attrs['version'] = self.version
        ## Save the output file
        if save:
            row.p2f_out.parent.mkdir(parents = True, exist_ok = True)
            ds.to_netcdf(row.p2f_out)
        ds.close()
        return ds