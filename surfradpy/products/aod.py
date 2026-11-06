import atmPy.aerosols.physics.column_optical_properties as atmcop
import socket
import atmPy.radiation.retrievals.spectral_irradiance as atmsi
import ibis
import productomator.worker as prowo
import xarray as xr
import pandas as pd
import pathlib as pl
import numpy as np


class AOD(prowo.Workplanner):
    """something"""
    version = '0.1'
    def __init__(self,
                p2fld_in=None,
                supplemental_data=None,
                path2absorption_correction_coeff_1625 = '/Volumes/grad/surfrad/products_level2/aod/1625nm_absorption_correction_coefficience.nc',
                path2v0database = None,
                p2fld_out=None,
                # database=None,
                output_file_format='productname_{date}.nc',
                file_name_format='*{date:%Y%m%d}*',
                start=None,
                end=None,
                days=10,
                input_directory_structure='yearly',
                # output_directory_structure=None,
                file_complete_check=True,
                reporter=None,
                verbose=False,
                site = None,
                ):
        super().__init__(
                        p2fld_in=p2fld_in,
                        supplemental_data=supplemental_data,
                        p2fld_out=p2fld_out,
                        database=None,
                        output_file_format=output_file_format,
                        file_name_format=file_name_format,
                        start=start,
                        end=end,
                        days=days,
                        input_directory_structure=input_directory_structure,
                        output_directory_structure=None,
                        file_complete_check=file_complete_check,
                        reporter=reporter,
                        verbose=verbose,
                        site = site,
                        version = self.version
                    )
        self.path2absorption_correction_coeff_1625 = path2absorption_correction_coeff_1625
        self.path2v0database = path2v0database


    def process_row(self, row=None, iloc=None, loc=None, save=True):
        if iloc is not None:
            row = self.workplan.iloc[iloc]
        elif loc is not None:
            row = self.workplan.loc[loc]
        self.tp_row = row
        
        out = {}
        ds = xr.open_dataset(row.p2f_in)
        
        drop_channels = [940]#, 1625]
        channels_sel = [int(i) for i in ds.channel if int(i) not in drop_channels]
        
        ds = ds.sel(channel = channels_sel, drop=True)
        
        si = atmsi.CombinedGlobalDiffuseDirect(ds.copy())
        
        db = ibis.sqlite.connect(self.path2v0database)
        dbv0 = db.table('gam_V0_prediction')
        calibration = dbv0.filter(dbv0.date == row.name).to_pandas()
        out['calibration'] = calibration
        calibration.index = calibration.channel
        calibration = calibration.drop(['date', 'channel'], axis=1).to_xarray()
        sical = si.apply_calibration_langley(calibration)
                
        # move to direct normal and work from there
        dn = sical.direct_normal_irradiation
        dn.dataset.attrs['serial_no'] = dn.dataset.attrs['sn_mfrsr']
        
        # dn.skip_1625_channel = True
        
        # load the ozon data and apply
        p2fmet = row.p2f_in_model #'/Volumes/stu3data2/Model_data/geos_fp/tavg1_2d_slv_Nx/surfrad/v0.1_wrongsurfacetemp/2026/surfrad_geosfp_tavg1_2d_slv_Nx_20261019.nc'
        p2f_ozone = pl.Path(p2fmet)
        
        dsmet = xr.open_dataset(p2fmet)
        dsmet = dsmet.sel(site = 'tbl', drop = True)
        
        dsmetrn = dsmet.rename({#'ps':'pressure',
                                #'ts': 'temperature', # !!!!! This is still the wrong temperature!!!
                                'to3': 'ozone',
                                'tqv': 'precipitable_water'})
        
        # dsmetrn['pressure'] = dsmetrn.pressure / 100
        # dsmetrn.pressure.attrs['units'] = 'hPa'
        
        # dsmetrn['temperature'] = dsmetrn.temperature - 273.15
        # dsmetrn.temperature.attrs['units'] = 'DegC'
        dsmetrn['precipitable_water'] = dsmetrn.precipitable_water /10
        dsmetrn.precipitable_water.attrs['units'] = 'cm'
        
        dsmetrn.ozone.attrs['units'] = 'DU'
        
        dn.ozone_data = dsmetrn
        
        # load the met data and apply
        p2f = row.p2f_in_met#'/Volumes/grad/surfrad/scaled/tbl/2026/gradobs-scaled.tbl.20261019.dat'
        p2f_met = pl.Path(p2f)
        
        df = pd.read_csv(p2f)
        df.index = pd.to_datetime(df.DateTime)
        df.index.name = df.index.name.lower()
        df = df.loc[:,['2m_t[DegC]', 'P[hPa]']]
        renameing = [[c] + [e.strip(']') for e in c.split('[')] for c in df.columns]
        df = df.rename({r[0]: r[1] for r in renameing}, axis = 1)
        dsmetsr = df.to_xarray()
        for rn in renameing:
            dsmetsr[rn[1]].attrs['units'] = rn[2]
        dsmetrs = dsmetsr.rename({"2m_t":'temperature', 'P': 'pressure'})
        
        dn.met_data = dsmetrs

        ds_abscorrcoeff1625 = xr.open_dataset(self.path2absorption_correction_coeff_1625)
        dn.absorption_correction_coeff_1625 = ds_abscorrcoeff1625
        
        dn.precipitable_water = dsmetrn
        
        dni = dn
        dni.aod #this will trigger the retrieval and populate the raw_data
        prod = dni.dataset.copy()
        # create cloud flag
        aod = atmcop.AOD_AOT(dni.dataset)
        cm = aod.cloudmask.cloudmask_michalsky
        cm = cm.where(~ np.isnan(cm), 2).astype('int8')
        cm.attrs.update({
            "long_name": "cloud detection flag",
            "standard_name": "status_flag",
            "units": "1",
            "flag_values": [0, 1, 2],
            "flag_meanings": "no_cloud_detected cloud_detected algorithm_not_run",
            "description": (
                "0 = no cloud detected; 1 = cloud detected; "
                "2 = algorithm could not be evaluated"
            ),
        })
        prod['cloud_flag'] = cm
        prod = prod.rename({'aod':'aod_unscreened'})
        prod['aod'] = prod.aod_unscreened.where(prod.cloud_flag == 0)
        prod.aod.attrs['long_name'] = 'cloud screened spectral aerosol optical depth'
        
        od1625 = dni.od_co2_ch4_h2o.sel(channel = 1625, drop = True)
        prod['od_1625nm_co2'] = od1625.co2
        prod['od_1625nm_ch4'] = od1625.ch4
        prod['od_1625nm_h2o'] = od1625.h2o_5cm
        
        #####
        # angstrom exp
        prod['angstrom_exponent'] = dni.angstrom_exponent.where(prod.cloud_flag == 0)
        prod.angstrom_exponent.attrs['long_name'] ="aerosol Ångström exponent from AOD at 500 and 870 nm"
        prod.angstrom_exponent.attrs['ancillary_variables'] = "angstrom_exponent_uncertainty cloud_flag"
        prod['angstrom_exponent_uncertainty'] = dni.angstrom_exponent_uncertainty.where(prod.cloud_flag == 0)
        prod['angstrom_turbidity_coefficient'] = dni.angstrom_turbidity_coefficient.where(prod.cloud_flag == 0)
        prod.angstrom_turbidity_coefficient.attrs['long_name'] = "Ångström turbidity coefficient from AOD at 500 and 870 nm"
        prod.angstrom_turbidity_coefficient.attrs.pop('ancillary_variables')


        #####
        # calibration related
        prod['langley_calibration_v0'] = calibration.V0
        prod.langley_calibration_v0.attrs.update({
            'long_name': 'Langley calibration extraterrestrial signal',
            'units': 'mV',
            'comment': 'Predicted V0 used to calibrate the MFRSR measurements.',
        })
        prod['days_since_last_valid_langley_calibration'] = calibration.days2valid_langley_before
        prod.days_since_last_valid_langley_calibration.attrs.update({
            'long_name': 'days since the most recent valid Langley calibration',
            'units': 'days',
            'comment': 'Elapsed time from the product date to the most recent valid Langley calibration.',
        })
        prod['irradiance_top_of_atmosphere'] = prod.toa_spectral_irradiance.mean(dim = 'datetime')
        prod.irradiance_top_of_atmosphere.attrs.update({
            'long_name': 'mean top-of-atmosphere spectral solar irradiance',
            'standard_name': 'solar_irradiance_per_unit_wavelength',
            'units': 'W m-2 nm-1',
            'cell_methods': 'datetime: mean',
            'source': 'ASTM E490-00a zero-air-mass solar spectral irradiance spectrum',
            'references': (
                'ASTM International, ASTM E490-00a, Standard Solar Constant and '
                'Zero Air Mass Solar Spectral Irradiance Tables, '
                'https://doi.org/10.1520/E0490-00A'
            ),
            'comment': (
                'Linearly interpolated to the MFRSR channel center wavelengths and '
                'adjusted for Earth-Sun distance.'
            ),
        })
        #####
        # order variables
        self.tp_prod_preorder = prod.copy()
        prod = prod[["aod",
                     "aod_unscreened",
                     "aod_uncertainty",
                     "angstrom_exponent",
                     "angstrom_exponent_uncertainty",
                     "angstrom_turbidity_coefficient", 
                     "cloud_flag",
                     "langley_calibration_v0",
                     "days_since_last_valid_langley_calibration",
                     "irradiance_top_of_atmosphere",
                    "precipitable_water",
                    "transmission",
                    "global_horizontal",
                    "diffuse_horizontal",
                    "direct_normal",
                    "irradiance_uncertainty",
                    "od_rayleigh",
                    "od_ozone",
                    "od_1625nm_co2",
                    "od_1625nm_ch4",
                    "od_1625nm_h2o",
                    "channel_wavelength",
                    # "cosine_calibraion_direct",
                    "solar_zenith_angle",
                    "solar_azimuth_angle",
                     "airmass",
                    "pressure",
                    # "temperature",
                    "ozone",
                    ]]
        
        
        # add some attributes
        attrs = dict(title = "MFRSR calibrated spectral irradiance and aerosol optical depth",#, and precipitable water vapor",
                     summary = f"Calibrated spectral irradiances and retrieval products from MFRSR serial number {dni.dataset.serial_no} which is permanently deployed at the SURFRAD {self.site} observatory.",
                     # description = 'Spectral AOD, downwelling fluxes, and precipitable water from NOAA/GRAD MFRSR instruments during 2025 FRC in Davos',
                     Conventions = "CF-1.11, ACDD-1.3 [todo: is this correct?]",
                     institution = "NOAA Global Monitoring Laboratory / Global Radiation, Aerosols, and Clouds Division",
                     project = "SURFRAD",
                     network = "SURFRAD",
                     # observatory_name = "Physikalisch-Meteorologisches Observatorium Davos / World Radiation Center",
                     # observatory_id = "PMOD-WRC",
                     site_id = self.site,
                     site_name = prod.attrs['site_name'],
                     site_longitude = f'{ds.attrs['site_longitude']:0.6f}',
                     site_latitude = f'{ds.attrs['site_latitude']:0.6f}',
                     site_elevation = f'{ds.attrs['site_elevation']:0.0f}',
                     instrument = "Multifilter Rotating Shadowband Radiometer",
                     # instrument_serial_no = prod.attrs['serial_no'],
                     # instrument = prod.attrs['instrument'],
                     input_files = ', '.join([row.p2f_in.as_posix(), 
                                              p2f_met.as_posix(),
                                              p2f_ozone.as_posix(),
                                              self.path2v0database,
                                             ]
                                            ),
                     product_version = self.version
                    )          
        attrs.update({
         'site_latitude': np.float64(40.12498),
         'site_longitude': np.float64(-105.2368),
         'site_elevation': np.float64(1689.0)})
        
        timestamp = pd.Timestamp.now(tz='UTC').isoformat()
        ds.attrs['processing_date'] = timestamp
        entry = f'{timestamp}: Generated MFRSR AOD product using surfrad_aod v{self.version}'
        ds.attrs['history'] = '\n'.join(
            filter(None, [ds.attrs.get('history'), entry])
        )
        attrs['processing_server'] = socket.gethostname()
        # attrs['processing_class'] = f"This file was generated using{self.__class__.__module__}.{self.__class__.__qualname__}"
        
        prod.attrs = attrs
        
        ###################################
        # Irradiance stuff
        
        irradiance_comment = (
            "Instrument calibration was performed using Langley-derived top-of-atmosphere "
            "(extraterrestrial) signals. Irradiance units were assigned by scaling the "
            "Langley calibration to bandpass-weighted top-of-atmosphere solar irradiances "
            "computed from the ASTM E490-00a AM0 reference spectrum."
        )
        
        attrs = dict(long_name = "spectral atmospheric transmission",
                     units = "1",
                     comment = "Direct-beam atmospheric transmission derived from langley calibrated MFRSR measurements."
                    )
        
        prod.transmission.attrs = attrs
        prod = prod.rename({'transmission': 'direct_beam_transmittance'})
        
        attrs = dict(long_name = "spectral global horizontal irradiance",
                     units = "W m-2 nm-1",
                     comment =irradiance_comment
                    )
        prod.global_horizontal.attrs = attrs
        prod.global_horizontal.attrs
        
        attrs = dict(long_name = "spectral diffuse horizontal irradiance",
                     units = "W m-2 nm-1",
                     comment = irradiance_comment
                    )
        prod.diffuse_horizontal.attrs = attrs
        prod.diffuse_horizontal.attrs
        
        attrs = dict(long_name = "spectral direct normal irradiance",
                     units = "W m-2 nm-1",
                     comment = irradiance_comment
                    )
        prod.direct_normal.attrs = attrs
        prod.direct_normal.attrs
        
        prod.irradiance_uncertainty.attrs['units'] = 'percent'
        
        ##############
        ### AOD related
        attrs = dict(long_name = "spectral aerosol optical depth",
                     standard_name = "atmosphere_optical_thickness_due_to_ambient_aerosol_particles",
                     units = "1",
                     ancillary_variables = "aod_uncertainty cloud_flag",
                     comment = "AOD after correction for Rayleigh scattering and gaseous absorption from O3, CH4, CO2, and H2O."
                    )
        prod.aod.attrs = attrs
        prod.aod.attrs
        
        attrs = dict(long_name = 'absolute uncertainty of aerosol optical depth',
                     units = "1",
                     comment = "One-sigma aerosol optical depth uncertainty estimated from the variability among the applicable Langley calibrations, with a small-sample correction applied.")
        prod.aod_uncertainty.attrs = attrs
        
        od_correction_comment = (
            "Modeled optical depth correction term calculated from auxiliary data and "
            "instrument channel properties; this variable is not a retrieved geophysical product."
        )
        
        attrs = dict(long_name = "Rayleigh optical depth correction",
                     units = "1",
                     comment = od_correction_comment)
        prod.od_rayleigh.attrs = attrs
        
        attrs = dict(long_name = "ozone optical depth correction",
                     units = "1",
                     comment = od_correction_comment)
        prod.od_ozone.attrs = attrs
        
        attrs = dict(long_name = "water vapor optical depth correction at 1625 nm",
                     units = "1",
                     comment = od_correction_comment)
        prod.od_1625nm_h2o.attrs = attrs
        
        attrs = dict(long_name = "carbon dioxide optical depth correction at 1625 nm",
                     units = "1",
                     comment = od_correction_comment)
        prod.od_1625nm_co2.attrs = attrs
        
        attrs = dict(long_name = "methane optical depth correction at 1625 nm",
                     units = "1",
                     comment = od_correction_comment)
        prod.od_1625nm_ch4.attrs = attrs
        
        attrs = dict(long_name = "effective center wavelength",
                     units = "nm",
                     comment = "Effective center wavelength for each MFRSR channel.",
                    )
        prod.channel_wavelength.attrs = attrs
        
        #################################################
        ### Coordinates
        attrs = dict(long_name = "nominal MFRSR channel wavelength",
                     units = "nm")
        prod.channel.attrs = attrs
        prod.channel.encoding['_FillValue'] = None
        
        prod.datetime.attrs.update({
            'standard_name': 'time',
            'long_name': 'Time of observation',
            'axis': 'T',
        })

        ####################
        ### day complete check
        if self.file_complete_check:
            prod.attrs['day_complete'] = row.day_complete.__str__()
        
        prod = prod.sel(angstrom_exponent_pair = '500_870', drop = True)
        out['prod'] = prod
        ds.close()
        dsmet.close()
        dsmetrn.close()
        dsmetrs.close()
        dsmetsr.close()
        ds_abscorrcoeff1625.close()
        if save:
            row.p2f_out.parent.mkdir(parents=True, exist_ok=True)
            prod.to_netcdf(row.p2f_out)
            prod.close()
        return out
