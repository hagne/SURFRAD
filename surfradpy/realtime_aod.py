#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 13:38:10 2023

@author: htelg
"""

import pandas as pd
import pathlib as pl
import xarray as xr
import numpy as np
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.surfrad.surfrad as atmsrf
import atmPy.aerosols.physics.column_optical_properties as atmcop
import traceback
import surfradpy.radiation as srfrad


class mfrsr_AOD_lev0(object):
    def __init__(self, site=['tbl',]):
        # product_version = 0.4
        self.p2fld_langley_in = '/nfs/grad/Inst/MFR/SURFRAD/{site}/mfrsr/ccc/'
        # p2fld_out = f'/mnt/telg/data/grad/surfrad/mfrsr/langleys/'
        self.p2fld_langley_out = pl.Path(
            '/export/htelg/data/grad/surfrad/mfrsr/langleys.0.4/')
        self.p2fld_langley_concat = pl.Path(
            '/export/htelg/data/grad/surfrad/mfrsr/langleys_concat.0.4/tbl/')

        self.product_version = 3.0
        self.file_type = 'tu'
        self.aod_lim_max = 2.8
        # p2fld_base = '/nfs/grad/Inst/MFR/SURFRAD/{site}/mfrsr/ccc/'
        self.p2fld_base = '/nfs/grad/Inst/MFR/SURFRAD/{site}/mfrsr/{file_type}/'
        self.p2fld_out = f'/export/htelg/data/grad/surfrad/aod_3/{self.product_version}'
        
        
        # settings
        # self.date_start = '2022-12-01'
        self.date_start = '2017-01-01'
        # self.date_end = '2023-03-02'
        self.date_end = None
        # self.sites = ['gwn', 'psu',
        #          #          # 'bld', 'fpk', 'bon' Those do not seam to be real surfrad sites
        #          'sxf', 'tbl', 'bnd',
        #          'fpe',  # This is the real fort Peck
        #          'dra']
        self.sites = site
        self.overwrite = False
        
        
        self.remove_old_predictions = False
        # i guess Gary and Allen are creating all those review files for the gaps?!?
        self.remove_review = True
        self.lands = None

        self.overwrite_langley = False
        self._workplan_langleys = None
        self._workplan_langleys_concat = None
        self._workplan = None
    
        self.processing_complete = False
        
        self.no_processed_error = 0
        self.no_processed_success = 0
        self.no_processed_warning = 0
        
    @property
    def no_of_files_that_need_processing(self):
        return self.workplan.shape[0]

    # @property
    # def no_of_files_unprocessed(self):
    #     if not self.processing_complete:
    #         self.process_all()
    #     workplan_rest = self.workplan[~(self.workplan.apply(lambda row: row.p2f_out.is_file(), axis=1))]
    #     return workplan_rest.shape[0]
        
    def process_all(self, test = False, verbose = True):
        
        #### radiation
        # this is needed for the met data
        # if verbose:
        #     print('producing netcdf from radiation data', end=' ... ')
        # srfrad.generate_netcdfs(gui=False)
        # if verbose:
        #     print('done')
        
        #### langleys
        
        # self.workplan_langleys
        self.process_langleys(raise_error=True, verbose=True)
        # self.workplan_langleys_concat
        ds = self.process_langley_concat(verbose=True)
                
        #### the product
        
        # self.workplan
        row = self.workit(
            raise_error=False,
            skip_error=False,
            print_error_and_return=True,
        )
        self.processing_complete = True
        if test:
            return row
        else:
            return True

    @property
    def workplan_langleys_concat(self):
        if isinstance(self._workplan_langleys_concat, type(None)):
            def _in2outname(row):
                nsplit = row.p2f.name.split('_')
                nsplit[4] = nsplit[4][:-2]
                p2out = self.p2fld_langley_concat.joinpath(
                    '_'.join(nsplit[:-1]) + '.nc')
                return p2out

            # assert(self.workplan_langleys.shape[0] == 0), 'There are unprocessed langleys!!! Take care of those first'
            # TODO: This needs to get generalized
            # p2fld_out = pl.Path('/export/htelg/data/grad/surfrad/mfrsr/langleys/tbl/')
            p2fld_langley_out = self.p2fld_langley_out.joinpath('tbl')
            df = pd.DataFrame(p2fld_langley_out.glob('*am*'), columns=['p2f'])
            df.index = df.apply(lambda row: pd.to_datetime(
                row.p2f.name.split('_')[4]), axis=1)
            df['ampm'] = 'am'
            df_am = df

            df = pd.DataFrame(p2fld_langley_out.glob('*pm*'), columns=['p2f'])
            df.index = df.apply(lambda row: pd.to_datetime(
                row.p2f.name.split('_')[4]) + pd.to_timedelta(12, 'hours'), axis=1)
            df['ampm'] = 'pm'
            df_pm = df

            workplan = pd.concat([df_am, df_pm])

            workplan.sort_index(inplace=True)
            p2fout = workplan.apply(lambda row: _in2outname(row), axis=1)
            workplan['p2f_out'] = p2fout

            if 0:  # this was used to remove the last month that should not be concatinated since more files are coming in ... i now do the concatintaion but don't save it .... see in process langley function
                last = workplan.index[-1]
                end = pd.to_datetime(
                    f'{last.year:04d}-{last.month:02d}-01 00:00:00') - pd.to_timedelta(1, 's')
                workplan = workplan.truncate(after=end)

            workplan = workplan[~(workplan.apply(
                lambda row: row.p2f_out.is_file(), axis=1))]
            self._workplan_langleys_concat = workplan
        return self._workplan_langleys_concat

    
    def process_langley_concat(self, verbose = False):
        if verbose:
            print('process langley')
            print('---------------')
            
        workplan = self.workplan_langleys_concat
        last = workplan.iloc[-1].p2f_out
        
        for p2fout, grp in workplan.groupby('p2f_out'):
            i = 0
            for idx, row in grp.iterrows():
                if verbose:
                    print(f'{row.p2f}')
                dst = xr.open_dataset(row.p2f)
                # in case its pm this will ensure that the indes has no duplicates
                dst['datetime'] = [idx,]
                # break
                # if bool(dst.status == 'fit failed'):
                # print(dst.status)
                # FIXME: The following is due to inconsistancies in how the variable is set. The mess below could be fixed by unifying this reporting during the langley generations.
                if type(dst.status) == xr.DataArray:
                    status = dst.status.values
                    if status.shape == ():
                        status = str(status)
                    else:
                        status = status[0]
                else:
                    status = dst.status

                assert (isinstance(status, str)
                        ), f'Status should be string ...is: {type(status)}'
                # return dst
                # print(str(dst.status))
                if 'fail' in status:
                    print('f', end='')
                    dst.close()
                    continue
                elif status == 'not convergent':
                    print('f', end='')
                    dst.close()
                    continue
                elif status == 'converged':
                    # ds.close()
                    # continue
                    pass
                else:
                    assert (False), f'New status: {status}'

                dst = dst.drop_dims('airmass')

                if i == 0:
                    ds = dst
                else:
                    ds = xr.concat([ds, dst], dim='datetime')
                i += 1

            for var in ds:
                if ds[var].dtype == np.float64:
                    ds[var] = ds[var].astype(np.float32)
                elif ds[var].dtype == np.int64:
                    ds[var] = ds[var].astype(np.int32)
                elif ds[var].dtype == object:
                    continue
                # status variable
                elif ds[var].dtype == '<U9':
                    continue
                else:
                    assert (
                        False), f'Variable {var} has unexpected type: {ds[var].dtype}'

            print('.', end='')

            p2fout_partial = p2fout.parent.joinpath(p2fout.name+'.part')
            if p2fout == last:
                p2fout = p2fout_partial
                print('p', end='')
            else:
                # clean up old partials
                if p2fout_partial.is_file():
                    p2fout_partial.unlink()
            ds.to_netcdf(p2fout)

        print('Done')
        self._workplan_langleys_concat = workplan[~(
            workplan.apply(lambda row: row.p2f_out.is_file(), axis=1))]
        return ds

    @property
    def workplan_langleys(self):
        if isinstance(self._workplan_langleys, type(None)):
            p2fld_out = pl.Path(self.p2fld_langley_out)
            workplan = pd.DataFrame()
            for site in self.sites:
                p2fld = pl.Path(self.p2fld_langley_in.format(site=site))
                wt = pd.DataFrame(p2fld.glob('**/*.ccc'), columns=['p2f',])
                wt.apply(lambda row: row.p2f.parent.name, axis=1).unique()
                wt = wt[~(wt.apply(lambda row: row.p2f.parent.name ==
                          '.canbedeleted', axis=1))]
                # wt.apply(lambda row: row.p2f.name.split('_')[0][:3], axis =1)
                wt['site'] = site
                workplan = pd.concat([workplan, wt])

            workplan.index = workplan.apply(lambda row: pd.to_datetime(
                row.p2f.name.split('_')[0][3:]), axis=1)
            workplan.sort_index(inplace=True)

            workplan = workplan.truncate(self.date_start, self.date_end)

            # workplan['serial_no'] = workplan.apply(lambda row:int(row.p2f.name.split('_')[-1].split('.')[0]), axis = 1)
            workplan['serial_no'] = workplan.apply(
                lambda row: row.p2f.name.split('_')[-1].split('.')[0], axis=1)

            # generate paths to the output files ... am + pm
            ampm = 'am'
            p2fld_out_site = p2fld_out.joinpath(site)
            workplan[f'p2f_out_{ampm}'] = workplan.apply(lambda row: p2fld_out_site.joinpath(
                f'srf_{site}_mfrsr_{row.p2f.name.split("_")[-1].split(".")[0]}_{row.name.year:04d}{row.name.month:02d}{row.name.day:02d}_langley_{ampm}.nc'), axis=1)
            ampm = 'pm'
            workplan[f'p2f_out_{ampm}'] = workplan.apply(lambda row: p2fld_out_site.joinpath(
                f'srf_{site}_mfrsr_{row.p2f.name.split("_")[-1].split(".")[0]}_{row.name.year:04d}{row.name.month:02d}{row.name.day:02d}_langley_{ampm}.nc'), axis=1)

            # for some reason there are those review folders ... ignore them for now
            workplan = workplan[workplan.apply(
                lambda row: 'Review' not in row.p2f.as_posix(), axis=1)]

            # remove from workplan if both am and pm path exist
            if not self.overwrite_langley:
                workplan = workplan[~(workplan.apply(
                    lambda row: row.p2f_out_am.is_file() and row.p2f_out_pm.is_file(), axis=1))]

            self._workplan_langleys = workplan
        return self._workplan_langleys

    @workplan_langleys.setter
    def workplan_langleys(self, value):
        self._workplan_langleys = value

    def process_langleys(self, raise_error=False, verbose = False):
        workplan = self.workplan_langleys
        for idx, row in workplan.iterrows():
            row.p2f_out_am.parent.mkdir(exist_ok=True)
            assert(row.p2f_out_am.parent.parent.is_dir()), f'not even the parent exists {row.p2f_out_am.parent.parent}'
            print('.', end='')
            sii = atmsrf.read_ccc(row.p2f)
            
            for ampm in ['am', 'pm']:
                print('>', end='')
                # there are still scenarios that cause errors, that need to be handled
                try:
                    # return sii
                    langleyinst = getattr(
                        sii.direct_normal_irradiation, f'langley_{ampm}')
                    if langleyinst.langleys.shape[0] == 0:
                        print('~', end='')
                        ds = xr.Dataset()
                        ds['status'] = 'fail, no data'
                    else:
                        result = langleyinst.clean(
                            use_channels=[415, 500, 670, 870])
                        spic = result['langley']
                        resstats = pd.DataFrame([spic.langley_fit_residual.mean(), spic.langley_fit_residual.median(
                        ), spic.langley_fit_residual.std()], index=['mean', 'median', 'sdt'])
                        resstats.index.name = 'resstats'
                        ds = xr.Dataset({
                            'langleys': spic.langleys,
                            'langley_fit_residual': spic.langley_fit_residual,
                            'langley_fitres': spic.langley_fitres,
                            'langley_residual_correlation_prop': spic.langley_residual_correlation_prop['determinant'],
                            'sp02_serial_no': sii.dataset.serial_no,
                            'valid_points': spic.langley_fit_residual.shape[0],
                            'residual_stats': resstats,
                        })
                        ds['cleaning_iterations'] = result['iterations']
                        ds['status'] = result['status']
                        ds = ds.expand_dims({'datetime': [row.name,]})
                    ds.to_netcdf(row[f'p2f_out_{ampm}'])
                except Exception as inst:
                    if inst.args[0] == 'attempt to get argmin of an empty sequence':
                        print('~', end='')
                        ds = xr.Dataset()
                        ds['status'] = 'fail, attempt to get argmin of an empty sequence.'
                        ds = ds.expand_dims({'datetime': [row.name,]})
                        ds.to_netcdf(row[f'p2f_out_{ampm}'])
                        continue
                    if raise_error:
                        # return sii
                        raise
                    else:
                        continue
            print('>', end='')
        print('Done')
        workplan = workplan[~(workplan.apply(
            lambda row: row.p2f_out_am.is_file() and row.p2f_out_pm.is_file(), axis=1))]
        self._workplan_langleys = workplan
        return None

    def make_product(self, path2file, product_version, verbose=False, langley_version='0.4', lands=None, aod_lim_max=2.8):
        # path2file = '/nfs/grad/Inst/MFR/SURFRAD/tbl/mfrsr/ccc/2021/tbl20210301_0660.ccc'
        p2f_ccc = pl.Path(path2file)
        sii = atmsrf.read_ccc(p2f_ccc)
        sii.direct_normal_irradiation.settings_calibration = 'atm_gam'
        sii.direct_normal_irradiation.path2absorption_correction_ceoff_1625 = '/export/htelg/products/grad/surfrad/aod_2/1625nm_absorption_correction_coefficience.nc'
        sii.direct_normal_irradiation._apply_calibration_atm_gam(p2fld=f'/export/htelg/data/grad/surfrad/mfrsr/langleys_concat.{langley_version}/tbl/',
                                                                 th=0.02,
                                                                 order_stderr=2,
                                                                 lam_overtime=2.5e4,
                                                                 ns_overtime=2,
                                                                 lam_season=1e4,
                                                                 ns_season=6,
                                                                 lands=lands)
        lands = sii.direct_normal_irradiation.lands

        ds = sii.direct_normal_irradiation.od_co2_ch4_h2o.transpose('datetime', 'channel')
        ds = ds.rename_vars({var: f'od_{var.split("_")[0]}' for var in ds})
        ds['aod'] = sii.direct_normal_irradiation.aod.transpose('datetime', 'channel')
        ds['od_rayleigh'] = sii.direct_normal_irradiation.od_rayleigh.transpose('datetime', 'channel')
        ds['od_total'] = sii.direct_normal_irradiation.od_total.transpose('datetime', 'channel')
        
        
        #### QF exceeds detectionlimit aod max
        dshighlim = (sii.direct_normal_irradiation.raw_data.global_horizontal_irradiation -
                     sii.direct_normal_irradiation.raw_data.diffuse_horizontal_irradiation) < aod_lim_max
        ds['flag_exceeds_aod_max'] = dshighlim
        ds.flag_exceeds_aod_max.attrs[
            'info'] = f'The difference between global and diffuse is smaller than {aod_lim_max}, which is ~4 standard deviations of the measurement noise.'

        # invalidate all date where certain flags are set
        ds = ds.where(~ds.flag_exceeds_aod_max, np.nan)

        #### angstrom exp
        if 1:
            aod = atmcop.AOD_AOT(ds.aod)
            angcombos = [(415, 670), (500, 870), (500, 1625), (870, 1625)]
            anglist = []
            # assert(False), 'aod is not defined anywhere!?!?' #something like: aodinst = atmcop.AOD_AOT(aod)
            for angcomb in angcombos:
                # col1 = 415
                # col2 = 673
                ang = aod.aod2angstrom_exponent(
                    column_1=angcomb[0], column_2=angcomb[1],)
                coordval = f'{angcomb[0]}_{angcomb[1]}'
                anglist.append(ang.expand_dims({'ang_channels': [coordval]}))
            ang = xr.concat(anglist, dim='ang_channels')
            ds['angstrom_exp'] = ang.transpose('datetime', 'ang_channels')
        
        #### legacy cloudmask
        ds['cloudmask_michalsky'] = aod.cloudmask.cloudmask_michalsky
        
        #### temporary
        ds['transmission'] = sii.direct_normal_irradiation.transmission
        ds['global_horizontal'] = sii.direct_normal_irradiation.raw_data.global_horizontal_irradiation
        ds['direct_normal'] = sii.direct_normal_irradiation.raw_data.direct_normal_irradiation
        ds['diffuse_horizontal'] = sii.direct_normal_irradiation.raw_data.diffuse_horizontal_irradiation

        sunpos = sii.direct_normal_irradiation.sun_position.copy()
        sunpos = sunpos.drop(['ampm'], axis=1)
        sunpos.columns.name = 'sun_params'
        ds['sun_position'] = sunpos

        # end temporary
        ds.attrs = {att: sii.direct_normal_irradiation.raw_data.attrs[att]
                    for att in sii.direct_normal_irradiation.raw_data.attrs if 'site' in att}
        ds.attrs['MFRSR_serial_no'] = sii.direct_normal_irradiation.raw_data.attrs['serial_no']
        ds.attrs['creation_timestamp'] = str(pd.Timestamp.now())
        ds.attrs['product_version'] = product_version
        # temporary variables for testing
        # is this the Langley claibration part, we might still need this
        if 0:
            V0df = sii.direct_normal_irradiation.tp_V0df
            V0df.index.name = 'channel'
            V0da = V0df.iloc[:, 0].to_xarray()
            V0da = V0da.expand_dims(
                {'date': [pd.to_datetime(pd.to_datetime(ds.datetime.values[0]).date())]})
        #
            ds['V0'] = V0da
            ds['sun_earth_distance'] = sii.direct_normal_irradiation.tp_sedistcorr.iloc[:, 0]
            # .iloc[:,0]
            ds['sun_earth_V0_correction'] = sii.direct_normal_irradiation.tp_calib_interp_secorr
        # output
        out = {'product': ds}
        out['sii'] = sii
        out['aod'] = aod
        if 0:
            out['VOdf'] = sii.direct_normal_irradiation.tp_V0df
        return out

    @property
    def workplan(self):
        if isinstance(self._workplan, type(None)):
            p2fld_out = pl.Path(self.p2fld_out)

            workplan = pd.DataFrame()
            for site in self.sites:
                assert(site == 'tbl'), 'Currently only works for TBL!!!'
                p2fld = pl.Path(self.p2fld_base.format(
                    site=site, file_type=self.file_type))

                wt = pd.DataFrame(p2fld.glob(
                    f'**/*.{self.file_type}'), columns=['p2f',])
                # break
                wt.apply(lambda row: row.p2f.parent.name, axis=1).unique()

                wt = wt[~(wt.apply(lambda row: row.p2f.parent.name ==
                          '.canbedeleted', axis=1))]

                # wt.apply(lambda row: row.p2f.name.split('_')[0][:3], axis =1)
                wt['site'] = site

                workplan = pd.concat([workplan, wt])

            if self.file_type == 'ccc':
                workplan.index = workplan.apply(lambda row: pd.to_datetime(
                    row.p2f.name.split('_')[0][3:]), axis=1)
            elif self.file_type == 'tu':
                workplan.index = workplan.apply(
                    lambda row: pd.to_datetime(row.p2f.name.split('_')[2]), axis=1)
            workplan.index.name = 'date'
            workplan.sort_index(inplace=True)

            workplan = workplan.truncate(self.date_start, self.date_end)

            workplan['p2f_out'] = workplan.apply(lambda row: p2fld_out.joinpath(row.site).joinpath(f'{row.name.year}').joinpath(
                f'srf_aod1625_{row.site}_{row.name.year:04d}{row.name.month:02d}{row.name.day:02d}.nc'), axis=1)
            
            #### site test
            sitefromname = workplan.apply(lambda row: row.p2f.name.split('_')[0], axis = 1)
            sitecheck = workplan.site == sitefromname
            
            if not sitecheck.all():
                print(f'There are {(~sitecheck).sum()} files where the site in the file does not match the site of the folder...they will be skipped!!')
                workplanex = workplan[~sitecheck]
                workplan = workplan[sitecheck]
                
            #### test if output files exist
            if not self.overwrite:
                workplan = workplan[~(workplan.apply(
                    lambda row: row.p2f_out.is_file(), axis=1))]

            if self.remove_review:
                workplan = workplan[~ workplan.apply(
                    lambda row: 'Review' in row.p2f.as_posix(), axis=1)]
            
            self._workplan = workplan
        return self._workplan
    
    @workplan.setter
    def workplan(self, value):
        self._workplan = value 
        
    def workit(self,
               raise_error = True,
               skip_error = False,
               print_error_and_return = False,):
        
        assert(sum([raise_error, skip_error, print_error_and_return]) == 1), 'Exactly one of the three kewargs, raise_error, skip_error, print_error_and_return, has to be True'
        finished = False
        
        for idx, row in self.workplan.iterrows():
            # if row.site != row.p2f.name.split('_')[0]:
            #     print('Rough file!! Site in file does not match folder!!.')
            #     continue
            print('|', end = '')
            try:
                out = self.make_product(row.p2f, product_version = self.product_version, lands = self.lands)
                ds = out['product']
                sii = out['sii']
                lands = sii.direct_normal_irradiation.lands
                self.lands = lands
                # reduce from 64 to 32 bit floats
                for var in ds:
                    ds[var] = ds[var].astype(np.float32)
        
                print('>', end = '')
        
                row.p2f_out.parent.parent.parent.mkdir(exist_ok=True)
                row.p2f_out.parent.parent.mkdir(exist_ok=True)
                row.p2f_out.parent.mkdir(exist_ok=True)
                
                # generate the path to the predicted data
                p2fout = row.p2f_out
                p2foutalt = p2fout.name.split('.')
                p2foutalt.insert(1, f'predict_{pd.Timestamp.now().date()}')
                p2foutalt = p2fout.parent.joinpath('.'.join(p2foutalt))
            
                # If after cutoff date for prediction make alternative path the actual path
                if row.name > sii.direct_normal_irradiation.lands.date_predict:
                    p2fout = p2foutalt
                    
                # cleanup... if there is still a predictive file after the cutoff date, remove it
                else:
                    if self.remove_old_predictions:
                        if p2foutalt.is_file():
                            p2foutalt.unlink()
                
                ds.to_netcdf(p2fout)
                self.no_processed_success +=1
            except FileNotFoundError:
                tb = traceback.format_exc()
                print(f'skip as file missing: {tb.split()[-1]}')
                self.no_processed_warning +=1
                continue
            except Exception:
                self.no_processed_error +=1
                if raise_error:
                    raise
                elif skip_error:
                    print(f'\nprocessing failed for file {row.p2f} ... skip')
                elif print_error_and_return:
                    print(traceback.format_exc())
                    return row
                else:
                    raise
                
            # print('', end = '')
        print('\ndone')
        finished = True
        return finished
    
    