import pandas as _pd
import numpy as _np
import xarray as _xr
import os as _os
import tarfile as _tarfile
from .surfrad import locations as _locations
import platform as _platform
import hashlib as _hashlib
from pathlib import Path

columntransrules = {'solar_zenith_angle': 'Z',
                    'distance_from_sun': 'AU',
                    'best_estimate_shortwave_irradiance': 'BESW',
                    'total_downwelling_shortwave_irradiance': 'GSW',
                    'diffuse_downwelling_shortwave_irradiance': 'DIF',
                    'direct_downwelling_shortwave_irradiance': 'DIR',
                    'upwelling_shortwave_irradiance': 'SWup',
                    'downwelling_longwave_irradiance': 'LWdn',
                    'upwelling_longwave_irradiance': 'LWup',
                    'atmospheric_temperature': 'Ta',
                    'relative_humidity': 'RH',
                    'surface_air_pressure': 'Prs',
                    'downwelling_longwave_case_temperature': 'LWdTc',
                    'downwelling_longwave_dome_temperature': 'LWdTd',
                    'upwelling_longwave_case_temperature': 'LWuTc',
                    'upwelling_longwave_dome_temperature': 'LWuTd',
                    'QC_flag_01_thresholds_global_shortwave': 'qc1',
                    'QC_flag_02_thresholds_diffuse_shortwave': 'qc2',
                    'QC_flag_03_thresholds_direct_shortwave': 'qc3',
                    'QC_flag_04_thresholds_upwelling_shortwave': 'qc4',
                    'QC_flag_05_thresholds_downwelling_longwave': 'qc5',
                    'QC_flag_06_thresholds_upwelling_longwave': 'qc6',
                    'QC_flag_07_global_shortwave_over_sum': 'qc7',
                    'QC_flag_08_diffuse_over_global_shortwave': 'qc8',
                    'QC_flag_09_upwelling_shortwave_versus_sum': 'qc9',
                    'QC_flag_10_downwelling_longwave_to_temperature': 'qc10',
                    'QC_flag_11_upwelling_longwave_to_temperature': 'qc11',
                    'QC_flag_12_downwelling_longwave_to_upwelling_longwave': 'qc12',
                    'QC_flag_13_downwelling_longwave_case_temperature_vs_air_temperature': 'qc13',
                    'QC_flag_14_downwelling_longwave_dome_temperature_vs_air_temperature': 'qc14',
                    'QC_flag_15_upwelling_longwave_case_temperature_vs_air_temperature': 'qc15',
                    'QC_flag_16_upwelling_longwave_dome_temperature_vs_air_temperature': 'qc16',
                    'QC_flag_17_downwelling_longwave_case_temperature_vs_dome_temperature': 'qc17',
                    'QC_flag_18_upwelling_longwave_case_temperature_vs_dome_temperature': 'qc18',
                    'QC_flag_19_air_temperature': 'qc19',
                    'global_shortwave_irloss_flag': 'gflg',
                    'diffuse_shortwave_irloss_flag': 'dflg',
                    'downwelling_shortwave_irradiance_assuming_clear_sky': 'ClrSW',
                    'ir_loss_correction_to_diffuse_shortwave': 'DifCorr',
                    'ir_loss_correction_to_global_shortwave': 'GSWCorr',
                    'downwelling_erythemal_uvb_irradiance': 'UVB_down',
                    'photosynthetically_active_radiation': 'PAR_dn',
                    'wind_speed': 'WindSp',
                    'wind_direction': 'WindDr',
                    }


def parse_CDL_file(fname = '../data/SURFRAD_QCrad_metadata.cdl'):

    cdlname = fname

    with open(cdlname, 'r') as rein:
        lines = rein.readlines()

    what_is_it = None
    varlist = []
    global_atts = {}
    for l in lines:
        if 'variables:' in l:
            what_is_it = 'variable'
            continue
        elif 'global Attributes:' in l:
            what_is_it = 'global_attributes'
            continue
        elif l.strip() == '':
            continue
        elif l.strip() == '}':
#             print('done')
            break

        if what_is_it == 'variable':
            # tests if header and if then which type the data will have
            if l.strip().split()[0] in ('int','float','byte'):
                if '(' in l:
                    name, index = l.strip().split()[1].split('(')
                    index = index.replace(');', '')
                else:
                    name, index = l.strip().split()[1].strip(';'), None
                var = {}
                varlist.append(var)
                var['name_nc'] = name
                var['index'] = index
                var['attributes'] = {}
            else:
                assert(name in l)
                attname, value = l.strip().replace(name + ':', '').split('=')
                attname = attname.strip()
                value = value.strip()
                value = value.split('//')[0].strip().strip(';').strip('"').strip()
                # this will convert string numbers to int ... unless when it is a unit
                if attname != 'units':
                    try:

                        value = eval(value)
                    except:
                        pass
                var['attributes'][attname.strip()] = value
        if what_is_it == 'global_attributes':
            attname, value = l.strip().strip(':').strip(';').split('=', maxsplit=1)
            global_atts[attname.strip()] = value.strip().strip('"')

    now = _pd.datetime.utcnow()
    global_atts['history'] = 'This file was created {:04d}-{:02d}-{:02d} {:02d}:{:02d}:{:02d} on {} by {}.'.format(now.year, now.month, now.day, now.hour, now.minute, now.second, _platform.node(), _os.environ['LOGNAME'])

    # The _NCProperties parameter seams to be somewhat protected so I pop it out for now
    bla = global_atts.pop('_NCProperties')
    out = dict(variable_list = varlist,
              global_atts_dict = global_atts)
    return out

def read_data(fname='../data/bon_19950510.qdat', missing_data=((-9999.0, -9999.9), _np.nan), pressure_mbar2pa=True):
    data = _pd.read_csv(fname, delim_whitespace=True)

    # column names changed over time ... fixin it here
    data.rename({'GSWcorr': 'GSWCorr'}, axis=1, inplace=True)

    # creating timestamp
    datetimestr = data.Date.apply(lambda x: '{0:0>8}'.format(x)) + data.Time.apply(
        lambda x: '{0:0>4}'.format(x)) + 'UTC'  # '+0000'

    data.index = _pd.to_datetime(datetimestr, format="%Y%m%d%H%M%Z")
    data.index.name = 'time'

    # removing obsolete columns
    data = data.drop(['Date', 'Time'], axis=1)

    # our missing data is np.nan
    if missing_data:
        for i in missing_data[0]:
            data[data == i] = missing_data[1]

    # pressure mbar to pa
    if pressure_mbar2pa:
        data.Prs *= 100

    # change column names according to columntransrules
    columntransrules_rev = {v: k for k, v in columntransrules.items()}
    data.rename(columntransrules_rev, axis=1, inplace=True)

    # adjust dtypes to save space
    for col in data:
        if 'float' in str(data[col].dtype):
            data[col] = data[col].astype(_np.float32)
        elif 'int' in str(data[col].dtype):
            data[col] = data[col].astype(_np.int8)

    return data


def qcrad2netcdf(data_path_in, nc_path_out, cdl_dict, verbose=False):
    def save2netcdf(ds, fname='../data/Bondville_IL_1995_May_10.new.nc', compression=False):
        # turns out the compression results in a larger file tstststs
        ## set encoding
        if compression:
            for var in ds.variables:
                enc = ds[var].encoding
                enc['chunksizes'] = ds[var].shape
                enc['contiguous'] = False
                enc['complevel'] = 9
                enc['zlib'] = False
            #     # incase hdf5 is used at some point:
            #     enc['compression'] = 'gzip'
            #     enc['compression_opts'] = 9
        # creates folders as required
        _os.makedirs(_os.path.dirname(fname), exist_ok=True)
        # save
        ds.to_netcdf(fname,
                     format='NETCDF4_CLASSIC',
                     )
    def df2ds(data, variable_list, global_atts, location):
        ds = _xr.Dataset(data, attrs=global_atts)

        # add variable attributes
        for var in ds.variables:
            if var == 'time':
                continue
            atts = [v for v in variable_list if v['name_nc'] == var][0]['attributes'].copy()
            if '_FillValue' in atts.keys():
                if atts['_FillValue'] == 'NaNf':
                    atts['_FillValue'] = _np.nan
            ds.variables[var].attrs = atts

        # set location variables
        # in case the lat, lon, alt is used instead of full names:
        if 'lat' in location.keys():
            location = dict(latitude=location['lat'], longitude=location['lon'], altitude=location['alt'])
        for lv in location:
            ds[lv] = location[lv]
        return ds
    if verbose:
        print('pcrad2netcdf: {} -> {}'.format(data_path_in, nc_path_out), end = '...')
    data = read_data(data_path_in)
    location = [e for e in _locations if data_path_in.name.split('_')[0] in e['abbriviations']][0]
    ds = df2ds(data, cdl_dict['variable_list'], cdl_dict['global_atts_dict'], location)
    save2netcdf(ds, nc_path_out)
    if verbose:
        print('done')
    return location

def tar_nc(pot, todo, manifest = True, verbose = False):
    def generate_md5_checksum(fname):
        with open(fname, "rb") as file_to_check:
            data = file_to_check.read()
            md5 = _hashlib.md5(data).hexdigest()

        f_size = _os.stat(fname).st_size

        fname_out = fname + '.mnf'

        file_content = [_os.path.split(fname)[-1], str(md5), str(f_size)]

        with open(fname_out, 'w') as fout:
            fout.write(','.join(file_content))
    tar_mode = 'w:gz'
    # make sure folder exists:
    _os.makedirs(_os.path.dirname(pot), exist_ok=True)
    if verbose:
        print('tar_nc: {}'.format(pot), end = '...')
    tar = _tarfile.open(pot, mode=tar_mode)

    for po in todo.path_out:
        tar.add(po, arcname = po.name)
    tar.close()
    if manifest:
        generate_md5_checksum(pot)
    if verbose:
        print('done')

def qcrad2ncei(folder_in = '/Volumes/HTelg_4TB_Backup/GRAD/SURFRAD/qcrad_v3/',
               folder_out= '/Volumes/HTelg_4TB_Backup/GRAD/SURFRAD/NCEI/',
               folder_out_tar= '/Volumes/HTelg_4TB_Backup/GRAD/SURFRAD/NCEI_tar/',
               station_abb = 'bon',
               year = 1995,
               month = 1,
               overwrite = False,
               do_qcrad2nc = True,
               do_tar = True,
               do_manifest = True,
               test = False,
               verbose = False
              ):

    fl_ext = '.tar.gz'

    # generate a DataFrame with all files in sub folder and populate with relevant data
    paths_in = list(Path(folder_in).rglob("*.qdat"))
    index = [_pd.to_datetime(i.name.split('_')[1].split('.')[0]) for i in paths_in]
    station = [i.name.split('_')[0] for i in paths_in]
    station_name = [[e for e in _locations if i in e['abbriviations']][0]['name'] for i in station]
    station_state = [[e for e in _locations if i in e['abbriviations']][0]['state'] for i in station]

    df = _pd.DataFrame({'path_in': paths_in,
    #                     'path_out': paths_out,
    #                     'path_out_exists': paths_out_exist,
                        'station': station,
                        'station_name': station_name,
                        'station_state': station_state},
                       index = index)

    df['year'] = [i.year for i in df.index]
    df['month'] = [i.month for i in df.index]
    df['day'] = [i.day for i in df.index]
    df['month_name'] = [i.month_name()[:3] for i in df.index]
    df['path_out'] = df.apply(lambda x: Path('{}{}/{}_{}_{}_{}_{:02d}{}'.format(folder_out, x.station, x.station_name, x.station_state, x.year, x.month_name, x.day, '.nc')), axis = 1)
    df['path_out_exists'] = df.path_out.apply(lambda i: i.is_file())
    df['path_out_tar'] = df.apply(lambda x: '{}{}/{}_{}_{}_{}{}'.format(folder_out_tar, x.station, x.station_name, x.station_state, x.year, x.month_name, fl_ext), axis = 1)

    # the full list is neaded in case a tar is produced ...
    # this way a new tar is created when a new file belongs into an existing tar
    df_full = df.copy()

    # check existing for overwrite
    if not overwrite:
        df = df[df.path_out_exists == False]

    # select particular things ... like stations or times

    if station_abb:
        df = df[df.station == station_abb]
    if year:
        df = df[df.year == year]
    if month:
        df = df[df.month == month]

    # qcrad2netcdf the remaining todos in df
    if test:
        print('Files to be processed:')
        for fn in df.path_in:
            print('\t{}'.format(fn.as_posix()))
        return df
    if do_qcrad2nc:
        cdl_dict = parse_CDL_file()
        for idx,line in df.iterrows():
            qcrad2netcdf(line.path_in, line.path_out, cdl_dict, verbose=verbose)
    else:
        if verbose:
            print('No NetCDF files created since do_qcrad2nc == False')
    # get all filles from the full file list that match the desired tar file
    if do_tar:
        for pot in _np.unique(df.path_out_tar):
            todo = df_full[df_full.path_out_tar == pot]
            tar_nc(pot, todo, manifest=do_manifest, verbose= verbose)
    else:
        if verbose:
            print('No archives created since do_tar == False')

