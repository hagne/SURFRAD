import pandas as _pd
from pandas import errors as _pderrors
import numpy as _np
import xarray as _xr
import os as _os
import tarfile as _tarfile
from .surfrad import locations as _locations
import platform as _platform
import hashlib as _hashlib
from pathlib import Path as _Path

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

                attname = attname.strip()
                if attname == 'flag_values':
                    value = _np.array([bi.replace('b', '').strip() for bi in value.split(',')]).astype(
                        _np.int8)
                var['attributes'][attname] = value
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

# def make_log_entry(fname_log, txt):
#     with open(fname_log, 'a') as log:
#         timestamp = _pd.Timestamp(_pd.datetime.now())
#         log.write('{} {} {}'.format(timestamp, txt, '\n'))

def qcrad2netcdf(data_path_in, nc_path_out, cdl_dict, messages = None, verbose=False):
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
    def clean_atts(atts):
        if '_FillValue' in atts.keys():
            if atts['_FillValue'] == 'NaNf':
                atts['_FillValue'] = _np.nan
        # Chunck size just creates trouble ... remove it
        try:
            atts.pop('_ChunkSize')
        except KeyError:
            pass
        return atts

    def df2ds(data, variable_list, global_atts, location):
        ds = _xr.Dataset(data, attrs=global_atts)

        # add variable attributes
        for var in ds.variables:
            # if var == 'time':
            #     continue
            atts = [v for v in variable_list if v['name_nc'] == var][0]['attributes'].copy()

            atts = clean_atts(atts)
            # if '_FillValue' in atts.keys():
            #     if atts['_FillValue'] == 'NaNf':
            #         atts['_FillValue'] = _np.nan

            ds.variables[var].attrs = atts

        # set location variables
        # in case the lat, lon, alt is used instead of full names:
        if 'lat' in location.keys():
            location = dict(latitude=location['lat'], longitude=location['lon'], altitude=location['alt'])
        for lv in location:
            ds[lv] = _np.float32(location[lv])
            atts = [v for v in variable_list if v['name_nc'] == lv][0]['attributes']
            atts = clean_atts(atts)
            # if '_FillValue' in atts.keys():
            #     if atts['_FillValue'] == 'NaNf':
            #         atts['_FillValue'] = _np.nan
            ds[lv].attrs = atts

        return ds

    txt = 'pcrad2netcdf: {} -> {}'.format(data_path_in, nc_path_out)
    if verbose:
        print(txt, end = '...')
    if messages:
        messages.append(txt)
    data = read_data(data_path_in)
    location = [e for e in _locations if data_path_in.name.split('_')[0] in e['abbriviations']][0]
    ds = df2ds(data, cdl_dict['variable_list'], cdl_dict['global_atts_dict'], location)
    save2netcdf(ds, nc_path_out)
    if verbose:
        print('done')
    if messages:
        messages[-1] += ' ... done'
    return location

def tar_nc(pot, todo, manifest = True, messages = None, errors = [], verbose = False):
    def generate_md5_checksum(fname):
        with open(fname, "rb") as file_to_check:
            data = file_to_check.read()
            md5 = _hashlib.md5(data).hexdigest()

        f_size = _os.stat(fname).st_size

        fname_out = _Path(fname.as_posix() + '.mnf')

        file_content = [_os.path.split(fname)[-1], str(md5), str(f_size)]

        with open(fname_out, 'w') as fout:
            fout.write(','.join(file_content))

    tar_mode = 'w:gz'
    # make sure folder exists:
    _os.makedirs(_os.path.dirname(pot), exist_ok=True)

    txt = 'tar_nc: {}'.format(pot)
    if verbose:
        print(txt, end = '...')
    if messages:
        messages.append(txt)

    tar = _tarfile.open(pot, mode=tar_mode)

    for po in todo.path_out:
        try:
            tar.add(po, arcname = po.name)
            txt_ext = ' ... done.'
            remove = False
        except FileNotFoundError:
            txt_ext = 'FileNotFoundError: {}'.format(po)
            errors.append('{}: {}'.format(pot, txt_ext))
            remove = True

    tar.close()

    # remove tar if error accured
    if remove:
        _os.remove(pot)
    else:
        if manifest:
            generate_md5_checksum(pot)

    if verbose:
        print(txt_ext)
    if messages:
        messages[-1] += txt_ext

def create_todo(folder_in, folder_out, folder_out_tar, overwrite = False, station_abb = None, year = None, month = None):
    paths_in = list(_Path(folder_in).rglob("*.qdat"))
    if len(paths_in) == 0:
        raise ValueError(
            'There are no valid qcrad data files (with extention .qdat) in this folder or its sub-folders. Make sure the input-folder is correct.')

    # select files with a valid format and station
    valid_stations = ["bon", "tbl", "dra", "fpk", "gwn", "psu", "sxf"]
    paths_in = [pi for pi in paths_in if pi.name.split('_')[0] in valid_stations]

    # create the full list
    index = [_pd.to_datetime(i.name.split('_')[1].split('.')[0]) for i in paths_in]
    station = [i.name.split('_')[0] for i in paths_in]
    station_name = [[e for e in _locations if i in e['abbriviations']][0]['name'].replace(' ', '') for i in station]
    station_state = [[e for e in _locations if i in e['abbriviations']][0]['state'] for i in station]

    df = _pd.DataFrame({'path_in': paths_in,
                        #                     'path_out': paths_out,
                        #                     'path_out_exists': paths_out_exist,
                        'station': station,
                        'station_name': station_name,
                        'station_state': station_state},
                       index=index)
    df.sort_index(inplace=True)
    df['year'] = [i.year for i in df.index]
    df['month'] = [i.month for i in df.index]
    df['day'] = [i.day for i in df.index]
    df['month_name'] = [i.month_name()[:3] for i in df.index]
    df['path_out'] = df.apply(lambda x: _Path(
        '{}{}/{}_{}_{}_{}_{:02d}{}'.format(folder_out, x.station, x.station_name, x.station_state, x.year,
                                           x.month_name, x.day, '.nc')), axis=1)
    # take care of tar-archive
    ## generate a temporary tar file name
    df['path_out_tar'] = df.apply(
        lambda x: '{}{}/{}_{}_{}_{}'.format(folder_out_tar, x.station, 'ESRL-GMD-GRAD_v1.0_SURFRADQCRAD',
                                            x.station.upper(), x.year,
                                            x.month_name), axis=1)
    df['old_tar'] = False
    ## for each station set last month to do not process
    df['do_tar'] = True
    df['do_process'] = True
    for st in _np.unique(df.station):
        last = df.loc[df.station == st, :].iloc[-1]
        tst = df.station == st
        ty = df.year == last.year
        tm = df.month == last.month
        df.loc[tst & ty & tm, 'do_tar'] = False

    ## for all other generate the final name
    now = _pd.datetime.now()
    for dtp in _np.unique(df.loc[df.do_tar].path_out_tar):
        ttp = df.path_out_tar == dtp
        tdf = df.loc[ttp]
        newname = '{fot}{st}/ESRL-GMD-GRAD_v{version}_SURFRADQCRAD_{stu}_s{year}{month:02d}{day_s:02d}_e{year}{month:02d}{day_e:02d}_c{year_c}{month_c:02d}{day_c:02d}.tar.gz'.format(
            fot=folder_out_tar,
            version=1.0,
            st=tdf.station[0],
            stu=tdf.station[0].upper(),
            year=tdf.year[0],
            month=tdf.month[0],
            day_s=tdf.day[0],
            day_e=tdf.day[-1],
            year_c=now.year,
            month_c=now.month,
            day_c=now.day)  # x.year[0], x.month[0], x.day[0], x.day[-1])
        df.loc[ttp, 'path_out_tar'] = _Path(newname)

    def tar_exists(i, return_fname = False):
        pattern = '{}*.tar.gz'.format(i.name.split('_c')[0])
        matches = list(i.parent.glob(pattern))
        if len(matches) == 1:
            if return_fname:
                return matches[0]
            else:
                return True
        elif len(matches) > 1:
            mt = [m.as_posix() for m in matches]
            txtt = '\n'.join(mt)
            txt = 'This should not be possible, there can not be more then one other version of this tar!\n{}'.format(txtt)
            raise ValueError(txt)
        else:
            return False

    # select particular things ... like stations or times
    if station_abb:
        tst = df.station == station_abb
        df.loc[~tst, ['do_process','do_tar']] = False
    if year:
        ty = df.year == year
        df.loc[~ty, ['do_process','do_tar']] = False
    if month:
        tm = df.month == month
        df.loc[~tm, ['do_process','do_tar']] = False

    if not overwrite:
        df.loc[df.do_process,'do_process'] = ~df.loc[df.do_process, 'path_out'].apply(lambda i: i.is_file())
        df.loc[df.do_tar, 'do_tar'] = ~df.loc[df.do_tar, 'path_out_tar'].apply(tar_exists)
    else:
        df.loc[df.do_tar, 'old_tar'] = df.loc[df.do_tar, 'path_out_tar'].apply(lambda i: tar_exists(i, return_fname=True))
    return df

def qcrad2ncei(folder_in = '/Volumes/HTelg_4TB_Backup/GRAD/SURFRAD/qcrad_v3/',
               folder_out= '/Volumes/HTelg_4TB_Backup/GRAD/SURFRAD/NCEI/',
               folder_out_tar= '/Volumes/HTelg_4TB_Backup/GRAD/SURFRAD/NCEI_tar/',
               # fname_cdl = '../data/SURFRAD_QCrad_metadata.cdl',
               messages = None,
               errors = [],
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
    fname_cdl = _os.path.join(_os.path.split(__file__)[0], 'SURFRAD_QCrad_metadata.cdl')#'../data/SURFRAD_QCrad_metadata.cdl'

    # generate a DataFrame with all files in sub folder and populate with relevant data
    if verbose:
        print('creating todo list', end = '....')
    df = create_todo(folder_in, folder_out, folder_out_tar, overwrite = overwrite, station_abb = station_abb, year = year, month = month)
    if verbose:
        print('done')

    # qcrad2netcdf the remaining todos in df
    if test:
        print('Files to be processed:')
        for fn in df.path_in[df.do_process]:
            print('\t{}'.format(fn.as_posix()))
        print('Files tar-archives to be created:')
        for fn in _np.unique(df.path_out_tar[df.do_tar]):
            print('\t{}'.format(fn.as_posix()))
        return df
    if do_qcrad2nc:
        cdl_dict = parse_CDL_file(fname = fname_cdl)#'../data/SURFRAD_QCrad_metadata.cdl')
        for idx,line in df[df.do_process].iterrows():
            try:
                qcrad2netcdf(line.path_in, line.path_out, cdl_dict, messages= messages, verbose=verbose)
            except _pderrors.ParserError:
                txt = 'Failed to read data! Probably badly shaped.'
                if messages:
                    messages.append(txt)
                if verbose:
                    print(txt)
                errors.append('{}: {}'.format(line.path_in.as_posix(), txt))
    else:
        if verbose:
            print('No NetCDF files created since do_qcrad2nc == False')
    # get all filles from the full file list that match the desired tar file
    if do_tar:
        for pot in _np.unique(df[df.do_tar].path_out_tar):
            #remove existing tars (which are defined in todo.tar_old
            for ot in df.old_tar:
                if ot:
                    if ot.is_file():
                        if verbose:
                            print('removing old tar and mnf')
                        _Path(ot.as_posix() + '.mnf').unlink()
                        ot.unlink()
            # generate the new tars
            todo = df[df.path_out_tar == pot]
            tar_nc(pot, todo, manifest=do_manifest, messages= messages, errors=errors, verbose= verbose)
    else:
        if verbose:
            print('No archives created since do_tar == False')

