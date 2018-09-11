import string as _string
import pandas as _pd
import numpy as _np
import xarray as _xr
import os as _os
import tarfile as _tarfile
from .surfrad import locations as _locations
import platform as _platform
import hashlib as _hashlib

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

# location_values = {'latitude': 40.05192,
#                    'longitude': 271.62692,
#                    'altitude': 230.}


def filename2info(filename):
    abbr, rest = filename.split('_')
    date = rest.split('.')[0]
    date_pd = _pd.to_datetime(date)
    return abbr, date_pd

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


# def generate_dataset(time_dimension=480):
#     template_file = '../data/template_{}.nc'.format(time_dimension)
#     ds = _xr.open_dataset(template_file, autoclose=True,
#                           chunks={'time': time_dimension}
#                           )
#
#     for var in ds.variables:
#         enc = ds[var].encoding
#         enc['chunksizes'] = enc['original_shape']
#         enc['contiguous'] = False
#         enc['complevel'] = 9
#         enc['zlib'] = True
#     return ds


# def populate_dataset(dataset, data):
#     """
#     TODO
#     ----
#     outsource the testing!
#     """
#
#     ds = dataset
#     ds['time'].values = data.index
#
#     invalid = _string.punctuation
#     for i in {'_', '-', '.', '+', '@'}:
#         invalid = invalid.replace(i, '')
#
#     e = 1
#     verbose = True
#     for varname, colname in columntransrules.items():
#         #     print(varname)
#         #     break
#         ds[varname].values = data.loc[:, colname].values.astype(ds[varname].dtype)
#         if 'flag_values' in ds[varname].attrs:
#             flag_values_len = len(ds[varname].attrs['flag_values'])
#             try:
#                 flag_meaning_len = len(ds[varname].attrs['flag_meanings'].split())
#             except KeyError:
#                 print('-----')
#                 print('{}) {}: no flag_meanings found'.format(e, varname))
#                 e += 1
#                 if verbose:
#                     print('\tflag_values: {}'.format(ds[varname].attrs['flag_values']))
#                 continue
#             if flag_values_len != flag_meaning_len:
#                 print('-----')
#                 print(('{}) {}: number of flag_values ({}) different than number',
#                        ' of flag_meanings ({}).'.format(e, varname, flag_values_len, flag_meaning_len)))
#                 e += 1
#                 if verbose:
#                     print('\tflag_values: {}'.format(ds[varname].attrs['flag_values']))
#                     print('\tflag_meanings: {}'.format(ds[varname].attrs['flag_meanings']))
#
#             invinstring = []
#             for char in ds[varname].attrs['flag_meanings']:
#                 if char in invalid:
#                     invinstring.append(char)
#
#                     #                 print('------')
#                     #                 print(ds['QC_flag_07_global_shortwave_over_sum'].attrs['flag_meanings'])
#                     #                 assert False
#
#             invinstring = set(invinstring)
#             #         if any(char in invalid for char in ds[varname].attrs['flag_meanings']):
#             if len(invinstring) > 0:
#                 print('-----')
#                 print('{}) {}: invalid character(s) {} found in flag_meaning'.format(e, varname, ''.join(invinstring)))
#                 e += 1
#                 if verbose:
#                     print('\tflag_meanings: {}'.format(ds[varname].attrs['flag_meanings']))
#
#         if 'ancillary_variables' in ds[varname].attrs:
#             #         raise ValueError()
#             for avar in ds[varname].attrs['ancillary_variables'].split():
#                 if avar not in ds.variables:
#                     print('-----')
#                     print(('{}) {}: the ancillary_variables {} is not a valid variable',
#                            '... spell check.'.format(e, varname, avar)))
#                     e += 1
#
#     for varname, value in location_values.items():
#         ds[varname].values = value


# def save_dataset2netcdf(dataset: _xr.Dataset, fname: str) -> bool:
#     ds = dataset
#     ds.to_netcdf(fname,
#                  format='NETCDF4_CLASSIC'
#                  )
#     return True


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

def save2netcdf(ds, fname = '../data/Bondville_IL_1995_May_10.new.nc', compression = False):
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

    ## save
    ds.to_netcdf(fname,
                 format = 'NETCDF4_CLASSIC',
                )


def qcrad2netcdf(fname, output_folder, path2CDL, overwrite=False, verbose=False):
    # folder_in = '/Volumes/HTelg_4TB_Backup/SURFRAD/qcrad_v3/bon/1995/'
    # folder_out = '/Volumes/HTelg_4TB_Backup/SURFRAD/NCEI/'
    # overwrite = False
    # verbose = True

    cdl_dict = parse_CDL_file(path2CDL)
    folder_in = fname
    folder_out = output_folder

    if _os.path.isfile(fname):
        folder_in, fn = _os.path.split(folder_in)
        folder_in += '/'
        fnames = [fn]
    elif _os.path.isdir(fname):
        fnames = _os.listdir(folder_in)
    else:
        raise ValueError('fname is neither a file nore a folder (fname = {})'.format(fname))
    # clean files
    fnames = [i for i in fnames if '.qdat' in i]

    # find files that need processing
    process_list = []
    for fn in fnames:
        # generate output folder and name
        ## get info from file name
        abbr, date_pd = filename2info(fn)

        # find location based on abbriviation
        # abbr = 'bon'
        loc = [i for i in _locations if abbr in i['abbriviations']][0]

        # generate file name
        file_out = '{name}_{state}_{year}_{month}_{day:02d}.nc'.format(name=loc['name'], state=loc['state'],
                                                                       year=date_pd.year,
                                                                       month=date_pd.month_name()[:3], day=date_pd.day)

        # generate folder name
        ## check if final folder is already the year
        if folder_out.strip('/').split('/')[-1] != str(date_pd.year):
            folder_out = folder_out + '{}/'.format(date_pd.year)

        # check if folder exists, if not create it
        if not _os.path.isdir(folder_out):
            _os.makedirs(folder_out)

        path_out = folder_out + file_out

        if not overwrite:
            if _os.path.isfile(path_out):
                continue
        out = {'path_in': folder_in + fn, 'path_out': path_out, 'location':loc}
        process_list.append(out)
    # print list of files that need to be processed
    if 0:
        print('\n'.join([i['path_in'] for i in process_list]))

    # process and save all designated files
    if len(process_list) == 0:
        print('process_list is empty ... nothing to do here')
    for todo in process_list:
        if verbose:
            print('processing {}\n ......'.format(todo['path_in']), end=' ')
        data = read_data(todo['path_in'])
        ds = df2ds(data, cdl_dict['variable_list'], cdl_dict['global_atts_dict'], todo['location'])
        if save2netcdf(ds, todo['path_out']):
            if verbose:
                print('done! netCDF saved to {}'.format(todo['path_out']))

def generate_md5_checksum(fname):
    with open(fname, "rb") as file_to_check:
        data = file_to_check.read()
        md5 = _hashlib.md5(data).hexdigest()

    f_size = _os.stat(fname).st_size

    fname_out = fname + '.mnf'

    file_content = [_os.path.split(fname)[-1], str(md5), str(f_size)]

    with open(fname_out, 'w') as fout:
        fout.write(','.join(file_content))

def tar_netcdf_files(input_folder, output_folder, compression=True, overwrite = False, test = True, manifest = True, verbose = False):
    # some definitions
    tar_mode = 'w:'
    if compression:
        f_ext = '_tar_gz'
        fl_ext = '.tar.gz'
        tar_mode += 'gz'
    else:
        f_ext = '_tar'
        fl_ext = '.tar'

    # generate output_folder name and create the folder if needed
    fnt = input_folder.strip('/').split('/')[-1]
    # year = fnt
    fnt += f_ext
    if output_folder.strip('/').split('/')[-1] != fnt:
        output_folder = output_folder + fnt + '/'

    if not _os.path.isdir(output_folder):
        _os.makedirs(output_folder)

    # get file names and provide attribute of month for grouping
    fnames = _os.listdir(input_folder)
    fnames = [fn for fn in fnames if '.nc' in fn]

    dft = _pd.DataFrame(_np.array([fnames]).transpose(), columns=['filename'])
    # df['month'] = [fn.split('_')[-2] for fn in df.filename]
    # df['year'] = [int(fn.split('_')[-3]) for fn in df.filename]

    df = dft.filename.apply(lambda x: _pd.Series(x.split('_'), index=['filename', 'state', 'year', 'month', 'day']))
    df['filename'] = dft.filename
    df['day'] = df.day.apply(lambda x: x.split('.')[0])
    df.index = df.apply(lambda x: _pd.to_datetime('{}-{}-{}'.format(x.year, x.month, x.day)), axis=1)
    df.sort_index(inplace=True)

    # remove all files from the last month -> a tar is only created if the next month has started
    last = df.index[-1]
    last = last - _pd.to_timedelta(last.day, 'd')
    df = df.truncate(after=last)

    # return df
    # write a tar of monthly files
    for month in _np.unique(df.month):
        fns = df[df.month == month].filename.values
        name, state, year, month = fns[0].split('_')[:-1]
        archive_name = '{}_{}_{}_{}{}'.format(name, state, year, month, fl_ext)
        if test:
            print('============')
            print('archive_name: {}'.format(output_folder + archive_name))
            for fn in fns:
                print('\t {}'.format(fn))
        else:
            archive_path = output_folder + archive_name
            if not overwrite:
                if _os.path.isfile(archive_path):
                    if verbose:
                        print('File {} already exists ... go to next'.format(archive_path))
                    continue
            tar = _tarfile.open(archive_path, mode=tar_mode)

            for fn in fns:
                tar.add(input_folder + fn, arcname = fn)

            tar.close()
            if manifest:
                generate_md5_checksum(archive_path)
    return df