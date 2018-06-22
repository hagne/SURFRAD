import string as _string
import pandas as _pd
import numpy as _np
import xarray as _xr
import os as _os
import tarfile as _tarfile
from .surfrad import locations as _locations

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

location_values = {'latitude': 40.05192,
                   'longitude': 271.62692,
                   'altitude': 230.}


def filename2info(filename):
    abbr, rest = filename.split('_')
    date = rest.split('.')[0]
    date_pd = _pd.to_datetime(date)
    return abbr, date_pd


def read_data(fname, missing_data=((-9999.0, -9999.9), _np.nan), pressure_mbar2pa=True):
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

    return data


def generate_dataset(time_dimension=480):
    template_file = '../data/template_{}.nc'.format(time_dimension)
    ds = _xr.open_dataset(template_file, autoclose=True,
                          chunks={'time': time_dimension}
                          )

    for var in ds.variables:
        enc = ds[var].encoding
        enc['chunksizes'] = enc['original_shape']
        enc['contiguous'] = False
        enc['complevel'] = 9
        enc['zlib'] = True
    return ds


def populate_dataset(dataset, data):
    """
    TODO
    ----
    outsource the testing!
    """

    ds = dataset
    ds['time'].values = data.index

    invalid = _string.punctuation
    for i in {'_', '-', '.', '+', '@'}:
        invalid = invalid.replace(i, '')

    e = 1
    verbose = True
    for varname, colname in columntransrules.items():
        #     print(varname)
        #     break
        ds[varname].values = data.loc[:, colname].values.astype(ds[varname].dtype)
        if 'flag_values' in ds[varname].attrs:
            flag_values_len = len(ds[varname].attrs['flag_values'])
            try:
                flag_meaning_len = len(ds[varname].attrs['flag_meanings'].split())
            except KeyError:
                print('-----')
                print('{}) {}: no flag_meanings found'.format(e, varname))
                e += 1
                if verbose:
                    print('\tflag_values: {}'.format(ds[varname].attrs['flag_values']))
                continue
            if flag_values_len != flag_meaning_len:
                print('-----')
                print(('{}) {}: number of flag_values ({}) different than number',
                       ' of flag_meanings ({}).'.format(e, varname, flag_values_len, flag_meaning_len)))
                e += 1
                if verbose:
                    print('\tflag_values: {}'.format(ds[varname].attrs['flag_values']))
                    print('\tflag_meanings: {}'.format(ds[varname].attrs['flag_meanings']))

            invinstring = []
            for char in ds[varname].attrs['flag_meanings']:
                if char in invalid:
                    invinstring.append(char)

                    #                 print('------')
                    #                 print(ds['QC_flag_07_global_shortwave_over_sum'].attrs['flag_meanings'])
                    #                 assert False

            invinstring = set(invinstring)
            #         if any(char in invalid for char in ds[varname].attrs['flag_meanings']):
            if len(invinstring) > 0:
                print('-----')
                print('{}) {}: invalid character(s) {} found in flag_meaning'.format(e, varname, ''.join(invinstring)))
                e += 1
                if verbose:
                    print('\tflag_meanings: {}'.format(ds[varname].attrs['flag_meanings']))

        if 'ancillary_variables' in ds[varname].attrs:
            #         raise ValueError()
            for avar in ds[varname].attrs['ancillary_variables'].split():
                if avar not in ds.variables:
                    print('-----')
                    print(('{}) {}: the ancillary_variables {} is not a valid variable',
                           '... spell check.'.format(e, varname, avar)))
                    e += 1

    for varname, value in location_values.items():
        ds[varname].values = value


def save_dataset2netcdf(dataset: _xr.Dataset, fname: str) -> bool:
    ds = dataset
    ds.to_netcdf(fname,
                 format='NETCDF4_CLASSIC'
                 )
    return True


def qdat2netcdf(fname, output_folder, overwrite=False, verbose=False):
    # folder_in = '/Volumes/HTelg_4TB_Backup/SURFRAD/qcrad_v3/bon/1995/'
    # folder_out = '/Volumes/HTelg_4TB_Backup/SURFRAD/NCEI/'
    # overwrite = False
    # verbose = True

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
        out = {'path_in': folder_in + fn, 'path_out': path_out}
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
        ds = generate_dataset(time_dimension=data.shape[0])
        populate_dataset(ds, data)
        if save_dataset2netcdf(ds, todo['path_out']):
            if verbose:
                print('done! netCDF saved to {}'.format(todo['path_out']))


def tar_netcdf_files(input_folder, output_folder, compression=True):
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
    year = fnt
    fnt += f_ext
    if output_folder.strip('/').split('/')[-1] != fnt:
        output_folder = output_folder + fnt + '/'

    if not _os.path.isdir(output_folder):
        _os.makedirs(output_folder)

    # get file names and provide attribute of month for grouping
    fnames = _os.listdir(input_folder)
    fnames = [fn for fn in fnames if '.nc' in fn]

    df = _pd.DataFrame(_np.array([fnames]).transpose(), columns=['filename'])
    df['month'] = [fn.split('_')[-2] for fn in df.filename]

    # write a tar of monthly files
    for month in _np.unique(df.month):
        fns = df[df.month == month].filename.values
        name, state, year, month = fns[0].split('_')[:-1]
        archive_name = '{}_{}_{}_{}{}'.format(name, state, year, month, fl_ext)

        tar = _tarfile.open(output_folder + archive_name, mode=tar_mode)

        for fn in fns:
            tar.add(input_folder + fn, arcname = fn)

        tar.close()
