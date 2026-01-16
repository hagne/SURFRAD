import pathlib as pl
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.surfrad.surfrad as atmsrf
import xarray as xr
import pandas as pd


class MfrsrRawToNetcdf:
    def __init__(self, 
                 path_in,
                 path_out,
                 name_pattern_netcdf, 
                 glob_pattern_raw = "**/*.xmd",
                 version = '0.1',
                 reporter = None,
                 verbose = True,
                 **kwargs,
                 ):
        """
        Conversion class which converts raw MFRSR files to daily NetCDF files.

        Parameters
        ----------
        path_in: str
            Path to the raw MFRSR files. Can contain placeholders for kwargs
            (e.g., '{serialnumber}').
            example: '/nfs/grad/Inst/MFR/Campaign/frc/2025/mfrsr/{serialnumber}',
        path_out: str
            Path where the daily NetCDF files will be stored. Can contain
            placeholders for kwargs (e.g., '{serialnumber}', '{version}').
            example: '/nfs/grad/Inst/MFR/Campaign/frc/2025/mfrsr/{serialnumber}.netcdf/v{version}/'
        name_pattern_netcdf: str
            Filename pattern for the output NetCDF files. Can contain
            placeholders for year, month, day, and kwargs
            (e.g., 'frc_{serialnumber}_{year}{month}{day}.nc').
            example:  'frc_{serialnumber}_{year}{month}{day}.nc'
        glob_pattern_raw: str, optional
            Glob pattern to find raw MFRSR files within `path_in`. Defaults to
            '**/*.xmd'.
        version: str, optional
            Version of the conversion script. Defaults to '0.1'.
        reporter: object, optional
            An object to report progress or issues. Not yet implemented.
        **kwargs:
            Arbitrary keyword arguments that can be used as placeholders in
            `path_in`, `path_out`, and `name_pattern_netcdf`.

        Returns
        -------
        None.

        """

        for k,v in kwargs.items():
            setattr(self, k, v)
        kwargs['version'] = version
        self.kwargs = kwargs

        self.version = version
        self.path_in = pl.Path(path_in.format(**self.kwargs))
        self.path_out = pl.Path(path_out.format(**self.kwargs))
        self.path_out.mkdir(exist_ok=True, parents=True)

        name_pattern_netcdf = name_pattern_netcdf.format(year = '{year}', month = '{month}', day = '{day}', **self.kwargs)
        self.name_pattern_netcdf = name_pattern_netcdf
        self.glob_pattern_raw = glob_pattern_raw

        self._masterplan_in = None 
        self._masterplan_out = None
        self._workplan = None
        self.reporter = reporter
        self.freq='d'
        self.verbose = verbose

    def _make_masterplan_out(self):
        # get all the files that should at the end exist
        # load very first file to get the earliest datapoint
        i=0
        mp_in = self._masterplan_in
        while 1:
            # print(i)
            if self.verbose:
                print(f'opening file {mp_in.iloc[i].p2f_in}', end = ' ... ')
            if i == 10:
                print('lets not do this')
                break
            try:
                opt = atmsrf.read_raw(mp_in.iloc[i].p2f_in)
                break
            except atmsrf.FileCorruptError:
                if self.verbose:
                    print('corrupt, try next')
                i += 1
                continue
        if self.verbose:
            print('done') 

        start = pd.to_datetime(opt.dataset.datetime.values[0])
        
        # Load the very last file to get the last datapoint
        i=1
        while 1:
            print(i)
            if i == 10:
                print('lets not do this')
                break
            try:
                opt = atmsrf.read_raw(mp_in.iloc[-i].p2f_in)
                break
            except atmsrf.FileCorruptError:
                i += 1
                continue

        end = pd.to_datetime(opt.dataset.datetime.values[-1])

        # data frame with dayly output files
        df = pd.DataFrame(index = pd.date_range(start.date(),end.date(), freq=self.freq), columns = ['p2out'])
        #TODO: below the pattern need to be defined somewhere
        df['p2out'] = df.apply(lambda row: self.path_out.joinpath(self.name_pattern_netcdf.format(year = f'{row.name.year:04d}',month = f'{row.name.month:02d}', day = f'{row.name.day:02d}')), axis = 1)
        return df
    
    def _make_masterplan_in(self):
        # Get all raw files
        df  = pd.DataFrame(self.path_in.glob(self.glob_pattern_raw), columns=['p2f_in'])
        df = df.sort_values('p2f_in')
        df['fname'] = df.apply(lambda row: row.p2f_in.name, axis = 1)
        return df

    @property
    def masterplan(self):
        if isinstance(self._masterplan_out, type(None)):
            self._masterplan_in = self._make_masterplan_in()
            self._masterplan_out = self._make_masterplan_out()
        return dict(raw_files = self._masterplan_in, processed_files = self._masterplan_out)
    
    @property
    def workplan(self):
        # what do we actually need to work on
        if isinstance(self._workplan, type(None)):
            mp_out = self.masterplan['processed_files']
            mp_in = self.masterplan['raw_files']
            mp_out_exist = mp_out[mp_out.apply(lambda row: row.p2out.is_file(), axis = 1)]

            if mp_out_exist.shape[0] == 0:
                print('No rawfiles have been processe yet, start from the beginning')
                wp_in = mp_in
            else:
                lastrow = mp_out_exist.iloc[-1]
                ds = xr.open_dataset(lastrow.p2out)

                # check if any new raw files have been produced, if not there is nothing to do
                last_used_rawfile = ds.parent_files.split(',')[-1].strip()
                last_used_rawfile = pl.Path(last_used_rawfile)
                if mp_in.iloc[-1].p2f_in == last_used_rawfile:
                    assert(False), 'No new raw files, nothing to do here'
                    
                if ds.day_complete == 'True':
                    print('Last file was complet')
                    #when the day was complet we still want the last file since that file reached into the next day and was truncated
                    start_at_this_file = last_used_rawfile
                else:
                    print('Last file was incomplet')
                    # when the day was not complete we want to start from scratch with this file and load all files that has been used in the last netcdf file
                    start_at_this_file = ds.parent_files.split(',')[0].strip()
                    start_at_this_file = pl.Path(start_at_this_file)

                rawlable = mp_in.index[mp_in.p2f_in == start_at_this_file][0]
                pos = mp_in.index.get_loc(rawlable)
                wp_in = mp_in.iloc[pos:]

            self._workplan = wp_in
        return self._workplan
    
    def process(self):
        # keep opening files until the start and end of the file are on a different day
        def open_next_row(wpiter, verbose = True):
            """ This opens the next readable (non-corrupt) file"""
            while 1:
                row_in = next(wpiter)[1]
                if verbose:
                    print(f'open file: {row_in.p2f_in}')
                try:
                    dsin = atmsrf.read_raw(row_in.p2f_in)
                    break
                except atmsrf.FileCorruptError:
                    print(f'Corrupt file encountered: {row_in.p2f_in.as_posix()}.')
                    continue
            return dsin 

        # get the start and end of the file, keep previous start if this is not the first file.
        def daysfromstart(dsin, start_file, whichend = 'end'):
            if whichend == 'end':
                test = pd.to_datetime(dsin.dataset.datetime.values[-1])
            elif whichend == 'start':
                test  = pd.to_datetime(dsin.dataset.datetime.values[0])
            return (test.date() - start_file.date()) / pd.to_timedelta('1d')

        def make_file_of_the_day(wpiter, active, start_file, verbose = True):
            ds_rawlist = []
            complete = True
            
            # open a file, but only if none is active from the last loop
            if not active:
                active = open_next_row(wpiter, verbose=verbose)
            
            ds_rawlist.append(active)
            
            # in first loop only
            if not start_file:
                start_file = pd.to_datetime(active.dataset.datetime.values[0].astype('datetime64[D]'))
            
            # check if the file ends on the same day. If so we want to check the next one and see if ends on the next day, while starting on the same day !!!
            noofdays = daysfromstart(active, start_file)
            
            while noofdays == 0:
                if verbose:
                    print('file ends the same day, get the next one')
                try:
                    active = open_next_row(wpiter, verbose = verbose)
                except StopIteration:
                    complete = False
                    if verbose:
                        print(f'Reached end of file list. This file will be incomplet')
                    break
                # Removed the blow: I don't thing it is that concerning when the file is a little longer, it gets concatinated anyway?
                # when the file starts on the next day it should not be included in the list
                # if daysfromstart(active, start_file, whichend='start') == 0:
                
                ds_rawlist.append(active)
                noofdays = daysfromstart(active, start_file)
                
            if verbose:     
                print('more than one day in files, concat and truncate')
            
            # ds_rawlist should now include all files that have data on this day
            # lets concatonate and truncate them
            dsout = xr.concat([i.dataset for i in ds_rawlist], 'datetime')
            
            end = start_file + pd.to_timedelta(1,self.freq) - pd.to_timedelta(1, 'ns')

            dsout = dsout.drop_duplicates('datetime', keep = 'last') # in rare cases the data is present in multiple files
            dsout = dsout.sel(datetime = slice(start_file,end))
            
            dsout.attrs = ds_rawlist[0].dataset.attrs
            dsout.attrs.pop('path2file')
            dsout.attrs['day_complete'] = complete.__str__()
            dsout.attrs['parent_files'] = ', '.join([i.dataset.attrs['path2file'].as_posix() for i in ds_rawlist])
            dsout.attrs['product_version'] = self.version
            dsout.attrs['info'] = 'Original raw files concatinated/truncated to daily files and converted to netcdf files'
            
            # create the pathname and save under that name
            
            dt = start_file
            p2f_out = self.path_out.joinpath(self.name_pattern_netcdf.format(year = f'{dt.year:04d}', month = f'{dt.month:02d}', day = f'{dt.day:02d}'))
            
            if dsout.datetime.shape[0] == 0:
                if verbose:
                    print('no data in dataset, skip saving')
            else:
                if p2f_out.is_file():
                    if verbose:
                        print('File exists, skip saving')
                else:
                    dsout.to_netcdf(p2f_out)
            return {'dsout': dsout, 'p2f_out': p2f_out, 'active': active, 'start_file': start_file, 'complete': complete}
        # Keep rocessing
        wp_in = self.workplan
        wpiter = wp_in.iterrows()
        active = False
        start_file = False
        complete = True
        # endsonsameday = False

        i = 0
        while complete:
            # print(f'make next file of the day ({i})')
            print('.', end = '')
            out = make_file_of_the_day(wpiter, active, start_file, verbose = False)
            active = out['active']
            complete = out['complete']
            start_file = out['start_file'] + pd.to_timedelta(1, 'd')
            dsout = out['dsout']
            i += 1
        return dsout



