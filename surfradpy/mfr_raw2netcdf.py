import pathlib as pl
import site
import warnings
# import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.surfrad.surfrad as atmsrf
import xarray as xr
import pandas as pd
import surfradpy.database as sfp_db
import socket
import surfradpy.file_io.mfrsr as srpmfrsrio


def files_between(root: pl.Path, start: pd.Timestamp, end: pd.Timestamp, globpattern: str = ""):
    """ Generator that yields all files between start and end dates (inclusive) in the given root directory.
    Parameters
    ----------
    root : pl.Path
        Root directory containing year subdirectories with files.
    start : pd.Timestamp
        Start date.
    end : pd.Timestamp
        End date.
    globpattern : str, optional
        Glob pattern to match files. The default is mainly covering the extension, set to '*.nc' for netcdf files."".
    Yields
    -------
    pl.Path
        Paths to files between start and end dates.
    """ 
    assert(end > start), f'End must come after start! (end: {end}, start{start})'
    root = root
    d = start
    while d <= end:
        year_dir = root / f"{d.year}"
        assert(year_dir.exists())
        yield from year_dir.glob(f"*{d:%Y%m%d}{globpattern}")
        d += pd.to_timedelta(1, 'D')



class MfrsrRawToNetcdf:
    def __init__(self, 
                 path_in,
                 path_out,
                 name_pattern_netcdf, 
                 glob_pattern_raw = "*.xmd",
                 start = None,
                 end = None,
                 site = None,
                 version = '0.3',
                 path2surfrad_database = None,
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
            example 1:  In subfolders per year: '{year}/tbl_{serialnumber}_{year}{month}{day}.nc'
            example 2:  All in one folder: 'frc_{serialnumber}_{year}{month}{day}.nc'
        site: str, optional
            Thre letter site code (e.g., 'tbl'). Used for validation against the SURFRAD database.
        glob_pattern_raw: str, optional
            Glob pattern to find raw MFRSR files within `path_in`. Defaults to
            '*.xmd'. Consider using start and end in addition to this pattern.
        start: str or pd.Timestamp, optional
            Start date for processing. Have to provide end as well. glob_pattern_raw is still needed to define extension.
        end: str or pd.Timestamp, optional
            See start.
        version: str, optional
            Version of the conversion script. Defaults to '0.1'.
        path2surfrad_database: str, optional
            Path to the SURFRAD database file. If provided, additional metadata and tests are performed. 
            E.g. if serialnumber and head_id correspond. If the instrumment was deployed at that time at that site etc.
        reporter: object, optional
            An object to report progress or issues. Not yet implemented.
        **kwargs:
            Arbitrary keyword arguments that can be used as placeholders in
            `path_in`, `path_out`, and `name_pattern_netcdf`.

        Returns
        -------
        None.

        """ 
        kwargs['site'] = site
        for k,v in kwargs.items():
            setattr(self, k, v)
        kwargs['version'] = version
        self.kwargs = kwargs

        self.version = version
        self.path_in = pl.Path(path_in.format(**self.kwargs))
        self.path_out = pl.Path(path_out.format(**self.kwargs))
        self.path_out.mkdir(exist_ok=True, parents=True)

        if path2surfrad_database is None:
            path2surfrad_database = sfp_db.get_default_db_path()
        if isinstance(path2surfrad_database, (str, pl.Path)):
            self.surfrad_db = sfp_db.SurfradDatabase(path2surfrad_database)
        else:
            self.surfrad_db = None


        try:
            name_pattern_netcdf = name_pattern_netcdf.format(year = '{year}', month = '{month}', day = '{day}', **self.kwargs)
        except KeyError as e:
            raise KeyError(f'Key: "{e.args[0]}" not found in kwargs. Please provide it when initializing the class MFRSRRawToNetcdf as an additional kwarg, e.g. site = "tbl".')
        
        self.name_pattern_netcdf = name_pattern_netcdf
        self.glob_pattern_raw = glob_pattern_raw

        self._masterplan_in = None 
        self._masterplan_out = None
        self._workplan = None
        self.reporter = reporter
        self.freq='d'
        self.verbose = verbose
        self._processing_start = start
        self._processing_end = end

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
                # opt = atmsrf.read_raw(mp_in.iloc[i].p2f_in)
                opt = srpmfrsrio.open_rsr(mp_in.iloc[i].p2f_in)  # test if the file is readable
                break
            # except atmsrf.FileCorruptError:
            #     if self.verbose:
            #         print('corrupt, try next')
            #     i += 1
            #     continue
            except Exception as e:
                raise e
                if self.verbose:
                    print(f'error {e}, try next')
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
                # opt = atmsrf.read_raw(mp_in.iloc[-i].p2f_in)
                opt = srpmfrsrio.open_rsr(mp_in.iloc[-i].p2f_in)  # test if the file is readable
                break
            # except atmsrf.FileCorruptError:
            #     i += 1
            #     continue
            except Exception as e:
                raise e
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
        if isinstance(self._processing_start, type(None)):
            if self.verbose:
                print(f'Get all files in {self.path_in} with glob pattern: {self.glob_pattern_raw}')
            gen = self.path_in.glob(self.glob_pattern_raw)
        else:
            start = pd.to_datetime(self._processing_start)
            end = pd.to_datetime(self._processing_end) if not isinstance(self._processing_end, type(None)) else pd.Timestamp.now()
            if self.verbose:
                print(f'Get all files in {self.path_in} with "files_between" function and start: {start}, end: {end} and glob pattern: {self.glob_pattern_raw}')
            gen = files_between(self.path_in, start, end, globpattern = self.glob_pattern_raw)
        df  = pd.DataFrame(gen, columns=['p2f_in'])
        assert(df.shape[0]>0), f'No raw files found in {self.path_in} between {self._processing_start} and {self._processing_end} with pattern {self.glob_pattern_raw}.'
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
                if mp_in.iloc[-1].p2f_in == last_used_rawfile: # no new raw files
                    print('No new raw files to process, workplan is empty')
                    wp_in = mp_in.iloc[:0]
                    self._workplan = wp_in
                    return self._workplan
                    # assert(False), 'No new raw files, nothing to do here'
                    
                if ds.day_complete == 'True':
                    print('Last file was complete')
                    #when the day was complete we still want the last file since that file reached into the next day and was truncated
                    start_at_this_file = last_used_rawfile
                else:
                    print('Last file was incomplete')
                    # when the day was not complete we want to start from scratch with this file and load all files that has been used in the last netcdf file
                    start_at_this_file = ds.parent_files.split(',')[0].strip()
                    start_at_this_file = pl.Path(start_at_this_file)

                rawlable = mp_in.index[mp_in.p2f_in == start_at_this_file][0]
                pos = mp_in.index.get_loc(rawlable)
                wp_in = mp_in.iloc[pos:]

            self._workplan = wp_in
        return self._workplan

    def process(self, justone = False, save = True, id_mismatch_error = True, verbose = True):
        """ Process the workplan
        Parameters
        ----------
        justone : bool, optional
            If True, a single file will be processed but not saved. The default is False.
        save: bool, optional
            If True, the processed files will be saved. The default is True.
        id_mismatch_error : bool, optional
            If True, a ValueError will be raised if head_id or logger_id mismatches are found
            between raw files being concatonated. The default is True.
        verbose : bool, optional
            If True, print progress messages. The default is True.
            
        """
        # keep opening files until the start and end of the file are on a different day
        def open_next_row(wpiter, verbose = True):
            """ This opens the next readable (non-corrupt) file"""
            while 1:
                row_in = next(wpiter)[1]
                if verbose:
                    print(f'open file: {row_in.p2f_in}')
                try:
                    # dsin = atmsrf.read_raw(row_in.p2f_in)
                    dsin = srpmfrsrio.open_rsr(row_in.p2f_in)  # test if the file is readable
                    dsin.dataset # to avoid linting error
                    self.tp_p2f_in = row_in.p2f_in
                    self.tp_dsin = dsin
                    # print(dsin.filename, end = '\t')
                    # print(dsin.head.logger_id, flush=True)
                    break
                # except atmsrf.FileCorruptError:
                #     print(f'Corrupt file encountered: {row_in.p2f_in.as_posix()}.')
                #     continue
                except Exception as e:
                    print(f'Error {e} encountered when opening file: {row_in.p2f_in.as_posix()}. Try next file.')
                    if not isinstance(self.reporter, type(None)):
                       self.reporter.errors_increment()
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
            # lets check if all head_ids and logger_ids are identical
            head_ids = [ds.dataset.head_id for ds in ds_rawlist]
            if len(set(head_ids))!=1:
                msg = f'Head IDs in the raw files do not match: {head_ids}'
                if id_mismatch_error:
                    raise ValueError(msg)
                warnings.warn(msg)

            logger_ids = [ds.dataset.logger_id for ds in ds_rawlist]
            if len(set(logger_ids))!=1:
                msg = f'Logger IDs in the raw files do not match: {logger_ids}'
                if id_mismatch_error:
                    raise ValueError(msg)
                warnings.warn(msg)

            # lets concatonate and truncate them
            dsout = xr.concat([i.dataset for i in ds_rawlist], 'datetime')
            
            end = start_file + pd.to_timedelta(1,self.freq) - pd.to_timedelta(1, 'ns')

            dsout = dsout.drop_duplicates('datetime', keep = 'last') # in rare cases the data is present in multiple files
            dsout = dsout.sel(datetime = slice(start_file,end))

            # check metadata consistency inculding database information
            ## get metadata from surfrad database
            query = f'SELECT * FROM instruments_mfrsr WHERE Logger_ID="${dsout.logger_id}"'
            db_meta = self.surfrad_db.execute_query(query)
            assert(db_meta.shape[0] == 1), f'No unique entry found in instruments_mfrsr for Logger_ID = {dsout.logger_id}. Found {db_meta.shape[0]} entries.'
            db_meta = db_meta.iloc[0]

            ## where was the instrument deployed at the time of interest

            query = f'SELECT * FROM deployments WHERE instrument_sn = "{db_meta.Instrument}"'
            dep_df = self.surfrad_db.execute_query(query)
            dep_df['Date_start'] = pd.to_datetime(dep_df.Date_start)
            dep_df = dep_df.sort_values('Date_start')
            dep_df = dep_df[dep_df.Date_start <= start_file ] # remove deployments after the file date
            assert(dep_df.shape[0]>0), f'No deployment found for instrument {db_meta.Instrument} before {start_file}.'
            deployment = dep_df.iloc[-1]     #last_install_before_file_date

            self.tp_deployment = deployment
            assert(deployment.Date_stop == 'present' or pd.to_datetime(deployment.Date_stop) >= start_file), f'Instrument {db_meta.Instrument} was not deployed at {start_file}. Deployment ended at {deployment.Date_stop}.'

            ## Double check that the deployment site matches the expected site

            assert(deployment.Location.lower() == self.site.lower()), f'Instrument {db_meta.Instrument} was not deployed at site {self.site} at {start_file}, but at {deployment.Location}.'

            ## get site metadata
            query = f'SELECT * FROM sites WHERE abb="{self.site}"'
            site_meta = self.surfrad_db.execute_query(query).iloc[0]

            self.tp_ds_rawlist = ds_rawlist
            self.tp_db_meta = db_meta
            self.tp_site_meta = site_meta

            # add attributes
            self.tp_ri = ds_rawlist[0]
            dsraw = ds_rawlist[0].dataset
            attrs = dsraw.attrs.copy()
            attrs['info'] = 'SURFRAD MFRSR raw data converted to netcdf format, concatonated and truncated to daily files in UTC time.'
            attrs['source'] = 'This netcdf file was created by surfradpy.mfr_raw2netcdf.MfrsrRawToNetcdf'
            attrs['processing_date'] = pd.Timestamp.now().isoformat()
            attrs['processing_server'] = socket.gethostname()
            attrs['product_version'] = self.version
            attrs['site'] = site_meta.abb
            attrs['site_name'] = site_meta['name']
            attrs['elevation'] = site_meta.elevation
            attrs['latitude'] = site_meta.latitude
            attrs['longitude'] = site_meta.longitude
            # attrs['band_on'] = dsraw.attrs['band_on']
            # attrs['band_on'] = dsraw.attrs['band_on']
            # attrs['avg_period'] = dsraw.attrs['avg_period']
            # attrs['sample_rate'] = dsraw.attrs['sample_rate']
            # attrs['instrument_type'] = dsraw.attrs['instrument_type']
            attrs['sn_mfrsr'] = db_meta.Instrument
            attrs['logger_id'] = db_meta.Logger_ID
            attrs['head_id'] = db_meta.Head_ID
            attrs['parent_files'] = ', '.join([i.dataset.attrs['path2file'] for i in ds_rawlist])
            attrs['day_complete'] = complete.__str__()
            dsout.attrs = attrs
            # create the pathname and save under that name
            
            dt = start_file
            p2f_out = self.path_out.joinpath(self.name_pattern_netcdf.format(year = f'{dt.year:04d}', month = f'{dt.month:02d}', day = f'{dt.day:02d}'))

            self.tp_p2fout = p2f_out
            self.tp_dsout = dsout
            self.tp_start_file = start_file
            # assert(False), 'debug stop'
            if justone:
                if verbose:
                    print('justone active, skip saving')
                return {'dsout': dsout, 'p2f_out': p2f_out, 'active': active, 'start_file': start_file, 'complete': complete}
            if dsout.datetime.shape[0] == 0:
                if verbose:
                    print('no data in dataset, skip saving')
            else:
                if p2f_out.is_file():
                    if verbose:
                        print('File exists, skip saving')
                elif not save:
                    # if verbose:
                    #     print('save is False, skip saving')
                    pass
                else:
                    p2f_out.parent.mkdir(parents=True, exist_ok=True)
                    dsout.to_netcdf(p2f_out)
                    if not isinstance(self.reporter, type(None)):
                        self.reporter.clean_increment()
                        self.reporter.log()
            return {'dsout': dsout, 'p2f_out': p2f_out, 'active': active, 'start_file': start_file, 'complete': complete}
        # Keep processing
        wp_in = self.workplan
        if wp_in.shape[0] == 0:
            print('Workplan is empty, nothing to do')
            return
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
            if justone:
                print('justone active, stop after first file')
                break
            active = out['active']
            complete = out['complete']
            start_file = out['start_file'] + pd.to_timedelta(1, 'd')
            dsout = out['dsout']
            i += 1
        return out


