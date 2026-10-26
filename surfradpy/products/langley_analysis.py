
import xarray as xr
import productomator.worker as prowo
import atmPy.radiation.retrievals.spectral_irradiance as atmspec
import atmPy.radiation.retrievals.langley_db as langley_db
import sqlite3
import pandas as pd
import atmPy.radiation.retrievals.langley_calibration as atmlc
import numpy as np
import pathlib as pl

class LangleyAnalysis(prowo.Workplanner):
    """
    changes
    -------
    0.1 - as PMOD-FRC 2025
    """

    def __init__(self,
                 p2fld_in,
                 p2database=None, #('test.db', 'langley_results', 'local_day', 'input_file'),
                #  database=('test.db', 'langley_results', 'local_day', 'input_file'),
                 file_name_format = '*{date:%Y%m%d}*',
                 start = None,
                 end = None,
                 days = None,
                 input_directory_structure='yearly',
                 reporter = None,
                 verbose = False,
                 langley_airmass_limits = (2, 6),
                 expected_instrument_noise = 0.001,
                 langley_clean_threshold = 2,
                 ):
        
        self.version = 0.1
        if p2database is not None:
            database = (p2database.format(version=self.version), 'langley_results', 'local_day', 'input_file')
            if verbose:
                print(f'Using database: {database}')
        super().__init__(
            p2fld_in = p2fld_in,
            supplemental_data=None, #{'name': 'ozon', 'path2folder': p2f_in_sup, 'file_name_format': '*{date:%Y%m%d}*','input_directory_structure': 'yearly'},
            p2fld_out=None,
            # database=(srfdb.get_default_db_path(), 'mfrsr_calibration_parameters', 'timestamp', 'input_file'),
            database=database,
            # date_from_name=lambda name: pd.to_datetime(name.split('_')[-1].strip('.nc')),
            file_name_format = file_name_format,
            output_file_format=None,
            # glob_pattern_in='*.nc',
            start = start,
            end= end,
            days = days,
            input_directory_structure=input_directory_structure,
            output_directory_structure=None,
            file_complete_check=False,
            reporter=reporter,
            verbose= verbose,
            )
        self.langley_airmass_limits = langley_airmass_limits
        self.expected_instrument_noise = expected_instrument_noise
        self.langley_clean_threshold = langley_clean_threshold

        return 
    
    def process_row(self, row = None, iloc = None, loc = None, save = True):
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
        out = {}
        if iloc is not None:
            row = self.workplan.iloc[iloc]
        elif loc is not None:
            row = self.workplan.loc[loc]
        self.tp_row = row

        #####
        # get last processed instance - usefull if processing depends on the previous day
        ######
        if self.verbose:
            print(f'Opening file: {row.p2f_in}')
        indict = self.open_local_day(row, lambda x: xr.open_dataset(x), cut_nighttime=True) # this will make sure you have a full local day, it will open the next day if needed to finish the current day
        if indict['status'] == "Day not done yet":
            return out
        elif indict['status'] == "File has no data":
            if save:
                db = langley_db.LangleyResultDatabase(self.database[0], version = self.version)
                for when in ['am', 'pm']:
                    for channel in ['415', '500', '670', '870', '940', '1625']:
                        db.add_empty_entry(row.name, path2file=row.p2f_in, when = when, channel = channel)
            return out
        elif indict['dataset'].datetime.shape[0] == 0: # the difference between this case and the above is that here the file is empty as a result of cutting it to the local day, while in the above case the file is empty frome the start
            if save:
                db = langley_db.LangleyResultDatabase(self.database[0], version = self.version)
                for when in ['am', 'pm']:
                    for channel in ['415', '500', '670', '870', '940', '1625']:
                        db.add_empty_entry(row.name, path2file=row.p2f_in, when = when, channel = channel)
            return out
        
        ds = indict['dataset']
        self.tp_indict = indict
        ## Do some processing here, e.g. add attributes, format the dataset, etc.
        #####

        si = atmspec.CombinedGlobalDiffuseDirect(ds)
        sir = si.direct_normal_irradiation
        sir.settings_langley_airmass_limits = self.langley_airmass_limits    



        clean_channels_individually = True 
        lang = sir.langley_am
        lc = lang.clean(threshold=self.langley_clean_threshold,
                        clean_channels_individually=clean_channels_individually, 
                        expected_instrument_noise=self.expected_instrument_noise,
                        verbose=self.verbose)
        langc_am = lc['langley']
        # print(lang.tp_cutoff)
        
        
        lang = sir.langley_pm
        lc = lang.clean(threshold=self.langley_clean_threshold, clean_channels_individually=clean_channels_individually,
                        expected_instrument_noise=self.expected_instrument_noise, verbose=self.verbose)
        langc_pm = lc['langley']
        if save:
            db = langc_am.save2database(self.database[0], path2file=row.p2f_in, version = self.version)
            db = langc_pm.save2database(self.database[0], path2file=row.p2f_in, version = self.version)
        else:
            db = None

        ds.close()
        out['dni_instance'] = sir
        out['lang_am'] = langc_am
        out['lang_pm'] = langc_pm
        out['database'] = db
        return out

class V0prediction(prowo.Workplanner):
    def __init__(self,
                 path2database_in, #currently this is site specific, so the site kwarg has no effect on this
                 database=('test.db', 'gam_V0_prediction','date',"None"),
                 start=None, #pd.to_datetime(pd.Timestamp.now().date()) - pd.to_timedelta(80, 'D'),
                 end=None, #pd.Timestamp.now(),
                 days=None,
                 reporter=None,
                 verbose=True,
                 site = 'tmp', #only needed if multiple sites will be saved in the database
                 ):
        super().__init__(
            database=database,
            start=start,
            end=end,
            days=days,
            reporter=reporter,
            verbose=verbose,
        )
        self.path2database_in = path2database_in
        self.site = site
        self.langley_analysis_version = '0.1'
        return

    def process_row(self, row = None, iloc = None, loc = None, save = True):
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
        out = {}
        if iloc is not None:
            row = self.workplan.iloc[iloc]
        elif loc is not None:
            row = self.workplan.loc[loc]
        self.tp_row = row
        
        
        #######
        ## test if input database contains the right data, this should probably go into the workplan property? Unless this is causing en error anyway?!?
        #######
        with sqlite3.connect(self.path2database_in) as con:
            db_max = pd.read_sql("SELECT MAX(local_day) FROM langley_results;", con)
            db_min = pd.read_sql("SELECT MIN(local_day) FROM langley_results;", con)
        
        db_max = pd.to_datetime(db_max.iloc[0,0])
        db_min = pd.to_datetime(db_min.iloc[0,0])
        
        # not sure if this should be prevented earlier? basically when the retrieval is called in the first place?!?
        daystodbend = ((db_max - row.name) / pd.to_timedelta(1, 'D')) 
        assert(daystodbend> 60), f'Final run should only be run 2 month later, last database entry is only {daystodbend} days old.'
        
        # determine the range to use for the gam analysis, the range that is taken from the input database
        daysingamanalysis = 3*356
        
        if daystodbend < (daysingamanalysis/2):
            end = db_max
            start = db_max - pd.to_timedelta(daysingamanalysis, 'D')
        else:
            end = row.name + pd.to_timedelta(np.floor(daysingamanalysis/2), 'D')
            start = row.name - pd.to_timedelta(np.ceil(daysingamanalysis/2), 'D')
            
        # some tests to make sure everthing is in order
        assert(start > db_min), f'This is currently not implement. WATCH OUT, It might be that you just need to run the langleys!!! The date is too close to the beginning such that the range for the gam anaylsis needs to be pushed closer to the end. Do so.'
        assert(end <= db_max), f'This really should not be possible!!assert({end} <= {db_max})'
        
        ####################
        ## do the processing
        ########
        
        lts = atmlc.read_langleys_from_database(self.path2database_in,start = start, end = end)
        out['lts'] = lts
        pred = lts.makev0prediction4date(row.name)
        
        bla= [{'date': row.name,
            'channel': p['channel'], 
            'V0': p['prediction'].prediction.data, 
            'V0_uncertainty': p['prediction'].prediction_confidence.data,
            'langley_analysis_version': self.langley_analysis_version,
            'input_db': pl.Path(self.path2database_in).name,
            'site': self.site,
            } 
            for p in pred]
        res_formatted = pd.DataFrame(bla)
        
        out_table_name = 'gam_V0_prediction'

        if save:
            with sqlite3.connect(self.database[0]) as con:
                res_formatted.to_sql(out_table_name, con, if_exists='append')
        out['res_formatted'] = res_formatted
        return out