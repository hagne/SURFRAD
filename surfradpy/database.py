"""
Docstring for surfradpy.database
created 2026-01-15

@author: hagen telg

This module provides database access functions for all SURFRAD data and instruments.
"""

import sqlite3
import pathlib as pl
import pandas as pd


class SurfradDatabase:
    """
    Class to handle SURFRAD database connections and queries.
    """
    def __init__(self, path2db: str, create_if_missing: bool = False):
        """
        Initialize the SurfradDatabase with the path to the database file.

        Parameters:
        -----------
        path2db: str
            Path to the SQLite database file.
        """
        if not create_if_missing:
            if not pl.Path(path2db).exists():
                raise FileNotFoundError(f"Database file {path2db} does not exist. Set create_if_missing=True to create a new database.")
        self.path2db = pl.Path(path2db)

    def snapshot(self, max_rows=20, include_schema=True):
        """Return a lightweight representation of the database."""

        with sqlite3.connect(self.path2db) as db:
            tables = [row[0] for row in db.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'").fetchall()]
            schema = {}
            preview = {}
            for tbl in tables:
                if include_schema:
                    schema_row = db.execute(
                        "SELECT sql FROM sqlite_master WHERE type='table' AND name=?",
                        (tbl,),
                    ).fetchone()
                    schema[tbl] = schema_row[0] if schema_row else ''
                preview[tbl] = pd.read_sql(
                    f"SELECT * FROM {tbl} LIMIT {int(max_rows)}", db)
        return {'tables': tables, 'schema': schema, 'preview': preview}

    def dump_table(self, tbl_name, index_col = None):
        qu = 'select * from {}'.format(tbl_name)
        with sqlite3.connect(self.path2db) as db:
            df_db = pd.read_sql(qu, db, index_col= index_col)
        return df_db
    
    def update_mfrsr_history(self, path2mfrsr_history = 'MFRSR_History.xlsx', replace: bool = False, dryrun: bool = False):
        """
        Update the MFRSR history information in the database from an Excel file that can be found on grad's google
        drive (search for MFRSR_History.xlsx).

        Parameters
        ----------
        path2mfrsr_history : str, optional
            Path to the MFRSR History Excel file.
        replace : bool, optional
            If True, replace existing data in the database. TODO: implement incremental updates.
        dryrun : bool, optional
            If True, do not write to the database, just return the parsed data.
        """
        def get_facing(row):
            if not isinstance(row['Table/Tower'], str):
                return 'unk'       
            tts = row['Table/Tower'].split('/')
            if len(tts) !=2:
                return 'unk'
            else:
                if k == tts[0]:
                    facing = 'up'
                elif k == tts[1]:
                    facing = 'down'
                else:
                    raise ValueError(f'Instrument {k} as an unconclusive "Table/Tower" value of {row['Table/Tower']} at Date_start: {row.Date_start}')
            return facing
        
        if_exists = 'replace' if replace else 'fail'
        out = {}

        # read the history file
        sheets = pd.read_excel(path2mfrsr_history, sheet_name=None)
        sheets.keys()

        df = sheets['Overview']

        # get rows where instrument is an integer, meaning it is actually an instrument with a valid serial number, which are all integers
        df = df[df.apply(lambda row:isinstance(row.Instrument, int), axis = 1)]
        instruments_mfrsr = df
        out['instruments_mfrsr'] = instruments_mfrsr

        # get all the tables
        instruments_keys = [k for k in list(sheets.keys()) if k.strip().isnumeric()]


        inst_list = []
        for k in instruments_keys:
            df = sheets[k]
            if df.columns[0] != 'Date_start':
                print(f'Instrument {k} has fist column with the name {df.columns[0]} not "Date_start" as expected. This indicates that the instrument was never put into action ... bought for spare parts?')
                continue
            # self.tp_df = df
            # print(type(df.Date_start))

            # assert(df.Date_start.__name__ == "Timestamp")
            assert(df.Date_start.dtypes.name == 'datetime64[us]'), f'Date_start column of the {k} sheet does not have a valid datetime format (has {df.Date_start.dtype.name}). Someone must have entered a non valid entry in that column of MFRSR_history.'
            
            # only get relevant columns
            df = df[['Date_start','Date_stop', 'Location', 'Table/Tower']]
            
            # remove any rows beyond the last time samp in Date_start
            df = df[~df.Date_start.isna()]
            
            # add a column with the instrument serial number
            df['instrument_type_id'] = 1
            df['instrument_type'] = 'mfrsr'
            df['instrument_sn'] = k
            
            # translate the Table/Tower meaning
            df['facing'] = df.apply(lambda row: get_facing(row), axis=1)
            
            df = df.drop('Table/Tower', axis=1)
            inst_list.append(df)

        deployments = pd.concat(inst_list)
        deployments['Date_stop'] = deployments.Date_stop.astype(str).str.lower() # this is a mixed format column ... also make sure 'present' is lowercase
        deployments['Location'] = deployments['Location'].str.lower() # make sure location is lowercase
        # update the table
        if not dryrun:
            with sqlite3.connect(self.path2db) as db:
                instruments_mfrsr.to_sql('instruments_mfrsr', db, if_exists='replace')
                deployments.to_sql('deployments', db, if_exists= if_exists)

        out['deployments'] = deployments
        return out
    
    def execute_query(self, query: str):
        """
        Execute a custom SQL query on the database and return the result as a pandas DataFrame.

        Parameters
        ----------
        query : str
            The SQL query to execute.

        Returns
        -------
        pd.DataFrame
            The result of the query.

        Examples
        --------
        >>> db = SurfradDatabase('surfrad.db')
        >>> df = db.execute_query('SELECT * FROM instruments_mfrsr WHERE instrument_sn = 648')
        """
        with sqlite3.connect(self.path2db) as db:
            cur = db.cursor()
            try:
                cur.execute(query)
                if cur.description is None:
                    db.commit()
                    return None
                rows = cur.fetchall()
                col_names = [desc[0] for desc in cur.description]
            finally:
                cur.close()
        return pd.DataFrame(rows, columns=col_names)
