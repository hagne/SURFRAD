"""
Docstring for surfradpy.database
created 2026-01-15

@author: hagen telg

This module provides database access functions for all SURFRAD data and instruments.
"""

import datetime as dt
import functools
import getpass
import inspect
import json
import sqlite3
import pathlib as pl
import pandas as pd
import surfradpy.config as sfp_config


def get_default_db_path() -> pl.Path | None:
    """
    Return the configured default database path, if available.
    """
    return sfp_config.get_db_path()


def _quote_identifier(name):
    if not name.isidentifier():
        raise ValueError(f"Invalid database identifier: {name}")
    return f'"{name}"'


def _filter_audit_arguments(func, args, kwargs, exclude_args):
    signature = inspect.signature(func)
    parameters = list(signature.parameters.values())
    if parameters and parameters[0].name == "self":
        parameters = parameters[1:]

    logged_args = []
    arg_index = 0
    for parameter in parameters:
        if parameter.kind == inspect.Parameter.VAR_POSITIONAL:
            if parameter.name not in exclude_args:
                logged_args.extend(args[arg_index:])
            arg_index = len(args)
            break
        if parameter.kind not in (
            inspect.Parameter.POSITIONAL_ONLY,
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
        ):
            continue
        if arg_index >= len(args):
            break
        if parameter.name not in exclude_args:
            logged_args.append(args[arg_index])
        arg_index += 1

    if arg_index < len(args):
        logged_args.extend(args[arg_index:])

    logged_kwargs = {
        key: value
        for key, value in kwargs.items()
        if key not in exclude_args
    }

    return logged_args, logged_kwargs


def audit_db_call(exclude_args=None, table_name="audit_log"):
    """
    Log decorated database method calls to an audit table.
    """
    exclude_args = set(exclude_args or [])
    quoted_table_name = _quote_identifier(table_name)

    def decorator(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            logged_args, logged_kwargs = _filter_audit_arguments(func, args, kwargs, exclude_args)

            result = func(self, *args, **kwargs)

            event = {
                "time_of_execution": dt.datetime.now(dt.UTC).isoformat(),
                "function": func.__qualname__,
                "args": json.dumps(
                    {"args": logged_args, "kwargs": logged_kwargs},
                    default=str,
                ),
                "user": getpass.getuser(),
            }

            self.ensure_audit_log_table(table_name=table_name)
            self.execute_query(
                f"""
                INSERT INTO {quoted_table_name} (
                    time_of_execution, function, args, user
                )
                VALUES (
                    :time_of_execution, :function, :args, :user
                )
                """,
                event,
            )

            return result

        return wrapper

    return decorator


class SurfradDatabase:
    """
    Class to handle SURFRAD database connections and queries.
    """
    def __init__(self, path2db: str | pl.Path | None, create_if_missing: bool = False):
        """
        Initialize the SurfradDatabase with the path to the database file.

        Parameters:
        -----------
        path2db: str
            Path to the SQLite database file.
        """
        if path2db is None:
            path2db = get_default_db_path()
        if path2db is None:
            cfg_path = sfp_config.get_config_path()
            raise FileNotFoundError(
                "No database path provided and no default configured. "
                f"Set SURFRAD_DB_PATH or add [database] path to {cfg_path}."
            )
        if not create_if_missing:
            if not pl.Path(path2db).exists():
                raise FileNotFoundError(
                    f"Database file {path2db} does not exist. "
                    "Set create_if_missing=True to create a new database."
                )
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

    
    def find_site_info(self, abb = None, site_id = None, name = None, allow_multiple = False):
        """
        Find site information for a given site code from the database.

        Parameters
        ----------
        site_code : str
            The site code to look up.

        Returns
        -------
        dict
            A dictionary containing site information such as name, latitude, and longitude.
        """

        # only one of the parameters should be provided
        if sum([abb is not None, site_id is not None, name is not None]) != 1:
            raise ValueError("Exactly one of abb, site_id, or name must be provided.")
        # which one is provided?
        if abb is not None:
            column = 'abb'
            value = abb
        elif site_id is not None:
            column = 'site_id'
            value = site_id
        elif name is not None:
            column = 'name'
            value = name
        else:
            raise ValueError("This should never happen, but exactly one of abb, site_id, or name must be provided.")
        
        query = f"SELECT * FROM sites WHERE {column} = ?"
        df = self.execute_query(query, (value,))
        if df.empty:
            raise ValueError(f"No site found with {column} {value}")
        elif len(df) > 1:
            if not allow_multiple:
                raise ValueError(f"Multiple sites found with {column} {value}, expected only one. Set allow_multiple=True to return the full dataframe.")
            else:
                return df
        row = df.iloc[0]
        return row
    
    def execute_query(self, query: str, *args):
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
                cur.execute(query, *args)
                if cur.description is None:
                    db.commit()
                    return None
                rows = cur.fetchall()
                col_names = [desc[0] for desc in cur.description]
            finally:
                cur.close()
        return pd.DataFrame(rows, columns=col_names)
    
    def ensure_audit_log_table(self, table_name="audit_log"):
        """
        Create the audit log table if it does not exist.
        """
        table_name = _quote_identifier(table_name)
        self.execute_query(
            f"""
            CREATE TABLE IF NOT EXISTS {table_name} (
                time_of_execution TEXT,
                function TEXT,
                args TEXT,
                user TEXT
            )
            """
        )

    @audit_db_call()
    def add_site(self, abb, name, elevation, latitude, longitude, network='', country=''):
        """
        Add a site entry to the sites table.
        Parameters
        ----------
        abb : str
            Site abbreviation (e.g., 'tbl').
        name : str
            Site name.
        elevation : float
            Site elevation in meter.
        latitude : float
            Site latitude in degrees north.
        longitude : float
            Site longitude in degrees east.
        network : str, optional
            Site network or Campaign.
        country : str, optional
            Site country.
        """

        idx = self.execute_query('SELECT COALESCE(MAX("index"), 0) + 1 FROM sites')
        idx = int(idx.iloc[0,0])

        site_id = self.execute_query('SELECT COALESCE(MAX("site_id"), 0) + 1 FROM sites')
        site_id = int(site_id.iloc[0,0])

        site = {
            "index": idx,
            "site_id": site_id,
            "abb": abb.lower(),
            "name": name,
            "elevation": elevation,
            "latitude": latitude,
            "longitude": longitude,
            "network": network,
            "country": country,
        }

        self.execute_query(
            """
            INSERT INTO sites ("index", site_id, abb, name, elevation, latitude, longitude, network, country)
            VALUES (:index, :site_id, :abb, :name, :elevation, :latitude, :longitude, :network, :country)
            """,
            site,
        )

    @audit_db_call()
    def remove_site(self, abb = None, site_id = None, name = None, allow_multiple = False):
        """
        Remove a site entry from the sites table.
        """
        site = self.find_site_info(
            abb=abb,
            site_id=site_id,
            name=name,
            allow_multiple=allow_multiple,
        )
        sites = site.to_dict("records") if isinstance(site, pd.DataFrame) else [site.to_dict()]

        for site in sites:
            columns = site.keys()

            where_clause = " AND ".join(
                f'"{column}" = :{column}'
                for column in columns
            )

            sql = f"""DELETE FROM sites WHERE {where_clause}"""

            self.execute_query(sql, site)

    @audit_db_call()
    def add_table(self, tbl_name):
        """
        Create a new table in the database.
        """
        tbl_name = _quote_identifier(tbl_name)
        self.execute_query(f"CREATE TABLE {tbl_name} (id INTEGER PRIMARY KEY)")

    @audit_db_call()
    def remove_table(self, tbl_name):
        """
        Delete a table from the database.
        """
        tbl_name = _quote_identifier(tbl_name)
        self.execute_query(f"DROP TABLE {tbl_name}")

    @audit_db_call()
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
    
