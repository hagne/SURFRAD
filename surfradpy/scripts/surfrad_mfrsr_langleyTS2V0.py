"""Generate MFRSR extraterrestrial calibration (V0) predictions.

For one or all SURFRAD sites, this script reads qualified AM/PM Langley fits
from a site SQLite database, fits a robust trend-and-seasonality GAM with
R/mgcv for each wavelength, and stores the predicted V0, its uncertainty, and
the distances to the nearest valid Langley fits in a per-site SQLite database.

The minimum Conda environment is::

    conda create -n [env_name] -c conda-forge \
        'python>=3.11' numpy pandas xarray scipy r-base r-mgcv

Install the local ``atmPy``, ``productomator``, and ``surfradpy`` packages into
that environment separately; add ``pip`` only if using it for those installs.
``argparse``, ``pathlib``, and ``sqlite3`` come with Python; NetCDF4, dask,
matplotlib, statsmodels, and pyGAM are not required for this V0-prediction
workflow. The configured ``atmPy`` Rscript path must point to the ``Rscript``
executable in the environment.
"""

import argparse

import productomator.lab as prolab
import surfradpy.products.langley_analysis  as srfla
import pandas as pd


def run(prefix = '/nfs',
        log_folder='/home/grad/htelg/.processlogs/',
        start = None,
        end = None,
        days = 90,
        realtime = False,
        site = None,
        test = 0,
        raise_errors = False,
        verbose = True,):

    out = {}
    if realtime:
        realtimeindicator = 'nearrealtime'
    else:
        realtimeindicator = 'final'
    logname = f'surfrad_mfrsr_langleyTS2V0_{realtimeindicator}'
    reporter = prolab.Reporter(logname, 
                                log_folder=log_folder,
                                reporting_frequency=(6, 'h'),
                            )
    
    if (start is None) and (end is None):
        if realtime:
            end = pd.to_datetime(pd.Timestamp().now().date())
        else:
            end = pd.to_datetime(pd.Timestamp().now().date()) - pd.to_timedelta(60, 'D')
        start = end-pd.to_timedelta(days, 'D')

            


    if site is not None:
        sites = [site]
    else:
        sites = ['tbl',
                'bnd',
                'dra',
                'gwn',
                'psu', 
                'sxf',
                # 'inl',
                'fpe',
        ]

    langley_analysis_version = '0.1'

    for site in sites:
        if verbose:
            print(site)
            print('-----')

        p2db_in = f'{prefix}/grad/surfrad/products_level2/aod/langley_analysis_{site}_v{langley_analysis_version}.db'
        p2db_out = f'{prefix}/grad/surfrad/products_level2/aod/v0_prediction_{site}_v{{version}}_{realtimeindicator}.db'
        li = srfla.V0prediction(
                 p2db_in,
                 database=[p2db_out, 'gam_V0_prediction', 'date', 'None'],
                 start = start,
                 end = end,
                 days = days,
                 realtime=realtime,
                 reporter = reporter,
                 verbose = verbose,
        )
        print(f'{site} workplan.shape: {li.workplan.shape}')
        if test == 1:
            last_processed = None
            break
        elif test == 2:
            last_processed = li.process_row(iloc = 1, save=False)
            break
        elif test == 3:
            last_processed = li.process_row(iloc = 1, save=True)
            break
        else:
            # try:

            last_processed = li.process(raise_errors = raise_errors)
    out['product_instance'] = li
    out['last_processed'] = last_processed
    reporter.wrapup(print_dagster_report=True)
    out['reporter'] = reporter
    if test > 0:
        return out
    return 

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate MFRSR V0 predictions from Langley-analysis results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--prefix",
        default="/nfs",
        help="Base filesystem prefix containing SURFRAD data products.",
    )
    parser.add_argument(
        "--log-folder",
        default="/home/grad/htelg/.processlogs/",
        help="Directory in which processing logs are stored.",
    )
    parser.add_argument(
        "--start",
        default=None,
        help="Start date or datetime for the processing window.",
    )
    parser.add_argument(
        "--end",
        default=None,
        help="End date or datetime for the processing window.",
    )
    parser.add_argument(
        "--days",
        type=int,
        default=90,
        help="Number of days to include in the processing window.",
    )
    parser.add_argument(
        "--site",
        default=None,
        help="SURFRAD site code to process; omit to process every site.",
    )
    parser.add_argument(
        "--test",
        type=int,
        choices=(0, 1, 2, 3),
        default=0,
        help=(
            "Test mode: 0=full run, 1=build workplan only, "
            "2=process one row without saving, 3=process one row and save."
        ),
    )
    parser.add_argument(
        "--raise-errors",
        action="store_true",
        help="Propagate processing errors instead of handling them internally.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        help="Print processing progress.",
    )
    parser.add_argument(
        "--no-verbose",
        action="store_false",
        dest="verbose",
        help="Suppress processing progress.",
    )
    parser.set_defaults(verbose=True)
    args = parser.parse_args(argv)
    return run(
        prefix=args.prefix,
        log_folder=args.log_folder,
        start=args.start,
        end=args.end,
        days=args.days,
        site=args.site,
        test=args.test,
        raise_errors=args.raise_errors,
        verbose=args.verbose,
    )

if __name__ == "__main__":
    main()
