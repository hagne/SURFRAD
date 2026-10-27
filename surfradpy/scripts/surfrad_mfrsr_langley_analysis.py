"""Run morning and afternoon Langley analyses for SURFRAD MFRSR data.

The :func:`run` workflow reads cosine-corrected MFRSR NetCDF files, cleans and
fits the Langley data for each requested site, stores the results in per-site
SQLite databases, and records processing status through ``productomator``.

Minimum Conda runtime packages (from ``conda-forge``)::

    python>=3.11
    numpy
    pandas
    xarray
    netcdf4
    scipy
    statsmodels
    pvlib

e.g.:
conda create -n [env_name] -c conda-forge \
    'python>=3.11' numpy pandas xarray scipy netcdf4 statsmodels pvlib

The local ``atmPy``, ``productomator``, and ``surfradpy`` projects must also be
installed. Add ``pip`` to the Conda environment only when using it to install
those projects; it is not a runtime dependency.
"""

import argparse

import productomator.lab as prolab
import surfradpy.products.langley_analysis  as srfla


def run(prefix = '/nfs',
        log_folder='/home/grad/htelg/.processlogs/',
        start = None,
        end = None,
        days = 90,
        site = None,
        test = 0,
        raise_errors = False,
        verbose = True,):

    out = {}
    logname = 'surfrad_mfrsr_langley_analysis'
    reporter = prolab.Reporter(logname, 
                                log_folder=log_folder,
                                reporting_frequency=(6, 'h'),
                            )

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

    version_in = '0.2'

    for site in sites:
        if verbose:
            print(site)
            print('-----')

        p2fld_in = f'{prefix}/grad/Inst/MFR/SURFRAD/{site}/mfrsr/cosine_corrected/v{version_in}/'
        p2db = f'{prefix}/grad/surfrad/products_level2/aod/langley_analysis_{site}_v{{version}}.db'
        li = srfla.LangleyAnalysis(
                 p2fld_in,
                 p2database=p2db,
                 file_name_format = '*{date:%Y%m%d}*',
                 start = start,
                 end = end,
                 days = days,
                 input_directory_structure='yearly',
                 reporter = reporter,
                 verbose = verbose,
                 langley_airmass_limits = (2, 6),
                 expected_instrument_noise = 0.01,
                 langley_clean_threshold = 2,
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
        description=(
            "Run morning and afternoon Langley analyses for SURFRAD MFRSR "
            "cosine-corrected data."
        ),
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
