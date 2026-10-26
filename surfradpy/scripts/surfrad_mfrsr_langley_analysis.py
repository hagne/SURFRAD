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

#todo: pupulate main with argument parsing and call run with those arguments. Add some documentation for meaningfull --help.
def main():
    return 

if __name__ == "__main__":
    main()
