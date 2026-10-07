"""
Dependencies
=============
xarray
netcdf4
pandas

"""

import argparse
import inspect
import warnings
import pathlib as pl
import surfradpy.products.radiation2netcdf as srfrad

warnings.simplefilter(action='ignore')

def run(prefix = '/nfs',
        start = None,
        end = None,
        days = 90,
        verbose = True,
        raise_errors = False,
        test = 0,
        ):
    """Convert SURFRAD radiation text products to netCDF files."""
    import productomator.lab as prolab

    reporter = prolab.Reporter(
                'clamps3_radiation2netcdf',
                log_folder='/home/grad/htelg/.processlogs/',
                verbose=True,
                reporting_frequency=(1, 'h'),
            )

    wi = srfrad.Clamps3Radiation2netcdf(
            site='clamps3',
            p2fld_in=f'{prefix}/grad/gradobs/raw/short_term/clamps3/radsys/scaled/csv/',
            p2fld_out=f'{prefix}/grad/gradobs/raw/short_term/clamps3/radsys/netcdf/radiation/v{{version}}/',
            file_name_format='*{date:%Y%m%d}*',
            output_file_format='clamps_rad_{date}.nc',
            start=start,
            end=end,
            days=days,
            input_directory_structure='yearly',
            reporter=reporter,
            verbose=verbose,
    )
    wi.combine_masterplan_duplicates()
    if test == 1:
        print(wi.workplan)
        return wi
    elif test == 2:
        out = wi.process_row(iloc = 0, save = False)
        return out
    elif test == 0:
        wi.process(raise_errors = raise_errors)

    reporter.wrapup(print_dagster_report=True)
    return


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=inspect.getdoc(run) or "",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--prefix', default='/nfs')
    parser.add_argument('--start', default=None)
    parser.add_argument('--end', default=None)
    parser.add_argument('--days', type=int, default=90)
    parser.add_argument('--test', type=int, choices=(0, 1, 2), default=0)
    parser.add_argument('--raise-errors', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose')
    parser.add_argument('--no-verbose', action='store_false', dest='verbose')
    parser.set_defaults(verbose=True)
    args = parser.parse_args(argv)
    return run(
        prefix=args.prefix,
        start=args.start,
        end=args.end,
        days=args.days,
        verbose=args.verbose,
        raise_errors=args.raise_errors,
        test=args.test,
    )


if __name__ == '__main__':
    main()
    
