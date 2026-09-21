import argparse
import inspect
import warnings
import pathlib as pl
warnings.simplefilter(action='ignore')

def run(prefix = '/nfs',
        start = None,
        end = None,
        days = 90,
        verbose = True,
        raise_errors = False,
        ):
    """Convert SURFRAD radiation text products to netCDF files."""
    import productomator.lab as prolab
    import surfradpy.products.radiation2netcdf as srfrad

    reporter = prolab.Reporter(
                'radiation2netcdf',
                log_folder='/home/grad/htelg/.processlogs/',
                verbose=True,
                reporting_frequency=(1, 'h'),
            )
    sites = ["inl",
            "bon",
            "dra",
            "gwn",
            "psu",
            "sxf",
            "tbl",
            "fpk",
    ]
    for site in sites:

        p2fld_in=f'{prefix}/aftp/data/radiation/surfrad/{site}'
        if not pl.Path(p2fld_in).exists():
            p2fld_in=f'{prefix}/iftp/aftp/data/radiation/surfrad/{site}'
            if not pl.Path(p2fld_in).exists():
                raise FileNotFoundError(f'Folder does not exist: {p2fld_in}.')

        wi = srfrad.SurfradRadiation2netcdf(
            site=site,
            p2fld_in=p2fld_in,
            p2fld_out=f'{prefix}/grad/surfrad/products_level1/radiation_netcdf/v{{version}}/{{site}}',
            file_name_format='*{date:%y%j}*',
            output_file_format='srf_rad_full_{site}_{date}.nc',
            start=start,
            end=end,
            days=days,
            input_directory_structure='yearly',
            reporter=reporter,
            verbose=verbose,
        )
        wi.process(raise_errors = raise_errors)

    reporter.wrapup()
    return

def run_deprecated():
    """Convert SURFRAD radiation text products to netCDF files."""
    import productomator.lab as prolab
    import surfradpy.products.radiation2netcdf as srfrad

    reporter = prolab.Reporter(
                'radiation2netcdf',
                log_folder='/home/grad/htelg/.processlogs/',
                verbose=True,
                reporting_frequency=(1, 'h'),
            )
    try:
        out = srfrad.generate_netcdfs(p2fld = '/nfs/iftp/aftp/data/radiation/surfrad/',
                                p2fldout = '/nfs/grad/surfrad/products_level1/radiation_netcdf/',
                                gui=False,
                                 verbose = False)
        reporter.clean_increment(out['numprocessed'])
        reporter.wrapup()
    except:
        reporter.errors_increment(7)
        reporter.wrapup()
        raise

    return


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=inspect.getdoc(run) or "",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--start', default=None)
    parser.add_argument('--end', default=None)
    parser.add_argument('--days', type=int, default=360)
    parser.add_argument('--raise-errors', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose')
    parser.add_argument('--no-verbose', action='store_false', dest='verbose')
    parser.set_defaults(verbose=True)
    args = parser.parse_args(argv)
    return run(
        start=args.start,
        end=args.end,
        days=args.days,
        verbose=args.verbose,
        raise_errors=args.raise_errors,
    )


if __name__ == '__main__':
    main()
    
