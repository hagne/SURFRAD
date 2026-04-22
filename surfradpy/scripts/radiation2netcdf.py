import argparse
import inspect
import warnings
warnings.simplefilter(action='ignore')

def run():
    """Convert SURFRAD radiation text products to netCDF files."""
    import productomator.lab as prolab
    import surfradpy.radiation as srfrad

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
    parser.parse_args(argv)
    return run()


if __name__ == '__main__':
    main()
    
