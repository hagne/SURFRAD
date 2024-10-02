import productomator.lab as prolab
import surfradpy.radiation as srfrad
import warnings
warnings.simplefilter(action='ignore')

def main():
    reporter = prolab.Reporter(
                'radiation2netcdf',
                log_folder='/home/grad/htelg/.processlogs/',
                verbose=True,
                reporting_frequency=(1, 'h'),
            )
    try:
        srfrad.generate_netcdfs(p2fld = '/nfs/iftp/aftp/data/radiation/surfrad/',
                                p2fldout = '/nfs/grad/surfrad/products_level1/radiation_netcdf/',
                                gui=False,
                                 verbose = False)
        reporter.clean_increment()
        reporter.wrapup()
    except:
        reporter.errors_increment()
        reporter.wrapup()
        raise

    return
    