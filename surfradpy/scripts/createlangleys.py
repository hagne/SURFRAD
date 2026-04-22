#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 13:29:47 2023

@author: hagen
"""
import argparse
import inspect
import warnings
warnings.simplefilter(action='ignore')


def run():
    """Create and concatenate SURFRAD MFRSR Langley calibration products."""
    import productomator.lab as prolab
    import surfradpy.realtime_aod as sufrtaod

    reporter = prolab.Reporter(
                'srf_langleys',
                log_folder='/home/grad/htelg/.processlogs/',
                verbose=True,
                reporting_frequency=(3, 'h'),
            )
    proi = sufrtaod.mfrsr_AOD_lev0()
        
    try:      
        # proi.process_all()
        out = proi.process_langleys(raise_error=True, verbose=True)
        reporter.clean_increment(out['numprocessed'])
    except:
        reporter.errors_increment(10)
        reporter.wrapup()
        raise
        
    try:
        proi.process_langley_concat(verbose=True)
    except:
        reporter.errors_increment(out['numprocessed'])
        reporter.wrapup()
        raise
    reporter.wrapup()    


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=inspect.getdoc(run) or "",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.parse_args(argv)
    return run()


if __name__ == '__main__':
    main()
