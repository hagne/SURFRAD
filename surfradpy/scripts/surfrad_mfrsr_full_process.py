#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2026-02-02

@author: Hagen Telg

This scirpt runs the complete processing chain of MFRSRs from raw files all the way to the highest level process. 
As most processes depends on the successfull completion of a previous step it does so successively. Each process can be run
independently.
"""

import argparse
import surfradpy.scripts.surfrad_mfrsr_raw2netcdf as raw2nc
import surfradpy.scripts.surfrad_mfrsr_cosinecalibration as cosine_cal

def run(prefix = '/nfs', verbose = False):
    """Run the full SURFRAD MFRSR processing chain. Including raw2netcdf, cosine calibration

    Changes
    -------
    2026-04-22
        - cosine correction added
        - timestamp adjustemnt takes affect

    Parameters
    ----------
    prefix : str, optional
        Base filesystem prefix used by downstream processing scripts, by default '/nfs'
    verbose : bool, optional
        Enable verbose output, by default False"""
    if verbose:
        print(f'Running MFRSR full process with prefix={prefix}')
    raw2nc.run(prefix=prefix, verbose=verbose)
    cosine_cal.run(prefix=prefix, verbose=verbose)
    # TODO: add cosine calibration script
    # TODO: add add langley calibration script, not sure if this is the way to go ... do same as PMOD
    # TODO: add aod script
    # TODO: add sizedistribution script…

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=run.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        '--prefix',
        default='/nfs',
        help='Base filesystem prefix used by downstream processing scripts.',
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output.',
    )
    args = parser.parse_args(argv)
    run(prefix=args.prefix, verbose=args.verbose)

if __name__ == '__main__':
    main()
