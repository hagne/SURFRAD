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
import surfradpy.scripts as surfscritp

def run(prefix = '/nfs/grad/', verbose = False):
    if verbose:
        print(f'Running MFRSR full process with prefix={prefix}')
    surfscritp.surfrad_mfrsr_raw2netcdf.run(prefix = prefix, verbose = verbose)
    # TODO: add cosine calibration script
    # TODO: add add langley calibration script, not sure if this is the way to go ... do same as PMOD
    # TODO: add aod script
    # TODO: add sizedistribution script…

def main():
    parser = argparse.ArgumentParser(
        description='Run the complete MFRSR processing chain.'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output.',
    )
    args = parser.parse_args()
    run(verbose=args.verbose)

if __name__ == '__main__':
    main()
