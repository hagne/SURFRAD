#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2026-02-02

@author: Hagen Telg

This scirpt runs the complete processing chain of MFRSRs from raw files all the way to the highest level process. 
As most processes depends on the successfull completion of a previous step it does so successively. Each process can be run
independently.
"""

import surfradpy.scripts as surfscritp

def run(prefix = '/nfs/grad/'):
    surfscritp.surfrad_mfrsr_raw2netcdf.run(prefix = prefix)
    # TODO: add cosine calibration script
    # TODO: add add langley calibration script, not sure if this is the way to go ... do same as PMOD
    # TODO: add aod script
    # TODO: add sizedistribution script…