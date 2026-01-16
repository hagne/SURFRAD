#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2026-01-15

@author: hagen telg

This script converts SURFRAD MFR raw (.xmd) data files to netcdf format. In the process it merges and truncates the data to daily files in UTC time.
"""

import surfradpy.mfr_raw2netcdf as mfr_r2nc

def run():
    prefix = '/Users/htelg/'

    ci = mfr_r2nc.MfrsrRawToNetcdf(f'{prefix}/nfs/grad/Inst/MFR/Campaign/frc/2025/mfrsr/{{serialnumber}}',
                        f'{prefix}/nfs/grad/Inst/MFR/Campaign/frc/2025/mfrsr/{{serialnumber}}.netcdf/v{{version}}/',
                        'frc_{serialnumber}_{year}{month}{day}.nc',
                        serialnumber = '648',
                        verbose = True
                            )