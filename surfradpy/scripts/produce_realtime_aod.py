#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 13:29:47 2023

@author: hagen
"""
import productomator.lab as prolab
import surfradpy.realtime_aod as sufrtaod
import warnings
warnings.simplefilter(action='ignore')


def routine():
    autothis = sufrtaod.mfrsr_AOD_lev0()
    autothis.process_all()

    #### TODO: update to reporter
    auto = prolab.Automation(autothis)
    auto.log()

def catchup():
    autothis = sufrtaod.mfrsr_AOD_lev0(site=['psu',])
    autothis.process_all()

    #### TODO: update to reporter
    auto = prolab.Automation(autothis, product_name='produce_aod_realtime_catchup')
    auto.log()