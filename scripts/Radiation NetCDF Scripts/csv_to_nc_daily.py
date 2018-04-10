"""
This module accepts .csv files and returns .nc files.
"""

import sys
from os.path import basename
import os
import nc_convert as ncc

def main(filesToProcess):
    numOfFiles = len(filesToProcess)
    
    if numOfFiles>1:
        for count,input in enumerate(filesToProcess):
            outName = set_outName(input)
            site,lat,lon,elev = ncc.get_testsite(input)
            ncc.normal_process(input, outName, site, ncc.get_year(input), ncc.get_month(input), ncc.get_day(input), lat, lon, elev)
    else:
        input = filesToProcess[0]
        outName = set_outName(input)
        site,lat,lon,elev = ncc.get_testsite(input)
        ncc.normal_process(input, outName, site, ncc.get_year(input), ncc.get_month(input), ncc.get_day(input), lat, lon, elev)

def set_outName(input):
    '''
    Sets the output name for daily files
    '''
    site,lat,lon,elev = ncc.get_testsite(input)
    year = ncc.get_year(input)
    month = ncc.get_month(input)
    day = ncc.get_day(input)
    #print("%s_%s_%s_%s.nc" % (site, year, month, day))
    #return "Data/nc/daily/%s/%s/%s_%s_%s_%s.nc" % (site, year, site, year, month, day)
    folderpath = os.path.dirname(os.path.abspath(input))
    print(folderpath+"/%s_%s_%s_%s_uncompressed.nc" % (site, year, month, day))
    return folderpath+"/%s_%s_%s_%s_uncompressed.nc" % (site, year, month, day)

if __name__ == '__main__':
    main(sys.argv[1:])
