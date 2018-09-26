"""
This module should accept .csv files and return .nc files.
"""

import sys
from os.path import basename
import os
import nc_convert as ncc

def main(filesToProcess):
    numOfFiles = len(filesToProcess)
    print ("number of files: "+str(numOfFiles))
    if numOfFiles>1:
        input = filesToProcess[0]
        outName = set_outName(input)
        ncc.normal_process(filesToProcess, outName, ncc.get_year(input), "", "")
    else:
        input = filesToProcess[0]
        outName = set_outName(input)
        #print("got here2")
        ncc.normal_process(input, outName, ncc.get_year(input), "", "")


def set_outName(input):
    '''
    Sets the output name for yearly files
    '''
    site = ncc.get_testsite(input)
    year = ncc.get_year(input)
    folderpath = os.path.dirname(os.path.abspath(input))
    #return folderpath+"/%s_%s_%s_%s_uncompressed.nc" % (site, year, month, day)
    return folderpath+"/%s_%s_uncompressed.nc" % (site, year)

if __name__ == '__main__':
    main(sys.argv[1:])
