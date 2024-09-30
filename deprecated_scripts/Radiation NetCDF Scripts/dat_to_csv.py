import sys
from os.path import basename
import os
import datetime
import numpy as np
import nc_convert as ncc

'''
EXTRACTING DATA FROM THE FILENAMES ------
Some testsites are prepended in the filename, others aren't.
Use this to add 3 letter code to filename:
for f in *.lw1; do mv "$f" "bon_$f"; done 

We've got 2 different filename formats in the data...
1) bon95001.dat
2) was_20160210.qdat
The functions get_day, get_month, and get_year
decide which, then parse it correctly.
'''

def main(filesToProcess):
	if len(filesToProcess)>1: # if list of 2 or more files
		for input in filesToProcess:
			dat_to_csv(input)
	else: # if just single file
		dat_to_csv(filesToProcess[0])

def dat_to_csv(input):

	headers = get_headers()
	filename = os.path.splitext(basename(input))[0]  # Returns the name of the file without the extension.

	month = get_month(filename)
	year = get_year(filename)
	day = get_day(filename)
	
	site,lat,lon,alt = ncc.get_testsite(input)

	# TO DEBUG
	#print("year: "+year)
	#print("month: "+month)
	#print("day: "+day)
	#print("input: "+input)

	site = get_site(input)
	#out_name = "Data/csv/%s/%s/%s/%s.csv" % (site, year, month, base)
	out_name = "%s.csv" % (input)
	
	# Clear all known extensions
	out_name = out_name.replace(".dat","").replace(".DAT","").replace(".qDAT","").replace(".qdat","").replace(".qadj","").replace(".lw1","")

	with open(input, 'r') as input_file:
		with open(out_name, 'w') as output_file:
			output_file.write(headers + "\n")
			for count, line in enumerate(input_file):
				if count < 2:
					pass
				else:
					if line.split()[0]=='\x1a':
						pass
					else:
						outLine = ",".join(line.split())

					    # Fix G-RAD time into proper minutes
					    # i.e. Time 105 = 65 seconds
						weirdseconds = outLine.split(',')[1].strip()
						hr = weirdseconds[0:-2]
						min = weirdseconds[-2:]
						if hr == '':
							hr = 0
						# add in corrected time column, and a column for year, month, and day and lat, lon, and altitude
						addToCSV = "," + str((int(hr) * 60) + int(min)) + "," + str(year) + "," + str(month) + "," + str(day) + "," + str(lat) + "," + str(lon) + "," + str(alt) + ","
						a = "," + str(weirdseconds) + ","
						thisrow = str.replace(outLine, a, addToCSV, 1)
						
						# Write output to file
						output_file.write(thisrow + "\n")

def get_headers(input="headers_csv.txt", get_csv=True):
    headers = []
    with open(input, 'r') as header_file:
        for line in header_file:
            headers = line.split()
    if get_csv:
        return ('%s') % ','.join(headers)
    else: # leave the option of returning an array of headers instead
        return headers

def get_year(input):
	'''
	Returns the year of the file
	For example, bon95002 returns '1995'
	'''
	if "_" not in input:    # file format 1  i.e. bon95001.dat
		date = input[3:]    # remove 3 characters from beginning
		date = datetime.datetime.strptime(date, '%y%j')
		return date.strftime('%Y')
	else:					# file format 2  i.e. was_20160210.qdat
		date = input[:-4]
		date = date[4:]     # remove 4 characters from beginning
		return date

	return date.strftime('%Y')
    
def get_month(input):
    '''
    Returns the month of the file
    For example, bon95002 returns '01'
    '''
    if "_" not in input:
    	date = input[3:]  # remove 3 characters from beginning
    	date = datetime.datetime.strptime(date, '%y%j')
    else:
    	date = input[4:]  # remove 4 characters from beginning
    	date = datetime.datetime.strptime(date, '%Y%m%d')

    return date.strftime('%m')
    
def get_day(input):
    '''
    Returns the day of the file
    For example, bon95002 returns '02'
    '''
    if "_" not in input:
    	date = input[3:]  # remove 3 characters from beginning
    	date = datetime.datetime.strptime(date, '%y%j')
    else:
    	date = input[4:]  # remove 4 characters from beginning
    	date = datetime.datetime.strptime(date, '%Y%m%d')
    return date.strftime('%d')
    
def get_site(input):
	return input[0:3]

if __name__ == '__main__':
    main(sys.argv[1:])
