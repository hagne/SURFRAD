from os.path import basename
import os
import datetime
import pandas as pd
import xarray as xr
import add_metadata as am

def normal_process(input, outName, site, yr, mo, dy, lat, lon, elev):
    '''
    Converts all all files in 'input' into a single NetCDF file with name 'outName'
    '''
    try:
    	basestring
    except NameError:
    	basestring = str
    if not isinstance(input, basestring):
        for count,day in enumerate(input):
            print ("%s into %s -- %i/%i" % (get_filename(day), outName, count+1, len(input)))
            df1 = csv_to_dataframe(day)
            df1 = set_headers(df1)
            if count == 0:
                df2 = df1
            else:
                df2 = pd.concat([df2, df1])
            del df1
        write_netcdf(df2, outName)
        am.main(outName, site, yr, mo, dy, lat, lon, elev)
    else:
        df1 = csv_to_dataframe(input)
        df1 = set_headers(df1)
        write_netcdf(df1, outName)
        am.main(outName, site, yr, mo, dy, lat, lon, elev)

def get_testsite(input):   
    '''
    Grabs the initials from the filename and returns the full name of the test site,
    the latitude, longitude, and elevation.
     see: https://www.esrl.noaa.gov/gmd/grad/surfrad/sitepage.html
     Lat and long have been converted to degrees East per CF Conventions.
    '''
    if os.path.isfile(input) == False:
    	exit()
    base = os.path.splitext(basename(input))[0]
    testsite = base[0:3].lower()
    if testsite == "bon":
        return "Bondville_IL","40.05192","271.62691","230"
    elif testsite == "tbl":
        return "Boulder_CO","40.12498","254.7632","1689"
    elif testsite == "dra":
        return "Desert_Rock_NV","36.62373","243.98053","1007"
    elif testsite == "fpk":
        return "Fort_Peck_MT","48.30783","254.8983","634"
    elif testsite == "gwn":
        return "Goodwin_Creek_MS","34.2547","270.1271","98"
    elif testsite == "psu":
        return "Penn_State_PA","40.72012","282.06915","376"
    elif testsite == "sxf":
        return "Sioux_Falls_SD","43.73403","263.37672","473"
    elif testsite == "bao":   # Boulder Atmospheric Observatory
        return "Erie_CO"  ,"","",""
    elif testsite == "slv":
        return "San_Luis_Valley_CO" ,"","" ,""
    elif testsite == "cdn":
        return "Condon_OR"  ,"","",""
    elif testsite == "red":
        return "Red_Lake_AZ","","",""
    elif testsite == "rut":
        return "Rutland_VT","","",""
    elif testsite == "ruf":
        return "Rufus_OR" ,"","",""
    elif testsite == "was":
        return "Wasco_OR","","",""
    else:
        print ("Testing site abbreviation ("+testsite+") does not match any options. Consult nc_convert.py")
        return "no_test_site","?","?","?"

def get_year(input):
	'''
	Returns the year of the file
	For example, bon95002 returns '1995'
	'''
	date = os.path.splitext(basename(input))[0]
	if "_" not in date:    # file format 1  i.e. bon95001.dat
		date = date[3:]    # remove 3 characters from beginning
		date = datetime.datetime.strptime(date, '%y%j')
		return date.strftime('%Y')
	else:				   # file format 2  i.e. was_20160210.qdat
		date = date[:-4]
		date = date[4:]    # remove 4 characters from beginning
		return date

	return date.strftime('%Y')
    
def get_month(input):
    '''
    Returns the month of the file
    For example, bon95002 returns '01'
    '''
    date = os.path.splitext(basename(input))[0]
    if "_" not in input:
    	date = date[3:]  # remove 3 characters from beginning
    	date = datetime.datetime.strptime(date, '%y%j')
    else:
    	date = date[4:]  # remove 4 characters from beginning
    	date = datetime.datetime.strptime(date, '%Y%m%d')

    return date.strftime('%m')
    
def get_day(input):
    '''
    Returns the month of the file
    For example, bon95002 returns '02'
    '''
    date = os.path.splitext(basename(input))[0]
    if "_" not in date:
    	date = date[3:]  # remove 3 characters from beginning
    	date = datetime.datetime.strptime(date, '%y%j')
    else:
    	date = date[4:]  # remove 4 characters from beginning
    	date = datetime.datetime.strptime(date, '%Y%m%d')
    return date.strftime('%d')

def write_netcdf(df1, name):
    '''
    Writes the dataframe into a NetCDF4 file
    'name' should be the full path: "Data/nc/daily/bon/1995/Bondville_IL_1995_01_01.nc"
    '''
    xds = xr.Dataset.from_dataframe(df1)
    
    # The "encoding" parameter allows us to set the datatype of the variables. if not specified, double is usually used.
    # If a variable does not need double precision, then use float32 (f4) to reduce file size.
    xds.to_netcdf(name, encoding={
    				#'date': {'dimension': 0},
    				'Time': {'dtype': 'i4'},
    				'year': {'dtype': 'i4'},
    				'month': {'dtype': 'i4'},
    				'day': {'dtype': 'i4'},
    				'altitude': {'dtype': 'i4'},
    				
    				'latitude': {'dtype': 'float32'},
    				'longitude': {'dtype': 'float32'},
    				
    				'distance_from_sun': {'dtype': 'float32'},
    				'best_estimate_shortwave': {'dtype': 'float32'},
    				'total_shortwave': {'dtype': 'float32'},
    				'diffuse_downwelling_shortwave_flux_in_air': {'dtype': 'float32'},
    				'direct_downwelling_shortwave_flux_in_air': {'dtype': 'float32'},
    				'upwelling_shortwave_flux_in_air': {'dtype': 'float32'},
    				'downwelling_longwave_flux_in_air': {'dtype': 'float32'},
    				'upwelling_longwave_flux_in_air': {'dtype': 'float32'},
    				'air_temperature': {'dtype': 'float32'},
    				'relative_humidity': {'dtype': 'float32'},
    				'surface_air_pressure': {'dtype': 'float32'},
    				'downwelling_longwave_case_temperature': {'dtype': 'float32'},
    				'downwelling_longwave_dome_temperature': {'dtype': 'float32'},
    				'upwelling_longwave_case_temperature': {'dtype': 'float32'},
    				'upwelling_longwave_dome_temperature': {'dtype': 'float32'},
    				'solar_zenith_angle': {'dtype': 'float32'},
    				'downwelling_shortwave_flux_in_air_assuming_clear_sky': {'dtype': 'float32'},
    				'ir_loss_correction_to_diffuse_shortwave': {'dtype': 'float32'},
    				'ir_loss_correction_to_global_shortwave': {'dtype': 'float32'},
    				'downwelling_erythemal_uvb_flux_in_air': {'dtype': 'float32'},
    				'surface_downwelling_photosynthetic_radiative_flux_in_air': {'dtype': 'float32'},
    				'wind_speed': {'dtype': 'float32'},
    				'wind_from_direction': {'dtype': 'float32'},
    				
    				'qc01': {'dtype': 'byte'},
    				'qc02': {'dtype': 'byte'},
    				'qc03': {'dtype': 'byte'},
    				'qc04': {'dtype': 'byte'},
    				'qc05': {'dtype': 'byte'},
    				'qc06': {'dtype': 'byte'},
    				'qc07': {'dtype': 'byte'},
    				'qc08': {'dtype': 'byte'},
    				'qc09': {'dtype': 'byte'},
    				'qc10': {'dtype': 'byte'},
    				'qc11': {'dtype': 'byte'},
    				'qc12': {'dtype': 'byte'},
    				'qc13': {'dtype': 'byte'},
    				'qc14': {'dtype': 'byte'},
    				'qc15': {'dtype': 'byte'},
    				'qc16': {'dtype': 'byte'},
    				'qc17': {'dtype': 'byte'},
    				'qc18': {'dtype': 'byte'},
    				'qc19': {'dtype': 'byte'},
    				
    				'g_flag': {'dtype': 'byte'},
    				'd_flag': {'dtype': 'byte'}
    				})
    	

def csv_to_dataframe(input):
    '''
    Loads a csv into a Pandas dataframe.
    Uses the Year, Julian Day, Hour, and Minute to create an index column of "Date".
    '''
    #parser = lambda time: pd.datetime.strptime(time, '%H%M')
    
    with open(input, 'r') as input_file:
        df1=pd.read_csv(input_file, 
        	encoding = 'utf-8',
            sep=",",
            parse_dates = {'Time': [1]},
            #date_parser = lambda x: pd.to_datetime(x, format="%H%M"),
            index_col = ['Time'])
    df1.drop(list(df1.filter(like="date")), axis=1, inplace=True)

    return df1

def set_headers(df1):
    '''
    Renames all the columns based off of "nc_headers.txt"
    '''
    input="headers_nc.txt"
    with open(input, 'r') as header_file:
        for line in header_file:
            headers = line.split()
    df1.columns = headers
    return df1



#--------------------------------------------------------    
# The following functions are for reference, and NOT actively used 
#-------------------------------------------------------- 
def replace_nan(df1):
    '''
    The value "-9999.9" is specified to be a placeholder for non-existent values. This replaces those values with "NaN"
    '''
    #A[A==-999] = np.nan
    df1.replace(to_replace="-9999.9",value="NaN", inplace=True)
    return df1
    
def add_date(df1):
    '''
    The value "-9999.9" is specified to be a placeholder for non-existent values. This replaces those values with "NaN"
    '''
    #A[A==-999] = np.nan
    df1.createVariable('thedate','c')
    return df1
    
def filter_qc(df1):
    '''
    Drops all the qc columns
    '''
    df1.drop(list(df1.filter(like="qc")), axis=1, inplace=True)
    return df1

def filter_dates(df1):
    '''
    Drops extra date columns
    '''
    df1.drop(df1.columns[[0,1,2]], axis=1, inplace=True)
    return df1

def splitTimeColumn(df1):
    #print (df1)
    # make string version of original column, call it 'col'
    df1['col'] = df1['Time'] #.astype(str)  #.zfill(3)
    #print (df1['col'])
    # make the new columns using string indexing
    df1['Hour'] = df1['col'].str[0:-2]
    df1['Minute'] = df1['col'].str[-2:]
    # get rid of the extra variable (if you want)
    df1.drop('col', axis=1, inplace=True)
    return df1

