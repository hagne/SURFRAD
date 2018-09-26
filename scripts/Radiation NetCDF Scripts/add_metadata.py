'''
Input: NetCDF file
Output: NetCDF file with Surfrad metadata
'''

import sys
import datetime
from netCDF4 import Dataset

def main(input, site, yr, mo, dy, lat, lon, elev):
    #print ("Adding metadata to %s" % (input))  # For debugging
    f = Dataset(input, "r+", format="NETCDF4")
    f = add_globals(f, input, site, yr, mo, dy, lat, lon, elev)
    f = add_var_attrs(f, yr, mo, dy)
    f.close()

def add_globals(f, input, site, yr, mo, dy, lat, lon, elev):
    f.Conventions='CF-1.6'
    f.title = 'QCRad V3 Dataset from NOAA/ESRL/GMD/GRAD Surfrad Radiation Archive'
    f.date = "Data recorded on "+str(mo).zfill(2)+'/'+str(dy)+'/'+str(yr)
    f.location = str(site)
    f.latitude = str(lat)+" N"
    f.longitude = str(lon)+" E"
    f.elevation = str(elev)+" meters"
    f.source = 'NOAA/ESRL/GMD/GRAD Surfrad Radiation Archive'
    f.references = 'Online: https://www.esrl.noaa.gov/gmd/grad/surfrad/index.html'
    f.history = '['+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")+'] Created compressed netCDF4 dataset'
    f.institution = 'National Oceanic and Atmospheric Administration (NOAA) - David Skaggs Research Center - Boulder, CO'
    f.comment = 'Surface radiation data saved in daily files by test site'
    return f

def add_var_attrs(f, yr, mo, dy):
	#f.createVariable('date','S4') # create scalar variable with the date
	#f.variables['date'].long_name = str(mo)+'/'+str(dy)+'/'+str(yr)
	
	f.variables["Time"].long_name = "Time (beginning at midnight)"
	f.variables["Time"].units = "minutes since "+yr+"-"+str(mo)+"-"+dy+" 00:00:00"
	f.variables["Time"].calendar = "none"
	
	# what UNITS should we use here?
	f.variables["year"].long_name = "Year that data was recorded"
	f.variables["year"].units = "Julian_year"
	# The CF compliance checker gives a warning for this unit 'year'
	f.variables["month"].long_name = "Month that data was recorded (1-12)"
	f.variables["month"].units = "month"
	# The CF compliance checker gives a warning for this unit 'month'
	f.variables["day"].long_name = "Day of month that data was recorded (1-31)"
	f.variables["day"].units = "day"
	
	f.variables["latitude"].units = "degrees_north"
	f.variables["latitude"].long_name = "Station latitude in degrees"
	f.variables["latitude"].standard_name = "latitude"
	
	f.variables["longitude"].units = "degrees_east"
	f.variables["longitude"].long_name = "Station longitude in degrees"
	f.variables["longitude"].standard_name = "longitude"
	
	f.variables["altitude"].units = "m"
	f.variables["altitude"].long_name = "Station altitude above mean sea level"
	f.variables["altitude"].standard_name = "altitude"
	f.variables["altitude"].positive = "up"
	
	f.variables["best_estimate_shortwave"].units = "W m-2"
	f.variables["best_estimate_shortwave"].long_name = "Best Estimate SW, sum of direct plus diffuse"
	f.variables["best_estimate_shortwave"].comments = "Sum of direct plus diffuse if both pass QC tests, else global SW if available"

	f.variables["total_shortwave"].units = "W m-2"
	f.variables["total_shortwave"].long_name = "Total SW"
	f.variables["total_shortwave"].comments = "Total (global) SW from unshaded pyranometer"
    
	f.variables["diffuse_downwelling_shortwave_flux_in_air"].units = "W m-2"
	f.variables["diffuse_downwelling_shortwave_flux_in_air"].long_name = "Measured downwelling diffuse SW"

	f.variables["direct_downwelling_shortwave_flux_in_air"].units = "W m-2"
	f.variables["direct_downwelling_shortwave_flux_in_air"].long_name = "Measured downwelling direct SW"

	f.variables["upwelling_shortwave_flux_in_air"].units = "W m-2"
	f.variables["upwelling_shortwave_flux_in_air"].long_name = "Upwelling SW from pyranometer"

	f.variables["downwelling_longwave_flux_in_air"].units = "W m-2"
	f.variables["downwelling_longwave_flux_in_air"].long_name = "Downwelling LW from pyrgeometer"

	f.variables["upwelling_longwave_flux_in_air"].units = "W m-2"
	f.variables["upwelling_longwave_flux_in_air"].long_name = "Upwelling LW from pyrgeometer"

	f.variables["air_temperature"].units = "K"
	f.variables["air_temperature"].long_name = "Air Temperature"

	#f.variables["water_vapor_partial_pressure_in_air"].units = "mb"
	#f.variables["water_vapor_partial_pressure_in_air"].long_name = "Vapor Pressure"

	f.variables["relative_humidity"].units = "percent"
	f.variables["relative_humidity"].long_name = "Relative Humidity"

	f.variables["surface_air_pressure"].units = "mb"
	f.variables["surface_air_pressure"].long_name = "Air Pressure"
	f.variables["surface_air_pressure"].comments = "Not adjusted to sea level"

	f.variables["downwelling_longwave_case_temperature"].units = "K"
	f.variables["downwelling_longwave_case_temperature"].long_name = "Downwelling LW pyrgeometer case temperature"

	f.variables["downwelling_longwave_dome_temperature"].units = "K"
	f.variables["downwelling_longwave_dome_temperature"].long_name = "Downwelling LW pyrgeometer dome temperature"

	f.variables["upwelling_longwave_case_temperature"].units = "K"
	f.variables["upwelling_longwave_case_temperature"].long_name = "Upwelling LW pyrgeometer case temperature"

	f.variables["upwelling_longwave_dome_temperature"].units = "K"
	f.variables["upwelling_longwave_dome_temperature"].long_name = "Upwelling LW pyrgeometer dome temperature"

	f.variables["solar_zenith_angle"].units = "degrees"
	f.variables["solar_zenith_angle"].long_name = "Angle between Sun and the vertical"

	f.variables["distance_from_sun"].units = "ua"
	f.variables["distance_from_sun"].long_name = "Earth-Sun distance in astronomical units"

	f.variables["qc01"].units = "1"
	f.variables["qc01"].long_name = "Is Global SW too high?"
	f.variables["qc02"].units = "1"
	f.variables["qc02"].long_name = "Is Dif too high?"
	f.variables["qc03"].units = "1"
	f.variables["qc03"].long_name = "Is DirN too high?"
	f.variables["qc04"].units = "1"
	f.variables["qc04"].long_name = "Is SWup too high?"
	f.variables["qc05"].units = "1"
	f.variables["qc05"].long_name = "Is LWdn too low?"
	f.variables["qc06"].units = "1"
	f.variables["qc06"].long_name = "Is Lwdn too high?"
	f.variables["qc07"].units = "1"
	f.variables["qc07"].long_name = "Is LWup too low?"
	f.variables["qc08"].units = "1"
	f.variables["qc08"].long_name = "Is LWup too high?"
	f.variables["qc09"].units = "1"
	f.variables["qc09"].long_name = "Is SWup too high (normal ground cover)?"
	f.variables["qc10"].units = "1"
	f.variables["qc10"].long_name = "Is SWup too high? (snow cover)"
	f.variables["qc11"].units = "1"
	f.variables["qc11"].long_name = "Is LWdn too low? (vs Ta)"
	f.variables["qc12"].units = "1"
	f.variables["qc12"].long_name = "Is LWdn too high? (vs Ta)"
	f.variables["qc13"].units = "1"
	f.variables["qc13"].long_name = "Is LWup too low? (vs Ta)"
	f.variables["qc14"].units = "1"
	f.variables["qc14"].long_name = "Is LWup too high? (vs Ta)"
	f.variables["qc15"].units = "1"
	f.variables["qc15"].long_name = "Is LWdn too low? (vs LWup - C_15)"
	f.variables["qc16"].units = "1"
	f.variables["qc16"].long_name = "Is LWdn too high? (vs LWup + C_16)"
	f.variables["qc17"].units = "1"
	f.variables["qc17"].long_name = "LW T_c or T_d fails vs t_a +/- Limit or T_avg +/- 15.0"
	f.variables["qc18"].units = "1"
	f.variables["qc18"].long_name = "LW (T_c - T_d) fails"
	f.variables["qc19"].units = "1"
	f.variables["qc19"].long_name = "T_a fails"  

	f.variables["g_flag"].units = "1"
	f.variables["g_flag"].long_name = "Type of IR loss correction applied to global SW"
	f.variables["g_flag"].comments = "none/full dry/full moist/detector dry/det moist (0/1/2/3/4)"

	f.variables["d_flag"].units = "1"
	f.variables["d_flag"].long_name = "Type of IR loss correction applied to diffuse SW"
	f.variables["d_flag"].comments = "none/full dry/full moist/detector dry/det moist (0/1/2/3/4)"
	
	f.variables["downwelling_shortwave_flux_in_air_assuming_clear_sky"].units = "W m-2"
	f.variables["downwelling_shortwave_flux_in_air_assuming_clear_sky"].long_name = "Estimated clear-sky total downwelling SW"
	f.variables["downwelling_shortwave_flux_in_air_assuming_clear_sky"].comments = "Estimated from user-set power law coefficient settings in configuration file"

	f.variables["ir_loss_correction_to_diffuse_shortwave"].units = "1"
	f.variables["ir_loss_correction_to_diffuse_shortwave"].long_name = "IR loss correction applied to diffuse (shaded) SW measurements"

	f.variables["ir_loss_correction_to_global_shortwave"].units = "1"
	f.variables["ir_loss_correction_to_global_shortwave"].long_name = "IR loss correction applied to global (unshaded) SW measurements"

	f.variables["downwelling_erythemal_uvb_flux_in_air"].units = "mW m-2"
	f.variables["downwelling_erythemal_uvb_flux_in_air"].long_name = "Erythemal UVB"
	f.variables["downwelling_erythemal_uvb_flux_in_air"].comments = "Instrument is weighted by the erythemal action spectrum"

	f.variables["surface_downwelling_photosynthetic_radiative_flux_in_air"].units = "W m-2"
	f.variables["surface_downwelling_photosynthetic_radiative_flux_in_air"].long_name = "Photosynthetically Active Radiation"
	f.variables["surface_downwelling_photosynthetic_radiative_flux_in_air"].comments = "Solar radiation from 400 to 700 nanometers"

	f.variables["wind_speed"].units = "m s-1"
	f.variables["wind_speed"].long_name = "Wind speed"

	f.variables["wind_from_direction"].units = "degree"
	f.variables["wind_from_direction"].long_name = "Direction of wind in degrees"
    
	return f

if __name__ == '__main__':
    main(sys.argv[1:])
