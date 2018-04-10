#!/usr/bin/env bash
# ----------------------------------------
# PIPELINE --- DAT to NetCDF
# 
# Parameters: 
# $1 - choose: daily monthly or yearly
# $2 - Folder containing DAT files
#
# example usage: bash completePipeline.sh daily ../bondville/
# 
# NOTE: add_metadata.py must be 
# updated to reflect the specific dataset
#
# Ezra Huscher, September 2017
# ezra.huscher@noaa.gov
# ----------------------------------------

start=`date +%s`
folderstruct=""

### resulting NetCDFs: daily, monthly, yearly
if [ $1 == 'daily' ]
then
  folderstruct="*/"
elif [ $1 == 'monthly' ]
then
  folderstruct="*/"
elif [ $1 == 'yearly' ]
then
  folderstruct=""
else
  echo 'First parameter must be daily, monthly, or yearly.'
  exit 0
fi

which_py_script="csv_to_nc_$1.py"

# How many files are we dealing with?
count=$(find $2 -type f -name '*.qDAT' -o -name '*.qdat' -o -name '*.dat' -o -name '*.qadj' -o -name '*.lw1' | wc -l)
count=$(echo "$count" | xargs) # trim white space

### Create nc directory, if it doesn't exit
#mkdir -p $2/nc
#rm $2/nc/*.nc

### Clean up file, if necessary
# Remove x commented lines from beginning
# tail -n +18 brw_F12_Day.dat > brw_F12_Day_rdy.dat

### Converts all DATs to CSVs
# all at once: bash _dat_to_csv.sh Data/Bondville_IL/1995/
# one at a time: python dat_to_csv.py ../brw_F12_Day_rdy.dat
# NOTE: Must adjust strptime for dat filename format
echo "Converting $count DAT files to CSV files..."
bash _dat_to_csv.sh $2

### Replace all csv files in Data folder: -9999.9 with NaN
# MacOS seems to require a backup file (.bak)
echo 'Replacing any -9999.9 values with NaN...'
find $2 -name '*.csv' -print0 | xargs -0 sed -i '.bak' 's/-9999.9000/NaN/g'
find $2 -name '*.csv' -print0 | xargs -0 sed -i '.bak' 's/-9999.900/NaN/g'
find $2 -name '*.csv' -print0 | xargs -0 sed -i '.bak' 's/-9999.90/NaN/g'
find $2 -name '*.csv' -print0 | xargs -0 sed -i '.bak' 's/-9999.9/NaN/g'
find $2 -name '*.csv' -print0 | xargs -0 sed -i '.bak' 's/-9999.0/NaN/g'
#rm $2/*.bak

# Delete .bak files in all subdirectories
find $2 -name '*.bak' -delete

### Convert to NetCDF format
# This will call add_metadata.py
echo "Converting $count CSV files to NetCDF files..."
python $which_py_script $2/*/*.csv

# Delete .csv files in all subdirectories
find $2 -name '*.csv' -delete

### Compress the NetCDFs to desired level
echo 'Compressing NetCDF files...'
# (d1-d9) is the compression level
i=0
for f in $2/*/*_uncompressed.nc
do
    nccopy -d9 $f ${f//_uncompressed/}  
    rm $f
    i=$((i+1))
done

echo "$i NetCDF files generated successfully!"

end=`date +%s`
runtime=$((end-start))
echo "Runtime was $runtime seconds."
exit 0

### run NCEI compliance checker or https://data.ioos.us/compliance/index.html
compliance-checker -t cf -v brw_F12_Day_rdy.nc
