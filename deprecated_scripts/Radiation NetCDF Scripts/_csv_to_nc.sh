#!/usr/bin/env bash
FILES=$(find $1 -type f -iname '*.csv')
python csv_to_nc_monthly.py $FILES
