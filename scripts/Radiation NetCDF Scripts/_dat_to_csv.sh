#!/usr/bin/env bash
FILES=$(find $1 -type f -iname '*.qDAT')
for f in $FILES
do
  python dat_to_csv.py $f
done

FILES=$(find $1 -type f -iname '*.qadj')
for f in $FILES
do
  python dat_to_csv.py $f
done

FILES=$(find $1 -type f -iname '*.qdat')
for f in $FILES
do
  python dat_to_csv.py $f
done

FILES=$(find $1 -type f -iname '*.lw1')
for f in $FILES
do
  python dat_to_csv.py $f
done
