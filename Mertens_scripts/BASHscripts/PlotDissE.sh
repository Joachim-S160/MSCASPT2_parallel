#!/usr/bin/env bash


#Transform data to usable array
echo -n "Give the molecule's name: "
read Molecule

FileD0="${Molecule}D0.txt"
FileDe="${Molecule}De.txt"

echo $FileDe

cat $FileD0 | sed 's/tpl//g' | sed 's/\/.*://g' | awk {'print $1 "\t" $3'} > D0vsTPL.txt
cat $FileDe | sed 's/tpl//g' | sed 's/\/.*://g' | awk {'print $1 "\t" $3'} > DevsTPL.txt

module load Python/3.5.1-intel-2016a
module load matplotlib/1.5.1-intel-2016a-Python-3.5.1

PlotD0vsTPL.py
