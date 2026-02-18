#!/usr/bin/env bash

module load Python/3.5.2-intel-2016b 
module load matplotlib/1.5.3-intel-2016b-Python-3.5.2


LIST=()

echo "Enter the name of the molecule: "
read Molecule

echo "Enter desired TPLnrs separated by space: "
read -a TPLS

echo "Enter outputname for pdf: "
read PDFname

PlotDiffTPL_chosen.py $Molecule $PDFname ${TPLS[@]}

