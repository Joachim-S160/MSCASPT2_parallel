#!/usr/bin/env bash

module load Python/3.5.2-intel-2016b 
module load matplotlib/1.5.3-intel-2016b-Python-3.5.2

CurrDir=$PWD

echo "Enter TPLnr of interest: "
read TPLNR

CIROOTS=()

for file in $CurrDir/CIROOT5*/*tpl${TPLNR}sorassiSORTED*
do
    Cut2=${file#*CIROOT*/}
    Molecule=${Cut2%tpl*}
    echo $Molecule
    Cut1=${file#$CurrDir/CIROOT}
    CIROOTnr=${Cut1%/*}
    echo $CIROOTnr
    CIROOTS+=(${CIROOTnr})
done

PlotDiffTPL_ImShift.py $Molecule $TPLNR ${CIROOTS[@]}

