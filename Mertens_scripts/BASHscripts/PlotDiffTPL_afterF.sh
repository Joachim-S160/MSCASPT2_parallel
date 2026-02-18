#!/usr/bin/env bash

module load Python/3.5.2-intel-2016b 
module load matplotlib/1.5.3-intel-2016b-Python-3.5.2


LIST=()

for file in *FsorassiSORTED*
do
    Cut1=${file#*tpl}
    TPLnr=${Cut1%FsorassiSORTED.txt}
    LIST+=(${TPLnr})
    Molecule=${file%tpl*}
done

PlotDiffTPL_afterF.py $Molecule ${LIST[@]}

