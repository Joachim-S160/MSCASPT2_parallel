#!/usr/bin/env bash

echo "Give the upper limit for the denominator: "
read DENOM
export DENOM

echo "Give the lower limit for the contribution: "
read CONTRI
export CONTRI

CurrDir=$PWD;

for i in $CurrDir/tpl*;
do
    echo "Getting intuders in $i | DENOM = ${DENOM} and CONTRI = ${CONTRI}"
    cd $i
    MakeCASPT2filesAllTrans_SC_WEAK_choice.sh
done

