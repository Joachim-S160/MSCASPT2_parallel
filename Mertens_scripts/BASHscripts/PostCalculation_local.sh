#!/usr/bin/env bash

module swap cluster/delcatty

CurrDir=$PWD

DissCurveDiffTPL_SC.sh
MakeVIBROTsh_SC_local.sh
./MakeVibrotAllTPL_SC_local.sh
PlotDiffTPL.sh

echo "Do you also want to determine the intruder states? (y/n)"
read INTRUDER


if [[ $INTRUDER == "y" ]]; then
    for i in $CurrDir/tpl*
    do
        cd $i
        echo -n 'Checking intuders in '
        echo $i
        MakeCASPT2filesAllTrans_SC.sh
    done
fi
