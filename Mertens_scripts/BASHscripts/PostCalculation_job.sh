#!/bin/sh
#
#
#PBS -N PostCalc
#PBS -l walltime=71:58:00
#PBS -l nodes=1:ppn=1
#PBS -o PostCalc.stdout
#PBS -e PostCalc.stderr

module swap cluster/delcatty
module load Python/3.5.2-intel-2016b
module load matplotlib/1.5.3-intel-2016b-Python-3.5.2

export CurrDir=$PBS_O_WORKDIR



DissCurveDiffTPL_SC.sh
MakeVIBROTsh_SC_local.sh
./MakeVibrotAllTPL_SC_local.sh
PlotDiffTPL.sh

INTRUDER="y"


if [[ $INTRUDER == "y" ]]; then
    for i in $CurrDir/tpl*
    do
        cd $i
        echo -n 'Checking intuders in '
        echo $i
        MakeCASPT2filesAllTrans_SC.sh
    done
fi
