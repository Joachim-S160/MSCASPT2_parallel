#!/usr/bin/env bash
#
#
#PBS -N INTRUDERS_MOL
#PBS -l walltime=71:58:00
#PBS -l nodes=1:ppn=1
#PBS -o job-INTRUDERS_MOL.stdout
#PBS -e job-INTRUDERS_MOL.stderr

export CurrDir=$PBS_O_WORKDIR

for i in $CurrDir/tpl*
do
    cd $i
    MakeCASPT2filesAllTrans_SC.sh
done

