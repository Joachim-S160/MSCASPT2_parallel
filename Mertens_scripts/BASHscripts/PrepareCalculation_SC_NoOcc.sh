#!/usr/bin/env bash
CurrentDir=$PWD
module load Python/3.5.2-intel-2016b
TpLmAkEr_SC_NoOcc.py
chmod 700 SubmitFiles.sh
chmod 700 StartCalc.sh
GetDuplicates.sh        
echo 'Please enter the name of the submit cluster: '
read Cluster
module swap cluster/${Cluster}
for i in tpl*
do
    cp $CurrentDir/SubmitFiles.sh $CurrentDir/$i
    cd $CurrentDir/$i
    ./SubmitFiles.sh
done



