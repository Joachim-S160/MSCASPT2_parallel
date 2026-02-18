#!/usr/bin/env bash

echo "Enter the name of the molecule: "
read Molecule

CurrDir=$PWD
for i in $CurrDir/tpl*
do
    cd $i;
    qsub ${Molecule}.sh
done
