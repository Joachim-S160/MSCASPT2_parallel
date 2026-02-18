#!/usr/bin/env bash

cp /user/scratch/gent/gvo000/gvo00003/vsc41682/PhD/MOLCAS_20161202/Amsterdam/BAZIS/BASHscripts/MakeVibrotAllTPL_SC_local.sh .

sed -i "s/ATOM1/${ATOM1}/g" MakeVibrotAllTPL_SC_local.sh
sed -i "s/ATOM2/${ATOM2}/g" MakeVibrotAllTPL_SC_local.sh
sed -i "s/MOLECULE/${Molecule}/g" MakeVibrotAllTPL_SC_local.sh
