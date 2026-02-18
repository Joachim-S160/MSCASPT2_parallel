#!/bin/sh

Atom1=ATOM1
Atom2=ATOM2
Molecule=MOLECULE

module load OpenMPI/1.8.1-GCC-4.8.3

module load Molcas/8.0sp1_OpenMPI-1.8.1_CentOS-6.6

module load Python/3.5.2-intel-2016b

#### Set MOLCAS location
export MOLCAS_MEM=2000
export MOLCAS_LICENSE=/user/scratch/gent/gvo000/gvo00003/vsc41682/Licenses

#### Set MOLCAS in running directory
export CurrDir=$PWD
export WorkDir=$CurrDir
export OutputDir=$CurrDir

rm $CurrDir/*/*VIBROT.inp

for file in $CurrDir/*sorassi.txt
do
	Cut1=${file%sorassi*}
	TPL=${Cut1#*tpl}
    echo $TPL
    echo $file
	sort -n -k1 $file > "${Molecule}tpl${TPL}sorassiSORTED.txt" 
	MakeVibRotAllTPL.py $Molecule $Atom1 $Atom2 $TPL $CurrDir/
done

for file in $CurrDir/*/*VIBROT.inp
do
	Out=${file%.inp}
	molcas $file > ${Out}.out 
done


