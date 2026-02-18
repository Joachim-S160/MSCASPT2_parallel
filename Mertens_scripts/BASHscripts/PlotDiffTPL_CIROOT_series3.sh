#!/usr/bin/env bash

module load Python/3.5.2-intel-2016b 
module load matplotlib/1.5.3-intel-2016b-Python-3.5.2

CurrDir=$PWD
echo "Enter the name of the molecule: "
read Molecule

for TPLNR in 17 18 2 20 22 24 49 50 51 52 54 56 6 60 64 1 19 59
do
    CIROOTS=()

    for file in $CurrDir/CIROOT*/*tpl${TPLNR}sorassiSORTED* # $CurrDir/CIROOT3S/*tpl${TPLNR}sorassiSORTED*
    do
#        Cut2=${file#*CIROOT*/}
#        Molecule=${Cut2%tpl*}
#        echo $Molecule
        Cut1=${file#$CurrDir/CIROOT}
        CIROOTnr=${Cut1%/*}
        echo $CIROOTnr
        CIROOTS+=(${CIROOTnr})
    done

    PlotDiffTPL_CIROOT.py $Molecule $TPLNR ${CIROOTS[@]}
done
: '

for TPLNR in 1 19 59
do
    CIROOTS=()

    for file in $CurrDir/CIROOT{1,2}/*tpl${TPLNR}sorassiSORTED*
    do
#        Cut2=${file#*CIROOT*/}
#        Molecule=${Cut2%tpl*}
#        echo $Molecule
        Cut1=${file#$CurrDir/CIROOT}
        CIROOTnr=${Cut1%/*}
        echo $CIROOTnr
        CIROOTS+=(${CIROOTnr})
    done

    PlotDiffTPL_CIROOT.py $Molecule $TPLNR ${CIROOTS[@]}
done

for TPLNR in 54
do
    CIROOTS=()

    for file in $CurrDir/CIROOT{1,2,3,4,5,6,7,8,9,10}/*tpl${TPLNR}sorassiSORTED*
    do
#        Cut2=${file#*CIROOT*/}
#        Molecule=${Cut2%tpl*}
#        echo $Molecule
        Cut1=${file#$CurrDir/CIROOT}
        CIROOTnr=${Cut1%/*}
        echo $CIROOTnr
        CIROOTS+=(${CIROOTnr})
    done

    PlotDiffTPL_CIROOT.py $Molecule $TPLNR ${CIROOTS[@]}
done

for TPLNR in 20
do
    CIROOTS=()

    for file in $CurrDir/CIROOT{6,7,8,10}/*tpl${TPLNR}sorassiSORTED*
    do
#        Cut2=${file#*CIROOT*/}
#        Molecule=${Cut2%tpl*}
#        echo $Molecule
        Cut1=${file#$CurrDir/CIROOT}
        CIROOTnr=${Cut1%/*}
        echo $CIROOTnr
        CIROOTS+=(${CIROOTnr})
    done

    PlotDiffTPL_CIROOT.py $Molecule $TPLNR ${CIROOTS[@]}
done




for TPLNR in 54
do
    CIROOTS=()

    for file in $CurrDir/CIROOT{6,7,8,9}/*tpl${TPLNR}sorassiSORTED*
    do
#        Cut2=${file#*CIROOT*/}
#        Molecule=${Cut2%tpl*}
#        echo $Molecule
        Cut1=${file#$CurrDir/CIROOT}
        CIROOTnr=${Cut1%/*}
        echo $CIROOTnr
        CIROOTS+=(${CIROOTnr})
    done

    PlotDiffTPL_CIROOT.py $Molecule $TPLNR ${CIROOTS[@]}
done
'
