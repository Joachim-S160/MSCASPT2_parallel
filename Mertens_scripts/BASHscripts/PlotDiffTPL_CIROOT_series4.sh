#!/usr/bin/env bash

module load Python/3.5.2-intel-2016b 
module load matplotlib/1.5.3-intel-2016b-Python-3.5.2

CurrDir=$PWD
echo "Enter the name of the molecule: "
read Molecule

for TPLNR in  18 21 22 23 24 50 52 54 55 56 5 6 63 64 #  1 17 18 19 2 20 22 24 50 51 52 54 56 59 6 64 #49 50 51 52 56 58 73 74 80
do
    CIROOTS=()
    
    for ciroot in CIROOT_FewStates CIROOT_AllSingleComb CIROOT_1Li3Po CIROOT3 CIROOT6 CIROOT1 CIROOT5 CIROOT2 CIROOT4 CIROOT7 CIROOT8 CIROOT9
    do
#        Cut2=${file#*CIROOT*/}
#        Molecule=${Cut2%tpl*}
#        echo $Molecule
        file=$CurrDir/${ciroot}/${Molecule}tpl${TPLNR}sorassiSORTED.txt
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
