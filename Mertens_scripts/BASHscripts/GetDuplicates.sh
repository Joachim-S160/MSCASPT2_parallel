#!/usr/bin/env bash

echo "Please enter the name of the molecule: "
read Molecule

echo "Enter 0 if number of electrons is even, 1 if number of electrons are odd"
read Elec

#NrOfTPL=$(cat ReadableOverviewTPL.txt | sort -n -k1 | tail -1 | awk '{print $1}')

if [ "$Elec" == "0" ]; then
    rm CleanOverview.txt
    cat ReadableOverviewTPL.txt | sort -n  -k6 -k7 -k8 -k9 -k10 -k11 -k12 > SortedOverview.txt
    mapfile -t OV < SortedOverview.txt
    SEQ=$(seq 1 ${#OV[@]})

    printf "%-5s%-3s%-12s%-3s%-12s%-9.3s%-9.5s%-9s%-6s%-9s%-21s%-4.4s" ${OV[0]}  >> CleanOverview.txt
    printf "\n" >> CleanOverview.txt
    echo "______________________________________________________________________________________________________________________________" >> CleanOverview.txt
    for i in $SEQ
    do
        j=$(($i + 1))
        A=$(echo ${OV[$i]} | awk '{print $4 $5 $6 $7 $8 $9 $10 $11 $12}')
        B=$(echo ${OV[$j]} | awk '{print $4 $5 $6 $7 $8 $9 $10 $11 $12}')
        if  [ "$A" == "$B" ]; then
            printf "%-5s%-15s%-15s%-3s%-6s%-3s%-6s%-3s%-6s%-6s%-3s%-6s%-3s%-3s%-3s%-3s%-3s%-6s%s" ${OV[$i]} >> CleanOverview.txt
            printf "\n" >> CleanOverview.txt
        else
            printf "%-5s%-15s%-15s%-3s%-6s%-3s%-6s%-3s%-6s%-6s%-3s%-6s%-3s%-3s%-3s%-3s%-3s%-6s%s" ${OV[$i]} >> CleanOverview.txt
            printf "\n" >> CleanOverview.txt
            echo "______________________________________________________________________________________________________________________________" >> CleanOverview.txt
        fi
    done

    cat ReadableOverviewTPL.txt  | sort -u -k4,12 | sed \$d | awk {'print $1'} | sort -n > UniqueTPL.txt
    mapfile -t UnTPL < UniqueTPL.txt
    for T in ${UnTPL[@]}
    do
	    mkdir tpl${T}
	    cp TPL${T} tpl${T}/template
	    cp ${Molecule}.sh tpl${T}
        sed -i "s/JOBNAME/${Molecule}TPL${T}/g" "tpl${T}/${Molecule}.sh"
    done
    rm -r tpl
fi

if [ "$Elec" == "1" ]; then
    rm CleanOverview.txt
    cat ReadableOverviewTPL.txt | sort -n  -k6 -k7 -k8 -k9 -k10 -k11 -k12 -k13 -k14 > SortedOverview.txt
    mapfile -t OV < SortedOverview.txt
    SEQ=$(seq 1 ${#OV[@]})

    printf "%-5s%-3s%-9s%-3s%-9s%-4s%-5s%-4s%-5s%-9.5s%-9s%-6s%-9s%-21s%-4.4s" ${OV[0]} >> CleanOverview.txt
    printf "\n" >> CleanOverview.txt
    echo "_________________________________________________________________________________________________________" >> CleanOverview.txt
    for i in $SEQ
    do
        j=$(($i + 1))
        A=$(echo ${OV[$i]} | awk '{print $4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14}')
        B=$(echo ${OV[$j]} | awk '{print $4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14}')
        if  [ "$A" == "$B" ]; then
            printf "%-5s%-12s%-12s%-3s%-6s%-3s%-6s%-3s%-6s%-3s%-6s%-6s%-3s%-6s%-3s%-3s%-3s%-3s%-3s%-6s%s" ${OV[$i]} >> CleanOverview.txt
            printf "\n" >> CleanOverview.txt
        else
            printf "%-5s%-12s%-12s%-3s%-6s%-3s%-6s%-3s%-6s%-3s%-6s%-6s%-3s%-6s%-3s%-3s%-3s%-3s%-3s%-6s%s" ${OV[$i]} >> CleanOverview.txt
            printf "\n" >> CleanOverview.txt
            echo "_________________________________________________________________________________________________________" >> CleanOverview.txt
        fi
    done

    cat ReadableOverviewTPL.txt  | sort -u -k4,14 | sed \$d | awk {'print $1'} | sort -n > UniqueTPL.txt
    mapfile -t UnTPL < UniqueTPL.txt
    for T in ${UnTPL[@]}
    do
        mkdir tpl${T}
        cp TPL${T} tpl${T}/template
        cp ${Molecule}.sh tpl${T}
        sed -i "s/JOBNAME/${Molecule}TPL${T}/g" "tpl${T}/${Molecule}.sh"
    done
    rm -r tpl
fi
