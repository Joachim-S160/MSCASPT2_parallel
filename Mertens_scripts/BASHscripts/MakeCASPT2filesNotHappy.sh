#!/usr/bin/env bash

module load Python/3.5.1-intel-2016a

OutputFiles=()

for i in out_*
do
    if [[ "$(grep 'executing module RASSI' $i)" != "" ]]
    then    
        OutputFiles+=(${i#*_})
    fi
done

#CASPT2Lines=()


for OUT in ${OutputFiles[@]}
do
    egrep -n "Start Module: casp|Stop Module:  casp" out_$OUT | cut -f1 -d: > SectionLines.txt
    #i=0
    #while [ $i -lt ${#LineNumbers[@]} ]
    #do
    #echo ${#CASPT2Lines[@]}
    egrep -n 'ATVX     1|ATVX     2|BVATP    1|BVATP    2|BVATM    1|BVATM    2' out_$OUT | sed $'s/:/\t/g' | sed -e 's/\<BVATM\>//g' | sed -e 's/\<BVATP\>//g' | awk -v OFS='\t' '{print $1, $6, $9}' > Denominators.txt
    #A=($(egrep -n 'ATVX     1|ATVX     2|BVATP    1|BVATP    2|BVATM    1|BVATM    2|CASE| Report on' out_$OUT | cut -f1 -d:))
    #SectionLines=()
    #SectionLines+=(${A[0]})
    #for i in "${!A[@]}"
    #do
    #     First=${A[$i]}
    #     Second=${A[$(($i+1))]}
    #     if [[ $(($Second-$First)) -gt 1 ]]
    #     then 
    #         SectionLines+=($First)
    #         SectionLines+=($Second)
    #     fi
    #done
    #SectionLines+=(${A[${#A[@]}-1]})
    
    grep '++ CASPT2' out_$OUT | awk '{print $4}' > RootNr.txt 
    egrep -n '\+\+ CASPT2|CASPT2 PROPERTY'  out_$OUT | cut -f1 -d: > RootLines.txt
    #CaseLines=($(egrep -n 'CASE' out_$OUT | cut -f1 -d:))
    grep 'State symmetry' out_$OUT | awk '{print $3}' > Symm.txt
    grep 'Spin quantum' out_$OUT | awk '{print $4}' > Spin.txt
    grep 'Number of root(s) av' out_$OUT | awk '{print $5}' > NumberOfRoots.txt
    echo $OUT > OutputFiles.txt
    #done
    GetIntruderLines.py
done

#        echo 'Making the TEST files'
#        sed -n "${LineNumbers[i]},${LineNumbers[i+1]}p" outNr_$OUT | grep 'Spin quantum' | sed 's/:/ /g' >> TEST_$OUT
#        sed -n 
#        sed -n "${LineNumbers[i]},${LineNumbers[i+1]}p" outNr_$OUT | grep 'State symmetry' | sed 's/:/ /g' >> TEST_$OUT
#        sed -n "${LineNumbers[i]},${LineNumbers[i+1]}p" outNr_$OUT | sed -n -e  '/CONTRIB/,/--/ p' | sed 's/:/ /g' >> TEST_$OUT
#        i=$((i+2))
#    done
#    rm outNr_$OUT
#    echo 'Start reading line by line'
#    while IFS='' read -r line || [[ -n "$line" ]]
#    do
# 
        #num1=1
        #num2=2
#        LINE=($line)
        #Case=${LINE[1]}
        #Denominator=${LINE[5]}
        #Contribution=${LINE[8]}
#        if [[ "${LINE[1]}" == "ATVX" ]] && (( $(echo "${LINE[5]} < 0.1" |bc -l) )) && (( $(echo "${LINE[8]} > 0.001" |bc -l) )) 
#        then
#            echo ${LINE[@]}
#        fi
#    done < TEST
#done
#COMMENT1
