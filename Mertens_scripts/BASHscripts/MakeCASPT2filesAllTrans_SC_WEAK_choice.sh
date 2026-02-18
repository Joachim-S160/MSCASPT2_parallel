#!/usr/bin/env bash

module load Python/3.5.1-intel-2016a

OutputFiles=()

for i in */out_*
do
    if [[ "$(grep 'executing module RASSI' $i)" != "" ]]
    then    
        OutputFiles+=(${i})
    fi
done

#CASPT2Lines=()


for OUT in ${OutputFiles[@]}
do
    egrep -n "Start Module: casp|Stop Module:  casp" $OUT | cut -f1 -d: > SectionLinesWEAK.txt
    #i=0
    #while [ $i -lt ${#LineNumbers[@]} ]
    #do
    #echo ${#CASPT2Lines[@]}
    egrep -n 'ATVX     1|ATVX     2|BVATP    1|BVATP    2|BVATM    1|BVATM    2|VJTU     1|VJTU     2|VJTIP    1|VJTIP    2|VJTIM    1|VJTIM    2|AIVX     1|AIVX     2|VJAIP    1|VJAIP    2|BJATP    1|BJATP    2|BJATM    1|BJATM    2|BJAIP    1|BJAIP    2|BJAIM    1|BJAIM    2' $OUT | sed $'s/:/\t/g' | sed -e 's/\<BVATM\>//g' | sed -e 's/\<BVATP\>//g' | sed -e 's/\<VJTIP\>//g' |sed -e 's/\<VJTIM\>//g' |sed -e 's/\<AIVX\>//g' |sed -e 's/\<VJAIP    1\>//g' |sed -e 's/\<VJAIP    2\>//g' |sed -e 's/\<VJAIM    1\>//g' |sed -e 's/\<VJAIM    2\>//g' |sed -e 's/\<BJATP    1\>//g' |sed -e 's/\<BJATP    2\>//g' |sed -e 's/\<BJATM    1\>//g' |sed -e 's/\<BJATM    2\>//g' |sed -e 's/\<BJAIP    1\>//g' | sed -e 's/\<BJAIP    2\>//g' |sed -e 's/\<BJAIM    1\>//g' |sed -e 's/\<BJAIM    2\>//g' | awk -v OFS='\t' '{print $1, $6, $9}' > DenominatorsWEAK.txt
    
    
    

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
    
    grep '++ CASPT2' $OUT | awk '{print $4}' > RootNrWEAK.txt 
    egrep -n '\+\+ CASPT2|CASPT2 PROPERTY'  $OUT | cut -f1 -d: > RootLinesWEAK.txt
    #CaseLines=($(egrep -n 'CASE' $OUT | cut -f1 -d:))
    grep 'State symmetry' $OUT | awk '{print $3}' > SymmWEAK.txt
    grep 'Spin quantum' $OUT | awk '{print $4}' > SpinWEAK.txt
    grep 'Number of root(s) av' $OUT | awk '{print $5}' > NumberOfRootsWEAK.txt
    echo ${OUT#*_} > OutputFilesWEAK.txt
    #done
    GetIntruderLinesWEAK_choice.py ${DENOM} ${CONTRI}
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
