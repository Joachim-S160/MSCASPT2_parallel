#!/usr/bin/env bash


rm TEMPCASSCF.txt TEMPCASPT2.txt TEMPMSCASPT2.txt TEMPRASSI.txt TEMPSORASSI.txt TEMPHF.txt

for DIR in $CurrDir/tpl*
do
	for i in $DIR/*/out*
	do
        	if grep -q 'executing module RASSI' $i
        	then
                echo -n $i | sed -e 's/.*out_//g' >> TEMPCASSCF.txt
                echo -ne ' \t ' >> TEMPCASSCF.txt
                grep 'RASSCF root number  1' $i | awk {'print $8'} | sort -n -k 1 | head -n 1  >> TEMPCASSCF.txt
                echo -n $i | sed -e 's/.*out_//g' >> TEMPCASPT2.txt
                echo -ne ' \t ' >> TEMPCASPT2.txt
                grep ' CASPT2 Root  1' $i | awk {'print $7'} | sort -n -k 1 | head -n 1  >> TEMPCASPT2.txt
                echo -n $i | sed -e 's/.*out_//g' >> TEMPMSCASPT2.txt
                echo -ne ' \t ' >> TEMPMSCASPT2.txt
                grep 'MS-CASPT2 Root  1' $i | awk {'print $7'} | sort -n -k 1 | head -n 1  >> TEMPMSCASPT2.txt
		echo -n $i | sed -e 's/.*out_//g' >> TEMPHF.txt
		echo -ne ' \t ' >> TEMPHF.txt
		grep 'Total SCF energy' $i | awk {'print $5'} >> TEMPHF.txt
    		fi

		    if grep -q 'Happy' $i
            then
                    echo -n $i | sed -e 's/.*out_//g' >> TEMPRASSI.txt
                    echo -ne ' \t ' >> TEMPRASSI.txt
                    grep ' RASSI State' $i | sed 's/.*-/-/g' | sort -n -k 1 | head -n 1 >> TEMPRASSI.txt
                    echo -n $i | sed -e 's/.*out_//g' >> TEMPSORASSI.txt
                    echo -ne ' \t ' >> TEMPSORASSI.txt
                    grep 'SO-RASSI State  1' $i | awk {'print $7'} >> TEMPSORASSI.txt
            fi

	done
	sort -n -r -k 1 < TEMPHF.txt > ${Molecule}${DIR#${CurrDir}/}HF.txt
	sort -n -r -k 1 < TEMPCASSCF.txt > ${Molecule}${DIR#${CurrDir}/}casscf.txt
	sort -n -r -k 1 < TEMPCASPT2.txt > ${Molecule}${DIR#${CurrDir}/}caspt2.txt
	sort -n -r -k 1 < TEMPMSCASPT2.txt > ${Molecule}${DIR#${CurrDir}/}mscaspt2.txt
	sort -n -r -k 1 < TEMPRASSI.txt > ${Molecule}${DIR#${CurrDir}/}rassi.txt
	sort -n -r -k 1 < TEMPSORASSI.txt > ${Molecule}${DIR#${CurrDir}/}sorassi.txt
	rm TEMPCASSCF.txt TEMPCASPT2.txt TEMPMSCASPT2.txt TEMPRASSI.txt TEMPSORASSI.txt TEMPHF.txt
done

