#!/bin/bash
for NR in {1..4} 
do
	echo "Generate CRNs with no species classes"
	../x64/Release/genCRN.exe -n3 crn_3_$NR.txt > expected3$NR.txt
	
	echo "Generate CRNs with species class [1;2]"
	../x64/Release/genCRN.exe -n3 -'s1;2' crn_3_$NR.txt > actual3$NR_12.txt
	../TestGenCRN/bin/Release/TestGenCRN.exe -s 1 2 -e expected3$NR.txt -a actual3$NR_12.txt
	
	if (($NR < 4))
	then
		echo "---"
		echo "Generate CRNs with species class [1;1;1]"
		../x64/Release/genCRN.exe -n3 -s'1;1;1' crn_3_$NR.txt > actual3$NR_111.txt	
		../TestGenCRN/bin/Release/TestGenCRN.exe -s 1 1 1 -e expected3$NR.txt -a actual3$NR_111.txt
	fi
	
	echo --------------------------------------
done