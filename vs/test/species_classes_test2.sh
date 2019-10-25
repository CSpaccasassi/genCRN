#!/bin/bash
for NR in {1..4} 
do
	echo "Generate CRNs with no species classes"
	../x64/Release/genCRN.exe -n2 crn_2_$NR.txt > expected2$NR.txt
	
	echo "Generate CRNs with species class [1;1]"
	../x64/Release/genCRN.exe -n2 -'s1;1' crn_2_$NR.txt > actual2$NR_111.txt	
	../TestGenCRN/bin/Release/TestGenCRN.exe -s 1 1 -e expected2$NR.txt -a actual2$NR_111.txt
	echo "--------------------------------------"
done