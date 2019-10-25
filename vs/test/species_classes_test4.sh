#!/bin/bash
for NR in {1..3}
do
	echo "Generate CRNs with no species classes"
	../x64/Release/genCRN.exe -n4 crn_4_$NR.txt > expected4$NR.txt
	
	echo "Generate CRNs with species class [2;2]"
	../x64/Release/genCRN.exe -n4 -'s2;2' crn_4_$NR.txt > actual4$NR_22.txt
	../TestGenCRN/bin/Release/TestGenCRN.exe -s 2 2 -e expected4$NR.txt -a actual4$NR_22.txt
	echo "---"
	
	echo "Generate CRNs with species class [1;3]"
	../x64/Release/genCRN.exe -n4 -'s1;3' crn_4_$NR.txt > actual4$NR_13.txt
	../TestGenCRN/bin/Release/TestGenCRN.exe -s 1 3 -e expected4$NR.txt -a actual4$NR_13.txt
	echo "---"
	
	echo "Generate CRNs with species class [1;1;2]"
	../x64/Release/genCRN.exe -n4 -'s1;1;2' crn_4_$NR.txt > actual4$NR_112.txt
	../TestGenCRN/bin/Release/TestGenCRN.exe -s 1 1 2 -e expected4$NR.txt -a actual4$NR_112.txt
	echo "---"
	
	echo "Generate CRNs with species class [1;1;1;1]"
	../x64/Release/genCRN.exe -n4 -'s1;1;1;1' crn_4_$NR.txt > actual4$NR_1111.txt	
	../TestGenCRN/bin/Release/TestGenCRN.exe -s 1 1 1 1 -e expected4$NR.txt -a actual4$NR_1111.txt	
	echo "--------------------------------------"
done