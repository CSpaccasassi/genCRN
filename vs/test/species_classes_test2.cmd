@echo off
FOR /L %%R IN (1,1,4) DO (
	echo Generate CRNs with no species classes
	..\x64\Release\genCRN.exe -n2 crn_2_%%R.txt > expected2%%R.txt
	
	echo Generate CRNs with species class [1;1]
	..\x64\Release\genCRN.exe -n2 -s1;1 crn_2_%%R.txt > actual2%%R_111.txt	
	..\TestGenCRN\bin\Release\TestGenCRN.exe -s 1 1 -e expected2%%R.txt -a actual2%%R_111.txt

	echo --------------------------------------
)