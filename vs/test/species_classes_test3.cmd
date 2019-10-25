@echo off
FOR /L %%R IN (1,1,4) DO (
	echo Generate CRNs with no species classes
	..\x64\Release\genCRN.exe -n3 crn_3_%%R.txt > expected3%%R.txt
	
	echo Generate CRNs with species class [1;2]
	..\x64\Release\genCRN.exe -n3 -s1;2 crn_3_%%R.txt > actual3%%R_12.txt
	..\TestGenCRN\bin\Release\TestGenCRN.exe -s 1 2 -e expected3%%R.txt -a actual3%%R_12.txt
	
	if %%R LSS 4 (
		echo ---
		echo Generate CRNs with species class [1;1;1]
		..\x64\Release\genCRN.exe -n3 -s1;1;1 crn_3_%%R.txt > actual3%%R_111.txt	
		..\TestGenCRN\bin\Release\TestGenCRN.exe -s 1 1 1 -e expected3%%R.txt -a actual3%%R_111.txt
	)
	
	echo --------------------------------------
)