# genCRN: Fast enumeration of non-isomorphic chemical reaction networks

## Building genCRN
1. Build NAUTY from this forked location: https://github.com/CSpaccasassi/nauty
2. Build genCRN. *More instructions to follow*

## Binaries
To get started with genCRN without compiling the code, you can use pre-built executables (for Windows or Unix), which are in the dist folder.

## Usage
CRNs are generated in two steps. First, non-isomorphic input graphs for genCRN are produced using `geng` and `directg` in NAUTY. Then, these input graphs are passed to genCRN to produce non-isomorphic CRNs using the *Complex-Species Graph* encoding.

### Step 1: Generate Input Graphs
The simple program genInputGraphs.exe can be used to produce a series of calls these NAUTY utilities that enables the concatenation of all input graphs for a specified number of species and reactions. 

#### Example
Produce all input graphs of 3 species and 3 reactions: `./genInputGraphs.exe 3 3`

It produces as output:

```
./geng.exe 2 1 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 2 2 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 2 3 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 3 1 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 3 2 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 3 3 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 4 1 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 4 2 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 4 3 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 5 1 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 5 2 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 5 3 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 6 1 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 6 2 -d1 | ./directg.exe -e3 >> crn_3_3.txt
./geng.exe 6 3 -d1 | ./directg.exe -e3 >> crn_3_3.txt
```

By *piping* these results to the shell (e.g. `./genInputGraphs.exe 3 3 | bash`, the results of calling `geng` and `directg` are appended to the file crn_3_3.txt.

N.B. genCRN relies on the *piping* results from one program to another, and so if using Windows, we recommend using the Linux sub-system for Windows 10. 

### Step 2: Generate non-isomorphic CRNs

*More to follow*

| Flag  | Category |
| ----- | -------- |
| -help | Display help |
| -n#   | Number of species |
| -z    | Do not include naught in the reactions (e.g. 0 -> A) |
| -c    | Do not produce CRNs containing independent sub-CRNs (e.g. A->B, C->D) |
| -t    | Only produce "dynamically non-trivial" CRNs (see https://arxiv.org/abs/1705.10820) |
| -m    | Only produce mass-conserving CRNs (all species are under some conservation law) |
| -l    | Only produce CRNs with at least a conservation law |
| -x    | Only produce CRNs with no conservation laws |
| -q    | Only count the number of CRNs |
