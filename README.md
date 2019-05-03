# genCRN: Fast enumeration of non-isomorphic chemical reaction networks

## Building genCRN
On Unix:
1. git clone NAUTY from this forked location: https://github.com/CSpaccasassi/nauty
2. git clone genCRN from this repository in the same root folder as NAUTY
3. run `./configure --enable-wordsize=32` in the NAUTY folder as NAUTY
4. run `make all` in the NAUTY folder, which produces the `geng` and `directg` executables
5. run `make` in the genCRN folder, which produces the `genCRN` and `genInputGraphs` executables

On Windows:
1. git clone NAUTY from this forked location: https://github.com/CSpaccasassi/nauty
2. git clone genCRN from this repository in the same root folder
3. Install Visual Studio 2019 (any edition will do, such as the Community Edition), making sure the "Desktop development for C++" workload is selected (this will install the MSVC compiler for C)
4. Start Visual Studio
5. Click "File" -> "Open" -> "Project/Solution..." and select the "genCRN.sln" file in the genCRN/vs folder
6. Click "Build" -> "Build solution" to create `genCRN.exe` and `genInputGraphs` in the genCRN/vs/Release folder

*More instructions to follow to compile NAUTY in Windows; pre-compiled binaries are available in the genCRN/dist folder*

## Binaries
To get started with genCRN without compiling the code, you can use pre-built executables (for Windows or Unix), which are in the dist folder.

## Usage
CRNs are generated in two steps. First, non-isomorphic input graphs for genCRN are produced using `geng` and `directg` in NAUTY. Then, these input graphs are passed to genCRN to produce non-isomorphic CRNs using the *Complex-Species Graph* encoding.

### Step 1: Generate Input Graphs
The simple program `genInputGraphs.exe` can be used to produce a series of calls to these NAUTY utilities, enabling the concatenation of all input graphs for a specified number of species and reactions. 

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

N.B. genCRN relies on the *piping* results from one program to another. On Windows, we recommend using Git Bash. On Unix, it is possible to trim out the `.exe` extensions with the command ``./genInputGraphs.exe 3 3 | sed -e 's/.exe//g' | bash``

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

#### Examples
Count the number of non-isomorphic CRNs with 3 species and 3 reactions:

`./genCRN -n3 -q crn_3_3.txt`

Print all non-isomorphic CRNs with 3 species and 3 reactions on console:

`./genCRN -n3 -q crn_3_3.txt`

Count all non-isomorphic mass-conserving CRNs with 3 species and 3 reactions:

`./genCRN -n3 -m -q crn_3_3.txt`

Count all non-isomorphic, non-trivial and connected CRNs with 3 species and 3 reactions:

`./genCRN -n3 -t -c -q crn_3_3.txt`

*More to follow*
