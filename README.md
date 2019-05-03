# genCRN: Fast enumeration of non-isomorphic chemical reaction networks

## Building genCRN
1. git clone NAUTY from this forked location: https://github.com/CSpaccasassi/nauty
2. git clone genCRN from this repository in the same root folder as NAUTY

#### On Unix:
3. run `./configure --enable-wordsize=32` in the NAUTY folder as NAUTY
4. run `make all` in the NAUTY folder, which produces the `geng` and `directg` executables
5. run `make` in the genCRN folder, which produces the `genCRN` and `genInputGraphs` executables

#### On Windows:
3. Install Visual Studio 2019 (any edition will do, such as the Community Edition), making sure the "Desktop development for C++" workload is selected (this will install the MSVC compiler for C)
4. Start Visual Studio
5. Click "File" -> "Open" -> "Project/Solution..." and select the "genCRN.sln" file in the genCRN/vs folder
6. Click "Build" -> "Build solution"/. This creates `geng.exe`, `directg.exe`, `genCRN.exe` and `genInputGraphs` in the genCRN/vs/x64/Release folder

N.B. Make sure that the build configuration in Visual Studio is set on  "Release" rather than "Debug". When compiled in "Debug" mode, `genCRN` may experience a severe penalty in performance. More information here: https://docs.microsoft.com/en-us/visualstudio/debugger/how-to-set-debug-and-release-configurations?view=vs-2019.

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

Non-isomorphic CRNs are produced from the input file generated above. GenCRN reads each digraph and produces all non-isomorphic CRNs. The output format is LBS, unless the flag `-q` is used, to produce a total count only.

N.B. If the input graph is undirected, it is considered a reversible network (each undirected edge is turned into two opposite directed edges).

Below is a table showing all possible flags, and examples of using them.

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

#### Reversible CRNs
To enumerate reversible CRNs, it is sufficient to remove the invokation of `directg.exe` from the output of `genInputGraphs`. This can be achieved by another call to `sed`, such as `./genInputGraphs 3 3 | sed -e 's/ | \.\/directg -e3//g'`. The resulting input undirected graphs can then be passed to `genCRN`, which will automatically interpret them as digraphs where each undirected graph is replaced by two edges with opposite directions.

#### Parallelisation
CRN enumeration can be conducted in parallel, by splitting the input digraphs into different files, and passing them as input to different instances of genCRN. This can be achieved easily in Bash (and Git Bash in Windows).

For example, suppose that we wanted to enumerate all non-isomorphic CRNs with 5 species and 6 reactions on 65 cores. 
We first create the input digraphs file using `genInputGraphs` as above. Then the following command splits the input file into 65 files, distributing each digraph in round-robin fashion:

`split -d --additional-suffix=.txt --number=r/65 crn_5_6.txt crn_5_6_part`

The command splits the input file `crn_5_6.txt` into 65 files called `crn_5_6_part00`, `crn_5_6_part01`, ..., `crn_5_6_part64`.

Next, we can run 65 instances of `genCRN` as follows:
`for i in `ls crn_5_6_part*`; do (time ../GenCRN.exe -q -n5 $i > result_$i) & done`

Each instance of `genCRN` creates a file called `result_crn_5_6_i`, where `i` ranges from 0 to 64, and populates it with CRN counts when the enumeration is finished. To check if the enumeration is still running, a useful command is `ps aux | grep genCRN`.
