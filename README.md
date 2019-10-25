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

## Citation
If you use this code or build upon it, please use the following (bibtex) citation:
```bibtex
@InProceedings{
	title = "Fast enumeration of non-isomorphic chemical reaction networks",
	author = "Carlo Spaccasassi and Boyan Yordanov and Andrew Phillips and Neil Dalchau",
	booktitle = "Computational Methods in Systems Biology (CMSB 2019)",
	year = "2019"
}
```

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
| -s#   | Assign special classes to species |

#### Examples
Count the number of non-isomorphic CRNs with 3 species and 3 reactions:

`./genCRN -n3 -q crn_3_3.txt`

Print all non-isomorphic CRNs with 3 species and 3 reactions on console:

`./genCRN -n3 -q crn_3_3.txt`

Count all non-isomorphic mass-conserving CRNs with 3 species and 3 reactions:

`./genCRN -n3 -m -q crn_3_3.txt`

Count all non-isomorphic, non-trivial and connected CRNs with 3 species and 3 reactions:

`./genCRN -n3 -t -c -q crn_3_3.txt`

Count all non-isomorphic, non-trivial and connected CRNs with 3 species and 3 reactions, where species A is not isomorphic to species B or C:

`./genCRN -n3 -t -c -q -'s1;2' crn_3_3.txt`

(the use of single quotes allows some terminals to escape the semicolon character.)

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

## Counts

The following tables show the counts for the number of *genuine* non-isomorphic CRNs (those that include all *n* species), the number of genuine non-isomorphic reversible CRNs, and the number of genuine non-isomorphic non-trivial CRNs. These counts now go beyond what is reported in our CMSB paper. 

#### Non-isomorphic CRNs

Species \ Reactions |  1  |  2  |  3  |  4  |  5  |  6  |  7  |
-------------------:| ---:| ---:| ---:| ---:| ---:| ---:| ---:|
**1** | 6 | 15 | 20 | 15 | 6 | 1 | 0 |
**2** | 10 | 210 | 2,024 | 13,740 | 71,338 | 297,114 | 1,018,264 |
**3** | 5 | 495 | 17,890 | 414,015 | 7,262,666 | 103,511,272 | 1,244,363,180 |
**4** | 1 | 451 | 47,323 | 2,900,934 | 128,328,834 | 4,518,901,463 | 133,379,120,523 |
**5** | 0 | 204 | 55,682 | 7,894,798 | 763,695,711 | 56,929,248,832 |
**6** | 0 | 54 | 35,678 | 10,704,289 | 2,069,783,947 |
**7** | 0 | 8 | 13,964 | 8,386,321 | 3,041,467,242 |
**8** | 0 | 1 | 3,594 | 4,182,295 | 2,715,774,734 |
**9** | 0 | 0 | 639 | 1,417,784 | 1,595,551,325 |
**10** | 0 | 0 | 83 | 618,885 | 653,346,685 |

#### Filtering reversible CRNs

Species \ Reactions |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |
-------------------:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:|
**1** | 3 | 6 | 7 | 4 | 4 | 1 | 1 | 0 |
**2** | 6 | 60 | 296 | 989 | 2,516 | 4,997 | 8,241 | 11,271 |
**3** | 3 | 138 | 4,788 | 26,988 | 230,595 | 1,589,808 | 9,161,056 | 45,107,712 |
**4** | 1 | 134 | 6,354 | 187,005 | 4,048,219 | 69,982,180 | 1,011,965,511 |
**5** | 0 | 65 | 7,677 | 513,036 | 24,186,053 | 888,323,405 |
**6** | 0 | 21 | 5,178 | 709,212 | 66,152,034 | 4,674,311,477 |
**7** | 0 | 4 | 2,188 | 572,058 | 98,576,689 |
**8** | 0 | 1 | 648 | 298,030 | 89,754,652 |

#### Filtering non-trivial CRNs (-t)

Species \ Reactions |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |
-------------------:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:|
**1** | 0 |  9 |  18 |     15 |         6 |           1 |              0 |             0 |
**2** | 0 | 19 | 304 |  5,016 |    41,500 |     221,728 |        871,330 |     2,700,277 |
**3** | 0 |  8 | 464 | 25,272 | 1,125,465 |  30,806,874 |    563,453,020 | 7,675,100,687 |
**4** | 0 |  2 | 223 | 28,052 | 3,279,132 | 321,921,288 | 20,669,624,467 |
**5** | 0 |  0 |  41 | 12,340 | 2,845,389 | 633,623,890 |
**6** | 0 |  0 |   5 |  2,606 | 1,127,294 |
**7** | 0 |  0 |   0 |    264 |   238,105 |
**8** | 0 |  0 |   0 |     17 |    28,191 |
**9** | 0 |  0 |   0 |      0 |     1,795 |
**10**| 0 |  0 |   0 |      0 |        60 |

## CRN binary encoding (-b)
CRNs are printed by default in LBS format, which is human readable but not very efficient space-wise. To enable compact CRNs storage, `genCRN` can encode CRN reactions in byte format using the `-b` flag. The encoding hashes a reaction to either a one byte integer, if the number of species is less than or equal to 4, or to two bytes. 


Let N be the total number of species, m be the maximum number of complexes, and C be a complex. The hash function to encode a reaction C1 -> C2 into an integer is the following:
```
hash(C_1 -> C_2) = hash(C_1) * m + hash(C_2)
hash(C) =  0                    if C = 0                    (naught)
	|  1+i                  if C = Ai                   (monomer)
        |  N+1+i                if C = 2Ai                  (monodimer)
	|  2N + 1 + hash(i,j)   if C = Ai + Aj and i < j    (heterodimer)
hash(i, j) = (N * i - i * (i+1)/ 2) + (j - i - 1)
```
(*Ai* stands for the *i*th species, where *i* is an index from  0 to N-1)

The reason why CRNs with up to 4 species can be represented in a single byte is the following. The total number of complexes in a bimolecular CRN with N species is (N+2)-choose-2, which in the case of 4 species is 15. A reaction can be considered as a pair of (not ) complexes, therefore there are 15 * 14 = 210 possible reactions in a 4 species CRN, which fits into a single byte (up to 255). 

## CRN class partitioning
In the above description, CRN enumeration assumes that all species are isomorphic to each other in a CRN. In some scenarios, this might be inappropriate, because some species might have significantly different properties than others. For example, in reaction-diffusion systems, species A might represent a molecule that diffuses much faster in a solution than species B and C. Reactions A + B -> 2C and B + C -> 2A would normally be considered isomorphic, but the faster diffusivity of A would makes the spatiotemporal dynamics of the former reaction very different from the latter one. Therefore, it can be desirable to consider those two CRNs as non-isomophic. Whereas, reactions A + B -> 2C and A + C -> 2B would still be considered isomorphic, since B and C don't play a special role.

To support this scenario, we can specify that species A belongs to one class, and that species B and C belong to another separate class. More generally, `genCRN` supports the partitioning of CRN species into N classes. Species are considered isomorphic if and only if they belong to the same class. Given a Complex-Species graph (see "Fast enumeration of non-isomorphic chemical reaction networks"), the enumeration of such CRNs can be encoded by adding to the CS graph a new layer of vertices that assigns a class to a set of species, in the same fashion that species nodes in a CS graph assign a species to complexes. The enumeration algorithm of classes is very similar to the enumeration algorithm of CS graphs.
