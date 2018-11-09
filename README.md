# getITD
version  
link to future paper

## Author / Support
T.J. Blätte  
tamara.blaette@uni-ulm.de

## Requirements
python3
- pandas
- numpy
- biopython

## Setup (Linux / Ubuntu & MacOS)
To install the necessary python3 modules, open the commandline and enter:
```
pip3 install --user numpy pandas biopython
```
##### For convenience: Change to the folder containing the _getitd.py_ program and _anno_ subfolder
Commands provided in this manual assume that you are working from within the folder that contains the _getitd.py_ script. Your current folder is shown at the beginning of the commandline. To go up in the folder hierarchy, enter:
```
cd ..
```
To enter a certain subfolder, type:
```
cd EXAMPLE
```
(Replace _EXAMPLE_ with the name of the folder you are trying to enter.)

Repeat these two steps to reach the desired destination.


Alternatively, provide full or relative paths to all files, including the _getitd.py_ script, input FASTQ files and reference and annotation files, whenever you use getITD. 

## Setup (Windows)
#### Python3
If necessary, download and install the latest python3 from https://www.python.org/downloads/windows/. Be sure to check the box to add it to your PATH. Check that python works by opening the Windows commandline (cmd), typing
```
py -3
```
and pressing ENTER. If python3 was correctly installed and added to your PATH, this will open the python interpreter. Exit by typing 
```
exit()
```
and pressing ENTER.

#### Visual Studio C++ Build Tools
Download and install Visual Studio C++ Build tools from https://visualstudio.microsoft.com/visual-cpp-build-tools (click _Download Build Tools_, scroll down to _All Downloads_ and select _Build Tools for Visual Studio_ under _Tools for Visual Studio_). 

For installation, select the following options in addition to the defaults:
- "C++/CLI support"
- "VC++ 2015.3 v14.00 (v140) toolset for desktop"

#### If applicable: Set proxy
In case you are working behind a firewall, set the appropriate proxy server (ask your IT administrator if you do not know which one this is) by entering the following commands in the Windws commandline (cmd):
```
set https_proxy=your.proxy.server:your_port
set http_proxy=your.proxy.server:your_port
```

#### Install python3 modules:
Open the Windows commandline (cmd) and install the required python3 modules:
```
pip3 install --user pandas numpy biopython
```
If this command fails with a connection error such as the one shown below, you are working behind a firewall and first need to set the correct proxy server as described above under _Set proxy_):
```
Retrying (Retry(total=4, connect=None, read=None, redirect=None, status=None)) after connection broken by 'NewConnectionError('<pip._vendor.urllib3.connection.VerifiedHTTPSConnection object at 0x00000218168F2A90>: Failed to establish a new connection: [WinError 10061] No connection could be made because the target machine actively refused it')': /simple/pandas/
```
#### For convenience: Change to the folder containing the _getitd.py_ program and _anno_ subfolder
Commands provided in this manual assume that you are working from within the folder that contains the _getitd.py_ script. Your current folder is shown at the beginning of the commandline. To go up in the folder hierarchy, enter:
```
cd ..
```
To enter a certain subfolder, type:
```
cd EXAMPLE
```
(Replace _EXAMPLE_ with the name of the folder you are trying to enter.)

Repeat these two steps to reach the desired destination.


Alternatively, provide full or relative paths to all files, including the _getitd.py_ script, input FASTQ files and reference and annotation files, whenever you use getITD. 

## Input
###### Required:
- Fastq files to process
- Sample ID / Output folder prefix  

###### Optional:
For a complete list of available commandline arguments and their default values, change to the getitd folder and, on Linux & MacOS, run
```
python3 getitd.py --help
```
On Windows, run
```
py -3 getitd.py --help
```

## Output 
##### Files
For each sample, an output directory will be created in the current directory, named using the provided Sample ID / Output folder prefix.
Inside, all generated output files reside:
- config.tsv contains parameters, date and time of the analysis
- out\_needle/ contains individual alignment files, needle\_\*.txt, of all the different reads processed
- itds\_collapsed-is-same\_is-similar\_is-close\_is-same-trailing\_hc.tsv contains filtered high-confidence (hc) ITDs, fully merged
- itds\_collapsed-is-same\_is-similar\_is-close\_is-same-trailing.tsv contains all ITDs, fully merged
- itds\_collapsed-is-same\_is-similar\_is-close.tsv contains all ITDs, having merged those of the same insertion length, nearby insertion sites (within one tandem length) and similar insert sequences 
- itds\_collapsed-is-same\_is-similar.tsv contains all ITDs, having merged those of the same insertion length and site and similar insert sequences
- itds\_collapsed-is-same.tsv contains all ITDs, having merged those that share the same insertion length, site and sequence
- insertions\*.tsv files are analogous to ITD files but list all insertions, regardless of whether these are also ITDs or not

##### On stdout
At each filtering step, the number of reads and insertions passing the specified requirements are printed.   
Currently, also the computation time of various steps is printed, but presumably this will change in the future. 

## Example
To analyze the provided test data, change to the getitd folder and, on Linux / MacOS, run:
```
python3 getitd.py test/test_R1.fastq test/test_R2.fastq test
```
On Windows, run:
```
py -3 getitd.py test\test_R1.fastq test\test_R2.fastq test
```

Expected output:  
```console
==== PROCESSING SAMPLE test ====
-- Reading FASTQ files --
Reading FASTQ files took 0.15001007914543152 s
Number of total reads: 5000
Number of total reads remainging after N-trimming: 5000 (100.0 %)
Mean read length after N-trimming: 251.0
Number of total reads with mean BQS >= 30: 4544 (90.88 %)
Getting unique reads took 0.04079139232635498 s

Number of unique reads with mean BQS >= 30: 987
Number of unique reads with at least 2 copies: 170
Total reads remaining for analysis: 3727 (74.54 %)

-- Aligning to Reference --
Alignment took 1.0301645826548338 s
Filtering 0 / 170 low quality alignments with a score < 50.0 % of max
Filtering 22 / 170 alignments with indels in primer bases
Total reads remaining for analysis: 3654 (73.08 %)
Calculating coverage took 0.11887545138597488 s

-- Looking for insertions & ITDs --
Collecting inserts took 0.020746199414134026 s
43 insertions were found
Filtering inserts for adapter sequences took 2.495385706424713e-05 s
0/43 insertions were part of adapters and filtered
Annotating coverage took 0.005261320620775223 s
Collecting ITDs took 0.3815697953104973 s
43 ITDs were found

-- Merging results --
3 insertions remain after merging
2 insertions remain after merging
1 insertions remain after merging
1 insertions remain after merging
3 itds remain after merging
2 itds remain after merging
1 itds remain after merging
1 itds remain after merging
Merging took 0.5721261408179998 s

-- Filtering --
Filtered 0 / 1 insertions based on the vaf
Filtered 0 / 1 insertions based on the number of unique supporting reads
Filtered 0 / 1 insertions based on the number of total supporting reads
1 insertions remain after filtering!
Filtered 0 / 1 itds based on the vaf
Filtered 0 / 1 itds based on the number of unique supporting reads
Filtered 0 / 1 itds based on the number of total supporting reads
1 itds remain after filtering!
Filtering took 0.03359447047114372 s
```

## How it works
1. Read in FASTQ files
    - Read in forward reads and BQS (R1)
    - Read in and reverse complement reverse reads and BQS (R2)
    - Trim trailing ambiguous N bases on all reads
2. Filter reads for minimum average base quality score (BQS)
    - Discard low quality reads
3. Filter unique reads
    - Assume that true & clinically relevant sequences will be present at least twice and discard unique reads
4. Align each read to the WT amplicon reference sequence using Needlemann-Wunsch alignment algorithm
5. Filter alignments, require
    - at least 50 % of the maximum possible alignment score
    - gap-free alignment to the primer sequence (when these are given)
6. Collect insertions within passing alignments, require
    - insert length of at least 6 bp
    - absence of ambiguous "N" bases within the actual insert sequence
    - in-frame insert (length divisible by 3) for inserts fully contained within a read
    - 3' or 5' trailing inserts not fully spanned by the sequenced reads are not required to be in-frame, since their exact length is not actually known then
7. Collect ITDs from passing insertions, require
    - insert sequence can be realigned to the WT amplicon sequence, again with at least 50 % of the maximum possible alignment score
    - for 3' and 5' trailing inserts that insert sequence is not an adapter artefact
8. Merge total insertions and ITDs independently (considering reads to describe the same event and adding up supporting counts), require
    - that they are actually the same: insert length, site and sequence are identical
    - that they are similar: insert length and site are the same, sequences are similar
    - that they are close: insert length is the same, sequences are similar and sites are within one insert length of each other
    - close trailing inserts with similar sequences are considered the same regardless of their (estimated) lengths
9. Filter insertions and ITDs independently, require
    - at least two different supporting reads
    - a minimum variant allele frequency (VAF) of 0.006 %
