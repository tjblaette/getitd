# getITD v1.5.8

getITD for FLT3-ITD-based MRD monitoring in AML
https://doi.org/10.1038/s41375-019-0483-z


## Author / Support
T.J. Luck, née Blätte


## Overview
getITD was developed for FLT3-ITD detection from Illumina sequencing data. Adaptations to data generated by similar assays, sequencing technologies or target regions should be possible. In fact, getITD was successfully used to analyze 454 data of different amplicons as presented at ASH 2018 (Rücker et al. "Prognostic Impact of Insertion Site in Acute Myeloid Leukemia (AML) with FLT3 Internal Tandem Duplication: Results from the Ratify Study (Alliance 10603)".).

The core getITD analysis program is contained within getitd.py.
In addition, two wrapper scripts, make_getitd_config.py and getitd_from_config_wrapper.py, which do not require use of the commandline are provided to facilitate getITD analysis.

Setup and usage of all programs is described in detail below.


## Requirements
python3
- pandas
- numpy
- biopython
- easygui
- matplotlib

The pandas, numpy and biopython modules are required for the core getITD program. The easygui and matplotlib modules are optional: easygui is only required by the make_getitd_config.py script, whereas matplotlib is required for coverage plotting (optionally activated using `-plot_coverage True`).


## Setup (Linux / Ubuntu & MacOS)
#### If applicable: Set proxy
In case you are working behind a firewall, set the appropriate proxy server (ask your IT administrator if you do not know which one this is) by opening the commandline and entering the following commands:
```
export https_proxy=your.proxy.server:your_port
export http_proxy=your.proxy.server:your_port
```
#### Python3 modules
To install the necessary python3 modules, open the commandline and enter:
```
pip3 install --user numpy pandas biopython easygui matplotlib
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
pip3 install --user pandas numpy biopython easygui matplotlib
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

## getITD core program
### Input
###### Required:
The only required inputs are a sample ID to name output files and folders and the actual FASTQ files containing the sequencing data to process. Supply these in the correct order (see also the examples below):

| Position | Description |
| --- | --- |
| 1st | *sample ID* which is used for naming  of output files and folders. This parameter is *required*. |
| 2nd | Path and filename of the FASTQ file containing forward reads to process, by convention the _R1_ FASTQ file. When processing single-end sequencing data, supply the single FASTQ file that was generated. This parameter is *required*. |
| 3rd | Path and filename of the FASTQ file containing reverse reads to process, by convention the _R2_ FASTQ file generated by paired-end sequencing. This parameter is *optional*. |


Note that FASTQ files may be provided as both compressed (gzipped) `*.fastq.gz` or uncompressed `*.fastq` files.
getITD will automatically detect file type and appropriately read them in; commands do not need to be changed.


###### Optional:
Many optional parameters are available to customize the analysis:

| Parameter | Description |
| --- | --- |
| `-reference X` | Path and filename of the reference to align to. This file should contain a single line of text, namely the WT sequence of the amplicon that was processed. getITD aligns all reads to this reference and detects insertions and ITDs relative to it. Default: _./anno/amplicon.txt_. |
| `-anno X` | Path and filename of the annotation file used to retrieve chromosomal, transcriptomic and proteomic coordinates and domains for each insertion and ITD identified by getITD. Columns should be tab-separated with header _amplicon_bp_, _region_, _chr_bp_, _transcript_bp_, _protein_as_. Default: _./anno/amplicon_kayser.tsv_. |
| `-forward_primer X` | Gene-specific sequence(s) of the forward primer(s) used to generate the sequenced amplicon. The sequence should be identical to the 5' end of supplied forward reads. When sequencing a pool of multiple amplicons, provide information on all primer pairs used. Separate individual sequences by space. If this is the last optional parameter used, pass also " -- " to signal that what follows are not more primer sequences. Default: _GCAATTTAGGTATGAAAGCCAGCTAC_. |
| `-reverse_primer X` | Gene-specific sequence(s) of the reverse primer(s) used. This sequence should be identical to the 5' end of supplied reverse reads. When sequencing a pool of multiple amplicons, provide information on all primer pairs used. Separate individual sequences by space. If this is the last optional parameter used, pass also " -- " to signal that what follows are not more primer sequences. Default: _CTTTCAGCATTTTGACGGCAACC_. |
| `-require_indel_free_primers X` | If True, discard reads containing insertions or deletions within the primer sequence, as these indicate low sequence fidelity. Note that, if set to _True_, this also filters reads not containing any primer sequence at all, so this must be set to _False_ in case primers have been trimmed. Note that adapters but not primers should be trimmed. Default: _True_. |
| `-forward_adapter X` | Sequencing adapter of the forward reads' primer as (potentially) present at the 5' end of the supplied forward reads, 5' of the gene-specific primer sequence. Default: _TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGA_. |
| `-reverse_adapter X` | Sequencing adapter of the reverse reads' primer as (potentially) present at the 5' end of the supplied reverse reads, 5' of the gene-specific primer sequence. Default: _GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGA_. |
| `-plot_coverage X` | If True, plot coverage distribution across the reference to a png file. Default: False. |
| `-technology X` | Sequencing technology used, options are _454_ and _Illumina_. Currently this has no effect for Illumina data. For 454 data, -infer_sense_from_alignment is set to True and -min_read_copies is set to 1, command line options passed to these parameters are then ignored. Default: _Illumina_. |
| `-infer_sense_from_alignment X` | If True, infer each read's sense from alignment by aligning both its actual sequence as well as its reverse-complement and keeping whichever achieves the higher mapping score to the reference. If False, infer each read's sense from the file of origin (R1 vs R2 FASTQ). Default: False. |
| `-nkern X` | The number of cores to use for parallel tasks. Set to 1 to disable parallelization. Default: _12_. |
| `-gap_open X` | Alignment cost of a gap opening. Default: _-36_ |
| `-gap_extend X` | Alignment cost of a gap extension. Default: _-0.5_ |
| `-match X` | Alignment cost of a base match. Default: _5_ |
| `-mismatch X` | Alignment cost of a base mismatch. Default: _-15_ |
| `-minscore_inserts X` | Fraction of the maximum possible alignment score required between inserts and respective WT tandems when deciding whether these are sufficiently similar to consider them ITDs and inserts supported by distinct reads when determining whether these are supporting the same mutation or not. Higher values will require sequences in both cases to be more similar. Default: _0.5_. |
| `-minscore_alignments X` | Fraction of the maximum possible alignment score required between each read and the reference. Reads that do not pass this filter are discarded to eliminate low fidelity sequences. Default: _0.4_. |
| `-min_bqs X` | Minimum average base quality score (BQS) required by each read. Reads that do not pass this filter are discarded to eliminate sequences likely to contain sequencing errors. Default: _30_. |
| `-min_insert_seq_length X` | Minimum number of sequenced insert bases required for getITD to call an insertion / ITD. Default: _6_. Increase to _9_ or higher to avoid false positives in lower quality samples. |
| `-min_read_length X` | Minimum read length required by each read. This filter is applied after trailing _N_ bases are trimmed at both the 5' and 3' end. Reads that do not pass this filter are discarded. Default: _100_. |
| `-min_read_copies X` | Minimum number of copies of each read required for processing. When set to _2_, all unique read sequences are discarded, based on the assumption that _true_ and _clinically relevant_ sequences will be present at least twice. Disable this filter, by setting it to _1_, when not all reads have the same length as it will otherwise incorrectly discard the majority of input reads. Default: _2_. |
| `-filter_ins_unique_reads X` | Minimum number of unique supporting reads required for an insertion / ITD to be considered _high confidence_. Disable this filter by setting it to _1_. Default: _2_. |
| `-filter_ins_total_reads X` | Minimum number of total supporting reads required for an insertion / ITD to be considered _high confidence_. Default: _1_ (disabled). |
| `-filter_ins_vaf X` | Minimum variant allele frequency (VAF), in percent, required for an insertion / ITD to be considered _high confidence_. Default: _0.006_. Set to _0.2_ for data multiplexed using single-unique barcodes. |




For a complete list of available commandline arguments, a short description and respective default values, change to the getitd folder and, on Linux & MacOS, run
```
python3 getitd.py --help
```
On Windows, run
```
py -3 getitd.py --help
```
This will print
```console
usage: getitd.py [-h] [-reference REFERENCE] [-anno ANNO]
                 [-forward_primer FORWARD_PRIMER [FORWARD_PRIMER ...]]
                 [-reverse_primer REVERSE_PRIMER [REVERSE_PRIMER ...]]
                 [-require_indel_free_primers REQUIRE_INDEL_FREE_PRIMERS]
                 [-forward_adapter FORWARD_ADAPTER]
                 [-reverse_adapter REVERSE_ADAPTER]
                 [-plot_coverage PLOT_COVERAGE] [-technology {Illumina,454}]
                 [-infer_sense_from_alignment INFER_SENSE_FROM_ALIGNMENT]
                 [-nkern NKERN] [-gap_open GAP_OPEN] [-gap_extend GAP_EXTEND]
                 [-match MATCH] [-mismatch MISMATCH]
                 [-max_trailing_bp MAX_TRAILING_BP]
                 [-minscore_inserts MINSCORE_INSERTS]
                 [-minscore_alignments MINSCORE_ALIGNMENTS] [-min_bqs MIN_BQS]
                 [-min_read_length MIN_READ_LENGTH]
                 [-min_read_copies MIN_READ_COPIES]
                 [-min_insert_seq_length MIN_INSERT_SEQ_LENGTH]
                 [-filter_ins_unique_reads FILTER_INS_UNIQUE_READS]
                 [-filter_ins_total_reads FILTER_INS_TOTAL_READS]
                 [-filter_ins_vaf FILTER_INS_VAF]
                 sampleID fastq1 [fastq2]

positional arguments:
  sampleID              sample ID used as output folder prefix (REQUIRED)
  fastq1                FASTQ file (optionally gzipped) of forward reads (REQUIRED)
  fastq2                FASTQ file (optionally gzipped) of reverse reads (optional)

optional arguments:
  -h, --help            show this help message and exit
  -reference REFERENCE  WT amplicon sequence as reference for read alignment
                        (default ./anno/amplicon.txt)
  -anno ANNO            WT amplicon sequence annotation (default
                        ./anno/amplicon_kayser.tsv)
  -forward_primer FORWARD_PRIMER [FORWARD_PRIMER ...]
                        Forward primer gene-specific sequence(s) as present at
                        the 5' end of supplied forward reads. Separate by
                        space when supplying more than one (default
                        GCAATTTAGGTATGAAAGCCAGCTAC)
  -reverse_primer REVERSE_PRIMER [REVERSE_PRIMER ...]
                        Reverse primer gene-specific sequence(s) as present at
                        the 5' end of supplied reverse reads. Separate by
                        space when supplying more than one (default
                        CTTTCAGCATTTTGACGGCAACC)
  -require_indel_free_primers REQUIRE_INDEL_FREE_PRIMERS
                        If True, discard i) reads containing insertions or
                        deletions within the primer sequence and ii) reads not
                        containing any primer sequence. Set to False if these
                        have been trimmed (default True)
  -forward_adapter FORWARD_ADAPTER
                        Sequencing adapter of the forward reads' primer as
                        (potentially) present at the 5' end of the supplied
                        forward reads, 5' of the gene-specific primer sequence
                        (default TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGA)
  -reverse_adapter REVERSE_ADAPTER
                        Sequencing adapter of the reverse reads' primer as
                        (potentially) present at the 5' end of the supplied
                        reverse reads, 5' of the gene-specific primer sequence
                        (default GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGA)
  -plot_coverage PLOT_COVERAGE
                        If True, plot read coverage across the reference to
                        'coverage.png' in the respective output folder
                        (default False)
  -technology {Illumina,454}
                        Sequencing technology used, options are '454' or
                        'Illumina' (default). '454' sets
                        -infer_sense_from_alignment to True and
                        -min_read_copies to 1, regardless of the respective
                        command line options used; 'Illumina' will instead use
                        these command line options or their respective
                        defaults.
  -infer_sense_from_alignment INFER_SENSE_FROM_ALIGNMENT
                        If True, infer each read's sense by aligning it as a
                        forward and reverse read and keeping the better
                        alignment (default False).
  -nkern NKERN          number of cores to use for parallel tasks (default 12)
  -gap_open GAP_OPEN    alignment cost of gap opening (default -36)
  -gap_extend GAP_EXTEND
                        alignment cost of gap extension (default -0.5)
  -match MATCH          alignment cost of base match (default 5)
  -mismatch MISMATCH    alignment cost of base mismatch (default -15)
  -max_trailing_bp MAX_TRAILING_BP
                        maximum number of aligned bp between the start / end
                        of an insertion and the start / end of the read to
                        consider the insertion 'trailing'. Trailing insertions
                        are not required to be in-frame and will be considered
                        ITDs even if the matching WT tandem is not directly
                        adjacent. Set this to 0 to disable (default 0).
  -minscore_inserts MINSCORE_INSERTS
                        fraction of max possible alignment score required for
                        ITD detection and insert collapsing (default 0.5)
  -minscore_alignments MINSCORE_ALIGNMENTS
                        fraction of max possible alignment score required for
                        a read to pass when aligning reads to amplicon
                        reference (default 0.4)
  -min_bqs MIN_BQS      minimum average base quality score (BQS) required by
                        each read (default 30)
  -min_read_length MIN_READ_LENGTH
                        minimum read length in bp required after N-trimming
                        (default 100)
  -min_read_copies MIN_READ_COPIES
                        minimum number of copies of each read required for
                        processing (1 to turn filter off, 2 (default) to
                        discard unique reads)
  -min_insert_seq_length MIN_INSERT_SEQ_LENGTH
                        minimum number of insert basepairs which must be
                        sequenced of each insert for it to be considered by
                        getITD. For non-trailing ITDs, this is the minimum
                        insert length; for trailing ITDs, it is the minimum
                        number of bp of a potentially longer ITD which have to
                        be sequenced (default 6).
  -filter_ins_unique_reads FILTER_INS_UNIQUE_READS
                        minimum number of unique reads required to support an
                        insertion for it to be considered 'high confidence'
                        (default 2)
  -filter_ins_total_reads FILTER_INS_TOTAL_READS
                        minimum number of total reads required to support an
                        insertion for it to be considered 'high confidence'
                        (default 1)
  -filter_ins_vaf FILTER_INS_VAF
                        minimum variant allele frequency (VAF) required for an
                        insertion to be considered 'high confidence' (default
                        0.006)

```

### Output
##### Files
For each sample, an output directory will be created in the current directory, named using the provided sample ID.
Inside, all of the generated output files reside:
- config.txt contains parameters, date and time of the analysis, as well as the getITD version that was used
- stats.txt contains the number of reads and insertions / ITDs processed and filtered at each step of the analysis
- coverage.txt contains the read coverage across the reference sequence
- coverage.png contains a plot of that same coverage data (optionally created using `-plot_coverage True`)
- out\_needle/ contains individual alignment files, needle\_\*.txt, of all the different reads processed
- itds\_collapsed-is-same\_is-similar\_is-close\_is-same-trailing\_hc.tsv contains filtered high-confidence (hc) ITDs, fully merged
- itds\_collapsed-is-same\_is-similar\_is-close\_is-same-trailing.tsv contains all ITDs, fully merged
- itds\_collapsed-is-same\_is-similar\_is-close.tsv contains all ITDs, having merged those of the same insertion length, nearby insertion sites (within one tandem length) and similar insert sequences
- itds\_collapsed-is-same\_is-similar.tsv contains all ITDs, having merged those of the same insertion length and site and similar insert sequences
- itds\_collapsed-is-same.tsv contains all ITDs, having merged those that share the same insertion length, site and sequence
- insertions\*.tsv files are analogous to ITD files but list all insertions, regardless of whether these are also ITDs or not

Note that `insertion*tsv` and `itd*tsv` files are only created if at least one insertion or ITD is found, respectively. Config and stats files should exist for all samples.

##### On stdout
At each filtering step, the number of reads and insertions passing the specified requirements are printed. This basically corresponds to the content of _stats.txt_.
Currently, also the computation time of various steps is printed, but presumably this will change in the future.

### Example
##### Using default parameters
To analyze the provided test data, change to the getitd folder and, on Linux / MacOS, run:
```
python3 getitd.py test test/test_R1.fastq test/test_R2.fastq
```
On Windows, run:
```
py -3 getitd.py test test\test_R1.fastq test\test_R2.fastq
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
Filtering 20 / 170 alignments with indels in primer bases
Total reads remaining for analysis: 3669 (73.38 %)
Calculating coverage took 0.11887545138597488 s

-- Looking for insertions & ITDs --
Collecting inserts took 0.020746199414134026 s
44 inserts >= 6 bp were found
0/44 insertions were part of adapters and filtered
Filtering inserts for adapter sequences took 2.495385706424713e-05 s
Annotating coverage took 0.005261320620775223 s
Collecting ITDs took 0.3815697953104973 s
44 ITDs were found

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

#### Changing optional parameters
To change, for example, the  number of cores used for the analysis to _2_, on Linux / MacOS instead run:
```
python3 getitd.py -nkern 2 test test/test_R1.fastq test/test_R2.fastq
```
On Windows, run:
```
py -3 getitd.py -nkern=2 test test\test_R1.fastq test\test_R2.fastq
```

#### Specifying primer sequences
Because an arbitrary number of primer sequences can be supplied following the `-forward_primer` and `-reverse_primer` parameters, an explicit signal is required to identify the end of the given list of primer sequences. This can be either any other optional argument which reads in only one word, such as `-nkern`, or the 'empty' argument ` -- `:
```
python3 getitd.py -forward_primer GCAATTTAGGTATGAAAGCCAGCTAC -nkern 2 test test/test_R1.fastq test/test_R2.fastq
```
or
```
python3 getitd.py -forward_primer GCAATTTAGGTATGAAAGCCAGCTAC -- test test/test_R1.fastq test/test_R2.fastq
```


### How it works
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
    - at least 40% of the maximum possible alignment score
    - gap-free alignment to the given primer sequence
6. Collect insertions within passing alignments, require
    - insert length of at least 6 bp
    - absence of ambiguous "N" bases within the actual insert sequence
    - in-frame insert (length divisible by 3) for inserts fully contained within a read
    - 3' or 5' trailing inserts not fully spanned by the sequenced reads are not required to be in-frame, since their exact length is not actually known
7. Collect ITDs from passing insertions, require
    - insert sequence can be realigned to the WT amplicon sequence, with at least 50% of the maximum possible alignment score
    - for 3' and 5' trailing inserts that insert sequence is not an adapter artefact
8. Merge total insertions and ITDs independently (considering reads to describe the same event and adding up supporting counts), require
    - that they are actually the same: insert length, site and sequence are identical
    - that they are similar: insert length and site are the same, sequences are similar
    - that they are close: insert length is the same, sequences are similar and sites are within one insert length of each other
    - close trailing inserts with similar sequences are considered the same regardless of their (estimated) lengths
9. Filter insertions and ITDs independently, require
    - at least two different supporting reads
    - a minimum variant allele frequency (VAF) of 0.006 %


## getITD wrappers - running getITD without the commandline
Two additional scripts are provided to facilitate getITD analysis by alleviating the need to use the commandline. The first, make_getitd_config.py, when double-clicked, creates a simple graphical interface and asks the user to select desired analysis parameters and annotation / reference files. These parameters are then saved to a config.txt file. When this config.txt is placed in a folder containing both getitd_from_config_wrapper.py and an arbitrary number of FASTQ files to analyze, double-clicking the wrapper script will initiate sequential analysis of all of the provided FASTQ files using parameters stored in config.txt. Though these wrapper scripts were originally intended primarily to facilitate analysis on Windows, they run also on Linux and Mac.

