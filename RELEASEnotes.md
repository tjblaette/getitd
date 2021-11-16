# getITD 1.5.14  2021-11-16

### Fix special cases of 3' trailing inserts & extend incomplete-wt-tandem.log
Special processing is required to properly recognize and annotate trailing
inserts. While this was previously done for 5' inserts, it had not yet been
applied to 3' trailing inserts. This update makes sure trailing inserts on
both ends of the read are handled in the same way.
Note that this would previously have affected the coverage estimate of
some 3' trailing inserts, especially those that were long and located in
the middle of the reference sequence. In rare cases, it could have caused
getITD to fail with an assertionError when calculting respective VAFs.

Inserts that map to a WT tandem which is only partially contained within
the read have previously been written to incomplete-wt-tandem.log with
their full read sequence only. I have now added read counts, to show the
number of times each sequence appeared in a given sample, and the insert
sequence itself, to allow comparison with results in the insertion*.tsv
output files. This should help users judge how relevant sequences with
these incompletely sequenced (potential) ITDs are.


# getITD 1.5.13  2021-04-30

### Clean up getITD output messages by rounding long decimals
Both computation time estimates and long fractions are now
rounded to make output on the console and in `stats.txt` more
legible.


# getITD 1.5.12  2021-04-23

### Clean up getITD output messages on incomplete WT tandems
getITD requires that an ITD's WT tandem is completely covered
by the read. If the insert maps beyond the read's boundaries, the
ITD is considered a likely false positive and not called. I'm not
sure it's a guaranteed false positive in all cases, for all protocols,
so I am keeping these reads saved to a separate file. What is new now
is that instead of printing an explanatory info message each time
such a read is found during the analysis, only one summarizing info
box is printed after the analysis, where such reads were present.
This should make the output more readable and concise.


# getITD 1.5.11  2021-04-16

### Print readable error message when reads are too short for analysis
By default, getITD trims off terminal N bases of all reads and then
filters by length to discard sequences too short for analysis. The
default cutoff is 100bp and shorter reads are discarded. While this
is reasonable for 150bp or longer reads, it will remove all reads
that start out being only 100bp long or shorter. Previously, getITD
failed with an incomprehensible numpy error. Now it prints an
informative error message and points users to the `-min_read_length`
parameter that they can adjust accordingly for shorter read protocols.


# getITD 1.5.10  2021-02-04

### Reintroduce explicit retraint that WT tandem is fully sequenced
getITD v1.5.4 removed the constraint that an ITD's WT has to be fully
sequenced. This filter does seem to be necessary, so it is reimplemented
with this release.


# getITD 1.5.9  2021-01-28

### Fix bug affecting samples with identical reads in R1 and R2 FASTQ files
There was previously a small bug in the way reads were handled which shared
identical sequences but originated from both R1 and R2 fastq files.
In rare cases, this could cause getITD to fail when calculating an
insert VAF.


# getITD 1.5.8  2021-01-28

### Add getITD version to config.txt output
To keep track of which version was used for which analysis, this info
is now printed to the `config.txt` output file, right at the top below
the date and time of the analysis.


# getITD 1.5.7  2021-01-28

### Add parameter to adjust an insertion/ITD's required minimum insert length
Previously, getITD called all insertions and ITDs where the insert was at least
6 bp long. Lower quality samples, however, can be prone to many spurious
false positive calls, which appear as variable 6-7 bp trailing inserts,
mapped to ITDs of varying lengths, but supported only at low VAFs by reverse
reads only. To allow filtering of such cases, an additional parameter
`-min_insert_seq_length` was introduced. The default value is `6`, to
recapitulate getITD's previous behaviour, but for low quality samples
it can now be set to `9` or higher.


# getITD 1.5.6  2021-01-22

### Output, and optionally plot, final read coverage across the reference
To better judge whether coverage across the reference was sufficient
to call ITDs, getITD will now generate 1-2 additional files per sample:

* `coverage.txt` in the getITD output folder will always be created and
list the coverage achieved at each position.
* if `-plot_coverage True` is given, coverage will additionally be plotted
to `coverage.png` in the same folder.

In both files, coverage is given separately for forward reads only, reverse
reads only, and all reads combined. Note that total coverage is not necessarily
the sum of forward and reverse read coverage, as paired-end reads of the same
physical DNA fragment are counted as one when they overlap.

Note also that i) positions are 0-based, so that the first WT base of the
reference is labelled `0` and the following bases are counted from there and ii)
all position coordinates for refer to the inter-bp coordinates implemented in
v1.5.2: Coordinate `0.5` therefore refers to the "space" between bases 0 and 1,
`1.5` refers to the space between 1 and 2, and so on. `-0.5` refers to the
space before the 5' start of the WT reference. Reads are counted as covering
a coordinate, if both adjacent bases are covered. Thus, if the inter-bp
coverage at coordinate `0.5` is `10`, this means that ten reads span across
bases 0 and 1. Similarly, read coverage at coordinate `-0.5` refers to reads
spanning across the 5' end of the reference with 5' trailing inserts.


# getITD 1.5.5  2021-01-19

### Merge trailing inserts/ITDs regardless of read sense
Previously, getITD required that trailing inserts could
only be merged when they were both trailing, shared a
similar insert sequence, had proximal WT tandems AND
were supported by reads of the same sense (either forward
or reverse) AND were trailing on the same end in those reads.

Depending on the assay design used, it is possible to obtain
overlapping reads of different sense, with trailing inserts
and different ends, describing the same insertion. To allow
merging in these cases, the restraint on sense and trailing
end has been removed for trailing insertions and ITDs.


# getITD 1.5.4  2021-01-19

### Remove explicit retraint that WT tandem is fully sequenced
getITD generally requires that at least one of the two tandems
of an ITD is fully sequenced. Previously, there was an explicit
check on this implemented, which filtered out all insertions
which did not fulfill it. This explicit check has now been removed,
as it did not handle all edge cases correctly, and standard ITDs
would not need it. It could filter some false positives though,
so may be reintroduced later.

To better evaluate this change, reads with incomplete WT tandems
are written to a separate log file in the getITD output folder,
named "incomplete-wt-tandem.log".


# getITD 1.5.3  2021-01-19

### Fix trailing ITD coordinate calculation
Previously, some 5' trailing ITDs did not have their coordinates
correctly calculated, which caused getITD to fail on some internal
checks.


# getITD 1.5.2  2021-01-19

### Fix coverage calculation
Previously, coverage calculation of 5' trailing insertions and ITDs
was not handled properly. While in most cases this should only lead
to slightly incorrect VAF estimates, in some extreme cases getITD
could fail as VAFs came out at > 100%.

This is now handled by using a completely new way of calculating
the coverage throughout the amplicon: Instead of counting the
number of sequencing reads aligned to each position, getITD counts
the number of reads spanning across each pair of bases. This makes
coverage calculation for insertions, which by definition have no
real coordinate within the WT reference, much more straight-forward
and avoids the above problem.


# getITD 1.5.1  2021-01-19

### Filter out ITDs with presumable insertion sites within the WT tandem
I found a sample where a short, 6bp trailing insertion generated a spurious
alignment during ITD detection, which overlapped the insertion site itself.
To avoid such alignments in the future, an additional filter step was added,
which should discard these (rare) false hits in the future.


# getITD 1.5.0  2020-11-24

### Change default `-minscore_alignments` from 0.5 to 0.4
The previous default value was too stringent for both 250 bp
and 300 bp reads: Some reads with long insertions were wrongly
filtered so that ITDs could be missed, especially when they
were long and located towards the middle of the read.

While this affected very few patients in our cohort, it is
recommended to rerun previously analyzed samples, especially
when there is an unusually large fraction of reads filtered at
this step of the analysis, as visible in getITD's output to the
terminal and `stats.txt` file. (That's how we noticed the problem.)


# getITD 1.4.1  2020-11-24

### Update help page to reflect optionality of gzipped FASTQs
This is only a minor fix to the output printed with `python getitd.py -h`
to be explicit about gzipped FASTQ files being optional: Uncompressed
FASTQ can still be used as input to getITD too.


# getITD 1.4.0  2020-11-24

### Allow (optionally) gzipped FASTQ files as input
Previously, all FASTQ files had to be uncompressed before they could
be analyzed by getITD. Thanks to @methylnick, getITD can now also
process compressed gzipped files directly.

Note that none of the commands change and it is still possible to
analyze raw FASTQ files. getITD will automatically detect whether
a gzipped `fastq.gz` file was used as input or, as before, an
uncompressed `fastq` file.


# getITD 1.3.0  2020-07-30

### Fix internal -minscore_alignments and -minscore_inserts mixup
The parameter `-minscore_alignments` is supposed to define the minimum
alignment score required between a read and the reference sequence,
whereas `-minscore_inserts` should define the minimum alignment
score required between two inserts or an insert and a potential WT
tandem when determining whether two insertions/itds should be merged
and whether an insertion qualifies as an ITD, respectively. (This is
also how it was/is described in the README and help pages of getITD.)
Previously, however, `minscore_alignments` had internally been used
instead of `-minscore_inserts` when filtering insert to WT tandem alignments.
This mixup is now fixed to match definitions in the README.

Note that this only affects getITD analyses where different, non-default
values were used for `-minscore_alignments` and `-minscore_inserts`, as
default values for both parameters are `0.5`, and thus the same.


# getITD 1.2.3  2020-06-02

### Fix getitd_from_config.py
The parameters `-infer_sense_from_alignment` and `-require_indel_free_primers`
were not properly read from the provided config file and ended up being
always set to `True`. This is now fixed, based on a previous fix already
implemented for getITD 1.0.2.


# getITD 1.2.2  2019-10-29

### Add recently added config params to make_getitd_config.py
The parameters `-infer_sense_from_alignment` and `-max_trailing_bp`
had been missing from the `make_getitd_config.py` script, causing
getITD to fail when run via `getitd_from_config_wrapper.py`.
Parameters are now included and the config wrappers are again
functional.


# getITD 1.2.1  2019-09-05

### Fix -infer_sense_from_alignment
A few things were not working previously when read sense was to be
determined for *two* FASTQ files. This is now fixed, with sense
truly only set once alignment is successful and based on that alignment,
so that with -infer_sense_from_alignment set to True, it no longer
matters in which order two (paired) FASTQ files are supplied or whether
a single merged file is supplied.


# getITD 1.2.0  2019-08-29

### Rename -filter_reads to -min_read_copies
The new name should be more intuitive and was even already
used in some parts of the program / documentation. Now, -min_read_copies
is used consistently throughout getITD.  Functionality is not changed.

### Convert all input primer and adapter sequences to upper case
Previously, only the reference sequence was input as uppercase but
adapter and primer sequences were taken exactly as provided.
Alignment of lower-case sequences would thus fail. To prevent this,
all sequences are converted to upper-case once read into getITD.

### Clarify primer sequence specification in README.md
Note that both `-forward_primer` and `-reverse_primer` take an
arbitrary number of primer sequences and therefore require an
explicit signal to stop recognizing everything that follows as
primer sequences. This can be either any other optional argument
which reads in only a fixed number of arguments or the empty
argument ` -- `.


# getITD 1.1.1  2019-08-29

### Fix 1.1.0: new option -infer_sense_from_alignment was set but not used
The new option is now used.


# getITD 1.1.0  2019-08-29

### Add / Change command line options -infer_sense_from_alignment & -technology
Previously, `-technology 454` caused the sense of each read to be inferred
from alignment, by aligning the actual sequence as well as the reverse-complement
of each read and keeping whichever achieved the higher mapping score relative
to the reference. `-technology Illumina` had the sense inferred from each read's
file of origin, depending on whether they came from the R1 or R2 FASTQ file.
This same functionality is now achieved by setting `-infer_sense_from_alignment True`,
which should be much more transparent.

The `-technology 454` option was simultaneously extended: It now overwrites the
`-infer_sense_from_alignment` option, setting it to True regardless of what is
otherwise passed to the getITD command line, and it sets `-filter_reads 1`,
which previously had to be done for 454 data when reads are not all of
the same length.

Thus, v1.0.2 `-technology 454` becomes `-infer_sense_from_alignment True` in v1.1.0
and v1.0.2 `-technology 454  -filter_reads 1` becomes `-technology 454`.


# getITD 1.0.2  2019-08-28

### Fix -require_indel_free_primers parameter
Previously, this paramter was set to 'True' regardless of the command
line argument given. It now correctly evaluates to `False` when
'False' or 'false' is given and to `True` when 'True' or 'true' is
provided. Any other case will raise an exception and abort the analysis.

# getITD 1.0.1  2019-07-29

### Fix unique sequence counting
Previously, read sequences were not counted correctly for fully
overlapping mates (i.e. reads with the same sequence but different
strands). This affected read statistics in `stats.txt` as well
as VAF and coverage estimates. This fix properly distributes
per-sequence counts by strand.

