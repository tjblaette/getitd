
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

