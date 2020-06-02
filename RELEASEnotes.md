
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

