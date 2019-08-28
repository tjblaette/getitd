
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

