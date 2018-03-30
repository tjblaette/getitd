#!/bin/bash


# make known vaf file
#head make_known_vaf.tsv 
#76-dx   41.7
#76-rl   6.6
#76-cr   0
#150-dx  16.3
#...
 
#while read line; do SAMPLE=$(echo "$line" | cut -f1); VAF=$(echo "$line" | cut -f2); echo $VAF > ${SAMPLE}_known_vaf.txt; done < make_known_vaf.tsv



# concordance between GeneScan and NGS - ITDs found
# ITDs
echo 'sample    length  vaf     ref_coverage    counts  tandem2_start   insert' > summary_known_itds_missed.tsv
echo 'sample    length  vaf     ref_coverage    counts  tandem2_start   insert' > summary_known_itds_found.tsv
echo 'sample    length  vaf     ref_coverage    counts  start   insert' > summary_known_ins_missed.tsv
echo 'sample    length  vaf     ref_coverage    counts  start   insert' > summary_known_ins_found.tsv

grep --no-filename 'NA' */*itds*known.tsv >> summary_known_itds_missed.tsv
grep -v --no-filename 'NA' */*itds*known.tsv | grep -v 'ref_coverage' >> summary_known_itds_found.tsv

# INS
grep --no-filename 'NA' */*ins*known.tsv >> summary_known_ins_missed.tsv
grep -v --no-filename 'NA' */*ins*known.tsv | grep -v 'ref_coverage' >> summary_known_ins_found.tsv


# concordance between GeneScan and NGS - VAF estimate
# ITDs 
echo 'sample    length  vaf     vaf_genescan    vaf_each        tandem2_start   ref_coverage    counts_each' > summary_known_itds_vaf.tsv
cat */*itd*known_col*tsv | grep -v 'ref_coverage' >> summary_known_itds_vaf.tsv
echo 'sample    length  vaf     vaf_genescan    vaf_each        start   ref_coverage    counts_each' > summary_known_ins_vaf.tsv
cat */*ins*known_col*tsv | grep -v 'ref_coverage' >> summary_known_ins_vaf.tsv

# for plotting
grep 'vaf' */*itd*hc_known_*  | sort -u > summary_known_itds_vaf_genescan_correlation.tsv
cut -f1,3,4 */*itd*hc_known_* | grep -v 'vaf' | sort -k2,2nr >>  summary_known_itds_vaf_genescan_correlation.tsv

# recurrence of ITDs 
cut -f2,4 */flt*itd*ful*filtered.tsv | sort -r | uniq -c | sort -n > summary_recurrence_itds.tsv
cut -f2,3 */flt*ins*ful*filtered.tsv | sort -r | uniq -c | sort -n > summary_recurrence_ins.tsv



# total found
echo 'sample  length  tandem2_start   vaf     ref_coverage    counts  insert' > summary_all_itd_filtered.tsv
echo 'sample  length  tandem2_start   vaf     ref_coverage    counts  insert' > summary_all_itd.tsv
echo 'sample  length  start   vaf     ref_coverage    counts  insert' > summary_all_ins_filtered.tsv
echo 'sample  length  start   vaf     ref_coverage    counts  insert' > summary_all_ins.tsv

cat */flt3_itds_collapsed_full_filtered.tsv | sort | grep -v 'ref_coverage' >> summary_all_itd_filtered.tsv
cat */flt3_ins_collapsed_full.tsv | sort | grep -v 'ref_coverage' >> summary_all_itd.tsv
cat */flt3_ins_collapsed_full_filtered.tsv | sort | grep -v 'ref_coverage' >> summary_all_ins_filtered.tsv
cat */flt3_itds_collapsed_full.tsv | sort | grep -v 'ref_coverage' >> summary_all_ins.tsv


# clones not detected
echo 'sample    length  vaf     ref_coverage    counts  tandem2_start   insert' > summary_clones_missed.tsv
tail -n +2 */*itd*known.tsv  | grep 'NA' >> summary_clones_missed.tsv
