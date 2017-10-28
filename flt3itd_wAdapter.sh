#!/bin/bash
# September 20th, 2017
# Tamara Jacqueline Blätte

## This script takes as input paired FASTQ files of Laura's FLT3 PCR product
## It filters FASTQ files by average base quality and aligns passing reads to the amplicon sequence using needle (Needleman-Wunsch alignment)
## It outputs alignments of each read in txt format which are used downstream to detect FLT3-ITDs

#########################################################
# set up variables and folder to run the analysis in

R1=$1
R2=$2
SAMPLE=$3
MIN_BQS=${4:-35}
DIR="${SAMPLE}_minBQS_${MIN_BQS}"

mkdir $DIR
cd $DIR

ln -s ../${R1} $R1 
ln -s ../${R2} $R2 
ln -s /NGS/known_sites/hg19/flt3-itd_anno/ anno


#########################################################
# quality filter and sort|uniq FASTQ files


perl /NGS/prinseq/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $R1 -min_qual_mean $MIN_BQS
perl /NGS/prinseq/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $R2 -min_qual_mean $MIN_BQS

mkdir other_fastq
mv *fastq other_fastq

# for reads that passed QC filtering, sort and count unique reads, use sed to get nice COUNT\tREAD format
# REMOVE SINGLE UNIQUE READS -> REQUIRE THAT EACH INSERT IS SUPPORTED BY AT LEAST 2 READS  ----> final awk
UNIQ_COUNTS_ALL=$(cat <(sed -n '2~4p' other_fastq/*_R1_*_good_*.fastq) <(sed -n '2~4p' other_fastq/*_R2_*_good_*.fastq | rev | tr 'ATCGatcg' 'TAGCtagc') | sort | uniq -c | sed -e 's/^ *//' -e 's/^\([0-9]*\) /\1\t/' | awk '$1 > 1')

echo "$UNIQ_COUNTS_ALL" > all.reads_counted
echo "$UNIQ_COUNTS_ALL" | cut -f2 > all.reads
echo "$UNIQ_COUNTS_ALL" | cut -f1 > all.readCounts


#########################################################
# run needle on a subset of all high-quality reads


# rev complement reference and keep reads unchanged or rev-complement reads of one direction and use one reference for both -> the latter will put start and end coords of ITDs into the same coord space (as long as they are generated relative to the reference and not the reads!)
i=0; while read line; do i=$((i+1)); needle -asequence <(echo $line) -bsequence anno/amplicon_wAdapter.txt -gapopen 20 -gapextend 0.5 -outfile all_${i}.needle; done < all.reads

# clean up files by copying everything into one file
# make sure to sort files so that lines match those of all.reads, all.readCounts and all.reads_counted above (= original sort order based on index i in file name as assigned by while loop above)
cat $(ls all_*.needle | sort -V ) > all.needle 

# extract alignments and alignment IDs only and save these to separate files as well
grep '^[ ]\+[0-9]\|outfile' all.needle | sed -e 's/ \+[0-9]*//g' -e 's/#-outfile/\n/' | tail -n +2 > all.alignments

# extract statistics of needle alignment per read and save to separate files for each statistic and forward/reverse read respectively: alignment score, total length of gaps in the alignment, number of matching bases, number of non-gap bases (matches + mismatches)
grep 'Score' all.needle | cut -f2 -d':' | awk '{print $1}' > all.scores
grep 'Gaps' all.needle | awk '{print $3}' | cut -f1 -d'/' > all.gaps
grep 'Identity' all.needle | awk '{print $3}' | cut -f1 -d'/' > all.identity
grep 'Similarity' all.needle | awk '{print $3}' | cut -f1 -d'/' > all.similarity

# move read alignment files to separate folder 
mkdir out_needle
mv all_*.needle out_needle


#########################################################

# plot results
python3 /media/data/tabl/laura_mrd_flt3/get_stats2.py "all.alignments"


# collect unique reads' alignment for each category of detected insertion/itd
rm -f all_insertions.txt all_insertions_single.txt all_itds_exact.txt  all_itds_nonexact_fail.txt  all_itds_nonexact.txt
while read file ; do cat out_needle/$file >> all_insertions.txt ; done < flt3_insertions.txt 
while read file ; do cat out_needle/$file >> all_itds_exact.txt ; done < flt3_insertions_single_itd_exact.txt 
while read file ; do cat out_needle/$file >> all_itds_nonexact.txt ; done < flt3_insertions_single_itd_nonexact.txt 
while read file ; do cat out_needle/$file >> all_itds_nonexact_fail.txt ; done < flt3_insertions_single_itd_nonexact_fail.txt 


cd ..




