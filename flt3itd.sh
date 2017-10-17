#!/bin/bash
# September 20th, 2017
# Tamara Jacqueline Blätte

## This script takes as input paired FASTQ files of Laura's FLT3 PCR product
## It filters FASTQ files stringently by average base quality scores and takes 100000 from these high quality reads and aligns these to the amplicon sequence using needle
## It outputs alignments of each read in txt format that are used downstream to detect FLT3-ITDs

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


#########################################################
# quality filter and sort|uniq FASTQ files


perl /NGS/prinseq/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $R1 -min_qual_mean $MIN_BQS
perl /NGS/prinseq/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $R2 -min_qual_mean $MIN_BQS

mkdir other_fastq
mv *fastq other_fastq

# for reads that passed QC filtering, sort and count unique reads, use sed to get nice COUNT\tREAD format
# REMOVE SINGLE UNIQUE READS -> REQUIRE THAT EACH INSERT IS SUPPORTED BY AT LEAST 2 READS  ----> final awk
#UNIQ_COUNTS_R1=$(sed -n '2~4p' other_fastq/*_R1_*_good_*.fastq | sort | uniq -c | sed -e 's/^ *//' -e 's/^\([0-9]*\) /\1\t/' | awk '$1 > 1')
#UNIQ_COUNTS_R2=$(sed -n '2~4p' other_fastq/*_R2_*_good_*.fastq | sort | uniq -c | sed -e 's/^ *//' -e 's/^\([0-9]*\) /\1\t/' | awk '$1 > 1')
UNIQ_COUNTS_ALL=$(cat <(sed -n '2~4p' other_fastq/*_R1_*_good_*.fastq) <(sed -n '2~4p' other_fastq/*_R2_*_good_*.fastq | rev | tr 'ATCGatcg' 'TAGCtagc') | sort | uniq -c | sed -e 's/^ *//' -e 's/^\([0-9]*\) /\1\t/' | awk '$1 > 1')
#UNIQ_COUNTS_ALL=$(sed -n '2~4p' other_fastq/*_R?_*_good_*.fastq | sort | uniq -c | sed -e 's/^ *//' -e 's/^\([0-9]*\) /\1\t/' | awk '$1 > 1')

#echo "$UNIQ_COUNTS_R1" | cut -f2 > r1.reads
#echo "$UNIQ_COUNTS_R2" | cut -f2 > r2.reads
echo "$UNIQ_COUNTS_ALL" | cut -f2 > all.reads

#echo "$UNIQ_COUNTS_R1" | cut -f1 > r1.readCounts
#echo "$UNIQ_COUNTS_R2" | cut -f1 > r2.readCounts
echo "$UNIQ_COUNTS_ALL" | cut -f1 > all.readCounts

#echo "$UNIQ_COUNTS_R1" > r1.reads_counted
#echo "$UNIQ_COUNTS_R2" > r2.reads_counted
echo "$UNIQ_COUNTS_ALL" > all.reads_counted


#cat -n r1.readCounts | sort | awk '{print $2}' > r1.readCounts_needleOrder
#cat -n r2.readCounts | sort | awk '{print $2}' > r2.readCounts_needleOrder
cat -n all.readCounts | sort | awk '{print $2}' > all.readCounts_needleOrder


ln -s /NGS/known_sites/hg19/flt3-itd_anno/ anno


#########################################################
# run needle on a subset of all high-quality reads


# run needle with adapter sequences joined to amplicon reference ends -> avoid ambiguous gap lengths of reads +/- a few adapter bases
i=0; while read line; do i=$((i+1)); needle -asequence <(echo $line) -bsequence anno/amplicon_wAdapter.txt -gapopen 20 -gapextend 0.5 -outfile all_${i}.needle; done < all.reads
# rev complement reference and keep reads unchanged or rev-complement reads of one direction and use one reference for both -> the latter will put start and end coords of ITDs into the same coord space (as long as they are generated relative to the reference and not the reads!)
#echo 'STARTING ALIGNMENT OF R1'
#i=0; while read line; do i=$((i+1)); needle -asequence <(echo $line) -bsequence anno/amplicon_wAdapter.txt -gapopen 20 -gapextend 0.5 -outfile r1_${i}.needle; done < r1.reads
#echo 'STARTING ALIGNMENT OF R2'
#i=0; while read line; do i=$((i+1)); needle -asequence <(echo $line | rev | tr 'ATCGatcg' 'TAGCtagc') -bsequence anno/amplicon_wAdapter.txt -gapopen 20 -gapextend 0.5 -outfile r2_${i}.needle; done < r2.reads
echo 'DONE'



# clean up files by copying everything into one file for r1 and r2 respectively
#cat r1_*.needle > r1.needle
#cat r2_*.needle > r2.needle
cat all_*.needle > all.needle

# extract alignments and alignment IDs only and save these to separate files as well
#grep '^[ ]\+[0-9]\|outfile' r1.needle | sed -e 's/ \+[0-9]*//g' -e 's/#-outfile/\n/' | tail -n +2 > r1.alignments
#grep '^[ ]\+[0-9]\|outfile' r2.needle | sed -e 's/ \+[0-9]*//g' -e 's/#-outfile/\n/' | tail -n +2 > r2.alignments
grep '^[ ]\+[0-9]\|outfile' all.needle | sed -e 's/ \+[0-9]*//g' -e 's/#-outfile/\n/' | tail -n +2 > all.alignments
#grep '^[ ]\+[0-9]' r2.needle | sed -e 's/ 1 -/\n 1 -/' -e 's/[ 0-9]*//g' > r2.alignments



# extract statistics of needle alignment per read and save to separate files for each statistic and forward/reverse read respectively
# alignment score
#grep 'Score' r1.needle | cut -f2 -d':' | awk '{print $1}' > r1.scores
#grep 'Score' r2.needle | cut -f2 -d':' | awk '{print $1}' > r2.scores
grep 'Score' all.needle | cut -f2 -d':' | awk '{print $1}' > all.scores

# total length of gaps in the alignment
#grep 'Gaps' r1.needle | awk '{print $3}' | cut -f1 -d'/' > r1.gaps
#grep 'Gaps' r2.needle | awk '{print $3}' | cut -f1 -d'/' > r2.gaps
grep 'Gaps' all.needle | awk '{print $3}' | cut -f1 -d'/' > all.gaps

# number of matching bases
#grep 'Identity' r1.needle | awk '{print $3}' | cut -f1 -d'/' > r1.identity
#grep 'Identity' r2.needle | awk '{print $3}' | cut -f1 -d'/' > r2.identity
grep 'Identity' all.needle | awk '{print $3}' | cut -f1 -d'/' > all.identity

# number of non-gap bases (matches + mismatches)
#grep 'Similarity' r1.needle | awk '{print $3}' | cut -f1 -d'/' > r1.similarity
#grep 'Similarity' r2.needle | awk '{print $3}' | cut -f1 -d'/' > r2.similarity
grep 'Similarity' all.needle | awk '{print $3}' | cut -f1 -d'/' > all.similarity

# move read alignment files to separate folder 
mkdir out_needle_100k_wAdapter
#mv r1_*.needle out_needle_100k_wAdapter
#mv r2_*.needle out_needle_100k_wAdapter
mv all_*.needle out_needle_100k_wAdapter


#########################################################
# run Rscript for plotting -> create histogram with read counts on y and number of unmapped bases per read on x

#Rscript /media/data/tabl/laura_mrd_flt3/plotGaps_wAdapters.R $(echo $R1 | cut -f1-2 -d'_')
python3 /media/data/tabl/laura_mrd_flt3/get_stats2.py "all.alignments"


# collect unique reads' alignment for each category of detected insertion/itd
rm -f all_insertions.txt all_insertions_single.txt all_itds_exact.txt  all_itds_nonexact_fail.txt  all_itds_nonexact.txt
while read file ; do cat out_needle_100k_wAdapter/$file >> all_insertions.txt ; done < flt3_insertions.txt 
while read file ; do cat out_needle_100k_wAdapter/$file >> all_insertions_single.txt ; done < flt3_insertions_single.txt 
while read file ; do cat out_needle_100k_wAdapter/$file >> all_itds_exact.txt ; done < flt3_insertions_single_itd_exact.txt 
while read file ; do cat out_needle_100k_wAdapter/$file >> all_itds_nonexact.txt ; done < flt3_insertions_single_itd_nonexact.txt 
while read file ; do cat out_needle_100k_wAdapter/$file >> all_itds_nonexact_fail.txt ; done < flt3_insertions_single_itd_nonexact_fail.txt 




cd ..




