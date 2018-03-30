#!/bin/bash


DIR=$1
BQS=$2


# min, max, average (all hc ITDs)
for file in *${DIR}*${BQS}/*itd*hc.tsv; do grep -c "${DIR}" $file; done | sort -n | head -n 1
for file in *${DIR}*${BQS}/*itd*hc.tsv; do grep -c "${DIR}" $file; done | sort -n | tail -n 1
for file in *${DIR}*${BQS}/*itd*hc.tsv; do grep -c "${DIR}" $file; done | awk '{ total += $1 } END { print total/NR }'

# min, max, average (non-trailing hc ITDs only)
#for file in *${DIR}*${BQS}/*itd*hc.tsv; do grep -c 'False' $file; done | sort -n | head -n 1
#for file in *${DIR}*${BQS}/*itd*hc.tsv; do grep -c 'False' $file; done | sort -n | tail -n 1
#for file in *${DIR}*${BQS}/*itd*hc.tsv; do grep -c 'False' $file; done | awk '{ total += $1 } END { print total/NR }'


# correlate coverage with number of ITDs found
# make sure FASTQ files and analysis ordners are sorted accordingly!
echo -e "pe_reads_sequenced\titds_found" > coverage_to_itds.tsv
paste <(for file in *R1*fastq ; do sed -n '2~4p' $file | wc -l ; done) <(for file in *${DIR}*${BQS}/*itd*hc.tsv; do grep -c "${DIR}" $file; done)  >> coverage_to_itds.tsv
