import Bio.pairwise2 as bio
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')  # required to use matplotlib without X (via ssh + screen)
import matplotlib.pyplot as plt
import datetime
import collections
import os
import multiprocessing
import argparse
import decimal as dc # required for accuracy of VAF (otherwise asserts fail for sum(counts_each) = counts due to floating point inaccuracies)
dc.getcontext().prec = 5 # number of digits to round decimals to


#######################################
## INITIALIZE VARIABLES

# prevent neg nkern/minBQS?
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
parser.add_argument("fastq1", help="FASTQ file of forward reads (REQUIRED)")
parser.add_argument("fastq2", help="FASTQ file of reverse reads (REQUIRED)")
parser.add_argument("sampleID", help="sample ID used as output folder prefix (REQUIRED)")
parser.add_argument("minBQS", help="minimum average base quality score (BQS) required by each read (default 30)", type=int, default=30, nargs='?')
parser.add_argument("-reference", help="WT amplicon sequence as reference for read alignment (default /NGS/known_sites/hg19/flt3-itd_anno/amplicon.txt)", default="/NGS/known_sites/hg19/flt3-itd_anno/amplicon.txt", type=str)
parser.add_argument("-anno", help="WT amplicon sequence annotation (default /NGS/known_sites/hg19/flt3-itd_anno/amplicon_kayser.tsv)", default="/NGS/known_sites/hg19/flt3-itd_anno/amplicon_kayser.tsv", type=str)
parser.add_argument('-nkern', help="number of cores to use for parallel tasks (default 14)", default="14", type=int)
parser.add_argument('-gap_open', help="alignment cost of gap opening (default -20)", default="-20", type=int)
parser.add_argument('-gap_extend', help="alignment cost of gap extension (default -0.5)", default="-0.5", type=float)
parser.add_argument('-match', help="alignment cost of base match (default 5)", default="5", type=int)
parser.add_argument('-mismatch', help="alignment cost of base mismatch (default -10)", default="-10", type=int)
parser.add_argument('-minscore_inserts', help="fraction of max possible alignment score required for ITD detection and insert collapsing (default 0.5)", default="0.5", type=float)
parser.add_argument('-minscore_alignments', help="fraction of max possible alignment score required for a read to pass when aligning reads to amplicon reference (default 0.5)", default="0.5", type=float)
parser.add_argument('-known_length', help="file with expected ITD length, one on each line (optional)")
group.add_argument('-known_vaf', help="file with total expected ITD VAF of all clones (optional)")
group.add_argument('-known_ar', help="file with total expected ITD allele ratio of all clones vs WT (optional)")
parser.add_argument('-filter_reads', help="minimum number of copies of each read required for processing (1 to turn filter off, 2 (default) to discard unique reads)", default="2", type=int)
parser.add_argument('-filter_ins_unique_reads', help="minimum number of unique reads required to support an insertion for it to be considered 'high confidence' (default 2)", default="2", type=int)
parser.add_argument('-filter_ins_total_reads', help="minimum number of total reads required to support an insertion for it to be considered 'high confidence' (default 1)", default="1", type=int)
parser.add_argument('-filter_ins_vaf', help="minimum variant allele frequency (VAF) required for an insertion to be considered 'high confidence' (default 0.001)", default="0.001", type=float)
cmd_args = parser.parse_args()

R1 = cmd_args.fastq1
R2 = cmd_args.fastq2
SAMPLE = cmd_args.sampleID
MIN_BQS = cmd_args.minBQS
REF = cmd_args.reference
ANNO = cmd_args.anno
NKERN = cmd_args.nkern
KNOWN_LENGTH_FILE = cmd_args.known_length
KNOWN_VAF_FILE = cmd_args.known_vaf
KNOWN_AR_FILE = cmd_args.known_ar
OUT_DIR = '_'.join([SAMPLE,'minBQS', str(MIN_BQS)])

COST_MATCH = cmd_args.match
COST_MISMATCH = cmd_args.mismatch
COST_GAPOPEN = cmd_args.gap_open
COST_GAPEXTEND = cmd_args.gap_extend
MIN_SCORE = cmd_args.minscore_inserts
MIN_SCORE_ALIGNMENTS = cmd_args.minscore_alignments

MIN_READ_COPIES = cmd_args.filter_reads
MIN_TOTAL_READS = cmd_args.filter_ins_total_reads
MIN_UNIQUE_READS = cmd_args.filter_ins_unique_reads
MIN_VAF = cmd_args.filter_ins_vaf


#######################################
## HELPER FUNCTIONS FOR PART 1 - PREPROCESSING

# read in txt file with known ITD length and VAF 
def read_known(filename, dtype=str):
    with open(filename) as f:
        return [dtype(x) for x in f.read().splitlines()]


# extract ITDs of known length from df  --> check that these are ints! (when best to fail if that's not the case?)
def get_known(df,known_length,sample=SAMPLE):
    df_found = df.ix[[x in known_length for x in df["length"]]]
    #
    # fill in available data on known ITDs/inserts that were missed (and not present in df)
    missed = [x for x in known_length if x not in list(df_found["length"])]
    df_missed = pd.DataFrame({"length": missed, "sample": [sample for x in range(len(missed))], "vaf": [0 for x in range(len(missed))], "counts": [0 for x in range(len(missed))]})
    #
    # concatenate known_found and known_missed
    df_known = pd.concat([df_found, df_missed])
    df_known[["length","counts"]] = df_known[["length","counts"]].astype(int)
    return df_known


# convert ITD allele ratio (AR) to variant allele frequency (VAF)
def ar_to_vaf(ar):
    return ar/(ar+1)


# read in FASTQ files, return reads and bqs
# --> possibly add saving to FASTQ file for all lines 
def read_fastq(filename):
    reads_and_bqs = []
    with open(filename,'r') as f:
        line = f.readline()
        while line:
            read_id = line
            read_seq = f.readline()
            read_desc = f.readline()
            read_bqs = f.readline()
            reads_and_bqs.append((read_seq.rstrip('\n'), read_bqs.rstrip('\n')))
            # read next line to see if there is one
            line = f.readline()
    return (reads_and_bqs)  


# filter reads based on average BQS 
def filter_bqs(args): 
    read, bqs = args[0]
    min_bqs = args[1]
    if sum([ord(x) -33 for x in bqs])/len(bqs) >= min_bqs:
        return read


# reverse complement and return a given read
def reverse_complement(read):
    return read.translate(str.maketrans('ATCGatcg','TAGCtagc'))[::-1]


# trim ambinguous N bases at reads' ends
def trim_n(read_and_bqs):
    read,bqs = read_and_bqs
    base_is_n = [x == 'n' or x == 'N' for x in read]
    n_start, n_end = 0,0
    while base_is_n.pop():
        n_end = n_end +1
    base_is_n.reverse()
    while base_is_n.pop():
        n_start = n_start +1
    return (read[n_start:len(read) - n_end], bqs[n_start:len(bqs) - n_end])


# read in wt reference for alignment
def get_reference(filename):
    ref = None
    with open(filename, 'r') as f:
        ref = f.read()
    ref = ref.splitlines()
    assert len(ref) == 1
    return ref[0]


# read in wt reference annotation (genomic / transcript / protein coordinates and exon/intron information)
def get_annotation(filename):
    return pd.read_csv(filename, sep='\t')


# actually annotate DataFrame
def annotate(df_data, df_anno, start_col):
    return pd.merge(df_data, df_anno, how='left', left_on=[start_col], right_on=['amplicon_bp']).drop('amplicon_bp', axis=1)


# get min score required to pass alignment score filter
def get_min_score(seq1, seq2, match_score=COST_MATCH, min_max_score_fraction=MIN_SCORE):
    return min(len(seq1),len(seq2)) * match_score * min_max_score_fraction


# callback function to calculate match and mismatch score for realignment of insert to its read
# insert is masked by 'Z' -> return max penalty (min score of -Inf) to prohibit realignment of insert to itself
def get_alignment_score(char1,char2):
  assert not (char1 == 'Z' and char2 == 'Z') # only ever one of the sequences chars are taken from should contain masking letter 'Z', i.e. the read sequence but not the ref
  if char1 == char2:
    return COST_MATCH
  elif char1 == 'Z' or char2 == 'Z':
    return -np.inf
  else:
    return COST_MISMATCH

# align two sequences
# args must be a tuple because multiprocessing.pool can only pass one argument to the parallelized function!
def align(args):
    seq1, seq2 = args
    min_score = get_min_score(seq1, seq2, min_max_score_fraction=MIN_SCORE_ALIGNMENTS)
    alignments = bio.align.globalcs(seq1, seq2, get_alignment_score, COST_GAPOPEN, COST_GAPEXTEND, penalize_end_gaps=False, one_alignment_only=True) # one alignment only until multiple ones are handled
    alignment_score = None
    if alignments != []:
        alignment_score = alignments[0][2]
        if alignment_score >= min_score:
            return alignments[0]
    return []


def connect_bases(char1, char2):
    if char1 == '-' or char2 == '-':
        return ' '
    if char1 == char2:
        return '|'
    return '.'


def connect_alignment(alignment):
    seq1, seq2 = alignment[0:2]
    return ''.join([connect_bases(char1,char2) for char1,char2 in zip(seq1,seq2)])


# count number of digits to align all lines in pretty alignment printout (need to know how many spaces to insert)
def get_number_of_digits(number):
    if number == 0: 
        return 1
    return int(np.log10(number)) +1


def print_alignment_connection(connection, pre_width,f):
    f.write(' ' * (pre_width +2))
    f.write(connection)
    f.write('\n')


def print_alignment_seq(seq, seq_coord, pre_width, post_width,f):
    f.write(' ' * (pre_width - get_number_of_digits(seq_coord) +1))
    f.write(str(seq_coord) + ' ')
    f.write(seq)
    seq_coord = seq_coord + len(seq) - seq.count('-') -1
    f.write(' ' * (post_width - get_number_of_digits(seq_coord)))
    f.write(str(seq_coord) + '\n')
    return seq_coord +1



# print pretty alignment, format inspired by EMBOSS needle output
def print_alignment(alignment, i, out_dir, command='bio.align.globalcs', command_seq='seq1', command_ref='seq2'):
    filename = 'needle_{}.txt'.format(i)
    seq, ref, score = alignment[0:3]
    al = connect_alignment(alignment)
    command_score_function = "get_alignment_score"
    cost_gapopen = COST_GAPOPEN
    cost_gapextend = COST_GAPEXTEND
    width = 50
    pre_width = 20
    post_width = 7
    score_width = 15
    #
    with open(os.path.join(out_dir,filename), 'w') as f:
        f.write('########################################\n')
        f.write('# Program: Biopython\n')
        f.write('# Rundate: {}\n'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%d")))
        f.write('# Commandline: {}(\n'.format(command))
        f.write('#    {},\n'.format(command_seq))
        f.write('#    {},\n'.format(command_ref))
        f.write('#    {},\n'.format(command_score_function))
        f.write('#    {},\n'.format(cost_gapopen))
        f.write('#    {})\n'.format(cost_gapextend))
        f.write('# Align_format: srspair\n')
        f.write('# Report_file: {}\n'.format(filename))
        f.write('########################################\n')
        f.write('\n')
        f.write('#=======================================\n')
        f.write('#\n')
        f.write('# Aligned_sequences: 2\n')
        f.write('# Sample: {}\n'.format(''.join([x for x in seq if x != '-'])))
        f.write('# Reference: {}\n'.format(''.join([x for x in ref if x != '-']).lower()))
        f.write('# Matrix: EDNAFULL\n')
        f.write('# Gap_penalty: {}\n'.format(cost_gapopen))
        f.write('# Extend_penalty: {}\n'.format(cost_gapextend))
        f.write('#\n')
        f.write('# Length: {}\n'.format(len(seq)))
        identity = '{}/{} ({}%)\n'.format(al.count('|'), len(seq), round(al.count('|')/len(seq) * 100,1))
        similarity = '{}/{} ({}%)\n'.format(al.count('|') + al.count('.'), len(seq), round((al.count('|') + al.count('.'))/len(seq) * 100,1))
        gaps = '{}/{} ({}%)\n'.format(len(seq) - al.count('|') - al.count('.'), len(seq), round((len(seq) - al.count('|') - al.count('.'))/len(seq) * 100,1))
        f.write('# Identity:     {}'.format(' ' * (score_width - len(identity)) + identity))
        f.write('# Similarity:   {}'.format(' ' * (score_width - len(similarity)) + similarity))
        f.write('# Gaps:         {}'.format(' ' * (score_width - len(gaps)) + gaps))
        f.write('# Score: {}\n'.format(score))
        f.write('#\n')
        f.write('#\n')
        f.write('#=======================================\n')
        f.write('\n')
        #
        # split alignment strings into per-line chunks for pretty printing
        alignment_chunks = [(seq[i:i+width],al[i:i+width],ref[i:i+width]) for i in range(0,len(seq),width)]
        seq_coord = 1
        ref_coord = 1
        for s,a,r in alignment_chunks:
            seq_coord = print_alignment_seq(s, seq_coord,pre_width,post_width,f) #better return string to write + seq_coord?
            print_alignment_connection(a, pre_width,f)
            ref_coord = print_alignment_seq(r, ref_coord,pre_width,post_width,f) 
            f.write('\n')
        #
        f.write('\n')
        f.write('#---------------------------------------\n')
        f.write('#---------------------------------------\n')


def parallelize(function, args, cores):
    with multiprocessing.Pool(cores) as p:
        return p.map(function, args)




#######################################
# MORE HELPER FUNCTIONS

# transform read coordinates to WT reference coords
def read_to_wt_coord(readn_coord, refn):
  wt_coord = readn_coord - sum(np.where(refn == '-')[0] < readn_coord)
  assert(wt_coord >= 0 and wt_coord < len(refn) - sum(refn == '-')) #wt_coord should be in the range of [0,len(wt_ref)[
  return wt_coord


# check that insert was realigned in one piece
def integral_insert_realignment(insert_alignment, insert_length):
    insert_idxs = [i for i in range(len(insert_alignment)) if insert_alignment[i] != '-']
    return insert_idxs[-1] - insert_idxs[0] +1 == insert_length


# check whether insert requires left normalization, i.e. has an ambiguous alignment and is not fully shifted to the left 
def left_normalize(readn, refn, insert_start, insert_end, i):
    if insert_start > 0 and insert_end < len(readn)-1 and refn[insert_start -1].lower() == readn[insert_end].lower():
        #print("LEFT NORMALIZE: {}".format(i))
        #print(''.join(readn))
        #print(''.join(refn))
        return True
    return False


# filter inserts from a df supported by less than n unique_reads 
def filter_number_unique_reads(df, min_unique_reads):
    fail = [len(x) < min_unique_reads for x in df["idx"]]
    print("Filtered {} / {} inserts supported by < {} distinct reads".format(sum(fail),df.shape[0],min_unique_reads))
    return df.loc[[not x for x in fail]]

# filter inserts from a df supported by less than n total reads
def filter_number_total_reads(df, min_total_reads):
    fail = [x < min_total_reads for x in df["counts"]]
    print("Filtered {} / {} inserts supported by < {} total reads".format(sum(fail),df.shape[0],min_total_reads))
    return df.loc[[not x for x in fail]]

# filter inserts from df that do not pass min_vaf
def filter_vaf(df, min_vaf):
    fail = [x < min_vaf for x in df["vaf"]]
    print("Filtered {} / {} inserts with VAF < {} %".format(sum(fail),df.shape[0],min_vaf))
    return df.loc[[not x for x in fail]]

# filter non-trailing inserts from df whose offset (= insert-tandem distance) does not equal their length --> means insert and tandem are not adjacent
def filter_offset(df):   
    fail = [df["trailing"][i] == False and df["offset"][i] != df["length"][i] for i in range(df.shape[0])]
    #print("Filtering {} / {} ITDs that were not adjacent to their second tandem\n".format(sum(fail),df.shape[0]))
    return df.loc[[not x for x in fail]]  

# filter ITDs
def filter_inserts(df):
    if MIN_UNIQUE_READS > 1:
        df = filter_number_unique_reads(df, MIN_UNIQUE_READS)
    else:
        print("Not filtered based on number of unique supporting reads!")
    if MIN_TOTAL_READS > 1:
        df = filter_number_total_reads(df, MIN_TOTAL_READS)
    if MIN_VAF > 0:
        df = filter_vaf(df, MIN_VAF)
    else:
        print("Not filtered based on VAF!")
    return df



def norm_start_col(df, start_col, ref_wt):
    return [min(max(0,x), len(ref_wt)-1) for x in df["start"]]


# collapse df  
# -> keep: columns containing the same value in all rows -> pick any 
# -> add: columns to sum up (such as total ITD VAF or counts)
# -> max_: columns to keep max of (such as offset or trailing)
# -> append: columns for which all rows should be collapsed to one entry with a single list
# ---> make this nicer, I don't like all of these exceptions for VAF column...
# ---> protect against empty df (will fail for col in append: ... [0] !!
def collapse(df,add=[],max_=[],append=[],keep=[]):
    df_collapsed = df[keep].drop_duplicates().sort_values(by=keep).reset_index(drop=True)
    # vaf is of type Decimal and will be omitted by sum() -> calculate sum separately for vaf column where required
    if "vaf" in add:
        df_collapsed["vaf"] = [sum(x[1]['vaf']) for x in df.groupby(by=keep, as_index=False)]
        add_other = [x for x in add if x != "vaf"]
        df_collapsed[add_other] = df.groupby(by=keep, as_index=False).sum()[add_other]
    else:
        df_collapsed[add] = df.groupby(by=keep, as_index=False).sum()[add]
    df_collapsed[max_] = df.groupby(by=keep, as_index=False).max()[max_]
    for col in append:
        df_collapsed[col] = df.groupby(by=keep, as_index=False)[col].apply(list).reset_index()[0]  #[0] extracts aggregated column / removes df[keep]
    # keep track of which ITD contributed each value of the "add" columns -> create "append" columns for these as well
    for col in add:
        df_collapsed[col + '_each'] = df.groupby(by=keep, as_index=False)[col].apply(list).reset_index()[0] # same command as for col in append
        if col == "vaf":
            assert [sum(x) for x in df_collapsed[col + '_each']] == [dc.Decimal(x) for x in df_collapsed[col]]
        else:
            assert [sum(x) for x in df_collapsed[col + '_each']] == [int(x) for x in df_collapsed[col]]
    assert not df_collapsed.empty
    return df_collapsed


# collect coverage 
def get_coverage(df, start_col, ref_coverage):
    return [ref_coverage[pos] for pos in df[start_col]]

# calculate VAF from insert-containing read counts and total coverage of each mutation
# move asserts to where I first set counts and coverage -> both must be defined ints >= 0
def get_vaf(df): 
    return [dc.Decimal(df.loc[i,"counts"].item()) / dc.Decimal(df.loc[i,"ref_coverage"].item()) * dc.Decimal(100) for i in range(df.shape[0])]


# create and return an empty dataframe
def empty_df(start_col):
    empty_df = pd.DataFrame(columns=['length', 'trailing', start_col, 'insert', 'idx', 'file', 'counts', 'ref_coverage', 'vaf', 'counts_each'])
    if start_col == "tandem2_start":
        empty_df["offset"] = np.nan
        empty_df["offset"] = empty_df["offset"].astype(int)
    empty_df[["length",start_col,"counts","ref_coverage"]] = empty_df[["length",start_col,"counts","ref_coverage"]].astype(int)
    empty_df["idx"] = []
    empty_df["file"] = []
    empty_df["counts_each"] = []
    return empty_df


def fix_dtypes(df):
    for col in ["length","start","tandem2_start","counts","ref_coverage","offset","trailing_end"]:
        if col in df:
            df[col] = df[col].astype(int)
    for col in ["trailing"]:
        if col in df:
            df[col] = df[col].astype(bool)
    for col in ["vaf"]:
        if col in df:
            df[col] = df[col].astype(dc.Decimal)
    return df


def collapse_complex(ii_row, i_row, start_col):
    df_collapsed = pd.Series()
    # some stats will always be the same --> just set these now
    df_collapsed["sample"] = i_row["sample"]
    if "trailing_end" in i_row:
        df_collapsed["trailing_end"] = i_row["trailing_end"]
    # add together some statistics
    df_collapsed["counts"] = ii_row["counts"] + i_row["counts"]
    df_collapsed["vaf"] = ii_row["vaf"] + i_row["vaf"]
    # max some
    df_collapsed["trailing"] = max(i_row["trailing"], ii_row["trailing"])
    #if i_row["trailing"] > 0 + ii_row["trailing"] == 8: # R1 and R2 reads are being merged here!!    # should never happen because of groupby trailing_end!
    #    print("--> CAVE: trailing inserts from R1 and R2 reads are being merged here!")
    #    print(i_row)
    #    print(ii_row)
    #df_collapsed["trailing_end"] = max(i_row["trailing_end"], ii_row["trailing_end"]) # <-- this is always the same because it is part of groupby group!
    #if 'offset' in i_row:
    #    df_collapsed["offset"] = max(i_row["offset"], ii_row["offset"])
    # append some
    df_collapsed["idx"] = i_row["idx"] + ii_row["idx"]
    df_collapsed["file"] = i_row["file"] + ii_row["file"]
    df_collapsed["counts_each"] = i_row["counts_each"] + ii_row["counts_each"]
    # pick one or the other for some more (selection based on number of counts)
    if i_row["counts"] > ii_row["counts"]:
        df_collapsed["length"] = i_row["length"] # will be the same for collapse similar inserts but not necessarily for collapse trailing!
        df_collapsed[start_col] = i_row[start_col] # will be the same for collapse similar inserts but not necessarily for collapse trailing!
        df_collapsed["insert"] = i_row["insert"]        
        df_collapsed["ref_coverage"] = i_row["ref_coverage"]
        if 'offset' in i_row:
            df_collapsed["offset"] = i_row["offset"]
    else:
        df_collapsed["length"] = ii_row["length"]        
        df_collapsed[start_col] = ii_row[start_col]
        df_collapsed["insert"] = ii_row["insert"]        
        df_collapsed["ref_coverage"] = ii_row["ref_coverage"]
        if 'offset' in i_row:
            df_collapsed["offset"] = ii_row["offset"]
    return df_collapsed


# collapse inserts that have the same length and start coordinate and a SIMILAR insert sequence (group = ["length",start_col])
# collapse trailing inserts that have similar length and sequence (group = ["trailing_end"], df should contain trailing mutations only!)
# for collapse trailing: Add column with max_offset collapsed! (keep length then as the one with max counts as usual)
def collapse_similar_inserts(df, start_col, group):
    df_collapsed = empty_df(start_col) # collect all collapsed inserts
    #
    for this_df in [group[1].sort_values(['length'])[::-1] for group in df.groupby(by=group, as_index=False)]: # loop over grouped dfs, sort by length (no effect when grouping by length... avoidable? For trailing muts it ensures that longer trailing sequences are considered first and shorter ones are collapsed onto these (prefers longer insert sequences over shorter ones, as the former are allowed to accumulate more (collapsed) counts over time
        this_df_collapsed = empty_df(start_col)  # collect inserts of the same length and start-coord in this tmp df
        min_score = get_min_score(this_df["insert"].iloc[0], this_df["insert"].iloc[0])
        #
        for i in range(this_df.shape[0]):
            i_idx = this_df.index[i]
            this_ins = this_df["insert"][i_idx]
            collapsed = False
            #
            for ii,ii_row in this_df_collapsed[::-1].iterrows(): #[::-1] to reverse df and speed up pos alignment
                other_ins = ii_row["insert"]
                #
                alignment_score = bio.align.globalcs(this_ins, other_ins, get_alignment_score, COST_GAPOPEN, COST_GAPEXTEND, one_alignment_only=True, penalize_end_gaps=(not (this_df.ix[i_idx]["trailing"] and ii_row["trailing"])), score_only=True)
                if alignment_score >= min_score:
                    # collapse
                    collapsed = True
                    this_df_collapsed.loc[ii] = collapse_complex(ii_row, this_df.ix[i_idx], start_col)
                    #this_df_collapsed.loc[ii, ['counts','vaf','trailing','idx','file','counts_each','insert','ref_coverage']] = collapse_complex(ii_row, this_df.ix[i_idx])
                    break # once collapsed, need not look further for this insert -> move on to the next one
            #
            if not collapsed:
                this_df_collapsed = this_df_collapsed.append(this_df.ix[i_idx], ignore_index=True)
        # append collapsed inserts of current length and start-coord to the overall collecting df
        df_collapsed = df_collapsed.append(this_df_collapsed, ignore_index=True)
    #
    # check that sum of "counts_each" (= read counts of each unique read) equals total counts in "counts"
    df_collapsed = fix_dtypes(df_collapsed)
    assert([sum(x) for x in df_collapsed["counts_each"]] == [int(x) for x in df_collapsed["counts"]])
    return df_collapsed



def index_of_max(lst):
    return lst.index(max(lst))


# collapse inserts that have the same length and SIMILAR insert sequence and are within one insert_length of each other
def collapse_close_inserts(df, start_col):
    df_collapsed = empty_df(start_col) # collect all collapsed inserts
    #
    for length in set(df["length"]):
        this_df = df.ix[df["length"] == length]
        this_df_collapsed = empty_df(start_col)  # collect inserts of the same length and start-coord in this tmp df -> use this one for collapsing!
        min_score = get_min_score(this_df["insert"].iloc[0], this_df["insert"].iloc[0])
        #
        for i in range(this_df.shape[0]):
            i_idx = this_df.index[i]
            this_max_counts = index_of_max(this_df["counts_each"][i_idx])
            this_read = all_reads[this_df["idx"][i_idx][this_max_counts]]
            this_ins_start = this_df[start_col][i_idx]
            collapsed = False
            #
            for ii,ii_row in this_df_collapsed[::-1].iterrows(): #[::-1] to reverse df and speed up pos alignment
                other_ins_start = ii_row[start_col]
                if abs(other_ins_start - this_ins_start) > length:
                    continue
                other_max_counts = index_of_max(ii_row["counts_each"])
                other_read = all_reads[ii_row["idx"][other_max_counts]]
                min_start = min(this_ins_start, other_ins_start)
                max_start = max(this_ins_start, other_ins_start)
                #
                alignment_score = bio.align.globalcs(this_read[min_start:max_start+length], other_read[min_start:max_start+length], get_alignment_score, COST_GAPOPEN, COST_GAPEXTEND, one_alignment_only=True, penalize_end_gaps=True, score_only=True)
                if alignment_score >= min_score:
                    # collapse
                    collapsed = True
                    this_df_collapsed.loc[ii] = collapse_complex(ii_row, this_df.ix[i_idx], start_col)
                    #this_df_collapsed.loc[ii, ['counts','vaf','trailing','idx','file','counts_each','insert','ref_coverage']] = collapse_complex(ii_row, this_df.ix[i_idx])
                    break # once collapsed, need not look further for this insert -> move on to the next one
            #
            if not collapsed:
                this_df_collapsed = this_df_collapsed.append(this_df.ix[i_idx], ignore_index=True)
        # append collapsed inserts of current length and start-coord to the overall collecting df
        df_collapsed = df_collapsed.append(this_df_collapsed, ignore_index=True)
    #
    # check that sum of "counts_each" (= read counts of each unique read) equals total counts in "counts"
    df_collapsed = fix_dtypes(df_collapsed)
    assert([sum(x) for x in df_collapsed["counts_each"]] == [int(x) for x in df_collapsed["counts"]])
    return df_collapsed



# update length of trailing ITDs to offset instead (should be max potential ITD length)
# -> update also insert sequence?
def fix_trailing_length(df):
    df_fixed = df.copy()
    df_fixed.loc[df["trailing"],"length"] = df_fixed.ix[df["trailing"],"offset"]
    return df_fixed.sort_values(['length','tandem2_start'])
    


#######################################
#######################################
## READ INPUT & CREATE OUTPUT FOLDERS
#######################################
#######################################

if __name__ == '__main__':
    print("==== PROCESSING SAMPLE {} ====".format(SAMPLE))
    #
    ## CREATE OUTPUT FOLDER
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    #
    #
    ## GET READS FROM FASTQ 
    print("-- Reading FASTQ files --")
    # collect reads that pass MIN_BQS filter -> reverse-complement R2 reads so that all reads can be aligned to the same reference
    reads_and_bqs = read_fastq(R1) + [(reverse_complement(r2_read),bqs) for r2_read,bqs in read_fastq(R2)]
    print("Number of total reads: {}".format(len(reads_and_bqs)))
    # trim ambiguous 'N' bases
    trimmed = parallelize(trim_n, reads_and_bqs, NKERN)
    #
    # filter based on BQS -> PASS returns read, FAIL returns None -> remove None from list!
    reads = None
    args = [(read_and_bqs, MIN_BQS) for read_and_bqs in trimmed]
    if MIN_BQS > 0:
        reads = [x for x in parallelize(filter_bqs, args, NKERN) if x is not None]
    else:
        reads = [x[0] for x in trimmed]
    #
    #
    # get unique reads and counts thereof
    print("Number of total reads with mean BQS >= {}: {}".format(MIN_BQS,len(reads)))
    tmp = collections.Counter(reads)
    unique_reads = list(tmp.keys())
    unique_reads_counts = list(tmp.values())
    assert len(unique_reads) == len(unique_reads_counts)
    print("Number of unique reads with mean BQS >= {}: {}".format(MIN_BQS,len(unique_reads)))
    #
    #
    # filter unique reads -> keep only reads that exist at least twice  -----> make this assignment nicer, it's like the same loop twice!!?
    if MIN_READ_COPIES == 1:
        print("Turned OFF unique reads filter!")
    else:
        unique_reads, unique_reads_counts = [read  for read,count in zip(unique_reads, unique_reads_counts) if count >= MIN_READ_COPIES], [count for read,count in zip(unique_reads, unique_reads_counts) if count >= MIN_READ_COPIES]
        print("Number of unique reads with at least {} copies: {}".format(MIN_READ_COPIES,len(unique_reads)))
    assert len(unique_reads) == len(unique_reads_counts)
    print("--> Number of total reads passing filters: {}".format(sum(unique_reads_counts)))
    #
    #
    ## GET REFERENCE and ANNOTATION
    # --> compare variable with those below (below it's ref_wt and list of chars instead of string?)
    wt_ref = get_reference(REF)
    wt_ref_upper = wt_ref.upper()
    anno = get_annotation(ANNO)
    #
    #
    #
    ## DO ALIGNMENTS & FILTER BASED ON ALIGNMENT SCORE
    print("\n-- Aligning to Reference --")
    all_alignments = None
    args = [(unique_read, wt_ref_upper) for unique_read in unique_reads]
    with multiprocessing.Pool(NKERN) as p:
            all_alignments = p.map(align, args)
    assert len(unique_reads) == len(all_alignments)
    #
    print("Filtering {} / {} low quality alignments with a score < {} % of max".format(all_alignments.count([]),len(all_alignments),  MIN_SCORE_ALIGNMENTS *100))
    #
    #
    all_readCounts  = [unique_reads_counts[i] for i in range(len(all_alignments)) if all_alignments[i] != []]
    all_alignments  = [x for x in all_alignments if x != []]
    all_reads       = [x[0] for x in all_alignments]
    all_refs        = [x[1] for x in all_alignments]
    all_scores      = [x[2] for x in all_alignments]
    all_files       = ['needle_{}.txt'.format(i) for i in range(len(all_alignments))]  ## --> do I even need this anymore?
    print("Total reads remaining for analysis: {}".format(sum(all_readCounts)))
    #
    #
    ## PRINT ALIGNMENTS
    # create output file directory for alignments print-outs
    needle_dir = os.path.join(OUT_DIR,'out_needle')
    if not os.path.exists(needle_dir):
        os.makedirs(needle_dir)
    #
    for i in range(len(all_alignments)):
        print_alignment(all_alignments[i], i, needle_dir, command='bio.align.globalcs', command_seq='unique_reads[i]', command_ref='wt_ref')
    #
    #
    # make sure there is a 1:1 matching between these files -> must all have the same length!  -> since I am reading out of all_alignments, this should always be the case?
    assert(len(all_reads) == len(all_readCounts))
    assert(len(all_reads) == len(all_alignments))
    assert(len(all_reads) == len(all_refs))
    assert(len(all_reads) == len(all_scores))
    assert(len(all_reads) == len(all_files))
    #
    #
    if len(all_reads) == 0:
        print("\nNO READS TO PROCESS!")
        quit()
    #
    #
    #######################################
    # EXTRACT INSERT SEQUENCE FROM READ
    print("\n-- Looking for insertions & ITDs --")
    #
    # check each alignment for insert/itd and save index in all_reads/all_refs/all_files to idx, insert/itd length to length and insert/itd start/stop position to start/end dicts based on insert/itd classification
    w_ins = {"idx": [], "file": [], "length": [], "start": [], "insert": [], "trailing": [], "trailing_end": []}
    w_itd_exact = {"idx": [], "file": [], "length": [], "start": [], "tandem2_start": [], "offset": [], "insert": [], "trailing": [], "trailing_end": []}
    w_itd_nonexact_fail = {"idx": [], "file": [], "length": [], "start": [], "insert": [], "trailing": [], "trailing_end": []}
    w_itd_nonexact = {"idx": [], "file": [], "length": [], "start": [], "tandem2_start": [], "offset": [], "insert": [], "trailing": [], "trailing_end": []}
    #
    ref_wt = [base for base in all_refs[0] if base != '-'] 
    ref_coverage = np.zeros(len(ref_wt)) # count number of reads covering each bp AND its successor (therefore do not calc coverage for last bp)
    #
    ambig_i = []
    ambig_als = []
    should_left_normalize = 0
    #
    # loop over all alignments, test for presence of an ITD
    for read,ref,score,counts,filename,i in zip(all_reads, all_refs, all_scores, all_readCounts, all_files, range(len(all_reads))):
        readn = np.array(list(read))
        refn = np.array(list(ref))
        assert(len(readn) == len(refn))
#
        # calculate ref_coverage
        readn_covered = np.where(readn[refn != '-'] != '-')
        readn_covered_range = np.arange(np.min(readn_covered), np.max(readn_covered)) # do not count last index -> read ending here holds no information on whether or not an ITD starts here
        ref_coverage[readn_covered_range] = ref_coverage[readn_covered_range] + counts
#	
        # get indeces of inserts in read (positions where reference has a gap and read does not)
        insert_idxs_all = np.arange(len(readn))[refn == '-']
        assert('-' not in readn[insert_idxs_all]) # two gaps should never align at the same pos!
#	
        # get indeces of each individual insert in each read
        insert_idxs_each = []
        insert_idxs_tmp = []
        i_prev = None
        for i_this in insert_idxs_all:
            if i_prev is None or i_prev == i_this -1: #start saving first/continue saving next insert index
                insert_idxs_tmp.append(i_this)
                i_prev = i_this
            else: #save current insert_idxs_tmp list and open up a new one for the next insert
                insert_idxs_each.append(insert_idxs_tmp)
                insert_idxs_tmp = [i_this]
                i_prev = i_this
        # save last insert as well
        insert_idxs_each.append(insert_idxs_tmp)
        assert(np.all(np.concatenate(insert_idxs_each) == insert_idxs_all))		
#
        for insert_idxs in insert_idxs_each:
            insert_length = len(insert_idxs)	
            if insert_length < 6 or "N" in readn[insert_idxs]: 
                continue
            ins = readn[insert_idxs]
            insert_start = insert_idxs[0]
            insert_end = insert_idxs[-1]
            #trailing = insert_end == sum(readn != '-')-1 # or insert_start == 0 --> there are no 5' trailing ITDs! (if they are shorter than read remainder, they will be contained within it (thus not trailing), otherwise alignment will switch so that 5' part will be aligned and shorter 3' end will be the insert
            trailing_end = 0
            if (readn[0:insert_start] == '-').all():
                trailing_end = 5
            if (readn[insert_end+1:] == '-').all():
                assert trailing_end != 5 # read would consist of insert only if this were the case
                trailing_end = 3
            trailing = insert_start != 0 and insert_end != len(readn)-1 and (trailing_end == 5 or trailing_end == 3) #(readn[0:insert_start] == '-').all() or (readn[insert_end+1:] == '-').all())
            # if there is an insert  --> require min 6 bp length, in-frame insert (except for trailing muts) and no "N"s within insert
            if(trailing or insert_length % 3 == 0):
                if insert_start > 0:
                    should_left_normalize = should_left_normalize + left_normalize(readn, refn, insert_start, insert_end, i)
#
                # relative to the reference, get coord of the first WT base before insert	
                insert_start_ref = insert_start - sum(refn[0:insert_start] == '-')  
                if insert_start == 0: 
                    insert_start_ref = insert_start_ref - insert_length
#	    
                w_ins["idx"].append(i)
                w_ins["file"].append(filename)
                w_ins["length"].append(insert_length)
                w_ins["start"].append(insert_start_ref)
                w_ins["insert"].append(''.join(ins))
                w_ins["trailing"].append(trailing)
                w_ins["trailing_end"].append(trailing_end)
#			
                # check whether the insert is contained within non-insert read a second time -> that makes it an ITD!
                readn_nonIns = np.delete(readn,insert_idxs)
                readn_maskedIns = readn.copy()
                readn_maskedIns[insert_idxs] = 'Z' # wild base for "no base" -> prevent alignment to already detected insert
#
                # search for nearest tandem before and after ITD
                tandem2_after = ''.join(readn_maskedIns).find(''.join(ins), insert_start + insert_length,len(readn_maskedIns))
                tandem2_before = ''.join(reversed(readn_maskedIns)).find(''.join(reversed(ins)), len(readn_maskedIns) -1 -insert_start +1, len(readn_maskedIns))
# 
                # take the one closest to the insert (should be relevant only for small ITDs that may be contained multiple times within a read
                tandem2_start = None
                if tandem2_after == -1 and tandem2_before == -1:
                    tandem2_start = -1 # not found --> no itd present
                elif tandem2_after == -1:
                    tandem2_start = len(readn_maskedIns) -1 -tandem2_before -insert_length +1  # convert coords back from reverse to forward sense
                elif tandem2_before == -1:
                    tandem2_start = tandem2_after
                elif tandem2_after < tandem2_before:
                    tandem2_start = tandem2_after
                elif tandem2_before < tandem2_after:
                    tandem2_start = len(readn_maskedIns) -1 -tandem2_before -insert_length +1  # convert coords back from reverse to forward sense
                assert tandem2_start is not None  # should be assigned something!
#                    
                offset = abs(tandem2_start - insert_start)
                if trailing and offset == insert_length:  # should test whether untrailed mut is in-frame at this point! (can I undo saving as insert if it is not?? Otherwise get rid of w_ins and create df_ins from df_itd_exact/nonexact/nonexactfail)
                    trailing = False ## should fix this in w_ins as well!!!
                    print("UNTRAILED") 
                    if insert_length % 3 != 0:
                        print("BUT NOT IN FRAME!!!")
                # save if an exact second tandem of the insert was found
                if tandem2_start != -1:   # ---> also check that index of second match is sufficiently close to insert! (for exact match and alignment approach!)
                    w_itd_exact["idx"].append(i)
                    w_itd_exact["file"].append(filename)
                    w_itd_exact["length"].append(insert_length)
                    w_itd_exact["start"].append(insert_start_ref)
                    w_itd_exact["tandem2_start"].append(read_to_wt_coord(tandem2_start, refn))
                    w_itd_exact["offset"].append(offset)
                    w_itd_exact["insert"].append(''.join(ins))
                    w_itd_exact["trailing"].append(trailing)
                    w_itd_exact["trailing_end"].append(trailing_end)
                    #if trailing:
                    #    w_itd_exact["length"][-1] = w_itd_exact["offset"][-1]
                else:
                    # otherwise search for sufficiently similar (> 90 % bases mapped) second tandem by realignment of the insert within the remainder of the read
                    max_score = len(ins) * COST_MATCH  # param?
                    min_score = max_score * MIN_SCORE
                    # arguments: seq1, seq2, match-score, mismatch-score, gapopen-score, gapextend-score --> match/mismatch from needle default (/usr/share/EMBOSS/data/EDNAFULL), gap as passed to needle in my script
                    # output: list of optimal alignments, each a list of seq1, seq2, score, start-idx, end-idx 
                    alignments = bio.align.localcs(''.join(ins), ''.join(readn_maskedIns), get_alignment_score, COST_GAPOPEN, COST_GAPEXTEND)
                    # filter alignments where insert cannot be realigned in one piece
                    alignments = [al for al in alignments if integral_insert_realignment(al[0],insert_length)]
                    alignment_score = None
                    #
                    if alignments == []:
                        alignment_score = -1
                    else:
                        alignment = alignments[0]
                        alignment_score, alignment_start, alignment_end = alignment[2:5]
                    # check that min alignment score is reach and that ref covers at least part of 2nd tandem (for an internal duplication, sequence must not only be duplicated but also contained in ref. I am realigning the insert to the remainder of the read (to not penalize true germline SNPs in the ITD) so I must make sure at this point that TD is in fact internal)
                    if alignment_score >= min_score and alignment_start < max(np.where(refn != '-')[0]):
                        offset = abs(alignment_start - insert_start)
                        w_itd_nonexact["idx"].append(i)
                        w_itd_nonexact["file"].append(filename)
                        w_itd_nonexact["length"].append(insert_length)
                        w_itd_nonexact["start"].append(insert_start_ref)
                        w_itd_nonexact["tandem2_start"].append(read_to_wt_coord(alignment_start, refn))
                        w_itd_nonexact["offset"].append(offset)
                        w_itd_nonexact["insert"].append(''.join(ins))
                        w_itd_nonexact["trailing"].append(trailing)
                        w_itd_nonexact["trailing_end"].append(trailing_end)
                    else:
                        w_itd_nonexact_fail["idx"].append(i)
                        w_itd_nonexact_fail["file"].append(filename)
                        w_itd_nonexact_fail["length"].append(insert_length)
                        w_itd_nonexact_fail["start"].append(insert_start_ref)
                        w_itd_nonexact_fail["insert"].append(''.join(ins))
                        w_itd_nonexact_fail["trailing"].append(trailing)
                        w_itd_nonexact_fail["trailing_end"].append(trailing_end)
                        #print(bio.format_alignment(*alignment))
                    if len(alignments) > 1:
                        ambig_i.append(i)
                        ambig_als.append(alignments)
    #
    # print number of ambiguous alignments (to see if this is sth I need to handle or not)
    #print("There were {} inserts that generated ambiguous alignments.".format(len(ambig_i)))
    #print("There were {} inserts whose alignment should have been left normalized.".format(should_left_normalize))
    #print()
    #print("There were {} trailing inserts.".format(sum(w_ins["trailing"])))
    #print("There were {} trailing exact ITDs.".format(sum(w_itd_exact["trailing"])))
    #print("There were {} trailing nonexact ITDs.".format(sum(w_itd_nonexact["trailing"])))
    #print("There were {} trailing nonexact ITDs failed.".format(sum(w_itd_nonexact_fail["trailing"])))
    #print()
    #
    #
    # fix ref_coverage -> coverage of last index is 0 since I am counting reads covering a position AND its successor -> for the final index, there is no successor but also this restraint is unnecessary (any trailing mut will be covered by any read covering the first/last base!)
    ref_coverage[-1] = ref_coverage[-2]
    ref_coverage = ref_coverage.astype(int)
    #
    #
    ########################################
    # COLLECT AND COLLAPSE ITDs
    print("-- Collecting results --")
    #
    #
    df_ins =  pd.DataFrame(w_ins)
    if not df_ins.empty:
        cols = ['sample','length', 'start', 'vaf', 'ref_coverage', 'counts', 'insert']
        df_ins["sample"] = [SAMPLE for i in range(df_ins.shape[0])]
        df_ins["counts"] = [all_readCounts[i] for i in df_ins["idx"]]
        #
        # collapse same inserts (same length, insert sequence and reference-based start coordinate)
        #df_ins_grouped = collapse(df_ins,keep=["sample","length","start","insert"],add=["counts"],max_=["trailing"],append=["idx","file"])
        df_ins_grouped = collapse(df_ins,keep=["sample","length","start","insert","trailing","trailing_end"],add=["counts"],max_=[],append=["idx","file"])
        assert sum(df_ins["counts"]) == sum(df_ins_grouped["counts"])
        #df_ins_grouped["norm_start"] = norm_start_col(df_ins_grouped, "start", ref_wt)
        df_ins_grouped["start"] = norm_start_col(df_ins_grouped, "start", ref_wt)
        df_ins_grouped["ref_coverage"] = get_coverage(df_ins_grouped, "start", ref_coverage) # should I be using norm_start at some other place as well?? Or is it just for coverage calculation?!
        #df_ins_grouped["ref_coverage"] = get_coverage(df_ins_grouped, "norm_start", ref_coverage) # should I be using norm_start at some other place as well?? Or is it just for coverage calculation?!
        df_ins_grouped["vaf"] = get_vaf(df_ins_grouped)
        df_ins_grouped[cols + ['counts_each', 'file']].to_csv(os.path.join(OUT_DIR,"flt3_ins.tsv"), index=False, float_format='%.2e', sep='\t')
        #
        # COLLAPSE #2
        # --> align inserts of same length and tandem2_start, collapse if they are sufficiently similar
        df_ins_collapsed = collapse_similar_inserts(df_ins_grouped, "start",group=["length","start"]).sort_values(['length','start'])
        assert sum(df_ins["counts"]) == sum(df_ins_collapsed["counts"])
        df_ins_collapsed[cols].to_csv(os.path.join(OUT_DIR,"flt3_ins_collapsed-similar.tsv"), index=False, float_format='%.2e', sep='\t')
        #
        # COLLAPSE #3
        # --> align inserts of same length,  collapse if they are sufficiently similar (and within one insert length of each other)
        df_ins_collapsed = collapse_close_inserts(df_ins_collapsed, "start").sort_values(['length','start'])
        assert sum(df_ins["counts"]) == sum(df_ins_collapsed["counts"])
        df_ins_collapsed[cols].to_csv(os.path.join(OUT_DIR,"flt3_ins_collapsed-similar-close.tsv"), index=False, float_format='%.2e', sep='\t')
        # 
        # COLLAPSE #4
        # --> collapse trailing
        df_ins_collapsed = pd.concat([df_ins_collapsed.ix[df_ins_collapsed["trailing"] == False], collapse_similar_inserts(df_ins_collapsed.ix[df_ins_collapsed["trailing"]==True], "start",group=["trailing_end"]).sort_values(['length','start'])])
        assert sum(df_ins["counts"]) == sum(df_ins_collapsed["counts"])
        annotate(df_ins_collapsed, anno, 'start')[cols + list(anno)[1:]].to_csv(os.path.join(OUT_DIR,"flt3_ins_collapsed-similar-close-trailing.tsv"), index=False, float_format='%.2e', sep='\t')
        # 
        # FILTER
        # --> filter inserts based on number of unique and total supporting reads
        #if 'cr' not in SAMPLE: # change this to some binary flag
        print("\n-- Filtering insertions --")
        df_ins_collapsed = filter_inserts(df_ins_collapsed).sort_values(['length','start'])
        annotate(df_ins_collapsed, anno, 'start')[cols + list(anno)[1:]].to_csv(os.path.join(OUT_DIR,"flt3_ins_collapsed-similar-close-trailing_hc.tsv"), index=False, float_format='%.2e', sep='\t')
        print("Found {} insertions!\n".format(df_ins_collapsed.shape[0]))
    else:
        print("Found 0 insertions!\n")
    #
    #
    #
    df_itd =  pd.concat([pd.DataFrame(w_itd_exact), pd.DataFrame(w_itd_nonexact)], ignore_index=True)
    df_itd[["idx","length","offset","start","tandem2_start","trailing_end"]] = df_itd[["idx","length","offset","start","tandem2_start","trailing_end"]].astype(int)
    df_itd["trailing"] = df_itd["trailing"].astype(bool) # wasn't there a function to fix dtypes??
    df_itd = filter_offset(df_itd)
    if not df_itd.empty:
        cols = ['sample','length', 'trailing', 'tandem2_start', 'vaf', 'ref_coverage', 'counts', 'insert']
        df_itd["sample"] = [SAMPLE for i in range(df_itd.shape[0])]
        df_itd["counts"] = [all_readCounts[i] for i in df_itd["idx"]]
        #
        # collapse same inserts (same length, insert sequence and reference-based start coordinate)
        # -> careful: losing "start" column here (only keeping tandem2_start) -> do I need it?
        #df_itd_grouped = collapse(df_itd,keep=["sample","length","tandem2_start","insert"],add=["counts"],max_=["trailing","offset"],append=["idx","file"])
        df_itd_grouped = collapse(df_itd,keep=["sample","length","tandem2_start","insert","trailing","trailing_end"],add=["counts"],max_=["offset"],append=["idx","file"])
        assert sum(df_itd["counts"]) == sum(df_itd_grouped["counts"])
        df_itd_grouped["ref_coverage"] = get_coverage(df_itd_grouped, "tandem2_start", ref_coverage)
        df_itd_grouped["vaf"] = get_vaf(df_itd_grouped)
        fix_trailing_length(df_itd_grouped)[cols + ['counts_each', 'file']].to_csv(os.path.join(OUT_DIR,"flt3_itds.tsv"), index=False, float_format='%.2e', sep='\t')
        #
        # COLLAPSE #2
        # --> align inserts of same length and tandem2_start, collapse if they are sufficiently similar
        df_itd_collapsed = collapse_similar_inserts(df_itd_grouped, "tandem2_start",group=["length","tandem2_start","trailing_end"]).sort_values(['length','tandem2_start'])
        assert sum(df_itd["counts"]) == sum(df_itd_collapsed["counts"])
        fix_trailing_length(df_itd_collapsed)[cols].to_csv(os.path.join(OUT_DIR,"flt3_itds_collapsed-similar.tsv"), index=False, float_format='%.2e', sep='\t')
        #
        # COLLAPSE #3
        # --> align inserts of same length,  collapse if they are sufficiently similar (and within one insert length of each other)
        df_itd_collapsed = collapse_close_inserts(df_itd_collapsed, "tandem2_start").sort_values(['length','tandem2_start'])
        assert sum(df_itd["counts"]) == sum(df_itd_collapsed["counts"])
        fix_trailing_length(df_itd_collapsed)[cols].to_csv(os.path.join(OUT_DIR,"flt3_itds_collapsed-similar-close.tsv"), index=False, float_format='%.2e', sep='\t')
        #
        # COLLAPSE #4
        # --> collapse trailing
        df_itd_collapsed = pd.concat([df_itd_collapsed.ix[df_itd_collapsed["trailing"] == False], collapse_similar_inserts(df_itd_collapsed.ix[df_itd_collapsed["trailing"]==True], "tandem2_start",group=["trailing_end"]).sort_values(['length','tandem2_start'])])
        assert sum(df_itd["counts"]) == sum(df_itd_collapsed["counts"])
        annotate(fix_trailing_length(df_itd_collapsed), anno, 'tandem2_start')[cols + list(anno)[1:]].to_csv(os.path.join(OUT_DIR,"flt3_itds_collapsed-similar-close-trailing.tsv"), index=False, float_format='%.2e', sep='\t')
        # 
        # FILTER
        # --> filter inserts based on number of unique and total supporting reads
        #if 'cr' not in SAMPLE: # change this to some binary flag
        print("-- Filtering ITDs --")
        df_itd_collapsed = filter_inserts(df_itd_collapsed).sort_values(['length','tandem2_start'])
        annotate(fix_trailing_length(df_itd_collapsed), anno, 'tandem2_start')[cols + list(anno)[1:]].to_csv(os.path.join(OUT_DIR,"flt3_itds_collapsed-similar-close-trailing_hc.tsv"), index=False, float_format='%.2e', sep='\t')
        print("Found {} ITDs!".format(df_itd_collapsed.shape[0]))
        print()
    else:
        print("Found 0 ITDs!\n")
    #   
    #
    #
    ########################################
    # PRINT SUMMARY STATISTICS on the number of reads in each category
    #
    #print("\nNumber of unique reads supporting each type of insert")
    #print("Insertions: {}".format(len(w_ins["idx"])))
    #print("Single exact ITD: {}".format(len(w_itd_exact["idx"])))
    #print("Single non-exact ITD: {}".format(len(w_itd_nonexact["idx"])))
    #print("Single insertion failed alignment: {}".format(len(w_itd_nonexact_fail["idx"])))
    assert len(w_ins["idx"]) == len(w_itd_exact["idx"]) + len(w_itd_nonexact["idx"]) + len(w_itd_nonexact_fail["idx"])
    #
    #
    #
    ########################################
    # GET KNOWN ITDs 
    # --> if length known is given, extract relevant ITDs 
    # --> if vaf known is given, compare numbers? (vaf is always known if length is known? is there a better way to supply this together?)
    #
    if KNOWN_LENGTH_FILE is not None:
        known_length = read_known(KNOWN_LENGTH_FILE,int)
       # df_itd_known=None
       # df_ins_known=None
        df_itd_known = get_known(fix_trailing_length(df_itd_collapsed), known_length)
        df_itd_known[['sample','length','vaf','ref_coverage','counts','tandem2_start','insert']].to_csv(os.path.join(OUT_DIR,"flt3_itds_collapsed-similar-close-trailing_hc_known.tsv"), index=False, float_format='%.2e', sep='\t', na_rep='NA')
        #
        df_ins_known = get_known(df_ins_collapsed, known_length)
        df_ins_known[['sample','length','vaf','ref_coverage','counts','start','insert']].to_csv(os.path.join(OUT_DIR,"flt3_ins_collapsed-similar-close-trailing_hc_known.tsv"), index=False, float_format='%.2e', sep='\t', na_rep='NA')
        #
        #
        # associate detected ITDs with expected VAF (if available)  -> useful to check correlation between VAF estimates of different experiments --> right now assuming there is only one VAF/AR for sum of all ITD clones!
        known_vaf = None
        if KNOWN_VAF_FILE is not None:
            known_vaf = read_known(KNOWN_VAF_FILE,dc.Decimal)[0]
        elif KNOWN_AR_FILE is not None:
            known_vaf = ar_to_vaf(read_known(KNOWN_AR_FILE,dc.Decimal)[0])
        #
        if known_vaf is not None:
            assert known_vaf <= 100 and known_vaf >= 0
        #
        # does this make sense with multiple inserts per read? counts/vaf would be messed up because counted twice, right? --> more accurate maybe: collect all supporting reads and count unique 
        df_itd_known_collapsed = collapse(df_itd_known,keep=["sample"],add=["counts","vaf"],append=["length","tandem2_start","ref_coverage"])
        df_itd_known_collapsed["vaf_genescan"] = known_vaf
        df_itd_known_collapsed[['sample','length','vaf','vaf_genescan','vaf_each','tandem2_start','ref_coverage','counts_each']].to_csv(os.path.join(OUT_DIR,"flt3_itds_collapsed-similar-close-trailing_hc_known_collapsed.tsv"), index=False, float_format='%.2e', sep='\t', na_rep='NA')
        # 
        #
        df_ins_known_collapsed = collapse(df_ins_known,keep=["sample"],add=["counts","vaf"],append=["length","start","ref_coverage"])
        df_ins_known_collapsed["vaf_genescan"] = known_vaf
        df_ins_known_collapsed[['sample','length','vaf','vaf_genescan','vaf_each','start','ref_coverage','counts_each']].to_csv(os.path.join(OUT_DIR,"flt3_ins_collapsed-similar-close-trailing_hc_known_collapsed.tsv"), index=False, float_format='%.2e', sep='\t', na_rep='NA')
        


