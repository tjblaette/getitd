import Bio.pairwise2 as bio
import timeit
import collections
import datetime
import multiprocessing
import argparse
import pandas as pd
import numpy as np
import decimal as dc
dc.getcontext().prec = 5
import pprint
import os
import copy

#######################################
## INITIALIZE VARIABLES

class Read(object):
    
    def __init__(self, seq, sense=1, bqs=None, counts=1, 
            al_score=None, al_seq=None, al_ref=None, al_file=None, al_index=None):
        self.seq = seq
        self.bqs = bqs
        self.length = len(seq)
        self.counts = counts
        self.sense = sense
        self.al_score = al_score
        self.al_seq = al_seq
        self.al_ref = al_ref
        self.al_file = al_file
        self.al_index = al_index # get rid of this? (or al_file)
        assert self.counts > 0
        if self.bqs is not None:
            assert len(self.seq) == len(self.bqs)
    
    def print(self):
        pprint.pprint(vars(self))
    
    # reverse complement and return a given read
    def reverse_complement(self):
        self.seq = self.seq.translate(str.maketrans('ATCGatcg','TAGCtagc'))[::-1]
        if self.bqs is not None:
            self.bqs = self.bqs[::-1]
        if not self.sense:
            self.sense = -1
        else:
            self.sense = self.sense * -1
        return self 
    
    # trim ambinguous N bases at reads' ends
    def trim_n(self): # rewrite using str.startswith / endswith?
        base_is_n = [x in 'nN' for x in self.seq]
        n_start,n_end = 0,0
        while base_is_n.pop():
            n_end = n_end + 1
        base_is_n.reverse()
        while base_is_n.pop():
            n_start = n_start + 1
        self.seq = self.seq[n_start:self.length - n_end]
        if self.bqs is not None:
            self.bqs = self.bqs[n_start:self.length - n_end]
            assert len(self.seq) == len(self.bqs)
        self.length = len(self.seq)
        return self
    
    # calculate average BQS
    def average_bqs(self):
        return sum([ord(x) - 33 for x in self.bqs]) / len(self.bqs)

    # filter reads based on average BQS
    def filter_bqs(self):
        global MIN_BQS
        # if bqs is None, don't try to filter
        if self.bqs is None or self.average_bqs() >= MIN_BQS:
            return self
        return None
    
    # align read to ref
    def align(self):
        global REF, TECH, COST_GAPOPEN, COST_GAPEXTEND
        # one_alignment_only until more are handled
        alignment = bio.align.globalcs(self.seq, REF, get_alignment_score,
            COST_GAPOPEN, COST_GAPEXTEND, penalize_end_gaps=False, one_alignment_only=True)
        if alignment:
            self.al_seq, self.al_ref, self.al_score = alignment[0][0:3]
        # for 454, I don't know which read is forward, which is reverse -> try both
        if TECH == '454':
            rev = self.reverse_complement()
            rev_alignment = bio.align.globalcs(rev.seq, REF, get_alignment_score,
                COST_GAPOPEN, COST_GAPEXTEND, penalize_end_gaps=False, one_alignment_only=True)
            if (rev_alignment and
                    (self.al_score is None or rev_alignment[0][2] > self.al_score)):
                rev.al_seq, rev.al_ref, rev.al_score = rev_alignment[0][0:3]
                return rev
        return self


class Insert(object):
    def __init__(self, seq, start, end, counts,
            trailing=None, trailing_end=None, reads=None, coverage=None, vaf=None):
        self.seq = seq
        self.length = len(seq)
        self.start = start
        self.end = end
        self.trailing = trailing
        self.trailing_end = trailing_end
        self.counts = counts
        self.coverage = coverage
        self.vaf = vaf
        if reads is None: # passing mutable objects in def() will have them shared between instances!
            reads = []
        self.reads = reads
    
    # function to call for sorting list of inserts by seq
    def get_seq(self):
        return self.seq

    def calc_vaf(self):
        self.vaf = dc.Decimal(self.counts) / self.coverage * 100
        assert self.vaf >= 0 and self.vaf <= 100
        return self

    def norm_start(self):
        global REF
        self.start = min(max(0,self.start), len(REF)-1)
        return self

    # inserts
    # norm start coords to [0,len(REF)[ before printing
    def prep_for_save(self):
        to_save = copy.deepcopy(self)
        to_save = to_save.norm_start()
        return to_save

    def annotate_domains(self, DOMAINS):
        domains = []
        for domain,start,end in DOMAINS:
            if self.start <= end and self.end >= start:
                domains.append(domain)
        return domains

    
    def print(self):
        # --> don't print reads, they only clutter the screen
        pprint.pprint({key: vars(self)[key] for key in vars(self).keys() if not key == 'reads'})

    def is_close_to(self, that):
        if hasattr(self, 'tandem2_start'):
            return abs(self.tandem2_start - that.tandem2_start) <= self.length
        return abs(self.start - that.start) <= 2 * self.length
    
    # align two inserts' seq
    # --> consider similar when al_score >= cutoff
    def is_similar_to(self, that):
        global COST_GAPOPEN, COST_GAPEXTEND, MIN_SCORE_INSERTS
        min_score = get_min_score(self.seq, that.seq, MIN_SCORE_INSERTS)
        al_score = bio.align.globalcs(self.seq, that.seq, get_alignment_score, COST_GAPOPEN, COST_GAPEXTEND, one_alignment_only=True, score_only=True, penalize_end_gaps=False)
        if al_score >= min_score:
            return True
        return False
    
    def should_merge(self, that, condition):
        if condition == 'is_same':
            return self.seq == that.seq and self.start == that.start
        if condition == 'is_similar':
            return self.length == that.length and self.start == that.start and self.is_similar_to(that)
        if condition == 'is_close':
            return self.length == that.length and self.is_similar_to(that) and self.is_close_to(that)
        if condition == 'is_same_trailing':
            return self.trailing and that.trailing and self.trailing_end == that.trailing_end and self.is_similar_to(that)
        if condition == 'any':
            return ((self.length == that.length and self.is_close_to(that)) or (self.trailing and that.trailing and self.trailing_end == that.trailing_end)) and self.is_similar_to(that)
        assert False # should never be reached!
    
    
    # filter based on number of unique supporting reads
    def filter_unique_supp_reads(self):
        global MIN_UNIQUE_READS
        return len(self.reads) >= MIN_UNIQUE_READS

    # filter based on number of total supporting reads
    def filter_total_supp_reads(self):
        global MIN_TOTAL_READS
        return self.counts >= MIN_TOTAL_READS

    def filter_vaf(self):
        global MIN_VAF
        return self.vaf >= MIN_VAF


class ITD(Insert):
    def __init__(self, insert, tandem2_start, offset):
        self.seq = insert.seq
        self.length = insert.length
        self.start = insert.start
        self.end = insert.end
        self.trailing = insert.trailing
        self.trailing_end = insert.trailing_end
        self.reads = insert.reads
        self.counts = insert.counts
        self.coverage = insert.coverage
        self.vaf = insert.vaf
        
        self.tandem2_start = tandem2_start
        self.offset = offset

    # update length of trailing ITDs to offset instead (should be max potential ITD length)
    # -> update also insert sequence?
    def fix_trailing_length(self):
        if self.trailing:
            self.length = self.offset
            return self
        else:
            return self
        
    # itds
    # for trailing itds, set length to offset
    # for all itds, set start to tandem2_start
    # and norm start coords to [0,len(REF)[ before printing
    def prep_for_save(self):
        to_save = copy.deepcopy(self)
        to_save = to_save.fix_trailing_length()
        to_save.start = to_save.tandem2_start
        to_save = to_save.norm_start()
        return to_save


# when merging inserts, put them into the same collection
class InsertCollection(object):
    def __init__(self, insert):
        self.inserts = [insert]
        self.rep = copy.deepcopy(insert)
        #self.counts = insert.counts
        #self.reads = list(insert.reads)
    
    def set_representative(self):
        # most abundant insert may represent the collection
        # --> how to handle trailing?
        # --> select more abundant or longer insert?
        # --> are ever trailing and non-trailing inserts merged? -> result is what?
        # --> are ever trailing inserts with different trailing_end merged? 
        # ----> result is what?? (should that be possible?? don't think so...)
        self.rep = copy.deepcopy(self.inserts[[insert.counts for insert in self.inserts].index(max([insert.counts for insert in self.inserts]))])
        # reads and counts must be summed for the representative -> overwrite these (that's why representative needs to be a copy!)
        self.rep.reads = flatten_list([insert.reads for insert in self.inserts])
        self.rep.counts = sum([insert.counts for insert in self.inserts])
        self.rep = self.rep.calc_vaf()
        return self

    def merge(self, insert):
        self.inserts = self.inserts + insert.inserts
        self = self.set_representative()
        return self

    def should_merge(self, that, condition):
        for insert in self.inserts:
            for this in that.inserts:
                if insert.should_merge(this, condition):
                    return True
        return False
            

def flatten_list(list_):
    return [item for sublist in list_ for item in sublist]


def connect_bases(char1, char2):
    if char1 == '-' or char2 == '-':
        return ' '
    if char1 == char2:
        return '|'
    return '.'

def connect_alignment(seq1, seq2):
    return ''.join([connect_bases(char1,char2) for char1,char2 in zip(seq1,seq2)])

# count number of digits to align all lines in pretty alignment printout
# --> (need to know how many spaces to insert)
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
def print_alignment(read, out_dir,
        command='bio.align.globalcs', command_seq='read.seq', command_ref='REF'):
    global COST_GAPOPEN, GAPEXTEND
    al = connect_alignment(read.al_seq, read.al_ref)
    al_len = len(read.al_seq)
    command_score_function = "get_alignment_score"
    width = 50
    pre_width = 20
    post_width = 7
    score_width = 15
    #
    with open(os.path.join(out_dir,read.al_file), 'w') as f:
        f.write('########################################\n')
        f.write('# Program: Biopython\n')
        f.write('# Rundate: {}\n'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%d")))
        f.write('# Commandline: {}(\n'.format(command))
        f.write('#    {},\n'.format(command_seq))
        f.write('#    {},\n'.format(command_ref))
        f.write('#    {},\n'.format(command_score_function))
        f.write('#    {},\n'.format(COST_GAPOPEN))
        f.write('#    {})\n'.format(COST_GAPEXTEND))
        f.write('# Align_format: srspair\n')
        f.write('# Report_file: {}\n'.format(read.al_file))
        f.write('########################################\n')
        f.write('\n')
        f.write('#=======================================\n')
        f.write('#\n')
        f.write('# Aligned_sequences: 2\n')
        f.write('# Sample: {}\n'.format(''.join([x for x in read.al_seq if x != '-'])))
        f.write('# Reference: {}\n'.format(''.join([x for x in read.al_ref if x != '-']).lower()))
        f.write('# Matrix: EDNAFULL\n')
        f.write('# Gap_penalty: {}\n'.format(COST_GAPOPEN))
        f.write('# Extend_penalty: {}\n'.format(COST_GAPEXTEND))
        f.write('#\n')
        f.write('# Length: {}\n'.format(al_len))
        identity = '{}/{} ({}%)\n'.format(
            al.count('|'),
            al_len,
            round(al.count('|') / al_len * 100, 1))
        similarity = '{}/{} ({}%)\n'.format(
            al.count('|') + al.count('.'),
            al_len,
            round((al.count('|') + al.count('.')) / al_len * 100, 1))
        gaps = '{}/{} ({}%)\n'.format(
            al_len - al.count('|') - al.count('.'),
            al_len,
            round((al_len - al.count('|') - al.count('.')) / al_len * 100,1))
        f.write('# Identity:     {}'.format(' ' * (score_width - len(identity)) + identity))
        f.write('# Similarity:   {}'.format(' ' * (score_width - len(similarity)) + similarity))
        f.write('# Gaps:         {}'.format(' ' * (score_width - len(gaps)) + gaps))
        f.write('# Score: {}\n'.format(read.al_score))
        f.write('#\n')
        f.write('#\n')
        f.write('#=======================================\n')
        f.write('\n')
        
        # split alignment strings into per-line chunks for pretty printing
        alignment_chunks = [(read.al_seq[i:i+width],al[i:i+width],read.al_ref[i:i+width]) 
            for i in range(0, al_len, width)]
        seq_coord = 1
        ref_coord = 1
        for s,a,r in alignment_chunks:
            seq_coord = print_alignment_seq(s, seq_coord,pre_width,post_width,f)
            print_alignment_connection(a, pre_width,f)
            ref_coord = print_alignment_seq(r, ref_coord,pre_width,post_width,f)
            f.write('\n')
        
        f.write('\n')
        f.write('#---------------------------------------\n')
        f.write('#---------------------------------------\n')


# callback function for align() to calc alignment score for 2 bases
# insert is masked by 'Z'
# --> return max penalty (min score of -Inf) to prohibit realignment of insert to itself
def get_alignment_score(char1,char2):
    global COST_MATCH, COST_MISMATCH
    # only ever one of the sequences chars are taken from should contain masking letter 'Z'
    # -->  i.e. the read sequence but not the ref
    assert not (char1 == 'Z' and char2 == 'Z')
    if char1 == char2:
        return COST_MATCH
    elif char1 == 'Z' or char2 == 'Z':
        return -np.inf
    else:
        return COST_MISMATCH

# get min score required to pass alignment score filter
def get_min_score(seq1, seq2, min_score):
    global COST_MATCH
    return min(len(seq1),len(seq2)) * COST_MATCH * min_score


def parallelize(function, args, cores):
    with multiprocessing.Pool(cores) as p:
        return p.map(function, args)

def read_fastq(filename):
# read in FASTQ files, init Read instances for each record
    reads = []
    with open(filename,'r') as f:
        line = f.readline()
        while line:
            read_id = line
            read_seq = f.readline().rstrip('\n')
            read_desc = f.readline()
            read_bqs = f.readline().rstrip('\n')
            assert len(read_seq) == len(read_bqs)
            reads.append(Read(seq=read_seq, bqs=read_bqs, sense=1))
            line = f.readline()
    return reads#[0:10000] # REMOVE THIS ###############################################################

# read in wt reference for alignment
def read_reference(filename):
    with open(filename, 'r') as f:
        ref = f.read()
    ref = ref.splitlines()
    assert len(ref) == 1
    return ref[0]

# read in wt reference annotation 
# --> genomic / transcript / protein coordinates and exon/intron information
def read_annotation(filename):
    try:
        return pd.read_csv(filename, sep='\t')
    except IOError as e:
        print("No annotation file given")
        return None

# extract domains with start / stop coords into list of tuples
def get_domains(ANNO):
    domains = []
    domain = start = end = None
    for i,row in ANNO.iterrows():
        if domain and domain == row["region"]:
            end = end + 1
        else:
            if domain:
                domains.append((domain,start,end))
            domain = row["region"]
            start = end = row["amplicon_bp"]
    return domains

# check that insert was realigned in one piece
def integral_insert_realignment(insert_alignment, insert_length):
    insert_idxs = [i for i in range(len(insert_alignment)) if insert_alignment[i] != '-']
    return insert_idxs[-1] - insert_idxs[0] +1 == insert_length

# add read.counts to coverage where read covers ref
# --> coords are 0-based and relative to WT ref
def update_coverage(coverage, read):
    refn = np.array(list(read.al_ref))
    readn = np.array(list(read.al_seq))
    wt_ref_covered_bp = np.where(readn[refn != '-'] != '-')
    wt_ref_covered_range = np.arange(
        np.min(wt_ref_covered_bp),
        np.max(wt_ref_covered_bp) + 1) # DO count last index
    coverage[wt_ref_covered_range] = coverage[wt_ref_covered_range] + read.counts
    return coverage

# get coverage at a certain insert position (for all pos in df["start_col"] 
# --> for VAF calculation
def get_coverage(df, start_col, ref_coverage):
    return [ref_coverage[pos-1] for pos in df[start_col]]

def annotate(df):
    global ANNO
    return pd.merge(df, ANNO, 
        how='left', left_on=['start'], right_on=['amplicon_bp']).drop('amplicon_bp', axis=1)

# read in txt file with known ITD length and VAF
def read_known(filename, dtype=str):
    with open(filename) as f:
        return [dtype(x) for x in f.read().splitlines()]

# extract ITDs of known length from df
def get_known(df,known_length):
    global SAMPLE
    df_found = df.ix[[x in known_length for x in df["length"]]]
    #
    # fill in available data on known ITDs/inserts that were missed (and not present in df)
    missed = [x for x in known_length if x not in list(df_found["length"])]
    df_missed = pd.DataFrame( {
        "length": missed, 
        "sample": [SAMPLE] * len(missed),
        "vaf": [0] * len(missed),
        "counts": [0] * len(missed)})
    #
    # concatenate known_found and known_missed
    df_known = pd.concat([df_found, df_missed])
    df_known[["length","counts"]] = df_known[["length","counts"]].astype(int)
    return df_known


# convert ITD allele ratio (AR) to variant allele frequency (VAF)
# --> AR = V-AF/WT-AF
# --> VAF = V-AF
# --> V-AF + WT-AF = 100 (%)
def ar_to_vaf(ar):
    return ar/(ar + 1) * 100 # * 100 because VAF is in %

# convert VAF to AR
# --> if VAF == 100, AR = 100/0??
def vaf_to_ar(vaf):
    if vaf == 100:
        return -1
    return vaf/(100 - vaf) 

# merge / collapse insert/itd records describing the same mutation
def merge(inserts, condition):
    merged = []
    for insert_collection in inserts:
        was_merged = False
        for minsert_collection in merged[::-1]:
            if minsert_collection.should_merge(insert_collection, condition):
                minsert_collection = minsert_collection.merge(insert_collection)
                was_merged = True
                break
        if not was_merged:
            merged.append(insert_collection)
    return merged

# save list of Insert objects to file
# --> create dict, convert to pandas df, print that
def save_to_file(inserts, filename):
    if inserts:
        global SAMPLE, ANNO, DOMAINS
        dict_ins = {}
        for key in vars(inserts[0]):
            dict_ins[key] = tuple(vars(insert)[key] for insert in [insert.prep_for_save() for insert in inserts])
        
        df_ins =  pd.DataFrame(dict_ins)
        df_ins["sample"] = [SAMPLE] * len(inserts)
        df_ins["ar"] = [vaf_to_ar(insert.vaf) for insert in inserts]
        df_ins["counts_each"] = [[read.counts for read in insert.reads] for insert in inserts]
        df_ins["file"] = [[read.al_file for read in insert.reads] for insert in inserts]
        
        cols = ['sample','length', 'start', 'vaf', 'ar', 'coverage', 'counts', 'trailing', 'seq'] 
        # print counts_each only when they contain fewer than X elements (i.e. unique reads)
        #cols = cols + [col for col in ['counts_each'] if max([len(x) for x in df_ins[col]]) <= 10]
        if ANNO is not None:
            # if annotation file exists, 
            # overwrite with annotated df
            # (same command as above!)
            df_ins = annotate(df_ins)
            df_ins["region"] = [insert.annotate_domains(DOMAINS) for insert in inserts]
            cols = cols + ["region", "chr13_bp", "transcript_bp", "protein_as"]
        cols = cols + ['file']
        df_ins[cols].to_csv(os.path.join(OUT_DIR,filename), index=False, float_format='%.2e', sep='\t') # move this command below if, delete the one before if (write csv once, add region column when possible)

def get_unique_reads(reads):
    tmp = collections.Counter([(read.seq,read.sense) for read in reads])
    unique_reads = list(tmp.keys())
    unique_reads_counts = list(tmp.values())
    assert len(unique_reads) == len(unique_reads_counts)
    assert sum(unique_reads_counts) == len(reads)
    reads = [Read(seq=this_seq, bqs=None, counts=this_count, sense=this_sense)
        for (this_seq,this_sense),this_count in zip(unique_reads, unique_reads_counts)]
    return reads

def filter_alignment_score(reads):
    # FILTER BASED ON ALIGNMENT SCORE (INCL FAILED ALIGNMENTS WITH read.al_score is None!
    global REF
    reads_filtered = [
        read for read in reads
        if read.al_score is not None and read.al_score >= get_min_score(
        read.seq, REF, MIN_SCORE_ALIGNMENTS)]
    print("Filtering {} / {} low quality alignments with a score < {} % of max".format(
        len(reads) - len(reads_filtered), len(reads), MIN_SCORE_ALIGNMENTS *100))
    return reads_filtered


########## MAIN ####################
if __name__ == '__main__':
    # prevent neg nkern/minBQS?
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("fastq1", help="FASTQ file of forward reads (REQUIRED)")
    parser.add_argument("fastq2", help="FASTQ file of reverse reads (REQUIRED)")
    parser.add_argument("sampleID", help="sample ID used as output folder prefix (REQUIRED)")
    parser.add_argument("minBQS", help="minimum average base quality score (BQS) required by each read (default 30)", type=int, default=30, nargs='?')
    parser.add_argument("-reference", help="WT amplicon sequence as reference for read alignment (default /NGS/known_sites/hg19/flt3-itd_anno/amplicon.txt)", default="/NGS/known_sites/hg19/flt3-itd_anno/amplicon.txt", type=str)
    parser.add_argument("-anno", help="WT amplicon sequence annotation (default /NGS/known_sites/hg19/flt3-itd_anno/amplicon_kayser.tsv)", default="/NGS/known_sites/hg19/flt3-itd_anno/amplicon_kayser.tsv", type=str)
    parser.add_argument("-technology", help="Sequencing technology used, options are '454' or 'Illumina' (default)", default="Illumina", type=str)
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
    REF_FILE = cmd_args.reference
    ANNO_FILE = cmd_args.anno
    TECH = cmd_args.technology
    NKERN = cmd_args.nkern
    KNOWN_LENGTH_FILE = cmd_args.known_length
    KNOWN_VAF_FILE = cmd_args.known_vaf
    KNOWN_AR_FILE = cmd_args.known_ar
    OUT_DIR = '_'.join([SAMPLE,'minBQS', str(MIN_BQS)])

    COST_MATCH = cmd_args.match
    COST_MISMATCH = cmd_args.mismatch
    COST_GAPOPEN = cmd_args.gap_open
    COST_GAPEXTEND = cmd_args.gap_extend
    MIN_SCORE_INSERTS = cmd_args.minscore_inserts
    MIN_SCORE_ALIGNMENTS = cmd_args.minscore_alignments

    MIN_READ_COPIES = cmd_args.filter_reads
    MIN_TOTAL_READS = cmd_args.filter_ins_total_reads
    MIN_UNIQUE_READS = cmd_args.filter_ins_unique_reads
    MIN_VAF = cmd_args.filter_ins_vaf


    print("==== PROCESSING SAMPLE {} ====".format(SAMPLE))

    ANNO = read_annotation(ANNO_FILE)
    DOMAINS = get_domains(ANNO)
    REF = read_reference(REF_FILE).upper()

    ## CREATE OUTPUT FOLDER
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

    print("-- Reading FASTQ files --")
    start_time = timeit.default_timer()
    reads = read_fastq(R1)

    # IF IT EXISTS:
    # --> reverse-complement R2 reads so that all reads can be aligned to the same reference
    if 'R2' in locals():# --> check how this works with argparse
        reads_rev = read_fastq(R2)
        reads_rev_rev = parallelize(Read.reverse_complement,reads_rev,NKERN)
        reads = reads + reads_rev_rev
    print("Reading FASTQ files took {} s".format(timeit.default_timer() - start_time))
    print("Number of total reads: {}".format(len(reads)))

    ## TRIM trailing AMBIGUOUS 'N's
    reads = parallelize(Read.trim_n, reads, NKERN)

    ## FILTER ON BQS
    if MIN_BQS > 0:
        reads = [x for x in parallelize(Read.filter_bqs, reads, NKERN) if x is not None]
    print("Number of total reads with mean BQS >= {}: {}".format(MIN_BQS,len(reads)))

    # get unique reads and counts thereof
    reads = get_unique_reads(reads)
    print("Number of unique reads with mean BQS >= {}: {}".format(MIN_BQS,len(reads)))

    # FILTER UNIQUE READS
    # --> keep only those that exist at least twice
    # --> assumption: if it's not reproducible, it's not (true and clinically relevant)
    if MIN_READ_COPIES == 1:
        print("Turned OFF unique reads filter!")
    else:
        reads = [read for read in reads if read.counts >= MIN_READ_COPIES ]
        print("Number of unique reads with at least {} copies: {}".format(MIN_READ_COPIES,len(reads)))
    print("Total reads remaining for analysis: {}".format(sum([read.counts for read in reads])))

    ## ALIGN TO REF
    print("\n-- Aligning to Reference --")
    start_time = timeit.default_timer()
    reads = parallelize(Read.align, reads, NKERN)
    print("Alignment took {} s".format(timeit.default_timer() - start_time))

    # FILTER BASED ON ALIGNMENT SCORE (INCL FAILED ALIGNMENTS WITH read.al_score is None!
    reads = filter_alignment_score(reads)

    # FILTER BASED ON MISALIGNED PRIMERS 
    # --> require that primers (26 bp forward / 23 bp reverse) are always aligned with max 3 gaps
    # --> only works for this specific MRD project! --> 454 has different and multiple primers (2 PCRs!)
    if TECH == "Illumina":
        rev_primer = 'GGTTGCCGTCAAAATGCTGAAAG'
        fwrd_primer = 'GCAATTTAGGTATGAAAGCCAGCTAC'
        primers_filtered = [read for read in reads if ( 
            (read.sense == -1 
            and read.al_seq.count('-', read.al_ref.find(rev_primer), read.al_ref.find(rev_primer) + len(rev_primer)) <= 0
            ) or 
                (read.sense == 1 
                and read.al_seq.count('-', read.al_ref.find(fwrd_primer), read.al_ref.find(fwrd_primer) + len(fwrd_primer)) <= 0))]
        primer_fail = [read for read in reads if not ( ### keep this for initial testing only!!!
            (read.sense == -1 
            and read.al_seq.count('-', read.al_ref.find(rev_primer), read.al_ref.find(rev_primer) + len(rev_primer)) <= 0
            ) or 
                (read.sense == 1 
                and read.al_seq.count('-', read.al_ref.find(fwrd_primer), read.al_ref.find(fwrd_primer) + len(fwrd_primer)) <= 0))]
        print("Filtering {} / {} alignments with more than 3 unaligned primer bases".format(
            len(reads) - len(primers_filtered), len(reads)))
        reads = primers_filtered

    # FINAL STATS
    print("Total reads remaining for analysis: {}".format(sum([read.counts for read in reads])))

    # PRINT PASSING ALIGNMENTS
    # create output file directory for alignments print-outs
    needle_dir = os.path.join(OUT_DIR,'out_needle')
    if not os.path.exists(needle_dir):
        os.makedirs(needle_dir)

    for i,read in enumerate(reads):
        reads[i].al_index = i
        reads[i].al_file = 'needle_{}.txt'.format(i)
        print_alignment(reads[i], needle_dir,
            command='bio.align.globalcs', command_seq='read.seq', command_ref='REF')

    if not reads:
        print("\nNO READS TO PROCESS!")
        quit()


    #######################################
    print("\n-- Looking for insertions & ITDs --")

    # COLLECT INSERTS, CALC COVERAGE
    inserts = []
    ref_wt = list(REF)
    ref_coverage = np.zeros(len(REF))

    start_time = timeit.default_timer()
    for read in reads:
        readn = np.array(list(read.al_seq))
        refn = np.array(list(read.al_ref))
        assert(len(readn) == len(refn))
        
        ref_coverage = update_coverage(ref_coverage, read)
        
        # if read contains insert
        insert_idxs_all = np.where(refn == '-')[0]
        if len(insert_idxs_all) > 0:
            # two gaps should never align at the same pos!
            assert('-' not in readn[insert_idxs_all])
            
            # get indeces of individual inserts
            insert_idxs_list = []
            insert_idxs= []
            i_prev = None
            for i_this in insert_idxs_all:
                #start saving first/continue saving next insert index
                if i_prev is None or i_prev == i_this -1:
                    insert_idxs.append(i_this)
                    i_prev = i_this
                #save current insert_idxs and open up a new one for the next insert
                else:
                    insert_idxs_list.append(insert_idxs)
                    insert_idxs = [i_this]
                    i_prev = i_this
            # save last insert as well
            insert_idxs_list.append(insert_idxs)
            assert np.all(np.concatenate(insert_idxs_list) == insert_idxs_all)
           
            # analyze largest insert per read only 
            # --> assume no two true inserts occur within the same read
            # --> assume the smaller one is more likely to be the false positive (not always true though!)
            #insert_idxs_list = [insert_idxs_list[[len(ins) for ins in insert_idxs_list].index(max([len(ins) for ins in insert_idxs_list]))]]
            # --> could try to discard all reads with multiple inserts!  (might be more accurate) 
            # --> would have to discard all mini-inserts (< 6 bp) first to allow alignment/minor sequencing errors
            for insert_idxs in insert_idxs_list:
                if len(insert_idxs) >= 6 and "N" not in readn[insert_idxs]:
                    insert_start = insert_idxs[0]
                    insert_end = insert_idxs[-1]
                    insert = Insert(
                        seq=read.al_seq[insert_start:insert_end+1], 
                        start=insert_start, 
                        end=insert_end, 
                        reads=[read], 
                        counts=read.counts)
                    assert insert.length == len(insert_idxs)
                    
                    if all(readn[0:insert.start] == '-'):
                        insert.trailing_end = 5
                    elif all(readn[insert.end+1:] == '-'):
                        insert.trailing_end = 3
                    else:
                        insert.trailing_end = 0
                    # inserts are considered trailing when it is unclear whether they were covered completely or not
                    # --> because primers guarantee that amplicon starts and ends with WT ref bases,
                    #       forward reads cannot have 5' and reverse reads cannot have 3' trailing insertions
                    #       (trailing_end would be set respectively but as insertions will in fact be fully contained
                    #       within the insert, trailing will be False nontheless)
                    # should I discard reads where the primer is not mapped? See lbseq:/media/data/tabl/laura*/mail/primer_unmapped.txt
                    insert.trailing = (read.sense == 1 and insert.trailing_end == 3) or (read.sense == -1 and insert.trailing_end == 5)
                    
                    if insert.trailing or insert.length % 3 == 0:
                        # change insert.start coord
                        #   from: 1st insert/gap bp in read-ref alignment
                        #   to: preceding bp in WT ref sequence
                        #   --> -sum(preceding gaps) -1
                        insert.start = insert.start - sum(refn[0:insert.start] == '-') -1
                        # distinguish insertions starting before WT base 0 (5' insertion) and those starting right after base 0 (0 would then be preceding WT base)
                        #   --> later negative/ 5' coords are reset to 0 to hide counterintuitive coords (should I set them to -1???)
                        if insert.start == 0:
                            insert.start = -insert.length
                        insert.end = insert.start + insert.length - 1
                        
                        # having passed all filters, save insert
                        inserts.append(insert)

    print("Collecting inserts took {} s".format(timeit.default_timer() - start_time))
    print("{} insertions were found".format(len(inserts)))


    start_time = timeit.default_timer()
    # ref_coverage was of type np.float because initialized as np.zeros 
    # --> (which was necessary to add read counts as I did) 
    # --> could change that by using a generator!
    ref_coverage = ref_coverage.astype(int) 

    # add coverage to inserts # -> delete timer, move up to insert block above
    for insert in inserts:
        # add coverage 
        # --> be sure to normalize start coord to [0,len(REF)[ first
        # --> negative start (-> 5' trailing_end) will result in
        #     coverage = ref_coverage[-X] which will silently report incorrect coverage!!
        insert.coverage = ref_coverage[copy.deepcopy(insert).norm_start().start]
        insert = insert.calc_vaf()
    print("Annotating coverage took {} s".format(timeit.default_timer() - start_time))


    # CHECK WHETHER INSERTIONS ARE ITDs -> place in method later and apply at each step of filtering?
    # Should I do exact searching before doing insert-ref alignment as I used to in the previous script?
    # --> put this in a method and use parallelize to speed things up!
    # --> (can I also do that for reads above when there are possibly multiple inserts per read?) -> yes: return [inserts found] per itd, remove None, flatten list
    itds = []
    start_time = timeit.default_timer()
    for insert in inserts:
        min_score = get_min_score(insert.seq, REF, MIN_SCORE_ALIGNMENTS)
        
        # arguments: seq1, seq2, match-score, mismatch-score, gapopen-score, gapextend-score 
        # output: list of optimal alignments, each a list of seq1, seq2, score, start-idx, end-idx
        alignments = bio.align.localcs(insert.seq, REF, get_alignment_score, COST_GAPOPEN, COST_GAPEXTEND)
        
        # filter alignments where insert cannot be realigned in one piece
        # --> is it possible to obtain alignments where inserts are in one piece and in multiple pieces but with the same al_score??
        # --> I don't think so... --> test below, see if this ever happens (otherwise remove following "if" block and alignments = reassignment)
        if not (not any([al for al in alignments if integral_insert_realignment(al[0],insert.length)]) or all([al for al in alignments if integral_insert_realignment(al[0],insert.length)])):
            print("Inserts broken and whole with same alignment score?!?!?")
            print(alignments)
            stop
        alignments = [al for al in alignments if integral_insert_realignment(al[0],insert.length)]
        if not alignments:
            alignment_score = -1
        else:
            # how often is there more than one alignment? 
            # --> (more than 1 => alignment is ambiguous) 
            # --> Can I choose a smart one somehow? Otherwise return only one in the first place...
            # ---> only if there can never be an integral and non-integral alignment!! 
            alignment = alignments[0]  
            alignment_score, alignment_start, alignment_end = alignment[2:5]
            #print(bio.format_alignment(*alignment))
        if alignment_score >= min_score:
            offset = abs(alignment_start - insert.start)
            # offset = 1 for adjacent insert-tandem2
            # offset = insert.length-1 for adjacent tandem2-insert
            # --> (for tandem2-insert: offset = abs((insert_start - insert.length +1) - insert_start))
            if insert.trailing or (offset == 1 or offset == insert.length - 1):
                # if by chance insert is completely contained within read in spite of it being trailing
                # (i.e. insert and tandem have the same length & are adjacent)
                # --> revert trailing to be able to apply more stringent filters of non-trailing inserts
                if insert.trailing and (offset == 1 or offset == insert.length - 1):
                    insert.trailing = False
                    print("UNTRAILED: {}".format(vars(read)))
                    if insert.length % 3 != 0:
                        print("BUT NOT IN FRAME!!!")
                        # remove from inserts list (see below loop)
                itd = ITD(
                    insert, 
                    offset=offset, 
                    tandem2_start=alignment_start)
                itds.append(itd)

    # in case any out-of-frame insert was untrailed: remove it from list of inserts
    inserts[:] = [insert for insert in inserts if insert.trailing or insert.length % 3 == 0] 
    inserts = sorted(inserts, key=Insert.get_seq)
    itds = sorted(itds, key=Insert.get_seq)
    print("Collecting ITDs took {} s".format(timeit.default_timer() - start_time))
    print("{} ITDs were found".format(len(itds)))


    ########################################
    # MERGE INSERTS
    print("\n-- Merging results --")

    merge_dic = {"insertions": inserts, "itds": itds}
    all_merged = {}
    for inserts_type,inserts_ in merge_dic.items():
        all_merged[inserts_type] = []
        suffix = ""
        # turn Insert objects into InsertCollection to keep merging methods simple and not have to distinguish between the two
        to_merge = [InsertCollection(insert) for insert in inserts_]
        for condition,abrev in [
                ("is_same",""), 
                ("is_similar","similar"), 
                ("is_close","close"), 
                ("is_same_trailing","trailing")]:
            to_merge = merge(to_merge, condition)
            all_merged[inserts_type].append(to_merge)
            print("{} {} remain after merging".format(len(to_merge), inserts_type))
            suffix = suffix + "-" + condition
            save_to_file([insert.rep for insert in to_merge], "flt3_" + inserts_type + "_collapsed" + suffix + ".tsv")


    # save as list of Insert(s) to process further as before (vs continuing with list of InsertCollection(s)
    final_merged = {}
    for inserts_type in merge_dic:
        final_merged[inserts_type] = [insert.rep for insert in all_merged[inserts_type][-1]]
    print("Merging took {} s".format(timeit.default_timer() - start_time))


    ########################################
    # FILTER INSERTS
    print("\n-- Filtering --")

    final_filtered = {}
    filter_dic = {
        "number of unique supporting reads": Insert.filter_unique_supp_reads, 
        "number of total supporting reads": Insert.filter_total_supp_reads,
        "vaf": Insert.filter_vaf}
    start_time = timeit.default_timer()
    for inserts_type,inserts_ in final_merged.items():
        filtered = copy.deepcopy(inserts_)
        for filter_type,filter_ in filter_dic.items():
            passed = [filter_(insert) for insert in filtered]
            filtered = [insert for (insert,pass_) in zip(filtered, passed) if pass_]
            print("Filtered {} / {} {} based on the {}".format(
                len(passed) - sum(passed), len(passed), inserts_type, filter_type))
        print("{} {} remain after filtering!".format(len(filtered), inserts_type))
        save_to_file(filtered, "flt3_" + inserts_type + "_collapsed" + suffix + "_hc.tsv")
        final_filtered[inserts_type] = filtered
    print("Filtering took {} s".format(timeit.default_timer() - start_time))



    # compare with get_itd to make sure I didn't forget anything

    ########################################
    # GET KNOWN ITDs
    # --> if length known is given, extract relevant ITDs
    # --> if vaf known is given, compare numbers? (vaf is always known if length is known? is there a better way to supply this together?)
    if False and KNOWN_LENGTH_FILE is not None:
        known_length = read_known(KNOWN_LENGTH_FILE,int)
        df_itds_known = get_known(fix_trailing_length(df_itds_collapsed), known_length)
        df_itds_known[['sample','length','vaf','coverage','counts','tandem2_start','seq']].to_csv(os.path.join(OUT_DIR,"flt3_itds_collapsed-similar-close-trailing_hc_known.tsv"), index=False, float_format='%.2e', sep='\t', na_rep='NA')
        #
        df_ins_known = get_known(df_ins_collapsed, known_length)
        df_ins_known[['sample','length','vaf','coverage','counts','start','seq']].to_csv(os.path.join(OUT_DIR,"flt3_ins_collapsed-similar-close-trailing_hc_known.tsv"), index=False, float_format='%.2e', sep='\t', na_rep='NA')
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
        df_itds_known_collapsed = collapse(df_itds_known,keep=["sample"],add=["counts","vaf"],append=["length","tandem2_start","coverage"])
        df_itds_known_collapsed["vaf_genescan"] = known_vaf
        df_itds_known_collapsed[['sample','length','vaf','vaf_genescan','vaf_each','tandem2_start','coverage','counts_each']].to_csv(os.path.join(OUT_DIR,"flt3_itds_collapsed-similar-close-trailing_hc_known_collapsed.tsv"), index=False, float_format='%.2e', sep='\t', na_rep='NA')
        #
        df_ins_known_collapsed = collapse(df_ins_known,keep=["sample"],add=["counts","vaf"],append=["length","start","coverage"])
        df_ins_known_collapsed["vaf_genescan"] = known_vaf
        df_ins_known_collapsed[['sample','length','vaf','vaf_genescan','vaf_each','start','coverage','counts_each']].to_csv(os.path.join(OUT_DIR,"flt3_ins_collapsed-similar-close-trailing_hc_known_collapsed.tsv"), index=False, float_format='%.2e', sep='\t', na_rep='NA')


