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


class Read(object):
    """
    Sequencing read.
    """ 
    def __init__(
                self, 
                seq, 
                index=None,
                sense=1, 
                bqs=None, 
                index_bqs=None, 
                counts=1,
                al_score=None, 
                al_seq=None, 
                al_ref=None,
                al_file=None):
        """
        Initialize Read.

        Args:
            seq (str): Base pair sequence of the read.
            index (int): Read index to identify paired reads
                for paired-end sequencing 
                --> will have the same index
            sense (int): 1 for forward reads, -1 for reverse.
            bqs (str): Base quality scores of the read.
            index_bqs (str): Base quality scores of the 
                corresponding index read.
            counts (int): Total number of reads with these
                exact attributes.
            al_score (int): Score of the read-to-reference
                alignment.
            al_seq (str): Read sequence aligned to the reference.
            al_ref (str): Reference sequence aligned to the read.
            al_file (str): Name of the file containing the 
                read-to-reference alignment.
        """
        self.seq = seq
        self.index = index
        self.bqs = bqs
        self.index_bqs = index_bqs
        self.length = len(seq)
        self.counts = counts
        self.sense = sense
        self.al_score = al_score
        self.al_seq = al_seq
        self.al_ref = al_ref
        self.al_file = al_file
        assert self.counts > 0
        if self.bqs is not None:
            assert len(self.seq) == len(self.bqs)
    
    def print(self):
        """
        Pretty print Read.
        """
        pprint.pprint(vars(self))
    
    def reverse_complement(self):
        """
        Reverse complement a given Read.

        Reverse complement the Read's sequence, 
        reverse its BQS and invert its sense.

        Returns:
            Reverse complemented Read.
        """
        rev = copy.deepcopy(self)
        rev.seq = self.seq.translate(str.maketrans('ATCGatcg','TAGCtagc'))[::-1]
        if self.bqs:
            rev.bqs = self.bqs[::-1]
        if not self.sense:
            rev.sense = -1
        else:
            rev.sense = self.sense * -1
        return rev
    
    def trim_n(self):
        """
        Trim trailing N's at Read's ends.

        Returns:
            Trimmed Read.
        """
        while self.seq.endswith('N'):
            self.seq = self.seq[:-1]
            if self.bqs:
                self.bqs = self.bqs[:-1]
        while self.seq.startswith('N'):
            self.seq = self.seq[1:]
            if self.bqs:
                self.bqs = self.bqs[1:]
        assert len(self.seq) == len(self.bqs)
        self.length = len(self.seq)
        if self.length < MIN_READ_LENGTH:
            return None
        return self
    

    def filter_bqs(self):
        """
        Filter Read based on its average BQS.

        Returns:
            Read when it passes the average BQS filter,
            None otherwise. 
            
            When no BQS is set, also return the Read.
        """
        if not self.bqs or average_bqs(self.bqs) >= MIN_BQS:
            if not self.index_bqs or average_bqs(self.index_bqs) >= MIN_BQS:
                return self
        return None
    
    def align(self):
        """
        Align a Read to the WT reference.
        
        For 454 data, I do not know which Read is forward, which is
        reverse - try both and keep the reverse Read if that succeeds
        instead of the forward.
        
        Consider:
            Is there a smart way to handle multiple alignments?

        Returns:
            Aligned Read. May be reversed for 454.
        """
        alignment = bio.align.globalcs(self.seq, REF, get_alignment_score,
            COST_GAPOPEN, COST_GAPEXTEND, penalize_end_gaps=False)
        if alignment:
            self.al_seq, self.al_ref, self.al_score = alignment[-1][0:3]

        if TECH == '454':
            rev = self.reverse_complement()
            rev_alignment = bio.align.globalcs(rev.seq, REF, get_alignment_score,
                COST_GAPOPEN, COST_GAPEXTEND, penalize_end_gaps=False)
            if (rev_alignment and
                    (self.al_score is None or rev_alignment[-1][2] > self.al_score)):
                rev.al_seq, rev.al_ref, rev.al_score = rev_alignment[-1][0:3]
                return rev
        return self

    def get_reference_range_covered(read):
        """
        Get most 3' and 5' reference base covered by read.
        Coordinates are 0-based and relative to the WT reference.

        Args:
            read (Read): Read for which to get range of WT reference spanned.

        Returns:
            Updated read.
        """
        refn = np.array(list(read.al_ref))
        readn = np.array(list(read.al_seq))
        ref_covered_bp = np.where(readn[refn != '-'] != '-')

        read.ref_span = [np.min(ref_covered_bp), np.max(ref_covered_bp)]
        return read


class Insert(object):
    """
    Insertion.
    """
    def __init__(
                self, 
                seq, 
                start, 
                end, 
                counts,
                trailing=None, 
                trailing_end=None, 
                reads=None, 
                coverage=None, 
                vaf=None):
        """
        Initialize Insert.

        Args:
            seq (str): Base pair sequence of the Insert.
            start (int): Start coordinate of the Insert. Relative to
                the WT reference, 0-based, refers to the last WT base
                before the insert.
            end (int): End coordinate. Should be start + length -1.
            counts (int): Total number of supporting reads.
            trailing (bool): True when insert is not fully spanned but
                instead occurs at the very beginning or end of supporting
                reads. In such a case, the insert's true length is
                unknown and only estimated using the position of the
                second tandem.
            trailing_end (int): 5 or 3, depending on the end at which
                the insert trails / is not encompassed by bases aligned
                to the reference.
            reads ([Read]): List of supporting reads.
            coverage (int): Total number of reads spanning start.
            vaf (dc.Decimal): Fraction of insert supporting reads in
                all spanning reads, i.e. counts / coverage.
        """
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
    
    def get_seq(self):
        """"
        Return Insert's sequence.
        Implemented to sort list of Inserts by their sequence.
        """
        return self.seq

    def calc_vaf(self):
        """
        Calculate and set Insert's variant allele frequency (VAF).

        Returns:
            Insert.
        """
        self.vaf = dc.Decimal(self.counts) / self.coverage * 100
        assert self.vaf >= 0 and self.vaf <= 100
        return self

    def norm_start(self):
        """
        Normalize Insert's start coordinate to [0, len(REF)[.

        Returns:
            Insert.
        """
        self.start = min(max(0,self.start), len(REF)-1)
        return self

    def prep_for_save(self):
        """
        Prepare Insert for saving to file by normalizing start
        coordinates to [0, len(REF)[.
        
        Returns:
            Prepared Insert. 
        """
        to_save = copy.deepcopy(self)
        to_save = to_save.norm_start()
        return to_save

    def annotate_domains(self, domains):
        """
        Get annotated domains of the Insert's sequence.

        Args:
            domains: List of tuples (domain_name, start_coord, end_coord)
            where coordinates refer to amplicon bp and the list contains
            all of its annotated domains.

        Returns:
            Subset of domains, containing only those affected by the Insert.
        """
        annotated = []
        for domain,start,end in domains:
            if self.start <= end and self.end >= start:
                annotated.append(domain)
        return annotated

    
    def print(self):
        """
        Pretty print Insert.
        Leave out supporting reads as they will clutter the screen.
        """
        pprint.pprint({key: vars(self)[key] for key in vars(self).keys() if not key == 'reads'})

    def is_close_to(self, that):
        """
        Test whether two Inserts are close and 
        located within one insert length of each other.
        
        Args:
            that: Insert to compare to.

        Returns:
            True when they are close, False otherwise.
        """
        if hasattr(self, 'tandem2_start'):
            return abs(self.tandem2_start - that.tandem2_start) <= self.length
        return abs(self.start - that.start) <= 2 * self.length
    
    def is_similar_to(self, that):
        """
        Test whether two Inserts' sequences are similar.

        Args:
            that: Insert to compare to.

        Returns:
            True when they are similar, False otherwise.
        """
        min_score = get_min_score(self.seq, that.seq, MIN_SCORE_INSERTS)
        al_score = bio.align.globalcs(self.seq, that.seq, get_alignment_score, COST_GAPOPEN, COST_GAPEXTEND, one_alignment_only=True, score_only=True, penalize_end_gaps=False)
        if al_score >= min_score:
            return True
        return False
    
    def should_merge(self, that, condition):
        """
        Test whether two Inserts should be merged.

        Args:
            that (Insert): Insert to potentially merge with.
            condition (str): Encodes condition to be met for merging,
                one of 'is-same', 'is-similar', 'is-close',
                'is-same-trailing' or 'any'.

        Returns:
            True when Inserts should be merged, False otherwise.
        """
        if condition == 'is-same':
            return self.seq == that.seq and self.start == that.start
        elif condition == 'is-similar':
            return self.length == that.length and self.start == that.start and self.is_similar_to(that)
        elif condition == 'is-close':
            return self.length == that.length and self.is_similar_to(that) and self.is_close_to(that)
        elif condition == 'is-same_trailing':
            return self.trailing and that.trailing and self.trailing_end == that.trailing_end and self.is_similar_to(that)
        elif condition == 'any':
            return ((self.length == that.length and self.is_close_to(that)) or (self.trailing and that.trailing and self.trailing_end == that.trailing_end)) and self.is_similar_to(that)
        else:
            print("\n---Undefined merging condition!---\n") 
            assert False
    
    
    def filter_unique_supp_reads(self):
        """
        Test whether Insert is supported by a given number of
        distinct supporting reads.

        Returns:
            True when that is the case, False otherwise. 
        """
        return len(self.reads) >= MIN_UNIQUE_READS

    def filter_total_supp_reads(self):
        """
        Test whether Insert is supported by a given minimum number 
        of total supporting reads.

        Returns:
            True when that is the case, False otherwise.
        """
        return self.counts >= MIN_TOTAL_READS

    def filter_vaf(self):
        """
        Test whether Insert has at least a given minimum VAF.
        
        Returns:
            True when that is the case, False otherwise.
        """
        return self.vaf >= MIN_VAF


class ITD(Insert):
    """
    Internal Tandem Duplication.
    """
    def __init__(self, insert, tandem2_start, offset, external_bp):
        """
        Initialize ITD.

        Args:
            insert (Insert): Insert that classifies as an ITD.
            tandem2_start (int): Start coordinate of the second tandem.
                Like start, this is 0-based and relative to the
                reference but refers to the actual first WT base of
                the tandem. For adjacent insert and tandem
                constellations therefore: | start - tandem2_start | = length
            offset (int): Number of bases between insert and tandem2
                start coordinates. For adjacent constellations therefore
                offset = length.
            external_bp (int): Number of external / non-WT bases inserted
                before and after the duplicated WT sequence. These are 
                included in insert.seq and length and must not be added on top!
        """
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
        self.external_bp = external_bp

    def fix_trailing_length(self):
        """
        For trailing ITDs, change length to offset.
        This should be the maximum potential ITD length.
        
        Consider: Should the insert sequence be updated as well?

        Returns:
           ITD, changed if trailing, unchanged if not. 
        """
        if self.trailing:
            self.length = self.offset
            return self
        else:
            return self
        
    def prep_for_save(self):
        """
        Prepare ITD for saving to file.

        Set the insert start coordinate to that of the second tandem.
        Normalize start coordinate to [0, len(REF)[ and 
        adjust end coordinate accordingly.
        For trailing ITDs, change the length to the offset, i.e. 
        distance between insert and second tandem. 
        
        Return:
            Prepared ITD.
        """
        to_save = copy.deepcopy(self)
        to_save = to_save.fix_trailing_length()
        to_save.start = to_save.tandem2_start
        to_save.end = to_save.start + to_save.length - 1
        assert to_save.end <= len(REF)
        to_save = to_save.norm_start()
        return to_save


class InsertCollection(object):
    """
    Bundle merged Inserts.
    """
    def __init__(self, insert):
        """
        Initialize InsertCollection.

        Args:
            inserts: List of inserts.
            rep (Insert): Single Insert to represent the entire
                collection. After each merge, will be set to the
                most abundant one.
        """
        self.inserts = [insert]
        self.rep = copy.deepcopy(insert)
    
    def set_representative(self):
        """
        Set InsertCollection's most abundant Insert as its representative.

        Consider: 
            How to handle trailing inserts? Should longest or most
                abundant insert be selected?
            Are ever trailing and non-trailing inserts merged? What
                would the result be?
            Are ever trailing inserts with different trailing_end merged?
                What would the result be? Should this even be possible?

        Returns:
            InsertCollection with set representative. 
        """
        self.rep = copy.deepcopy(self.inserts[[insert.counts for insert in self.inserts].index(max([insert.counts for insert in self.inserts]))])
        # reads and counts must be summed for the representative -> overwrite these (that's why representative needs to be a copy!)
        self.rep.reads = flatten_list([insert.reads for insert in self.inserts])
        self.rep.counts = len(set(flatten_list([read.index for insert in self.inserts for read in insert.reads])))
        self.rep = self.rep.calc_vaf()
        return self

    def merge(self, insert):
        """
        Merge two InsertCollections.

        Args:
            insert: InsertCollection to merge with.

        Returns:
            Merged InsertCollection.
        """
        self.inserts = self.inserts + insert.inserts
        self = self.set_representative()
        return self

    def should_merge(self, that, condition):
        """
        Test whether two InsertCollections should be merged.

        Args:
            that (InsertCollection): InsertCollection to potentially
                merge with.
            condition (str): String encoding condition to be fullfilled
                for merging.
    
        Returns:
            True when condition is met and InsertCollections should be
            merged, False otherwise.
        """
        for insert in self.inserts:
            for this in that.inserts:
                if insert.should_merge(this, condition):
                    return True
        return False
            

def flatten_list(list_):
    """
    Turn list of lists into list of sublists' elements.

    Args:
        list_ (list): List to flatten.

    Returns:
        Flattened list.
    """
    return [item for sublist in list_ for item in sublist]

def average_bqs(bqs):
    """
    Calculate the mean BQS of a given string of quality scores.
    Assumes BQS are in Sanger format, encoded as Phred +33.

    Args:
        bqs (str): String of base quality scores.

    Returns:
        Mean BQS of that string.
    """
    return sum([ord(x) - 33 for x in bqs]) / len(bqs)

def connect_bases(char1, char2):
    """
    For two aligned bases, get connecting symbol.
    Bases on whether bases match, mismatch or are part of a gap.

    Args:
        char1, char2 (str): Two aligned bases.
    """
    if char1 == '-' or char2 == '-':
        return ' '
    if char1 == char2:
        return '|'
    return '.'

def connect_alignment(seq1, seq2):
    """
    For two aligned sequences, get connecting symbols.
    Based on whether bases match, mismatch or are part of a gap.

    Args:
        seq1, seq2 (str): Two aligned sequences.
    """
    return ''.join([connect_bases(char1,char2) for char1,char2 in zip(seq1,seq2)])

def get_number_of_digits(number):
    """
    Count the number of digits in a given number. 
    
    Use this, to print that many less spaces in front of each 
    part of the well-formatted read-to-reference alignment.

    Args:
        number (int): Number whose digits to count.

    Returns:
        The number of digits.
    """
    if number == 0:
        return 1
    return int(np.log10(number)) +1

def print_alignment_connection(connection, pre_width, f):
    """
    Print the symbols connecting aligned read and reference, encoding
    match, mismatch and gap at each bp.

    Args: 
        connection (str): String of connecting symbols.
        pre_width (int): Number of spaces printed before the alignment to keep everything formatted.
        f (file_object): Output file to contain the alignment. 
    """
    f.write(' ' * (pre_width +2))
    f.write(connection)
    f.write('\n')

def print_alignment_seq(seq, seq_coord, pre_width, post_width, f):
    """
    Print part of an alignment.

    Args:
        seq (str): The part of the read's sequence to be printed.
        seq_coord (int): Start coordinate of the part of the alignment printed.
        pre_width (int): Number of spaces printed before the alignment to keep everything formatted.
        post_width (int): Number of spaces printed after it.
        f (file_object): Output file to contain the alignment as a whole.

    Return:
       Start coordinate of the next part of this alignment to be printed.
    """
    f.write(' ' * (pre_width - get_number_of_digits(seq_coord) +1))
    f.write(str(seq_coord) + ' ')
    f.write(seq)
    seq_coord = seq_coord + len(seq) - seq.count('-') -1
    f.write(' ' * (post_width - get_number_of_digits(seq_coord)))
    f.write(str(seq_coord) + '\n')
    return seq_coord +1

def print_alignment(read, out_dir):
    """
    Print read-to-reference alignment in a nice format, inspired by EMBOSS needle output.

    Args:
        read (Read): Read whose alignment will be printed.
        out_dir (str): Name of output directory to save alignment files to.
    """
    al = connect_alignment(read.al_seq, read.al_ref)
    al_len = len(read.al_seq)

    command = 'bio.align.globalcs'
    command_seq = 'read.seq'
    command_ref = 'REF'
    command_score_function = "get_alignment_score"

    width = 50
    pre_width = 20
    post_width = 7
    score_width = 15
    
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
        alignment_chunks = [(
                read.al_seq[i:i+width],
                al[i:i+width],
                read.al_ref[i:i+width])
                for i in range(0, al_len, width)]
        seq_coord = 1
        ref_coord = 1
        for s, a, r in alignment_chunks:
            seq_coord = print_alignment_seq(s, seq_coord,pre_width,post_width,f)
            print_alignment_connection(a, pre_width, f)
            ref_coord = print_alignment_seq(r, ref_coord,pre_width,post_width,f)
            f.write('\n')
        
        f.write('\n')
        f.write('#---------------------------------------\n')
        f.write('#---------------------------------------\n')


def get_alignment_score(char1,char2):
    """
    Calculate the alignment score of two aligned bases.

    When realigning an insert to the WT reference, the actual insert
    is masked by 'Z'. Return the maximum penalty (- np.inf) to probibit
    realignment of the insert to itself.

    Args:
        char1, char2 (str): The two aligned bases.

    Returns:
        Alignment score (float).
    """
    if char1 == char2:
        return COST_MATCH
    elif char1 == 'Z' or char2 == 'Z':
        return -np.inf
    else:
        return COST_MISMATCH

def get_min_score(seq1, seq2, min_score):
    """
    For two sequences, calculate the minimum alignment score required
    to pass the respective alignment filter, which is a percentage of the
    maximum possible alignment score between these sequences - and 
    therefore dependent on the two aligned sequences' length.

    Args:
        seq1, seq2 (str): The two aligned sequences.

        min_score (float): Proportion of the maximum possible alignment score
                that is required to pass the alignment score filter.

    Returns:
        Minimum required alignment score (float).
    """
    return min(len(seq1),len(seq2)) * COST_MATCH * min_score


def parallelize(function, args, cores):
    """
    Parallelize a given function across a given number of cores.

    Args:
        function (function): Function or method to parallelize.

        args (tuple): Tuple of function's arguments.

        cores (int): Number of cores to utilize.

    Returns:
        List of function's outputs.
    """
    with multiprocessing.Pool(cores) as p:
        return p.map(function, args)

def read_fastq(fastq_file):
    """
    Read sequence fastq file and extract sequences and BQS.

    Args:
        fastq_file: Name of the fastq file to read, R1 or R2.        

    Returns:
        List of Read() objects.
    """
    reads = []
    read_index = 0
    try:
        with open(fastq_file,'r') as f:
            line = f.readline()
            while line:
                read_id = line
                read_seq = f.readline().rstrip('\n')
                read_desc = f.readline()
                read_bqs = f.readline().rstrip('\n')
                assert len(read_seq) == len(read_bqs)
                reads.append(Read(seq=read_seq, index=read_index, bqs=read_bqs, sense=1))
                line = f.readline()
                read_index += 1
    # catch missing file or permissions
    except IOError as e:
        print("---\nCould not read fastq file {}!\n---".format(fastq_file))
    return reads


def read_index_bqs(index_file):
    """
    Read index BQS. 

    Args:
        index_file: Name of the index fastq file, I1 or I2.

    Returns:
        List of index BQS.
    """
    try:
        with open(index_file, "r") as f:
            all_lines = f.readlines()
            index_bqs = all_lines[3::4]
            index_bqs = [bqs.rstrip('\n') for bqs in index_bqs]
        return index_bqs
    except IOError as e:
        print("---\nCould not read index file {}!\n---".format(index_file))


def read_reference(filename):
    """
    Read in WT reference sequence.

    Args:
        filename (str): Name of the file to be read.

    Returns: 
        Reference sequence, stripped of trailing newlines.
    """
    with open(filename, 'r') as f:
        ref = f.read()
    ref = ref.splitlines()
    assert len(ref) == 1
    return ref[0]

# add column names!
def read_annotation(filename):
    """
    Read in WT reference annotation file.
    
    For each bp of the WT reference, provides genomic, transcriptomic
    and proteomic coordinate, exon/intron annotation and the respective
    reference bp.

    Args:
        filename (str): Name of the file to be read.

    Returns:
        pd.DataFrame of the annotation.
    """
    try:
        return pd.read_csv(filename, sep='\t')
    except IOError as e:
        print("No annotation file given")
        return None

def get_domains(anno):
    """
    Extract start and stop coordinates of annotated domains.

    Args:
       anno (pd.DataFrame): Annotation supplied by user with one row
                per reference bp and at least two columns: 
                "region" contains the domain name and 
                "amplicon_bp" the annotated coordinates.
     
    Returns:
        List of tuples (domain_name, start_coord, end_coord) where
        coordinates refer to amplicon bp.
    """
    domains = []
    domain = start = end = None
    for i,row in anno.iterrows():
        if domain and domain == row["region"]:
            end = end + 1
        else:
            if domain:
                domains.append((domain,start,end))
            domain = row["region"]
            start = end = row["amplicon_bp"]
    return domains

def integral_insert_realignment(insert_alignment, insert_length):
    """
    Check whether insert realigned without gaps to the reference.

    Inserts are only considered as ITDs if they realign to the reference
    in one piece (to their respective second tandem).

    Args:
        insert_alignment (str): Alignment string of insert relative to
                WT reference, output by bio.align.localcs().
        insert_length:  Length / number of bp of the insert.

    Returns:
        bool, True when alignment contains one or more gaps, False otherwise.
    """
    insert_idxs = [i for i in range(len(insert_alignment)) if insert_alignment[i] != '-']
    return insert_idxs[-1] - insert_idxs[0] + 1 == insert_length


def annotate(df):
    """
    Annotate insertions with chromosomal regions and coordinates.

    Add genomic, transcriptomic and proteomic coordinates and
    genomic region (exonic, intronic, splicing).

    Args:
        df (pd.DataFrame): DataFrame with insertions.

    Returns:
        Annotated df.
    """
    df = pd.merge(df, ANNO,
        how='left', left_on=['end'], right_on=['amplicon_bp']).drop(['amplicon_bp', 'region'], axis=1)
    df = df.rename(columns={"chr13_bp": "end_chr13_bp", "transcript_bp": "end_transcript_bp", "protein_as": "end_protein_as"})
    df = pd.merge(df, ANNO,
        how='left', left_on=['insertion_site'], right_on=['amplicon_bp']).drop(['amplicon_bp', 'region'], axis=1)
    df = df.rename(columns={"chr13_bp": "insertion_site_chr13_bp", "transcript_bp": "insertion_site_transcript_bp", "protein_as": "insertion_site_protein_as"})
    df = pd.merge(df, ANNO,
        how='left', left_on=['start'], right_on=['amplicon_bp']).drop('amplicon_bp', axis=1)
    df = df.rename(columns={"chr13_bp": "start_chr13_bp", "transcript_bp": "start_transcript_bp", "protein_as": "start_protein_as"})
    return df


def ar_to_vaf(ar):
    """
    Convert AR to VAF.

    VAF (variant allele frequency) = V-AF
    AR (allele ratio) = V-AF / WT-AF
    V-AF + WT-AF = 100 (%)

    Args:
        ar (float): AR to convert.
    
    Returns:
        VAF (float)
    
    """
    return ar/(ar + 1) * 100 # * 100 because VAF is in %

def vaf_to_ar(vaf):
    """
    Convert VAF to AR.

    VAF (variant allele frequency) = V-AF
    AR (allele ratio) = V-AF / WT-AF
    V-AF + WT-AF = 100 (%)
    
    Note:
        if VAF == 100:
            AR = -1
            (instead of 100 / 0)

    Args:
        vaf (dc.Decimal): VAF to convert.

    Returns:
        AR (dc.Decimal)
    """
    if vaf == 100:
        return -1
    return vaf/(100 - vaf)

def merge(inserts, condition):
    """
    Merge insertions describing the same mutation.

    Args:
        inserts ([InsertCollection]): List of insertions to merge.
        condition (str): Encodes condition to determine whether
                insertions do describe the same mutation. One of 
                "is-same", "is-similar", "is-close", "is-same_trailing"

    Returns:
        InsertCollection of merged insertions.
    """
    still_need_to_merge = True
    while still_need_to_merge:
        merged = []
        still_need_to_merge = False
        for insert_collection in inserts:
            was_merged = False
            for minsert_collection in merged[::-1]:
                if minsert_collection.should_merge(insert_collection, condition):
                    if not was_merged:
                        minsert_collection = minsert_collection.merge(insert_collection)
                        was_merged = True
                    else:
                        if not still_need_to_merge:
                            print("need another round of merging!")
                        still_need_to_merge = True
                        break
            if not was_merged:
                merged.append(insert_collection)
        inserts = merged
    return merged

def save_to_file(inserts, filename):
    """
    Write insertions detected to TSV file.
    
    Add additional columns with sample ID, actual coordinate of the
    insertion site, allelic ratio (MUT : WT) and counts and alignment 
    filenames of each of the distinct supporting reads. If annotation 
    was supplied, add that as well. 

    Args:
        inserts ([InsertCollection]): List of inserts to save.
        filename (str): Name of the file, will be saved in specified OUT_DIR.
    """
    if inserts:
        dict_ins = {}
        for key in vars(inserts[0]):
            dict_ins[key] = tuple(vars(insert)[key] for insert in [insert.prep_for_save() for insert in inserts])
        
        df_ins =  pd.DataFrame(dict_ins)
        df_ins["sample"] = [SAMPLE] * len(inserts)
        df_ins["insertion_site"] = df_ins["end"] + 3 # insertion site = WT AS after insert --> +3 to make sure the next AS is annotated instead of the current one
        df_ins["ar"] = [vaf_to_ar(insert.vaf) for insert in inserts]
        df_ins["counts_each"] = [[read.counts for read in insert.reads] for insert in inserts]
        df_ins["file"] = [[read.al_file for read in insert.reads] for insert in inserts]
        
        cols = ['sample','length', 'start', 'end', 'vaf', 'ar', 'coverage', 'counts', 'trailing', 'seq']
        if 'external_bp' in df_ins:
            cols = cols + ['external_bp']

        # print counts_each only when they contain fewer than X elements (i.e. unique reads)
        #cols = cols + [col for col in ['counts_each'] if max([len(x) for x in df_ins[col]]) <= 10]
        if ANNO is not None:
            # if annotation file exists,
            # overwrite with annotated df
            # (same command as above!)
            df_ins = annotate(df_ins)
            df_ins["region"] = [insert.annotate_domains(DOMAINS) for insert in inserts]
            cols = cols + ["region", "start_chr13_bp", "start_transcript_bp", "start_protein_as", "end_chr13_bp", "end_transcript_bp", "end_protein_as", "insertion_site_protein_as"]
        cols = cols + ['file']
        df_ins[cols].sort_values(by=['length','start','vaf']).to_csv(os.path.join(OUT_DIR,filename), index=False, float_format='%.2e', sep='\t')

def get_unique_reads(reads):
    """
    Merge reads with identical read sequences.
    
    Create a Read() object and sum up supporting read counts 
    in Read.counts for each Read.seq, keep reads of different
    orientation distinct, i.e. create one Read for each combination
    of Read.seq and Read.sense. 

    This is really slow. Come up with something better!
    Would it help to process forward/reverse reads separately?

    Args:
        reads ([Read]): Reads to merge, may share the same Read.seq.
    
    Returns:
        Merged reads ([Read]), each with a unique Read.seq.
        
    """
    seqs = [read.seq for read in reads]
    unique_seqs, inverse_indices = np.unique(seqs, return_inverse=True)
    nreads = np.array(reads)
    unique_reads = []
    for inverse_index, seq in enumerate(unique_seqs):
        list_reads = nreads[inverse_indices == inverse_index]
        list_reads_index = [read.index for read in list_reads]
        list_reads_sense = set([read.sense for read in list_reads])
        for sense in list_reads_sense:
            unique_reads.append(
                    Read(seq=seq, sense=sense, bqs=None, counts=len(list_reads), index=list_reads_index))
    return unique_reads


def filter_alignment_score(reads):
    """
    Filter reads based on alignment score.

    Keep reads with an alignment score above the specified fraction of
    the max achievable alignment score. Failed alignments' score is None.

    Args:
        reads ([Read]): Reads to filter.
    Returns:
        Passing reads ([Read]).
    """
    reads_filtered = [
        read for read in reads
        if read.al_score is not None and read.al_score >= get_min_score(
        read.seq, REF, MIN_SCORE_ALIGNMENTS)]
    save_stats("Filtering {} / {} low quality alignments with a score < {} % of max".format(
            len(reads) - len(reads_filtered), len(reads), MIN_SCORE_ALIGNMENTS *100), STATS_FILE)
    return reads_filtered


def save_config(cmd_args, filename):
    """
    Write timestamp and commandline arguments to file.

    Args:
        cmd_args (argparse.Namespace): Commandline arguments and values to write.
        filename (str): Name of the file to write to.
    """
    with open(filename, "w") as f:
        f.write("Commandline_argument\tValue\n")
        f.write("Time\t{}\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%d")))
        #for arg in vars(cmd_args):
        # write arguments in alphabetical order
        for arg in sorted(list(vars(cmd_args).keys())):
            f.write("{}\t{}\n".format(arg, vars(cmd_args)[arg]))

def save_stats(stat, filename):
    """
    Write statistics to file.

    Args:
        stat (str): Statistic to save.
        filename (str): Name of the file to write to.
    """
    print(stat)
    with open(filename, "a") as f:
        f.write(stat + "\n")


########## MAIN ####################
if __name__ == '__main__':
    # prevent neg nkern/minBQS?
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("fastq1", help="FASTQ file of forward reads (REQUIRED)")
    parser.add_argument("fastq2", help="FASTQ file of reverse reads (REQUIRED)")
    parser.add_argument("sampleID", help="sample ID used as output folder prefix (REQUIRED)")
    parser.add_argument("minBQS", help="minimum average base quality score (BQS) required by each read (default 30)", type=int, default=30, nargs='?')
    parser.add_argument("-index1", help="FASTQ file of I1 index (may be omitted)", default=None)
    parser.add_argument("-index2", help="FASTQ file of I2 index (may be omitted)", default=None)
    parser.add_argument("-reference", help="WT amplicon sequence as reference for read alignment (default ./anno/amplicon.txt)", default="./anno/amplicon.txt", type=str)
    parser.add_argument("-anno", help="WT amplicon sequence annotation (default ./anno/amplicon_kayser.tsv)", default="./anno/amplicon_kayser.tsv", type=str)
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
    parser.add_argument('-min_read_length', help="minimum read length in bp required after N-trimming (default 100)", default="100", type=int)
    parser.add_argument('-filter_reads', help="minimum number of copies of each read required for processing (1 to turn filter off, 2 (default) to discard unique reads)", default="2", type=int)
    parser.add_argument('-filter_ins_unique_reads', help="minimum number of unique reads required to support an insertion for it to be considered 'high confidence' (default 2)", default="2", type=int)
    parser.add_argument('-filter_ins_total_reads', help="minimum number of total reads required to support an insertion for it to be considered 'high confidence' (default 1)", default="1", type=int)
    parser.add_argument('-filter_ins_vaf', help="minimum variant allele frequency (VAF) required for an insertion to be considered 'high confidence' (default 0.001)", default="0.001", type=float)
    cmd_args = parser.parse_args()

    R1 = cmd_args.fastq1
    R2 = cmd_args.fastq2
    I1 = cmd_args.index1
    I2 = cmd_args.index2
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
    STATS_FILE = os.path.join(OUT_DIR, "stats.txt")
    CONFIG_FILE = os.path.join(OUT_DIR, "config.txt")

    COST_MATCH = cmd_args.match
    COST_MISMATCH = cmd_args.mismatch
    COST_GAPOPEN = cmd_args.gap_open
    COST_GAPEXTEND = cmd_args.gap_extend
    MIN_SCORE_INSERTS = cmd_args.minscore_inserts
    MIN_SCORE_ALIGNMENTS = cmd_args.minscore_alignments

    MIN_READ_LENGTH = cmd_args.min_read_length
    MIN_READ_COPIES = cmd_args.filter_reads
    MIN_TOTAL_READS = cmd_args.filter_ins_total_reads
    MIN_UNIQUE_READS = cmd_args.filter_ins_unique_reads
    MIN_VAF = cmd_args.filter_ins_vaf


    ANNO = read_annotation(ANNO_FILE)
    DOMAINS = get_domains(ANNO)
    REF = read_reference(REF_FILE).upper()

    ## CREATE OUTPUT FOLDER
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

    save_config(cmd_args, CONFIG_FILE)

    # stats are appended, remove previous output prior to new analysis
    try:
        os.remove(STATS_FILE)
    except OSError:
        pass
    save_stats("==== PROCESSING SAMPLE {} ====".format(SAMPLE), STATS_FILE)

    save_stats("-- Reading FASTQ files --", STATS_FILE)
    start_time = timeit.default_timer()
    reads = read_fastq(R1)

    # IF IT EXISTS:
    # --> reverse-complement R2 reads so that all reads can be aligned to the same reference
    if R2:
        reads_rev = read_fastq(R2)
        reads_rev_rev = parallelize(Read.reverse_complement,reads_rev,NKERN)
        reads = reads + reads_rev_rev
    print("Reading FASTQ files took {} s".format(timeit.default_timer() - start_time))
    save_stats("Number of total reads: {}".format(len(reads)), STATS_FILE)
    TOTAL_READS = len(reads)

    ## IF GIVEN, GET AND FILTER ON INDEX BQS
    start_time = timeit.default_timer()
    if I1 or I2:
        if I1:
            indices1_bqs = read_index_bqs(I1)
            indices_bqs = indices1_bqs
        if I2:
            indices2_bqs = read_index_bqs(I2)
            indices_bqs = indices2_bqs
        if I1 and I2:
            merged_indices_bqs = [index1_bqs + index2_bqs for index1_bqs,index2_bqs in zip(indices1_bqs, indices2_bqs)]
            indices_bqs = merged_indices_bqs
        reads = [read for read,index_bqs in zip(reads, indices_bqs) if average_bqs(index_bqs) >= MIN_BQS]
        print("Reading and filtering index BQS took {} s".format(timeit.default_timer() - start_time))
        save_stats("Number of total reads with index BQS >= {}: {} ({} %)".format(MIN_BQS, len(reads), len(reads) * 100 / TOTAL_READS), STATS_FILE)

    ## TRIM trailing AMBIGUOUS 'N's
    reads = [x for x in parallelize(Read.trim_n, reads, NKERN) if x is not None]
    save_stats("Number of total reads remainging after N-trimming: {} ({} %)".format(len(reads), len(reads) * 100 / TOTAL_READS), STATS_FILE)
    save_stats("Mean read length after N-trimming: {}".format(np.mean([read.length for read in reads])), STATS_FILE)

    ## FILTER ON BQS
    if MIN_BQS > 0:
        reads = [x for x in parallelize(Read.filter_bqs, reads, NKERN) if x is not None]
    save_stats("Number of total reads with mean BQS >= {}: {} ({} %)".format(MIN_BQS, len(reads), len(reads) * 100 / TOTAL_READS), STATS_FILE)

    # get unique reads and counts thereof
    start_time = timeit.default_timer()
    reads = get_unique_reads(reads)
    print("Getting unique reads took {} s\n".format(timeit.default_timer() - start_time))
    save_stats("Number of unique reads with mean BQS >= {}: {}".format(MIN_BQS,len(reads)), STATS_FILE)

    # FILTER UNIQUE READS
    # --> keep only those that exist at least twice
    # --> assumption: if it's not reproducible, it's not (true and clinically relevant)
    if MIN_READ_COPIES == 1:
        save_stats("Turned OFF unique reads filter!", STATS_FILE)
    else:
        reads = [read for read in reads if read.counts >= MIN_READ_COPIES ]
        save_stats("Number of unique reads with at least {} copies: {}".format(MIN_READ_COPIES,len(reads)), STATS_FILE)
    save_stats("Total reads remaining for analysis: {} ({} %)".format(sum([read.counts for read in reads]), sum([read.counts for read in reads]) * 100 / TOTAL_READS), STATS_FILE)

    ## ALIGN TO REF
    save_stats("\n-- Aligning to Reference --", STATS_FILE)
    start_time = timeit.default_timer()
    reads = parallelize(Read.align, reads, NKERN)
    print("Alignment took {} s".format(timeit.default_timer() - start_time))

    # FILTER BASED ON ALIGNMENT SCORE (INCL FAILED ALIGNMENTS WITH read.al_score is None!
    reads = filter_alignment_score(reads)
    reads = parallelize(Read.get_reference_range_covered, reads, NKERN)

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
        save_stats("Filtering {} / {} alignments with more than 3 unaligned primer bases".format( len(reads) - len(primers_filtered), len(reads)), STATS_FILE)
        reads = primers_filtered

    # FINAL STATS
    save_stats("Total reads remaining for analysis: {} ({} %)".format(sum([read.counts for read in reads]), sum([read.counts for read in reads]) * 100 / TOTAL_READS), STATS_FILE)

    # PRINT PASSING ALIGNMENTS
    # create output file directory for alignments print-outs
    needle_dir = os.path.join(OUT_DIR,'out_needle')
    if not os.path.exists(needle_dir):
        os.makedirs(needle_dir)

    for i,read in enumerate(reads):
        reads[i].al_file = 'needle_{}.txt'.format(i)
        print_alignment(reads[i], needle_dir)

    if not reads:
        save_stats("\nNO READS TO PROCESS!", STATS_FILE)
        quit()


    #######################################
    # CALCULATE COVERAGE
    start_time = timeit.default_timer()
    ref_coverage = []
    for coord, bp  in enumerate(REF):
        spanning_reads = [read for read in reads if coord >= read.ref_span[0] and coord <= read.ref_span[1]]
        spanning_reads_index = flatten_list([read.index for read in spanning_reads])
        ref_coverage.append(len(set(spanning_reads_index)))
    print("Calculating coverage took {} s".format(timeit.default_timer() - start_time))


    #######################################
    # COLLECT INSERTS
    save_stats("\n-- Looking for insertions & ITDs --", STATS_FILE)

    inserts = []
    start_time = timeit.default_timer()
    for read in reads:
        readn = np.array(list(read.al_seq))
        refn = np.array(list(read.al_ref))
        assert(len(readn) == len(refn))
        
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
    save_stats("{} insertions were found".format(len(inserts)), STATS_FILE)


    start_time = timeit.default_timer()
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
        alignments = [al for al in alignments if integral_insert_realignment(al[0],insert.length)]
        if not alignments:
            alignment_score = -1
        else:
            # how often is there more than one alignment?
            # --> (more than 1 => alignment is ambiguous)
            # --> Can I choose a smart one somehow? Otherwise return only one in the first place...
            alignment = alignments[-1]
            alignment_score, alignment_start, alignment_end = alignment[2:5]
            #print(bio.format_alignment(*alignment))
        if alignment_score >= min_score:
            #offset = abs(alignment_start - insert.start)
            tandem2_start = [i for i,bp in enumerate(alignment[0]) if bp != '-'][0]
            offset = abs(tandem2_start - insert.start)
            # offset = 1 for adjacent insert-tandem2
            # offset = insert.length-1 for adjacent tandem2-insert
            # --> (for tandem2-insert: offset = abs((insert_start - insert.length +1) - insert_start))
            if (offset == 1 or offset == insert.length - 1) or (insert.trailing and (insert.trailing_end == 3 and alignment_start < insert.start) or (insert.trailing_end == 5 and alignment_start > insert.start)):
                # do not allow gaps in tandems of trailing ITDs
                # -> also filters tandems not covered by read at all 
                #    such as small trailing inserts by chance also found in some other part of the the reference
                #    (can never be the case for non-trailing ITDs anyway)
                # --> careful: alignment_end is exclusive coord, i.e. the index of the first bp after the alignment!
                #if alignment_start >= insert.reads[0].ref_span[0] or alignment_end-1 <= insert.reads[1].ref_span[1]:
                if alignment_start >= insert.reads[0].ref_span[0] and alignment_end-1 <= insert.reads[0].ref_span[1]:
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
                        tandem2_start=tandem2_start,
                        #tandem2_start=alignment_start,
                        external_bp=abs(tandem2_start - alignment_start)
                        #external_bp=insert.length - (alignment_end - alignment_start)
                        )
                    itds.append(itd)
                else:
                    print("ITD's tandem not covered by read")
                    insert.print()
                    insert.reads[0].print()
                    print(bio.format_alignment(*alignment))
                    print(alignment_start)
                    print(alignment_end)

    # in case any out-of-frame insert was untrailed: remove it from list of inserts
    inserts[:] = [insert for insert in inserts if insert.trailing or insert.length % 3 == 0]
    inserts = sorted(inserts, key=Insert.get_seq)
    itds = sorted(itds, key=Insert.get_seq)
    print("Collecting ITDs took {} s".format(timeit.default_timer() - start_time))
    save_stats("{} ITDs were found".format(len(itds)), STATS_FILE)


    ########################################
    # MERGE INSERTS
    save_stats("\n-- Merging results --", STATS_FILE)

    merge_dic = {"insertions": inserts, "itds": itds}
    all_merged = {}
    for inserts_type,inserts_ in merge_dic.items():
        all_merged[inserts_type] = []
        suffix = ""
        # turn Insert objects into InsertCollection to keep merging methods simple and not have to distinguish between the two
        to_merge = [InsertCollection(insert) for insert in inserts_]
        for condition,abrev in [
                ("is-same",""),
                ("is-similar","similar"),
                ("is-close","close"),
                ("is-same_trailing","trailing")]:
            to_merge = merge(to_merge, condition)
            all_merged[inserts_type].append(to_merge)
            save_stats("{} {} remain after merging".format(len(to_merge), inserts_type), STATS_FILE)
            suffix = suffix + condition
            save_to_file([insert.rep for insert in to_merge], inserts_type + "_collapsed-" + suffix + ".tsv")
            suffix = suffix + "_"


    # save as list of Insert(s) to process further as before (vs continuing with list of InsertCollection(s)
    final_merged = {}
    for inserts_type in merge_dic:
        final_merged[inserts_type] = [insert.rep for insert in all_merged[inserts_type][-1]]
    print("Merging took {} s".format(timeit.default_timer() - start_time))


    ########################################
    # FILTER INSERTS
    save_stats("\n-- Filtering --", STATS_FILE)

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
            save_stats("Filtered {} / {} {} based on the {}".format(
                len(passed) - sum(passed), len(passed), inserts_type, filter_type), STATS_FILE)
        save_stats("{} {} remain after filtering!".format(len(filtered), inserts_type), STATS_FILE)
        save_to_file(filtered, inserts_type + "_collapsed-" + suffix + "hc.tsv")
        final_filtered[inserts_type] = filtered
    print("Filtering took {} s".format(timeit.default_timer() - start_time))


