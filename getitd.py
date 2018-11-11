import Bio.pairwise2 as bio
import timeit
import collections
import itertools
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


def save_config(config, filename):
    """
    Write timestamp and commandline arguments to file.

    Args:
        config (dict): Config parameters and values to write.
        filename (str): Name of the file to write to.
    """
    with open(filename, "w") as f:
        f.write("Commandline_argument\tValue\n")
        f.write("Time\t{}\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%d")))
        for param in sorted(config.keys()):
            if param not in ["ANNO", "DOMAINS"]:
                f.write("{}\t{}\n".format(param, config[param]))

def load_config(filename):
    """
    Load config parameters from file.

    Args:
        filename (str): Name of the file to read config from.

    Returns:
        Dictionary with config parameter - value pairs.
    """
    config = {}
    with open(filename, "r") as f:
        for line in f:
            key, val = line.strip("\n").split("\t")
            try:
                config[key] = float(val)
            except:
                config[key] = val
    return config



# child processes spawned on Windows by multiprocessing do not
# receive variables set in __main__ of parent process
# -->  they cannot access config {} values set in __main__
# --> to circumvent this, __main__ saves config and children
#     spawned by multiprocessing load it from file
if __name__ == '__mp_main__':
    try:
        current_dir = os.getcwd()
        config = load_config(os.path.join(current_dir, "config.txt"))
    except OSError:
        print("NO CONFIG FOUND")
else:
    # mimic Windows style process spawning on Linux:
    # multiprocessing.set_start_method("spawn")
    config = {}




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

    def trim_n(self, config=config):
        """
        Trim trailing N's at Read's ends.

        Returns:
            Trimmed Read.
        """
        self.seq = self.seq.rstrip('N')
        if self.bqs:
            self.bqs = self.bqs[:len(self.seq)]
        self.seq = self.seq.lstrip('N')
        if self.bqs:
            self.bqs = self.bqs[::-1][:len(self.seq)][::-1]
        assert len(self.seq) == len(self.bqs)
        self.length = len(self.seq)
        if self.length < config["MIN_READ_LENGTH"]:
            return None
        return self


    def filter_bqs(self, config=config):
        """
        Filter Read based on its average BQS.

        Returns:
            Read when it passes the average BQS filter,
            None otherwise.

            When no BQS is set, also return the Read.
        """
        if not self.bqs or average_bqs(self.bqs) >= config["MIN_BQS"]:
            if not self.index_bqs or average_bqs(self.index_bqs) >= config["MIN_BQS"]:
                return self
        return None

    def align(self, config=config):
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
        alignment = bio.align.globalcs(self.seq, config["REF"], get_alignment_score,
            config["COST_GAPOPEN"], config["COST_GAPEXTEND"], penalize_end_gaps=False)
        if alignment:
            self.al_seq, self.al_ref, self.al_score = alignment[-1][0:3]

        if config["TECH"] == '454':
            rev = self.reverse_complement()
            rev_alignment = bio.align.globalcs(rev.seq, config["REF"], get_alignment_score,
                config["COST_GAPOPEN"], config["COST_GAPEXTEND"], penalize_end_gaps=False)
            if (rev_alignment and
                    (self.al_score is None or rev_alignment[-1][2] > self.al_score)):
                rev.al_seq, rev.al_ref, rev.al_score = rev_alignment[-1][0:3]
                return rev
        return self

    def get_ref_span(self):
        """
        Get most 3' and 5' reference base covered by Read.
        Coordinates are 0-based and relative to the WT reference.

        Returns:
            Updated read.
        """
        refn = np.array(list(self.al_ref))
        readn = np.array(list(self.al_seq))
        ref_covered_bp = np.where(readn[refn != '-'] != '-')
        self.ref_span = [np.min(ref_covered_bp), np.max(ref_covered_bp)]

        # extend ref_span for trailing inserts
        if refn[-1] == '-':
            assert readn[-1] != '-'
            self.ref_span[-1] = len(refn) -1
        if refn[0] == '-':
            assert readn[0] != '-'
            self.ref_span[0] = 0
        return self

    def reorder_trailing_inserts(self):
        """
        Change the alignment of

        Ref     ----xxxxxxXXXXXXXX
        Sample  xxxx------XXXXXXXX

        to

        Ref     xxxx------XXXXXXXX
        Sample  ----xxxxxxXXXXXXXX

        and analoguously treat 3' trailing inserts to guarantee
        that reference seq is more terminal than sample insert.
        This is required to guarantee correct ref_span calculation:
        ref_span (i.e. reference coordinates spanned by sample read)
        corresponds to most terminal reference coordinates that
        reside within sample read alignment. In the above example,
        the 5' ref_span is fixed from 0 to 4 to reflect that
        the first four reference bases are not really part of the
        sample read.
        """
        if self.al_ref[0] == '-':
            insert_idxs = get_gaps(self.al_ref)
            if self.al_seq[insert_idxs[0][-1] + 1] == '-':
                del_idxs = get_gaps(self.al_seq)
                self.print()
                print(insert_idxs)
                print(del_idxs)
                self.al_seq = self.al_seq[del_idxs[0][0] : del_idxs[0][-1] + 1] \
                        + self.al_seq[insert_idxs[0][0] : insert_idxs[0][-1] + 1] \
                        + self.al_seq[del_idxs[0][-1] + 1 : ]
                self.al_ref = self.al_ref[del_idxs[0][0] : del_idxs[0][-1] + 1] \
                        + self.al_ref[insert_idxs[0][0] : insert_idxs[0][-1] + 1] \
                        + self.al_ref[del_idxs[0][-1] + 1 : ]
                self.print()
        return self


    def get_inserts(self):
        """
        Collect all inserts contained within a given Read.

        Returns:
            List of collected Insert objects.
        """
        inserts = []

        # if read contains insert
        if '-' in self.al_ref:
            readn = np.array(list(self.al_seq))
            refn = np.array(list(self.al_ref))
            assert(len(readn) == len(refn))

            # collect all inserts' alignment coords
            insert_idxs_list = get_gaps(self.al_ref)

            # check & process each insert found
            for insert_idxs in insert_idxs_list:
                if len(insert_idxs) >= 6 and "N" not in readn[insert_idxs]:
                    insert_start = insert_idxs[0]
                    insert_end = insert_idxs[-1]
                    insert = Insert(
                        seq=self.al_seq[insert_start:insert_end+1],
                        start=insert_start,
                        end=insert_end,
                        reads=[self],
                        counts=self.counts)
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
                    insert.trailing = (self.sense == 1 and insert.trailing_end == 3) or (self.sense == -1 and insert.trailing_end == 5)

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
                        inserts.append(insert)
        return inserts


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
        self.sense = set([read.sense for read in reads])

    def get_seq(self):
        """"
        Return Insert's sequence.
        Implemented to sort list of Inserts by their sequence.
        """
        return self.seq

    def set_insertion_site(self):
        """
        Get and set Insertion's insertion_site, i.e. the WT bp
        following the insert.
        """
        self.insertion_site = self.end + 1
        return self

    def set_sense(self):
        """
        Get and set Insert's supporting reads' sense.
        """
        self.sense = set([read.sense for read in self.reads])
        return self

    def set_specific_sense(self, sense):
        """
        Set Insert's sense manually.

        Args:
            sense (set): sense to manually set to.
        """
        self.sense = sense
        return self

    def set_coverage(self, ref_coverage, ref_coverage_frwd, ref_coverage_rev):
        """
        Set Insert's coverage based on supporting reads' sense.
        """
        if self.sense == {1}:
            self.coverage = ref_coverage_frwd[copy.deepcopy(self).norm_start().start]
        elif self.sense == {-1}:
            self.coverage = ref_coverage_rev[copy.deepcopy(self).norm_start().start]
        elif self.sense == {1,-1}:
            self.coverage = ref_coverage[copy.deepcopy(self).norm_start().start]
        return self

    def calc_vaf(self):
        """
        Calculate and set Insert's variant allele frequency (VAF).

        Returns:
            Insert.
        """
        self.vaf = dc.Decimal(self.counts) / self.coverage * 100
        assert self.vaf >= 0 and self.vaf <= 100
        return self

    def get_itd(self, config=config):
        """
        Check whether Insert qualifies as an ITD.

        Args:
            config (dict): Dict containing analysis parameters.

        Returns:
            ITD if the Insert qualifies as one.
            None otherwise.
        """
        alignments = bio.align.localcs(self.seq, config["REF"], get_alignment_score, config["COST_GAPOPEN"], config["COST_GAPEXTEND"])
        alignments = [al for al in alignments if integral_insert_realignment(al[0],self.length)]
        if not alignments:
            return None
        # bio.align arguments: seq1, seq2, match-score, mismatch-score, gapopen-score, gapextend-score
        # output: list of optimal alignments, each a list of seq1, seq2, score, start-idx, end-idx
        # to print: print(bio.format_alignment(*alignment))
        alignment = alignments[-1]
        alignment_score, alignment_start, alignment_end = alignment[2:5]
        if alignment_score >= get_min_score(self.seq, config["REF"], config["MIN_SCORE_ALIGNMENTS"]):
            tandem2_start = [i for i,bp in enumerate(alignment[0]) if bp != '-'][0]
            offset = abs(tandem2_start - self.start)
            # offset = 1 for adjacent insert-tandem2
            # offset = insert.length-1 for adjacent tandem2-insert
            # --> (for tandem2-insert: offset = abs((insert_start - insert.length +1) - insert_start))
            if (offset == 1 or offset == self.length - 1) or (self.trailing and (self.trailing_end == 3 and alignment_start < self.start) or (self.trailing_end == 5 and alignment_start > self.start)):
                # do not allow gaps in tandems of trailing ITDs
                # -> also filters tandems not covered by read at all
                #    such as small trailing inserts by chance also found in some other part of the the reference
                #    (can never be the case for non-trailing ITDs anyway)
                # --> careful: alignment_end is exclusive coord, i.e. the index of the first bp after the alignment!
                if alignment_start >= self.reads[0].ref_span[0] and alignment_end-1 <= self.reads[0].ref_span[1] and tandem2_start + offset <= self.reads[0].ref_span[1]:
                    # if by chance insert is completely contained within read in spite of it being trailing
                    # (i.e. insert and tandem have the same length & are adjacent)
                    # --> revert trailing to be able to apply more stringent filters of non-trailing inserts
                    if self.trailing and (offset == 1 or offset == self.length - 1):
                        self.trailing = False
                        print("UNTRAILED: {}".format(vars(read)))
                        if self.length % 3 != 0:
                            print("BUT NOT IN FRAME!!!")
                            return None
                    return ITD(
                        self,
                        offset=offset,
                        tandem2_start=tandem2_start,
                        external_bp=abs(tandem2_start - alignment_start)
                        )
                else:
                    if tandem2_start + offset > self.reads[0].ref_span[1]:
                        print("Too long!")
                    print("ITD's tandem not covered by read and therefore discarded")
                    self.print()
                    self.reads[0].print()
                    print(bio.format_alignment(*alignment))
                    print("alignment_start: {}".format(alignment_start))
                    print("alignment_end: {}".format(alignment_end))
        return None


    def norm_start(self, config=config):
        """
        Normalize Insert's start coordinate to [0, len(REF)[.

        Returns:
            Insert.
        """
        self.start = min(max(0,self.start), len(config["REF"])-1)
        return self

    def prep_for_save(self):
        """
        Prepare Insert for saving to file by normalizing start
        coordinates to [0, len(REF)[.
        Convert coordinates from 0- to 1-based to match
        needle output files and annotation table.
        Set end coordinate equal to start -> insertion_site
        will be the WT base after the insertion, i.e. start +1
        (start + length has no meaning for non-ITD insertions)

        Returns:
            Prepared Insert.
        """
        to_save = copy.deepcopy(self)
        to_save = to_save.norm_start()
        to_save.start += 1
        to_save.end = to_save.start
        return to_save

    def annotate(self, coord, to_annotate, config=config):
        """
        Add specific annotation to Insertion.

        Args:
            coord (str): Name of the Insertion's attribute to be
                    annotated, i.e. self.coord is the coordinate
                    that will be annotated and one of {"start", "end",
                    "insertion_site"}
            to_annotate (str): Name of the annotation to add, must
                    be one of {"chr13_bp", "transcript_bp", "protein_as"}.
            config (dict): Dictionary where config["ANNO"] contains
                    pandas DataFrame with annotation to retrieve.
        Returns:
            Annotated Insertion.
        """
        anno = config["ANNO"]
        annotated = anno.loc[anno["amplicon_bp"] == getattr(self, coord), to_annotate].values[0]
        setattr(self, coord + "_" + to_annotate, annotated)
        return self


    def annotate_domains(self, domains):
        """
        Get annotated domains of Insert's flanking bp.
        --> self.start is WT base preceding insert
        --> self.start +1 is WT base after insert

        Args:
            domains: List of tuples (domain_name, start_coord, end_coord)
            where coordinates refer to amplicon bp and the list contains
            all of its annotated domains.

        Returns:
            Subset of domains, containing only those of the Insert's
            flanking base pairs.
        """
        annotated = []
        for domain,start,end in domains:
            if self.start <= end and self.start + 1 >= start:
                annotated.append(domain)
        self.domains = annotated
        return self


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
            return abs(self.tandem2_start - that.tandem2_start) <= max(self.length, that.length)
        return abs(self.start - that.start) <= self.length + that.length ## <---- what's the eqivalent for ITDs?? Does this even apply for ITDs? Can trailing ITDs have different length and still describe the same mutation?

    def is_similar_to(self, that, config=config):
        """
        Test whether two Inserts' sequences are similar.

        Args:
            that: Insert to compare to.

        Returns:
            True when they are similar, False otherwise.
        """
        min_score = get_min_score(self.seq, that.seq, config["MIN_SCORE_INSERTS"])
        al_score = bio.align.globalcs(self.seq, that.seq, get_alignment_score, config["COST_GAPOPEN"], config["COST_GAPEXTEND"], one_alignment_only=True, score_only=True, penalize_end_gaps=False)
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
            return self.trailing and that.trailing and self.trailing_end == that.trailing_end and self.sense.intersection(that.sense) and self.is_similar_to(that) and self.is_close_to(that)
        assert False


    def filter_unique_supp_reads(self, config=config):
        """
        Test whether Insert is supported by a given number of
        distinct supporting reads.

        Returns:
            True when that is the case, False otherwise.
        """
        return len(self.reads) >= config["MIN_UNIQUE_READS"]

    def filter_total_supp_reads(self, config=config):
        """
        Test whether Insert is supported by a given minimum number
        of total supporting reads.

        Returns:
            True when that is the case, False otherwise.
        """
        return self.counts >= config["MIN_TOTAL_READS"]

    def filter_vaf(self, config=config):
        """
        Test whether Insert has at least a given minimum VAF.

        Returns:
            True when that is the case, False otherwise.
        """
        return self.vaf >= config["MIN_VAF"]


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
        self.sense = insert.sense
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

    def annotate_domains(self, domains):
        """
        Get annotated domains of the ITD's WT tandem.

        Args:
            domains: List of tuples (domain_name, start_coord, end_coord)
            where coordinates refer to amplicon bp and the list contains
            all of its annotated domains.

        Returns:
            Subset of domains, containing only those duplicated by the ITD.
        """
        annotated = []
        for domain,start,end in domains:
            if self.start <= end and self.end >= start:
                annotated.append(domain)
        self.domains = annotated
        return self

    def prep_for_save(self, config=config):
        """
        Prepare ITD for saving to file.

        Set the insert start coordinate to that of the second tandem.
        Normalize start coordinate to [0, len(REF)[ and
        adjust end coordinate accordingly.
        Convert coordinates from 0- to 1-based to match
        needle output files and annotation table.
        For trailing ITDs, change the length to the offset, i.e.
        distance between insert and second tandem.

        Return:
            Prepared ITD.
        """
        to_save = copy.deepcopy(self)
        to_save = to_save.fix_trailing_length()
        to_save.start = to_save.tandem2_start
        to_save.end = to_save.start + to_save.length - 1
        if to_save.end > len(config["REF"]):
            self.print()
            self.reads[0].print()
            to_save.print()
        assert to_save.end <= len(config["REF"])
        to_save = to_save.norm_start()
        to_save.start += 1
        to_save.end += 1
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
        self.rep = self.rep.set_sense()
        self.rep.coverage = max([insert.coverage for insert in self.inserts])
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
        if self.rep.sense != insert.rep.sense:
            new_sense = self.rep.sense.union(insert.rep.sense)
            self.inserts = [insert.set_specific_sense(new_sense).set_coverage() for insert in self.inserts]
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

def get_gaps(seq):
    """
    Extract gap indices from alignment string.

    Args:
        seq (str): Alignment string of one of the aligned
                sequences.

    Returns:
        List of lists, each containing indices of consecutive gaps.

    """
    gap_idxs_sep = []

    if '-' in seq:
        seqn = np.array(list(seq))
        gap_idxs_all = np.where(seqn == '-')[0]
        # x = i,e from enumerate(), lambda return diff between e and i, groupby breaks up
        #   list whenever this difference changes (changes between non-consecutive indices)
        for key, group in itertools.groupby(enumerate(gap_idxs_all), lambda x: x[1]-x[0]):
            gap_idxs_sep.append([e for i,e in group])
        assert np.all(np.concatenate(gap_idxs_sep) == gap_idxs_all)
    return gap_idxs_sep


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
    seq_coord += int(len(seq) > seq.count('-'))
    f.write(' ' * (pre_width - get_number_of_digits(seq_coord) +1))
    f.write(str(seq_coord) + ' ')
    f.write(seq)
    seq_coord = seq_coord + len(seq) - seq.count('-') -1 +int(len(seq) == seq.count('-'))
    f.write(' ' * (post_width - get_number_of_digits(seq_coord)))
    f.write(str(seq_coord) + '\n')
    return seq_coord

def print_alignment(read, out_dir, config=config):
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
        f.write('#    {},\n'.format(config["COST_GAPOPEN"]))
        f.write('#    {})\n'.format(config["COST_GAPEXTEND"]))
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
        f.write('# Gap_penalty: {}\n'.format(config["COST_GAPOPEN"]))
        f.write('# Extend_penalty: {}\n'.format(config["COST_GAPEXTEND"]))
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
        seq_coord = 0
        ref_coord = 0
        for s, a, r in alignment_chunks:
            seq_coord = print_alignment_seq(s, seq_coord,pre_width,post_width,f)
            print_alignment_connection(a, pre_width, f)
            ref_coord = print_alignment_seq(r, ref_coord,pre_width,post_width,f)
            f.write('\n')

        f.write('\n')
        f.write('#---------------------------------------\n')
        f.write('#---------------------------------------\n')


def get_alignment_score(char1,char2, config=config):
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
        return config["COST_MATCH"]
    elif char1 == 'Z' or char2 == 'Z':
        return -np.inf
    else:
        return config["COST_MISMATCH"]

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
    return (min(len(seq1),len(seq2)) * config["COST_MATCH"]) * min_score
    # min score with penalize_trailing_gaps=True:
    #return (min(len(seq1),len(seq2)) * config["COST_MATCH"] + (abs(len(seq1) - len(seq2)) -1) * config["COST_GAPEXTEND"] + config["COST_GAPOPEN"]) * min_score


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
                read_seq = f.readline().rstrip(os.linesep)
                read_desc = f.readline()
                read_bqs = f.readline().rstrip(os.linesep)
                assert len(read_seq) == len(read_bqs)
                reads.append(Read(seq=read_seq, index=read_index, bqs=read_bqs, sense=1))
                line = f.readline()
                read_index += 1
    # catch missing file or permissions
    except IOError as e:
        print("---\nCould not read fastq file {}!\n---".format(fastq_file))
    return reads



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
        still_need_to_merge = False
        merged = []
        for insert_collection in inserts:
            was_merged = False
            for minsert_collection in merged[::-1]:
                if minsert_collection.should_merge(insert_collection, condition):
                    minsert_collection = minsert_collection.merge(insert_collection)
                    was_merged = True
                    still_need_to_merge = True
                    break
            if not was_merged:
                merged.append(insert_collection)
        inserts = merged
    return merged

def save_to_file(inserts, filename, config=config):
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
        to_save = [insert.prep_for_save() for insert in inserts]
        if config["ANNO"] is not None:
            for insert in to_save:
                insert = insert.set_insertion_site()
                insert = insert.annotate_domains(config["DOMAINS"])
                cols = ["domains"]
                for to_annotate in ["start", "end", "insertion_site"]:
                    for coord in ["chr13_bp", "transcript_bp", "protein_as"]:
                        insert = insert.annotate(to_annotate, coord)
                        cols.append(to_annotate + "_" + coord)
                insert = insert.annotate("insertion_site", "region")
                cols.append("insertion_site_domain") # rename in df...
                #cols = ["domains", "start_chr13_bp", "start_transcript_bp", "start_protein_as", "end_chr13_bp", "end_transcript_bp", "end_protein_as", "insertion_site_protein_as"]

        dict_ins = {}
        for key in vars(to_save[0]):
            dict_ins[key] = tuple(vars(insert)[key] for insert in to_save)

        df_ins =  pd.DataFrame(dict_ins)
        df_ins["sample"] = [config["SAMPLE"]] * len(to_save)
        df_ins["ar"] = [vaf_to_ar(insert.vaf) for insert in to_save]
        df_ins["counts_each"] = [[read.counts for read in insert.reads] for insert in to_save]
        df_ins["file"] = [[read.al_file for read in insert.reads] for insert in to_save]

        if 'external_bp' in df_ins:
            cols = ['external_bp'] + cols

        cols = ['sample','length', 'start', 'vaf', 'ar', 'coverage', 'counts', 'trailing', 'seq', 'sense'] + cols + ['file']
        df_ins = df_ins.rename(index=str, columns={"insertion_site_region": "insertion_site_domain"})
        df_ins[cols].sort_values(by=['length','start','vaf']).to_csv(os.path.join(config["OUT_DIR"],filename), index=False, float_format='%.2e', sep='\t')

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


def filter_alignment_score(reads, config=config):
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
                read.seq, config["REF"], config["MIN_SCORE_ALIGNMENTS"])]
    save_stats("Filtering {} / {} low quality alignments with a score < {} % of max".format(
            len(reads) - len(reads_filtered), len(reads), config["MIN_SCORE_ALIGNMENTS"] *100), config["STATS_FILE"])
    return reads_filtered

def filter_alignments_with_gaps_in_primers(reads):
    """
    Filter reads that are not aligned fully to at least one primer.

    Args:
        reads: List of Read objects to filter

    Returns:
        List of passing reads
    """
    rev_primer = 'GGTTGCCGTCAAAATGCTGAAAG'
    fwrd_primer = 'GCAATTTAGGTATGAAAGCCAGCTAC'
    filtered = [read for read in reads if (
            (read.sense == -1
            and rev_primer in read.al_ref
            and read.al_seq.count('-', read.al_ref.find(rev_primer), read.al_ref.find(rev_primer) + len(rev_primer)) <= 0
            ) or
            (read.sense == 1
            and fwrd_primer in read.al_ref
            and read.al_seq.count('-', read.al_ref.find(fwrd_primer), read.al_ref.find(fwrd_primer) + len(fwrd_primer)) <= 0))]
    save_stats("Filtering {} / {} alignments with indels in primer bases".format( len(reads) - len(filtered), len(reads)), config["STATS_FILE"])
    return filtered


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

def parse_config_from_cmdline(config=config):
    """
    Get analysis parameters from commandline.

    Args:
        config (dict): Dict to save parameters to

    Returns:
        Filled config dict
    """
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("fastq1", help="FASTQ file of forward reads (REQUIRED)")
    parser.add_argument("fastq2", help="FASTQ file of reverse reads (REQUIRED)")
    parser.add_argument("sampleID", help="sample ID used as output folder prefix (REQUIRED)")
    parser.add_argument("minBQS", help="minimum average base quality score (BQS) required by each read (default 30)", type=int, default=30, nargs='?')
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
    parser.add_argument('-min_read_length', help="minimum read length in bp required after N-trimming (default 100)", default="100", type=int)
    parser.add_argument('-filter_reads', help="minimum number of copies of each read required for processing (1 to turn filter off, 2 (default) to discard unique reads)", default="2", type=int)
    parser.add_argument('-filter_ins_unique_reads', help="minimum number of unique reads required to support an insertion for it to be considered 'high confidence' (default 2)", default="2", type=int)
    parser.add_argument('-filter_ins_total_reads', help="minimum number of total reads required to support an insertion for it to be considered 'high confidence' (default 1)", default="1", type=int)
    parser.add_argument('-filter_ins_vaf', help="minimum variant allele frequency (VAF) required for an insertion to be considered 'high confidence' (default 0.006)", default="0.006", type=float)
    cmd_args = parser.parse_args()

    config["R1"] = cmd_args.fastq1
    config["R2"] = cmd_args.fastq2
    config["SAMPLE"] = cmd_args.sampleID
    config["MIN_BQS"] = cmd_args.minBQS
    config["REF_FILE"] = cmd_args.reference
    config["ANNO_FILE"] = cmd_args.anno
    config["TECH"] = cmd_args.technology
    config["NKERN"] = cmd_args.nkern
    config["OUT_DIR"] = '_'.join([config["SAMPLE"],'minBQS', str(config["MIN_BQS"])])
    config["STATS_FILE"] = os.path.join(config["OUT_DIR"], "stats.txt")
    config["CONFIG_FILE"] = os.path.join(config["OUT_DIR"], "config.txt")

    config["COST_MATCH"] = cmd_args.match
    config["COST_MISMATCH"] = -abs(cmd_args.mismatch)
    config["COST_GAPOPEN"] = -abs(cmd_args.gap_open)
    config["COST_GAPEXTEND"] = -abs(cmd_args.gap_extend)
    config["MIN_SCORE_INSERTS"] = cmd_args.minscore_inserts
    config["MIN_SCORE_ALIGNMENTS"] = cmd_args.minscore_alignments

    config["MIN_READ_LENGTH"] = cmd_args.min_read_length
    config["MIN_READ_COPIES"] = cmd_args.filter_reads
    config["MIN_TOTAL_READS"] = cmd_args.filter_ins_total_reads
    config["MIN_UNIQUE_READS"] = cmd_args.filter_ins_unique_reads
    config["MIN_VAF"] = cmd_args.filter_ins_vaf

    # make all input & output file / folder names absolute paths
    for path in ["R1", "R2", "REF_FILE", "ANNO_FILE", "OUT_DIR", "STATS_FILE", "CONFIG_FILE"]:
        if config[path] and not os.path.isabs(config[path]):
            config[path] = os.path.join(os.getcwd(), config[path])
    return config


def get_reads(config=config):
    """
    Read in FASTQ files.

    Args:
        config (dict): Dict containing analysis parameters.

    Returns:
        List of Read objects, one for each read from input FASTQ files.
    """
    save_stats("-- Reading FASTQ files --", config["STATS_FILE"])
    start_time = timeit.default_timer()
    reads = read_fastq(config["R1"])

    # IF IT EXISTS:
    # --> reverse-complement R2 reads so that all reads can be aligned to the same reference
    if config["R2"]:
        reads_rev = read_fastq(config["R2"])
        reads_rev_rev = parallelize(Read.reverse_complement, reads_rev, config["NKERN"])
        reads = reads + reads_rev_rev
    print("Reading FASTQ files took {} s".format(timeit.default_timer() - start_time))
    save_stats("Number of total reads: {}".format(len(reads)), config["STATS_FILE"])
    return reads


def filter_sequencing_adapter_artefacts(inserts):
    """
    Filter Inserts to remove false positives originating from
    sequencing adapters partially contained within the read.
    (Alternatively, one could trim reads in advance but using
    the generated alignment instead should ensure that really
    no unnecessary sequence is removed.)

    Args:
        inserts: List of Inserts to filter.

    Returns:
        List of passing Inserts.
    """
    fwrd_adapter = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGA"
    rev_adapter = "AGACAGAGAATATGTGTAGAGGCTCGGGTGCTCTG".translate(str.maketrans('ATCGatcg','TAGCtagc'))[::-1]
    non_adapter_inserts = [insert for insert in inserts if
            not insert.trailing
            or (insert.trailing_end == 5 and not insert.is_similar_to(Insert(seq=fwrd_adapter[-insert.length:], start=0, end=0, counts=0)))
            or (insert.trailing_end == 3 and not insert.is_similar_to(Insert(seq=rev_adapter[-insert.length:], start=0, end=0, counts=0)))
            ]

    save_stats("{}/{} insertions were part of adapters and filtered".format(len(inserts) - len(non_adapter_inserts), len(inserts)), config["STATS_FILE"])
    return non_adapter_inserts


def get_merged_inserts(inserts, type_, config):
    """
    Merge Inserts based on different conditions.

    Args:
        inserts: List of Inserts to merge.
        type_ (str): "insertions" or "itds",
                     description for output files.
        config (dict): Dict containing analysis parameters.

    Returns:
        List of merged Inserts.
    """
    save_stats("\n-- Merging {} --".format(type_), config["STATS_FILE"])

    merged = []
    # turn Insert objects into InsertCollection to keep merging methods simple and not have to distinguish between the two
    to_merge = [InsertCollection(insert) for insert in inserts]
    suffix = ""
    for condition in [
            "is-same",
            "is-similar",
            "is-close",
            "is-same_trailing"]:
        to_merge = merge(to_merge, condition)
        merged.append(to_merge)
        save_stats("{} {} remain after merging".format(len(to_merge), type_), config["STATS_FILE"])
        suffix = suffix + condition
        save_to_file([insert.rep for insert in to_merge], type_ + "_collapsed-" + suffix + ".tsv")
        suffix = suffix + "_"

    # convert InsertCollection back to list of (representative) Inserts
    return [insert.rep for insert in merged[-1]]


def get_hc_inserts(inserts, type_, suffix="", config=config):
    """
    Filter for high-confidence (hc) inserts only.

    Args:
        inserts: List of Insertion() objects to filter

        type_ (str): Description of type of inserts,
                     "insertions" or "itds",
                     for output files.

        config (dict): Dict containing analysis parameters.

    Returns:
        Dict with hc inserts, same format as input dict
    """
    save_stats("\n-- Filtering {} --".format(type_), config["STATS_FILE"])

    filter_dic = {
        "number of unique supporting reads": Insert.filter_unique_supp_reads,
        "number of total supporting reads": Insert.filter_total_supp_reads,
        "vaf": Insert.filter_vaf}

    filtered = copy.deepcopy(inserts)
    for filter_type, filter_ in filter_dic.items():
        passed = [filter_(insert) for insert in filtered]
        filtered = [insert for (insert, pass_) in zip(filtered, passed) if pass_]

        save_stats("Filtered {} / {} {} based on the {}".format(
                len(passed) - sum(passed), len(passed), type_, filter_type), config["STATS_FILE"])
    save_stats("{} {} remain after filtering!".format(len(filtered), type_), config["STATS_FILE"])
    save_to_file(filtered, type_ + suffix + ".tsv")

    return filtered

def main(config=config):
    config["ANNO"] = read_annotation(config["ANNO_FILE"])
    config["DOMAINS"] = get_domains(config["ANNO"])
    config["REF"] = read_reference(config["REF_FILE"]).upper()


    ## CREATE OUTPUT FOLDER
    if not os.path.exists(config["OUT_DIR"]):
        os.makedirs(config["OUT_DIR"])

    ## CHANGE TO OUTPUT FOLDER
    #  this is required for parallel child processes to retrieve
    #  the correct config.txt file later on despite static / constant filename
    os.chdir(config["OUT_DIR"])
    save_config(config, config["CONFIG_FILE"])

    ## REMOVE OLD STATS FILE & START CREATING A NEW ONE
    try:
        os.remove(config["STATS_FILE"])
    except OSError:
        pass
    save_stats("==== PROCESSING SAMPLE {} ====".format(config["SAMPLE"]), config["STATS_FILE"])

    ## READ FASTQ FILES
    reads = get_reads(config)
    TOTAL_READS = len(reads)

    ## TRIM trailing AMBIGUOUS 'N's
    reads = [x for x in parallelize(Read.trim_n, reads, config["NKERN"]) if x is not None]
    save_stats("Number of total reads remainging after N-trimming: {} ({} %)".format(len(reads), len(reads) * 100 / TOTAL_READS), config["STATS_FILE"])
    save_stats("Mean read length after N-trimming: {}".format(np.mean([read.length for read in reads])), config["STATS_FILE"])

    ## FILTER ON BQS
    if config["MIN_BQS"] > 0:
        reads = [x for x in parallelize(Read.filter_bqs, reads, config["NKERN"]) if x is not None]
    save_stats("Number of total reads with mean BQS >= {}: {} ({} %)".format(config["MIN_BQS"], len(reads), len(reads) * 100 / TOTAL_READS), config["STATS_FILE"])

    ## GET UNIQUE READS AND COUNTS THEREOF
    start_time = timeit.default_timer()
    reads = get_unique_reads(reads)
    print("Getting unique reads took {} s\n".format(timeit.default_timer() - start_time))
    save_stats("Number of unique reads with mean BQS >= {}: {}".format(config["MIN_BQS"],len(reads)), config["STATS_FILE"])

    # FILTER UNIQUE READS
    # --> keep only those that exist at least twice
    # --> assumption: if it's not reproducible, it's not (true and clinically relevant)
    if config["MIN_READ_COPIES"] == 1:
        save_stats("Turned OFF unique reads filter!", config["STATS_FILE"])
    else:
        reads = [read for read in reads if read.counts >= config["MIN_READ_COPIES"] ]
        save_stats("Number of unique reads with at least {} copies: {}".format(config["MIN_READ_COPIES"],len(reads)), config["STATS_FILE"])
    save_stats("Total reads remaining for analysis: {} ({} %)".format(sum([read.counts for read in reads]), sum([read.counts for read in reads]) * 100 / TOTAL_READS), config["STATS_FILE"])

    ## ALIGN TO REF
    save_stats("\n-- Aligning to Reference --", config["STATS_FILE"])
    start_time = timeit.default_timer()
    reads = parallelize(Read.align, reads, config["NKERN"])
    print("Alignment took {} s".format(timeit.default_timer() - start_time))

    # FILTER BASED ON ALIGNMENT SCORE (INCL FAILED ALIGNMENTS WITH read.al_score is None)
    reads = filter_alignment_score(reads)

    # FILTER BASED ON UNALIGNED PRIMERS
    # --> require that primers are always aligned without gaps / indels
    # --> do allow mismatches
    if config["TECH"] == "Illumina":
        reads = filter_alignments_with_gaps_in_primers(reads)

    # FINAL STATS
    save_stats("Total reads remaining for analysis: {} ({} %)".format(sum([read.counts for read in reads]), sum([read.counts for read in reads]) * 100 / TOTAL_READS), config["STATS_FILE"])

    # REORDER TRAILING INSERTS TO GUARANTEE REF SEQ MORE TRAILING THAN INSERT SEQ AND CORRECT REF_SPAN
    reads = parallelize(Read.reorder_trailing_inserts, reads, config["NKERN"])
    reads = parallelize(Read.get_ref_span, reads, config["NKERN"])

    # PRINT PASSING ALIGNMENTS
    # create output file directory for alignments print-outs
    needle_dir = os.path.join(config["OUT_DIR"],'out_needle')
    if not os.path.exists(needle_dir):
        os.makedirs(needle_dir)

    for i,read in enumerate(reads):
        reads[i].al_file = 'needle_{}.txt'.format(i)
        print_alignment(reads[i], needle_dir)

    if not reads:
        save_stats("\nNO READS TO PROCESS!", config["STATS_FILE"])
        quit()


    #######################################
    # CALCULATE COVERAGE
    start_time = timeit.default_timer()
    ref_coverage = []
    ref_coverage_frwd = []
    ref_coverage_rev = []
    for coord, bp  in enumerate(config["REF"]):
        spanning_reads = [read for read in reads if coord >= read.ref_span[0] and coord <= read.ref_span[1]]
        spanning_reads_index = flatten_list([read.index for read in spanning_reads])
        ref_coverage.append(len(set(spanning_reads_index)))
        ref_coverage_frwd.append(sum([read.counts for read in reads if coord >= read.ref_span[0] and coord <= read.ref_span[1] and read.sense == 1]))
        ref_coverage_rev.append(sum([read.counts for read in reads if coord >= read.ref_span[0] and coord <= read.ref_span[1] and read.sense == -1]))
    print("Calculating coverage took {} s".format(timeit.default_timer() - start_time))


    #######################################
    # COLLECT INSERTS
    save_stats("\n-- Looking for insertions & ITDs --", config["STATS_FILE"])

    inserts = []
    start_time = timeit.default_timer()
    for read in reads:
        inserts.append(read.get_inserts())
    inserts = flatten_list([insert for insert in inserts if insert is not None])
    print("Collecting inserts took {} s".format(timeit.default_timer() - start_time))
    save_stats("{} insertions were found".format(len(inserts)), config["STATS_FILE"])

    # filter inserts that are actually adapter sequences
    # (instead of trimming adapters in advance)
    if config["TECH"] == "Illumina":
        start_time = timeit.default_timer()
        inserts = filter_sequencing_adapter_artefacts(inserts)
        print("Filtering inserts for adapter sequences took {} s".format(timeit.default_timer() - start_time))

    # add coverage to inserts
    start_time = timeit.default_timer()
    for insert in inserts:
        # add coverage
        # --> be sure to normalize start coord to [0,len(REF)[ first
        # --> negative start (-> 5' trailing_end) will result in
        #     coverage = ref_coverage[-X] which will silently report incorrect coverage!!
        insert = insert.set_sense()
        insert = insert.set_coverage(ref_coverage, ref_coverage_frwd, ref_coverage_rev)
        insert = insert.calc_vaf()
    print("Annotating coverage took {} s".format(timeit.default_timer() - start_time))


    #######################################
    # COLLECT ITDs
    # --> put this in a method and use parallelize to speed things up!
    # --> (can I also do that for reads above when there are possibly multiple inserts per read?) -> yes: return [inserts found] per itd, remove None, flatten list
    itds = []
    start_time = timeit.default_timer()
    for insert in inserts:
        itds.append(insert.get_itd(config))
    itds = [itd for itd in itds if itd is not None]

    inserts = sorted(inserts, key=Insert.get_seq)
    itds = sorted(itds, key=Insert.get_seq)
    print("Collecting ITDs took {} s".format(timeit.default_timer() - start_time))
    save_stats("{} ITDs were found".format(len(itds)), config["STATS_FILE"])


    ########################################
    # MERGE INSERTS
    ins_and_itds = {"insertions": inserts, "itds": itds}

    merged_ins_and_itds = {}
    start_time = timeit.default_timer()
    for type_, inserts_ in ins_and_itds.items():
        merged_ins_and_itds[type_] = get_merged_inserts(inserts_, type_, config)
    print("Merging took {} s".format(timeit.default_timer() - start_time))


    ########################################
    # FILTER INSERTS
    filtered_ins_and_itds = {}
    start_time = timeit.default_timer()
    for type_, inserts_ in merged_ins_and_itds.items():
        filtered_ins_and_itds[type_] = get_hc_inserts(inserts_, type_, "_collapsed-is-same_is-similar_is-close_is-same_trailing_hc", config)
    print("Filtering took {} s".format(timeit.default_timer() - start_time))


########## MAIN ####################
if __name__ == '__main__':

    config = parse_config_from_cmdline(config)
    main(config)
