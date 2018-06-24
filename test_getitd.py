from getitd import *
import random

config["SAMPLE"] = "testtesttest"
config["MIN_BQS"] = 30
config["REF"] = "GCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGGTTGCCGTCAAAATGCTGAAAG"
config["TECH"] = "ILLUMINA"
config["NKERN"] = 14

config["COST_MATCH"] = 5
config["COST_MISMATCH"] = -10
config["COST_GAPOPEN"] = -20
config["COST_GAPEXTEND"] = -0.5
config["MIN_SCORE_INSERTS"] = 0.5
config["MIN_SCORE_ALIGNMENTS"] = 0.5

config["MIN_READ_LENGTH"] = 100
config["MIN_READ_COPIES"] = 2
config["MIN_TOTAL_READS"] = 1
config["MIN_UNIQUE_READS"] = 2
config["MIN_VAF"] = 0.001


def simulate_read(ref, length, sense):
    return Read(
            seq = config["REF"][::sense][0:length][::sense],
            sense = sense
            )

def add_substitution(read):
    site = random.randrange(read.length)
    read.seq[site] = random.choice(["A","T","C","G"])
    return read
    
def add_itd(seq, site, length):
    print(seq)
    print(seq[0:site] + 2* seq[site:site+length] + seq[site+length:])



def test_vaf_to_ar_50():
    assert vaf_to_ar(50) == 1

def test_vaf_to_ar_20():
    assert vaf_to_ar(20) == 1/4

def test_vaf_to_ar_0():
    assert vaf_to_ar(0) == 0

def test_vaf_to_ar_100():
    assert vaf_to_ar(100) == -1



def test_ar_to_vaf_0():
    assert ar_to_vaf(0) == 0

def test_ar_to_vaf_1():
    assert ar_to_vaf(1) == 50

def test_ar_to_vaf_one_fourth():
    assert ar_to_vaf(1/4) == 20



def test_get_number_of_digits_0():
    assert get_number_of_digits(0) == 1

def test_get_number_of_digits_2():
    assert get_number_of_digits(2) == 1

def test_get_number_of_digits_1340():
    assert get_number_of_digits(1340) == 4



def test_flatten_list():
    assert flatten_list([[1],[2,3]]) == [1,2,3]
    assert flatten_list([[[1]], [2,3]]) == [[1], 2,3]



def test_connect_bases_same():
    assert connect_bases("A","A") == "|"
    assert connect_bases("T","T") == "|"
    assert connect_bases("C","C") == "|"
    assert connect_bases("G","G") == "|"

def test_connect_bases_gap():
    assert connect_bases("A","-") == " "
    assert connect_bases("T","-") == " "
    assert connect_bases("-","C") == " "
    assert connect_bases("-","G") == " "

def test_connect_bases_different():
    assert connect_bases("A","T") == "."
    assert connect_bases("T","C") == "."
    assert connect_bases("C","G") == "."
    assert connect_bases("G","A") == "."
   
   
   
def test_align_ref():
    read = Read(config["REF"])
    read_aligned = Read.align(read, config=config) 
    assert read_aligned.al_seq == config["REF"]
    assert read_aligned.al_ref == config["REF"]
    assert read_aligned.al_score == len(config["REF"]) * config["COST_MATCH"]

def test_align_wt_R1():
    read = Read(config["REF"][0:250])
    read_aligned = Read.align(read, config=config) 
    assert read_aligned.al_seq == read.seq + "-" * (len(config["REF"]) - read.length)
    assert read_aligned.al_ref == config["REF"]
    assert read_aligned.al_score == read.length * config["COST_MATCH"]


def test_align_wt_R1_with_5prime_insert():
    insert = "AA"
    read = Read(seq = insert + config["REF"][0:(250 - len(insert))])
    read_aligned = Read.align(read, config=config) 
    assert read_aligned.al_seq == read.seq + "-" * (len(config["REF"]) - read.length + len(insert))
    assert read_aligned.al_ref == "--" + config["REF"]
    assert read_aligned.al_score == (len(read.seq) - len(insert)) * config["COST_MATCH"]

def test_molm_21bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGG").align().get_reference_range_covered().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd()
    assert itd is not None


def test_molm_21bp_02():
    read = Read(seq="AAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGG").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd()
    assert itd is not None



def test_1610_9263_90bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTT").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 90
    itd = insert.get_itd()
    assert itd is not None



def test_1610_9263_78bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAATGGGCTGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTC").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 78
    itd = insert.get_itd()
    assert itd is not None


def test_1610_9181_174():
    read = Read(seq="AACATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCACCTTCTGATT").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert (insert.length - 174) <= 6
    itd = insert.get_itd()
    assert (itd.fix_trailing_length().length - 174 <= 6)
    assert itd is not None


def test_1610_9181_45bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTGGGGTGGAACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTAC").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 45
    itd = insert.get_itd()
    assert itd is not None


def test_1610_9181_27bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCCCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAACCAGAAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 27
    itd = insert.get_itd()
    assert itd is not None


def test_1610_6948_54bp():
    read = Read(seq="AAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAAT").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 54
    itd = insert.get_itd()
    assert itd is not None


def test_1610_8230_9bp():
    read = Read(seq="AACAATCTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATGGCTTCATTATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 9
    itd = insert.get_itd()
    assert itd is None


def test_1610_111_60bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCT").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 60
    itd = insert.get_itd()
    assert itd is not None


def test_1610_111_39bp():
    read = Read(seq="AGCAATTTAGGTATGAAAGCCAGCTACAGATGGCACAGGCGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 39
    itd = insert.get_itd()
    assert itd is not None


def test_1610_14_21bp():
    read = Read(seq="AAGCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd()
    assert itd is not None


def test_1610_14_21bp_02():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd()
    assert itd is not None


def test_1610_150_93bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 93
    itd = insert.get_itd()
    assert itd is not None


def test_1610_150_36bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAGGGCGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTT").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 36
    itd = insert.get_itd()
    assert itd is not None


def test_1610_189_87bp():
    read = Read(seq="AGAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATCCCACCGGGTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGATTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTG").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 87
    itd = insert.get_itd()
    assert itd is not None


def test_1610_189_87bp_02():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATCCCACCGGGTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTG").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 87
    itd = insert.get_itd()
    assert itd is not None


def test_1610_189_87bp_03():
    read = Read(seq="AGAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATCCCACCGGGTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGATTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 87
    itd = insert.get_itd()
    assert itd is not None


def test_1610_232_42bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTGGTCCCTTCCTTAGTGAAGGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 42
    itd = insert.get_itd()
    assert itd is not None


def test_1610_232_60bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAGCCACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCACCATTTCTTTTCCATTGGAAAATCT").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 60
    itd = insert.get_itd()
    assert itd is not None


def test_1610_232_60bp_02():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAGCCACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCT").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 60
    itd = insert.get_itd()
    assert itd is not None


def test_1610_264_198bp():
    read = Read(seq="ACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGCCACAGGTGACCGGCTCCTCAA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    itd = insert.get_itd()
    assert itd is not None
    assert abs(itd.fix_trailing_length().length - 198) <= 6


def test_1610_264_198bp_02():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGCCACAGGTGACCGGCTCCTCA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    itd = insert.get_itd()
    assert itd is not None
    assert abs(itd.fix_trailing_length().length - 198) <= 6


def test_1610_38_42bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAACCTTTTCCACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGAACTCA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 42
    itd = insert.get_itd()
    assert itd is not None


def test_1610_38_42bp_02():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAACCTTTTCCACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 42
    itd = insert.get_itd()
    assert itd is not None


def test_1610_38_42bp_03():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAACCTTTTCCACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGAACTCA").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 42
    itd = insert.get_itd()
    assert itd is not None


def test_1610_7_72bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCC").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 72
    itd = insert.get_itd()
    assert itd is not None


def test_1610_7_72_02():
    read = Read(seq="AGAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCC").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 72
    itd = insert.get_itd()
    assert itd is not None


def test_1610_76_21bp():
    read = Read(seq="AAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACCCACCATTTGTCTTTCCAGGGAAGG").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd()
    assert itd is not None


def test_1610_76_69bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATT").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 69
    itd = insert.get_itd()
    assert itd is not None


def test_1610_76_21bp():
    read = Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAG").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd()
    assert itd is not None


def test_pl21_126bp():
    read = Read(seq="ACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTTCAAATCGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATT").align().get_reference_range_covered()
    inserts = read.get_inserts()
    assert len(inserts) == 1
    insert = inserts[0]
    itd = insert.get_itd()
    assert itd is not None
    assert abs(itd.fix_trailing_length().length - 126) <= 3
