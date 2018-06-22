from getitd import *

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
    REF = "GCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGGTTGCCGTCAAAATGCTGAAAG"
    COST_MATCH = 5
    read = Read(REF)
    read_aligned = Read.align(read, config=config) 
    assert read_aligned.al_seq == REF
    assert read_aligned.al_ref == REF
    assert read_aligned.al_score == len(REF) * COST_MATCH
