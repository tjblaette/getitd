import getitd
import random
import itertools

config = {}
config["SAMPLE"] = "testtesttest"
config["OUT_DIR"] = "."
config["OUT_NEEDLE"] = "."

config["MIN_BQS"] = 30
config["REF"] = "GCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGGTTGCCGTCAAAATGCTGAAAG"
config["TECH"] = "ILLUMINA"
config["NKERN"] = 14

config["COST_MATCH"] = 5
config["COST_MISMATCH"] = -15
config["COST_ALIGNED"] = {(c1, c2): getitd.get_alignment_score(c1, c2, config) for c1, c2 in itertools.combinations_with_replacement(["A","T","G","C","Z","N"], 2)}
config["COST_GAPOPEN"] = -36
#config["COST_GAPOPEN"] = -41
config["COST_GAPEXTEND"] = -0.5
config["MIN_SCORE_INSERTS"] = 0.5
config["MIN_SCORE_ALIGNMENTS"] = 0.5

config["MIN_READ_LENGTH"] = 100
config["MIN_READ_COPIES"] = 2
config["MIN_TOTAL_READS"] = 1
config["MIN_UNIQUE_READS"] = 2
config["MIN_VAF"] = 0.001

config["ANNO_FILE"] = "./anno/amplicon_kayser.tsv"
config["ANNO"] = getitd.read_annotation(config["ANNO_FILE"])
config["DOMAINS"] = getitd.get_domains(config["ANNO"])

config["MAX_TRAILING_BP"] = 0

getitd.config = {}
for key in config:
    getitd.config[key] = config[key]

print(config)
print(getitd.config)

def simulate_read(ref, length, sense):
    return getitd.Read(
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
    assert getitd.vaf_to_ar(50) == 1

def test_vaf_to_ar_20():
    assert getitd.vaf_to_ar(20) == 1/4

def test_vaf_to_ar_0():
    assert getitd.vaf_to_ar(0) == 0

def test_vaf_to_ar_100():
    assert getitd.vaf_to_ar(100) == -1



def test_ar_to_vaf_0():
    assert getitd.ar_to_vaf(0) == 0

def test_ar_to_vaf_1():
    assert getitd.ar_to_vaf(1) == 50

def test_ar_to_vaf_one_fourth():
    assert getitd.ar_to_vaf(1/4) == 20



def test_get_number_of_digits_0():
    assert getitd.get_number_of_digits(0) == 1

def test_get_number_of_digits_2():
    assert getitd.get_number_of_digits(2) == 1

def test_get_number_of_digits_1340():
    assert getitd.get_number_of_digits(1340) == 4



def test_flatten_list():
    assert getitd.flatten_list([[1],[2,3]]) == [1,2,3]
    assert getitd.flatten_list([[[1]], [2,3]]) == [[1], 2,3]



def test_connect_bases_same():
    assert getitd.connect_bases("A","A") == "|"
    assert getitd.connect_bases("T","T") == "|"
    assert getitd.connect_bases("C","C") == "|"
    assert getitd.connect_bases("G","G") == "|"

def test_connect_bases_gap():
    assert getitd.connect_bases("A","-") == " "
    assert getitd.connect_bases("T","-") == " "
    assert getitd.connect_bases("-","C") == " "
    assert getitd.connect_bases("-","G") == " "

def test_connect_bases_different():
    assert getitd.connect_bases("A","T") == "."
    assert getitd.connect_bases("T","C") == "."
    assert getitd.connect_bases("C","G") == "."
    assert getitd.connect_bases("G","A") == "."
   
   
   
def test_align_ref():
    read = getitd.Read(config["REF"])
    read_aligned = getitd.Read.align(read, config=config) 
    assert read_aligned.al_seq == config["REF"]
    assert read_aligned.al_ref == config["REF"]
    assert read_aligned.al_score == getitd.get_min_score(read.seq, config["REF"], 1)

def test_align_wt_R1():
    read = getitd.Read(config["REF"][0:250])
    read_aligned = getitd.Read.align(read, config=config) 
    assert read_aligned.al_seq == read.seq + "-" * (len(config["REF"]) - read.length)
    assert read_aligned.al_ref == config["REF"]
    assert read_aligned.al_score == getitd.get_min_score(read.seq, config["REF"], 1)


def test_align_wt_R1_with_5prime_insert():
    insert = "AAAA"
    read = getitd.Read(seq = insert + config["REF"][0:(250 - len(insert))])
    read_aligned = getitd.Read.align(read, config=config) 
    assert read_aligned.al_seq == read.seq + "-" * (len(config["REF"]) - read.length + len(insert))
    assert read_aligned.al_ref == "-" * len(insert) + config["REF"]
    assert read_aligned.al_score == getitd.get_min_score(read.seq[len(insert):], config["REF"], 1)




##############################
##  RUN 1

def test_molm_21bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGG").align(config).get_ref_span().get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["598", "599"]
    itd.annotate_domains(config["DOMAINS"])
    assert itd.domains[-1] == "exon14_JMD_zipperMotif"

def test_molm_21bp_02():
    read = getitd.Read(seq="AAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGG").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["598", "599"]
    itd.annotate_domains(config["DOMAINS"])
    assert itd.domains[-1] == "exon14_JMD_zipperMotif"


def test_1610_9263_90bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 90
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in [".", "612", "613"]
    itd.annotate_domains(config["DOMAINS"])
    assert "TKD1_beta1Sheet" in itd.domains[-1] or "intron14" in itd.domains[-1]




def test_1610_9263_78bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAATGGGCTGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTC").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 78
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "602"
    itd.annotate_domains(config["DOMAINS"])
    assert itd.domains[-1] == "exon14_JMD_zipperMotif"



def test_1610_9181_174bp():
    read = getitd.Read(seq="AACATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCACCTTCTGATT").align(config).get_ref_span()
    read.print()
    inserts = read.get_inserts(config)
    print(inserts)
    assert len(inserts) == 1
    insert = inserts[0]
    insert.print()
    assert (insert.length - 174) <= 4
    itd = insert.get_itd(config)
    assert (itd.fix_trailing_length().length - 174 <= 4)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "618" or itd.insertion_site_protein_as == "619"
    itd.annotate_domains(config["DOMAINS"])
    assert itd.domains[-1] == "exon15_TKD1_nucleotideBindingLoop"


def test_1610_9181_174bp_02():
    read = getitd.Read(seq="GCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCACCTTCTGATT").align(config).get_ref_span()
    read.print()
    inserts = read.get_inserts(config)
    print(inserts)
    assert len(inserts) == 1
    insert = inserts[0]
    insert.print()
    assert (insert.length - 174) <= 4
    itd = insert.get_itd(config)
    assert (itd.fix_trailing_length().length - 174 <= 4)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "618" or itd.insertion_site_protein_as == "619"
    itd.annotate_domains(config["DOMAINS"])
    assert itd.domains[-1] == "exon15_TKD1_nucleotideBindingLoop"



def test_1610_9181_45bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTGGGGTGGAACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTAC").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 45
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "594"
    itd.annotate_domains(config["DOMAINS"])
    assert itd.domains[-1] == "exon14_JMD_zipperMotif"



def test_1610_9181_27bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCCCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAACCAGAAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 27
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "603"
    itd.annotate_domains(config["DOMAINS"])
    assert itd.domains[-1] == "exon14_JMD_zipperMotif"



def test_1610_6948_54bp():
    read = getitd.Read(seq="AAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAAT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 54
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "594"
    itd.annotate_domains(config["DOMAINS"])
    assert itd.domains[-1] == "exon14_JMD_zipperMotif"


def test_1610_8230_9bp():
    read = getitd.Read(seq="AACAATCTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATGGCTTCATTATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 9
    itd = insert.get_itd(config)
    assert itd is None



##############################
##  RUN 2

def test_1610_111_60bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 60
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["613", "."]


def test_1610_111_39bp():
    read = getitd.Read(seq="AGCAATTTAGGTATGAAAGCCAGCTACAGATGGCACAGGCGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 39
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "596"


def test_1610_14_21bp():
    read = getitd.Read(seq="AAGCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["600", "601"]


def test_1610_14_21bp_02():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["600", "601"]


def test_1610_150_93bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 93
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["611", "612"]


def test_1610_150_36bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAGGGCGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 36
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "596"


def test_1610_189_87bp():
    read = getitd.Read(seq="AGAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATCCCACCGGGTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGATTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTG").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 87
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["609", "610"]


def test_1610_189_87bp_02():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATCCCACCGGGTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTG").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 87
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["609", "610"]


def test_1610_189_87bp_03():
    read = getitd.Read(seq="AGAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATCCCACCGGGTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGATTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 87
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["609", "610"]


def test_1610_232_42bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTGGTCCCTTCCTTAGTGAAGGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 42
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "590"


def test_1610_232_60bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAGCCACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCACCATTTCTTTTCCATTGGAAAATCT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 60
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "598"


def test_1610_232_60bp_02():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAGCCACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 60
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "598"


def test_1610_264_198bp():
    read = getitd.Read(seq="ACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGCCACAGGTGACCGGCTCCTCAA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    itd = insert.get_itd(config)
    assert itd is not None
    assert abs(itd.fix_trailing_length().length - 198) <= 1
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["614", "615"]


def test_1610_264_198bp_02():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGCCACAGGTGACCGGCTCCTCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    itd = insert.get_itd(config)
    assert itd is not None
    assert abs(itd.fix_trailing_length().length - 198) <= 1
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["614", "615"]

def test_1610_264_198bp_03():
    read = getitd.Read(seq="AGCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGCACTAGGGACCGGCGCCTTA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.trailing == True
    itd = insert.get_itd(config)
    assert itd is not None
    assert abs(itd.fix_trailing_length().length - 198) <= 1
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["614", "615"]


def test_1610_38_42bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAACCTTTTCCACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGAACTCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 42
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "602"


def test_1610_38_42bp_02():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAACCTTTTCCACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 42
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "602"


def test_1610_38_42bp_03():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAACCTTTTCCACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGAACTCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 42
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "602"


def test_1610_7_72bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCC").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 72
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["606", "607"]


def test_1610_7_72_02():
    read = getitd.Read(seq="AGAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCC").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 72
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["606", "607"]


def test_1610_76_21bp():
    read = getitd.Read(seq="AAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACCCACCATTTGTCTTTCCAGGGAAGG").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["601", "602"]


def test_1610_76_69bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 69
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["596", "597"]


def test_1610_76_21bp_02():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAG").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["601", "602"]


##############################
##  RUN 3

def test_1610_145_27bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 27
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "608"


def test_1610_42_27bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATCCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 27
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "599"

def test_1610_54_36bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 36
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["608", "609"]

def test_1610_89_45bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGACCTCCTCGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTAC").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 45
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "598"



##############################
##  RUN 6/7

def test_pl21_126bp():
    read = getitd.Read(seq="ACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTTCAAATCGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATT").align(config).get_ref_span()
    read.print()
    read.al_file = "pl21.txt"
    read.print_alignment(config)
    inserts = read.get_inserts(config)
    print(inserts)
    assert len(inserts) == 1
    insert = inserts[0]
    insert.print()
    itd = insert.get_itd(config)
    itd.print()
    itd = itd.prep_for_save(config)
    itd.print()
    itd.reads[0].print()
    assert itd is not None
    assert abs(itd.fix_trailing_length().length - 126) <= 1

def maybe_test_pl21_126bp_artificial():
    read = getitd.Read(seq="ACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTTCAAATCGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTXXXXXXGGAGTTTCCAAGAGAAAATT").align(config).get_ref_span()
    read.print()
    read.al_file = "pl21_artificial.txt"
    read.print_alignment(config)
    inserts = read.get_inserts(config)
    print(inserts)
    assert len(inserts) == 1
    insert = inserts[0]
    insert.print()
    itd = insert.get_itd(config)
    itd.print()
    itd = itd.prep_for_save(config)
    itd.print()
    itd.reads[0].print()
    assert itd is not None
    assert abs(itd.fix_trailing_length().length - 126) <= 1


##############################
##  RUN 4/5

def test_1610_109_75bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGGTGACCGGCTCGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGGTTCTGCAGCATTTCTTT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 75
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["600", "601", "602"]


def test_1610_115_24bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGG").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 24
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["601", "602"]


def test_1610_133_66bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAATGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 66
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "604"


def test_1610_15_21bp():
    read = getitd.Read(seq="ACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGG").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "596"

def test_1610_15_180bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGCTCCTCAGATAATGAGTACTTCTAC").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert abs(itd.length - 180) <= 1
    assert itd.insertion_site_protein_as == "613"


def test_1610_164_21bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAG").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 21
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "602"


def test_1610_164_99bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGACGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGC").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 99
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in [".", "613"]


def test_1610_17_39bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 39
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "593"


def test_1610_197_39bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 39
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["599","600"]


def test_1610_197_66bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAATGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 66
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in [".", "613"]


def test_1610_209_54bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 54
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in [".", "613"]


def test_1610_21_57bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 57
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in [".", "613"]


def test_1610_224_48bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGGGGAATTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACG").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 48
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "608"


def test_1610_240_75bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 75
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in [".", "613"]


def test_1610_281_54bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATAGGAGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 54
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["599", "600"]


def test_1610_49_45bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTAC").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 45
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in [".", "613"]


def test_1610_52_54bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 54
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "609"


def test_1610_6_54bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAAGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 54
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["598", "599"]


def test_1610_86_45bp():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCCGACCGGAAAAATGGTCGGTCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTAC").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 45
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "591"


##############################
##  RUN 8

def test_1610_111_39bp_02():
    read = getitd.Read(seq="ACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCAT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 39
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "596"


def test_1610_111_60bp_02():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 60
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["613","."]


def test_1610_189_87bp_04():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATCCCACCGGGTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTG").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 87
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as in ["609", "610"]


def test_1610_232_42bp_02():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTGGTCCCTTCCTTAGTGAAGGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 42
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "590"


def test_1610_232_60bp_03():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAGCCACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCT").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    assert insert.length == 60
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert itd.insertion_site_protein_as == "598"


def test_1610_264_198bp_03():
    read = getitd.Read(seq="AACAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGCCACAGGTGACCGGCTCCTCA").align(config).get_ref_span()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    itd = insert.get_itd(config)
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert abs(itd.length - 198) <= 1
    assert itd.insertion_site_protein_as in ["614", "615"]

def test_1610_264_198bp_04_nonTrailingThoughShouldBeTrailingITD():
    read = getitd.Read(seq="GAAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGCCACAGGTGACCGGCTCCTCAG").align(config).get_ref_span()
    read.print()
    inserts = read.get_inserts(config)
    assert len(inserts) == 1
    insert = inserts[0]
    insert.print()
    itd = insert.get_itd(config)
    itd.print()
    assert itd is not None
    itd = itd.prep_for_save(config)
    itd.set_insertion_site()
    itd.annotate("insertion_site", "protein_as", config)
    assert abs(itd.length - 198) <= 1
    assert itd.insertion_site_protein_as in ["614", "615"]

