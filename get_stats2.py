import Bio.pairwise2 as bio
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')  # required to use matplotlib without X (via ssh + screen)
import matplotlib.pyplot as plt
import sys

print("###\n###\n###")

#######################################
## OPEN FILE

#INPUT_FILE = sys.argv[1]
#f = open(INPUT_FILE,"r")


#######################################
# READ IN ALIGNMENTS

# "all.alignments" file should contain aligned sequences only with different alignments separated by an empty line
next_ = [] # save read and ref sequences of the next alignment
all_reads = [] # save all read sequences in file, each element is one read stored as one string 
all_refs = [] # equivalent to all_reads but for ref sequences -> all_reads and all_refs are matched by index (all_reads[0] was aligned to all_refs[0]
all_files = [] # store the file name of each read's needle alignment to be able to manually check these later
all_scores = [] # store alignment scores for these needle alignments --> filter very poor alignments as these are due to really weird reads somehow passing QC filters

def is_number(string):
  try:
    float(string)
    return True
  except ValueError:
    return False


# loop over file to fill all_reads and all_refs
with open("all.alignments", "r") as f:
  for line in f:
    clean = line.rstrip('\n')
    if '.needle' in clean:
      all_files.append(clean)
    elif is_number(clean):
      all_scores.append(float(clean))
    else:
      # append all sequences of the next alignment to next_
      next_.append(line.rstrip('\n'))
      # once the empty newline is encountered that separates distinct alignments, separate read and ref sequences of next_ and store as one string element of all_reads and all_refs respectively
      if(line == '\n' or line == ''):
        all_reads.append(''.join(next_[0:len(next_):2]))
        all_refs.append(''.join(next_[1:len(next_):2]))
        next_ = []

# readline() returns '' for file end, for-loop doesn't -> compensate by running this code once more (think of a better way later...)
all_reads.append(''.join(next_[0:len(next_):2]))
all_refs.append(''.join(next_[1:len(next_):2]))

# alignments are of unique reads only, read in separately the number of each of these reads in the original (qc-filtered) FASTQ file 
with open("all.readCounts", "r") as counts:
  all_readCounts = [int(readCount) for readCount in counts.read().splitlines()]


# make sure there is a 1:1 matching between these files -> must all have the same length!
assert(len(all_reads) == len(all_refs))
assert(len(all_reads) == len(all_readCounts))
assert(len(all_reads) == len(all_scores))


# filter alignments by their score
max_score = max(all_scores) # should be wt, but doesn't have to be!  --> use len(read
min_score = max_score * 0.5
print("Filtering {} / {} low quality alignments with a score < {}".format(sum(np.array(all_scores) < min_score),len(all_scores),  min_score))
print("--------------------")

all_reads = [read for read,score in zip(all_reads,all_scores) if score >= min_score]
all_refs = [ref for ref,score in zip(all_refs,all_scores) if score >= min_score]
all_readCounts = [readCount for readCount,score in zip(all_readCounts,all_scores) if score >= min_score]
all_files = [file_ for file_,score in zip(all_files,all_scores) if score >= min_score]
all_scores = [score for score in all_scores if score >= min_score]

# again make sure there is a 1:1 matching between these files -> must all have the same length!
assert(len(all_reads) == len(all_refs))
assert(len(all_reads) == len(all_readCounts))
assert(len(all_reads) == len(all_scores))
assert(len(all_reads) == len(all_files))


#######################################
# HELPER FUNCTIONS

# callback function to calculate match and mismatch score for realignment of insert to its read
# insert is masked by 'Z' -> return max penalty (min score of -Inf) to prohibit realignment of insert to itself
def get_alignment_score(char1,char2):
  assert(not (char1 == 'Z' and char2 == 'Z')) # only ever one of the sequences chars are taken from should contain masking letter 'Z', i.e. the read sequence but not the ref
  if char1 == char2:
    return 5
  elif char1 == 'Z' or char2 == 'Z':
    return -np.inf
  else:
    return -4


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
def left_normalize(refn, readn, insert_start, insert_end, i):
    if insert_start > 0 and refn[insert_start -1].lower() == readn[insert_end].lower():
        print("LEFT NORMALIZE: {}".format(i))
        return True
    return False


#######################################
# EXTRACT INSERT SEQUENCE FROM READ

# check each alignment for insert/itd and save index in all_reads/all_refs/all_files to idx, insert/itd length to length and insert/itd start/stop position to start/end dicts based on insert/itd classification
w_ins = {"idx": [], "file": [], "length": [], "start": [], "insert": []}
w_itd_exact = {"idx": [], "file": [], "length": [], "start": [], "tandem2_start": [], "offset": [], "insert": []}
w_itd_nonexact = {"idx": [], "file": [], "length": [], "start": [], "tandem2_start": [], "offset": [], "insert": []}
w_itd_nonexact_fail = {"idx": [], "file": [], "length": [], "start": [], "offset": [], "insert": []}

ref_wt = [base for base in all_refs[0] if base != '-'] 
ref_coverage = np.zeros(len(ref_wt)) # count number of reads covering each bp AND its successor (therefore do not calc coverage for last bp)

ambig_i = []
ambig_als = []
should_left_normalize = 0

trailing5 = []
trailing3 = []

# loop over all alignments, test for presence of an ITD
test = False
#i_of_interest = 10
#test = True
for read,ref,score,counts,filename,i in zip(all_reads, all_refs, all_scores, all_readCounts, all_files, range(len(all_reads))):
	readn = np.array(list(read))
	refn = np.array(list(ref))
	assert(len(readn) == len(refn))
#
	readn_onRef = readn[refn != '-'] ## compare readn_nonIns below
	readn_onRef_covered = np.where(readn_onRef != '-')
	readn_onRef_covered_range = np.arange(np.min(readn_onRef_covered), np.max(readn_onRef_covered)) # do not count last index -> read ending here holds no information on whether or not an ITD starts here --> but what about a read that covers the first 5 bases of an ITD > 5bp, I still wouldn't know...
	ref_coverage[readn_onRef_covered_range] = ref_coverage[readn_onRef_covered_range] + counts
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
	    # if there is an insert  --> require min 6 bp length, in-frame insert and no "N"s within insert
	    if(insert_length >= 6 and insert_length % 3 == 0 and "N" not in readn[insert_idxs]):
		    ins = readn[insert_idxs]
		    insert_start = insert_idxs[0]
		    insert_end = insert_idxs[-1]
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
		    # take the one closest to the insert (should be relevant only for small ITDs that may be contained multiple time within a read
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
		    # save if an exact second tandem of the insert was found
		    if tandem2_start != -1:   # ---> also check that index of second match is sufficiently close to insert! (for exact match and alignment approach!)
			    w_itd_exact["idx"].append(i)
			    w_itd_exact["file"].append(filename)
			    w_itd_exact["length"].append(insert_length)
			    w_itd_exact["start"].append(insert_start_ref)
			    w_itd_exact["tandem2_start"].append(read_to_wt_coord(tandem2_start, refn))
			    w_itd_exact["offset"].append(abs(tandem2_start - insert_start))
			    w_itd_exact["insert"].append(''.join(ins))
		    else:
			    # otherwise search for sufficiently similar (> 90 % bases mapped) second tandem by realignment of the insert within the remainder of the read
			    max_score = len(ins) * 5  # +5 score per matching base
			    min_score = max_score * 0.9
			    # arguments: seq1, seq2, match-score, mismatch-score, gapopen-score, gapextend-score --> match/mismatch from needle default (/usr/share/EMBOSS/data/EDNAFULL), gap as passed to needle in my script
			    # output: list of optimal alignments, each a list of seq1, seq2, score, start-idx, end-idx 
			    alignments = bio.align.localcs(''.join(ins), ''.join(readn_maskedIns), get_alignment_score, -20, -0.05)
			    # filter alignments where insert cannot be realigned in one piece
			    alignments = [al for al in alignments if integral_insert_realignment(al[0],21)]
			    alignment_score = None
			    #
			    if alignments == []:
				    alignment_score = -1
			    else:
				    alignment = alignments[0]
				    alignment_score, alignment_start, alignment_end = alignment[2:5]
#			
			    if alignment_score >= min_score:
				    w_itd_nonexact["idx"].append(i)
				    w_itd_nonexact["file"].append(filename)
				    w_itd_nonexact["length"].append(insert_length)
				    w_itd_nonexact["start"].append(insert_start_ref)
				    w_itd_nonexact["tandem2_start"].append(read_to_wt_coord(alignment_start, refn))
				    w_itd_nonexact["offset"].append(abs(alignment_start - insert_start))
				    w_itd_nonexact["insert"].append(''.join(ins))
			    else:
				    w_itd_nonexact_fail["idx"].append(i)
				    w_itd_nonexact_fail["file"].append(filename)
				    w_itd_nonexact_fail["length"].append(insert_length)
				    w_itd_nonexact_fail["start"].append(insert_start_ref)
				    w_itd_nonexact_fail["insert"].append(''.join(ins))
				    #print(bio.format_alignment(*alignment))
				    if insert_length > insert_start:
					    trailing5.append(i)
				    if insert_length > len(readn) -1 -insert_end:
					    trailing3.append(i)
			    if len(alignments) > 1:
				    ambig_i.append(i)
				    ambig_als.append(alignments)
	if test == True and i == i_of_interest:
		break

# print number of ambiguous alignments (to see if this is sth I need to handle or not)
print("There were {} inserts that generated ambiguous alignments.".format(len(ambig_i)))
print("--------------------")

print("There were {} inserts whose alignment should have been left normalized.".format(should_left_normalize))
print("--------------------")

#w_itd = {"idx": w_itd_exact["idx"] + w_itd_nonexact["idx"], "file": w_itd_exact["file"] + w_itd_nonexact["file"], "length": w_itd_exact["length"] + w_itd_nonexact["length"], "tandem2_start": w_itd_exact["tandem2_start"] + w_itd_nonexact["tandem2_start"], "insert": w_itd_exact["insert"] + w_itd_nonexact["insert"]}
w_itd = {"idx": w_itd_exact["idx"] + w_itd_nonexact["idx"], "file": w_itd_exact["file"] + w_itd_nonexact["file"], "length": w_itd_exact["length"] + w_itd_nonexact["length"], "tandem2_start": w_itd_exact["tandem2_start"] + w_itd_nonexact["tandem2_start"], "insert": w_itd_exact["insert"] + w_itd_nonexact["insert"], 'offset': w_itd_exact["offset"] + w_itd_nonexact["offset"]}


########################################
# COLLECT AND COLLAPSE ITDs

df_itd = pd.DataFrame(w_itd)
df_itd["counts"] = [all_readCounts[i] for i in df_itd["idx"]]

# require that insert offset == insert length --> means they are adjacent   -> test this much earlier already, when saving inserts above! (left shift + extend first)
#df_itd = pd.concat([df_itd.ix[df_itd["offset"] == df_itd["length"]], df_itd.ix[df_itd["offset"] != df_itd["length"]]])
#df_itd = pd.concat([df_itd.ix[df_itd["offset"] != df_itd["length"]], df_itd.ix[df_itd["offset"] == df_itd["length"]]])
df_itd = df_itd.ix[df_itd["offset"] == df_itd["length"]][['file', 'idx', 'insert', 'length', 'tandem2_start', 'counts']]

df_itd_grouped = df_itd.groupby(by=["length","tandem2_start","insert"], as_index=False).sum()
df_itd_grouped["ref_coverage"] = [ref_coverage[pos] for pos in df_itd_grouped["tandem2_start"]]
df_itd_grouped["vaf"] = (df_itd_grouped["counts"]/df_itd_grouped["ref_coverage"] * 100).round(5)
df_itd_grouped["file"] = np.zeros(len(df_itd_grouped)) 
df_itd_grouped["counts_each"] = np.zeros(len(df_itd_grouped)) 
df_itd_grouped[["idx","file","counts_each"]] = df_itd_grouped[["idx","file","counts_each"]].astype("object")

for i in range(len(df_itd_grouped)):
	this_itd = df_itd[np.array(df_itd["length"] == df_itd_grouped.ix[i,"length"]) * np.array(df_itd["tandem2_start"] == df_itd_grouped.ix[i,"tandem2_start"]) * np.array(df_itd["insert"] == df_itd_grouped.ix[i,"insert"])]
	df_itd_grouped.set_value(i,"idx",this_itd["idx"].tolist())
	df_itd_grouped.set_value(i,"file",this_itd["file"].tolist())
	df_itd_grouped.set_value(i,"counts_each",[np.int(x) for x in this_itd["counts"]])

# check that sum of "counts_each" (= read counts of each unique read) equals total counts in "counts"
assert([sum(x) for x in df_itd_grouped["counts_each"]] == [int(x) for x in df_itd_grouped["counts"]])

df_itd_grouped[['length', 'tandem2_start', 'vaf', 'ref_coverage', 'counts', 'counts_each', 'file']].to_csv("flt3_itds.csv", index=False)



########################################
# COLLAPSE ITDs #2 
# --> align inserts of same length and tandem2_start, collapse if they are sufficiently similar


df_itd_collapsed = pd.DataFrame(columns=['length', 'tandem2_start', 'insert', 'idx', 'file', 'counts', 'ref_coverage', 'vaf', 'counts_each'])
df_itd_collapsed["idx"] = []
df_itd_collapsed["file"] = []
df_itd_collapsed["counts_each"] = []

for length in set(df_itd_grouped["length"]):
    this_df_length = df_itd_grouped.ix[df_itd_grouped["length"] == length]
    #
    for tandem2_start in set(this_df_length["tandem2_start"]):
        this_df = this_df_length.ix[this_df_length["tandem2_start"] == tandem2_start]
        #
        max_score = length * 5
        min_score = max_score * 0.5
        #
        for i in range(this_df.shape[0]):
            i_idx = this_df.index[i]
            this_ins = this_df["insert"][i_idx]
            collapsed = False
            #
            for ii,iirow in df_itd_collapsed[::-1].iterrows(): #[::-1] to reverse df and speed up pos alignment
                other_ins = iirow["insert"]
                #
                alignment = bio.align.globalcs(this_ins, other_ins, get_alignment_score, -20, -0.05)[0]
                alignment_score = alignment[2]
                #
                if alignment_score >= min_score:
                    collapsed = True
                    # collapse
                    # add together some statistics
                    df_itd_collapsed.at[ii,"counts"] = df_itd_collapsed["counts"][ii] + this_df["counts"][i_idx]
                    df_itd_collapsed.at[ii, "vaf"] = df_itd_collapsed["vaf"][ii] + this_df["vaf"][i_idx]
                    #
                    # pick one or the other for the others OR keep both but in specific order (first list for picked insert) -> go for the most abundant one (or the one closest to reference?!)
                    if this_df["counts"][i_idx] > df_itd_collapsed["counts"][ii]:
                        df_itd_collapsed.at[ii, "insert"] = this_df["insert"][i_idx]
                        df_itd_collapsed.at[ii, "ref_coverage"] = this_df["ref_coverage"][i_idx]
                        #
                        df_itd_collapsed.at[ii, "idx"] = this_df["idx"][i_idx] + df_itd_collapsed["idx"][ii]
                        df_itd_collapsed.at[ii, "file"] = this_df["file"][i_idx] + df_itd_collapsed["file"][ii]
                        df_itd_collapsed.at[ii, "counts_each"] = this_df["counts_each"][i_idx] + df_itd_collapsed["counts_each"][ii]
                    else:
                        df_itd_collapsed.at[ii, "idx"] = df_itd_collapsed["idx"][ii] + this_df["idx"][i_idx]
                        df_itd_collapsed.at[ii, "file"] = df_itd_collapsed["file"][ii] + this_df["file"][i_idx]
                        df_itd_collapsed.at[ii, "counts_each"] = df_itd_collapsed["counts_each"][ii] + this_df["counts_each"][i_idx]
                    break
            #
            if not collapsed:
                df_itd_collapsed = df_itd_collapsed.append(this_df.ix[i_idx], ignore_index=True)


# check that sum of "counts_each" (= read counts of each unique read) equals total counts in "counts"
assert([sum(x) for x in df_itd_collapsed["counts_each"]] == [int(x) for x in df_itd_collapsed["counts"]])

df_itd_collapsed[['length', 'tandem2_start', 'vaf', 'ref_coverage', 'counts', 'insert']].to_csv("flt3_itds_collapsed.csv", index=False)



########################################
# COLLAPSE ITDs #3
# --> align inserts of same length with masked flanking sequence, collapse if they are sufficiently similar



########################################
# APPLY SOME MORE FILTERS

# filter low support ITDs
#print(df_itd_collapsed)
#df_itd_collapsed = df_itd_collapsed.ix[df_itd_collapsed["counts"] > 10]

# filter duplications that change codons instead of simply duplicating them (check with other samples, if this is always a valid filter & ask the others)
#print(df_itd_collapsed)
df_itd_collapsed_noShift = df_itd_collapsed.ix[(df_itd_collapsed["tandem2_start"] +1 ) % 3 == 0]

df_itd_collapsed_noShift[['length', 'tandem2_start', 'vaf', 'ref_coverage', 'counts', 'insert']].to_csv("flt3_itds_collapsed_noShift.csv", index=False)


########################################
# FIND MOST ABUNDANT CLONE PER ITD LENGTH

df_itd_maxClone = df_itd_collapsed.groupby(by=["length"], as_index=False).max()[["length","counts"]]
df_itd_maxClone["tandem2_start"] = [list(df_itd_collapsed.ix[np.array(df_itd_collapsed["length"] == this_length) * np.array(df_itd_collapsed["counts"] == max_counts)]["tandem2_start"]) for this_length,max_counts in zip(df_itd_maxClone["length"],df_itd_maxClone["counts"])]
df_itd_maxClone["vaf"] = [list(df_itd_collapsed.ix[np.array(df_itd_collapsed["length"] == this_length) * np.array(df_itd_collapsed["counts"] == max_counts)]["vaf"]) for this_length,max_counts in zip(df_itd_maxClone["length"],df_itd_maxClone["counts"])]
df_itd_maxClone.to_csv("flt3_itds_mostFrequentClonePerLength.csv", index=False)



########################################
# PRINT FILENAMES OF EACH CATEGORY TO FILE
			
out = open("flt3_insertions.txt","w")
out.write("\n".join(np.array(all_files)[w_ins["idx"]]) + "\n")
out.close()

out = open("flt3_itd_exact.txt","w")
out.write("\n".join(np.array(all_files)[w_itd_exact["idx"]]) + "\n")
out.close()

out = open("flt3_itd_nonexact.txt","w")
out.write("\n".join(np.array(all_files)[w_itd_nonexact["idx"]]) + "\n")
out.close()

out = open("flt3_itd_nonexact_fail.txt","w")
out.write("\n".join(np.array(all_files)[w_itd_nonexact_fail["idx"]]) + "\n")
out.close()


########################################
# PRINT SUMMARY STATISTICS on the number of reads in each category

print("Unique reads supporting each type of insert -> NOT number of distinct inserts!")
print("Insertions: {}".format(len(w_ins["idx"])))
print("Single exact ITD: {}".format(len(w_itd_exact["idx"])))
print("Single non-exact ITD: {}".format(len(w_itd_nonexact["idx"])))
print("Single insertion failed alignment: {}".format(len(w_itd_nonexact_fail["idx"])))
print("--------------------")


#########################################
# COLLECT DATA TO PLOT MOLM DILUTION SERIES' VAF AND READ COUNTS OF KNOWN ITD
# ----> LATER EXTEND THIS TO SUPPORT PLOTTING OF ANY/ALL KNOWN ITDs SUPPLIED?!

def get_known(df, klength, kstart):
	kstat = df.ix[np.bitwise_and(df["length"] == klength, df["tandem2_start"] == kstart)][["counts","vaf"]]
	kindex = "_".join([str(x) for x in [klength,kstart]])
#	
	known = pd.DataFrame({"itd": kindex, "counts": kstat["counts"], "vaf": kstat["vaf"]})[["itd","counts","vaf"]]
	print(known)
	known.to_csv("known_itds.csv", index=False)
	


# get known ITD counts for MOLM14 cellline
get_known(df_itd_collapsed,21,71) 




########################################
# PLOT SEPARATE HISTOGRAMS FOR COUNTS OF INSERT/ITD LENGTH/START/STOP/START-STOP PER CATEGORY

def plot_hist(title_, xlab_, top=True, bottom=True):
	ax.bar(counts_table_value,vafs, align='center', color="lightgrey")
	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	ax.tick_params(direction='out')
	ax2.tick_params(direction='out')
	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax2.spines['left'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.spines['bottom'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	# Only show ticks on the left and bottom spines
	ax.yaxis.set_ticks_position('left')
	ax2.yaxis.set_ticks_position('right')
	# set axes limits
	vaf_ylim = 0.001 * 10**(i-1)
	counts_ylim = vaf_ylim/100 * sum(all_readCounts)
	ax.set_ylim([0, vaf_ylim])
	ax2.set_ylim([0, counts_ylim])
	ax.set_xlim([0,len(counts_table)])
	# set axes labels and tick width
	ax.set_ylabel("% VAF")
	ax2.set_ylabel("Counts")
	ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
	ax.yaxis.major.formatter._useMathText = True
	ax2.yaxis.major.formatter._useMathText = True
	if(top):
		ax.set_title(title_)
	if(bottom):
		ax.xaxis.set_ticks_position('bottom')
		ax.set_xlabel(xlab_)
	else:
		plt.tick_params(axis='x',which='both',bottom='on',top='off',labelbottom='on') 
		# changes apply to the x-axis, both major and minor ticks are affected, ticks along the bottom edge are on, ticks along the top edge are off, labels along the bottom edge are on

for ins_type,ins_filename,title_ in zip([w_itd_exact, w_itd_nonexact, w_itd],["w_itd_exact", "w_itd_nonexact", "w_itd"],["ITDs - exact matching", "ITDs - non-exact matching", "ITDs - all"]):
	for stat,xlab_ in zip(["length","tandem2_start"], ["Insert length (bp)","Tandem2 start position"]):
		if not ins_type[stat]: # prevent failure for empty lists (happened for small test sets)
			continue
		
		counts_table = None
		counts_table_value = None
		counts_table_value_toIndexOffset = 0
		if stat == "length":
			counts_table = np.zeros(max(ins_type["length"]) +1)
			counts_table_value = np.arange(0,len(counts_table))
		elif stat == "tandem2_start":
			min_val = min(0,min(ins_type[stat]))
			max_val = max(ins_type[stat])
			counts_table = np.zeros(max_val + abs(min_val) + 1)
			counts_table_value = np.arange(min_val, max_val +1) 
			assert(len(counts_table) == len(counts_table_value))
			counts_table_value_toIndexOffset = abs(min_val)
		assert(counts_table is not None) # should be assigned sth
		
		# sum up counts of insert/itds with the same lengths/statistic (length, start or stop)
		for i in range(len(ins_type[stat])):
			i_stat = ins_type[stat][i] + counts_table_value_toIndexOffset
			i_idx = ins_type["idx"][i]
			counts_table[i_stat] = counts_table[i_stat] + all_readCounts[i_idx]

		vafs = counts_table/sum(all_readCounts) * 100   # in percent
		# PLOT
		fig = plt.figure(figsize=(8.27, 11.69)) #A4 DIN size
		n_plots = 6
		for i in range(1,n_plots+1):
			ax = fig.add_subplot(n_plots,1,n_plots -i +1)
			ax2 = ax.twinx() # Create another axis that shares the same x-axis as ax.	
			top_plt = False
			bottom_plt = False
			if(i == n_plots):
				top_plt = True
			if(i == 1):
				bottom_plt = True
			plot_hist(title_, xlab_, top_plt, bottom_plt)	
	
		fig.tight_layout() # required to not have subplots overlap
		# save counts as CSV table and as histogram 
		pd.DataFrame(counts_table, index=counts_table_value, columns=["counts"]).to_csv("table_" + stat + "_" + ins_filename + ".csv", index=False)
		plt.savefig("plot_" + stat + "_" + ins_filename + ".pdf")




