import Bio.pairwise2 as bio
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')  # required to use matplotlib without X (via ssh + screen)
import matplotlib.pyplot as plt
import sys

#######################################
## OPEN FILE

#INPUT_FILE = sys.argv[1]
#f = open(INPUT_FILE,"r")

# file should contain aligned sequences only with different alignments separated by an empty line
f = open("./all.alignments","r")
#f = open("./test.alignments","r")

# alignments are of unique reads only, read in separately the number of each of these reads in the original (qc-filtered) FASTQ file 
#cfile = open("all.readCounts_needleOrder","r")
cfile = open("all.readCounts","r")
#cfile = open("test.readCounts_needleOrder","r")
all_readCounts = [int(readCount) for readCount in cfile.read().splitlines()]

#######################################
# READ IN ALIGNMENTS

next_ = [] # save read and ref sequences of the next alignment
all_reads = [] # save all read sequences in f, each element is one read stored as one string 
all_refs = [] # equivalent to all_reads but for ref sequences -> all_reads and all_refs are matched by index (all_reads[0] was aligned to all_refs[0]
all_files = [] # store the file name of each read's needle alignment


# loop over file to fill all_reads and all_refs
for line in f:
	if '.needle' in line:
		all_files.append(line.rstrip('\n'))
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

# close the file!
f.close()
cfile.close()

# make sure there is a 1:1 matching between these files -> must all have the same length!
assert(len(all_reads) == len(all_refs))
assert(len(all_reads) == len(all_readCounts))


#######################################
# EXTRACT INSERT SEQUENCE FROM READ


# check each alignment for insert/itd and save index in all_reads/all_refs/all_files to idx, insert/itd length to length and insert/itd start/stop position to start/end dicts based on insert/itd classification
w_ins = {"idx": [], "length": [], "start": []}
w_itd_exact = {"idx": [], "length": [], "start": [], "tandem2_start": []}
w_itd_nonexact = {"idx": [], "length": [], "start": [], "tandem2_start": []}
w_itd_nonexact_fail = {"idx": [], "length": [], "start": []}

ref_wt = [base for base in all_refs[0] if base != '-'] 
ref_coverage = np.zeros(len(ref_wt)) # count number of reads covering each bp AND its successor (therefore do not calc coverage for last bp)

# loop over all alignments, test for presence of an ITD
for read,ref,counts,i in zip(all_reads, all_refs, all_readCounts, range(len(all_reads))):
#this_i = 12432
#if True:
#	read = all_reads[this_i]
#	ref = all_refs[this_i]
#	counts = all_readCounts[this_i]
#	i = this_i
	readn = np.array(list(read))
	refn = np.array(list(ref))
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
		    insert_start = insert_idxs[0]
		    insert_end = insert_idxs[-1]
#
		    # relative to the reference, get coord of the first WT base before insert	
		    insert_start_ref = insert_start - sum(refn[0:insert_start] == '-')  
		    if insert_start == 0: 
			    insert_start_ref = insert_start_ref - insert_length
#	    
		    w_ins["idx"].append(i)
		    w_ins["length"].append(insert_length)
		    w_ins["start"].append(insert_start_ref)
#			
		    # check whether the insert is contained within non-insert read a second time -> that makes it an ITD!
		    readn_nonIns = np.delete(readn,insert_idxs)
		    ins = readn[insert_idxs]
#
		    # search for nearest tandem before and after ITD
		    tandem2_after = ''.join(ref_wt).find(''.join(ins).lower(), insert_start_ref,len(ref_wt))
		    tandem2_before = ''.join(reversed(ref_wt)).find(''.join(reversed(ins)).lower(), len(ref_wt) -1 -insert_start_ref, len(ref_wt))
		    tandem2_before_converted = len(ref_wt) -1 -tandem2_before -insert_length +1  # convert coords back from reverse to forward sense
#
		    # take the one closest to the insert (should be relevant only for small ITDs that may be contained multiple time within a read
		    tandem2_start_ref = None
		    if tandem2_after == -1 and tandem2_before == -1:
			    tandem2_start_ref = -1 # not found --> no itd present
		    elif tandem2_after != -1:
			    tandem2_start_ref = tandem2_after
		    elif tandem2_before != -1:
			    tandem2_start_ref = tandem2_before_converted
		    elif tandem2_after < tandem2_before:
			    tandem2_start_ref = tandem2_after
		    elif tandem2_before < tandem2_after:
			    tandem2_start_ref = tandem2_before_converted
		    assert tandem2_start_ref is not None  # should be assigned something!
#
		    # save if an exact second tandem of the insert was found
		    if tandem2_start_ref != -1:   # ---> also check that index of second match is sufficiently close to insert! (for exact match and alignment approach!)
			    w_itd_exact["idx"].append(i)
			    w_itd_exact["length"].append(insert_length)
			    w_itd_exact["start"].append(insert_start_ref)
			    w_itd_exact["tandem2_start"].append(tandem2_start_ref)
		    else:
			    # otherwise search for sufficiently similar (> 90 % bases mapped) second tandem by realignment of the insert within the remainder of the read
			    max_score = len(ins) * 5  # +5 score per matching base
			    min_score = max_score * 0.9
			    # arguments: seq1, seq2, match-score, mismatch-score, gapopen-score, gapextend-score --> match/mismatch from needle default (/usr/share/EMBOSS/data/EDNAFULL), gap as passed to needle in my script
			    # output: list of optimal alignments, each a list of seq1, seq2, score, start-idx, end-idx 
			    alignment = bio.align.localms(''.join(ins), ''.join(readn_nonIns), 5, -4, -20, -0.05)[0]
			    alignment_score, alignment_start, alignment_end = alignment[2:5]
#			
			    if alignment_score >= min_score:
				    w_itd_nonexact["idx"].append(i)
				    w_itd_nonexact["length"].append(insert_length)
				    w_itd_nonexact["start"].append(insert_start_ref)
				    w_itd_nonexact["tandem2_start"].append(alignment_start)
			    else:
				    w_itd_nonexact_fail["idx"].append(i)
				    w_itd_nonexact_fail["length"].append(insert_length)
				    w_itd_nonexact_fail["start"].append(insert_start_ref)
				    #print(bio.format_alignment(*alignment))

w_itd = {"idx": w_itd_exact["idx"] + w_itd_nonexact["idx"], "length": w_itd_exact["length"] + w_itd_nonexact["length"], "start": w_itd_exact["start"] + w_itd_nonexact["start"], "tandem2_start": w_itd_exact["tandem2_start"] + w_itd_nonexact["tandem2_start"]}
print("-----------------------------------")

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

df_itd = pd.DataFrame(w_itd)
df_itd["counts"] = np.zeros(len(df_itd))

for ins_type,ins_filename,title_ in zip([w_ins, w_itd_exact, w_itd_nonexact, w_itd, w_itd_nonexact_fail],["w_ins", "w_itd_exact", "w_itd_nonexact", "w_itd", "w_itd_nonexact_fail"],["Single insertions", "Single ITDs - exact matching", "Single ITDs - non-exact matching", "Single ITDs - all", "Single insertions - failed ITD matching"]):
	#print(ins_type)
	for stat,xlab_ in zip(["length","start"], ["Insert length (bp)","Insert start position"]):
		#print(stat)
		if not ins_type[stat]: # prevent failure for empty lists (happened for small test sets)
			continue
		
		counts_table = None
		counts_table_value = None
		counts_table_value_toIndexOffset = 0
		if stat == "length":
			counts_table = np.zeros(max(ins_type["length"]) +1)
			counts_table_value = np.arange(0,len(counts_table))
		elif stat == "start":
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
			if((ins_filename=="w_itd_exact" or ins_filename=="w_itd_nonexact") and stat=="length"):
				df_itd.ix[df_itd["idx"] == i_idx, "counts"] = all_readCounts[i_idx]

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
		pd.DataFrame(counts_table, index=counts_table_value, columns=["counts"]).to_csv("table_" + stat + "_" + ins_filename + ".csv")
		plt.savefig("plot_" + stat + "_" + ins_filename + ".pdf")


########################################
# COLLECT AND SAVE TO FILE ALL ITDs DETECTED ALONG WITH LENGTH, START, STOP, TOTAL COUNTS AND INDICES AND COUNTS OF INDIVIDUAL UNIQUE READS

df_itd_grouped = df_itd.groupby(by=["length","start","tandem2_start"], as_index=False).sum()
df_itd_grouped["counts_each"] = np.zeros(len(df_itd_grouped))
df_itd_grouped[["idx","counts_each"]] = df_itd_grouped[["idx","counts_each"]].astype("object")

for i in range(len(df_itd_grouped)):
	this_itd = df_itd[np.array(df_itd["length"] == df_itd_grouped.ix[i,"length"]) * np.array(df_itd["start"] == df_itd_grouped.ix[i,"start"]) * np.array(df_itd["tandem2_start"] == df_itd_grouped.ix[i,"tandem2_start"])]
	
	df_itd_grouped.set_value(i,"idx",np.array(all_files)[this_itd["idx"]].tolist())
	df_itd_grouped.set_value(i,"counts_each",[np.int(x) for x in this_itd["counts"]])

# check that sum of "counts_each" (= read counts of each unique read) equals total counts in "counts"
assert([sum(x) for x in df_itd_grouped["counts_each"]] == [int(x) for x in df_itd_grouped["counts"]])
df_itd_grouped.to_csv("flt3_itds.csv")

df_itd_maxClone = df_itd_grouped.groupby(by=["length"], as_index=False).max()[["length","counts"]]
df_itd_maxClone["vaf"] = (df_itd_maxClone["counts"]/sum(all_readCounts) * 100).round(5)   # in percent
df_itd_maxClone.to_csv("flt3_itds_mostFrequentClonePerLength.csv")

########################################
# PRINT FILENAMES OF EACH CATEGORY TO FILE
			
out = open("flt3_insertions.txt","w")
out.write("\n".join(np.array(all_files)[w_ins["idx"]]) + "\n")
out.close()

out = open("flt3_insertions_single_itd_exact.txt","w")
out.write("\n".join(np.array(all_files)[w_itd_exact["idx"]]) + "\n")
out.close()

out = open("flt3_insertions_single_itd_nonexact.txt","w")
out.write("\n".join(np.array(all_files)[w_itd_nonexact["idx"]]) + "\n")
out.close()

out = open("flt3_insertions_single_itd_nonexact_fail.txt","w")
out.write("\n".join(np.array(all_files)[w_itd_nonexact_fail["idx"]]) + "\n")
out.close()


########################################
# PRINT SUMMARY STATISTICS on the number of reads in each category

print("--------------------")
print("Unique reads supporting each type of insert -> NOT number of distinct inserts!")
print("Insertions: {}".format(len(w_ins["idx"])))
print("Single exact ITD: {}".format(len(w_itd_exact["idx"])))
print("Single non-exact ITD: {}".format(len(w_itd_nonexact["idx"])))
print("Single insertion failed alignment: {}".format(len(w_itd_nonexact_fail["idx"])))


#########################################
# COLLECT DATA TO PLOT MOLM DILUTION SERIES' VAF AND READ COUNTS OF KNOWN ITD
# ----> LATER EXTEND THIS TO SUPPORT PLOTTING OF ANY/ALL KNOWN ITDs SUPPLIED?!

def get_known(klength, kstart, kstart2):
	kcounts = df_itd_grouped.ix[np.array(df_itd_grouped["length"] == klength) * np.array(df_itd_grouped["start"] == kstart) * np.array(df_itd_grouped["tandem2_start"] == kstart2)]["counts"]
	kvaf = (kcounts/sum(all_readCounts) * 100).round(5)   # in percent
	kindex = "_".join([str(x) for x in [klength,kstart,kstart2]])
#	
	known = pd.DataFrame({"itd": kindex, "counts": kcounts, "vaf": kvaf})
	print(known)
	known.to_csv("known_itds.csv")
	


# get known ITD counts for MOLM14 cellline
get_known(21,71,71) 
#get_known(21,77,77) 



