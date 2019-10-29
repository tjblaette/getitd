import getitd
import easygui


if __name__ == "__main__":
    msg = "Enter the parameters to save to 'config.txt' for your next getITD analysis.\nWhen replacing default values, be sure to replace decimal values with decimals and whole numbers by whole numbers, respectively."
    title = "getITD parameters"

    ref_file = easygui.fileopenbox(msg="Select the text (.txt) file containing the WT sequence to use as a reference.", title="getITD WT reference sequence", default='*', filetypes=["*.txt"], multiple=False)
    anno_file = easygui.fileopenbox(msg="Select the tab-separated (.tsv) file containing the genomic, transcriptomic and proteomic annotation for each WT base of the previously selected reference.", title="getITD reference sequenceannotation", default='*', filetypes=["*.tsv"], multiple=False)

    fieldNames = [
#            "R1",
#            "R2",
#            "SAMPLE",
            "REF_FILE",
            "ANNO_FILE",
            "NKERN",
            "COST_GAPEXTEND",
            "COST_GAPOPEN",
            "COST_MATCH",
            "COST_MISMATCH",
            "MIN_BQS",
            "MIN_READ_COPIES",
            "MIN_READ_LENGTH",
            "MIN_SCORE_ALIGNMENTS",
            "MIN_SCORE_INSERTS",
            "MIN_TOTAL_READS",
            "MIN_UNIQUE_READS",
            "MIN_VAF",
            "TECH",
            "INFER_SENSE_FROM_ALIGNMENT",
            "FORWARD_PRIMERS",
            "REVERSE_PRIMERS",
            "REQUIRE_INDEL_FREE_PRIMERS",
            "MAX_TRAILING_BP",
            "FORWARD ADAPTER",
            "REVERSE ADAPTER"
            ]
    fieldValues = [
#            "",
#            "",
#            "",
            ref_file,
            anno_file,
            "4",
            "-0.5",
            "-20",
            "5",
            "-10",
            "30",
            "2",
            "100",
            "0.5",
            "0.5",
            "1",
            "2",
            "0.006",
            "Illumina",
            "False",
            "GCAATTTAGGTATGAAAGCCAGCTAC",
            "CTTTCAGCATTTTGACGGCAACC",
            "True",
            "0",
            "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGA",
            "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGA"
            ]
    fieldValues = easygui.multenterbox(msg, title, fieldNames, fieldValues)

    # make sure that none of the fields was left blank
    valid_input = False
    while True:
        errmsg = ""
        for key, val in zip(fieldNames, fieldValues):
          if val.strip() == "":
            errmsg = errmsg + '{} is a required field.\n\n'.format(key)
        if not errmsg:
            valid_input = True
            break
        fieldValues = easygui.multenterbox(errmsg, title, fieldNames, fieldValues)

    config = {}
    if fieldValues:
        for key, val in zip(fieldNames, fieldValues):
            config[key] = val
    print(config)
    
    if valid_input:
        getitd.save_config(config, "config.txt")
        
