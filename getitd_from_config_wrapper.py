import getitd
import easygui
import os


if __name__ == "__main__":
    # get analysis config
    config = getitd.load_config("config.txt")
    for key in config:
        getitd.config[key] = config[key]

    # convert strings to bools
    # (otherwise `"False"` is evaluated as `True`)
    getitd.config["REQUIRE_INDEL_FREE_PRIMERS"] = getitd.str_to_bool(
        getitd.config["REQUIRE_INDEL_FREE_PRIMERS"])
    getitd.config["INFER_SENSE_FROM_ALIGNMENT"] = getitd.str_to_bool(
        getitd.config["INFER_SENSE_FROM_ALIGNMENT"])


    # get FASTQ files to analyze
    fastq_files = []
    for file_ in os.listdir("."):
        if file_.endswith(".fastq"):
            fastq_files.append(file_)
    fastq_files.sort()
    print(fastq_files)

    # from FASTQ file names, extract sample names
    samples = []
    for file_ in fastq_files:
        sample = file_.split("R")[0][:-1]
        if sample not in samples:
            samples.append(sample)
    samples.sort()
    print(samples)
            
    # match samples and files to start getITD
    for sample in samples:
        for file_ in fastq_files:
            print(file_)
            if file_.startswith(sample + "_R1"):
                getitd.config["R1"] = file_
                continue
            if file_.startswith(sample + "_R2"):
                getitd.config["R2"] = file_
                break
        getitd.config["SAMPLE"] = sample

        getitd.main(getitd.config)
        

