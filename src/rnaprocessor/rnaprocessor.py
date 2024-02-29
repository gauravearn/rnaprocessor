import pandas as pd
import polars as pl
import seaborn as sns


def preparemiRNAsML(transcripts_file,\
                        psRNATargetfile,\
                             expectation_value,
                                     upstream, 
                                          downstream):
    """summary_line
    a machine learning preparation for the microRNAs
    after the target predictions. It will extract the 
    microRNA binding sites and also the upstream and the
    dowstream to prepare for the motif machine learning. 
    Keyword arguments:
    argument -- description
    transcripts_file: microRNA transcripts_file used in
    target predictions. 
    psRNATargetfile: psRNAtarget outcome. 
    expectation_value: expectation_value filter in the range of
    0.0 to 1.0
    upstream: how many bases upstream is needed for neural networks
    downstream: how many bases downstream is needed for neural networks
    Return: return_description
    """
    transcripts = list(filter(None,[x.strip() for x in open(transcripts_file).readlines()]))
    transcripts_read = {}
    for i in transcripts:   
        if i.startswith(">"):
            transcript_path = i.strip()
            if i not in transcripts_read:
                transcripts_read[i] = ""
                continue
        transcripts_read[transcript_path] += i.strip()
    with open(psRNATargetfile, "r") as psRNA:
        with open(psRNATargetfile + "processed.txt", "w") as psRNAprocessed:
            for line in psRNA.readlines():
                if line.startswith("#"):
                    continue
                psRNAprocessed.write(line)
    upstream_read = int(upstream)
    downstream_read = int(downstream)
    expectation = int(expectation_value)                                     
    fasta_keys = list(map(lambda n: n.split()[0].replace(">", ""),list(transcripts_read.keys())))
    fasta_sequences = list(transcripts_read.values())
    fasta_keys_sequences = [(i,j) for i,j in zip(fasta_keys,fasta_sequences)]
    miRNAs = pd.read_csv(psRNATargetfile + "processed.txt", sep= "\t")
    read_miRNAs = miRNAs[miRNAs["Expectation"] > expectation]
    read_miRNAs_target = list(read_miRNAs["Target_Acc."])
    read_miRNAs_accession = list(read_miRNAs["miRNA_Acc."])
    read_miRNAs_start = list(read_miRNAs["Target_start"])
    read_miRNAs_end = list(read_miRNAs["Target_end"])
    storing_miRNA_start_stop_coordinates = [(i,j-upstream_read,k+downstream_read) for i, j, \
                                          k in zip(read_miRNAs_target, read_miRNAs_start, read_miRNAs_end)]
    storing_matches = [(k,v) for k,v in fasta_keys_sequences for i in read_miRNAs_target if k == i]
    spliced_transcripts = [(storing_matches[j][0],storing_matches[j][1][storing_miRNA_start_stop_coordinates[i][1] \
                                      : storing_miRNA_start_stop_coordinates[i][2]]) for \
                                        i in range(len(storing_miRNA_start_stop_coordinates)) for \
                                                 j in range(len(storing_matches)) if intermediate[j][0] ==  \
                                                                                 storing_miRNA_start_stop_coordinates[i][0]]
    return spliced_transcriptsdef preparemiRNAsML(transcripts_file,\
                        psRNATargetfile,\
                             expectation_value,
                                     upstream, 
                                          downstream):
    """summary_line
    a machine learning preparation for the microRNAs
    after the target predictions. It will extract the 
    microRNA binding sites and also the upstream and the
    dowstream to prepare for the motif machine learning. 
    Keyword arguments:
    argument -- description
    transcripts_file: microRNA transcripts_file used in
    target predictions. 
    psRNATargetfile: psRNAtarget outcome. 
    expectation_value: expectation_value filter in the range of
    0.0 to 1.0
    upstream: how many bases upstream is needed for neural networks
    downstream: how many bases downstream is needed for neural networks
    Return: return_description
    """
    transcripts = list(filter(None,[x.strip() for x in open(transcripts_file).readlines()]))
    transcripts_read = {}
    for i in transcripts:   
        if i.startswith(">"):
            transcript_path = i.strip()
            if i not in transcripts_read:
                transcripts_read[i] = ""
                continue
        transcripts_read[transcript_path] += i.strip()
    with open(psRNATargetfile, "r") as psRNA:
        with open(psRNATargetfile + "processed.txt", "w") as psRNAprocessed:
            for line in psRNA.readlines():
                if line.startswith("#"):
                    continue
                psRNAprocessed.write(line)
    upstream_read = int(upstream)
    downstream_read = int(downstream)
    expectation = int(expectation_value)                                     
    fasta_keys = list(map(lambda n: n.split()[0].replace(">", ""),list(transcripts_read.keys())))
    fasta_sequences = list(transcripts_read.values())
    fasta_keys_sequences = [(i,j) for i,j in zip(fasta_keys,fasta_sequences)]
    miRNAs = pd.read_csv(psRNATargetfile + "processed.txt", sep= "\t")
    read_miRNAs = miRNAs[miRNAs["Expectation"] > expectation]
    read_miRNAs_target = list(read_miRNAs["Target_Acc."])
    read_miRNAs_accession = list(read_miRNAs["miRNA_Acc."])
    read_miRNAs_start = list(read_miRNAs["Target_start"])
    read_miRNAs_end = list(read_miRNAs["Target_end"])
    storing_miRNA_start_stop_coordinates = [(i,j-upstream_read,k+downstream_read) for i, j, \
                                          k in zip(read_miRNAs_target, read_miRNAs_start, read_miRNAs_end)]
    storing_matches = [(k,v) for k,v in fasta_keys_sequences for i in read_miRNAs_target if k == i]
    spliced_transcripts = [(storing_matches[j][0],storing_matches[j][1][storing_miRNA_start_stop_coordinates[i][1] \
                                      : storing_miRNA_start_stop_coordinates[i][2]]) for \
                                        i in range(len(storing_miRNA_start_stop_coordinates)) for \
                                                 j in range(len(storing_matches)) if intermediate[j][0] ==  \
                                                                                 storing_miRNA_start_stop_coordinates[i][0]]
    return spliced_transcripts
