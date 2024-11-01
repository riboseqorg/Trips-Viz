from typing import Dict, List, Tuple, Union

import polars as pl

from tripsSplice import (get_protein_coding_transcript_ids,
                         get_reads_per_transcript_location,
                         get_start_stop_codon_positions)


def get_counts_meanFLD(transcripts: List[str],
                       read_file: str) -> Tuple[Dict[str, int], int]:
    """ 
    Get counts and mean FLD for a list of transcripts. 


    Parameters:
    transcripts (List[str]): list of transcript ids
    read_file (str): name of the sqlite file which contains information in dictionary.

    Returns:

    Example:
    >>> get_counts_meanFLD(["transcript1", "transcript2", "transcript3"],"read_file")
    ({'transcript1': 1, 'transcript2': 1, 'transcript3': 1}, if length not in length_fre

    read_file is the name of the sqlite file which contains information in dictionary. 

    read_file dictionary is as follows:
    {'transcript1': {'read_file1': 1, 'read_file2': 1, 'read_file3': 1}}, 
    
    """

    all_dfs = []
    for transcript in set(transcripts):
        #print "transcript, read_file", transcript, read_file
        reads = get_reads_per_transcript_location(transcript, read_file).with_columns(
            transctipt=pl.lit(transcript)
        )

        if not reads:
            continue
        all_dfs.append(reads)
    if all_dfs:

        reads = pl.concat(all_dfs)
        transcript_counts = reads.groupby("transctipt").agg(
            pl.sum("count").alias("count")
            
        )
        mean_fld = reads.groupby("length").agg( # Need to correct it
            pl.sum("count").alias("count")
    ).with_columns(
        pl.col("length") * pl.col("count") / pl.col("count").sum()
    )
        return transcript_counts, mean_fld


    else:
        return None, 1


def transcript_reads_per_kilobase(transcript_counts: Dict[str, int],
                                  cds_lengths: Dict[str, int],
                                  meanFLD: float) -> Dict[str, float]:
    """ 
    Calculate transcript per kilobase for a given read type.

    Parameters:
    - transcript_counts (Dict[str, int]): transcript counts
    - cds_lengths (Dict[str, int]): cds lengths
    - meanFLD (float): mean FLD

    Returns:

    Example:

    """
    cds_lengths = pl.DataFrame({'transcript': list(cds_lengths.keys()),'length': list(cds_lengths.values())}).with_columns(
        length = /(pl.col("length") -  meanFLD + 1)/1000. # Effective length per kb
    )
    transcript_counts = pl.DataFrame({'transcript': list(transcript_counts.keys()),'count': list(transcript_counts.values())})

    transcript_counts = cds_lengths.join(transcript_counts, on='transcript')

    return transcript_counts.with_columns(
        RPK = pl.col("count") / pl.col("length")
    )

def TPM(gene: str, sqlite_path_organism: str, sqlite_path_reads: List[str],
        seq_type: str) -> Dict[str, float]:
    """
    Calculate transcript per million for a given read type

    Parameters:
    - gene (str): The name of the gene
    - sqlite_path_organism (str): The path to the sqlite database file
    - sqlite_path_reads (List[str]): The path to the sqlite database file

    Returns:

    Example:

    """
    transcripts = []
    lengths = {}
    gene_lengths = get_protein_coding_transcript_ids(gene, sqlite_path_organism)
    if seq_type == "ribo":
        lengths = gene_lengths.filter(pl.col("tran_type") == 1).with_columns(lengths = pl.col("cds_stop") - pl.col("cds_start"))['lengths']

        

    else:  #if seq_type == "rna":
        lengths = gene_lengths["length"]

    all_TPMs = {transcript: [] for transcript in transcripts}

    for read_file in sqlite_path_reads:
        counts, meanFLD = get_counts_meanFLD(transcripts, read_file)

        if not counts:
            continue

        RPK = transcript_reads_per_kilobase(counts, lengths, meanFLD)
        per_million_scaling_factor = sum(RPK.values()) / 1000000.
        TPM = {}
        for transcript in transcripts:
            try:
                TPM[transcript] = round(
                    RPK[transcript] / per_million_scaling_factor, 2)
            except ZeroDivisionError:
                TPM[transcript] = 0
        for transcript in TPM:
            all_TPMs[transcript].append(TPM[transcript])

    avg_TPM = {
        transcript: sum(all_TPMs[transcript]) / len(all_TPMs[transcript])
        for transcript in all_TPMs
    }

    return avg_TPM
