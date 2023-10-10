from typing import Dict, List, Tuple, Union
import numpy as np
from tripsSplice import get_reads_per_transcript_location
from tripsSplice import get_protein_coding_transcript_ids
from tripsSplice import get_start_stop_codon_positions


def get_counts_meanFLD(transcripts: List[str],
                       read_file: str) -> Tuple[Dict[str, int], int]:
    ''' 
    Get counts and mean FLD for a list of transcripts. 
    >>> get_counts_meanFLD(["transcript1", "transcript2", "transcript3"],"read_file")
    ({'transcript1': 1, 'transcript2': 1, 'transcript3': 1}, if length not in length_fre

    read_file is the name of the sqlite file which contains information in dictionary. 

    read_file dictionary is as follows:
    {'transcript1': {'read_file1': 1, 'read_file2': 1, 'read_file3': 1}}, 
    
    '''

    length_freq = {}
    transcript_counts = {}
    for transcript in set(transcripts):
        transcript_counts[transcript] = 0
        #print "transcript, read_file", transcript, read_file
        reads = get_reads_per_transcript_location(transcript, read_file)
        if not reads:
            continue
        for length, positions in reads.items():
            if length not in length_freq:
                length_freq[length] = 0
            for position in positions:
                length_freq[length] += reads[length][position]
                transcript_counts[transcript] += reads[length][position]
    lengths = np.array(list(length_freq.keys()))
    freqs = np.array(list(length_freq.values()))
    total_count = np.sum(lengths * freqs)
    count = np.sum(freqs)
    mean_fld = 1 if not count else round(total_count * 1. / count)
    return transcript_counts, mean_fld


def transcript_reads_per_kilobase(transcript_counts: Dict[str, int],
                                  cds_lengths: Dict[str, int],
                                  meanFLD: float) -> Dict[str, float]:
    RPK = {}
    for transcript in cds_lengths:
        effective_length = cds_lengths[transcript] - meanFLD + 1
        effective_length_per_kilobase = effective_length / 1000

        try:
            RPK[transcript] = transcript_counts[
                transcript] / effective_length_per_kilobase
        except ZeroDivisionError:
            RPK[transcript] = 0
    return RPK


def TPM(gene: str, sqlite_path_organism: str, sqlite_path_reads: List[str],
        seq_type: str) -> Dict[str, float]:
    """
    Calculate transcript per million for a given read type

    """
    transcripts = []
    lengths = {}
    if seq_type == "ribo":
        transcripts = get_protein_coding_transcript_ids(
            gene, sqlite_path_organism)
        start_stops = {
            transcript:
            get_start_stop_codon_positions(transcript, sqlite_path_organism)
            for transcript in transcripts
        }
        lengths = {
            transcript: start_stop[1] - start_stop[0]
            for transcript, start_stop in start_stops.items()
        }

    else:  #if seq_type == "rna":
        transcripts = [
            transcript[0]
            for transcript in get_gene_info(gene, sqlite_path_organism)
        ]
        lengths = get_transcript_length(
            gene, sqlite_path_organism)  # TODO: Look for this function

    all_TPMs = {transcript: [] for transcript in transcripts}

    for read_file in sqlite_path_reads:
        counts, meanFLD = get_counts_meanFLD(transcripts, read_file)

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
