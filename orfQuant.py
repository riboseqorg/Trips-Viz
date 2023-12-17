#!/usr/bin/env python
from typing import Dict, List
import pandas as pd
from pandas.core.frame import DataFrame
from numpy.typing import NDArray

from sqlitedict import SqliteDict

from tripsCount import count_read_supporting_regions_per_transcript
from tripsSplice import (
    genomic_exon_coordinate_ranges,
    genomic_junction_positions,
    genomic_junction_scores,
    get_protein_coding_transcript_ids,
    get_reads_per_genomic_location,
    get_start_stop_codon_positions,
)


def classify_regions_shared_unique(
        regions: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """

    Parameters:
    - regions: dictionary of genomic coordinates that make up the features for each transcript

    Returns:

    Example:


    """
    # Takes in dictionary of genomic coordinates that make up the features for each transcript
    # for each region consider each transcript and count the number of occurences of this exact region in that transcript
    # If this count is > 1 it is considered not unique. Otherwise unique

    classified = {}
    for tran1 in regions:
        if tran1 not in classified:
            classified[tran1] = {"unique": [], "shared": []}
        for region in regions[tran1]:
            counter = 0
            for tran2 in regions:
                counter += regions[tran2].count(region)

            if counter > 1:
                classified[tran1]["shared"].append(region)
            else:
                classified[tran1]["unique"].append(region)
    return classified


def region_coverage(
        regions: Dict[str, List[str]],
        counts: Dict[str, List[int]]) -> Dict[str, Dict[str, float]]:
    # Returns a dictioanry of coverage of each coding region. Calculated as counts/length rounded to 4 places
    region_counts = {}
    for transcript in regions:
        region_counts_transcript = zip(regions[transcript], counts[transcript])
        if transcript not in region_counts:
            region_counts[transcript] = region_counts_transcript

    coverage_region_transcript = {}
    for transcript in region_counts:
        if transcript not in coverage_region_transcript:
            coverage_region_transcript[transcript] = {}

        for item in region_counts[transcript]:
            length = item[0][1] - item[0][0]
            coverage = round(float(item[1]) / float(length), 4)
            coverage_region_transcript[transcript][item[0]] = coverage

    return coverage_region_transcript


def coverage_junction_transcript(
    junction_scores: Dict[str, Dict[str,
                                    float]]) -> Dict[str, Dict[str, float]]:
    # function to calculate the normalised coverage of junction features. returns a nested dict.
    # 60 = 30 + 30 which is the sum of two riboseq read lengths.
    coverage_junction = {}
    for transcript in junction_scores:
        if transcript not in coverage_junction:
            coverage_junction[transcript] = {}
        for junction in junction_scores[transcript]:
            coverage = float(junction_scores[transcript][junction]) / float(60)
            coverage_junction[transcript][junction] = coverage

    return coverage_junction


def average_unique_coverage(
        unique_shared_exons: Dict[str, Dict[str, Dict[str, int]]],
        unique_shared_junctions: Dict[str, Dict[str, Dict[str, int]]],
        coverage_junction: Dict[str, Dict[str, float]],
        coverage_exons: Dict[str, Dict[str, float]]) -> Dict[str, float]:
    # Return the average coverage of all unique features. Sum of coverages over number of unique features

    average_features = {}

    for transcript in unique_shared_exons:
        sum = 0
        count = 0
        for item in unique_shared_exons[transcript]['unique']:
            sum += coverage_exons[transcript][item]
            count += 1

        for junction in unique_shared_junctions[transcript]['unique']:
            sum += coverage_junction[transcript][junction]
            count += 1

        if transcript not in average_features:
            if unique_shared_exons[transcript][
                    'unique'] == [] and unique_shared_junctions[transcript][
                        'unique'] == []:
                average_features[transcript] = 0
            else:
                average_features[transcript] = float(sum) / float(count)

    return average_features


def all_feature_average(
        coverage_junction: Dict[str, Dict[str, float]],
        coverage_exons: Dict[str, Dict[str, float]]) -> Dict[str, float]:
    # Average of the coverage of all features. Sum of coverages over number of features
    average_features = {}

    for transcript in coverage_exons:
        sum = 0
        count = 0
        for exon in coverage_exons[transcript]:
            sum += coverage_exons[transcript][exon]
            count += 1

        for junction in coverage_junction[transcript]:
            sum += coverage_junction[transcript][junction]
            count += 1

        if transcript not in average_features:
            average_features[transcript] = float(sum) / float(count)
    return average_features


def cORF_ratio(average_unique: Dict[str, float],
               average_all: Dict[str, float]) -> Dict[str, float]:
    """Calculate the ratio of average unique coverage to average all coverage."""
    cORF = {}
    for transcript in average_all:
        cORF[transcript] = average_unique[transcript] / average_all[transcript]

    return cORF


def adjusted_coverage_for_non_unique_orfs(
        transcript: str, coverage_exons: Dict[str, Dict[str, float]],
        coverage_junctions: Dict[str, Dict[str, float]],
        average_coverage_all: Dict[str, float]) -> float:
    """
    This is coverage of exons - coverage of exons * sum of all coverages from that transcript
    Calculate an adjusted value for cases where a transcript has no unique features
    Average is calculated of these adjusted values
    """

    adjusted_coverage_exons = {transcript: {}}
    for exon in coverage_exons[transcript]:
        signal = 0
        for trans in coverage_exons:
            if exon in coverage_exons[trans]:
                signal += coverage_exons[trans][exon]

        if exon not in adjusted_coverage_exons[transcript]:
            adjusted_coverage_exons[transcript][exon] = coverage_exons[
                transcript][exon] - (coverage_exons[transcript][exon] * signal)

    adjusted_coverage_junctions = {transcript: {}}
    for junction in coverage_junctions[transcript]:
        signal = 0
        for trans in coverage_junctions:
            if junction in coverage_junctions[trans]:
                signal += coverage_junctions[trans][junction]

        if junction not in adjusted_coverage_junctions[transcript]:
            adjusted_coverage_junctions[transcript][
                junction] = coverage_junctions[transcript][junction] - (
                    coverage_junctions[transcript][junction] * signal)
    sum = 0
    count = 0
    for region in adjusted_coverage_exons[transcript]:
        sum += adjusted_coverage_exons[transcript][region]
        count += 1

    for junction in adjusted_coverage_junctions[transcript]:
        sum += adjusted_coverage_junctions[transcript][junction]
        count += 1
    average_adjusted_coverage = float(sum) / float(count)

    cORF_transcript = average_adjusted_coverage / average_coverage_all[
        transcript]
    return cORF_transcript


def shared_coverage(
    coverage_exons: Dict[str, Dict[str, float]],
    coverage_junctions: Dict[str, Dict[str, float]]
) -> Dict[str, Dict[str, float]]:
    # Calculate an adjusted average value where all features are shared between transcripts.

    shared_exon_coverage = {}
    for transcript in coverage_exons:
        if transcript not in shared_exon_coverage:
            shared_exon_coverage[transcript] = {}
        for exon in coverage_exons[transcript]:
            number_ORF_over_F = 0
            for tran in coverage_exons:
                for exon2 in coverage_exons[tran]:
                    if transcript == tran:
                        continue

                    if exon == exon2:
                        number_ORF_over_F += 1
            number_ORF_over_F += 1

            shared_exon_coverage[transcript][
                exon] = coverage_exons[transcript][exon] / number_ORF_over_F

    shared_junction_coverage = {}
    for transcript in coverage_junctions:
        if transcript not in shared_junction_coverage:
            shared_junction_coverage[transcript] = {}
        for junction in coverage_junctions[transcript]:
            number_ORF_over_F = 0
            for tran in coverage_junctions:
                for junction2 in coverage_junctions[tran]:
                    if transcript == tran:
                        continue

                    if junction == junction2:
                        number_ORF_over_F += 1
            number_ORF_over_F += 1
            shared_junction_coverage[transcript][
                junction] = coverage_junctions[transcript][
                    junction] / number_ORF_over_F

    shared_cORF = {}
    for transcript in shared_exon_coverage:
        sum = 0
        count = 0
        for exon in shared_exon_coverage[transcript]:
            sum += shared_exon_coverage[transcript][exon]
            count += 1

        for junction in shared_junction_coverage[transcript]:
            sum += shared_junction_coverage[transcript][junction]
            count += 1

        if transcript not in shared_cORF:
            shared_cORF[transcript] = float(sum) / float(count)

    return shared_cORF


def aORF(cORF: DataFrame, counts: DataFrame) -> DataFrame:
    # Returns a dictionary of the number of a sites found in each transcript. Coverage by the sum of counts.
    cORF = cORF.drop_duplicates('transcript')
    counts = counts.drop_duplicates('transcript')
    cORF_count = cORF.merge(counts, on='transcript')
    cORF_count['a_sites_perORF'] = cORF_count.apply(
        lambda x: sum(x['counts']) * x['cORF'], axis=1)
    return cORF_count


def lORF(coding: List[str], sqlite_path_organism: str) -> DataFrame:
    # Returns the lengths of each open reading frame in a dictionary. Transcript as keys, lengths as values

    lORFs = get_start_stop_codon_positions(
        coding, sqlite_path_organism).drop_duplicates('transcript')
    lORFs['length'] = lORFs['stop'] - lORFs['start'] + 1
    return lORFs[['transcript', 'length']]


def ORFs_per_million(aORF: Dict[str, float],
                     lORF: Dict[str, int]) -> Dict[str, float]:
    # Return the TPM like value OPM. Normalise coverage by length of each orf.
    # One million divided by scaling factor is the sum of all normalised coverages

    aORF_lORF = {}
    full_set = 0
    for transcript in aORF:
        if transcript not in aORF_lORF:
            aORF_lORF[transcript] = aORF[transcript] / lORF[transcript]
            full_set += aORF[transcript] / lORF[transcript]

    orfs_per_million = {}
    for transcript in aORF_lORF:
        if transcript not in orfs_per_million:
            orfs_per_million[transcript] = aORF_lORF[transcript] * (10**6 /
                                                                    full_set)
    return orfs_per_million


def orfQuant(sqlite_path_organism: str, sqlite_path_reads: str,
             supported: List[str], counts: Dict[str, List],
             exons: Dict[str, Dict[str, float]]) -> Dict[str, float]:
    # Main function for implementing method. Executes above functions and determines if all features are shared
    # Returms OPM values from ORFs_per_million function
    junctions = genomic_junction_positions(exons)
    junction_scores = genomic_junction_scores(sqlite_path_organism,
                                              sqlite_path_reads, supported,
                                              junctions)

    unique_shared_exons = classify_regions_shared_unique(exons)
    coverage_exons = region_coverage(exons, counts)

    unique_shared_junctions = classify_regions_shared_unique(junctions)
    coverage_junctions = coverage_junction_transcript(junction_scores)

    average_unique = average_unique_coverage(unique_shared_exons,
                                             unique_shared_junctions,
                                             coverage_junctions,
                                             coverage_exons)
    average_all = all_feature_average(coverage_junctions, coverage_exons)

    cORF = cORF_ratio(average_unique, average_all)

    all_shared = True
    for transcript in supported:
        if unique_shared_junctions[transcript][
                "unique"] == [] and unique_shared_exons[transcript][
                    "unique"] == []:
            cORF[transcript] = adjusted_coverage_for_non_unique_orfs(
                transcript, coverage_exons, coverage_junctions, average_all)
        else:
            all_shared = False

    if all_shared:
        cORF = shared_coverage(coverage_exons, coverage_junctions)
    else:
        cORF = cORF_ratio(average_unique, average_all)
    # -- Need to0 recheck
    cORF = pd.DataFrame.from_dict(cORF, orient="index").rename(columns={
        0: "cORF"
    }).drop_duplicates("cORF")
    # counts = pd.DataFrame.from_dict(counts, orient="index").rename(columns={
    # 0: "counts"
    # }).drop_duplicates("counts")

    adjusted_a_sites = aORF(cORF, counts)
    orf_lengths = lORF(supported, sqlite_path_organism)

    orfs_per_million = ORFs_per_million(adjusted_a_sites, orf_lengths)

    return orfs_per_million


def incl_OPM_run_orfQuant(gene: str, sqlite_path_organism: str,
                          sqlite_path_reads: NDArray) -> Dict[str, float]:
    #WARNING: Find what is OPM

    # Run the ORFquant algorithm on the files that do not have a OPM value stored for the required trancripts
    # Store calculated OPMs for a file back in the read sqlitedict under "OPM"
    # Returns the average OPMs of the selected files for this locus
    # print "gene", gene
    # print "sqlitepath", sqlite_path_organism
    # print "sqlite_path_reads",sqlite_path_reads
    coding = get_protein_coding_transcript_ids(gene, sqlite_path_organism)
    exons = genomic_exon_coordinate_ranges(gene, sqlite_path_organism, True)

    transcript_OPMs = {}
    read_files = []
    for file in sqlite_path_reads:
        infile = SqliteDict(file)
        for transcript in coding:
            try:
                if transcript not in transcript_OPMs:
                    transcript_OPMs[transcript] = []
                transcript_OPMs[transcript].append(infile[transcript]["OPM"])
            except KeyError:
                read_files.append(file)
        infile.close()

    for file in set(read_files):
        genomic_read_positions = get_reads_per_genomic_location(
            gene, [file],
            sqlite_path_organism,
            coding,
            exons,
            filte_r=True,
            site="asite")

        counts = count_read_supporting_regions_per_transcript(
            exons, genomic_read_positions)

        orfQuant_res = orfQuant(sqlite_path_organism, [file], coding, counts,
                                exons)
        infile = SqliteDict(file, autocommit=False)
        for transcript in orfQuant_res:
            transcript_dict = infile[transcript].copy()
            transcript_dict['OPM'] = orfQuant_res[transcript]
            infile[transcript] = transcript_dict

            if transcript not in transcript_OPMs:
                transcript_OPMs[transcript] = [orfQuant_res[transcript]]
            else:
                transcript_OPMs[transcript].append(orfQuant_res[transcript])
        infile.commit()
        infile.close()

    avg_opms_per_transcript = {}
    for transcript in transcript_OPMs:
        avg_opms_per_transcript[transcript] = round(
            sum(transcript_OPMs[transcript]) /
            float(len(transcript_OPMs[transcript])), 2)
    return avg_opms_per_transcript


# if __name__ == "__main__":
#     gene = "phpt1"
#     # sqlite_path_organism = "/home/jackt/Tools/Trips-Viz/Trips-Viz-master/trips_annotations/homo_sapiens/old_homo_sapiens.Gencode_v25.sqlite"
#     sqlite_path_organism = "/home/jackt/Tools/Trips-Viz/Trips-Viz-master/trips_annotations/homo_sapiens/homo_sapiens.Gencode_v25.sqlite"

#     # sqlite_path_reads = ["/home/jackt/Tools/Trips-Viz/Trips-Viz-master/trips_shelves/rnaseq/homo_sapiens/Park16/SRR3306574.sqlite"]
#     sqlite_path_reads = [ '/home/jackt/Tools/Trips-Viz/Trips-Viz-master/trips_shelves/riboseq/homo_sapiens/Park16/SRR3306586.sqlite', '/home/jackt/Tools/Trips-Viz/Trips-Viz-master/trips_shelves/riboseq/homo_sapiens/Park16/SRR3306587.sqlite']
#     orfQuant_res = incl_OPM_run_orfQuant(gene, sqlite_path_organism, sqlite_path_reads)

#     print orfQuant_res
