from typing import Dict, List, Literal, Tuple
from tripsSplice import get_reads_per_genomic_location

from tripsSplice import genomic_exon_coordinate_ranges

from tripsSplice import get_protein_coding_transcript_ids
from tripsSplice import genomic_orf_coordinate_ranges
import pandas as pd


def get_unique_regions(
    genomic_exon_coordinates: Dict[str, List[Tuple[int, int]]]
) -> Dict[str, List[Tuple[int, int]]]:

    unique_regions = {}
    exon_coordinates = []

    for transcript in genomic_exon_coordinates:
        for exon in genomic_exon_coordinates[transcript]:
            exon_coordinates.append(exon)

    for transcript in genomic_exon_coordinates:
        if transcript not in unique_regions:
            unique_regions[transcript] = []
        for exon in genomic_exon_coordinates[transcript]:
            unique_exon = exon
            exact_match_counter = 0

            for coordinates in exon_coordinates:
                if exon == coordinates:
                    exact_match_counter += 1
                    if exact_match_counter > 1:
                        unique_exon = (exon[1], exon[1])
                    continue
                if (unique_exon[0] in range(
                        coordinates[0],
                        coordinates[1] + 1)) and (unique_exon[1] in range(
                            coordinates[0], coordinates[1] + 1)):
                    unique_exon = (unique_exon[1], unique_exon[1])
                    continue

                if (unique_exon[0] in range(coordinates[0],
                                            coordinates[1] + 1)):
                    unique_exon = (coordinates[1], unique_exon[1])

                if (unique_exon[1] in range(coordinates[0],
                                            coordinates[1] + 1)):
                    unique_exon = (unique_exon[0], coordinates[0])

                if (unique_exon == coordinates):
                    unique_exon = (exon[1], exon[1])

            unique_regions[transcript].append(unique_exon)
    return unique_regions


def count_read_supporting_regions_per_transcript(
        regions: Dict[str, List[Tuple[int, int]]],
        genomic_read_positions: List[int]) -> Dict[str, List[int]]:
    exons_counts = {}
    for read in genomic_read_positions:
        for transcript in regions:
            if transcript not in exons_counts:
                exons_counts[transcript] = [0] * len(regions[transcript])

            exon_num = 0
            for exon in regions[transcript]:
                exon_num += 1
                if exon[0] == exon[1]:
                    continue
                if (read in range(exon[0], exon[1] + 1)):
                    exons_counts[transcript][exon_num] += 1
    return exons_counts


def get_coverage_per_region(
        regions: Dict[str, List[Tuple[int, int]]],
        counts: Dict[str, List[int]]) -> Dict[str, List[float]]:

    region_lengths = {}
    for transcript in regions:
        if transcript not in region_lengths:
            region_lengths[transcript] = []

        for region in regions[transcript]:
            region_lengths[transcript].append(region[1] - region[0])

    #coverage = reads/length of region
    region_coverage = {}
    for transcript in counts:
        if transcript not in region_coverage:
            region_coverage[transcript] = []

        for region, _ in enumerate(counts[transcript]):
            if region_lengths[transcript][region] == 0:
                region_coverage[transcript].append(None)
            else:
                region_coverage[transcript].append(
                    round(
                        float(counts[transcript][region]) /
                        float(region_lengths[transcript][region]), 4))
    return region_coverage


def average_coverage_per_transcript(
        region_coverage: Dict[str, List[float]]) -> Dict[str, float | None]:

    # TODO: Check what values are being passed here
    average_coverage = {}
    for transcript in region_coverage:
        total = 0
        count = 0
        if (not region_coverage[transcript]) and all(
                region_coverage[transcript]):
            # print(f"transcript {transcript} had no unique regions")
            if transcript not in average_coverage:
                average_coverage[transcript] = None
        else:
            for region in region_coverage[transcript]:
                if region:
                    count += 1
                    total += region

        if transcript not in average_coverage:
            average_coverage[transcript] = total / count

    return average_coverage


def ribo_seq_read_counting(
        gene: str,
        sqlite_path_organism: str,
        sqlite_path_reads: List[str],
        count_type: Literal["range", "fiveprime", "asite"] = "range",
        unique: bool = True) -> Dict[str, float | None] | str:
    supported = get_protein_coding_transcript_ids(gene, sqlite_path_organism)
    exons = genomic_exon_coordinate_ranges(gene, sqlite_path_organism,
                                           supported)

    orf_regions = genomic_orf_coordinate_ranges(sqlite_path_organism,
                                                supported, exons)
    if unique:
        orf_regions = get_unique_regions(orf_regions)

    if count_type not in ["range", "fiveprime", "asite"]:
        return "ERROR"

    genomic_read_positions = get_reads_per_genomic_location(
        gene,
        sqlite_path_reads,
        sqlite_path_organism,
        supported,
        exons,
        filte_r=True,
        site=count_type)
    counts = count_read_supporting_regions_per_transcript(
        orf_regions, list(genomic_read_positions.keys()))

    region_coverage = get_coverage_per_region(orf_regions, counts)
    coverage = average_coverage_per_transcript(region_coverage)
    return coverage


if __name__ == "__main__":

    gene = "phpt1"
    sqlite_path_organism = "/home/jackt/Tools/Trips-Viz/Trips-Viz-master/trips_annotations/homo_sapiens/homo_sapiens.Gencode_v25.sqlite"
    sqlite_path_reads = [
        "/home/jackt/Tools/Trips-Viz/Trips-Viz-master/trips_shelves/rnaseq/homo_sapiens/Park16/SRR3306574.sqlite"
    ]
    rankings = ribo_seq_read_counting(gene,
                                      sqlite_path_organism,
                                      sqlite_path_reads,
                                      count_type="asite",
                                      unique=True)
    print(rankings)
    # transcripts = rna_seq_read_counting(gene, sqlite_path_organism, sqlite_path_reads, exclude=False, count_type="range")
    # print "To explain all reads in the selected files must include: ", transcripts
