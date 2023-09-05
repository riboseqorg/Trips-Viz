from typing import Dict, List
from tripsSplice import get_reads_per_genomic_location

from tripsSplice import genomic_exon_coordinate_ranges

from tripsSplice import get_protein_coding_transcript_ids
from tripsSplice import genomic_orf_coordinate_ranges


def get_unique_regions(
        genomic_exon_coordinates: Dict[str,
                                       List[int]]) -> Dict[str, List[int]]:

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


def count_readranges_supporting_exons_per_transcript(
        regions: Dict[str, List[int]],
        genomic_read_ranges: List[List[int]]) -> Dict[str, List[int]]:
    # Count the number of reads that overlap with each exon
    exons_counts = {}
    for read in genomic_read_ranges:
        for transcript in regions:
            exons_counts[transcript] = [0] * len(regions[transcript])

            for exon_num, exon in enumerate(regions[transcript]):

                exon_range = range(exon[0], exon[1] + 1)

                if (read[0] in exon_range) or (read[1] in exon_range):
                    exons_counts[transcript][exon_num] += genomic_read_ranges[
                        read]
    return exons_counts


def count_read_supporting_regions_per_transcript(
        regions: Dict[str, List[int]],
        genomic_read_positions: List[int]) -> Dict[str, List[int]]:
    exons_counts = {}
    for read in genomic_read_positions:
        for transcript in regions:
            if transcript not in exons_counts:
                exons_counts[transcript] = [
                    0 for i in range(len(regions[transcript]))
                ]

            exon_num = 0
            for exon in regions[transcript]:
                if exon[0] == exon[1]:
                    exon_num += 1
                    continue
                if (read in range(exon[0], exon[1] + 1)):
                    exons_counts[transcript][exon_num] += 1
                exon_num += 1
    return exons_counts


def get_coverage_per_region(
        regions: Dict[str, List[int]],
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

        for region, i in enumerate(counts[transcript]):
            if region_lengths[transcript][region] == 0:
                region_coverage[transcript].append(None)
            else:
                region_coverage[transcript].append(
                    round(
                        float(counts[transcript][region]) /
                        float(region_lengths[transcript][region]), 4))
    return region_coverage


def average_coverage_per_transcript(
        region_coverage: Dict[str, List[float]]) -> Dict[str, float]:
    average_coverage = {}
    for transcript in region_coverage:
        sum = 0
        count = 0
        if (not region_coverage[transcript]) and all(
                region_coverage[transcript]):
            print("transcript {transcript} had no unique regions".format(
                transcript=transcript))
            if transcript not in average_coverage:
                average_coverage[transcript] = None
        else:
            for region in region_coverage[transcript]:
                if region != None:
                    count += 1
                    sum += region

        if transcript not in average_coverage:
            average_coverage[transcript] = sum / count

    return average_coverage


def ribo_seq_read_counting(gene: str,
                           sqlite_path_organism: str,
                           sqlite_path_reads: List[str],
                           count_type: str = "range",
                           unique: bool = True) -> Dict[str, List[int]]:
    supported = get_protein_coding_transcript_ids(gene, sqlite_path_organism)
    exons = genomic_exon_coordinate_ranges(gene, sqlite_path_organism,
                                           supported)

    exclude = True
    orf_regions = genomic_orf_coordinate_ranges(gene, sqlite_path_organism,
                                                supported, exons)
    if unique:
        orf_regions = get_unique_regions(orf_regions)

    try:
        ["range", "fiveprime", "asite"].index(count_type)
    except ValueError:
        print(
            "The count type must be one of 'range', 'fiveprime' or 'asite'. "
            "count_type refers to the part of the read that is used in the feature counting process"
        )
        return "ERROR"

    if count_type == "range":
        genomic_read_ranges = get_reads_per_genomic_location(
            gene,
            sqlite_path_reads,
            sqlite_path_organism,
            supported,
            exons,
            filte_r=exclude,
            site="range")
        counts = count_readranges_supporting_exons_per_transcript(
            orf_regions, genomic_read_ranges)

    if count_type == "fiveprime":
        genomic_read_positions = get_reads_per_genomic_location(
            gene,
            sqlite_path_reads,
            sqlite_path_organism,
            supported,
            exons,
            filte_r=exclude,
            site="5prime")
        counts = count_read_supporting_regions_per_transcript(
            orf_regions, genomic_read_positions)

    if count_type == "asite":
        genomic_read_positions = get_reads_per_genomic_location(
            gene,
            sqlite_path_reads,
            sqlite_path_organism,
            supported,
            exons,
            filte_r=exclude,
            site="asite")
        counts = count_read_supporting_regions_per_transcript(
            orf_regions, genomic_read_positions)

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
    print rankings
    # transcripts = rna_seq_read_counting(gene, sqlite_path_organism, sqlite_path_reads, exclude=False, count_type="range")
    # print "To explain all reads in the selected files must include: ", transcripts
