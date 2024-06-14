from typing import Dict, List, Literal, Tuple


def count_read_supporting_regions_per_transcript(
        regions: Dict[str, List[Tuple[int, int]]],
        genomic_read_positions: List[int]) -> Dict[str, List[int]]:
    """

    Parameters:
    - regions (Dict[str, List[Tuple[int, int]]]): dictionary
    with transcript ids as keys and list of coordinates as values
    - genomic_read_positions (List[int]): list of read positions

    Returns:

    Example:
    """
    exons_counts = {}
    for read in genomic_read_positions:
        for transcript in regions:
            if transcript not in exons_counts:
                exons_counts[transcript] = [0] * len(regions[transcript])

            for exon_num, exon in enumerate(regions[transcript], 1):
                if (exon[0] != exon[1]) and (exon[0] <= read <= exon[1]):
                    exons_counts[transcript][exon_num] += 1
    return exons_counts
