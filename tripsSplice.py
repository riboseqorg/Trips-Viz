from typing import Dict, Tuple, List, Union
from typing_extensions import Literal
from pandas.core.frame import DataFrame
from sqlqueries import sqlquery
from sqlitedict import SqliteDict
import pandas as pd


class TripsSplice:

    def __init__(self, sqlite_path_organism: str) -> None:
        self.transcript_table = sqlquery(sqlite_path_organism, 'transcripts')

    def string_to_integer_list(self, lst: List[str]) -> List[int]:
        return [int(i) if i else 0 for i in lst]

    def get_gene_info(self, gene_id: str) -> DataFrame:
        return self.transcript_table.loc[
            self.transcript_table.gene == gene_id,
            ['transcript', 'exon_junctions', 'sequence']]

    def get_transcript_length(self, gene: str) -> DataFrame:
        return self.transcript_table.loc[self.transcript_table.gene == gene,
                                         ["transcript", "length"]]

    def get_genes_principal(self, gene_id: str) -> str:
        return self.transcript_table.loc[self.transcript_table.gene == gene_id
                                         & self.transcript_table.principle,
                                         "transcript"].values[0]

    def get_transcript_info(self, transcript_id: str) -> DataFrame:
        return self.transcript_table[self.transcript_table.transcript ==
                                     transcript_id,
                                     ['exon_junctions', 'sequence']].iloc[0]

    def _exon_coordinates_list(self, lst: List[int],
                               transcript_length: int) -> List[List[int]]:
        # Return the list of integers from a list of strings. If the list is empty then return 0
        if not lst:
            return [[0, transcript_length]]
        lst = [-1] + [int(i) for i in lst] + [transcript_length]
        return [[lst[i] + 1, lst[i + 1]] for i in range(len(lst) - 1)]

    def get_exon_coordinates_for_orf(self,
                                     transcript_id: str) -> List[List[int]]:
        transcript_info = self.get_transcript_info(transcript_id)
        exon_junctions = self.string_to_integer_list(
            transcript_info['exon_junctions'].split(","))

        exon_coordinates = self._exon_coordinates_list(
            exon_junctions, len(transcript_info['sequence']))
        return exon_coordinates

    def get_protein_coding_transcript_ids(self, gene: str) -> DataFrame:
        return self.transcript_table[self.transcript_table.gene == gene]

    def get_start_stop_codon_positions(self, transcript_id: str) -> DataFrame:
        return self.transcript_table[self.transcript_table.transcript ==
                                     transcript_id]

    def get_orf_exon_coordinates(self, transcript_id: str) -> DataFrame:
        return self.transcript_table[self.transcript_table.transcript ==
                                     transcript_id]


def get_3prime_exon(junction_list: List[int],
                    sequence: str) -> Tuple[str, str]:
    # Part of the process of producing the sequences of all exons in a transcript.
    # this function slices the 3' exon sequence from the transcript sequence returning
    # both the exon sequence and the remaining sequence of the transcript
    cut_off = junction_list[-1]
    exon = sequence[cut_off - 1:]
    seq_less_exon = sequence[:cut_off - 1]
    return exon, seq_less_exon


def get_exon_coordinates_for_orf(transcript_id: str,
                                 sqlite_path_organism: str) -> List[List[int]]:
    # return the coordinates of the exons in a given transcript. 0 -> first junction, first junction -> second junction
    # last junction -> end
    tripsplice = TripsSplice(sqlite_path_organism)
    transcript_info = tripsplice.get_transcript_info(transcript_id)
    exon_junctions = list(map(int, transcript_info[0][0].split(",")))
    transcript_end = len(transcript_info[0][1])
    exon_coordinates = []
    for junction in exon_junctions:
        if exon_junctions.index(junction) == 0:
            exon_coordinates.append((0, junction))
        else:
            exon_coordinates.append(
                (exon_junctions[exon_junctions.index(junction) - 1] + 1,
                 junction))

    exon_coordinates.append((exon_junctions[-1] + 1, transcript_end))

    return exon_coordinates


def get_protein_coding_transcript_ids(gene: str,
                                      sqlite_path_organism: str) -> List[str]:
    """gets the transcript IDs for protein coding genes for the given gene"""
    gene_id = gene.upper()
    transcripts = sqlquery(sqlite_path_organism, "transcripts")
    protein_coding_transcripts = transcripts.loc[transcripts.gene == gene_id
                                                 & transcripts.tran_type == 1,
                                                 "transcript"].values

    return protein_coding_transcripts


def get_start_stop_codon_positions(
        transcript_id: str, sqlite_path_organism: str) -> Tuple[int, int]:
    # Get start and stop codon positions from annotation sqlite
    transcripts = get_table("transcripts", sqlite_path_organism)
    return transcripts.loc[transcripts.transcript == transcript_id,
                           ['cds_start', 'cds_stop']].values[0]


def get_orf_exon_structure(
        start_stop: Tuple[int, int],
        exon_coordinates: List[List[int]]) -> List[List[int]]:
    # determine the stucture of each ORF. Returns coordinates in the form:
    # Initiation site to junction, junction to junction, junction to translation stop
    start, stop = start_stop
    # TODO: Sort exon coordinates and select last
    exon_coordinates_df = DataFrame(exon_coordinates,
                                    columns=["start", "stop"])
    exon_coordinates_df = exon_coordinates_df[
        exon_coordinates_df.stop > start
        & exon_coordinates_df.start < stop]
    exon_coordinates_df.loc[exon_coordinates_df.start <= start,
                            "start"] = start
    exon_coordinates_df.loc[exon_coordinates_df.stop >= stop, "stop"] = stop
    return exon_coordinates_df.values

    # started = False
    # structure = []
    # for exon in exon_coordinates:
    # if start in range(exon[0], exon[1]):
    # started = True
    # structure.append((start, exon[1]))

    # elif started and (stop in range(exon[0], exon[1])):
    # structure.append((exon[0], stop))
    # elif started:  # NOTE:What if another exon is after the stop codon? It will insert wrong exon
    # structure.append(exon)
    # return structure


def get_coding_regions_for_genes_transcripts(
        gene: str, sqlite_path_organism: str) -> List[Tuple[int, int]]:
    # Returns the coding region coordinates annotated on each transcript for a given gene
    gene_id = gene.upper()
    transcripts = get_table("transcripts", sqlite_path_organism)
    transcripts = transcripts.loc[transcripts.gene == gene_id,
                                  "transcript"].values
    coding_regions = get_table("coding_regions", sqlite_path_organism)
    coding_regions = coding_regions[coding_regions.transcript.isin(
        transcripts)]
    return coding_regions

    # conn = sqlite3.connect(sqlite_path_organism)
    # c = conn.cursor()
    # c.execute(
    # 'SELECT * FROM coding_regions WHERE transcript IN (SELECT transcript FROM transcripts WHERE gene = (:gene));',
    # {'gene': gene_id})
    # coding_regions = c.fetchall()
    # return coding_regions


def exons_of_transcript(transcript_id: str,
                        sqlite_path_organism: str) -> List[str]:
    # For a given transcript return the exon sequences in a list in 5' to 3' direction
    # WARNING: This function doesn't make sense
    exon_lst = []
    tripsplice = TripsSplice(sqlite_path_organism)
    trans_info = tripsplice.get_transcript_info(transcript_id)
    exon_junct_int = list(map(int, trans_info[0][0].split(",")))
    sequence = trans_info[0][1]

    while len(exon_junct_int) != 0:
        exon, sequence = get_3prime_exon(exon_junct_int, sequence)
        exon_lst.append(exon)

        exon_junct_int.pop(-1)
    exon_lst.append(sequence)
    return exon_lst[::-1]


def get_exon_coordinate_ranges(sequence: str, exons: List[str],
                               junctions: List[int]) -> List[Tuple[int, int]]:
    """Return list of transcript coordinate ranges (start position, stop position)."""
    end = len(sequence) - 1
    junctions.append(end)

    ranges = []
    for i, exon in enumerate(exons):
        ranges.append((sequence.find(exon), junctions[i] - 1))
    return ranges


def get_exon_info(gene: str,
                  sqlite_path_organism: str,
                  supported_transcripts: List[str],
                  filte_r: bool = True) -> DataFrame:
    # Return the exon starts,stops and transcript for a given gene
    gene = gene.upper()
    exons = sqlquery(sqlite_path_organism, "exons")
    transcripts = sqlquery(sqlite_path_organism, "transcripts")
    exons = exons.loc[exons.transcript.isin(transcripts.loc[
        transcripts.gene == gene, "transcript"].values),
                      ["transcript", "exon_start", "exon_stop"]]
    if filte_r:
        return exons[exons.transcript.isin(supported_transcripts)]
    return exons


def genomic_exon_coordinate_ranges(gene: str,
                                   sqlite_path_organism: str,
                                   supported_transcripts: List[str],
                                   filte_r: bool = True) -> DataFrame:
    # create a dictionary with transcript_ids as keys and exon coordinates in tuple (start, stop) as values
    # subtract the minumum start codon position of any exon for the gene
    exon_info = get_exon_info(gene,
                              sqlite_path_organism,
                              supported_transcripts,
                              filte_r=filte_r)
    minimum = min(exon_info.exon_start)
    exon_info["exon_start"] = list(
        map(lambda x: x - minimum, exon_info.exon_start))
    exon_info["exon_stop"] = list(
        map(lambda x: x - minimum, exon_info.exon_stop))
    return exon_info


def genomic_orf_coordinate_ranges(
    # gene,
    sqlite_path_organism: str,
    supported_transcripts: List[str],
    genomic_coordinates: Dict[str, List[Tuple[int, int]]],
    # filter=True
) -> Dict[str, List[Tuple[int, int]]]:
    # Translate the orf structure to genomic coordinates using genomic coordinates from the exons table of sqlite file
    orf_structures = {}
    genomic_orf_structures = {}
    for transcript in supported_transcripts:
        start_stop = get_start_stop_codon_positions(transcript,
                                                    sqlite_path_organism)
        transcript_coordinates = get_exon_coordinates_for_orf(
            transcript, sqlite_path_organism)
        orf_structures[transcript] = get_orf_exon_structure(
            start_stop, transcript_coordinates)

        for region in orf_structures[transcript]:
            for exon in transcript_coordinates:
                if region[1] == exon[1]:
                    difference = region[1] - region[0]

                    genomic_exon = genomic_coordinates[transcript][
                        transcript_coordinates.index(exon)]
                    genomic_orf = (genomic_exon[1] - difference,
                                   genomic_exon[1])
                    if transcript in genomic_orf_structures:
                        genomic_orf_structures[transcript].append(genomic_orf)
                    else:
                        genomic_orf_structures[transcript] = [genomic_orf]

                elif region[0] == exon[0]:
                    difference = region[1] - region[0]

                    genomic_exon = genomic_coordinates[transcript][
                        transcript_coordinates.index(exon)]
                    genomic_orf = (genomic_exon[0],
                                   genomic_exon[0] + difference)
                    if transcript in genomic_orf_structures:
                        genomic_orf_structures[transcript].append(genomic_orf)
                    else:
                        genomic_orf_structures[transcript] = [genomic_orf]

    return genomic_orf_structures


def genomic_junction_positions(
    # gene,
    # sqlite_path_organism,
    # supported_transcripts,
    exons: Dict[str, List[Tuple[int, int]]],
    # filter=True
) -> Dict[str, List[Tuple[int, int]]]:
    genomic_junctions = {}
    for transcript in exons:
        if transcript not in genomic_junctions:
            genomic_junctions[transcript] = []

        number_of_exons = len(exons[transcript])
        for index, _ in enumerate(exons[transcript]):
            if index < number_of_exons - 1:
                genomic_junctions[transcript].append(
                    (exons[transcript][index][1],
                     exons[transcript][index + 1][0]))
    return genomic_junctions


def genomic_junction_scores(
    # gene_name,
    sqlite_path_organism: str,
    sqlite_path_reads: str,
    supported_transcripts: List[str],
    # exons,
    genomic_junctions: Dict[str, List[Tuple[int, int]]],
    # filter=True
) -> Dict[str, Dict[str, List[Tuple[int, int]]]]:

    scores = get_scores_per_exonjunction_for_gene(sqlite_path_organism,
                                                  sqlite_path_reads,
                                                  supported_transcripts)
    paired_dictionary = {}
    for transcript in genomic_junctions:
        sorted_scores = sorted(scores[transcript].keys())
        paired = zip(genomic_junctions[transcript], sorted_scores)

        paired_dictionary[transcript] = paired
    genomic_junctions_scores = {}
    for transcript in paired_dictionary:
        if transcript not in genomic_junctions_scores:
            genomic_junctions_scores[transcript] = {}

        for junction in paired_dictionary[transcript]:
            genomic_junctions_scores[transcript][
                junction[0]] = scores[transcript][junction[1]]

    return genomic_junctions_scores


def get_reads_per_transcript_location(
        transcript_id: str,
        sqlite_path_reads: str) -> Union[None, Dict[int, List[int]]]:
    infile = SqliteDict(sqlite_path_reads)
    return infile[transcript_id]["unambig"] if transcript_id in infile else None
    # print("No unambiguous reads support this gene " + transcript_id)


def get_reads_per_genomic_location(
    gene: str,
    sqlite_path_reads: List[str],
    sqlite_path_organism: str,
    supported_transcripts: List[str],
    genomic_exon_coordinates: DataFrame,
    filte_r: bool = True,
    site: Literal[
        'asite', 'fiveprime',
        'range'] = 'asite'  # NOTE: fiveprime is not required, but keep for reference
) -> Dict[int, Dict[str, List[Tuple[int, int]]]]:
    # get the number of reads supporting each genomic position to be used in the display of support of the
    # supertranscript model. This function takes the reads mapped for each transcript of a gene and uses a combination
    # of genomic and transcriptomic ranges to translate each transcript position to a genomic one.
    trips_splice = TripsSplice(sqlite_path_organism)

    gene_info = trips_splice.get_gene_info(gene)
    gene_info['exon_junctions'] = gene_info['exon_junctions'].apply(
        lambda x: list(map(int, x.split(","))))
    genomic_read_dictionary = {}
    counted_reads = []
    for read_file in sqlite_path_reads:
        infile = SqliteDict(read_file)
        t_gene_info = gene_info[gene_info.transcript.isin(infile.keys())]
        if filte_r:
            t_gene_info = pd.concat([
                t_gene_info,
                gene_info[gene_info.transcript.isin(supported_transcripts)]
            ])
            # ------

        for transcript in gene_info:  # gene_info is a data frame ['transcript', 'exon_junctions', 'sequence']
            if ((filte_r and transcript[0] not in supported_transcripts)
                    or (transcript[0] not in infile)):
                continue
            transcript_read_dictionary = infile[transcript[0]]["unambig"]
            genomic_ranges = genomic_exon_coordinates.loc[
                genomic_exon_coordinayes.gene == transcript[0],
                ['exon_start', 'exon_end']].sort_values('exon_start')

            exon_junctions = list(map(int, transcript[1].split(",")))
            sequence = transcript[2]
            exons = exons_of_transcript(transcript[0], sqlite_path_organism)
            transcript_ranges = get_exon_coordinate_ranges(
                sequence, exons, exon_junctions)

            for length in transcript_read_dictionary:
                for location in transcript_read_dictionary[length]:
                    position = location + (
                        infile["offsets"]["fiveprime"]["offsets"][length]
                        if site == 'asite' else 0)

                    range_counter = 0

                    for exon in transcript_ranges:
                        if position in range(exon[0], exon[1]):
                            difference_between_read_position_and_exon_site = (
                                position - exon[0])

                            genomic_site = genomic_ranges[range_counter][
                                0] + difference_between_read_position_and_exon_site
                            if site == "range":
                                genomic_site = (genomic_site,
                                                genomic_site + length)

                            if genomic_site not in genomic_read_dictionary:
                                genomic_read_dictionary[
                                    genomic_site] = transcript_read_dictionary[
                                        length][location]
                            else:
                                genomic_read_dictionary[
                                    genomic_site] += transcript_read_dictionary[
                                        length][location]
                            if (length, genomic_site) not in counted_reads:
                                counted_reads.append((length, genomic_site))
                        range_counter += 1
    return genomic_read_dictionary


def get_exonjunction_pileup_for_transcript(
        transcript_id: str, sqlite_path_organism: str,
        sqlite_path_reads: str) -> Dict[str, int]:
    # count the number of reads in the read file that span each exon-exon junction. for a given transcript
    # returns a dictionary with junctions as keys and counts as values d
    trips_splice = TripsSplice(sqlite_path_organism)
    transcript_info = trips_splice.get_transcript_info(transcript_id)
    exon_junctions = list(map(int, transcript_info[0][0].split(",")))
    counts = {}
    for read_file in sqlite_path_reads:
        reads = get_reads_per_transcript_location(transcript_id, read_file)
        if not reads:
            continue
        for junction in exon_junctions:
            if junction not in counts:
                counts[junction] = 0
            for read_length in reads:
                for position in reads[read_length]:
                    if position <= junction < position + read_length:
                        counts[junction] += reads[read_length][position]
    return counts


def get_scores_per_exonjunction_for_gene(
        sqlite_path_organism: str, sqlite_path_reads: str,
        supported: List[str]) -> Dict[str, Dict[str, List[Tuple[int, int]]]]:
    # count the reads in the reads file whos p sites lite within the exon sequence
    # returns a dictionary with all unique exons in the gene as keys and counts as values
    pileup = {}
    for trans in supported:
        pileup[trans] = get_exonjunction_pileup_for_transcript(
            trans, sqlite_path_organism, sqlite_path_reads)
    return pileup
