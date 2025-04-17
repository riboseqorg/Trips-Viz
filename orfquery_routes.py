import json
import numpy as np
from time import time
import logging
import os
import random
import sqlite3
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from functools import partial
from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
import polars as pl
from flask import Blueprint, render_template, request
from flask_login import current_user
from sqlitedict import SqliteDict

import config
from core_functions import (
    build_profile,
    build_proteomics_profile,
    fetch_file_paths,
    fetch_files,
    fetch_studies,
    fetch_study_info,
    fetch_user,
    form_filler,
    generate_short_code,
    nuc_to_aa,
)
from sqlqueries_2 import get_table, get_user_id, sqlquery

HEADER = "Gene,transcript,start,stop,start_codon,type,sru,coverage,median_diff,read_density,split,first_diff\n"
MIN_READ_DENSITY = 3
MIN_COVERAGE = 0.33
MAX_ATTEMPTS = 20
CASES_PER_TRAN = 2

# This page is used to detect translated open reading frames
translated_orf_blueprint = Blueprint(
    "orf_translationpage", __name__, template_folder="templates"
)


@translated_orf_blueprint.route("/<organism>/<transcriptome>/orf_translation/")
def orf_translationpage(organism: str, transcriptome: str) -> str:
    # ip = request.environ['REMOTE_ADDR']
    data = form_filler(organism, transcriptome)
    accepted_studies = fetch_studies(data["gwips_info"][0, "organism_id"])
    data["files"] = fetch_files(accepted_studies).to_pandas()
    return render_template("orf_translation.html", template_dict=data)


def tran_to_genome(
    transcriptome_info_dict: Dict[str, Dict[str, Union[str, List[Tuple[int, int]]]]],
    tran: str,
    pos: int,
) -> Union[str, Tuple[str, int]]:
    if tran not in transcriptome_info_dict:
        return None, 0  # This test can be done before coming here
    traninfo = transcriptome_info_dict[tran]
    chrom = traninfo["chrom"]
    strand = traninfo["strand"]
    exons = traninfo["exons"]
    # logging.debug(exons)
    if strand == "+":
        exon_start = 0
        exonlen = 0
        for tup in exons:
            exon_start = tup[0]
            exonlen = tup[1] - tup[0]
            if pos > exonlen:
                pos = (pos - exonlen) - 1
            else:
                break
        genomic_pos = (exon_start + pos) - 1
    else:
        exon_start = 0
        for tup in exons[::-1]:
            exon_start = tup[1]
            exonlen = tup[1] - tup[0]
            if pos > exonlen:
                pos = (pos - exonlen) - 1
            else:
                break
        genomic_pos = (exon_start - pos) + 1
    return "{}_{}".format(chrom, genomic_pos)


def create_aggregate(file_paths_dict: Dict):
    ambig = "unambig"
    file_count = 0
    profile_dict = {}
    file_list = []
    for file_id in file_paths_dict[seq_type]:
        file_list.append(file_id)
        file_count += 1
        sqlite_dict = SqliteDict(
            f"{file_paths_dict[seq_type][file_id]}", autocommit=False, decode=my_decoder
        )
        # sqlite_dict = SqliteDict(file_paths_dict[seq_type][file_id])
        sqlite_db = dict(sqlite_dict)
        sqlite_dict.close()
        # offsets = offset_dict[file_id]
        offsets = {}
        scores = {}
        if seq_type == "riboseq":
            try:
                offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
                scores = sqlite_db["offsets"]["fiveprime"]["read_scores"]
            except Exception:
                pass

        for transcript in sqlite_db:
            if transcript not in profile_dict:
                profile_dict[transcript] = {"riboseq": {}, "proteomics": {}}
            try:
                counts = sqlite_db[transcript]
            except Exception:
                continue
            if seq_type == "riboseq":
                subprofile = build_profile(
                    counts, offsets, ambig, minscore=0.5, scores=scores
                )
            elif seq_type == "proteomics":
                subprofile = build_proteomics_profile(counts, ambig)
            for pos in subprofile:
                try:
                    profile_dict[transcript][seq_type][pos] += subprofile[pos]
                except Exception:
                    profile_dict[transcript][seq_type][pos] = subprofile[pos]
        logging.debug("{} files read".format(file_count))
    outfile = SqliteDict("{}/aggregate_0.5_{}.sqlite".format(study_path, seq_type))
    for tran in profile_dict:
        outfile[tran] = profile_dict[tran]
    outfile["file_list"] = file_list
    outfile.commit()
    outfile.close()


def create_profiles(
    file_paths_dict,
    region,
    ambig,
    # total_files,
    minscore,
):
    # logging.debug(file_paths_dict)
    # If
    # This will be populated with the users chosen file_ids and passed to the table, so that the trips link can use these files aswell.

    profiles = {}
    # print(file_paths_dict, region)

    for row in file_paths_dict.iter_rows(named=True):

        sqlite_db = SqliteDict(
            row["path"],
            autocommit=False,
            #        decode=my_decoder,
        )
        print(row)
        if row["file_type"] == "riboseq":
            offsets_5p_offsetsNscores = sqlite_db["offsets"]["fiveprime"]
            if minscore:
                offsets_5p_offsetsNscores = offsets_5p_offsetsNscores.filter(
                    pl.col("read_scores") >= int(minscore)
                )

        for transcript in set(region["transcript"]) & set(sqlite_db):
            # print(transcript, row)
            counts = sqlite_db[transcript]
            # print(offsets_5p_offsetsNscores, "==========", counts)
            subprofile = pl.DataFrame()
            if row["file_type"] == "riboseq":
                subprofile = build_profile(counts, offsets_5p_offsetsNscores, ambig)
            elif row["file_type"] == "proteomics":
                subprofile = build_proteomics_profile(counts)
            if subprofile.shape[0]:
                profiles[transcript] = {row["file_type"]: subprofile}

    # print(profiles)
    return profiles


def extract_features(start, stop, profile):
    # selected_range = np.array(range(start - 9, stop + 12))
    selected_range = np.array(range(start - 10, stop + 13))
    profile = (
        pd.DataFrame({"pos": list(profile.keys()), "counts": list(profile.values())})
        .merge(pd.DataFrame({"pos": selected_range}), how="right")
        .fillna(0)
    )
    inframe_values = profile.loc[profile.pos % 3 == 0, "counts"].values
    minusone_values = profile.loc[profile.pos % 3 == 2, "counts"].values
    plusone_values = profile.loc[profile.pos % 3 == 1, "counts"].values
    oof_val = np.maximum(minusone_values, plusone_values)[
        4:-4
    ]  # Pairwise compration and trimming
    if_val = inframe_values[4:-4]
    diff_list = if_val / (oof_val + if_val + 1)
    standard_cov = sum(if_val > 0)
    if_count = sum(if_val > (oof_val + 1))
    # Till here

    if_cov = 0.0
    if_len = stop - start

    max_inframe = max(
        inframe_values[5:-5]
    )  # What if sequences is smaller that 10 nuc. I know it is not possible
    for i in range(5, len(inframe_values) - 5):
        if inframe_values[i] == max_inframe:
            inframe_values[i] = 0
            break
    sru = (sum(inframe_values[4:8])) / (sum(inframe_values[:8]) + 4.0)
    minusone_sum = sum(minusone_values[4:-4]) * 1.0
    plusone_sum = sum(plusone_values[4:-4]) * 1.0
    if_count = 0.0
    standard_cov = 0.0
    first_diff = False
    diff_list = []
    for i in range(4, len(inframe_values) - 4):
        if_val = inframe_values[i]
        oof_val = max(minusone_values[i], plusone_values[i])
        diff = if_val * 1.0 / (oof_val + if_val + 1)
        if not first_diff:
            first_diff = diff
        diff_list.append(diff)
        if if_val > 0:
            standard_cov += 1
        if if_val > (oof_val + 1):
            if_count += 1
    if_cov = if_count / (len(inframe_values) - 8)
    if_values = inframe_values[4:-4]
    if_values_sum = float(sum(if_values))
    read_density = if_values_sum / float(if_len / 3)
    mod_inframe_sum = sum(sorted(if_values)[:-1])
    standard_coverage = standard_cov / float(if_len / 3)
    median_diff = mod_inframe_sum / (
        (max(minusone_sum, plusone_sum) + 1) + mod_inframe_sum
    )
    if_values_len = len(if_values) / 2
    split = sum(if_values[:if_values_len]) / (if_values_sum + 1)
    return [
        sru,
        if_cov,
        median_diff,
        read_density,
        standard_coverage,
        split,
        first_diff,
    ]


def extract_values(
    regions, data, tran_gene_dict, selected_seq_types, profile_dict, all_cases
):
    logging.debug("extract values called")
    best_high_frame = 0
    best_low_frame = 0
    best_start_score = 0
    best_stop_score = 0
    best_inframe_cov = 0
    all_values = []
    tot_loc = 0
    regions = regions.filter(pl.col("transcript").is_in(profile_dict))
    best_values = {}

    for row in regions.iter_rows(named=True):
        locus = row["transcript"]
        stop = row["stop"]
        start = row["start"]
        best_values[locus] = {
            "start": -1,
            "high_frame_count": 1,
            "low_frame_count": 1,
            "start_score": -10,
            "stop_score": 1,
            "final_score": -10000,
            "coverage": 0,
            "length": 0,
            "transcript": locus,
            "stop": 0,
            "proteomics_count": 0,
            "ratio": 0,
            "inframe_count": 0,
            "start_ratio": 0,
            "stop_ratio": 0,
            "high_ratio": 0,
            "low_ratio": 0,
            "stop_geo": 0,
            "start_geo": 0,
        }
        best_values["transcript"] = row["transcript"]
        if best_values["start"] == -1:
            best_values["start"] = start
            best_values["stop"] = stop
        gene = tran_gene_dict[row["transcript"]]
        if "riboseq" in selected_seq_types and "proteomics" in selected_seq_types:
            profile = profile_dict[row["transcript"]]["riboseq"]
            prot_profile = profile_dict[row["transcript"]]["proteomics"]
        elif "proteomics" not in selected_seq_types:
            profile = profile_dict[row["transcript"]]["riboseq"]
            prot_profile = {}
        else:
            profile = profile_dict[row["transcript"]]["proteomics"]
            prot_profile = {}
        transcriptome_stop = regions[locus][stop][start]["stop"]
        orftype = regions[locus][stop][start]["orftype"]
        length = (transcriptome_stop - start) + 1

        if_cov = 0.0
        if_len = (stop - start) / 3
        for p in set(np.arange(start + 1, transcriptome_stop, 3) - 1) & set(
            prot_profile
        ):
            proteomics_count += prot_profile[p - 1]

        if_cov = 0.0
        if_cov_total = 0
        inframe_values = []
        minusone_values = []
        plusone_values = []
        for i in range(start - 9, stop + 12, 3):
            if i - 1 in profile:
                minusone_values.append(profile[i - 1])
            else:
                minusone_values.append(0)
            if i in profile:
                inframe_values.append(profile[i])
            else:
                inframe_values.append(0)
            if i + 1 in profile:
                plusone_values.append(profile[i + 1])
            else:
                plusone_values.append(0)
        for x in range(4, len(inframe_values) - 4):
            if max(inframe_values[x], minusone_values[x], plusone_values[x]) > 0:
                if_cov_total += 1
                min_max = np.sort([minusone_values[x], plusone_values[x]])

                if "highest_frame_diff_check" in data:
                    if inframe_values[x] > min_max[1]:
                        if_cov += 1
                else:
                    if inframe_values[x] > min_max[0]:
                        if_cov += 1

        if_cov *= (100.0 / if_cov_total) if if_cov_total else 0

        start_score_raw = (sum(inframe_values[4:8]) + 1.0) / (
            sum(inframe_values[:4]) + 1
        )
        start_score = np.log(start_score_raw)
        stop_score_raw = (sum(inframe_values[-8:-4]) + 1.0) / (
            sum(inframe_values[-4:]) + 1
        )
        stop_score = np.log(stop_score_raw)
        inframe_sum = sum(inframe_values[4:-4])
        lowframe__highframe_sum = np.sort(
            sum(minusone_values[4:-4]), sum(plusone_values[4:-4])
        )
        non_zero_counts = 0.0
        for i in inframe_values[4:-4]:
            if i:
                non_zero_counts += 1
        read_density = non_zero_counts / (length / 3.0)

        high_frame_count = inframe_sum - lowframe__highframe_sum[1]
        lowest_frame_count = inframe_sum - lowframe__highframe_sum[0]

        # Ratios
        start_ratio = sum(inframe_values[4:8]) / (sum(inframe_values[:4]) + 1.0)
        stop_ratio = sum(inframe_values[-8:-4]) / (sum(inframe_values[-4:]) + 1.0)
        low__high_ratio = inframe_sum / (lowframe__highframe_sum + 1.0)
        high_ratio = inframe_sum / (highframe_sum + 1.0)
        low_ratio = inframe_sum / (lowframe_sum + 1.0)
        final_score_values = []
        if "start_increase_check" in data:
            final_score_values.append(start_score)

        if "stop_decrease_check" in data:
            final_score_values.append(stop_score)

        if "lowest_frame_diff_check" in data:
            final_score_values.append(lowest_frame_count)

        if "highest_frame_diff_check" in data:
            final_score_values.append(high_frame_count)
        if "coverage_check" in data:
            final_score_values.append(if_cov)
        final_score_values.append(float(inframe_sum) / float(if_len))

        # TO DO: Instead of just summing these, they should be normalised over the current best value and then compared. The way it works now, coverage
        # has very little affect on the final outcome as it's a number below one. Instead normalise everything to whatever is highest between the current
        # value and the best value, e.g if current start score is 50 and best start score is 100, current start becomes 0.5 and best becomes 1,
        # another e.g if current coverage is 0.8 and best coverage is 0.2, current coverage becomes 1, best coverage becomes 0.25
        final_score = sum(final_score_values)
        if not all_cases:
            if start_score > best_values["start_score"] or (
                start_score == best_values["start_score"]
                and start < best_values["start"]
            ):
                best_values["start"] = start
                best_values["stop"] = transcriptome_stop
                best_values["transcript"] = row["transcript"]
                best_values["length"] = length
                best_values["high_frame_count"] = high_frame_count
                best_values["low_frame_count"] = lowest_frame_count
                best_values["start_score"] = start_score
                best_values["stop_score"] = stop_score
                best_values["final_score"] = final_score
                best_values["coverage"] = if_cov
                best_values["proteomics_count"] = proteomics_count
                best_values["read_density"] = read_density
                best_values["inframe_count"] = float(inframe_sum) / float(if_len)
                best_values["start_ratio"] = start_ratio
                best_values["stop_ratio"] = stop_ratio
                best_values["high_ratio"] = high_ratio
                best_values["low_ratio"] = low_ratio
                if high_frame_count > best_high_frame:
                    best_high_frame = high_frame_count
                if lowest_frame_count > best_low_frame:
                    best_low_frame = lowest_frame_count
                if stop_score > best_stop_score:
                    best_stop_score = stop_score
                if start_score > best_start_score:
                    best_start_score = start_score
                if if_cov > best_inframe_cov:
                    best_inframe_cov = if_cov
        else:
            all_values.append(
                [
                    gene,
                    row["transcript"],
                    start,
                    transcriptome_stop,
                    length,
                    high_frame_count,
                    lowest_frame_count,
                    stop_score,
                    start_score,
                    if_cov,
                    float(inframe_sum) / float(if_len),
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    orftype,
                    proteomics_count,
                    read_density,
                    start_ratio,
                    stop_ratio,
                    high_ratio,
                    low_ratio,
                ]
            )
    if not all_cases:
        all_values.append(
            [
                gene,
                best_values["transcript"],
                best_values["start"],
                best_values["stop"],
                best_values["length"],
                best_values["high_frame_count"],
                best_values["low_frame_count"],
                best_values["stop_score"],
                best_values["start_score"],
                best_values["coverage"],
                best_values["inframe_count"],
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                orftype,
                best_values["proteomics_count"],
                best_values["read_density"],
                best_values["start_ratio"],
                best_values["stop_ratio"],
                best_values["high_ratio"],
                best_values["low_ratio"],
            ]
        )
    len_all_rows = float(len(all_values))
    logging.debug("LENGTH OF ALL ROWS {}".format(len_all_rows))
    if len_all_rows == 0:
        return None

    sorted_all_values = sorted(all_values, key=lambda x: x[5], reverse=True)
    rank = 1
    prev_value = sorted_all_values[0][5]
    for row in sorted_all_values:
        curr_value = row[5]
        if curr_value != prev_value:
            rank += 1
            prev_value = curr_value
        row[11] = rank
    sorted_all_values = sorted(all_values, key=lambda x: x[6], reverse=True)
    # FDR correction
    # sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    rank = 1
    prev_value = sorted_all_values[0][6]
    for row in sorted_all_values:
        curr_value = row[6]
        if curr_value != prev_value:
            rank += 1
            prev_value = curr_value
        row[12] = rank

    sorted_all_values = sorted(all_values, key=lambda x: x[7], reverse=True)
    # FDR correction
    # sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    rank = 1
    prev_value = sorted_all_values[0][7]
    for row in sorted_all_values:
        curr_value = row[7]
        if curr_value != prev_value:
            rank += 1
            prev_value = curr_value
        row[13] = rank

    sorted_all_values = sorted(all_values, key=lambda x: x[8], reverse=True)
    # FDR correction
    # sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    rank = 1
    prev_value = sorted_all_values[0][8]
    for row in sorted_all_values:
        curr_value = row[8]
        if curr_value != prev_value:
            rank += 1
            prev_value = curr_value
        row[14] = rank

    sorted_all_values = sorted(all_values, key=lambda x: x[9], reverse=True)
    rank = 1
    prev_value = sorted_all_values[0][9]
    for row in sorted_all_values:
        curr_value = row[9]
        if curr_value != prev_value:
            rank += 1
            prev_value = curr_value
        row[15] = rank

    sorted_all_values = sorted(all_values, key=lambda x: x[10], reverse=True)
    rank = 1
    prev_value = sorted_all_values[0][10]
    for row in sorted_all_values:
        curr_value = row[10]
        if curr_value != prev_value:
            rank += 1
            prev_value = curr_value
        row[16] = rank

    for row in all_values:
        normalised_score_values = []
        # inframe count
        normalised_score_values.append(row[16])
        if "coverage_check" in data:
            normalised_score_values.append(row[15])
        if "start_increase_check" in data:
            # if row[13] >0.05:
            # continue
            normalised_score_values.append(row[14])
        if "stop_decrease_check" in data:
            # if row[12] >0.05:
            # continue
            normalised_score_values.append(row[13])
        if "lowest_frame_diff_check" in data:
            # if row[11] >0.05:
            # continue
            normalised_score_values.append(row[12])
        if "highest_frame_diff_check" in data:
            # if row[10] >0.05:
            # continue
            normalised_score_values.append(row[11])
        normalised_score = sum(normalised_score_values)
        # logging.debug("normalised_score_values, normalised score", normalised_score_values, normalised_score)
        row.append(round(normalised_score, 2))
    logging.debug("sorting all values")
    sorted_all_values = sorted(all_values, key=lambda x: x[-1], reverse=False)
    logging.debug("LENGTH OF ALL SORTED ROWS {}".format(len(sorted_all_values)))
    final_rank = 1
    for tup in sorted_all_values:
        tup[17] = final_rank
        final_rank += 1
    return sorted_all_values


def write_to_file(
    sorted_all_values, filename, sequence_dict, organism, transcriptome, file_string
):
    returnstr = "Table|"
    tmp_result_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename), "w")
    tmp_result_file.write(
        "Gene,Tran,Start,Stop,Length,Global_Rank,Type,Trips-viz link,Start Codon,Highframe rank,Highframe value,Lowframe rank,Lowframe value,Stop rank,Stop value,Start rank,Start value,Coverage rank,Coverage value,Inframe Count Rank,Inframe Count Value,Amino acid sequence,Proteomics count,Read density\n"
    )
    tup_count = 0
    # logging.debug("writing to file",len(sorted_all_values))
    for tup in sorted_all_values:
        gene = tup[0]
        transcript = tup[1]
        start = tup[2]
        stop = tup[3]
        length = tup[4]
        global_rank = tup[17]
        orftype = tup[18]
        proteomics_count = tup[19]
        read_density = tup[20]
        try:
            seq = sequence_dict[transcript][start - 1 : stop]
        except Exception:
            seq = ""
        start_codon = seq[:3]
        # Get amino acid sequence of this ORF, excluding the stop codon.
        while len(seq) % 3 != 0:
            seq = seq[:-1]

        aa_seq = nuc_to_aa(seq[:-3])
        trips_link = (
            '<a href="https://trips.ucc.ie/'
            + organism
            + "/"
            + transcriptome
            + "/interactive_plot/?tran="
            + transcript
            + "&hili="
            + str(start)
            + "_"
            + str(stop)
            + "&files="
            + file_string
            + '" target="_blank_" >View on trips-viz</a>'
        )
        tmp_result_file.write(
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
                gene,
                transcript,
                start,
                stop,
                length,
                global_rank,
                orftype,
                trips_link,
                start_codon,
                tup[11],
                tup[5],
                tup[12],
                tup[6],
                tup[13],
                tup[7],
                tup[14],
                tup[8],
                tup[15],
                tup[9],
                tup[16],
                tup[10],
                aa_seq,
                proteomics_count,
                read_density,
            )
        )
        if tup_count < 1000:
            returnstr += "{},{},{},{},{},{},{},NULL,NULL,NULL,NULL,{}.,/".format(
                gene, transcript, start, stop, length, global_rank, orftype, trips_link
            )
        tup_count += 1
    return returnstr


def orfquery(data):
    # data = request.args.to_dict()
    t1 = time()
    user, logged_in = fetch_user()

    # print(data)
    logging.debug("orfquery called")
    owner = get_table("organisms").filter(
        (pl.col("organism_name") == data["organism"])
        & (pl.col("transcriptome_list") == data["transcriptome"])
    )[0, "owner"]
    # TODO: Check transcriptome list option

    files = []
    for key in data:
        if key.startswith(
            "file_{}".format(data["file_type"])
        ):  # NOTE:Might need to apply for proteomics
            # print(key)
            files.append(int(key.split("__")[1]))
    data["file_ids"] = files
    # print(files)

    total_files = len(data["file_ids"])
    if total_files > 100:
        return "A maximum of 100 files can be selected on this page, currently there are {} selected".format(
            total_files
        )

    file_paths_dict = fetch_file_paths(data)
    # print(file_paths_dict)
    # Find out which studies have all files of a specific sequence type selected (to create aggregates)

    # full_studies = []

    # for study_path in aggregate_dict:
    #     create_aggregate(
    #         file_paths_dict
    #     )

    # logging.debug("Full studies {}".format(full_studies))
    # return str(full_studies)
    # return str(all_study_ids)

    # html_args = data["html_args"]
    returnstr = ""

    start_codons = []

    # Used to extract columns from csv for neural neural_net
    for s_codon in ["aug", "cug", "gug"]:  # , "none"
        if "sc_" + s_codon in data:
            start_codons.append(s_codon.upper())
    if not start_codons:
        return "Error no start codon types selected"
    data["start_codons"] = start_codons
    logging.debug("start codons {}".format(start_codons))

    all_cases = False
    if int(data["max_len"]) > 10_000:
        all_cases = True
        data["max_len"] = 10_000

    user_defined_transcripts = []
    # tran_list is a radio button, user can choose between principal, all or a custom list
    # custom_tran_list is the actual comma seperated list of transcripts that user would enter should they choose the custom option in tranlist

    filtered_transcripts = pl.DataFrame()

    # if "saved_check" in data and current_user.is_authenticated:
    #     user_id = get_user_id(current_user.name)
    #     filtered_transcripts = (
    #         get_table("users_saved_cases")
    #         .filter(
    #             (pl.col("user_id") == user_id)
    #             & (pl.col("organism") == data["organism"])
    #         )
    #         .select("tran", "stop")
    #     )

    feature_list = ["Type"]
    if "start_increase_check" in data:
        feature_list.append("Start value")

    if "stop_decrease_check" in data:
        feature_list.append("Stop value")

    if "coverage_check" in data:
        feature_list.append("Coverage value")

    if "lowest_frame_diff_check" in data:
        feature_list.append("Lowframe value")

    if "highest_frame_diff_check" in data:
        feature_list.append("Highframe value")
    # feature_list.append("Inframe Count Value")
    # if html_args["user_short"] == "None":
    #     short_code = generate_short_code(
    #         data, organism, data["transcriptome"], "orf_translation"
    #     )
    # else:
    #     short_code = html_args["user_short"]
    short_code = "anml"

    filename = short_code + ".csv"
    filepath = "{}/static/tmp/{}".format(config.SCRIPT_LOC, filename)
    # if os.path.isfile(filepath):
    #     logging.debug(
    #         "File exists {}/static/tmp/{}".format(config.SCRIPT_LOC, filename)
    #     )
    #     returnstr = (
    #         pd.read_csv(filepath)
    #         .head(1000)
    #         .to_html(index=False, classes="table table-striped", table_id="table")
    #     )
    #     return returnstr  # TODO: add user short
    # else:
    #     logging.debug(
    #         "File does not exists {}/static/tmp/{}".format(config.SCRIPT_LOC, filename)
    #     )

    if data["custom_tran_list"]:
        data["custom_tran_list"] = data["custom_tran_list"].split(",")

    tran_gene_dict = {}
    # structure of orf dict is transcript[stop][start] = {"length":x,"score":0,"cds_cov":0} each stop can have multiple starts
    accepted_orf_dict = {}

    # print("jasgdjasgfjgsdfjgsdjfgshj")
    if owner:
        sqlite_path = "{0}/{1}/{2}/{2}.{3}.sqlite".format(
            config.SCRIPT_LOC,
            config.ANNOTATION_DIR,
            data["organism"],
            data["transcriptome"],
        )
        if not os.path.isfile(sqlite_path):
            return "Cannot find annotation file {}.{}.sqlite".format(
                data["organism"], data["transcriptome"]
            )
    else:
        sqlite_path = "{0}/transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
            config.UPLOADS_DIR, owner, data["organism"], data["transcriptome"]
        )
    transcripts = sqlquery(sqlite_path, "transcripts")

    # traninfo_dict = {}
    # if data["output"] == "nnet":
    #     pass

    if data["tranlist"] == "prin_trans":
        transcripts = transcripts.filter(pl.col("principal") == 1)
    elif data["tranlist"] == "custom_trans":
        transcripts = transcripts.filter(
            pl.col("transcript").is_in(data["custom_tran_list"])
        )
    ambig = False

    # exons = (
    #     sqlquery(sqlite_path, "exons")
    #     .filter(pl.col("transcript").is_in(transcripts["transcript"]))
    #     .join(transcripts, on="transcript")
    # ).select("transcript", "exon_start", "exon_stop", "strand")

    # tmp_exon = []
    # for transcript in exons["transcript"].unique():
    #     t_exon = exons.filter(pl.col("transcript") == transcript)
    #     strand = t_exon[0, "strand"]
    #     if strand == "+":
    #         t_exon = t_exon.sort("exon_start")
    #     else:
    #         t_exon = t_exon.sort("exon_start", descending=True)
    #     begin = 0
    #     for tc in t_exon.iter_rows(named=True):
    #         exon_len = tc["exon_stop"] - tc["exon_start"]
    #         tmp_exon.append(
    #             (
    #                 transcript,
    #                 tc["exon_start"],
    #                 tc["exon_stop"],
    #                 strand,
    #                 begin,
    #                 begin + exon_len,
    #             )
    #         )
    #         begin += exon_len
    # tmp_exon = pl.DataFrame(
    #     tmp_exon, columns=["transcript", "start", "stop", "strand", "begin", "end"]
    # )

    # exons = tmp_exon.copy()
    # print(tmp_exon)

    logging.debug("building transcriptom info dict")
    # Holds a list of all transcripts in accepted_orf_dict
    region = sqlquery(sqlite_path, data["region"]).filter(
        (pl.col("cds_coverage") >= int(data["min_cds"]))
        & (pl.col("cds_coverage") <= float(data["max_cds"]))
        & (pl.col("length") >= int(data["min_len"]))
        & (pl.col("length") <= int(data["max_len"]))
        & pl.col("transcript").is_in(transcripts["transcript"])
    )
    # transcripts = transcripts.filter(pl.col("transcript").is_in(table["transcript"]))
    # print(region)
    # print(transcripts)
    # tran2gnome = partial(tran_to_genome, transcriptome_info_dict)
    # table = table.filter(pl.col("start_codon").is_in(start_codons)).with_columns(
    #     # Apply tran to dict
    #     stopx=pl.col("stop_codon").apply(tran2gnome)
    # )

    logging.debug("for row in result")
    # rows = 0
    #     for row in result:
    #         rows += 1
    #         if rows % 1000 == 0:
    #             logging.debug("Rows: {}".format(rows))
    #             prog_count = 100 + (100 * (rows / total_trans))
    #         # logging.debug("row", row)
    #         transcript = str(row[0])
    #         start_codon = str(row[1])
    #         length = row[2]
    #         cds_cov = row[3]
    #         start = row[4]
    #         # Orfs with multiple potential starts will be grouped by start codon (so only one potential start is reported)
    #         # If the user has selected all transcripts then instead group by the genomic stop codon co-ordinates (so only one ORF will be reported,
    #         # even if it occurs on multiple transcript isoforms)
    #         stop = tran_to_genome(transcript, row[5], transcriptome_info_dict)
    #         locus = tran_gene_dict[transcript]
    #         if transcript in filtered_transcripts:
    #             if stop in filtered_transcripts[transcript]:
    #                 continue
    #         transcriptome_stop = row[5]
    #         if locus not in accepted_orf_dict:
    #             accepted_orf_dict[locus] = {}
    #         if stop not in accepted_orf_dict[locus]:
    #             accepted_orf_dict[locus][stop] = {}
    #         if transcript not in accepted_transcript_list:
    #             accepted_transcript_list.append(transcript)
    #         accepted_orf_dict[locus][stop][start] = {
    #             "length": length,
    #             "score": 0,
    #             "cds_cov": cds_cov,
    #             "start_codon": start_codon,
    #             "orftype": table_name,
    #             "stop": transcriptome_stop,
    #             "transcript": transcript,
    #         }

    # logging.debug("accepted orf dict", accepted_orf_dict)
    logging.debug("accepted orf dict built")
    # Now build a profile for every transcript in accepted_transcripts

    if file_paths_dict.is_empty():
        return "Error no files selected"
    if file_paths_dict.filter(pl.col("file_type") == "riboseq").is_empty() and (
        "te_check" in data
    ):
        del data["te_check"]

    total_files = file_paths_dict.shape[0]
    selected_seq_types = set(file_paths_dict["file_type"])

    # if data["output"] == "nnet":
    #     for transcript in traninfo_dict:
    #         if transcript not in accepted_transcript_list:
    #             accepted_transcript_list.append(transcript)

    profile_dict = create_profiles(file_paths_dict, region, ambig, data["minscore"])
    print(time() - t1)
    print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    logging.debug("profile dict built\nextracting values")
    if data["output"] != "nnet":
        sorted_all_values = extract_values(
            region,
            data,
            tran_gene_dict,
            selected_seq_types,
            profile_dict,
            all_cases,
        )
        if not sorted_all_values:
            return "No results, try making filters less restrictive"
    # else:
    #     create_training_set(profile_dict, traninfo_dict, short_code, data["min_len"])
    #     create_test_set(
    #         profile_dict,
    #         accepted_orf_dict,
    #         traninfo_dict,
    #         short_code,
    #         data["min_len"],
    #         sequence_dict,
    #         region,
    #     )
    #     returnstr = neural_net(
    #         short_code,
    #         ["type", "coverage", "median_diff", "first_diff"],
    #         organism,
    #         transcriptome,
    #         file_string,
    #     )

    if data["output"] != "nnet":
        # NOTE: It will never come here
        logging.debug("Writing to file")
        returnstr = write_to_file(
            sorted_all_values,
            filename,
            sequence_dict,
            organism,
            transcriptome,
            file_string,
        )

    logging.debug("creating returnstr")
    returnstr += "|"
    returnstr += "Min,0,0,0,0,0,0.,/"
    returnstr += "10th_percentile,0,0,0,0,0,0.,/"
    returnstr += "25th_percentile,0,0,0,0,0,0.,/"
    returnstr += "50th_percentile,0,0,0,0,0,0.,/"
    returnstr += "Max,0,0,0,0,0,0.,/"
    returnstr += "|{}".format(filename)
    returnstr += "|{}".format(
        str(data["file_list"]).strip("[]").replace("'", "").replace(" ", "")
    )  # file_list is empty after passin through generate short_Code need to make a copy of it beforehand
    returnstr += "|{}".format(short_code)
    total_time = time.time() - start_time
    logging.debug("sending email")
    # If the job was > 5 minutes and the user is using an email address, send an email to say the job is done
    if current_user.is_authenticated:
        if total_time > 300 and logged_in:
            try:
                fromaddr = "ribopipe@gmail.com"
                toaddr = user
                msg = MIMEMultipart()
                msg["From"] = fromaddr
                msg["To"] = toaddr
                msg["Subject"] = "Trips-Viz job completion"
                msg.attach(
                    MIMEText(
                        "Your Trips-Viz job is complete: https:trips.ucc.ie/short/{}".format(
                            short_code
                        )
                    )
                )
                server = smtplib.SMTP("smtp.gmail.com", 587)
                server.starttls()
                # TODO, move this to the config file
                server.login(fromaddr, "Ribosome")
                text = msg.as_string()
                logging.debug("sending now")
                server.sendmail(fromaddr, toaddr, text)
                server.quit()
            except Exception:
                pass
    logging.debug("Returning result")
    return returnstr
