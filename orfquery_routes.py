from typing import Dict, Tuple, List, Union
from flask import Blueprint, render_template, request
import sqlite3
import polars as pl
from sqlitedict import SqliteDict
import os
import time
import logging
import config
from core_functions import (fetch_studies, fetch_files, fetch_study_info,
                            fetch_file_paths, generate_short_code, form_filler,
                            build_profile, build_proteomics_profile, nuc_to_aa,
                            fetch_user)
from sqlqueries_2 import get_table, get_user_id
import json
import pandas as pd
from flask_login import current_user
import random

import numpy as np
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

HEADER = "Gene,transcript,start,stop,start_codon,type,sru,coverage,median_diff,read_density,split,first_diff\n"
MIN_READ_DENSITY = 3
MIN_COVERAGE = 0.33
MAX_ATTEMPTS = 20
CASES_PER_TRAN = 2

# This page is used to detect translated open reading frames
translated_orf_blueprint = Blueprint("orf_translationpage",
                                     __name__,
                                     template_folder="templates")


@translated_orf_blueprint.route('/<organism>/<transcriptome>/orf_translation/')
def orf_translationpage(organism: str, transcriptome: str) -> str:
    # ip = request.environ['REMOTE_ADDR']
    data = form_filler(organism, transcriptome)
    accepted_studies = fetch_studies(data["gwips_info"][0, "organism_id"])
    data["files"] = fetch_files(accepted_studies).to_pandas()
    return render_template('orf_translation.html', template_dict=data)


def tran_to_genome(
    tran: str, pos: int,
    transcriptome_info_dict: Dict[str, Dict[str, Union[str, List[Tuple[int,
                                                                       int]]]]]
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


def create_aggregate(file_paths_dict: Dict, study_path, seq_type):
    ambig = "unambig"
    file_count = 0
    profile_dict = {}
    file_list = []
    for file_id in file_paths_dict[seq_type]:
        file_list.append(file_id)
        file_count += 1
        sqlite_dict = SqliteDict(f"{file_paths_dict[seq_type][file_id]}",
                                 autocommit=False,
                                 decode=my_decoder)
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
                subprofile = build_profile(counts,
                                           offsets,
                                           ambig,
                                           minscore=0.5,
                                           scores=scores)
            elif seq_type == "proteomics":
                subprofile = build_proteomics_profile(counts, ambig)
            for pos in subprofile:
                try:
                    profile_dict[transcript][seq_type][pos] += subprofile[pos]
                except Exception:
                    profile_dict[transcript][seq_type][pos] = subprofile[pos]
        logging.debug("{} files read".format(file_count))
    outfile = SqliteDict("{}/aggregate_0.5_{}.sqlite".format(
        study_path, seq_type))
    for tran in profile_dict:
        outfile[tran] = profile_dict[tran]
    outfile["file_list"] = file_list
    outfile.commit()
    outfile.close()


def create_profiles(
        file_paths_dict,
        accepted_transcript_list,
        ambig,
        # total_files,
        minscore):
    # logging.debug(file_paths_dict)
    # If
    file_count = 0
    # This will be populated with the users chosen file_ids and passed to the table, so that the trips link can use these files aswell.
    file_string = ""
    ribo_study_string = "&ribo_studies="
    proteomics_study_string = "&proteomics_studies="
    seq_types = ["riboseq", "proteomics"]

    tot_file_ids = 0.0
    file_count = 0
    for seq_type in seq_types:
        tot_file_ids += len(file_paths_dict[seq_type])

    for seq_type in seq_types:
        if seq_type not in file_paths_dict:
            continue
        for file_id in file_paths_dict[seq_type]:
            file_count += 1

            sqlite_db = SqliteDict(f"{file_paths_dict[seq_type][file_id]}",
                                   autocommit=False,
                                   decode=my_decoder)
            if type(file_id) == str:
                if "STUDY" in file_id:
                    if seq_type == "riboseq":
                        ribo_study_string += "{};".format(
                            file_id.replace("STUDY_", ""))
                    else:
                        proteomics_study_string += "{};".format(
                            file_id.replace("STUDY_", ""))
                    for transcript in accepted_transcript_list:
                        if transcript not in profile_dict:
                            profile_dict[transcript] = {
                                "riboseq": {},
                                "proteomics": {}
                            }
                        if transcript in sqlite_db:
                            subprofile = sqlite_db[transcript][seq_type]
                            for pos in subprofile:
                                try:
                                    profile_dict[transcript][seq_type][
                                        pos] += subprofile[pos]
                                except Exception:
                                    profile_dict[transcript][seq_type][
                                        pos] = subprofile[pos]
            else:
                file_string += "{};".format(file_id)
                offsets = {}
                scores = {}
                if seq_type == "riboseq":
                    try:
                        offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
                    except Exception:
                        pass
                    try:
                        scores = sqlite_db["offsets"]["fiveprime"][
                            "read_scores"]
                    except Exception:
                        pass

                for transcript in accepted_transcript_list:
                    if transcript not in profile_dict:
                        profile_dict[transcript] = {
                            "riboseq": {},
                            "proteomics": {}
                        }
                    try:
                        counts = sqlite_db[transcript]
                    except Exception:
                        continue
                    if seq_type == "riboseq":
                        subprofile = build_profile(counts, offsets, ambig,
                                                   minscore, scores)
                    elif seq_type == "proteomics":
                        subprofile = build_proteomics_profile(counts, ambig)
                    for pos in subprofile:
                        try:
                            profile_dict[transcript][seq_type][
                                pos] += subprofile[pos]
                        except Exception:
                            profile_dict[transcript][seq_type][
                                pos] = subprofile[pos]
            sqlite_db.close()
    # logging.debug("profile_dict", profile_dict["ENST00000233893"])
    file_string += ribo_study_string
    file_string += proteomics_study_string
    return (profile_dict, file_string)


def extract_features(start, stop, profile):
    # selected_range = np.array(range(start - 9, stop + 12))
    selected_range = np.array(range(start - 10, stop + 13))
    profile = pd.DataFrame({
        "pos": list(profile.keys()),
        "counts": list(profile.values())
    }).merge(pd.DataFrame({"pos": selected_range}), how="right").fillna(0)
    inframe_values = profile.loc[profile.pos % 3 == 0, "counts"].values
    minusone_values = profile.loc[profile.pos % 3 == 2, "counts"].values
    plusone_values = profile.loc[profile.pos % 3 == 1, "counts"].values
    oof_val = np.maximum(
        minusone_values,
        plusone_values)[4:-4]  # Pairwise compration and trimming
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
    minusone_sum = sum(minusone_values[4:-4]) * 1.
    plusone_sum = sum(plusone_values[4:-4]) * 1.
    if_count = 0.0
    standard_cov = 0.0
    first_diff = False
    diff_list = []
    for i in range(4, len(inframe_values) - 4):
        if_val = inframe_values[i]
        oof_val = max(minusone_values[i], plusone_values[i])
        diff = if_val * 1. / (oof_val + if_val + 1)
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
        (max(minusone_sum, plusone_sum) + 1) + mod_inframe_sum)
    if_values_len = len(if_values) / 2
    split = sum(if_values[:if_values_len]) / (if_values_sum + 1)
    return [
        sru, if_cov, median_diff, read_density, standard_coverage, split,
        first_diff
    ]


def create_training_set(profile_dict, traninfo_dict, short_code, MINLEN):
    filename = short_code + ".training.csv"
    outfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename), "w")
    outfile.write(HEADER)
    # Positive cases
    acceptable_trans = []
    tot_positive = 0
    tot_negative = 0
    # Positives
    for transcript in traninfo_dict:
        if not traninfo_dict[transcript]["cds_start"]:
            continue
        acceptable_cases = 0
        attempts = 0
        while acceptable_cases < CASES_PER_TRAN and attempts < MAX_ATTEMPTS and tot_positive < 30000:
            attempts += 1
            gene = traninfo_dict[transcript]["gene"]
            cds_start = traninfo_dict[transcript]["cds_start"]
            cds_stop = traninfo_dict[transcript]["cds_stop"]
            if cds_stop - cds_start < ((MINLEN * 2) + 100):
                continue
            start = random.randint(cds_start + 100, cds_stop - MINLEN)
            while start % 3 != cds_start % 3:
                start = random.randint(cds_start + 100, cds_stop - MINLEN)
            stop = start + MINLEN
            [sru, if_cov, median_diff, read_density, _, split, first_diff
             ] = extract_features(start, stop,
                                  profile_dict[transcript]["riboseq"])
            if read_density > MIN_READ_DENSITY:
                outfile.write("{},{},{},{},ATG,1,{},{},{},{},{},{}\n".format(
                    gene, transcript, start, stop, sru, if_cov, median_diff,
                    read_density, split, first_diff))
                acceptable_trans.append(transcript)
                acceptable_cases += 1
                tot_positive += 1

    # Negative cases, wrong frame (min highrame)
    for transcript in acceptable_trans:
        if not traninfo_dict[transcript]["cds_start"]:
            continue
        acceptable_cases = 0
        attempts = 0
        while acceptable_cases < (
                CASES_PER_TRAN /
                2) and attempts < MAX_ATTEMPTS and tot_negative < tot_positive:
            attempts += 1
            gene = traninfo_dict[transcript]["gene"]
            cds_start = traninfo_dict[transcript]["cds_start"]
            cds_stop = traninfo_dict[transcript]["cds_stop"]
            if cds_stop - cds_start < ((MINLEN * 2) + 100):
                continue
            start = random.randint(cds_start + 100, cds_stop)
            while start % 3 == cds_start % 3:
                start = random.randint(cds_start + 100, cds_stop)
            stop = start + MINLEN
            [sru, if_cov, median_diff, read_density, _, split, first_diff
             ] = extract_features(start, stop,
                                  profile_dict[transcript]["riboseq"])
            if read_density > MIN_READ_DENSITY:
                outfile.write("{},{},{},{},NULL,0,{},{},{},{},{},{}\n".format(
                    gene, transcript, start, stop, sru, if_cov, median_diff,
                    read_density, split, first_diff))
                acceptable_cases += 1
                tot_negative += 1

    # negative cases 3'trailer
    for transcript in acceptable_trans:
        if not traninfo_dict[transcript]["cds_start"]:
            continue
        acceptable_cases = 0
        attempts = 0
        while (acceptable_cases < CASES_PER_TRAN and attempts < MAX_ATTEMPTS
               and tot_negative < tot_positive):
            attempts += 1
            gene = traninfo_dict[transcript]["gene"]
            try:
                cds_start = random.randint(
                    traninfo_dict[transcript]["cds_stop"] + 12,
                    traninfo_dict[transcript]["length"] - MINLEN)
            except Exception:
                continue
            cds_stop = cds_start + MINLEN
            # Try 20 times to find a position with at least 2 inframe reads

            cds_start = random.randint(
                traninfo_dict[transcript]["cds_stop"] + 12,
                traninfo_dict[transcript]["length"] - MINLEN)
            cds_stop = cds_start + MINLEN
            [sru, if_cov, median_diff, read_density, _, split, first_diff
             ] = extract_features(cds_start, cds_stop,
                                  profile_dict[transcript]["riboseq"])
            outfile.write("{},{},{},{},NULL,0,{},{},{},{},{},{}\n".format(
                gene, transcript, cds_start, cds_stop, sru, if_cov,
                median_diff, read_density, split, first_diff))
            acceptable_cases += 1
            tot_negative += 1
    outfile.close()


def create_test_set(
        profile_dict,
        orf_dict,
        traninfo_dict,
        short_code,  # MINLEN,
        sequence_dict,
        region):
    filename = short_code + ".test.csv"
    outfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename), "w")
    # outfile = open("test.csv","w")
    outfile.write(HEADER)
    # Positive cases
    for transcript in orf_dict:
        try:
            gene = traninfo_dict[transcript]["gene"]
        except Exception:
            gene = transcript
        # if transcript != "ENST00000262030":
        # continue
        for stop in orf_dict[transcript]:
            # if stop != 1239:
            # continue
            for start in orf_dict[transcript][stop]:
                [sru, if_cov, median_diff, read_density, _, split, first_diff
                 ] = extract_features(start, stop,
                                      profile_dict[transcript]["riboseq"])
                if read_density > MIN_READ_DENSITY:
                    try:
                        seq = sequence_dict[transcript][start - 1:stop]
                        start_codon = seq[:3]
                    except Exception:
                        seq = ""
                        start_codon = None

                    outfile.write(
                        "{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
                            gene, transcript, start, stop, start_codon, region,
                            sru, if_cov, median_diff, read_density, split,
                            first_diff))
    outfile.close()


def create_model(columns):
    # create model
    model = Sequential()
    model.add(
        BatchNormalization(axis=-1,
                           momentum=0.99,
                           epsilon=0.001,
                           center=True,
                           scale=True))
    model.add(Dropout(0.33, input_shape=(columns - 1, )))
    # model.add(Dropout(0.33))
    # model.add(BatchNormalization(axis=-1,momentum=0.99,epsilon=0.001,center=True,scale=True))
    model.add(Dense(2, kernel_initializer='uniform', activation='relu'))
    # model.add(Dropout(0.5))
    # model.add(Dense(8, kernel_initializer='uniform', activation='relu'))
    # model.add(Dropout(0.5))
    # model.add(BatchNormalization(axis=-1,momentum=0.99,epsilon=0.001,center=True,scale=True))
    # model.add(Dense(10, kernel_initializer='uniform', activation='relu'))
    # model.add(Dropout(0.1))
    # model.add(Dense(10, kernel_initializer='uniform', activation='relu'))
    model.add(Dense(1, kernel_initializer='uniform', activation='sigmoid'))
    opt = optimizers.Adam(learning_rate=0.001)
    # Compile model
    model.compile(loss='binary_crossentropy',
                  optimizer=opt,
                  metrics=['accuracy', metrics.AUC()])
    return model


def neural_net(short_code, feature_list, organism, transcriptome, file_string):
    training_filename = short_code + ".training.csv"
    test_filename = short_code + ".test.csv"
    prediction_filename = short_code + ".predictions.csv"
    testfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,
                                              test_filename)).readlines()
    outfile = open(
        "{}/static/tmp/{}".format(config.SCRIPT_LOC, prediction_filename), "w")
    outfile.write("{},prediction\n".format(testfile[0].replace("\n", "")))
    returnstr = "Table|"
    columns = len(feature_list)
    model = create_model(columns)

    openfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,
                                              training_filename)).readlines()
    headers = openfile[0].replace(" ", "").replace("\n", "").split(",")

    index_list = []

    for feature in feature_list:
        if feature in headers:
            ind = headers.index(feature)
            index_list.append(ind)
        else:

            returnstr = "Could not find feature {} in training dataset: {}".format(
                feature, training_filename)
            return returnstr
    dataset = np.loadtxt("{}/static/tmp/{}".format(config.SCRIPT_LOC,
                                                   training_filename),
                         delimiter=",",
                         skiprows=1,
                         usecols=index_list)
    # split into input (X) and output (Y) variables
    X = dataset[:, 1:columns]
    Y = dataset[:, 0]
    randomize = np.arange(len(X))
    np.random.shuffle(randomize)
    X = X[randomize]
    Y = Y[randomize]
    # Fit the model
    model.fit(X,
              Y,
              nb_epoch=600,
              batch_size=50,
              validation_split=0.3,
              verbose=0)
    # evaluate the model
    _ = model.evaluate(X, Y)

    testset = np.loadtxt("{}/static/tmp/{}".format(config.SCRIPT_LOC,
                                                   test_filename),
                         delimiter=",",
                         skiprows=1,
                         usecols=index_list[1:])
    X = testset[:, 0:columns]

    predictions = model.predict(X)
    tup_count = 0
    prev_row = []
    for row in testfile[1:]:
        row = row.replace("\n", "")
        splitrow = row.split(",")
        gene = splitrow[0]
        tran = splitrow[1]
        start = splitrow[2]
        stop = splitrow[3]

        # continue
        prediction = predictions[tup_count][0]
        tup_count += 1
        splitrow.append(prediction)

        if prev_row == []:
            prev_row = splitrow
            continue
        if tran == prev_row[1] and stop == prev_row[3]:
            if prediction > prev_row[-1]:
                prev_row = splitrow
        else:
            gene = prev_row[0]
            transcript = prev_row[1]
            start = prev_row[2]
            stop = prev_row[3]
            orftype = prev_row[5]
            outfile.write("{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
                prev_row[0], prev_row[1], prev_row[2], prev_row[3],
                prev_row[4], prev_row[5], prev_row[6], prev_row[7],
                prev_row[8], prev_row[9], prev_row[10], prev_row[11],
                prev_row[12]))
            prev_row = splitrow
    outfile.close()

    # Open outfile again, read in all lines with probablity > 0.5 and sort
    result_list = []
    results = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,
                                             prediction_filename)).readlines()
    for line in results[1:]:
        splitline = line.split(",")
        gene = splitline[0]
        transcript = splitline[1]
        start = splitline[2]
        stop = splitline[3]
        orftype = splitline[5]
        length = int(stop) - int(start)
        if prediction > 0.5:
            trips_link = '<a href="https://trips.ucc.ie/' + organism + '/' + transcriptome + '/interactive_plot/?tran=' + transcript + '&hili=' + str(
                start
            ) + "_" + str(
                stop
            ) + '&files=' + file_string + '&rs=0.5" target="_blank_" >View on trips-viz</a>'
            tup = [
                gene, transcript, start, stop, length, prediction, orftype,
                trips_link
            ]
            result_list.append(tup)
    sorted_result_list = sorted(result_list, key=lambda tup: tup[5])
    for tup in sorted_result_list[-1000:]:
        returnstr += "{},{},{},{},{},{},{},NULL,NULL,NULL,NULL,{}.,/".format(
            tup[0], tup[1], tup[2], tup[3], tup[4], tup[5], tup[6], tup[7])

    return returnstr


def extract_values(accepted_orf_dict, data, tran_gene_dict, selected_seq_types,
                   profile_dict, all_cases):
    logging.debug("extract values called")
    best_high_frame = 0
    best_low_frame = 0
    best_start_score = 0
    best_stop_score = 0
    best_inframe_cov = 0
    all_values = []
    tot_loc = 0
    for locus in accepted_orf_dict:
        tot_loc += 1
        if tot_loc % 100 == 0:
            print("total transcripts {}".format(tot_loc))
        for stop in accepted_orf_dict[locus]:
            best_values = {
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
                "start_geo": 0
            }
            for start in accepted_orf_dict[locus][stop]:
                transcript = accepted_orf_dict[locus][stop][start][
                    "transcript"]
                best_values["transcript"] = transcript
                if best_values["start"] == -1:
                    best_values["start"] = start
                    best_values["stop"] = stop
                gene = tran_gene_dict[transcript]
                if "riboseq" in selected_seq_types and "proteomics" in selected_seq_types:
                    profile = profile_dict[transcript]["riboseq"]
                    prot_profile = profile_dict[transcript]["proteomics"]
                elif "proteomics" not in selected_seq_types:
                    profile = profile_dict[transcript]["riboseq"]
                    prot_profile = {}
                else:
                    profile = profile_dict[transcript]["proteomics"]
                    prot_profile = {}
                transcriptome_stop = accepted_orf_dict[locus][stop][start][
                    "stop"]
                orftype = accepted_orf_dict[locus][stop][start]["orftype"]
                length = (transcriptome_stop - start) + 1

                if_cov = 0.0
                if_len = (stop - start) / 3
                for p in range(start + 1, transcriptome_stop, 3):
                    if p - 1 in prot_profile:
                        proteomics_count += prot_profile[p - 1]

                if_cov = 0.0
                if_cov_total = 0.0
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
                    if max(inframe_values[x], minusone_values[x],
                           plusone_values[x]) > 0:
                        if_cov_total += 1
                        if "highest_frame_diff_check" in data:
                            if inframe_values[x] > max(minusone_values[x],
                                                       plusone_values[x]):
                                if_cov += 1
                        else:
                            if inframe_values[x] > min(minusone_values[x],
                                                       plusone_values[x]):
                                if_cov += 1

                if if_cov_total > 0:
                    if_cov = if_cov / if_cov_total
                    if_cov = if_cov * 100
                else:
                    if_cov = 0
                start_score_raw = (float(sum(inframe_values[4:8])) +
                                   1) / (float(sum(inframe_values[:4])) + 1)
                start_score = np.log(start_score_raw)
                stop_score_raw = (float(sum(inframe_values[-8:-4])) +
                                  1) / (float(sum(inframe_values[-4:])) + 1)
                stop_score = np.log(stop_score_raw)
                inframe_sum = float(sum(inframe_values[4:-4]))
                highframe_sum = float(
                    max(sum(minusone_values[4:-4]), sum(plusone_values[4:-4])))
                lowframe_sum = float(
                    min(sum(minusone_values[4:-4]), sum(plusone_values[4:-4])))
                non_zero_counts = 0.0
                for i in inframe_values[4:-4]:
                    if i != 0:
                        non_zero_counts += 1
                read_density = non_zero_counts / (float(length) / 3)

                high_frame_count = inframe_sum - highframe_sum
                lowest_frame_count = inframe_sum - lowframe_sum

                # Ratios
                start_ratio = float(sum(inframe_values[4:8])) / float(
                    (sum(inframe_values[:4]) + 1))
                stop_ratio = float(sum(inframe_values[-8:-4])) / float(
                    (sum(inframe_values[-4:]) + 1))
                high_ratio = inframe_sum / (highframe_sum + 1)
                low_ratio = inframe_sum / (lowframe_sum + 1)
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
                            and start < best_values["start"]):
                        best_values["start"] = start
                        best_values["stop"] = transcriptome_stop
                        best_values["transcript"] = transcript
                        best_values["length"] = length
                        best_values["high_frame_count"] = high_frame_count
                        best_values["low_frame_count"] = lowest_frame_count
                        best_values["start_score"] = start_score
                        best_values["stop_score"] = stop_score
                        best_values["final_score"] = final_score
                        best_values["coverage"] = if_cov
                        best_values["proteomics_count"] = proteomics_count
                        best_values["read_density"] = read_density
                        best_values["inframe_count"] = float(
                            inframe_sum) / float(if_len)
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
                    all_values.append([
                        gene, transcript, start, transcriptome_stop, length,
                        high_frame_count, lowest_frame_count, stop_score,
                        start_score, if_cov,
                        float(inframe_sum) / float(if_len), 0, 0, 0, 0, 0, 0,
                        0, orftype, proteomics_count, read_density,
                        start_ratio, stop_ratio, high_ratio, low_ratio
                    ])
            if not all_cases:
                all_values.append([
                    gene, best_values["transcript"], best_values["start"],
                    best_values["stop"], best_values["length"],
                    best_values["high_frame_count"],
                    best_values["low_frame_count"], best_values["stop_score"],
                    best_values["start_score"], best_values["coverage"],
                    best_values["inframe_count"], 0, 0, 0, 0, 0, 0, 0, orftype,
                    best_values["proteomics_count"],
                    best_values["read_density"], best_values["start_ratio"],
                    best_values["stop_ratio"], best_values["high_ratio"],
                    best_values["low_ratio"]
                ])
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
    logging.debug("LENGTH OF ALL SORTED ROWS {}".format(
        len(sorted_all_values)))
    final_rank = 1
    for tup in sorted_all_values:
        tup[17] = final_rank
        final_rank += 1
    return sorted_all_values


def write_to_file(sorted_all_values, filename, sequence_dict, organism,
                  transcriptome, file_string):
    returnstr = "Table|"
    tmp_result_file = open(
        "{}/static/tmp/{}".format(config.SCRIPT_LOC, filename), "w")
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
            seq = sequence_dict[transcript][start - 1:stop]
        except Exception:
            seq = ""
        start_codon = seq[:3]
        # Get amino acid sequence of this ORF, excluding the stop codon.
        while len(seq) % 3 != 0:
            seq = seq[:-1]

        aa_seq = nuc_to_aa(seq[:-3])
        trips_link = '<a href="https://trips.ucc.ie/' + organism + '/' + transcriptome + '/interactive_plot/?tran=' + transcript + '&hili=' + str(
            start
        ) + "_" + str(
            stop
        ) + '&files=' + file_string + '" target="_blank_" >View on trips-viz</a>'
        tmp_result_file.write(
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n"
            .format(gene, transcript, start, stop, length, global_rank,
                    orftype, trips_link, start_codon, tup[11], tup[5], tup[12],
                    tup[6], tup[13], tup[7], tup[14], tup[8], tup[15], tup[9],
                    tup[16], tup[10], aa_seq, proteomics_count, read_density))
        if tup_count < 1000:
            returnstr += "{},{},{},{},{},{},{},NULL,NULL,NULL,NULL,{}.,/".format(
                gene, transcript, start, stop, length, global_rank, orftype,
                trips_link)
        tup_count += 1
    return returnstr


def find_orfs(data, user, logged_in):
    logging.debug("orfquery called")
    owner = get_table("organisms").filter(
        (pl.col("organism_name") == data["organism"])
        & (pl.col("transcriptome_list") == data["transcriptome"]))[0, "owner"]
    # TODO: Check transcriptome list option

    start_time = time.time()

    prog_count = 0
    if "nnet_check" in data:
        data["output"] = "nnet"
    else:
        pass

    total_files = len(data["file_list"])
    if total_files > 100:
        return "A maximum of 100 files can be selected on this page, currently there are {} selected".format(
            total_files)

    file_paths_dict = fetch_file_paths(data["file_list"], organism)
    # Find out which studies have all files of a specific sequence type selected (to create aggregates)
    files = get_table("files")

    aggregate_dict = {}
    full_studies = []
    if data["minscore"] == 0:
        for seq_type in file_paths_dict:
            all_file_ids = list(file_paths_dict[seq_type].keys())
            all_study_ids = files.loc[(files["file_id"].isin(all_file_ids)
                                       & files["file_type"] == seq_type),
                                      ["study_id", "file_id"]].unique()
            prog_count += 100.
            # NOTE: Till here

            for study_id in all_study_ids:
                cursor.execute(
                    "SELECT file_id FROM files WHERE study_id = {} AND file_type = '{}'"
                    .format(study_id, seq_type))
                result = cursor.fetchall()
                study_file_ids = []
                for row in result:
                    study_file_ids.append(row[0])
                logging.debug("study id {}".format(study_id))
                # Do not create an aggregate if there is only one file in the study:
                if len(study_file_ids) > 1:
                    logging.debug("length over 1")
                    # Check if all the file_ids in this study are present in all_file_ids, if they are this study is "full" and an aggreagate will be made
                    if set(study_file_ids).issubset(all_file_ids):
                        # This is a full study, create an aggregate for it if non exits
                        # Get path to study directory
                        study_path = "/".join(
                            (file_paths_dict[seq_type][study_file_ids[0]]
                             ).split("/")[:-1])
                        logging.debug("study_path {}".format(study_path))
                        full_studies.append(study_id)
                        # Check if aggregate file exists
                        agg_filename = "{}/aggregate_0.5_{}.sqlite".format(
                            study_path, seq_type)
                        if not os.path.isfile(agg_filename):
                            sub_file_paths_dict = {seq_type: {}}
                            for file_id in study_file_ids:
                                sub_file_paths_dict[seq_type][
                                    file_id] = file_paths_dict[seq_type][
                                        file_id]
                            # create_aggregate(sub_file_paths_dict,study_path, seq_type)
                            aggregate_dict[study_path] = [
                                sub_file_paths_dict, seq_type
                            ]
                            # logging.debug("aggregate created")
                        # Remove all file_ids associated with this study from file_paths_dict
                        for file_id in study_file_ids:
                            del file_paths_dict[seq_type][file_id]
                        # Add the aggregate to file_paths_dict
                        file_paths_dict[seq_type]["STUDY_" + str(
                            study_id)] = "{}/aggregate_0.5_{}.sqlite".format(
                                study_path, seq_type)
    cursor.close()
    connection.close()

    for study_path in aggregate_dict:
        create_aggregate(aggregate_dict[study_path][0], study_path,
                         aggregate_dict[study_path][1])

    logging.debug("Full studies {}".format(full_studies))
    # return str(full_studies)
    # return str(all_study_ids)

    html_args = data["html_args"]
    returnstr = ""

    start_codons = []

    # Used to extract columns from csv for neural neural_net
    feature_list = ["Type"]
    for s_codon in ["aug", "cug", "gug", "none"]:
        if 'sc' + s_codon in data:
            start_codons.append(s_codon.upper())
    logging.debug("start codons {}".format(start_codons))

    all_cases = False
    if data["max_len"] == "all_cases":
        all_cases = True
        max_len = 10000
    else:
        max_len = float(data["max_len"])

    min_avg = float(data["min_avg"])
    max_avg = float(data["max_avg"])

    accepted_orftypes = []

    region = data["region"]
    accepted_orftypes.append(region)

    logging.debug("accepted orftypes {}".format(accepted_orftypes))

    try:
        cons_score = data["cons_score"].strip()
    except Exception:
        cons_score = ""
    user_defined_transcripts = []
    # tran_list is a radio button, user can choose between principal, all or a custom list
    # custom_tran_list is the actual comma seperated list of transcripts that user would enter should they choose the custom option in tranlist

    ambig = True if "ambig_check" in data else False

    filtered_transcripts = {}

    if "saved_check" in data:
        if current_user.is_authenticated:
            user_id = get_user_id(current_user.name)
            user_saved_cases = get_table('users_saved_cases').filter(
                (pl.col('user_id') == user_id)
                & (pl.col('organism') == organism))[["tran", "stop"]]
            "SELECT tran,stop FROM users_saved_cases WHERE user_id = '{}' and organism = '{}';"

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
    if html_args["user_short"] == "None":
        short_code = generate_short_code(data, organism, data["transcriptome"],
                                         "orf_translation")
    else:
        short_code = html_args["user_short"]
        user_short_passed = True

    filename = short_code + ".csv"
    if os.path.isfile("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename)):
        logging.debug("File exists {}/static/tmp/{}".format(
            config.SCRIPT_LOC, filename))
        tmp_result_file = open(
            "{}/static/tmp/{}".format(config.SCRIPT_LOC, filename), "r")
        total_rows = 0
        returnstr = "Table|"
        for row in tmp_result_file:
            splitrow = row.split(",")
            if len(splitrow) <= 1:
                continue
            tran = splitrow[1].strip(" ")
            try:
                stop = int(splitrow[3].strip(" "))
            except Exception:
                continue
            if tran in filtered_transcripts:
                if stop in filtered_transcripts[tran]:
                    continue

            if total_rows < 1001:
                returnstr += "{},{},{},{},{},{},{},NULL,NULL,NULL,NULL,{}.,/".format(
                    splitrow[0], splitrow[1], splitrow[2], splitrow[3],
                    splitrow[4], splitrow[5], splitrow[6], splitrow[7])
                total_rows += 1
            else:
                break
        returnstr += "|"
        returnstr += "Min,NULL,NULL,NULL,NULL,NULL,NULL.,/"
        returnstr += "10th_percentile,NULL,NULL,NULL,NULL,NULL,NULL.,/"
        returnstr += "25th_percentile,NULL,NULL,NULL,NULL,NULL,NULL.,/"
        returnstr += "50th_percentile,NULL,NULL,NULL,NULL,NULL,NULL.,/"
        returnstr += "Max,NULL,NULL,NULL,NULL,NULL,NULL.,/"
        returnstr += "|{}".format(filename)
        returnstr += "|{}".format(
            str(data["file_list"]).strip("[]").replace("'",
                                                       "").replace(" ", ""))
        returnstr += "|{}".format(data["html_args"]["user_short"])
        return returnstr
    else:
        logging.debug("File does not exists {}/static/tmp/{}".format(
            config.SCRIPT_LOC, filename))

    if data["custom_tran_list"]:
        data["custom_tran_list"] = data["custom_tran_list"].split(",")
        for item in data["custom_tran_list"]:
            user_defined_transcripts.append(item)

    if not accepted_orftypes:
        return "Error no ORF type selected"
    tran_gene_dict = {}
    # structure of orf dict is transcript[stop][start] = {"length":x,"score":0,"cds_cov":0} each stop can have multiple starts
    accepted_orf_dict = {}

    if not start_codons:
        return "Error no start codon types selected"

    if owner:
        sqlite_path = "{0}/{1}/{2}/{2}.{3}.sqlite".format(
            config.SCRIPT_LOC, config.ANNOTATION_DIR, organism, transcriptome)
        if not os.path.isfile(sqlite_path):
            return "Cannot find annotation file {}.{}.sqlite".format(
                organism, transcriptome)
    else:
        sqlite_path = "{0}/transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
            config.UPLOADS_DIR, owner, organism, transcriptome)
    transcripts = sqlquery(sqlite_path, "transcripts")
    principle_transcripts = transcripts.filter(
        pl.col("transcript").is_in(principal_transcripts))

    traninfo_dict = {}
    if data["output"] == "nnet":
        pass

    principal_transcripts = []
    if data["tran_list"] == "prin_trans":
        transcripts = transcripts.filter(pl.col("principal") == 1)
    elif data["tran_list"] == "all_trans":
        pass
    elif data["tran_list"] == "custom_trans":
        transcripts = transcripts.filter(
            pl.col("transcript").is_in(data["custom_tran_list"]))

    exons = sqlquery(sqlite_path, "exons").filter(
        pl.col("transcript").is_in(principal_transcripts)).join(
            principle_transcripts, on="transcript")

    logging.debug("building transcriptom info dict")
    # Holds a list of all transcripts in accepted_orf_dict
    accepted_transcript_list = []
    for table_name in accepted_orftypes:
        # logging.debug("table_name", table_name)
        # logging.debug("start codons", start_codons)
        table = sqlquery(sqlite_path, table_name).filter(
            (pl.col("cds_coverage") >= data["min_cds"]) & (pl.col("cds_coverage") <= data["max_cds"]) & (pl.col("length") >= data["min_len"]) & (pl.col("length") <= data["max_len"]) & (pl.col("transcript"").is_in(principle_transcript))
        )
        if "any" in start_codons:
            logging.debug("selecting any start")
            pass
        else:
            logging.debug("selecting aug,cug,gug")
            table = table.filter(pl.col("start_codon").is_in(start_codons))
        logging.debug("for row in result")
        rows = 0
        for row in result:
            rows += 1
            if rows % 1000 == 0:
                logging.debug("Rows: {}".format(rows))
                prog_count = 100 + (100 * (rows / total_trans))
            # logging.debug("row", row)
            transcript = str(row[0])
            start_codon = str(row[1])
            length = row[2]
            cds_cov = row[3]
            start = row[4]
            # Orfs with multiple potential starts will be grouped by start codon (so only one potential start is reported)
            # If the user has selected all transcripts then instead group by the genomic stop codon co-ordinates (so only one ORF will be reported,
            # even if it occurs on multiple transcript isoforms)
            if data["tran_list"] != "all_trans2":
                stop = row[5]
                locus = transcript
            else:
                stop = tran_to_genome(transcript, row[5],
                                      transcriptome_info_dict)
                locus = tran_gene_dict[transcript]
            if transcript in filtered_transcripts:
                if stop in filtered_transcripts[transcript]:
                    continue
            transcriptome_stop = row[5]
            if locus not in accepted_orf_dict:
                accepted_orf_dict[locus] = {}
            if stop not in accepted_orf_dict[locus]:
                accepted_orf_dict[locus][stop] = {}
            if transcript not in accepted_transcript_list:
                accepted_transcript_list.append(transcript)
            accepted_orf_dict[locus][stop][start] = {
                "length": length,
                "score": 0,
                "cds_cov": cds_cov,
                "start_codon": start_codon,
                "orftype": table_name,
                "stop": transcriptome_stop,
                "transcript": transcript
            }

    # logging.debug("accepted orf dict", accepted_orf_dict)
    logging.debug("accepted orf dict built")
    # Now build a profile for every transcript in accepted_transcripts

    if file_paths_dict["rnaseq"] == {}:
        if "te_check" in data:
            del data["te_check"]

    if not file_paths_dict["riboseq"] and not file_paths_dict[
            "proteomics"]:
        return "Error no files selected"

    total_files = 0
    selected_seq_types = []
    if "riboseq" in file_paths_dict:
        total_files += len(file_paths_dict["riboseq"])
        if "riboseq" not in selected_seq_types:
            selected_seq_types.append("riboseq")
    if "proteomics" in file_paths_dict:
        total_files += len(file_paths_dict["proteomics"])
        if "proteomics" not in selected_seq_types:
            selected_seq_types.append("proteomics")

    total_trans = len(accepted_transcript_list)
    if data["output"] == "nnet":
        for transcript in traninfo_dict:
            if transcript not in accepted_transcript_list:
                accepted_transcript_list.append(transcript)
    logging.debug("total trans".format(total_trans))

    profile_dict, file_string = create_profiles(file_paths_dict,
                                                accepted_transcript_list,
                                                ambig, total_files,
                                                data["minscore"])

    logging.debug("profile dict built\nextracting values")
    if data["output"] != "nnet":
        sorted_all_values = extract_values(accepted_orf_dict, data,
                                           tran_gene_dict, selected_seq_types,
                                           profile_dict, all_cases)
        if not sorted_all_values:
            return "No results, try making filters less restrictive"
    else:
        create_training_set(profile_dict, traninfo_dict, short_code,
                            data["min_len"])
        create_test_set(profile_dict, accepted_orf_dict, traninfo_dict,
                        short_code, data["min_len"], sequence_dict, region)
        returnstr = neural_net(
            short_code, ["type", "coverage", "median_diff", "first_diff"],
            organism, transcriptome, file_string)

    if data["output"] != "nnet":
        logging.debug("Writing to file")
    returnstr = write_to_file(sorted_all_values, filename, sequence_dict,
                              organism, transcriptome, file_string)

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
                msg['From'] = fromaddr
                msg['To'] = toaddr
                msg['Subject'] = "Trips-Viz job completion"
                msg.attach(
                    MIMEText(
                        "Your Trips-Viz job is complete: https:trips.ucc.ie/short/{}"
                        .format(short_code)))
                server = smtplib.SMTP('smtp.gmail.com', 587)
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


# Returns a table with ranked orf scores
orfquery_blueprint = Blueprint("orfquery",
                               __name__,
                               template_folder="templates")


@orfquery_blueprint.route('/orfquery', methods=['POST'])
def orfquery():
    data = request.args.to_dict()
    user, logged_in = fetch_user()
    return find_orfs(data, user, logged_in)
