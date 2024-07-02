from typing import Dict, List, Tuple, Union
from flask import Blueprint, render_template, request
from sqlitedict import SqliteDict
from sqlqueries_2 import sqlquery, get_table
import os
import logging
import config
import subprocess
from core_functions import (fetch_studies, fetch_files, fetch_study_info,
                            fetch_file_paths, generate_short_code, form_filler,
                            build_profile, build_proteomics_profile,
                            fetch_user)
from fixed_values import my_decoder
import polars as pl

# This page is used to detect pauses
pause_detection_blueprint = Blueprint("pause_detection_page",
                                      __name__,
                                      template_folder="templates")


@pause_detection_blueprint.route(
    '/<organism>/<transcriptome>/pause_detection/')
def pause_detection_page(organism: str, transcriptome: str) -> str:
    data = form_filler(organism, transcriptome)
    accepted_studies = fetch_studies(data["gwips_info"][0, "organism_id"])
    data['files'] = fetch_files(accepted_studies).to_pandas()
    return render_template('pause_detection.html', template_dict=data)


def create_profiles(file_paths_dict, accepted_transcript_list, total_files,
                    min_read_length, max_read_length):
    ambig = False
    file_count = 0
    # This will be populated with the users chosen file_ids and passed to the table, so that the trips link can use these files aswell.
    file_string = ""
    label_string = "&labels="
    seq_types = ["riboseq", "proteomics"]
    color_list = [
        "#ff0000", "#2bff00", "#0004ff", "#ffa200", "#c800ff",
        "#000000", "#969696", "#fa00f2"
    ]
    color_ind = 0
    profile_dict = {}

    tot_file_ids = 0.0
    file_count = 0
    for seq_type in seq_types:
        tot_file_ids += len(file_paths_dict[seq_type])

    for seq_type in seq_types:
        if seq_type not in file_paths_dict:
            continue
        for file_id in file_paths_dict[seq_type]:
            profile_dict[file_id] = {}
            file_count += 1
            file_name = (
                file_paths_dict[seq_type][file_id].split("/")[-1]).replace(
                    ".sqlite", "").replace("_", " ")
            sqlite_db = SqliteDict(f"{file_paths_dict[seq_type][file_id]}",
                                   autocommit=False,
                                   decode=my_decoder)
            file_string += "{};{}_".format(file_id, color_list[color_ind])
            label_string += "{};{}_".format(file_name, color_list[color_ind])
            color_ind += 1
            offsets = {}
            scores = {}
            if seq_type == "riboseq":
                try:
                    offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
                except Exception:
                    pass
                try:
                    scores = sqlite_db["offsets"]["fiveprime"]["read_scores"]
                except Exception:
                    pass
            for transcript in accepted_transcript_list:
                if transcript not in profile_dict[file_id]:
                    profile_dict[file_id][transcript] = {
                        "riboseq": {},
                        "proteomics": {}
                    }
                try:
                    counts = sqlite_db[transcript]
                except Exception:
                    continue
                if seq_type == "riboseq":
                    subprofile = build_profile(counts,
                                               offsets,
                                               ambig,
                                               minscore=None,
                                               scores=scores)
                elif seq_type == "proteomics":
                    subprofile = build_proteomics_profile(counts, ambig)
                for pos in subprofile:
                    try:
                        profile_dict[file_id][transcript][seq_type][
                            pos] += subprofile[pos]
                    except Exception:
                        profile_dict[file_id][transcript][seq_type][
                            pos] = subprofile[pos]
    file_string = file_string[:-1]
    label_string = label_string[:-1]
    return (profile_dict, file_string, label_string)


def extract_values(traninfo_dict, data, tran_gene_dict, selected_seq_types,
                   profile_dict, min_fold_change, window, min_coverage,
                   nuc_output):
    step = window // 2
    all_values_dict = {}
    file_output_dict = {}
    file_number = len(list(profile_dict))
    for file_id in profile_dict:
        print("file id", file_id)
        for tran in profile_dict[file_id]:
            if tran not in all_values_dict:
                all_values_dict[tran] = {}

            score_dict = {}
            cov_dict = {}
            count_dict = {}
            region_dict = {}
            tranlen = int(traninfo_dict[tran]["length"])
            seq = traninfo_dict[tran]["seq"].replace("T", "U")
            gene = traninfo_dict[tran]["gene"]
            cds_start = traninfo_dict[tran]["cds_start"]
            cds_stop = traninfo_dict[tran]["cds_stop"]
            profile = profile_dict[file_id][tran]["riboseq"]
            region = "non-coding"
            for i in range(0, tranlen, step):
                count = 0.0
                cov_count = 0.0
                for x in range(i, i + window + 15):
                    if x in profile:
                        count += profile[x]
                        cov_count += 1
                cov = cov_count / window
                avg = count / window
                if cov >= min_coverage and count > 0:
                    for x in range(i, i + window):
                        if x in profile:
                            pos_count = profile[x]
                            score = pos_count / avg
                            if score > min_fold_change:
                                if cds_start:
                                    if x < cds_start - 3:
                                        region = "5' Leader"
                                    elif x >= cds_start - 3 and x <= cds_start + 3:
                                        region = "CDS Start"
                                    elif x > cds_start + 3 and x < cds_stop - 3:
                                        region = "CDS"
                                    elif x >= cds_stop - 3 and x <= cds_stop + 3:
                                        region = "CDS Stop"
                                    elif x > cds_stop + 3:
                                        region = "3' Trailer"
                                if x not in score_dict:
                                    score_dict[x] = []
                                    cov_dict[x] = []
                                    count_dict[x] = 0
                                    region_dict[x] = region
                                score_dict[x].append(score)
                                cov_dict[x].append(cov)
                                count_dict[x] = pos_count
            for pos in score_dict:
                all_over_min = True
                for score in score_dict[pos]:
                    if score < min_fold_change:
                        all_over_min = False
                if all_over_min:
                    avg_score = sum(score_dict[pos]) / len(score_dict[pos])
                    avg_cov = sum(cov_dict[pos]) / len(cov_dict[pos])
                    count = count_dict[pos]
                    region = region_dict[pos]
                    if pos not in all_values_dict[tran]:
                        all_values_dict[tran][pos] = {}
                    if file_id not in all_values_dict[tran][pos]:
                        all_values_dict[tran][pos][file_id] = {}
                    all_values_dict[tran][pos][file_id] = [
                        gene, tran, pos, avg_score, seq[pos - nuc_output:pos],
                        seq[pos:pos + nuc_output], avg_cov, count, region
                    ]
    all_values_list = []
    for tran in all_values_dict:
        for pos in all_values_dict[tran]:
            if len(all_values_dict[tran][pos]) == file_number:
                tot_score = 0.0
                tot_cov = 0.0
                tot_count = 0.0
                for file_id in all_values_dict[tran][pos]:
                    if file_id not in file_output_dict:
                        file_output_dict[file_id] = []

                    gene = all_values_dict[tran][pos][file_id][0]
                    tran = all_values_dict[tran][pos][file_id][1]
                    tot_score += all_values_dict[tran][pos][file_id][3]
                    useq = all_values_dict[tran][pos][file_id][4]
                    dseq = all_values_dict[tran][pos][file_id][5]
                    tot_cov += all_values_dict[tran][pos][file_id][6]
                    tot_count += all_values_dict[tran][pos][file_id][7]
                    region = all_values_dict[tran][pos][file_id][8]
                    file_output_dict[file_id].append([
                        gene, tran, pos,
                        all_values_dict[tran][pos][file_id][3], useq, dseq,
                        all_values_dict[tran][pos][file_id][6],
                        all_values_dict[tran][pos][file_id][7], region
                    ])
                avg_score = tot_score / file_number
                avg_cov = tot_cov / file_number
                avg_count = tot_count / file_number
                all_values_list.append([
                    gene, tran, pos, avg_score, useq, dseq, avg_cov, avg_count,
                    region
                ])

    sorted_all_values = sorted(all_values_list,
                               key=lambda x: x[3],
                               reverse=True)
    return (sorted_all_values, file_output_dict)


def write_to_file(sorted_all_values, file_output_dict, sequence_dict, organism,
                  transcriptome, file_string, label_string, short_code):
    # TODO: Write only when the number of results are more than 1000
    # logging.debug("all sorted all values", sorted_all_values)
    print("writing to file")
    returnstr = "Table|"
    tmp_filepath = "{}/static/tmp/{}.csv".format(config.SCRIPT_LOC, short_code)
    all_filepaths = tmp_filepath

    for file_id in file_output_dict:
        file_name = get_table("files").filter(pl.col("file_id") == file_id)[0,
        logging.debug(file_name)
        filepath= "{}/static/tmp/{}_pauses.csv".format(
            config.SCRIPT_LOC, file_name)
        outfile= open(filepath, "w")
        all_filepaths += " {}".format(filepath)
        for line in file_output_dict[file_id]:
            outfile.write("{},{},{},{},{},{},{},{},{}\n".format(
                line[0], line[1], line[2], line[3], line[4], line[5], line[6],
                line[7], line[8]))
        outfile.close()

    tmp_result_file= open(tmp_filepath, "w")
    print("tmp filepath", tmp_filepath)
    tmp_result_file.write(
        "Gene,Tran,Position,Region, Coverage,Pause Score,Upstream_sequence, Downstream_sequence,Count,Link\n"
    )
    tup_count= 0

    # logging.debug("writing to file",len(sorted_all_values))
    for tup in sorted_all_values:
        # logging.debug("tup", tup)
        gene= tup[0]
        transcript= tup[1]
        position= tup[2]
        pause_score= round(tup[3], 2)
        upstream_seq= tup[4]
        downstream_seq= tup[5]
        cov= tup[6]
        count= round(tup[7], 2)
        region= tup[8]

        comparison_url= "/{}/{}/comparison/?files={}{}&transcript={}&normalize=F&cov=T&ambig=F&minread=25&maxread=150&hili_start={}&hili_stop={}".format(
            organism, transcriptome, file_string, label_string, transcript,
            position - 15, position + 15)
        ebc_link= '<a href="{}" target="_blank_" >View</a>'.format(
            comparison_url)

        tmp_result_file.write("{},{},{},{},{},{},{},{},{},{}\n".format(
            gene, transcript, position, region, cov, pause_score, upstream_seq,
            downstream_seq, count, ebc_link))
        if tup_count < 1000:
            returnstr += "{},{},{},{},{},{},{},{},{}.,/".format(
                gene, transcript, position, region, pause_score, upstream_seq,
                downstream_seq, count, ebc_link)
        tup_count += 1
    tmp_result_file.close()
    # Create a zip file of all output files
    print("zip -j {}/static/tmp/{}.zip {}".format(config.SCRIPT_LOC,
                                                  short_code, all_filepaths))
    subprocess.call("zip -j {}/static/tmp/{}.zip {}".format(
        config.SCRIPT_LOC, short_code, all_filepaths),
                    shell=True)
    return returnstr


def find_pauses(data, user, logged_in):
    logging.debug("pause query called")

    print("organism, transcriptome", organism, transcriptome)
    owner= get_table("organisms").filter(
        pl.col("organism_name") == data["organism"])[0, "organism_owner"]

    file_paths_dict= fetch_file_paths(data["file_list"], organism)
    # Find out which studies have all files of a specific sequence type selected (to create aggregates)

    full_studies= []

    logging.debug("Full studies {}".format(full_studies))

    min_coverage= data["min_coverage"] / 100.

    # feature_list.append("Inframe Count Value")
    if not html_args["user_short"]:
        short_code= generate_short_code(data)
    else:
        short_code= html_args["user_short"]
        user_short_passed= True

    if data['tranlist'] == "custom_trans":
        data['custom_tran_list']= data['custom_tran_list'].split(',')

    # structure of orf dict is transcript[stop][start] = {"length":x,"score":0,"cds_cov":0} each stop can have multiple starts

    if owner == 1:
        sqlfile= "{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC,
                                                      config.ANNOTATION_DIR,
                                                      data['organism'],
                                                      data['transcriptome'])
        if not os.path.isfile(sqlfile):
            return "Cannot find annotation file {}.{}.sqlite".format(
                data['organism'], data['transcriptome'])
    else:
        sqlfile= "{0}/transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
            config.UPLOADS_DIR, owner, data['organism'], data['transcriptome'])
    traninfo= sqlquery(sqlfile, "transcripts")
    tran_gene_dict= {}

    principal_transcripts= []
    if data['tranlist'] == "prin_trans":
        traninfo= traninfo.filter(pl.col("principal") == 1)
    elif data['tranlist'] == "custom_trans":
        traninfo= traninfo.filter(
            pl.col("transcript").is_in(data['custom_tran_list']))
    else:
        pass
    tran_gene= traninfo[["transcript", "gene"]]
    tran_gene.gene= tran_gene.gene.apply(lambda x: x.replace(",", "_"))

    transcriptome_info_dict= traninfo[["transcript", "strand", "chrom"]]
    exons= sqlquery(sqlfile, "exons").filter(
        pl.col("transcript").is_in(traninfo['transcript']))
    transcriptome_info_dict= transcriptome_info_dict.merge(exons,
                                                            on="transcript")

    # logging.debug("accepted orf dict", accepted_orf_dict)
    logging.debug("accepted orf dict built")
    # Now build a profile for every transcript in accepted_transcripts

    if (not file_paths_dict["rnaseq"]) and ("te_check" in data):
        del data["te_check"]

    if not file_paths_dict["riboseq"] and not file_paths_dict["proteomics"]:
        return "Error no files selected"

    total_files= 0
    selected_seq_types= []
    if "riboseq" in file_paths_dict:
        total_files += len(file_paths_dict["riboseq"])
        if "riboseq" not in selected_seq_types:
            selected_seq_types.append("riboseq")
    if "proteomics" in file_paths_dict:
        total_files += len(file_paths_dict["proteomics"])
        if "proteomics" not in selected_seq_types:
            selected_seq_types.append("proteomics")

    profile_dict, file_string, label_string= create_profiles(
        file_paths_dict, principal_transcripts, total_files, min_read_length,
        max_read_length)
    sorted_all_values, file_output_dict= extract_values(
        traninfo_dict, data, tran_gene_dict, selected_seq_types, profile_dict,
        min_fold_change, window, min_coverage, nuc_output)
    if sorted_all_values:
        return "No results, try making filters less restrictive"

    # TODO change extension to csv if only one file
    filename= short_code + ".zip"
    returnstr= write_to_file(sorted_all_values, file_output_dict,
                              sequence_dict, organism, transcriptome,
                              file_string, label_string, short_code)

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

    logging.debug("returning result")
    print("return returnstr")
    return returnstr


# Returns a table with ranked orf scores
pausequery_blueprint= Blueprint("pausequery",
                                 __name__,
                                 template_folder="templates")


@ pausequery_blueprint.route('/pausequery', methods=['POST'])
def pausequery():

    data= request.args.to_dict()
    user, logged_in= fetch_user()
    return find_pauses(data, user, logged_in)
