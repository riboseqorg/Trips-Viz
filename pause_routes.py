from typing import Dict, List, Tuple, Union
from flask import Blueprint, render_template, request
import sqlite3
from sqlitedict import SqliteDict
import os
import logging
import config
import subprocess
from core_functions import (fetch_studies, fetch_files, fetch_study_info,
                            fetch_file_paths, generate_short_code,
                            build_profile, build_proteomics_profile,
                            fetch_user, fetch_filename_file_id)
import json
from fixed_values import my_decoder

# This page is used to detect pauses
pause_detection_blueprint = Blueprint("pause_detection_page",
                                      __name__,
                                      template_folder="templates")


@pause_detection_blueprint.route('/<organism>/<transcriptome>/pause_detection/'
                                 )
def pause_detection_page(organism: str, transcriptome: str) -> str:
    # ip = request.environ['REMOTE_ADDR']
    organism = str(organism)
    user, logged_in = fetch_user()
    accepted_studies = fetch_studies(user, organism, transcriptome)
    file_id_to_name_dict, accepted_studies, accepted_files, seq_types = fetch_files(
        accepted_studies)
    advanced = False

    # holds all values the user could possibly pass in the url (keywords are after request.args.get), anything not passed by user will be a string: "None"
    html_args = request.args.to_dict()

    user_files = request.args.get('files')
    try:
        user_files = user_files.split(",")
        html_args["files"] = [str(x) for x in html_args["files"].split(",")]
    except KeyError:
        html_args["files"] = []

    # If ever need to use other seq types than these modify fetch_files to return a proper list
    seq_types = ["riboseq", "proteomics"]
    # If seq type not riboseq remove it from accepted files as orf translation is only appicable to riboseq
    del_types = []
    for seq_type in accepted_files:
        if seq_type not in seq_types:
            del_types.append(seq_type)
    for seq_type in del_types:
        del accepted_files[seq_type]
    studyinfo_dict = fetch_study_info(organism)

    print("html args", html_args)
    return render_template('pause_detection.html',
                           studies_dict=accepted_studies,
                           accepted_files=accepted_files,
                           user=user,
                           default_tran="",
                           advanced=advanced,
                           seq_types=seq_types,
                           studyinfo_dict=studyinfo_dict,
                           html_args=html_args,
                           organism=organism,
                           transcriptome=transcriptome)


def tran_to_genome(tran: str, pos: int, transcriptome_info_dict: Dict) -> str:
    if tran in transcriptome_info_dict:
        traninfo = transcriptome_info_dict[tran]
    else:
        return None, 0  # Fix same as another file
    exon_start = traninfo["exons"][0][0]
    if traninfo["strand"] == "-":
        traninfo["exons"] = traninfo["exons"][::-1]
        exon_start = traninfo["exons"][0][1]
    # logging.debug(exons)
    for tup in traninfo["exons"]:
        exonlen = tup[1] - tup[0]
        if pos > exonlen:
            pos = (pos - exonlen) - 1
        else:
            break
    if traninfo["strand"] == "+":
        genomic_pos = (exon_start + pos) - 1
    else:
        genomic_pos = (exon_start - pos) + 1
    return "{}_{}".format(traninfo["chrom"], genomic_pos)


def create_profiles(file_paths_dict, accepted_transcript_list, total_files,
                    min_read_length, max_read_length):
    ambig = False
    file_count = 0
    # This will be populated with the users chosen file_ids and passed to the table, so that the trips link can use these files aswell.
    file_string = ""
    label_string = "&labels="
    seq_types = ["riboseq", "proteomics"]
    color_list = [
        "%23ff0000", "%232bff00", "%230004ff", "%23ffa200", "%23c800ff",
        "%23000000", "%23969696", "%23fa00f2"
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
        sqlite_db.close()
    file_string = file_string[:-1]
    label_string = label_string[:-1]
    return (profile_dict, file_string, label_string)


def extract_values(traninfo_dict, data, tran_gene_dict, selected_seq_types,
                   profile_dict, min_fold_change, window, min_coverage,
                   nuc_output):
    step = int(window / 2)
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
    # logging.debug("all sorted all values", sorted_all_values)
    print("writing to file")
    returnstr = "Table|"
    tmp_filepath = "{}/static/tmp/{}.csv".format(config.SCRIPT_LOC, short_code)
    all_filepaths = tmp_filepath

    for file_id in file_output_dict:
        file_name = fetch_filename_file_id(file_id)
        logging.debug(file_name)
        filepath = "{}/static/tmp/{}_pauses.csv".format(
            config.SCRIPT_LOC, file_name)
        outfile = open(filepath, "w")
        all_filepaths += " {}".format(filepath)
        for line in file_output_dict[file_id]:
            outfile.write("{},{},{},{},{},{},{},{},{}\n".format(
                line[0], line[1], line[2], line[3], line[4], line[5], line[6],
                line[7], line[8]))
        outfile.close()

    tmp_result_file = open(tmp_filepath, "w")
    print("tmp filepath", tmp_filepath)
    tmp_result_file.write(
        "Gene,Tran,Position,Region, Coverage,Pause Score,Upstream_sequence, Downstream_sequence,Count,Link\n"
    )
    tup_count = 0

    # logging.debug("writing to file",len(sorted_all_values))
    for tup in sorted_all_values:
        # logging.debug("tup", tup)
        gene = tup[0]
        transcript = tup[1]
        position = tup[2]
        pause_score = round(tup[3], 2)
        upstream_seq = tup[4]
        downstream_seq = tup[5]
        cov = tup[6]
        count = round(tup[7], 2)
        region = tup[8]

        comparison_url = "/{}/{}/comparison/?files={}{}&transcript={}&normalize=F&cov=T&ambig=F&minread=25&maxread=150&hili_start={}&hili_stop={}".format(
            organism, transcriptome, file_string, label_string, transcript,
            position - 15, position + 15)
        ebc_link = '<a href="{}" target="_blank_" >View</a>'.format(
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
    global user_short_passed
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    cursor = connection.cursor()

    organism = data["organism"]
    transcriptome = data["transcriptome"]
    print("organism, transcriptome", organism, transcriptome)
    cursor.execute(
        "SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';"
        .format(organism, transcriptome))
    owner = (cursor.fetchone())[0]

    file_paths_dict = fetch_file_paths(data["file_list"], organism)
    # Find out which studies have all files of a specific sequence type selected (to create aggregates)

    aggregate_dict = {}
    full_studies = []

    logging.debug("Full studies {}".format(full_studies))
    # return str(full_studies)
    # return str(all_study_ids)

    html_args = data["html_args"]
    returnstr = ""

    min_fold_change = int(data["min_fold_change"])
    min_read_length = int(data["min_read_length"])
    max_read_length = int(data["max_read_length"])
    window = int(data["window_size"])
    min_coverage = float(data["min_coverage"]) / 100
    nuc_output = int(data["nuc_output"])

    user_defined_transcripts = []
    # tran_list is a radio button, user can choose between principal, all or a custom list
    tranlist = data["tran_list"]
    # custom_tran_list is the actual comma seperated list of transcripts that user would enter should they choose the custom option in tranlist
    custom_tran_list = data["custom_tran_list"]

    # feature_list.append("Inframe Count Value")
    if not html_args["user_short"]:
        short_code = generate_short_code(data, organism, data["transcriptome"],
                                         "pause_pred")
    else:
        short_code = html_args["user_short"]
        user_short_passed = True

    if custom_tran_list != "":
        custom_tran_list = custom_tran_list.replace(",", " ")
        for item in custom_tran_list.split(" "):
            user_defined_transcripts.append(item)

    tran_gene_dict = {}
    # structure of orf dict is transcript[stop][start] = {"length":x,"score":0,"cds_cov":0} each stop can have multiple starts

    if owner == 1:
        if os.path.isfile("{0}/{1}/{2}/{2}.{3}.sqlite".format(
                config.SCRIPT_LOC, config.ANNOTATION_DIR, organism,
                transcriptome)):
            traninfo_connection = sqlite3.connect(
                "{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC,
                                                    config.ANNOTATION_DIR,
                                                    organism, transcriptome))
        else:
            return prepare_return_str(
                "Cannot find annotation file {}.{}.sqlite".format(
                    organism, transcriptome))
    else:
        traninfo_connection = sqlite3.connect(
            "{0}/transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
                config.UPLOADS_DIR, owner, organism, transcriptome))

    # traninfo_connection = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips_annotations/{0}/{0}.{1}.sqlite".format(organism,transcriptome))
    traninfo_cursor = traninfo_connection.cursor()

    traninfo_dict = {}
    traninfo_cursor.execute(
        "SELECT transcript, cds_start, cds_stop,length,gene,sequence FROM transcripts;"
    )
    result = traninfo_cursor.fetchall()
    for row in result:
        tran = str(row[0])
        try:
            cds_start = int(row[1])
            cds_stop = int(row[2])
        except Exception:
            cds_start = row[1]
            cds_stop = row[2]
        length = int(row[3])
        gene = str(row[4])
        seq = str(row[5])
        traninfo_dict[tran] = {
            "cds_start": cds_start,
            "cds_stop": cds_stop,
            "length": length,
            "gene": gene,
            "coding_regions": [],
            "seq": seq
        }

    principal_transcripts = []
    if tranlist == "prin_trans":
        traninfo_cursor.execute(
            "SELECT transcript,gene FROM transcripts WHERE principal = 1;")
        result = traninfo_cursor.fetchall()
        for row in result:
            principal_transcripts.append(str(row[0]))
            tran_gene_dict[row[0]] = row[1].replace(",", "_")
    elif tranlist == "all_trans":
        traninfo_cursor.execute("SELECT transcript,gene FROM transcripts;")
        result = traninfo_cursor.fetchall()
        for row in result:
            principal_transcripts.append(str(row[0]))
            tran_gene_dict[row[0]] = row[1].replace(",", "_")
    elif tranlist == "custom_trans":
        principal_transcripts = user_defined_transcripts
        traninfo_cursor.execute(
            "SELECT transcript,gene FROM transcripts WHERE transcript IN ({});"
            .format(str(principal_transcripts).strip("[]")))
        result = traninfo_cursor.fetchall()
        for row in result:
            tran_gene_dict[row[0]] = row[1].replace(",", "_")

    transcriptome_info_dict = {}
    traninfo_cursor.execute(
        "SELECT transcript,strand,chrom from transcripts WHERE transcript IN ({});"
        .format(str(principal_transcripts).strip("[]")))
    result = traninfo_cursor.fetchall()
    for row in result:
        transcriptome_info_dict[str(row[0])] = {
            "strand": row[1],
            "chrom": row[2],
            "exons": []
        }
    traninfo_cursor.execute(
        "SELECT * from exons WHERE transcript IN ({});".format(
            str(principal_transcripts).strip("[]")))
    result = traninfo_cursor.fetchall()
    for row in result:
        transcriptome_info_dict[str(row[0])]["exons"].append((row[1], row[2]))
    logging.debug("building transcriptom info dict")
    sequence_dict = {}
    traninfo_cursor.execute(
        "SELECT transcript,sequence FROM transcripts WHERE transcript IN ({})".
        format(str(principal_transcripts).strip("[]")))
    result = traninfo_cursor.fetchall()
    for row in result:
        sequence_dict[row[0]] = row[1]
    # logging.debug("sequence dict keys",sequence_dict.keys())
    # Holds a list of all transcripts in accepted_orf_dict

    # logging.debug("accepted orf dict", accepted_orf_dict)
    logging.debug("accepted orf dict built")
    # Now build a profile for every transcript in accepted_transcripts

    if file_paths_dict["rnaseq"] == {}:
        if "te_check" in data:
            del data["te_check"]

    if file_paths_dict["riboseq"] == {} and file_paths_dict[
            "proteomics"] == {}:
        returnstr = "Error no files selected"
        return returnstr

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

    profile_dict, file_string, label_string = create_profiles(
        file_paths_dict, principal_transcripts, total_files, min_read_length,
        max_read_length)
    sorted_all_values, file_output_dict = extract_values(
        traninfo_dict, data, tran_gene_dict, selected_seq_types, profile_dict,
        min_fold_change, window, min_coverage, nuc_output)
    if sorted_all_values:
        return ("No results, try making filters less restrictive")

    # TODO change extension to csv if only one file
    filename = short_code + ".zip"
    returnstr = write_to_file(sorted_all_values, file_output_dict,
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
pausequery_blueprint = Blueprint("pausequery",
                                 __name__,
                                 template_folder="templates")


@pausequery_blueprint.route('/pausequery', methods=['POST'])
def pausequery():

    data = request.args.to_dict()
    print("pausequery called", data)
    user, logged_in = fetch_user()
    return find_pauses(data, user, logged_in)
