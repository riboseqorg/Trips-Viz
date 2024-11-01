import collections
import json
import os
import re  # TODO: replace it with re2
import subprocess
import time
from typing import Dict, List, Union

import polars as pl
from flask import Blueprint, render_template, request
from flask_login import current_user
from sqlitedict import SqliteDict

import config
import fixed_values
from core_functions import (build_profile, fetch_file_paths, fetch_files,
                            fetch_studies, fetch_study_info, fetch_user,
                            form_filler, generate_short_code)
from fixed_values import my_decoder
from sqlqueries_2 import get_table, get_user_id, sqlquery, table2dict


def get_nuc_comp_reads(sqlite_db: SqliteDict,
                       nuccomp_reads: str, data) -> Union[str, pl.DataFrame]:
    transhelve = "{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC,
                                                     config.ANNOTATION_DIR,
                                                     data["organism"], data["transcriptome"])
    if not os.path.isfile(transhelve):
        return "Cannot find annotation file {0}.{1}.sqlite".format(
            data["organism"], data["transcriptome"])

    if nuccomp_reads in sqlite_db["nuc_counts"]:
        return pl.DataFrame(sqlite_db["nuc_counts"][nuccomp_reads], schema=["readlen", "pos", "nuc", "count"])
    transcripts = sqlquery(transhelve, "transcripts").filter(
        pl.col('principle') & pl.col('cds_start') & pl.col('cds_stop')
        & pl.col('transcript').is_in(sqlite_db)).select(
            'transcript', 'cds_start', 'cds_stop', 'sequence').with_columns(
                pl.col('sequence').apply(lambda x: x.replace("T", "U"))

    )
    # TODO: Check some files for nan values in the table and remove them

    master_df = []
    offsets: Dict = sqlite_db["offsets"]["fiveprime"]["offsets"]
    for transcript in transcripts.iterrows(names=True):
        counts = sqlite_db[transcript["transcript"]]["unambig"]
        for readlen in counts:
            for pos, count in counts[readlen].items():
                offset_pos = (pos +
                              offsets[readlen] if readlen in offsets else pos + 15)
                if transcript["cds_start"] <= offset_pos <= transcript["cds_stop"]:
                    readframe = offset_pos % 3
                    cds_frame = (transcript["cds_start"] + 2) % 3

                    inframe = True if readframe == cds_frame else False

                    if not inframe or (nuccomp_reads == "offrame"):
                        continue
                    for i, c in enumerate(transcript["sequence"][pos:pos + readlen]):
                        master_df.append([readlen, i, c, count])

    try:
        sqlite_db["nuc_counts"][nuccomp_reads] = master_df
    except Exception:
        sqlite_db["nuc_counts"] = {nuccomp_reads: master_df}
    # save results so they won't have to be computed again later

    sqlite_db.commit()
    return pl.DataFrame(master_df, schema=[
        "readlen", "pos", "nuc", "count"]).groupby("readlen", "pos", "nuc").sum("count")


metainfo_plotpage_blueprint = Blueprint("metainfo_plotpage",
                                        __name__,
                                        template_folder="templates")


@metainfo_plotpage_blueprint.route(
    "/<organism>/<transcriptome>/metainfo_plot/")
def metainfo_plotpage(organism: str, transcriptome: str):
    # global user_short_passed

    data: Dict = form_filler(organism, transcriptome)
    accepted_studies = fetch_studies(data["gwips_info"][0, "organism_id"])
    data['files'] = fetch_files(accepted_studies).to_pandas()

    return render_template("metainfo_index.html", template_dict=data)


# Used to create custom metagene plots on the metainformation plot page
def create_custom_metagene(
    data: Dict,
    sqlite_db: Dict[str, Dict[str, Dict[str, str]]],
    metagene_frame: str,
    coverage: bool = False,
) -> str:
    data["custom_seq_list"] = data["custom_seq_list"].upper().replace(" ",
                                                                      "").replace("T", "U")
    custom_metagene_id = "_".join([
        "cmgc",
        data["custom_seq_list"].replace(",", "_"),
        data["custom_search_region"],
        data["exclude_first"],
        data["exclude_last"],
        data["include_first"],
        data["include_last"],
        data["excclude_first_val"],
        data["exclude_last_val"],
        data["include_first_val"],
        data["include_last_val"],
        data["metagene_tranlist"]]
    )
    offsets = sqlite_db["offsets"]
    transcripts = get_table('transcripts')
    # Apply filters
    transcripts['sequence'] = transcripts['sequence'].apply(
        lambda x: x.replace('T', 'U'))
    if not data["metagene_tranlist"]:
        transcripts = transcripts.filter(pl.col('principle'))
    else:
        data["metagene_tranlist"] = data["metagene_tranlist"].split(",")
        if "coverage" in data["metagene_tranlist"]:
            data["metagene_tranlist"].remove("coverage")
        transcripts = transcripts.filter(
            pl.col('transcript').is_in(data["metagene_tranlist"]))

    result = transcripts
    mgc = {"fiveprime": {}, "threeprime": {}}
    iupac_dict = fixed_values.iupac_dict.copy()

    subseq_list = []
    for subseq in data["custom_seq_list"].split(","):
        for ambig_nuc in fixed_values.ambig_nucs:
            subseq = subseq.replace(ambig_nuc, f"[{iupac_dict[ambig_nuc]}]")
        subseq_list.append(subseq)
    for row in result:
        tran = row[0]
        cds_start, cds_stop, min_pos, max_pos = [0] * 4
        if data["custom_search_region"] != "whole_gene":
            try:
                cds_start = int(row[1])
                cds_stop = int(row[2])
            except Exception:
                continue
        seq = row[3].replace("T", "U")
        if data["custom_search_region"] == "whole_gene":
            min_pos = 0
            max_pos = len(seq)
        elif data["custom_search_region"] == "five_leader":
            min_pos = 0
            max_pos = cds_start
        elif data["custom_search_region"] == "cds":
            min_pos = cds_start
            max_pos = cds_stop
        elif data["custom_search_region"] == "three_trailer":
            min_pos = cds_stop
            max_pos = len(seq)
        if data["include_first"]:
            max_pos = min_pos + data["include_first_val"]
        if data["include_last"]:
            min_pos = max_pos - data["include_last_val"]
        if data["exclude_first"]:
            min_pos = min_pos + data["excclude_first_val"]
        if data["exclude_last"]:
            max_pos = max_pos - data["exclude_last_val"]

        seq_positions = []
        if cds_start:
            for subseq in subseq_list:
                pattern = re.compile(r"{}".format(subseq))
                m = True
                search_pos = min_pos - 1

                while m:
                    m = pattern.search(seq, search_pos, max_pos)
                    if m:
                        position = m.span()[0]
                        seq_positions.append(position)
                        search_pos = position + 1

        if seq_positions:
            profile = {"fiveprime": {}, "threeprime": {}}
            try:
                tran_reads = sqlite_db[tran]["unambig"]
            except Exception:
                tran_reads = {"unambig": {}}

            for readlen in tran_reads:
                if readlen in offsets["fiveprime"]["offsets"]:
                    offset = offsets["fiveprime"]["offsets"][readlen]
                else:
                    offset = 15
                profile["fiveprime"][readlen] = {}
                profile["threeprime"][readlen] = {}
                for pos in tran_reads[readlen]:
                    relative_frame = ((pos + offset) - cds_start) % 3
                    if ((metagene_frame == "minus_frame"
                         and relative_frame != 1) or
                        (metagene_frame == "in_frame" and relative_frame != 2)
                            or (metagene_frame == "plus_frame"
                                and relative_frame != 0)):
                        continue
                    three_pos = pos + readlen
                    count = tran_reads[readlen][pos]
                    if not coverage:
                        try:
                            profile["fiveprime"][readlen][pos] += count
                        except Exception:
                            profile["fiveprime"][readlen][pos] = 0
                            profile["fiveprime"][readlen][pos] += count
                        try:
                            profile["threeprime"][readlen][three_pos] += count
                        except Exception:
                            profile["threeprime"][readlen][three_pos] = 0
                            profile["threeprime"][readlen][three_pos] += count
                    else:
                        for x in range(pos, pos + readlen):
                            try:
                                profile["fiveprime"][readlen][x] += count
                            except Exception:
                                profile["fiveprime"][readlen][x] = 0
                                profile["fiveprime"][readlen][x] += count
                            try:
                                profile["threeprime"][readlen][x] += count
                            except Exception:
                                profile["threeprime"][readlen][x] = 0
                                profile["threeprime"][readlen][x] += count

            for readlen in profile["fiveprime"]:
                if readlen not in mgc["fiveprime"]:
                    mgc["fiveprime"][readlen] = {}
                    for i in range(-600, 601):
                        mgc["fiveprime"][readlen][i] = 0
                for pos in profile["fiveprime"][readlen]:
                    count = profile["fiveprime"][readlen][pos]
                    for seq_position in seq_positions:
                        if seq_position >= pos - 600 and seq_position <= pos + 600:
                            relative_seq_position = pos - seq_position
                            mgc["fiveprime"][readlen][
                                relative_seq_position] += count

            for readlen in profile["threeprime"]:
                if readlen not in mgc["threeprime"]:
                    mgc["threeprime"][readlen] = {}
                    for i in range(-600, 601):
                        mgc["threeprime"][readlen][i] = 0
                for pos in profile["threeprime"][readlen]:
                    count = profile["threeprime"][readlen][pos]
                    for seq_position in seq_positions:
                        if seq_position >= pos - 600 and seq_position <= pos + 600:
                            relative_seq_position = pos - seq_position
                            mgc["threeprime"][readlen][
                                relative_seq_position] += count
        else:
            pass
    sqlite_db[custom_metagene_id] = mgc
    sqlite_db.commit()
    return mgc


def redo_periodicity_plots(file_path: str) -> None:
    traninfo_dict = {}
    transcripts = get_table("transcripts").filter(
        pl.col('principle') & pl.col('tran_type')).select(
            "transcript", "cds_start", "cds_stop", "length")
    for row in result:
        traninfo_dict[str(row[0])] = [int(row[1]), int(row[2]), int(row[3])]
    trip_periodicity_reads = 0
    master_read_dict = SqliteDict(f"{file_path}",
                                  autocommit=False,
                                  decode=my_decoder)
    master_trip_dict = {"threeprime": {}, "fiveprime": {}}

    for tran in traninfo_dict:
        if tran not in master_read_dict:
            continue
        cds_start = traninfo_dict[tran][0]
        cds_stop = traninfo_dict[tran][1]
        length = traninfo_dict[tran][2]
        if cds_start > 1 and cds_stop < length:
            for primetype in ["fiveprime", "threeprime"]:
                for readlength in master_read_dict[tran]["unambig"]:
                    # for each fiveprime postion for this readlength within
                    #    this transcript
                    for raw_pos in master_read_dict[tran]["unambig"][
                            readlength]:
                        trip_periodicity_reads += 1
                        # get the five prime postion minus the cds start postion
                        real_pos = (raw_pos -
                                    cds_start if primetype == "fiveprime" else
                                    (raw_pos + readlength) - cds_start)
                        if real_pos >= cds_start and real_pos <= cds_stop:
                            readcount = master_read_dict[tran]["unambig"][
                                readlength][raw_pos]
                            frame = real_pos % 3
                            if readlength in master_trip_dict[primetype]:
                                master_trip_dict[primetype][readlength][str(
                                    frame)] += readcount
                            else:
                                master_trip_dict[primetype][readlength] = {
                                    "0": 0.0,
                                    "1": 0.0,
                                    "2": 0.0,
                                }
                                master_trip_dict[primetype][readlength][str(
                                    frame)] += readcount
    master_read_dict["trip_periodicity"] = master_trip_dict
    master_read_dict.commit()


def metainfoquery(data) -> str | tuple:
    # global user_short_passed
    gene_dict = {}

    if not data["file_list"]:
        return "No files selected"
    file_paths = fetch_file_paths(data)
    if file_paths.is_empty():
        return "Files not found, please contact the site admin"
    # TODO: Select files with ids
    # TODO:: Discuss with pash tha sequence read should be sequence type specific or combined
    traninfo_connection = (
        "/home/DATA/www/tripsviz/tripsviz/trips_annotations/{0}/{0}.{1}.sqlite".format(data["organism"], data["transcriptome"]))

    transcripts = sqlquery(traninfo_connection, "transcripts")
    if data["plottype"] == "readlen_dist":
        master_dict = []
        if data["metagene_list"]:
            data["metagene_list"] = data["metagene_list"].split(",")
        if not (not data["metagene_list"] and data["custom_search_region"] == "whole_gene"):
            if not data["metagene_list"]:
                transcripts = transcripts.filter(
                    pl.col("principal") == 1
                )
            else:
                transcripts = transcripts.filter(
                    pl.col("transcript").is_in(data["metagene_list"])
                )
        for filepath in file_paths["path"]:
            sqlite_db = SqliteDict(filepath, autocommit=False)
            # If no transcripts given and no region specified, get the precomputed read lengths (all transcripts, entire gene)
            if not data["metagene_list"] and data["custom_search_region"] == "whole_gene":
                if data["readlen_ambig"] and "read_lengths" in sqlite_db:
                    read_lengths = sqlite_db["read_lengths"]
                elif "unambig_read_lengths" in sqlite_db:
                    read_lengths = sqlite_db["unambig_read_lengths"]

                else:
                    sqlite_db.close()
                    return ("No readlength distribution data for this file, please report this to tripsvizsite@gmail.com or via the contact page.")

                sqlite_db.close()
                master_dict = pl.DataFrame(
                    {'pos': read_lengths.keys(), 'count': read_lengths.values()})
            else:
                transcripts = transcripts.filter(
                    pl.col("transcript").is_in(sqlite_db)
                )
                master_dict = []
                for transcript in transcripts & set(sqlite_db):
                    counts = sqlite_db[transcript]["unambig"]
                    master_dict.append(pl.DataFrame([(pos, counts)
                                                     for pos_dict in counts.values()
                                                     for pos, counts in pos_dict.items()
                                                     ], schema=["pos", "counts"]))

                master_dict = pl.concat(master_dict).groupby("readlen").agg([
                    pl.col("counts").sum()
                ]).sort("readlen")
                if data["custom_search_region"] == "five_leader":
                    master_dict = master_dict.filter(pl.col("pos") < cds_start)
                elif data["custom_search_region"] == "three_trailer":
                    master_dict = master_dict.filter(pl.col("pos") > cds_stop)
                elif data["custom_search_region"] == "cds":
                    master_dict = master_dict.filter(
                        (pl.col("pos") > cds_start) & (pl.col("pos") < cds_stop))
                else:
                    pass

        return metainfo_plots.readlen_dist(master_dict)

    if data["trip_minreadlen"] == "redo":
        redo_periodicity = True
        trip_minreadlen = 25
    else:
        trip_minreadlen = int(data["trip_minreadlen"])
        redo_periodicity = False
    color_palette = data["color_palette"]
    minimum_reads = int(data["minimum_reads"])
    raw_te_tranlist = data["te_tranlist"]
    raw_te_tranlist = raw_te_tranlist.replace(" ", ",")
    te_tranlist = []
    for item in raw_te_tranlist.split(","):
        te_tranlist.append(item)

    user_settings = config.DEFAULT_USER_SETTINGS.copy()

    # get a list of organism id's this user can access
    if current_user.is_authenticated:
        # get user_id
        user_name = current_user.name
        user_id = get_user_id(user_name)
        user_settings = get_table('user_settings').filter(
            pl.col('user_id') == user_id
        )

    user_short_passed = True

    if data["html_args"]["user_short"] == "None" or user_short_passed:
        short_code = generate_short_code(data)
    else:
        short_code = data["html_args"]["user_short"]
        user_short_passed = True

    owner = get_table("organism").filter((pl.col("organism_name") == data["organism"]) & (
        pl.col("transcriptome_list") == data["transcriptome"]))[0, "owner"]
    sqlfile = "{0}transcriptomes/{1}/{2}/{2}_{3}.sqlite".format(
        config.UPLOADS_DIR, owner, data["organism"], data["transcriptome"])
    if owner:
        sqlfile = "{0}/{1}/{2}/{2}.{3}.sqlite".format(
            config.SCRIPT_LOC, config.ANNOTATION_DIR, data["organism"], data["transcriptome"])
        if not os.path.isfile(sqlfile):
            return ("Cannot find annotation file {}.{}.sqlite".format(data["organism"], data["transcriptome"]))
    transcripts = sqlquery(sqlfile, "transcripts")
    if data["plottype"] == "te":

        if te_tranlist == "" or te_tranlist == ['']:
            te_tranlist = None
        if count_type != "tpm":
            if not te_tranlist:
                transcripts = transcripts.filter(pl.col("principal"))
                if region != "all":
                    transcripts = transcripts.filter(pl.col("tran_type") == 1)
            else:
                transcripts = transcripts.filter(
                    pl.col("transcript").is_in(te_tranlist)
                )

        else:
            if len(data["file_list"]) > 200:  # TODO: This can be chnage to higher value
                return ("Error: A maximum of 200 files can be used when TPM is selected.")

        if not te_tranlist:
            if len(data["file_list"]) > 30 and not count_agg:
                return ("Error: A maximum of 30 files can be selected if not aggregating counts and using all transcripts. Reduce number of selected files or click the 'Aggregate counts' checkbox at the top right of the page. Alternatively input a list of transcripts in the 'Transcript list' box.")

                # User may have passed list of genes instead of transcripts
        if not transcript_list and te_tranlist:
            transcripts = transcripts.filter(
                pl.col("gene").is_in(te_tranlist) & pl.col("principal")
            )

        traninfo_dict = {}
        for result in transcript_list:
            if result[0]:
                traninfo_dict[result[0]] = {"transcript": result[0], "gene": result[1].replace(",", "_").replace(";", "_"), "length": result[2], "cds_start": result[3], "cds_stop": result[4], "seq": result[5].upper(),
                                            "strand": result[6], "stop_list": result[7].split(","), "start_list": result[8].split(","), "exon_junctions": result[9].split(","),
                                            "tran_type": result[10], "principal": result[11]}

        longest_tran_list = traninfo_dict.keys()
        if count_agg:
            table_str = aggregate_counts(file_paths_dict, traninfo_dict, longest_tran_list, region,
                                         organism, all_seq_types, te_minimum_reads, data["html_args"], count_type, te_tranlist)
        else:
            table_str = sample_counts(file_paths_dict, traninfo_dict, longest_tran_list, region,
                                      organism, all_seq_types, te_minimum_reads, data["html_args"], count_type, te_tranlist)

        return metainfo_plots.te_table(table_str)

    elif data["plottype"] == "mrna_dist":
        longest_tran_list = []
        cds_dict = {}

        transcripts = transcripts.filter(
            pl.col("principal") & pl.col("tran_type"))
        # "SELECT transcript,cds_start,cds_stop from transcripts where principal = 1 and tran_type = 1;")
        mrna_dist_dict = {}
        for filepath in file_paths["path"]:
            filename = filepath.split("/")[-1].split(".sqlite")[0]
            mrna_dist_dict[filename] = {"5_leader": 0,
                                        "start_codon": 0,
                                        "cds": 0,
                                        "stop_codon": 0,
                                        "3_trailer": 0,
                                        "total": 0}
            if os.path.isfile(filepath):
                sqlite_db = SqliteDict(filepath, autocommit=False)
            else:
                return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath))
            if "mrna_dist_dict" in sqlite_db:
                mrna_dist_dict[filename]["5_leader"] = sqlite_db["mrna_dist_dict"]["5_leader"]
                mrna_dist_dict[filename]["start_codon"] = sqlite_db["mrna_dist_dict"]["start_codon"]
                mrna_dist_dict[filename]["cds"] = sqlite_db["mrna_dist_dict"]["cds"]
                mrna_dist_dict[filename]["stop_codon"] = sqlite_db["mrna_dist_dict"]["stop_codon"]
                mrna_dist_dict[filename]["3_trailer"] = sqlite_db["mrna_dist_dict"]["3_trailer"]
                mrna_dist_dict[filename]["total"] = float(sqlite_db["mrna_dist_dict"]["5_leader"]+sqlite_db["mrna_dist_dict"]["start_codon"] +
                                                          sqlite_db["mrna_dist_dict"]["cds"]+sqlite_db["mrna_dist_dict"]["stop_codon"]+sqlite_db["mrna_dist_dict"]["3_trailer"])
            if not mrna_dist_dict[filename]["total"]:
                for transcript in longest_tran_list:
                    try:
                        transcript_dict = sqlite_db[transcript]["unambig"]
                    except KeyError:
                        continue
                    try:
                        cds_start = cds_dict[transcript]["cds_start"]
                        cds_stop = cds_dict[transcript]["cds_stop"]
                    except KeyError:
                        continue
                    for readlen in transcript_dict:
                        for five_pos in transcript_dict[readlen]:
                            three_pos = five_pos+readlen
                            if three_pos <= cds_start+3:
                                mrna_dist_dict[filename]["5_leader"] += transcript_dict[readlen][five_pos]
                            elif five_pos <= cds_start-4 and three_pos >= cds_start+4:
                                mrna_dist_dict[filename]["start_codon"] += transcript_dict[readlen][five_pos]
                            elif five_pos >= cds_start-3 and three_pos <= cds_stop-2:
                                mrna_dist_dict[filename]["cds"] += transcript_dict[readlen][five_pos]
                            elif five_pos <= cds_stop-9 and three_pos >= cds_stop-1:
                                mrna_dist_dict[filename]["stop_codon"] += transcript_dict[readlen][five_pos]
                            elif five_pos >= cds_stop-9:
                                mrna_dist_dict[filename]["3_trailer"] += transcript_dict[readlen][five_pos]
                        sqlite_db["mrna_dist_dict"] = {"5_leader": mrna_dist_dict[filename]["5_leader"],
                                                       "start_codon": mrna_dist_dict[filename]["start_codon"],
                                                       "cds": mrna_dist_dict[filename]["cds"],
                                                       "stop_codon": mrna_dist_dict[filename]["stop_codon"],
                                                       "3_trailer": mrna_dist_dict[filename]["3_trailer"],
                                                       "total": (mrna_dist_dict[filename]["5_leader"]+mrna_dist_dict[filename]["start_codon"]+mrna_dist_dict[filename]["cds"]+mrna_dist_dict[filename]["stop_codon"]+mrna_dist_dict[filename]["3_trailer"])}
                    sqlite_db.commit()
            sqlite_db.close()

        return metainfo_plots.mrna_dist(mrna_dist_dict, short_code, background_col, title_size, axis_label_size, subheading_size, marker_size, mrna_dist_per, md_start, md_stop, legend_size)

    if data["plottype"] == "mrna_dist_readlen":
        minreadlen = 15
        maxreadlen = 100
        if owner == 1:
            traninfo_dict = SqliteDict("{0}/{1}/{2}/{2}.{3}.sqlite".format(
                config.SCRIPT_LOC, config.ANNOTATION_DIR, data["organism"], data["transcriptome"]), autocommit=False)
        else:
            traninfo_dict = SqliteDict("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
                config.UPLOADS_DIR, owner, data["organism"], data["transcriptome"]), autocommit=False)
        longest_tran_list = []
        cds_dict = {}

        transcripts = transcripts.filter(
            pl.col("principal") & pl.col("tran_type"))
        for row in transcripts:
            longest_tran_list.append(str(row[0]))
            cds_dict[str(row[0])] = {"cds_start": int(
                row[1]), "cds_stop": int(row[2])}
        mrna_dist_dict = {"5_leader": collections.OrderedDict(),
                          "start_codon": collections.OrderedDict(),
                          "cds": collections.OrderedDict(),
                          "stop_codon": collections.OrderedDict(),
                          "3_trailer": collections.OrderedDict(),
                          "total": collections.OrderedDict()}

        for filetype in file_paths_dict:
            for file_id in file_paths_dict[filetype]:
                filepath = file_paths_dict[filetype][file_id]
                try:
                    sqlite_db = SqliteDict(filepath, autocommit=False)
                except:
                    return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath))

                if "mrna_dist_readlen_dict" in sqlite_db:
                    for readlen in range(minreadlen, maxreadlen+1):
                        if readlen not in mrna_dist_dict["5_leader"]:
                            mrna_dist_dict['5_leader'][readlen] = 0
                            mrna_dist_dict['start_codon'][readlen] = 0
                            mrna_dist_dict['cds'][readlen] = 0
                            mrna_dist_dict['stop_codon'][readlen] = 0
                            mrna_dist_dict['3_trailer'][readlen] = 0
                            mrna_dist_dict['total'][readlen] = 0
                            mrna_dist_dict['5_leader'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["5_leader"][readlen]
                            mrna_dist_dict['start_codon'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["start_codon"][readlen]
                            mrna_dist_dict['cds'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["cds"][readlen]
                            mrna_dist_dict['stop_codon'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["stop_codon"][readlen]
                            mrna_dist_dict['3_trailer'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["3_trailer"][readlen]
                            mrna_dist_dict['total'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["total"][readlen]
                else:
                    file_specific_mrna_dist_dict = {"5_leader": {},
                                                    "start_codon": {},
                                                    "cds": {},
                                                    "stop_codon": {},
                                                    "3_trailer": {},
                                                    "total": {}}
                    for readlen in range(minreadlen, maxreadlen+1):
                        file_specific_mrna_dist_dict['5_leader'][readlen] = 0
                        file_specific_mrna_dist_dict['start_codon'][readlen] = 0
                        file_specific_mrna_dist_dict['cds'][readlen] = 0
                        file_specific_mrna_dist_dict['stop_codon'][readlen] = 0
                        file_specific_mrna_dist_dict['3_trailer'][readlen] = 0
                        file_specific_mrna_dist_dict['total'][readlen] = 0
                    for transcript in longest_tran_list:
                        try:
                            transcript_dict = sqlite_db[transcript]["unambig"]
                        except KeyError:
                            continue
                        try:
                            cds_start = cds_dict[transcript]["cds_start"]
                            cds_stop = cds_dict[transcript]["cds_stop"]
                        except KeyError:
                            continue
                    for readlen in range(minreadlen, maxreadlen+1):
                        if readlen not in mrna_dist_dict["5_leader"]:
                            mrna_dist_dict['5_leader'][readlen] = 0
                            mrna_dist_dict['start_codon'][readlen] = 0
                            mrna_dist_dict['cds'][readlen] = 0
                            mrna_dist_dict['stop_codon'][readlen] = 0
                            mrna_dist_dict['3_trailer'][readlen] = 0
                            mrna_dist_dict['total'][readlen] = 0
                        if readlen not in transcript_dict:
                            continue
                        for five_pos in transcript_dict[readlen]:
                            three_pos = five_pos+readlen
                            if three_pos <= cds_start+3:
                                mrna_dist_dict["5_leader"][readlen] += transcript_dict[readlen][five_pos]
                                file_specific_mrna_dist_dict["5_leader"][readlen] += transcript_dict[readlen][five_pos]
                            elif five_pos <= cds_start-4 and three_pos >= cds_start+4:
                                mrna_dist_dict["start_codon"][readlen] += transcript_dict[readlen][five_pos]
                                file_specific_mrna_dist_dict["start_codon"][readlen] += transcript_dict[readlen][five_pos]
                            elif five_pos >= cds_start-3 and three_pos <= cds_stop-2:
                                mrna_dist_dict["cds"][readlen] += transcript_dict[readlen][five_pos]
                                file_specific_mrna_dist_dict["cds"][readlen] += transcript_dict[readlen][five_pos]
                            elif five_pos <= cds_stop-9 and three_pos >= cds_stop-1:
                                mrna_dist_dict["stop_codon"][readlen] += transcript_dict[readlen][five_pos]
                                file_specific_mrna_dist_dict["stop_codon"][readlen] += transcript_dict[readlen][five_pos]
                            elif five_pos >= cds_stop-9:
                                mrna_dist_dict["3_trailer"][readlen] += transcript_dict[readlen][five_pos]
                                file_specific_mrna_dist_dict["3_trailer"][readlen] += transcript_dict[readlen][five_pos]
                    sqlite_db["mrna_dist_readlen_dict"] = file_specific_mrna_dist_dict
                    sqlite_db.commit()
                    sqlite_db.close()

        # If mrna_readlen_per is true normalize everything over the max in it's category
        if mrna_readlen_per:
            # For each category find the max value, max value will default to 1 if there is no values in that category.
            max_five = max(1, float(max(mrna_dist_dict["5_leader"].values())))
            max_start = max(
                1, float(max(mrna_dist_dict["start_codon"].values())))
            max_cds = max(1, float(max(mrna_dist_dict["cds"].values())))
            max_stop = max(
                1, float(max(mrna_dist_dict["stop_codon"].values())))
            max_three = max(
                1, float(max(mrna_dist_dict["3_trailer"].values())))
            for readlen in mrna_dist_dict["5_leader"]:
                mrna_dist_dict["5_leader"][readlen] = (
                    float(mrna_dist_dict["5_leader"][readlen])/max_five)*100
            for readlen in mrna_dist_dict["start_codon"]:
                mrna_dist_dict["start_codon"][readlen] = (
                    float(mrna_dist_dict["start_codon"][readlen])/max_start)*100
            for readlen in mrna_dist_dict["cds"]:
                mrna_dist_dict["cds"][readlen] = (
                    float(mrna_dist_dict["cds"][readlen])/max_cds)*100
            for readlen in mrna_dist_dict["stop_codon"]:
                mrna_dist_dict["stop_codon"][readlen] = (
                    float(mrna_dist_dict["stop_codon"][readlen])/max_stop)*100
            for readlen in mrna_dist_dict["3_trailer"]:
                mrna_dist_dict["3_trailer"][readlen] = (
                    float(mrna_dist_dict["3_trailer"][readlen])/max_three)*100

        return metainfo_plots.mrna_dist_readlen(mrna_dist_dict, mrna_readlen_per, short_code, background_col, title_size, axis_label_size, subheading_size, marker_size, legend_size)

    if data["plottype"] == "replicate_comp":
        prin_tran_list = [row[0] for row in result if row[1] == 1]
        minimum_reads = int(minimum_reads)
        mapped_reads_dict = {}
        factor_dict = {}
        if data["normalise"]:
            for filetype in file_paths_dict:
                for file_id in file_paths_dict[filetype]:
                    filepath = file_paths_dict[filetype][file_id]
                    sqlite_db = SqliteDict(filepath, autocommit=False)
                    mapped_reads = sqlite_db["coding_counts"]
                    if mapped_reads == None or mapped_reads == 0:
                        connection.close()
                        return ("Error: Mapped reads info missing for one or more files, cannot normalise")
                    mapped_reads_dict[file_id] = float(mapped_reads)
            minval = min(mapped_reads_dict.values())
            for file_id in mapped_reads_dict:
                factor = minval/mapped_reads_dict[file_id]
                factor_dict[file_id] = factor
        if minimum_reads > 0:
            min_log_val = log(minimum_reads, 2)
        else:
            min_log_val = 0
        labels = []
        transcript_dict = {}
        for filetype in file_paths_dict:
            for file_id in file_paths_dict[filetype]:
                cursor.execute(
                    "SELECT file_description from files where file_id = '{}';".format(file_id))
                result = cursor.fetchone()
                label = result[0]
                lbl_tag = 0
                while label in labels:
                    lbl_tag += 1
                    label = result[0] + str(lbl_tag)
                labels.append(label)
                filepath = file_paths_dict[filetype][file_id]

                if os.path.isfile(filepath):
                    sqlite_db = SqliteDict(filepath, autocommit=False)
                    opendict = sqlite_db["unambiguous_all_totals"]
                    sqlite_db.close()
                else:
                    return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath))
                for transcript in prin_tran_list:
                    if transcript not in transcript_dict:
                        transcript_dict[transcript] = {}
                    if transcript in opendict:
                        try:
                            count = opendict[transcript]
                            if normalise:
                                count = float(count)*factor_dict[file_id]
                            if count >= min_log_val:
                                transcript_dict[transcript][label] = log(
                                    count, 2)
                        except:
                            pass
        del_list = []
        for transcript in transcript_dict:
            if len(transcript_dict[transcript]) != len(labels):
                del_list.append(transcript)
        for tran in del_list:
            del transcript_dict[tran]

        return metainfo_plots.replicate_comp(labels, transcript_dict, min_log_val, short_code, background_col, str(title_size)+"pt", str(axis_label_size)+"pt", str(subheading_size)+"pt", str(marker_size)+"pt", data["coor_type"])
    if data["plottype"] == "nuc_comp":
        master_count_dict = {"A": collections.OrderedDict(), "T": collections.OrderedDict(
        ), "G": collections.OrderedDict(), "C": collections.OrderedDict()}
        for filetype in file_paths_dict:
            for file_id in file_paths_dict[filetype]:
                filepath = file_paths_dict[filetype][file_id]
                if os.path.isfile(filepath):
                    sqlite_db = SqliteDict(filepath, autocommit=False)
                else:
                    return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath))
                if "nuc_counts" not in sqlite_db:
                    return ("No nucleotide counts data for this file, please report this to tripsvizsite@gmail.com or via the contact page.")
                if data["nuc_comp_direction"] == "nuc_comp_five":
                    if "nuc_counts" in sqlite_db:
                        if data["nuccomp_reads"] in sqlite_db["nuc_counts"]:
                            nuc_counts = sqlite_db["nuc_counts"][data["nuccomp_reads"]]
                        else:
                            nuc_counts = get_nuc_comp_reads(
                                sqlite_db, data["nuccomp_reads"], data)
                    else:
                        nuc_counts = get_nuc_comp_reads(
                            sqlite_db, data["nuccomp_reads"], data)
                elif data["nuc_comp_direction"] == "nuc_comp_three":
                    if "threeprime_nuc_counts" in sqlite_db:
                        if data["nuccomp_reads"] in sqlite_db["threeprime_nuc_counts"]:
                            nuc_counts = sqlite_db["threeprime_nuc_counts"][data["nuccomp_reads"]]
                        else:
                            nuc_counts = get_nuc_comp_reads(
                                sqlite_db, data["nuccomp_reads"], data)
                    else:
                        nuc_counts = get_nuc_comp_reads(
                            sqlite_db, data["nuccomp_reads"], data)
                # return str(nuccomp_reads)
                if data["nuc_comp_direction"] == "nuc_comp_five":
                    for readlen in range(nuc_minreadlen, nuc_maxreadlen+1):
                        for i in range(0, readlen+1):
                            for nuc in "ATGC":
                                if i not in master_count_dict[nuc]:
                                    master_count_dict[nuc][i] = 0
                                if readlen in nuc_counts:
                                    if i in nuc_counts[readlen]:
                                        master_count_dict[nuc][i] += nuc_counts[readlen][i][nuc]
                elif data["nuc_comp_direction"] == "nuc_comp_three":
                    for readlen in range(nuc_minreadlen, nuc_maxreadlen+1):
                        for i in range(-1, -(nuc_maxreadlen), -1):
                            for nuc in "ATGC":
                                if i not in master_count_dict[nuc]:
                                    master_count_dict[nuc][i] = 0
                                if readlen in nuc_counts:
                                    if i in nuc_counts[readlen]:
                                        master_count_dict[nuc][i] += nuc_counts[readlen][i][nuc]
                sqlite_db.commit()
                sqlite_db.close()
        master_dict = {"A": collections.OrderedDict(),
                       "T": collections.OrderedDict(),
                       "G": collections.OrderedDict(),
                       "C": collections.OrderedDict()}
        master_dict["A"][0] = 0
        master_dict["T"][0] = 0
        master_dict["G"][0] = 0
        master_dict["C"][0] = 0
        if data["nuc_comp_direction"] == "nuc_comp_five":
            for nuc in "ATGC":
                for i in range(0, nuc_maxreadlen):
                    if i in master_count_dict[nuc]:
                        thiscount = master_count_dict[nuc][i]
                        othercount = 0.01
                        for subnuc in "ATGC":
                            othercount += master_count_dict[subnuc][i]
                        if data["nuc_comp_type"] == "nuc_comp_per":
                            master_dict[nuc][i] = (
                                float(thiscount)/float(othercount))*100
                        elif data["nuc_comp_type"] == "nuc_comp_count":
                            master_dict[nuc][i] = float(thiscount)
        elif data["nuc_comp_direction"] == "nuc_comp_three":
            for nuc in "ATGC":
                for i in range(-1, -(nuc_maxreadlen), -1):
                    if i in master_count_dict[nuc]:
                        thiscount = master_count_dict[nuc][i]
                        othercount = 0.01
                        for subnuc in "ATGC":
                            othercount += master_count_dict[subnuc][i]
                        if data["nuc_comp_type"] == "nuc_comp_per":
                            master_dict[nuc][i] = (
                                float(thiscount)/float(othercount))*100
                        elif data["nuc_comp_type"] == "nuc_comp_count":
                            master_dict[nuc][i] = float(thiscount)
        title = "Nucleotide composition"

        return metainfo_plots.nuc_comp(master_dict, nuc_maxreadlen, title, data["nuc_comp_type"], data["nuc_comp_direction"], short_code, background_col, a_col, t_col, g_col, c_col, title_size, axis_label_size, subheading_size, marker_size, legend_size)

    elif data["plottype"] == "dinuc_bias":
        master_count_dict = collections.OrderedDict([("AA", 0), ("AT", 0), ("AG", 0), ("AC", 0),
                                                     ("TA", 0), ("TT",
                                                                 0), ("TG", 0), ("TC", 0),
                                                     ("GA", 0), ("GT",
                                                                 0), ("GG", 0), ("GC", 0),
                                                     ("CA", 0), ("CT", 0), ("CG", 0), ("CC", 0)])
        master_count_dict = []
        for filetype in file_paths_dict:
            for file_id in file_paths_dict[filetype]:
                filepath = file_paths_dict[filetype][file_id]
                if os.path.isfile(filepath):
                    sqlite_db = SqliteDict(filepath, autocommit=False)
                else:
                    return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath))
                dinuc_counts = sqlite_db["dinuc_counts"]
                # TODO: put dinuc count as list instead of dict
                for readlen in dinuc_counts:
                    
                    for dinuc,count in dinuc_counts[readlen].items():
                        master_count_dict.append([dinuc, count]) 

            return metainfo_plots.dinuc_bias(master_count_dict, short_code, background_col, title_size, axis_label_size, subheading_size, marker_size)

        else:
            return ("No fastq_screen file available for this dataset")
    elif data["plottype"] == "metagene_plot":

        minpos = -300
        maxpos = 300
        pos_list = []
        mapped_reads_dict = {}
        for i in range(minpos, maxpos+1):
            pos_list.append(i)
        count_dict = {"fiveprime": {}, "threeprime": {}}
        for primetype in ["fiveprime", "threeprime"]:
            for i in range(minpos, maxpos+1):
                count_dict[primetype][i] = 0

        if metagene_aggregate:
            fiveprime_counts = []
            threeprime_counts = []
        else:
            if len(data["file_list"]) > 6:
                return ("Can only choose a maximum of 6 files if not using aggregate option")
            fiveprime_counts = {}
            threeprime_counts = {}

        for filetype in file_paths_dict:
            for file_id in file_paths_dict[filetype]:
                cursor.execute(
                    "SELECT file_description from files where file_id = '{}';".format(file_id))
                file_desc = cursor.fetchone()[0]
                filepath = file_paths_dict[filetype][file_id]
                if os.path.isfile(filepath):
                    sqlite_db = SqliteDict(filepath, autocommit=False)
                    if metagene_normalise:
                        mapped_reads = sqlite_db["coding_counts"] + \
                            sqlite_db["noncoding_counts"]
                        mapped_reads_dict[file_desc] = float(mapped_reads)
                else:
                    connection.close()
                    return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath))
                if data["metagene_type"] != "metagene_custom":
                    if "metagene_counts" not in sqlite_db:
                        connection.close()
                        return ("No metagene counts data for this file, please report this to tripsvizsite@gmail.com or via the contact page.")
                if data["metagene_type"] == "metagene_start":
                    if data["metagene_list"] == "":
                        mgc = sqlite_db["metagene_counts"]
                    else:
                        if "coverage" in data["metagene_list"]:
                            mgc = create_custom_metagene("AUG", 0, 0, 3, 0, "cds", False, False, True, False, sqlite_db,
                                                         organism, data["metagene_list"], "All", transcriptome, transhelve, coverage=True)
                        else:
                            mgc = create_custom_metagene("AUG", 0, 0, 3, 0, "cds", False, False, True, False,
                                                         sqlite_db, organism, data["metagene_list"], "All", transcriptome, transhelve)
                elif data["metagene_type"] == "metagene_stop":
                    if data["metagene_list"] == "":
                        mgc = sqlite_db["stop_metagene_counts"]
                    else:
                        if "coverage" in data["metagene_list"]:
                            mgc = create_custom_metagene("UAG,UAA,UGA", 0, 0, 0, 3, "cds", False, False, False, True,
                                                         sqlite_db, organism, data["metagene_list"], "All", transcriptome, transhelve, coverage=True)
                        else:
                            mgc = create_custom_metagene("UAG,UAA,UGA", 0, 0, 0, 3, "cds", False, False, False,
                                                         True, sqlite_db, organism, data["metagene_list"], "All", transcriptome, transhelve)
                elif data["metagene_type"] == "metagene_second_aug":
                    mgc = sqlite_db["secondary_metagene_counts"]
                elif data["metagene_type"] == "metagene_custom":
                    mgc = create_custom_metagene(custom_seq_list, exclude_first_val, exclude_last_val, include_first_val, include_last_val,
                                                 data["custom_search_region"], exclude_first, exclude_last, include_first, include_last, sqlite_db, organism, data["metagene_list"], data["metagene_frame"], transcriptome, transhelve)
                    if custom_seq_list == "AUG" and data["custom_search_region"] == "cds" and include_first_val == 3 and data["metagene_list"] == "":
                        mod_mgc = {}
                        for key in mgc:
                            if key != "unambig":
                                mod_mgc[key] = mgc[key]
                        sqlite_db["metagene_counts"] = mod_mgc
                        sqlite_db.commit()
                if metagene_offsets:
                    new_mgc = {"fiveprime": {}, "threeprime": {}}
                    offsets = sqlite_db["offsets"]
                    for readlen in mgc["fiveprime"]:
                        new_mgc["fiveprime"][readlen] = {}
                        if readlen in offsets["fiveprime"]["offsets"]:
                            offset = offsets["fiveprime"]["offsets"][readlen]
                        else:
                            offset = 15
                        for pos in mgc["fiveprime"][readlen]:
                            count = mgc["fiveprime"][readlen][pos]
                            new_pos = pos+offset
                            if new_pos not in new_mgc["fiveprime"][readlen]:
                                new_mgc["fiveprime"][readlen][new_pos] = 0
                            new_mgc["fiveprime"][readlen][new_pos] += count
                    for readlen in mgc["threeprime"]:
                        new_mgc["threeprime"][readlen] = {}
                        if readlen in offsets["threeprime"]["offsets"]:
                            offset = offsets["threeprime"]["offsets"][readlen]
                        else:
                            offset = -12
                        for pos in mgc["threeprime"][readlen]:
                            count = mgc["threeprime"][readlen][pos]
                            new_pos = pos+offset
                            if new_pos not in new_mgc["threeprime"][readlen]:
                                new_mgc["threeprime"][readlen][new_pos] = 0
                            new_mgc["threeprime"][readlen][new_pos] += count
                    mgc = new_mgc
                sqlite_db.close()

                for primetype in ["fiveprime", "threeprime"]:
                    for i in pos_list:
                        for readlen in range(minreadlen, maxreadlen+1):
                            if readlen in mgc[primetype]:
                                if i in mgc[primetype][readlen]:
                                    count_dict[primetype][i] += (
                                        mgc[primetype][readlen][i])

                if not metagene_aggregate:

                    fiveprime_counts[file_desc] = []
                    threeprime_counts[file_desc] = []
                    for i in pos_list:
                        fiveprime_counts[file_desc].append(
                            count_dict["fiveprime"][i])
                        threeprime_counts[file_desc].append(
                            count_dict["threeprime"][i])
                    # reset the count dict so that counts are not aggregated
                    for primetype in ["fiveprime", "threeprime"]:
                        for i in range(minpos, maxpos+1):
                            count_dict[primetype][i] = 0

        if metagene_normalise and not metagene_aggregate:
            min_mapped_reads = 1000000000
            for file_desc in mapped_reads_dict:
                if mapped_reads_dict[file_desc] < min_mapped_reads:
                    min_mapped_reads = mapped_reads_dict[file_desc]
            for file_desc in mapped_reads_dict:
                try:
                    factor = float(min_mapped_reads /
                                   mapped_reads_dict[file_desc])
                except:
                    return ("Error, missing mapped reads value for one of the files so cannot normalize")
                mapped_reads_dict[file_desc] = factor
            for file_desc in fiveprime_counts:
                norm_counts = []
                for count in fiveprime_counts[file_desc]:
                    factor = mapped_reads_dict[file_desc]
                    normalised_count = count*factor
                    norm_counts.append(normalised_count)
                fiveprime_counts[file_desc] = norm_counts
            for file_desc in threeprime_counts:
                norm_counts = []
                for count in threeprime_counts[file_desc]:
                    norm_counts.append(count*mapped_reads_dict[file_desc])
                threeprime_counts[file_desc] = norm_counts

        if metagene_aggregate:
            for i in pos_list:
                fiveprime_counts.append(count_dict["fiveprime"][i])
                threeprime_counts.append(count_dict["threeprime"][i])
        title = "Metagene profile"

        return metainfo_plots.metagene_plot(pos_list, fiveprime_counts, threeprime_counts, data["metagene_type"], title, minpos, maxpos, short_code, background_col, metagene_fiveprime_col, metagene_threeprime_col, title_size, axis_label_size, subheading_size, marker_size, data["metagene_end"], metagene_aggregate)

    elif data["plottype"] == "trip_periodicity":
        read_dict = {"readlengths": [],
                     "frame1": [],
                     "frame2": [],
                     "frame3": []}
        if trip_maxreadlen < trip_minreadlen:
            connection.close()
            return ("Error: max read length less than min read length, increase max read length using the input at the top of the page.")
        for i in range(trip_minreadlen, trip_maxreadlen+1):
            read_dict["readlengths"].append(i)
            read_dict["frame1"].append(0)
            read_dict["frame2"].append(0)
            read_dict["frame3"].append(0)
        for filetype in file_paths_dict:
            for file_id in file_paths_dict[filetype]:
                filepath = file_paths_dict[filetype][file_id]
                if os.path.isfile(filepath):
                    sqlite_db = SqliteDict(filepath, autocommit=False)
                else:
                    return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath))
                if "trip_periodicity" not in sqlite_db or redo_periodicity:

                    redo_periodicity_plots(transhelve, filepath)
                    # return "No triplet periodicity data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
                trip_periodicity_dict = sqlite_db["trip_periodicity"]
                sqlite_db.close()
                readlen_index = 0
                for readlength in read_dict["readlengths"]:
                    if readlength in trip_periodicity_dict["fiveprime"]:
                        read_dict["frame1"][readlen_index] += trip_periodicity_dict["fiveprime"][readlength]["0"]
                        read_dict["frame2"][readlen_index] += trip_periodicity_dict["fiveprime"][readlength]["1"]
                        read_dict["frame3"][readlen_index] += trip_periodicity_dict["fiveprime"][readlength]["2"]
                    readlen_index += 1
        title = "Triplet periodicity"
        # logging.warn("Location is {}".format(url_for('taskstatus',task_id=task.id)) )
        return metainfo_plots.trip_periodicity_plot(read_dict, title, short_code, background_col, title_size, axis_label_size, subheading_size, marker_size, legend_size)

    elif data["plottype"] == "mapped_reads_plot":
        labels = [""]
        unmapped = [0]
        mapped_coding = [0]
        mapped_noncoding = [0]
        ambiguous = [0]
        cutadapt_removed = [0]
        rrna_removed = [0]
        pcr_duplicates = [0]

        for filetype in file_paths_dict:
            for file_id in file_paths_dict[filetype]:
                filepath = file_paths_dict[filetype][file_id]

                cursor.execute(
                    "SELECT file_name,file_description from files where file_id = '{}';".format(file_id))
                result = cursor.fetchone()
                file_name = (result[0]).replace(".shelf", "")
                labels.append(result[1])

                if os.path.isfile(filepath):
                    sqlite_db = SqliteDict(filepath, autocommit=False)
                else:
                    return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath))

                if "unmapped_reads" in sqlite_db:
                    unmapped.append(sqlite_db["unmapped_reads"])
                else:
                    unmapped.append(0)

                if "coding_counts" in sqlite_db:
                    mapped_coding.append(sqlite_db["coding_counts"])
                else:
                    mapped_coding.append(0)

                if "noncoding_counts" in sqlite_db:
                    mapped_noncoding.append(sqlite_db["noncoding_counts"])
                else:
                    mapped_noncoding.append(0)

                if "ambiguous_counts" in sqlite_db:
                    ambiguous.append(sqlite_db["ambiguous_counts"])
                else:
                    ambiguous.append(0)

                if "cutadapt_removed" in sqlite_db:
                    if sqlite_db["cutadapt_removed"] != "NULL":
                        cutadapt_removed.append(sqlite_db["cutadapt_removed"])
                    else:
                        cutadapt_removed.append(0)
                else:
                    cutadapt_removed.append(0)

                if "rrna_removed" in sqlite_db:
                    rrna_removed.append(sqlite_db["rrna_removed"])
                else:
                    rrna_removed.append(0)

                if "pcr_duplicates" in sqlite_db:
                    pcr_duplicates.append(sqlite_db["pcr_duplicates"])
                else:
                    pcr_duplicates.append(0)
                sqlite_db.close()

        labels.append("")

        # Append a 0 to the end of every list so that there will be an empty space on the plot at the right hand side
        for listname in [unmapped, mapped_coding, mapped_noncoding, ambiguous, cutadapt_removed, rrna_removed, pcr_duplicates]:
            listname.append(0)

        return metainfo_plots.mapped_reads_plot(unmapped, mapped_coding, mapped_noncoding, labels, ambiguous, cutadapt_removed, rrna_removed, short_code, background_col, title_size, axis_label_size, subheading_size, marker_size, breakdown_per, pcr_duplicates, legend_size)

    elif data["plottype"] == "heatmap":

        min_readlen = heatmap_minreadlen
        max_readlen = heatmap_maxreadlen
        min_pos = heatmap_startpos
        max_pos = heatmap_endpos
        master_count_list = []

        for filetype in file_paths_dict:
            for file_id in file_paths_dict[filetype]:
                positions = []
                readlengths = []
                count_list = []
                filepath = file_paths_dict[filetype][file_id]
                if os.path.isfile(filepath):
                    sqlite_db = SqliteDict(filepath, autocommit=False)
                else:
                    return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath))
                if "metagene_counts" not in sqlite_db:
                    return ("No metagene counts data for this file, please report this to tripsvizsite@gmail.com or via the contact page.")
                if data["heatmap_metagene_type"] == "metagene_start":
                    if data["metagene_list"] == "":
                        mgc = sqlite_db["metagene_counts"]
                    else:
                        mgc = create_custom_metagene("AUG", 0, 0, 3, 0, "cds", False, False, True, False,
                                                     sqlite_db, organism, data["metagene_list"], "All", transcriptome, transhelve)
                elif data["heatmap_metagene_type"] == "metagene_stop":
                    if data["metagene_list"] == "":
                        mgc = sqlite_db["stop_metagene_counts"]
                    else:
                        mgc = create_custom_metagene("UAG,UAA,UGA", 0, 0, 0, 3, "cds", False, False, False,
                                                     True, sqlite_db, organism, data["metagene_list"], "All", transcriptome, transhelve)
                elif data["heatmap_metagene_type"] == "metagene_second_aug":
                    mgc = sqlite_db["secondary_metagene_counts"]
                sqlite_db.close()

                for primetype in [data["heatmap_direction"]]:
                    for readlen in range(max_readlen, min_readlen-1, -1):
                        for i in range(min_pos, max_pos+1):
                            if readlen in mgc[primetype]:
                                if i in mgc[primetype][readlen]:
                                    if mgc[primetype][readlen][i] != 0:
                                        count_list.append(
                                            mgc[primetype][readlen][i])
                                    else:
                                        count_list.append(0)
                                else:
                                    count_list.append(0)
                                readlengths.append(readlen)
                                positions.append(i)
                            else:
                                readlengths.append(readlen)
                                positions.append(i)
                                count_list.append(0)
                if master_count_list == []:
                    master_count_list = count_list
                else:
                    for i in range(0, len(count_list)):
                        master_count_list[i] = master_count_list[i] + \
                            count_list[i]

        title = "Heatmap"
        # Replace 0 values with None in master_count_list and apply log transformation if applicable
        fixed_master_count_list = []
        for i in range(0, len(master_count_list)):
            if master_count_list[i] == 0:
                fixed_master_count_list.append(None)
            else:
                if log_scale:
                    fixed_master_count_list.append(
                        log(master_count_list[i], 2))
                else:
                    fixed_master_count_list.append(master_count_list[i])

        return metainfo_plots.heatplot(min_readlen, max_readlen, min_pos, max_pos, positions, readlengths, fixed_master_count_list, data["heatmap_metagene_type"], title, reverse_scale, color_palette, short_code, background_col, maxscaleval, str(title_size)+"pt", str(axis_label_size)+"pt", str(subheading_size)+"pt", str(marker_size)+"pt")

    else:
        plt_type = data["plottype"].strip(" ").replace("\n", "")
        if plt_type != "replicate_comp":
            print("unknown plot type", data["plottype"])
        return ("Error, unknown plot type selected: {}".format(data["plottype"]))


# Groups together counts from different filepaths for the metainformation counts table

def aggregate_counts(
    file_paths_dict: Dict[str, List[str]],
    traninfo_dict: Dict[str, List[str]],
    longest_tran_list: List[str],
    region: str,
    organism: str,
    all_seq_types: List[str],
    te_minimum_reads: float,
    html_args: Dict[str, str],
    count_type: str,
    te_tranlist: str,
) -> str:
    # Remove seq types that aren't used
    seq_list = []
    for seq_type in all_seq_types:
        if seq_type in file_paths_dict and len(file_paths_dict[seq_type]):
            seq_list.append(seq_type)
    all_seq_types = seq_list

    file_list = ""
    transcript_dict = {}
    table_str = ""
    header_str = (
        "<thead><tr><th>Filename</th> <th>Gene</th> <th>Transcript</th> <th>Region</th>"
    )
    filename = organism + "_translation_efficiencies_" + str(
        time.time()) + ".csv"
    table_str += filename + "?~"
    mapped_reads = 0.0
    if "riboseq" in all_seq_types:
        header_str += "<th>Ribo-Seq</th>"
    if "rnaseq" in all_seq_types:
        header_str += "<th>RNA-Seq</th>"
    if "riboseq" in all_seq_types and "rnaseq" in all_seq_types:
        header_str += "<th>Translation Efficency</th>"
    for seq_type in file_paths_dict:
        if seq_type != "riboseq" and seq_type != "rnaseq":
            if len(file_paths_dict[seq_type]) != 0:
                header_str += "<th>{}</th>".format(seq_type)
        for file_id in file_paths_dict[seq_type]:
            file_list += "{},".format(file_id)
            filepath = file_paths_dict[seq_type][file_id]
            if os.path.isfile(filepath):
                sqlite_db = SqliteDict(f"{filepath}",
                                       autocommit=False,
                                       decode=my_decoder)
            else:
                return ("File not found, please report this to"
                        " tripsvizsite@gmail.com or via the contact page.")
            try:
                mapped_reads += sqlite_db["coding_counts"] + \
                    sqlite_db["noncoding_counts"]
            except Exception:
                pass
            opendict = sqlite_db[f"unambiguous_{region}_totals"]
            sqlite_db.close()
            # print"opendict", opendict
            for transcript in longest_tran_list:
                if transcript not in transcript_dict:
                    transcript_dict[transcript] = {}
                if seq_type not in transcript_dict[transcript]:
                    transcript_dict[transcript][seq_type] = 0
                if transcript in opendict:
                    transcript_dict[transcript][seq_type] += opendict[
                        transcript]

    # If count_type is rpkm or tpm, transform the raw counts
    for transcript in longest_tran_list:
        tranlength = float(traninfo_dict[transcript]["length"])
        if traninfo_dict[transcript]["cds_start"]:
            if region == "fiveprime":
                tranlength = float(traninfo_dict[transcript]["cds_start"] -
                                   1)
            elif region == "cds":
                tranlength = float(
                    (traninfo_dict[transcript]["cds_stop"] + 3) -
                    (traninfo_dict[transcript]["cds_start"] - 1))
            else:
                tranlength = float(traninfo_dict[transcript]["length"] -
                                   (traninfo_dict[transcript]["cds_stop"] +
                                   3))
        tranlength /= 1000.
    if count_type == "rpkm":
        if not mapped_reads:
            return "Cannot calculate RPKM. No mapped reads data available"
        # per million scaling factor
        pmsf = mapped_reads / 1000000.

        # get tranlength in kilobases
        # print"AGGREGATE,pmsf, tranlength",pmsf, tranlength
        for seq_type in transcript_dict[transcript]:
            # print"COUNT ", transcript_dict[transcript][seq_type]
            transcript_dict[transcript][seq_type] = round(
                (transcript_dict[transcript][seq_type] / pmsf) /
                tranlength, 2)
    # If count_type is rpkm or tpm, transform the raw counts
    if count_type == "tpm":
        total_rpk = 0.0

        for transcript in longest_tran_list:
            tranlength = float(traninfo_dict[transcript]["length"])
            if traninfo_dict[transcript]["cds_start"]:
                if region == "fiveprime":
                    tranlength = float(traninfo_dict[transcript]["cds_start"] -
                                       1)
                elif region == "cds":
                    tranlength = float(
                        (traninfo_dict[transcript]["cds_stop"] + 3) -
                        (traninfo_dict[transcript]["cds_start"] - 1))
                else:
                    tranlength = float(traninfo_dict[transcript]["length"] -
                                       traninfo_dict[transcript]["cds_stop"] +
                                       3)
            # get tranlength in kilobases
            tranlength /= 1000
            for seq_type in transcript_dict[transcript]:
                rpk = transcript_dict[transcript][seq_type] / tranlength
                transcript_dict[transcript][seq_type] = rpk
                total_rpk += rpk
        pmsf = total_rpk / 1000000
        for transcript in transcript_dict:
            for seq_type in transcript_dict[transcript]:
                transcript_dict[transcript][seq_type] = round(
                    transcript_dict[transcript][seq_type] / pmsf, 2)

    total_rows = 0
    tmp_te_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename),
                       "w")
    tmp_te_file.write(
        "Filename,Gene,Transcript,Region,Riboseq count,Rnaseq count,"
        "Translation efficiency")
    for seq_type in all_seq_types:
        if seq_type not in ["riboseq", "rnaseq"]:
            tmp_te_file.write(",{}".format(seq_type))
    tmp_te_file.write("\n")
    all_rows = []
    if te_tranlist:
        final_tranlist = te_tranlist
    else:
        final_tranlist = transcript_dict.keys()
    for transcript in final_tranlist:
        seq_count_dict = {}
        for seq_type in transcript_dict[transcript]:
            try:
                gene = ((traninfo_dict[transcript]["gene"]).replace(
                    ",", "_").replace(";", "_"))
            except Exception:
                gene = "Unknown"
            if seq_type not in seq_count_dict:
                seq_count_dict[seq_type] = transcript_dict[transcript][
                    seq_type]
        if "riboseq" in seq_count_dict:
            riboseq_count = seq_count_dict["riboseq"]
        else:
            riboseq_count = 0
        if "rnaseq" in seq_count_dict:
            rnaseq_count = seq_count_dict["rnaseq"]
        else:
            rnaseq_count = 0
        if rnaseq_count < te_minimum_reads or riboseq_count < te_minimum_reads:
            continue
        if rnaseq_count == 0 or riboseq_count == 0:
            te = 0
        else:
            te = float(transcript_dict[transcript]["riboseq"]) / float(
                transcript_dict[transcript]["rnaseq"])
            # te = round(te,2)
        tmp_te_file.write("Aggregate,{},{},{},{},{},{}".format(
            gene, transcript, region, riboseq_count, rnaseq_count, te))
        input_list = [
            "Aggregate",
            gene,
            transcript,
            region,
        ]  # , riboseq_count, rnaseq_count, te]
        if "riboseq" in all_seq_types:
            input_list.append(riboseq_count)
        if "rnaseq" in all_seq_types:
            input_list.append(rnaseq_count)
        if "riboseq" in all_seq_types and "rnaseq" in all_seq_types:
            input_list.append(te)
        for seq_type in all_seq_types:
            if seq_type != "riboseq" and seq_type != "rnaseq":
                if seq_type in seq_count_dict:
                    input_list.append(seq_count_dict[seq_type])
                    tmp_te_file.write(",{}".format(seq_count_dict[seq_type]))
                else:
                    if seq_type in file_paths_dict:
                        input_list.append(0)
                        tmp_te_file.write(",0")
        input_list.append("<a href='http://trips.ucc.ie/" + organism + "/" +
                          html_args["transcriptome"] +
                          "/interactive_plot/?tran=" + transcript + "&files=" +
                          file_list + "' target='_blank_' >View plot</a>")
        all_rows.append(input_list)
        tmp_te_file.write("\n")
    tmp_te_file.close()
    os.chmod("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename), 0o660)
    # if both rnaseq and riboseq files, sort by te, else sort by the relevant count
    anyfile = False
    for seq_type in set(all_seq_types) & set(file_paths_dict):
        if len(file_paths_dict[seq_type]):
            anyfile = True
            break
    if not anyfile:
        return (
            "No files selected. Select a file by clicking on a study name in the"
            " studies section. Then select at least one of the files that"
            " appear in the files section.")

    both = len(file_paths_dict["riboseq"]) and len(file_paths_dict["rnaseq"])
    all_sorted_rows = sorted(
        all_rows, key=lambda x: x[6 if both else 4], reverse=True)
    for row in all_sorted_rows:
        total_rows += 1
        if total_rows <= 1000:
            input_str = ""
            for item in row:
                input_str += "{}.;".format(item)
            input_str += "?~"
            table_str += input_str
    header_str += "<th>View Plot</th></tr></thead>"
    table_str = "TE?~" + str(total_rows) + "?~" + header_str + "?~" + table_str
    return table_str


def sample_counts(
    file_paths_dict: Dict[str, List[str]],
    traninfo_dict: Dict[str, Dict[str, str]],
    longest_tran_list: List[str],
    region: str,
    organism: str,
    all_seq_types: List[str],
    te_minimum_reads: int,
    html_args: Dict[str, str],
    count_type: str,
    te_tranlist: List[str],
) -> str:
    # Remove seq types that aren't used
    del_list = []
    for seq_type in all_seq_types:
        if seq_type not in ["riboseq", "rnaseq"]:
            if seq_type in file_paths_dict:
                if not len(file_paths_dict[seq_type]):
                    del_list.append(seq_type)
    for seq_type in del_list:
        all_seq_types.remove(seq_type)
    file_list = ""
    transcript_dict = {}
    table_str = ""
    filename = organism + "_translation_efficiencies_" + str(
        time.time()) + ".csv"
    table_str += filename + "?~"
    mapped_reads_dict = {}
    for seq_type in file_paths_dict:
        for file_id in file_paths_dict[seq_type]:
            file_list += "{},".format(file_id)
            filepath = file_paths_dict[seq_type][file_id]
            inputfilename = (filepath.split("/")[-1]).replace(".sqlite", "")
            if inputfilename not in mapped_reads_dict:
                mapped_reads_dict[inputfilename] = 0

            if os.path.isfile(filepath):
                sqlite_db = SqliteDict(f"{filepath}",
                                       autocommit=False,
                                       decode=my_decoder)
            else:
                return ("File not found, please report this to"
                        " tripsvizsite@gmail.com or via the contact page.")
            try:
                mapped_reads_dict[inputfilename] += sqlite_db["coding_counts"] + sqlite_db[
                    "noncoding_counts"]
            except Exception:
                pass
            if region in ["all", "cds", "fiveprime", "threeprime"]:
                try:
                    opendict = sqlite_db[f"unambiguous_{region}_totals"]
                except Exception:
                    return "No unambigous totals availabel for file: {}".format(
                        "/".join(filepath.split("/")[-2:]))

            for transcript in longest_tran_list:
                if transcript not in transcript_dict:
                    transcript_dict[transcript] = {}
                if seq_type not in transcript_dict[transcript]:
                    transcript_dict[transcript][seq_type] = {}
                if transcript in opendict:
                    transcript_dict[transcript][seq_type][inputfilename] = {
                        "count": opendict[transcript],
                        "file_id": str(file_id),
                    }
                else:
                    transcript_dict[transcript][seq_type][inputfilename] = {
                        "count": 0,
                        "file_id": str(file_id),
                    }

    # If count_type is rpkm or tpm, transform the raw counts
    if count_type == "rpkm":  # TODO: This should be merged with next for loop
        for transcript in longest_tran_list:
            tranlength = float(traninfo_dict[transcript]["length"])
            if traninfo_dict[transcript]["cds_start"]:
                if region == "fiveprime":
                    tranlength = float(traninfo_dict[transcript]["cds_start"] -
                                       1)
                elif region == "cds":
                    tranlength = float(
                        (traninfo_dict[transcript]["cds_stop"] + 3) -
                        (traninfo_dict[transcript]["cds_start"] -
                         1))  # NOTE: Why Float?
                elif region == "threeprime":
                    tranlength = float(traninfo_dict[transcript]["length"] -
                                       (traninfo_dict[transcript]["cds_stop"] +
                                       3))
            # get tranlength in kilobases
            tranlength = tranlength / 1000
            for seq_type in transcript_dict[transcript]:
                for inputfilename in transcript_dict[transcript][seq_type]:
                    mapped_reads = float(mapped_reads_dict[inputfilename])
                    if mapped_reads == 0.0:
                        return "No mapped reads data for file {}".format(
                            inputfilename)
                    # per million scaling factor
                    pmsf = mapped_reads / 1000000
                    transcript_dict[transcript][seq_type][inputfilename][
                        "count"] = round(
                            (transcript_dict[transcript][seq_type]
                             [inputfilename]["count"] / pmsf) / tranlength,
                            2,
                    )
    if count_type == "tpm":
        total_rpk_dict = {}
        for transcript in longest_tran_list:
            tranlength = float(traninfo_dict[transcript]["length"])
            if traninfo_dict[transcript]["cds_start"]:
                if region == "fiveprime":
                    tranlength = float(traninfo_dict[transcript]["cds_start"] -
                                       1)
                elif region == "cds":
                    tranlength = float(
                        (traninfo_dict[transcript]["cds_stop"] + 3) -
                        (traninfo_dict[transcript]["cds_start"] - 1))
                elif region == "threeprime":
                    tranlength = float(traninfo_dict[transcript]["length"] -
                                       (traninfo_dict[transcript]["cds_stop"] +
                                       3))
            # get tranlength in kilobases
            tranlength /= 1000
            for seq_type in transcript_dict[transcript]:
                for inputfilename in transcript_dict[transcript][seq_type]:
                    if inputfilename not in total_rpk_dict:
                        total_rpk_dict[inputfilename] = 0
                    rpk = (transcript_dict[transcript][seq_type][inputfilename]
                           ["count"] / tranlength)
                    transcript_dict[transcript][seq_type][inputfilename][
                        "count"] = rpk
                    total_rpk_dict[inputfilename] += rpk

        for transcript in transcript_dict:
            for seq_type in transcript_dict[transcript]:
                for inputfilename in transcript_dict[transcript][seq_type]:
                    pmsf = total_rpk_dict[inputfilename] / 1000000
                    transcript_dict[transcript][seq_type][inputfilename][
                        "count"] = round(
                            transcript_dict[transcript][seq_type]
                            [inputfilename]["count"] / pmsf,
                            2,
                    )
    total_rows = 0
    tmp_te_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename),
                       "w")
    tmp_te_file.write("Gene,Transcript,Region,")
    # write out headers to tmp_te_file
    transcript = transcript_dict.keys()[0]
    for seq_type in transcript_dict[transcript]:
        for inputfilename in transcript_dict[transcript][seq_type]:
            # if inputfilename not in skipped_files:
            tmp_te_file.write("{}_{},".format(inputfilename, seq_type))
    tmp_te_file.write("\n")

    all_rows = []

    final_tranlist = te_tranlist if te_tranlist else transcript_dict.keys()
    for transcript in final_tranlist:
        try:
            gene = traninfo_dict[transcript]["gene"]
        except KeyError:
            gene = "Unknown"
        tmp_te_file.write("{},{},{},".format(gene, transcript, region))
        for seq_type in transcript_dict[transcript]:
            for inputfilename in transcript_dict[transcript][seq_type]:
                # if inputfilename not in skipped_files:
                count = transcript_dict[transcript][seq_type][inputfilename][
                    "count"]
                tmp_te_file.write("{},".format(count))
        tmp_te_file.write("\n")

        for seq_type in transcript_dict[transcript]:
            if seq_type not in file_paths_dict:
                continue
            for inputfilename in transcript_dict[transcript][seq_type]:
                # if inputfilename not in skipped_files:
                count = transcript_dict[transcript][seq_type][inputfilename][
                    "count"]
                file_id = transcript_dict[transcript][seq_type][inputfilename][
                    "file_id"]
                if count < te_minimum_reads:
                    continue

                riboseq_count = 0.001
                rnaseq_count = 0.001
                if seq_type == "riboseq":
                    riboseq_count = count
                elif seq_type == "rnaseq":
                    rnaseq_count = count if count else 0.001
                try:
                    te = riboseq_count / rnaseq_count
                except Exception:
                    te = 0
                # tmp_te_file.write("{},{},{},{},{},{},{}".format(inputfilename,gene,transcript,region, riboseq_count, rnaseq_count, te))
                input_list = [
                    inputfilename,
                    gene,
                    transcript,
                    region,
                    riboseq_count,
                    rnaseq_count,
                    te,
                ]
                if seq_type not in ["riboseq", "rnaseq"]:
                    input_list.append(count)
                input_list.append("<a href='/" + organism +
                                  "/" + html_args["transcriptome"] +
                                  "/interactive_plot/?tran=" + transcript +
                                  "&files=" + file_id +
                                  "' target='_blank_' >View plot</a>")
                all_rows.append(input_list)
                # tmp_te_file.write("\n")
    # tmp_te_file.close()
    os.chmod("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename), 0o660)
    # if both rnaseq and riboseq files, sort by te, else sort by the relevant count
    anyfile = False
    if file_paths_dict["riboseq"] and file_paths_dict["rnaseq"]:
        all_sorted_rows = sorted(all_rows, key=lambda x: x[6], reverse=True)
    elif file_paths_dict["riboseq"]:  # WARNING: how to you priotize riboseq?
        all_sorted_rows = sorted(all_rows, key=lambda x: x[4], reverse=True)
    elif file_paths_dict["rnaseq"]:
        all_sorted_rows = sorted(all_rows, key=lambda x: x[5], reverse=True)
    else:
        for seq_type in all_seq_types:
            if seq_type in file_paths_dict and file_paths_dict[seq_type]:
                anyfile = True
                all_sorted_rows = sorted(all_rows,
                                         key=lambda x: x[7],
                                         reverse=True)
                break
        if not anyfile:
            return ("No files selected, please report this to "
                    "tripsvizsite@gmail.com or via the contact page.")
    for row in all_sorted_rows:
        total_rows += 1
        if total_rows <= 1000:
            input_str = ""
            for item in row:
                input_str += "{}.;".format(item)
            input_str += "?~"
            table_str += input_str
    table_str = "TE?~" + str(total_rows) + "?~" + table_str
    return table_str
