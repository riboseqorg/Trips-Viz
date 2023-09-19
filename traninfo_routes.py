from flask import Blueprint, render_template, request
import sqlite3
import os
import time
import config
import pandas as pd
from core_functions import (fetch_studies, fetch_files, fetch_study_info,
                            generate_short_code, fetch_user)
import traninfo_plots
from flask_login import current_user
import json
import fixed_values
from sqlqueries import get_table, sqlquery, get_user_id, table2dict

traninfo_plotpage_blueprint = Blueprint("traninfo_plotpage",
                                        __name__,
                                        template_folder="templates")


def calc_gc(seq):
    count_dict = {"A": 0.0, "T": 0.0, "G": 0.0, "C": 0.0, "N": 0.0}
    for char in seq:
        count_dict[char] += 1
    return count_dict


@traninfo_plotpage_blueprint.route('/<organism>/<transcriptome>/traninfo_plot/'
                                   )
def traninfo_plotpage(organism: str, transcriptome: str) -> str:
    # global user_short_passed
    global local
    try:
        print(local)
    except Exception:
        local = False
    print(organism, transcriptome)

    user = fetch_user()[0]
    accepted_studies = fetch_studies(user, organism, transcriptome)
    _, accepted_studies, accepted_files, seq_types = fetch_files(
        accepted_studies)
    organisms = get_table("organisms")

    result = organisms.loc[organism.organism_name == organism, [
        "gwips_clade", "gwips_organism", "gwips_database", "default_transcript"
    ]].iloc[0]

    studyinfo_dict = fetch_study_info(organism)
    gwips_clade = result[0]
    gwips_org = result[1]
    gwips_db = result[2]

    default_tran = result[3]
    # holds all values the user could possibly pass in the url
    # (keywords are after request.args.get), anything not passed
    # by user will be a string: "None"
    html_args = {
        "user_short": str(request.args.get('short')),
        "user_plot_type": str(request.args.get('plot')),
        "gc_tranlist": str(request.args.get('gc_tranlist')),
        "gc_tranlist2": str(request.args.get('gc_tranlist2')),
        "gc_tranlist3": str(request.args.get('gc_tranlist3')),
        "gc_tranlist4": str(request.args.get('gc_tranlist4')),
        "plot_type": str(request.args.get('plot_type')),
        "gc_location": str(request.args.get('gc_location')),
        "nucleotide": str(request.args.get('nucleotide')),
        "nuc_freq_plot_anchor": str(request.args.get('nuc_freq_plot_anchor')),
        "nuc_freq_plot_window": str(request.args.get('nuc_freq_plot_window')),
        "nuc_freq_plot_tranlist":
        str(request.args.get('nuc_freq_plot_tranlist')),
        "orftype": str(request.args.get('orftype')),
        "transcriptome": str(transcriptome),
        "metagene_tranlist": str(request.args.get('metagene_tranlist'))
    }
    return render_template('traninfo_index.html',
                           gwips_clade=gwips_clade,
                           gwips_org=gwips_org,
                           gwips_db=gwips_db,
                           transcriptome=transcriptome,
                           organism=organism,
                           default_tran=default_tran,
                           current_username=user,
                           local=local,
                           studies_dict=accepted_studies,
                           accepted_files=accepted_files,
                           html_args=html_args,
                           studyinfo_dict=studyinfo_dict,
                           seq_types=seq_types)


# Used to create custom metagene plots on the traninformation plot page

traninfoquery_blueprint = Blueprint("traninfoquery",
                                    __name__,
                                    template_folder="templates")


@traninfoquery_blueprint.route('/traninfoquery', methods=['POST'])
def traninfoquery() -> str:
    # global user_short_passed
    user_short_passed = True
    data = json.loads(request.data)
    plottype = data["plottype"]

    exclude_first_val = int(data["exclude_first_val"])
    exclude_last_val = int(data["exclude_last_val"])
    include_first_val = int(data["include_first_val"])
    include_last_val = int(data["include_last_val"])
    metagene_tranlist = (data["metagene_tranlist"].strip(" ")).replace(
        " ", ",")

    gc_tranlist = data['gc_tranlist'].upper()
    gc_tranlist2 = data['gc_tranlist2'].upper()
    gc_tranlist3 = data['gc_tranlist3'].upper()
    gc_tranlist4 = data['gc_tranlist4'].upper()
    nuc_freq_plot_tranlist = data['nuc_freq_plot_tranlist']
    nuc_freq_plot_window = int(data['nuc_freq_plot_window'])
    maxscaleval = data['maxscaleval']

    if maxscaleval != "None" and maxscaleval != "":
        try:
            maxscaleval = int(maxscaleval)
        except Exception:
            maxscaleval = "None"

    organism = data['organism']
    transcriptome = data['transcriptome']
    gc_location = data["gc_location"]
    nucleotide = data["nucleotide"]
    plot_type = data["plot_type"]
    nuc_freq_plot_anchor = data["nuc_freq_plot_anchor"]
    orftype = data["orftype"]
    html_args = data["html_args"]

    principal = "principal" in data
    exons = "exons" in data

    user_settings = config.DEFAULT_USER_SETTINGS.copy()

    # get a list of organism id's this user can access
    if current_user.is_authenticated:
        # get user_id
        user_name = current_user.name
        user_id = get_user_id(user_name)
        user_settings = get_table('user_settings')
        user_settings = user_settings[user_settings.user_id == user_id]
        user_settings = table2dict(user_settings,
                                   '')  # TODO: Correct for empty

    if html_args["user_short"] == "None" or user_short_passed:
        short_code = generate_short_code(data, organism,
                                         html_args["transcriptome"],
                                         "traninfo_plot")
    else:
        short_code = html_args["user_short"]
        user_short_passed = True
    if plottype == "fetch_seq":
        filename = organism + "_sequences_" + str(time.time()) + ".fa"
        table_str = filename + "?~"
        splitlist = (gc_tranlist.replace(" ", ",")).split(",")
        tmp_fa_file = open(
            "{}/static/tmp/{}".format(config.SCRIPT_LOC, filename), "w")
        organisms = get_table("organisms")
        owner = organisms.loc[organisms.organism_name == organism
                              & organisms.transcriptome_list == transcriptome,
                              "owner"].values[0]
        if owner:
            sqlite_path_organism = "{0}/{1}/{2}/{2}.{3}.sqlite".format(
                config.SCRIPT_LOC, config.ANNOTATION_DIR, organism,
                transcriptome)
            if not os.path.isfile(sqlite_path_organism):
                return_str = "Cannot find annotation file {}.{}.sqlite".format(
                    organism, transcriptome)
                return return_str
        else:
            sqlite_path_organism = ("{0}/transcriptomes/{1}/{2}/{3}/{2}"
                                    "_{3}.sqlite".format(
                                        config.UPLOADS_DIR, owner, organism,
                                        transcriptome))
        transhelve = sqlite3.connect(sqlite_path_organism)
        cursor = transhelve.cursor()
        transcripts = sqlquery(sqlite_path_organism, "transcripts")
        transcripts = transcripts.loc[
            transcripts.transcript.isin(splitlist),
            ["gene", "transcript", "sequence", "cds_start", "cds_stop"
             ]].drop_duplicates("transcript")
        transcript_without_start = transcripts.loc[
            pd.isna(transcripts["cds_start"]), "transcript"].value

        if gc_location != "all":
            if transcript_without_start:
                return "Error: {} is non-coding, region cannot be {}".format(
                    ", ".join(transcript_without_start), gc_location)
            transcripts["cds_start"] -= 1

            # gc_location,exclude_first_val,exclude_last_val,include_first_val,include_last_val,gc_tranlist
            # keep track of transcript that have been written to file to avoid duplicates
        if gc_location == "five":
            transcripts.seq = transcripts.seq.apply(lambda x: x[:cds_start])
        elif gc_location == "start":
            transcripts.seq = transcripts.seq.apply(lambda x: x[
                cds_start - include_first_val:cds_start + include_last_val])
        elif gc_location == "cds":
            transcripts.seq = transcripts.seq.apply(
                lambda x: x[cds_start:cds_end])
        elif gc_location == "stop":
            transcripts.seq = transcripts.seq.apply(
                lambda x: x[:cds_stop + include_last_val])
        else:  # gc_location == "three":
            transcripts.seq = transcripts.seq.apply(lambda x: x[cds_stop:])

        subseq = ""  # TODO: This need to be corrected here
        if sum([
                include_first_val, include_last_val, exclude_first_val,
                exclude_last_val
        ]) == 0 or (gc_location in ["start", "stop"]):
            subseq = seq
        else:
            if include_first_val != 0:
                subseq += seq[:include_first_val]
            if include_last_val != 0:
                subseq += seq[(include_last_val * -1):]
            if exclude_first_val != 0:
                seqlen = len(seq)
                subseq += seq[(seqlen - exclude_first_val) * -1:]
            if exclude_last_val != 0:
                seqlen = len(seq)
                subseq += seq[:(seqlen - exclude_last_val)]
        tmp_fa_file.write(">{}_{}\n{}\n".format(gene, tran, subseq))
        tmp_fa_file.close()
        total_rows = transcripts.shape[0]
        table_str = "FA?~" + str(total_rows) + "?~" + table_str

        return table_str
    if plottype == "nuc_freq_plot":
        splitlist = (nuc_freq_plot_tranlist.replace(" ", ",")).split(",")
        filename = "Sequences_{}.fa".format(time.time())
        outfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename),
                       "w")
        organisms = get_table("organisms")
        owner = organisms.loc[organisms.organism_name == organism
                              & organisms.transcriptome_list == transcriptome,
                              "owner"].values[0]

        if owner == 1:
            transhelve = ("{0}/{1}/{2}/{2}.{3}.sqlite".format(
                config.SCRIPT_LOC, config.ANNOTATION_DIR, organism,
                transcriptome))
        else:
            transhelve = (
                "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
                    config.UPLOADS_DIR, owner, organism, transcriptome))
        transcripts = sqlquery(transhelve, "transcrips")

        if splitlist:
            transcripts = transcripts[transcripts.transcript.isin(splitlist)]
        result = transcrpts[transcript, cds_start, cds_stop, sequence]
        master_dict = {}
        for i in range(nuc_freq_plot_window * -1, nuc_freq_plot_window):
            master_dict[i] = {"A": 0, "T": 0, "G": 0, "C": 0}
        for row in result:
            tran = row[0]
            try:
                cds_start = int(row[1])
                cds_stop = int(row[2])
            except Exception:
                cds_start = None
                cds_stop = None
            seq = row[3]
            seqlen = len(seq)
            if nuc_freq_plot_anchor == "tss":
                outfile.write(">{}\n{}\n".format(tran,
                                                 seq[0:nuc_freq_plot_window]))
                for i in range(0, nuc_freq_plot_window):
                    try:
                        master_dict[i][seq[i]] += 1
                    except Exception:
                        pass
            if nuc_freq_plot_anchor == "cds_start":
                if cds_start != None:
                    outfile.write(">{}\n{}\n".format(
                        tran,
                        seq[cds_start + (nuc_freq_plot_window * -1):cds_start +
                            nuc_freq_plot_window]))
                    for i in range(nuc_freq_plot_window * -1,
                                   nuc_freq_plot_window):
                        try:
                            seqpos = cds_start + i
                            if seqpos >= 0:
                                master_dict[i][seq[seqpos]] += 1
                        except Exception:
                            pass
            if nuc_freq_plot_anchor == "cds_stop":
                if cds_stop:
                    outfile.write(">{}\n{}\n".format(
                        tran,
                        seq[cds_stop + (nuc_freq_plot_window * -1):cds_stop +
                            nuc_freq_plot_window]))
                    for i in range(nuc_freq_plot_window * -1,
                                   nuc_freq_plot_window):
                        seqpos = cds_stop + i
                        try:
                            master_dict[i][seq[seqpos]] += 1
                        except Exception:
                            pass
            if nuc_freq_plot_anchor == "tts":
                outfile.write(">{}\n{}\n".format(
                    tran, seq[-1 * nuc_freq_plot_window:]))
                for i in range(nuc_freq_plot_window * -1, 0):
                    try:
                        master_dict[i][seq[seqlen + i]] += 1
                    except Exception:
                        pass

        title = nuc_freq_plot_anchor
        return traninfo_plots.nuc_freq_plot(master_dict, title, short_code,
                                            config.DEFAULT_USER_SETTING,
                                            filename)

    # Nucleotide composition (single transcript)
    if plottype == "nuc_comp_single":
        master_dict = {}
        metagene_tranlist = metagene_tranlist.split(",")
        if len(metagene_tranlist) == 1:
            tran = metagene_tranlist[0]
            organisms = get_table("organisms")
            owner = organisms.loc[
                organisms.organism_name == organism
                # NOTE: Something is not right here, please recheck
                & organisms.transcriptome_list == transcriptome,
                "owner"].values[0]

            if owner == 1:
                transhelve = ("{0}/{1}/{2}/{2}.{3}.sqlite".format(
                    config.SCRIPT_LOC, config.ANNOTATION_DIR, organism,
                    transcriptome))
            else:
                transhelve = (
                    "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
                        config.UPLOADS_DIR, owner, organism, transcriptome))
            result = sqlquery(transhelve, "transcript").iloc[0]
            traninfo = {
                "transcript": result[0],
                "gene": result[1],
                "length": result[2],
                "cds_start": result[3],
                "cds_stop": result[4],
                "seq": result[5],
                "strand": result[6],
                "stop_list": result[7].split(","),
                "start_list": result[8].split(","),
                "exon_junctions": result[9].split(","),
                "tran_type": result[10],
                "principal": result[11]
            }

            title = tran
            return traninfo_plots.nuc_comp_single(tran, master_dict, title,
                                                  short_code,
                                                  config.DEFAULT_USER_SETTING,
                                                  traninfo)
        else:
            organisms = get_table("organisms")
            owner = organisms.loc[
                organisms.organism_name == organism
                & organisms.transcriptome_list == transcriptome,
                "owner"].values[0]

            if owner == 1:
                transhelve = ("{0}/{1}/{2}/{2}.{3}.sqlite".format(
                    config.SCRIPT_LOC, config.ANNOTATION_DIR, organism,
                    transcriptome))
            else:
                transhelve = (
                    "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
                        config.UPLOADS_DIR, owner, organism, transcriptome))
            transcripts = sqlquery(transhelve, "transcipts")
            if metagene_tranlist:
                transcrips = transcripts[transcripts.transcript.isin(
                    metagene_tranlist)]
            else:
                transcrips = transcripts[transcripts.principle
                                         & transcripts.tran_type]
            traninfo = transcripts[[
                "transcript", "cds_start", "cds_stop", "sequence"
            ]].values.tolist()
            title = "GC metagene of {} genes".format(len(traninfo))
            return traninfo_plots.gc_metagene(title, short_code,
                                              config.DEFAULT_USER_SETTING,
                                              traninfo)

    elif plottype == "orfstats":
        filename = organism + "_orfstats_" + str(time.time()) + ".csv"
        table_str = filename + "?~"
        tmp_te_file = open(
            "{}/static/tmp/{}".format(config.SCRIPT_LOC, filename), "w")
        organisms = get_table("organisms")
        owner = organisms.loc[organisms.organism_name == organism
                              & organisms.transcriptome_list == transcriptome,
                              "owner"].values[0]

        if owner == 1:
            transhelve = "{0}/{1}/{2}/{2}.{3}.sqlite".format(
                config.SCRIPT_LOC, config.ANNOTATION_DIR, organism,
                transcriptome)
        else:
            transhelve = "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
                config.UPLOADS_DIR, owner, organism, transcriptome)
        tmp_te_file.write(
            "Gene,Tran,ORF Type,Start Codon,Length,Start position,"
            "Stop position,CDS coverage\n")
        transcripts = sqlquery(transhelve, "transcripts")
        if principal:
            transcripts = transcripts[transcripts.principal]
        orf_type = sqlquery(transhelve, orftype)
        transcripts = transcripts.merge(orf_type, on="transcript")[[
            "gene", "transcript", "start_codon", "length", "start", "stop",
            "cds_coverage"
        ]]

        input_list = []
        result = transcripts
        for row in result:  # TODO: All these part should be cleaned
            gene = row[0]
            tran = row[1]
            start_codon = row[2]
            length = row[3]
            start = row[4]
            stop = row[5]
            cds_cov = row[6]
            input_list.append([
                gene, tran, orftype, start_codon, length, start, stop, cds_cov
            ])
            tmp_te_file.write("{},{},{},{},{},{},{},{}\n".format(
                gene, tran, orftype, start_codon, length, start, stop,
                cds_cov))
        tmp_te_file.close()
        total_rows = 0
        all_sorted_rows = input_list
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

    # Nucleotide composition (multiple transcripts)
    elif plottype == "nuc_comp_multi":
        organisms = get_table("organisms")
        owner = organisms.loc[organisms.organism_name == organism
                              & organisms.transcriptome_list == transcriptome,
                              "owner"].values[0]

        if owner == 1:
            transhelve = "{0}/{1}/{2}/{2}.{3}.sqlite".format(
                config.SCRIPT_LOC, config.ANNOTATION_DIR, organism,
                transcriptome)
        else:
            transhelve = "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
                config.UPLOADS_DIR, owner, organism, transcriptome)
        filename = "{}_{}_nucleotide_comp_{}.csv".format(
            organism, gc_location, str(time.time()))
        outfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename),
                       "w")

        if plot_type == "scatter":
            master_dict = {
                1: {
                    "trans": [],
                    "lengths": []
                },
                2: {
                    "trans": [],
                    "lengths": []
                },
                3: {
                    "trans": [],
                    "lengths": []
                },
                4: {
                    "trans": [],
                    "lengths": []
                }
            }
            gc_dict = {}
            transcripts = sqlquery(transhelve, "transcripts")
            if gc_tranlist != "":
                gc_tranlist = gc_tranlist.upper()
                splitlist = (gc_tranlist.replace(" ", ",")).split(",")

                if gc_tranlist2 != "":
                    splitlist2 = (gc_tranlist2.replace(" ", ",")).split(",")
                    for item in splitlist2:
                        splitlist.append(item)
                if gc_tranlist3 != "":
                    splitlist3 = (gc_tranlist3.replace(" ", ",")).split(",")
                    for item in splitlist3:
                        splitlist.append(item)
                if gc_tranlist4 != "":
                    splitlist4 = (gc_tranlist4.replace(" ", ",")).split(",")
                    for item in splitlist4:
                        splitlist.append(item)
                strlist = str(splitlist).strip("[]")
                transcripts = transcripts[transcripts.transcript.isin(
                    strlist.upper())]

            else:
                if gc_location == "all":
                    transcripts = transcripts[transcripts.principal]
                else:
                    transcripts = transcripts[transcripts.principal
                                              & transcripts.tran_type]
            result = transcripts[[
                "transcript", "sequence", "cds_start", "cds_stop"
            ]]
            for row in result:
                tran = row[0]
                seq = row[1]
                if gc_location != "all":
                    cds_start = int(row[2])
                    cds_stop = int(row[3])
                if gc_location == "all":
                    nuc_dict = calc_gc(seq)
                elif gc_location == "five":
                    nuc_dict = calc_gc(seq[:cds_start - 1])
                elif gc_location == "cds":
                    nuc_dict = calc_gc(seq[cds_start - 1:cds_stop + 3])
                elif gc_location == "three":
                    nuc_dict = calc_gc(seq[cds_stop + 3:])
                tot = sum([
                    nuc_dict["G"], nuc_dict["C"], nuc_dict["T"], nuc_dict["A"]
                ])
                if tot > 0:
                    if nucleotide == "GC":
                        gc = ((nuc_dict["G"] + nuc_dict["C"]) / tot)
                        gc = gc * 100
                        gc_dict[tran] = gc
                    elif nucleotide == "A":
                        a = ((nuc_dict["A"]) / tot)
                        a = a * 100
                        gc_dict[tran] = a
                    elif nucleotide == "T":
                        t = ((nuc_dict["T"]) / tot)
                        t = t * 100
                        gc_dict[tran] = t
                    elif nucleotide == "G":
                        g = ((nuc_dict["G"]) / tot)
                        g = g * 100
                        gc_dict[tran] = g
                    elif nucleotide == "C":
                        c = ((nuc_dict["C"]) / tot)
                        c = c * 100
                        gc_dict[tran] = c
            if gc_tranlist == "":
                for tran in gc_dict:
                    master_dict[1]["trans"].append(tran)
                    master_dict[1]["lengths"].append(gc_dict[tran])
            else:
                splitlist = (gc_tranlist.replace(" ", ",")).split(",")
                for tran in splitlist:
                    if tran in gc_dict:
                        master_dict[1]["trans"].append(tran)
                        master_dict[1]["lengths"].append(gc_dict[tran])
                if gc_tranlist2 != "":
                    splitlist = (gc_tranlist2.replace(" ", ",")).split(",")
                    for tran in splitlist2:
                        if tran in gc_dict:
                            master_dict[2]["trans"].append(tran)
                            master_dict[2]["lengths"].append(gc_dict[tran])
                if gc_tranlist3 != "":
                    splitlist = (gc_tranlist3.replace(" ", ",")).split(",")
                    for tran in splitlist3:
                        if tran in gc_dict:
                            master_dict[3]["trans"].append(tran)
                            master_dict[3]["lengths"].append(gc_dict[tran])
                if gc_tranlist4 != "":
                    splitlist = (gc_tranlist4.replace(" ", ",")).split(",")
                    for tran in splitlist4:
                        if tran in gc_dict:
                            master_dict[4]["trans"].append(tran)
                            master_dict[4]["lengths"].append(gc_dict[tran])

            outfile.write("Group , GC\n")
            for group in master_dict:
                for gc in master_dict[group]["lengths"]:
                    outfile.write("{}.{}\n".format(group, gc))

            outfile.close()
            return traninfo_plots.nuc_comp_scatter(master_dict, filename,
                                                   str(title_size) + "pt",
                                                   str(axis_label_size) + "pt",
                                                   str(marker_size) + "pt",
                                                   nucleotide, short_code)
        elif plot_type == "box":
            transcripts = sqlquery(transhelve, "transcripts")
            master_dict = {
                1: {
                    "trans": [],
                    "gc": [],
                    "a": [],
                    "t": [],
                    "g": [],
                    "c": []
                },
                2: {
                    "trans": [],
                    "gc": [],
                    "a": [],
                    "t": [],
                    "g": [],
                    "c": []
                },
                3: {
                    "trans": [],
                    "gc": [],
                    "a": [],
                    "t": [],
                    "g": [],
                    "c": []
                },
                4: {
                    "trans": [],
                    "gc": [],
                    "a": [],
                    "t": [],
                    "g": [],
                    "c": []
                }
            }
            gc_dict = {}
            if gc_tranlist != "":
                splitlist = (gc_tranlist.replace(" ", ",")).split(",")

                if gc_tranlist2 != "":
                    splitlist2 = (gc_tranlist2.replace(" ", ",")).split(",")
                    for item in splitlist2:
                        splitlist.append(item)
                if gc_tranlist3 != "":
                    splitlist3 = (gc_tranlist3.replace(" ", ",")).split(",")
                    for item in splitlist3:
                        splitlist.append(item)
                if gc_tranlist4 != "":
                    splitlist4 = (gc_tranlist4.replace(" ", ",")).split(",")
                    for item in splitlist4:
                        splitlist.append(item)
                strlist = str(splitlist).strip("[]")
                transcripts = transcripts[transcripts.transcript.isin(
                    strlist.upper())]
            else:
                if gc_location == "all":
                    transcripts = transcripts[transcripts.principal]
                else:
                    transcripts = transcripts[transcripts.principal
                                              & transcripts.tran_type]
            result = transcripts[[
                "transcript", "seq", "cds_start", "cds_stop"
            ]]
            for row in result:
                tran = row[0]
                seq = row[1]
                if gc_location != "all":
                    try:
                        cds_start = int(row[2])
                        cds_stop = int(row[3])
                    except Exception:
                        return (
                            f"Error: {tran} is non coding, remove this from"
                            " the transcript list or change the region to 'All'"
                        )
                if gc_location == "all":
                    nuc_dict = calc_gc(seq)
                    # outfile.write(">{}\n{}\n".format(tran,seq))
                elif gc_location == "five":
                    nuc_dict = calc_gc(seq[:cds_start - 1])
                    outfile.write(">{}\n{}\n".format(tran,
                                                     seq[:cds_start - 1]))
                elif gc_location == "cds":
                    nuc_dict = calc_gc(seq[cds_start - 1:cds_stop + 3])
                elif gc_location == "three":
                    nuc_dict = calc_gc(seq[cds_stop + 3:])
                tot = sum([
                    nuc_dict["G"], nuc_dict["C"], nuc_dict["T"], nuc_dict["A"]
                ])
                if tot > 0:
                    gc = ((nuc_dict["G"] + nuc_dict["C"]) / tot)
                    a = ((nuc_dict["A"]) / tot)
                    t = ((nuc_dict["T"]) / tot)
                    g = ((nuc_dict["G"]) / tot)
                    c = ((nuc_dict["C"]) / tot)

                    gc = gc * 100
                    a = a * 100
                    t = t * 100
                    g = g * 100
                    c = c * 100

                    gc_dict[tran] = {"gc": gc, "a": a, "t": t, "g": g, "c": c}
            if gc_tranlist == "":
                for tran in gc_dict:
                    master_dict[1]["trans"].append(tran)
                    master_dict[1]["gc"].append(gc_dict[tran]["gc"])
                    master_dict[1]["a"].append(gc_dict[tran]["a"])
                    master_dict[1]["t"].append(gc_dict[tran]["t"])
                    master_dict[1]["g"].append(gc_dict[tran]["g"])
                    master_dict[1]["c"].append(gc_dict[tran]["c"])
            else:
                splitlist = (gc_tranlist.replace(" ", ",")).split(",")
                for tran in splitlist:
                    if tran in gc_dict:
                        master_dict[1]["trans"].append(tran)
                        master_dict[1]["gc"].append(gc_dict[tran]["gc"])
                        master_dict[1]["a"].append(gc_dict[tran]["a"])
                        master_dict[1]["t"].append(gc_dict[tran]["t"])
                        master_dict[1]["g"].append(gc_dict[tran]["g"])
                        master_dict[1]["c"].append(gc_dict[tran]["c"])
                if gc_tranlist2 != "":
                    splitlist = (gc_tranlist2.replace(" ", ",")).split(",")
                    for tran in splitlist2:
                        if tran in gc_dict:
                            master_dict[2]["trans"].append(tran)
                            master_dict[2]["gc"].append(gc_dict[tran]["gc"])
                            master_dict[2]["a"].append(gc_dict[tran]["a"])
                            master_dict[2]["t"].append(gc_dict[tran]["t"])
                            master_dict[2]["g"].append(gc_dict[tran]["g"])
                            master_dict[2]["c"].append(gc_dict[tran]["c"])
                if gc_tranlist3 != "":
                    splitlist = (gc_tranlist3.replace(" ", ",")).split(",")
                    for tran in splitlist3:
                        if tran in gc_dict:
                            master_dict[3]["trans"].append(tran)
                            master_dict[3]["gc"].append(gc_dict[tran]["gc"])
                            master_dict[3]["a"].append(gc_dict[tran]["a"])
                            master_dict[3]["t"].append(gc_dict[tran]["t"])
                            master_dict[3]["g"].append(gc_dict[tran]["g"])
                            master_dict[3]["c"].append(gc_dict[tran]["c"])
                if gc_tranlist4 != "":
                    splitlist = (gc_tranlist4.replace(" ", ",")).split(",")
                    for tran in splitlist4:
                        if tran in gc_dict:
                            master_dict[4]["trans"].append(tran)
                            master_dict[4]["gc"].append(gc_dict[tran]["gc"])
                            master_dict[4]["a"].append(gc_dict[tran]["a"])
                            master_dict[4]["t"].append(gc_dict[tran]["t"])
                            master_dict[4]["g"].append(gc_dict[tran]["g"])
                            master_dict[4]["c"].append(gc_dict[tran]["c"])
            outfile.write("Group ,GC\n")
            for group in master_dict:
                for gc in master_dict[group]["gc"]:
                    outfile.write("{},{}\n".format(group, gc))
            outfile.close()
            return traninfo_plots.nuc_comp_box(master_dict, filename,
                                               nucleotide,
                                               str(title_size) + "pt",
                                               config.A_COL,
                                               str(axis_label_size) + "pt",
                                               str(marker_size) + "pt",
                                               short_code)

    elif plottype == "codon_usage":
        aa_dict = fixed_values.codon_aa_full.copy()
        filename = organism + "_codon_usage_" + str(time.time()) + ".csv"
        cu_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename),
                       "w")
        cursor.execute(
            "SELECT owner FROM organisms WHERE organism_name = '{}' and"
            " transcriptome_list = '{}';".format(organism, transcriptome))
        owner = (cursor.fetchone())[0]
        if owner == 1:
            transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(
                config.SCRIPT_LOC, config.ANNOTATION_DIR, organism,
                transcriptome))
        else:
            transhelve = sqlite3.connect(
                "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
                    config.UPLOADS_DIR, owner, organism, transcriptome))

        traninfo_cursor = transhelve.cursor()
        codon_dict = {}
        if gc_tranlist != "":
            splitlist = (gc_tranlist.replace(" ", ",")).split(",")
            strlist = str(splitlist).strip("[]")
            traninfo_cursor.execute(
                "SELECT transcript,sequence,cds_start,cds_stop from"
                " transcripts WHERE transcript IN ({});".format(
                    strlist.upper()))
        else:
            traninfo_cursor.execute(
                "SELECT transcript,sequence,cds_start,cds_stop from"
                " transcripts WHERE principal = 1 and tran_type = 1;")
        result = traninfo_cursor.fetchall()
        total_trans = 0
        total_length = 0
        for row in result:
            tran = row[0]
            seq = row[1]
            total_trans += 1
            if not row[2] or not row[3]:
                return "Enter coding transcripts only, {} is noncoding".format(
                    row[0])
            start_pos = int(row[2]) - 1
            end_pos = int(row[3]) + 2
            total_length += (end_pos - start_pos) + 1
            for i in range(start_pos, end_pos, 3):
                try:
                    codon = seq[i:i + 3]
                except Exception:
                    pass
                if len(codon) != 3:
                    continue
                if codon not in codon_dict:
                    codon_dict[codon] = 0
                codon_dict[codon] += 1
        cu_file.write("Total_transcripts,{}\n".format(total_trans))
        cu_file.write("Total_length_(nts),{}\n".format(total_length))
        for codon in sorted(codon_dict.keys()):
            cu_file.write("{},{},{}\n".format(codon, aa_dict[codon],
                                              codon_dict[codon]))

        return fixed_values.codon_usage(codon_dict, short_code,
                                        str(title_size) + "pt",
                                        str(axis_label_size) + "pt",
                                        str(marker_size) + "pt", filename)
    # This is the lengths plot
    elif plottype == "lengths_plot":
        cursor.execute(
            "SELECT owner FROM organisms WHERE organism_name = '{}' and"
            " transcriptome_list = '{}';".format(organism, transcriptome))
        owner = (cursor.fetchone())[0]
        if owner == 1:
            transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(
                config.SCRIPT_LOC, config.ANNOTATION_DIR, organism,
                transcriptome))
        else:
            transhelve = sqlite3.connect(
                "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
                    config.UPLOADS_DIR, owner, organism, transcriptome))

        trancursor = transhelve.cursor()

        if plot_type == "box":
            master_dict = {
                1: {
                    "trans": [],
                    "lengths": []
                },
                2: {
                    "trans": [],
                    "lengths": []
                },
                3: {
                    "trans": [],
                    "lengths": []
                },
                4: {
                    "trans": [],
                    "lengths": []
                }
            }
            gc_dict = {}
            if gc_tranlist != "":
                splitlist = (gc_tranlist.replace(" ", ",")).split(",")

                if gc_tranlist2 != "":
                    splitlist2 = (gc_tranlist2.replace(" ", ",")).split(",")
                    for item in splitlist2:
                        splitlist.append(item)
                if gc_tranlist3 != "":
                    splitlist3 = (gc_tranlist3.replace(" ", ",")).split(",")
                    for item in splitlist3:
                        splitlist.append(item)
                if gc_tranlist4 != "":
                    splitlist4 = (gc_tranlist4.replace(" ", ",")).split(",")
                    for item in splitlist4:
                        splitlist.append(item)
                strlist = str(splitlist).strip("[]")
                trancursor.execute(
                    "SELECT transcript,length,cds_start,cds_stop,"
                    "exon_junctions from transcripts WHERE transcript IN ({});"
                    .format(strlist.upper()))
            else:
                if gc_location == "all":
                    trancursor.execute(
                        "SELECT transcript,length,cds_start,cds_stop,"
                        "exon_junctions from transcripts WHERE principal = 1;")
                else:
                    trancursor.execute(
                        "SELECT transcript,length,cds_start,cds_stop,"
                        "exon_junctions from transcripts WHERE principal = 1"
                        " and tran_type = 1;")
            result = trancursor.fetchall()
            if result == []:
                return "Could not find any info on given transcript list"
            for row in result:
                tran = row[0]
                tranlen = int(row[1])
                exon_junctions_raw = row[4].split(",")
                exon_junctions = []
                for item in exon_junctions_raw:
                    if item != '':
                        exon_junctions.append(int(item))
                if gc_location != "all":
                    cds_start = int(row[2])
                    cds_stop = int(row[3])

                if gc_location == "all":
                    length = tranlen
                    exon_no = length / (len(exon_junctions) + 1)
                elif gc_location == "five":
                    length = cds_start
                    local_exons = sum(i < cds_start for i in exon_junctions)
                    exon_no = length / (local_exons + 1)
                elif gc_location == "cds":
                    length = cds_stop - cds_start
                    local_exons = sum((i > cds_start and i < cds_stop)
                                      for i in exon_junctions)
                    exon_no = length / (local_exons + 1)
                elif gc_location == "three":
                    length = tranlen - cds_stop
                    local_exons = sum(i > cds_stop for i in exon_junctions)
                    exon_no = length / (local_exons + 1)
                if not exons:
                    gc_dict[tran] = length
                else:
                    gc_dict[tran] = exon_no
            filename = "Lengths_{}.csv".format(time.time())
            outfile = open(
                "{}/static/tmp/{}".format(config.SCRIPT_LOC, filename), "w")
            if gc_tranlist == "":
                for tran in gc_dict:
                    master_dict[1]["trans"].append(tran)
                    master_dict[1]["lengths"].append(gc_dict[tran])
            else:
                splitlist = (gc_tranlist.replace(" ", ",")).split(",")
                for tran in splitlist:
                    if tran in gc_dict:
                        master_dict[1]["trans"].append(tran)
                        master_dict[1]["lengths"].append(gc_dict[tran])
                        outfile.write("Group1,{},{}\n".format(
                            tran, gc_dict[tran]))
                    else:
                        return tran
                if gc_tranlist2 != "":
                    splitlist = (gc_tranlist2.replace(" ", ",")).split(",")
                    for tran in splitlist2:
                        if tran in gc_dict:
                            master_dict[2]["trans"].append(tran)
                            master_dict[2]["lengths"].append(gc_dict[tran])
                            outfile.write("Group2,{},{}\n".format(
                                tran, gc_dict[tran]))
                if gc_tranlist3 != "":
                    splitlist = (gc_tranlist3.replace(" ", ",")).split(",")
                    for tran in splitlist3:
                        if tran in gc_dict:
                            master_dict[3]["trans"].append(tran)
                            master_dict[3]["lengths"].append(gc_dict[tran])
                            outfile.write("Group3,{},{}\n".format(
                                tran, gc_dict[tran]))
                if gc_tranlist4 != "":
                    splitlist = (gc_tranlist4.replace(" ", ",")).split(",")
                    for tran in splitlist4:
                        if tran in gc_dict:
                            master_dict[4]["trans"].append(tran)
                            master_dict[4]["lengths"].append(gc_dict[tran])
                            outfile.write("Group4,{},{}\n".format(
                                tran, gc_dict[tran]))
            return traninfo_plots.lengths_box(master_dict, filename,
                                              config.A_COL, short_code,
                                              str(title_size) + "pt",
                                              str(marker_size) + "pt",
                                              str(axis_label_size) + "pt")
        elif plot_type == "scatter":
            master_dict = {
                1: {
                    "trans": [],
                    "lengths": []
                },
                2: {
                    "trans": [],
                    "lengths": []
                },
                3: {
                    "trans": [],
                    "lengths": []
                },
                4: {
                    "trans": [],
                    "lengths": []
                }
            }
            gc_dict = {}
            if gc_tranlist != "":
                splitlist = (gc_tranlist.replace(" ", ",")).split(",")

                if gc_tranlist2 != "":
                    splitlist2 = (gc_tranlist2.replace(" ", ",")).split(",")
                    for item in splitlist2:
                        splitlist.append(item)
                if gc_tranlist3 != "":
                    splitlist3 = (gc_tranlist3.replace(" ", ",")).split(",")
                    for item in splitlist3:
                        splitlist.append(item)
                if gc_tranlist4 != "":
                    splitlist4 = (gc_tranlist4.replace(" ", ",")).split(",")
                    for item in splitlist4:
                        splitlist.append(item)
                strlist = str(splitlist).strip("[]")
                trancursor.execute(
                    "SELECT transcript,length,cds_start,cds_stop from"
                    " transcripts WHERE transcript IN ({});".format(
                        strlist.upper()))
            else:
                if gc_location == "all":
                    trancursor.execute(
                        "SELECT transcript,length,cds_start,cds_stop from"
                        " transcripts WHERE principal = 1;")
                else:
                    trancursor.execute(
                        "SELECT transcript,length,cds_start,cds_stop from"
                        " transcripts WHERE principal = 1 and tran_type = 1;")
            result = trancursor.fetchall()
            for row in result:
                tran = row[0]
                tranlen = int(row[1])
                if gc_location != "all":
                    cds_start = int(row[2])
                    cds_stop = int(row[3])
                if gc_location == "all":
                    length = tranlen
                elif gc_location == "five":
                    length = cds_start
                elif gc_location == "cds":
                    length = cds_stop - cds_start
                elif gc_location == "three":
                    length = tranlen - cds_stop
                gc_dict[tran] = length
            filename = "Lengths_{}.csv".format(time.time())
            outfile = open(
                "{}/static/tmp/{}".format(config.SCRIPT_LOC, filename), "w")

            if gc_tranlist == "":
                for tran in gc_dict:
                    master_dict[1]["trans"].append(tran)
                    master_dict[1]["lengths"].append(gc_dict[tran])
            else:
                splitlist = (gc_tranlist.replace(" ", ",")).split(",")
                for tran in splitlist:
                    if tran in gc_dict:
                        master_dict[1]["trans"].append(tran)
                        master_dict[1]["lengths"].append(gc_dict[tran])
                if gc_tranlist2 != "":
                    splitlist = (gc_tranlist2.replace(" ", ",")).split(",")
                    for tran in splitlist2:
                        if tran in gc_dict:
                            master_dict[2]["trans"].append(tran)
                            master_dict[2]["lengths"].append(gc_dict[tran])
                if gc_tranlist3 != "":
                    splitlist = (gc_tranlist3.replace(" ", ",")).split(",")
                    for tran in splitlist3:
                        if tran in gc_dict:
                            master_dict[3]["trans"].append(tran)
                            master_dict[3]["lengths"].append(gc_dict[tran])
                if gc_tranlist4 != "":
                    splitlist = (gc_tranlist4.replace(" ", ",")).split(",")
                    for tran in splitlist4:
                        if tran in gc_dict:
                            master_dict[4]["trans"].append(tran)
                            master_dict[4]["lengths"].append(gc_dict[tran])
            outfile.close()
            return traninfo_plots.lengths_scatter(master_dict, filename,
                                                  str(title_size) + "pt",
                                                  str(axis_label_size) + "pt",
                                                  str(marker_size) + "pt",
                                                  short_code)

    elif plottype == "gene_count":
        organisms = get_table("organisms")
        owner = organisms.loc[organisms.organism_name == organism
                              & organisms.transcriptome_list == transcriptome,
                              "owner"].values[0]

        if owner == 1:
            transhelve = "{0}/{1}/{2}/{2}.{3}.sqlite".format(
                config.SCRIPT_LOC, config.ANNOTATION_DIR, organism,
                transcriptome)
        else:
            transhelve = "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
                config.UPLOADS_DIR, owner, organism, transcriptome)
        transcripts = sqlquery(transhelve, "transcripts")
        all_transcripts = transcripts.shape[0]
        coding_transcripts = transcripts[transcripts.tran_type].shape[0]
        all_genes = len(transcripts.gene.unique())
        coding_genes = len(
            transcripts[transcripts.gene_type == 1].gene.unique())
        coding = [1, coding_genes, coding_transcripts, 1]
        noncoding = [
            1, all_genes - coding_genes, all_transcripts - coding_transcripts,
            1
        ]
        return traninfo_plots.gene_count(short_code, background_col,
                                         title_size, axis_label_size,
                                         subheading_size, marker_size, coding,
                                         noncoding)

    else:
        if (plottype.strip(" ").replace("\n", "")) not in ["replicate_comp"]:
            print("Unknown plottype", plottype)
    return "Error, unknown plot type selected: {}".format(plottype)
