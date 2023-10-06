from flask import Blueprint, render_template, request, make_response, Response
from flask import current_app as app
from typing import Text
import sqlite3
import os
import config
from core_functions import (fetch_studies, fetch_files, fetch_study_info,
                            fetch_file_paths, generate_short_code, fetch_user)
import riboflask
from flask_login import current_user
import logging
import json
from sqlqueries import sqlquery, get_table
try:
    from orfQuant import incl_OPM_run_orfQuant
    from tripsTPM import TPM
except Exception:
    pass

# This is the single transcript plot page, user chooses gene, files and other settings
single_transcript_plotpage_blueprint = Blueprint("interactiveplotpage",
                                                 __name__,
                                                 template_folder="templates")


@single_transcript_plotpage_blueprint.route(
    '/<organism>/<transcriptome>/single_transcript_plot/')
def interactiveplotpage(organism: str, transcriptome: str) -> Response | Text:

    template_dict = request.args.to_dict()
    organism_id, accepted_studies = fetch_studies(organism, transcriptome)
    template_dict['studies_and_files'] = fetch_files(accepted_studies)
    gwips = get_table("organisms")
    print(accepted_studies)
    gwips_info = gwips.loc[(gwips.organism_id == organism_id)
                           & (gwips.transcriptome_list == transcriptome), [
                               "gwips_clade", "gwips_organism",
                               "gwips_database", "default_transcript"
                           ]].iloc[0]

    template_dict['transcript'] = gwips_info['default_transcript']
    template_dict['gwips_info'] = gwips_info
    template_dict['studyinfo_dict'] = fetch_study_info(organism_id)
    template_dict['organism'] = organism
    template_dict['transcriptome'] = transcriptome

    user_hili_starts = []
    user_hili_stops = []
    try:
        for item in user_hili.split(","):
            user_hili_starts.append(int(item.split("_")[0]))
            user_hili_stops.append(int(item.split("_")[1]))
    except Exception:
        pass

    try:
        user_minread = int(user_minread)
        user_maxread = int(user_maxread)
    except Exception:
        user_minread = None
        user_maxread = None
    advanced = 'True'
    consent = request.cookies.get("cookieconsent_status")
    rendered_template = render_template('single_transcript_plot.html',
                                        template_dict=template_dict)
    if consent == "deny":
        rendered_template = make_response(rendered_template)
        for cookie_name in request.cookies:
            if cookie_name != "cookieconsent_status":
                rendered_template.delete_cookie(cookie_name)
    return rendered_template


# Creates and serves the plots for the single transcript plot page
single_transcript_query_blueprint = Blueprint("query",
                                              __name__,
                                              template_folder="templates")


@single_transcript_query_blueprint.route('/query', methods=['POST'])
def query() -> str:
    # global user_short_passed
    try:
        user = current_user.name
    except Exception:
        user = None
    data = request.get_json()
    print(data, 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
    file_list = []
    file_ids = []
    study_ids = []
    for key, value in data.items():
        if key.startswith('file') and value:
            file_list.append(value)
            file_ids.append(int(key.split('_')[-1]))
            study_ids.append(int(key.split('_')[2]))
    data["file_list"] = file_list
    data["file_ids"] = file_ids
    data["study_ids"] = study_ids
    print(file_ids)

    # TODO: set the file count limit on js side

    total_files = len(data["file_list"])
    if total_files > 1500:
        return "A maximum of 1500 files can be selected on this page, currently there are {} selected".format(
            total_files)

    file_paths_dict = fetch_file_paths(data)

    # user_short = data["user_short"]

    owner = get_table('organisms')
    print(data['transcript'])
    owner = owner.loc[(owner.organism_name == data["organism"]) &
                      (owner.transcriptome_list == data["transcriptome"]),
                      "owner"].values[0]
    print(owner)

    user, logged_in = fetch_user()

    if owner == 1:
        sql_path = "{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC,
                                                       config.ANNOTATION_DIR,
                                                       data["organism"],
                                                       data["transcriptome"])
        if not os.path.isfile(sql_path):
            return_str = "Cannot find annotation file {}.{}.sqlite".format(
                data["organism"], data["transcriptome"])
            if app.debug == True:
                return return_str, "NO_CELERY", {'Location': None}
            else:
                return jsonify({
                    'current': 400,
                    'total': 100,
                    'status': 'return_str',
                    'result': return_str
                }), 200, {
                    'Location': ""
                }
    else:
        sql_path = "{0}/transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
            config.UPLOADS_DIR, owner, data["organism"], data["transcriptome"])
    transcripts = sqlquery(sql_path, "transcripts")
    transcripts = transcripts[transcripts.transcript ==
                              data["transcript"]].iloc[0]
    inputtran = True

    if result:
        newtran = result[0]
    else:
        inputtran = False

    if not inputtran:
        cursor.execute(
            "SELECT * from transcripts WHERE gene = '{}'".format(tran))
        result = cursor.fetchall()

        if result != []:
            if len(result) == 1:
                tran = str(result[0][0])
            else:
                return_str = "TRANSCRIPTS"
                if user == "test":
                    return_str = "QUANT_TRANSCRIPTS"
                    if len(file_paths_dict["riboseq"].values()) > 0:
                        pre_orfQuant_res = incl_OPM_run_orfQuant(
                            tran, sqlite_path_organism,
                            file_paths_dict["riboseq"].values())
                        pre_TPM_Ribo = TPM(tran, sqlite_path_organism,
                                           file_paths_dict["riboseq"].values(),
                                           "ribo")

                        max_TPM_Ribo = max(pre_TPM_Ribo.values())
                        TPM_Ribo = {
                            transcript:
                            round((pre_TPM_Ribo[transcript] / max_TPM_Ribo) *
                                  100, 2)
                            for transcript in pre_TPM_Ribo
                        }

                        max_orf = max(pre_orfQuant_res.values())
                        orfQuant_res = {
                            transcript:
                            round(
                                (pre_orfQuant_res[transcript] / max_orf) * 100,
                                2)
                            for transcript in pre_orfQuant_res
                        }

                    else:
                        orfQuant_res = {
                            transcript[0]: "Null"
                            for transcript in result
                        }
                        TPM_Ribo = {
                            transcript[0]: "Null"
                            for transcript in result
                        }

                    if len(file_paths_dict["rnaseq"].values()) > 0:
                        pre_TPM_RNA = TPM(tran, sqlite_path_organism,
                                          file_paths_dict["rnaseq"].values(),
                                          "rna")
                        max_TPM_RNA = max(pre_TPM_RNA.values())
                        TPM_RNA = {
                            transcript:
                            round(
                                (pre_TPM_RNA[transcript] / max_TPM_RNA) * 100,
                                2)
                            for transcript in pre_TPM_RNA
                        }

                    else:
                        TPM_RNA = {
                            transcript[0]: "Null"
                            for transcript in result
                        }

                for transcript in result:
                    cursor.execute(
                        "SELECT length,cds_start,cds_stop,principal,version from transcripts WHERE transcript = '{}'"
                        .format(transcript[0]))
                    tran_result = cursor.fetchone()
                    tranlen = tran_result[0]
                    cds_start = tran_result[1]
                    cds_stop = tran_result[2]
                    if str(tran_result[3]) == "1":
                        principal = "principal"
                    else:
                        principal = ""
                    version = tran_result[4]
                    if cds_start == "NULL" or cds_start == None:
                        cdslen = "NULL"
                        threeutrlen = "NULL"
                    else:
                        cdslen = cds_stop - cds_start
                        threeutrlen = tranlen - cds_stop
                    if user == "test":
                        if transcript[0] in orfQuant_res:
                            OPM_coverage = orfQuant_res[transcript[0]]
                        else:
                            OPM_coverage = "NULL"

                        if transcript[0] in TPM_RNA:
                            RNA_coverage = TPM_RNA[transcript[0]]
                        else:
                            RNA_coverage = "NULL"

                        if transcript[0] in TPM_Ribo:
                            ribo_coverage = TPM_Ribo[transcript[0]]
                        else:
                            ribo_coverage = "NULL"

                        return_str += (":{},{},{},{},{},{},{},{},{}".format(
                            transcript[0], version, tranlen, cds_start, cdslen,
                            threeutrlen, OPM_coverage, ribo_coverage,
                            RNA_coverage))
                    else:
                        return_str += (":{},{},{},{},{},{},{}".format(
                            transcript[0], version, tranlen, cds_start, cdslen,
                            threeutrlen, principal))
                return return_str

        else:
            return_str = "ERROR! Could not find any gene or transcript corresponding to {}".format(
                tran)
            logging.debug(return_str)
            return return_str
    transhelve.close()

    lite = "y" if 'varlite' in data else "n"
    preprocess = True if 'preprocess' in data else False
    uga_diff = True if 'uga_diff' in data else False
    color_readlen_dist = True if 'color_readlen_dist' in data else False
    ribo_coverage = True if 'ribo_coverage' in data else False
    nucseq = True if 'nucseq' in data else False
    mismatches = True if 'mismatches' in data else False
    ambiguous = True if 'ambiguous' in data else False
    pcr = True if 'pcr' in data else False
    noisered = True if 'noisered' in data else False
    mismatch = True if 'mismatch' in data else False

    user_short_passed = False
    if data["user_short"] == "None" or user_short_passed == True:
        short_code = generate_short_code(data, organism, data["transcriptome"],
                                         "interactive_plot")
    else:
        short_code = data["user_short"]
        user_short_passed = True

    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    connection.text_factory = str
    cursor = connection.cursor()
    # Put any publicly available seq types (apart from riboseq and rnaseq) here
    seq_rules = {
        "proteomics": {
            "frame_breakdown": 1
        },
        "conservation": {
            "frame_breakdown": 1
        },
        "tcpseq": {
            "frame_breakdown": 0
        }
    }

    # get user_id
    if current_user.is_authenticated:
        user_name = current_user.name
        cursor.execute(
            "SELECT user_id from users WHERE username = '{}';".format(
                user_name))
        result = (cursor.fetchone())
        user_id = result[0]
        print("current user id is", user_id)
        # get a list of organism id's this user can access
        cursor.execute(
            "SELECT background_col,uga_col,uag_col,uaa_col,title_size,subheading_size,axis_label_size,marker_size,cds_marker_width,cds_marker_colour,legend_size,ribo_linewidth from user_settings WHERE user_id = '{}';"
            .format(user_id))
        result = (cursor.fetchone())
        background_col = result[0]
        uga_col = result[1]
        uag_col = result[2]
        uaa_col = result[3]
        title_size = result[4]
        subheading_size = result[5]
        axis_label_size = result[6]
        marker_size = result[7]
        cds_marker_size = result[8]
        cds_marker_colour = result[9]
        legend_size = result[10]
        ribo_linewidth = result[11]
        # get rules for all custom seq types
        cursor.execute(
            "SELECT * from seq_rules WHERE user_id = {};".format(user_id))
        result = (cursor.fetchall())
        for row in result:
            seq_name = row[1]
            frame_breakdown = row[2]
            seq_rules[seq_name] = {"frame_breakdown": frame_breakdown}
        connection.close()
    if tran != "":
        return riboflask.generate_plot(
            tran, ambiguous, minread, maxread, lite, ribocoverage, organism,
            readscore, noisered, primetype, minfiles, nucseq, user_hili_starts,
            user_hili_stops, uga_diff, file_paths_dict, short_code,
            color_readlen_dist, background_col, uga_col, uag_col, uaa_col,
            advanced, seqhili, seq_rules, title_size, subheading_size,
            axis_label_size, marker_size, transcriptome, config.UPLOADS_DIR,
            cds_marker_size, cds_marker_colour, legend_size, ribo_linewidth,
            secondary_readscore, pcr, mismatches, hili_start, hili_stop)

    else:
        return "ERROR! Could not find any transcript or gene corresponding to {}".format(
            tran)
