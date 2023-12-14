from flask import (
    Blueprint,
    render_template,
    request,
    make_response,
    Response,
    jsonify,
)
from flask import current_app as app
from typing import Text
from sqlqueries import get_user_id
import os
import config
from core_functions import (fetch_studies, fetch_files, fetch_study_info,
                            fetch_file_paths, generate_short_code, fetch_user)
import riboflask
from flask_login import current_user
import logging
import json
from sqlqueries import sqlquery, get_table
from orfQuant import incl_OPM_run_orfQuant
from tripsTPM import TPM

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
    # print(accepted_studies)
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
        for item in template_dict['user_hili'].split(","):
            user_hili_starts.append(int(item.split("_")[0]))
            user_hili_stops.append(int(item.split("_")[1]))
    except Exception:
        pass

    advanced = True
    consent = request.cookies.get("cookieconsent_status")
    rendered_template = render_template('single_transcript_plot.html',
                                        template_dict=template_dict)
    if consent == "deny":
        # TODO: Convert this according to django
        rendered_template = make_response(rendered_template)
        for cookie_name in request.cookies:
            if cookie_name != "cookieconsent_status":
                rendered_template.delete_cookie(cookie_name)
    return rendered_template


# Creates and serves the plots for the single transcript plot page
single_transcript_query_blueprint = Blueprint("query",
                                              __name__,
                                              template_folder="templates")


@single_transcript_query_blueprint.route('/query', methods=['POST',"GET"])
def query():  #TODO: add return type
    # global user_short_passed
    data = request.form.to_dict()
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

    file_paths_dict = fetch_file_paths(data)
    print(file_paths_dict, "Kiran")

    # user_short = data["user_short"]

    owner = get_table('organisms')
    owner = owner.loc[(owner.organism_name == data["organism"]) &
                      (owner.transcriptome_list == data["transcriptome"]),
                      "owner"].values[0]

    user = fetch_user()[0]

    if owner == 1:
        sql_path = "{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC,
                                                       config.ANNOTATION_DIR,
                                                       data["organism"],
                                                       data["transcriptome"])
        if not os.path.isfile(sql_path):
            return_str = "Cannot find annotation file {}.{}.sqlite".format(
                data["organism"], data["transcriptome"])
            if app.debug:
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
    transcripts_full = sqlquery(sql_path, "transcripts")
    transcripts = transcripts_full[
        (transcripts_full.transcript == data["transcript"].upper()) |
        (transcripts_full.gene == data["transcript"].upper())]

    if transcripts.empty:
        return_str = "ERROR! Could not find any gene or transcript corresponding to {}".format(
            data['transcript'])
        logging.debug(return_str)
        return return_str

    if data['transcript'].upper() not in transcripts.transcript.values:
        return_str = "TRANSCRIPTS"
        if user == "test":
            return_str = "QUANT_TRANSCRIPTS"
            try:  #  Riboseq
                pre_orfQuant_res = incl_OPM_run_orfQuant(
                    transcripts.transcript[0], sql_path,
                    file_paths_dict.loc[file_paths_dict["file_type"] ==
                                        "riboseq", "path"].values)
                pre_TPM_Ribo = TPM(
                    transcripts.transcript[0], sql_path,
                    file_paths_dict.loc[file_paths_dict["file_type"] ==
                                        "riboseq", "path"].values, "ribo")

                max_TPM_Ribo = max(pre_TPM_Ribo.values())
                TPM_Ribo = {
                    transcript:
                    round((pre_TPM_Ribo[transcript] * 100. / max_TPM_Ribo), 2)
                    for transcript in pre_TPM_Ribo
                }

                max_orf = max(pre_orfQuant_res.values())
                orfQuant_res = {
                    transcript:
                    round((pre_orfQuant_res[transcript] / max_orf) * 100, 2)
                    for transcript in pre_orfQuant_res
                }

            except KeyError:
                orfQuant_res = {
                    transcript: None
                    for transcript in transcripts.transcript
                }
                TPM_Ribo = orfQuant_res.copy()

            try:  #RNA Seq
                pre_TPM_RNA = TPM(
                    transcripts.transcript[0], sql_path,
                    file_paths_dict.loc[file_paths_dict["file_type"] ==
                                        "rnaseq", "path"].values, "rna")
                max_TPM_RNA = max(pre_TPM_RNA.values())
                TPM_RNA = {
                    transcript:
                    round((pre_TPM_RNA[transcript] / max_TPM_RNA) * 100, 2)
                    for transcript in pre_TPM_RNA
                }

            except KeyError:
                TPM_RNA = {
                    transcript: None
                    for transcript in transcripts.transcript
                }

        for _, transcript in transcripts.iterrows(
        ):  # TODO: Replace with iter tuple
            if not transcript.cds_start:
                cdslen = None
                three_utr_len = None
            else:
                cdslen = transcript.cds_stop - transcript.cds_start
                three_utr_len = transcript.length - transcript.cds_stop
            if user == "test":
                try:
                    OPM_coverage = orfQuant_res[transcript.length]
                except KeyError:
                    OPM_coverage = None
                try:
                    RNA_coverage = TPM_RNA[transcript.length]
                except KeyError:
                    RNA_coverage = None
                try:
                    ribo_coverage = TPM_Ribo[transcript.length]
                except KeyError:
                    ribo_coverage = None
                return_str += (":{},{},{},{},{},{},{},{},{}".format(
                    transcript.transcript, transcript.version,
                    transcript.length, transcript.cds_start, cdslen,
                    three_utr_len, OPM_coverage, ribo_coverage, RNA_coverage))

            else:
                return_str += (":{},{},{},{},{},{},{}".format(
                    transcript.transcript, transcript.version,
                    transcript.length, transcript.cds_start, cdslen,
                    three_utr_len, transcript.principal))
        print(return_str)
        return return_str

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
    settings = config.DEFAULT_USER_SETTINGS.copy()
    if current_user.is_authenticated:

        user_name = current_user.name
        user_id = get_user_id(user_name)
        user_settings = get_table('user_settings')
        user_settings = user_settings.loc[
            user_settings.user_id == user_id]  # TODO: Push to default settings
        for key in settings:
            settings[key] = user_settings[key]
        sequence_rule = get_table('seq_rules')
        sequence_rule = sequence_rule.loc[sequence_rule.user_id == user_id]
    return riboflask.generate_plot({
        'user_settings': settings,
        'seq_rules': seq_rules,
        'data': data
    })
