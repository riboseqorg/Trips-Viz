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
from sqlqueries_2 import get_user_id, sqlquery, get_table
import os
import polars as pl
import config
from core_functions import (fetch_studies, fetch_files, fetch_study_info,
                            fetch_file_paths, generate_short_code, fetch_user,
                            string2other, form_filler)
import riboflask
from flask_login import current_user
import logging
from orfQuant import incl_OPM_run_orfQuant
from tripsTPM import TPM

# This is the single transcript plot page, user chooses gene, files and other settings
single_transcript_plotpage_blueprint = Blueprint("interactiveplotpage",
                                                 __name__,
                                                 template_folder="templates")


@single_transcript_plotpage_blueprint.route(
    '/<organism>/<transcriptome>/single_transcript_plot/')
def interactiveplotpage(organism: str, transcriptome: str) -> Response | Text:
    """
    Single transcript plot page.

    Parameters:
    - organism (str): name of the organism
    - transcriptome (str): name of the transcript

    Returns:
    - html page

    Example:
    """

    # data = request.args.to_dict()
    data = form_filler(organism, transcriptome)
    organism_id, accepted_studies = fetch_studies(organism, transcriptome)
    # Accepted_studies is a DataFrame of (study_id, study_name)
    data['studies_and_files'] = fetch_files(accepted_studies).to_pandas()

    data['studyinfo_dict'] = fetch_study_info(organism_id)

    consent = request.cookies.get("cookieconsent_status")
    print(data)
    rendered_template = render_template('single_transcript_plot.html',
                                        template_dict=data)
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


@single_transcript_query_blueprint.route('/query', methods=['POST', "GET"])
def query():  #TODO: add return type
    """
    jquery route for single transcript plot.

    Parameters:
    - request

    Returns:
    """
    # global user_short_passed
    data = string2other(request.form.to_dict())
    data["transcript"] = data["transcript"].upper()
    print(data, "Anmol")
    data["file_ids"] = []

    for key in data:
        if key.startswith(f'file_{data["file_type"]}'):
            file_id = int(key.split('__')[-1])
            data["file_ids"].append(file_id)
    file_paths_dict = fetch_file_paths(data)
    print(file_paths_dict, "Kiran")

    # user_short = data["user_short"]

    owner = get_table('organisms').filter(
        (pl.col('organism_name') == data["organism"])
        & (pl.col('transcriptome_list') == data['transcriptome']))[0, 'owner']
    data['owner'] = owner

    print(owner, "Anmol")

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
    print(sql_path, "I am here")
    transcripts = sqlquery(
        sql_path,
        "transcripts").filter((pl.col('transcript') == data['transcript'])
                              | (pl.col('gene') == data['transcript']))

    if transcripts.is_empty():
        return_str = "ERROR! Could not find any gene or transcript corresponding to {}".format(
            data['transcript'])
        logging.debug(return_str)
        return return_str
    print(data['transcript'] in transcripts['transcript'], "yyyyyyyyyyyyy")

    if data['transcript'] not in transcripts['transcript']:
        return_str = "TRANSCRIPTS"
        if user == "test":
            return_str = "QUANT_TRANSCRIPTS"
            try:  #  Riboseq
                pre_orfQuant_res = incl_OPM_run_orfQuant(
                    transcripts[0, 'transcript'], sql_path,
                    file_paths_dict.loc[file_paths_dict["file_type"] ==
                                        "riboseq", "path"].values)
                pre_TPM_Ribo = TPM(
                    transcripts[0, 'transcript'], sql_path,
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
                    transcripts[0, 'transcript'], sql_path,
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
                    for transcript in transcripts["transcript"]
                }

        for transcript in transcripts.iter_rows(
                named=True):  # TODO: Replace with iter tuple
            if not transcript['cds_start']:
                cdslen = None
                three_utr_len = None
            else:
                cdslen = transcript['cds_stop'] - transcript['cds_start']
                three_utr_len = transcript['length'] - transcript['cds_stop']
            if user == "test":
                try:
                    OPM_coverage = orfQuant_res[transcript['length']]
                except KeyError:
                    OPM_coverage = None
                try:
                    RNA_coverage = TPM_RNA[transcript['length']]
                except KeyError:
                    RNA_coverage = None
                try:
                    ribo_coverage = TPM_Ribo[transcript['length']]
                except KeyError:
                    ribo_coverage = None
                return_str += (":{},{},{},{},{},{},{},{},{}".format(
                    transcript['transcript'], transcript['version'],
                    transcript['length'], transcript['cds_start'], cdslen,
                    three_utr_len, OPM_coverage, ribo_coverage, RNA_coverage))

            else:
                return_str += (":{},{},{},{},{},{},{}".format(
                    transcript['transcript'], transcript['version'],
                    transcript['length'], transcript['cds_start'], cdslen,
                    three_utr_len, transcript['principle']))
        print(return_str)
        return return_str
    # NOTE: Till here

    # get user_id
    settings = config.DEFAULT_USER_SETTINGS.copy()
    print(current_user.is_authenticated, "ZZZZZZZZZZZ")
    if current_user.is_authenticated:

        user_name = current_user.name
        user_id = get_user_id(user_name)
        user_settings = get_table('user_settings').filter(
            pl.col('user_id') == user_id)
        for key in settings:
            settings[key] = user_settings[0, key]
        data['sequence_rule'] = get_table('seq_rules').filter(
            pl.col('user_id') == user_id)
    # return ""
    data['file_paths_dict'] = file_paths_dict

    return riboflask.generate_plot(data, settings)
