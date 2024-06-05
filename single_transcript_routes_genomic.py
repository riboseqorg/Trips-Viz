from flask import Blueprint, render_template, request, make_response
from flask import current_app as app
import sqlite3
import os
import config
from core_functions import fetch_studies, fetch_files, fetch_study_info, fetch_file_paths, generate_short_code, fetch_user
import riboflask
from flask_login import current_user
import logging
import json
try:
    from orfQuant import incl_OPM_run_orfQuant
    from tripsTPM import TPM
except Exception:
    pass

#This is the single transcript plot page, user chooses gene, files and other settings
single_transcript_plotpage_genomic_blueprint = Blueprint(
    "interactiveplotpage_genomic", __name__, template_folder="templates")


@single_transcript_plotpage_genomic_blueprint.route(
    '/<organism>/<transcriptome>/interactive_plot_genomic/')
def interactiveplotpage_genomic(organism: str, transcriptome: str) -> str:
    """
    Interactive plot page.

    Parameters:
    - organism (str): name of the organism
    - transcriptome (str): name of the transcript

    Returns:
    - html page

    Example:

    """
    #global user_short_passed
    data = form_filler(organism, transcriptome)

    accepted_studies = fetch_studies(organism, transcriptome)
    file_id_to_name_dict, accepted_studies, accepted_files, seq_types = fetch_files(
        accepted_studies)

    if user_files != None:
        user_files = user_files.split(",")
        user_files = [str(x) for x in user_files]
    else:
        user_files = []

    user_ribo_studies = request.args.get('ribo_studies')
    if user_ribo_studies != None:
        user_ribo_studies = user_ribo_studies.split(",")
        user_ribo_studies = [str(x) for x in user_ribo_studies]
    else:
        user_ribo_studies = []
    user_proteomics_studies = request.args.get('proteomics_studies')
    if user_proteomics_studies != None:
        user_proteomics_studies = user_proteomics_studies.split(",")
        user_proteomics_studies = [str(x) for x in user_proteomics_studies]
    else:
        user_proteomics_studies = []

    user_rna_studies = request.args.get('rna_studies')
    if user_rna_studies != None:
        user_rna_studies = user_rna_studies.split(",")
        user_rna_studies = [str(x) for x in user_rna_studies]
    else:
        user_rna_studies = []

    if user_generate_shorturl == "F":
        user_generate_shorturl = False
    else:
        user_generate_shorturl = True

        user_maxread = None
    advanced = 'True'
    connection.close()
    consent = request.cookies.get("cookieconsent_status")
    if consent == "deny":
        resp = make_response(render_template('index.html', template_dict=data))
        for cookie_name in request.cookies:
            if cookie_name != "cookieconsent_status":
                resp.delete_cookie(cookie_name)
        return resp
    return render_template('index_genomic.html', template_dict=data)


# Creates and serves the plots for the single transcript plot page
single_transcript_query_genomic_blueprint = Blueprint(
    "query_genomic", __name__, template_folder="templates")


@single_transcript_query_genomic_blueprint.route('/query_genomic',
                                                 methods=['POST'])
def query_genomic() -> str:
    """

    """
    #global user_short_passed
    try:
        user = current_user.name
    except Exception:
        user = None
    #print "user", user
    data = request.data.to_dict()

    advanced = data["advanced"]
    logging.debug("FILE LIST")
    logging.debug(str(data["file_list"]))
    #logging.warn(len(data["file_list"]))
    #logging.debug("Length of alt file list is"+ len(data["alt_file_list"]))
    # Send file_list (a list of integers intentionally encoded as strings due to javascript), to be converted to a dictionary with riboseq/rnaseq lists of file paths.

    total_files = len(data["file_list"])
    if total_files > 1500:
        return "A maximum of 1500 files can be selected on this page, currently there are {} selected".format(
            total_files)

    file_paths_dict = fetch_file_paths(data["file_list"], organism)
    organisms = get_tables("organisms")
    owner = organisms.loc[(organisms.organism_name == organism) &
                          (organisms.transcriptome_list == transcriptome),
                          "owner"].values[0]

    user = fetch_user()[0]

    if owner == 1:
        sqlpath = "{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC,
                                                      config.ANNOTATION_DIR,
                                                      organism, transcriptome)
        if os.path.isfile(sqlpath):
            print("sqlite_path_organism", sqlpath)
        else:
            return_str = "Cannot find annotation file {}.{}.sqlite".format(
                organism, transcriptome)
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
        sqlpath = "{0}/transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
            config.UPLOADS_DIR, owner, organism, transcriptome)
        print("sqlite_path_organism", sqlpath)

    inputtran = True

    try:
        newtran = sqlquery(sqlpath,
                           "transcripts").iloc[0]["transcript"]  # 0:transcript
    except Exception:
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
                            transcript[0]: None
                            for transcript in result
                        }
                        TPM_Ribo = {
                            transcript[0]: None
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
                            transcript[0]: None
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
                    if not cds_start:
                        cdslen = None
                        threeutrlen = None
                    else:
                        cdslen = cds_stop - cds_start
                        threeutrlen = tranlen - cds_stop
                    if user == "test":
                        if transcript[0] in orfQuant_res:
                            OPM_coverage = orfQuant_res[transcript[0]]
                        else:
                            OPM_coverage = None

                        if transcript[0] in TPM_RNA:
                            RNA_coverage = TPM_RNA[transcript[0]]
                        else:
                            RNA_coverage = None

                        if transcript[0] in TPM_Ribo:
                            ribo_coverage = TPM_Ribo[transcript[0]]
                        else:
                            ribo_coverage = None

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

    user_short_passed = False
    if not data["user_short"] or user_short_passed == True:
        short_code = generate_short_code(data, organism, data["transcriptome"],
                                         "interactive_plot")
    else:
        short_code = data["user_short"]
        user_short_passed = True

    #Put any publicly available seq types (apart from riboseq and rnaseq) here
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

    new_conf = config.copy()

    #get user_id
    if current_user.is_authenticated:
        user_name = current_user.name
        user_id = get_user_id(current_user.name)
        user_conf = get_user_settings(user_id)
        for key, value in user_conf.items():
            new_conf[key] = value

        #get rules for all custom seq types
        cursor.execute(
            "SELECT * from seq_rules WHERE user_id = {};".format(user_id))
        result = (cursor.fetchall())
        for row in result:
            seq_name = row[1]
            frame_breakdown = row[2]
            seq_rules[seq_name] = {"frame_breakdown": frame_breakdown}
        connection.close()
    if tran:
        return riboflask.generate_plot(
            tran, minread, maxread, lite, ribocoverage, organism, readscore,
            noisered, primetype, minfiles, nucseq, user_hili_starts,
            user_hili_stops, uga_diff, file_paths_dict, short_code,
            color_readlen_dist, background_col, uga_col, uag_col, uaa_col,
            advanced, seqhili, seq_rules, title_size, subheading_size,
            axis_label_size, marker_size, transcriptome, config.UPLOADS_DIR,
            cds_marker_size, cds_marker_colour, legend_size, ribo_linewidth,
            secondary_readscore, pcr, mismatches, hili_start, hili_stop)

    else:
        return "ERROR! Could not find any transcript or gene corresponding to {}".format(
            tran)
