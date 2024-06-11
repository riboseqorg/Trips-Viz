from typing import Tuple
from flask import Blueprint, render_template, request, jsonify
from flask import current_app as app
from sqlitedict import SqliteDict
import os
import polars as pl
import config
from core_functions import (fetch_studies, fetch_files, fetch_study_info,
                            fetch_file_paths, generate_short_code, form_filler)
# import riboflask_compare
from flask_login import current_user
from fixed_values import my_decoder

from sqlqueries_2 import get_table, sqlquery, get_user_id
# Single transcript comparison page, user chooses a gene and groups of files to display
comparison_plotpage_blueprint = Blueprint("comparisonpage",
                                          __name__,
                                          template_folder="templates")


@comparison_plotpage_blueprint.route('/<organism>/<transcriptome>/comparison/')
def comparisonpage(organism: str, transcriptome: str) -> str:
    """
    Parameters:
    - organism (str): name of the organism
    - transcriptome (str): name of the transcript

    Returns:
    - html page
    """
    # global user_short_passed
    data = form_filler(organism, transcriptome)

    data['studyinfo_dict'] = fetch_study_info(
        data["gwips_info"][0, "organism_id"])

    data['box_colors'] = config.BOX_COLORS

    accepted_studies = fetch_studies(data["gwips_info"][0, "organism_id"])
    data['files'] = fetch_files(accepted_studies).to_pandas()
    print(data['files'])
    return render_template('index_compare.html', template_dict=data)


# Creates/serves the comparison plots
comparisonquery_blueprint = Blueprint("comparequery",
                                      __name__,
                                      template_folder="templates")


@comparisonquery_blueprint.route('/comparequery', methods=['POST'])
def comparequery() -> str | Tuple:
    """
    Parameters:
    - None

    Returns:
    - html page
    """
    # global user_short_passed
    data = request.data.to_dict()
    if not data["master_file_dict"]:
        return "Error: No files in the File list box. To add files to the file list box click on a study in the studies section above. This will populate the Ribo-seq and RNA-Seq sections with a list of files. Click on one of the files and then press the  Add button in the studies section. This will add the file to the File list box. Selecting another file and clicking Add again will add the new file to the same group in the File list. Alternatively to add a new group simply change the selected colour (by clicking on the coloured box in the studies section) and then click the Add file button."
    user_short_passed = False
    owner = get_table("organisms").filter(
        (pl.col('organism_name') == data['organism'])
        & (pl.col('transcriptome_list') == data["transcriptome"]))[0, "owner"]

    if owner:
        transhelve = "{0}/{1}/{2}/{2}.{3}.sqlite".format(
            config.SCRIPT_LOC, config.ANNOTATION_DIR, data['organism'],
            data['transcriptome'])
        if not os.path.isfile(transhelve):
            return f"Cannot find annotation file {data['organism']}.{data['transcriptome']}.sqlite"
    else:
        transhelve = "{0}/transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
            config.UPLOADS_DIR, owner, data['organism'], data['transcriptome'])
    transcripts = sqlquery(transhelve, "transcripts").filter(
        (pl.col('transcript') == data['transcriptome'])
        | (pl.col('gene') == data['transcriptome'])).unique(
            subset=['transcript'])
    if transcripts.is_empty():
        return f"ERROR! Could not find any transcript corresponding to {data['transcript']}"
    if transcripts.shape[0] > 1:
        # TODO: Add input element to select transcript
        transcripts = transcripts.with_columns(
            cdslen=pl.col("cds_stop") - pl.col("cds_start"),
            threeutrlen=pl.col("length") - pl.col("cds_stop"),
        ).select([
            "transcript", "length", "cds_start", "cds_stop", "principle",
            "cdslen", "threeutrlen"
        ]).to_pandas().to_html()

    for color in data["master_file_dict"]:
        files_ids = data["master_file_dict"][color]["file_ids"]
        files_infos = get_table("files").filter(
            pl.col("file_id").is_in(files_ids)).unique(subset=["file_id"])
        file_paths = fetch_file_paths(data)
        # TODO: Continue here

    master_filepath_dict = {}

    # This section is purely to sort by label alphabetically
    for color in master_file_dict:
        master_filepath_dict[color] = {
            "filepaths": [],
            "file_ids": [],
            "file_names": [],
            "file_descs": [],
            "mapped_reads": 0,
            "minread": minread,
            "maxread": maxread
        }

        for file_id in master_file_dict[color]["file_ids"]:
            trips_cursor.execute(
                "SELECT file_name,file_description,file_type from files WHERE file_id = {};"
                .format(file_id))
            result = (trips_cursor.fetchone())
            file_name = master_file_dict[color]["label"]
            file_paths = fetch_file_paths([file_id], data['organism'])

            for filetype in file_paths:
                for file_id in file_paths[filetype]:
                    filepath = file_paths[filetype][file_id]
                    if os.path.isfile(filepath):
                        sqlite_db = SqliteDict(f"{filepath}",
                                               autocommit=False,
                                               decode=my_decoder)
                    else:
                        return_str = "File not found, please report this to tripsvizsite@gmail.com or via the contact page."
                        if app.debug:
                            return return_str, "NO_CELERY", {'Location': None}
                        else:
                            return jsonify({
                                'current': 100,
                                'total': 100,
                                'status': 'return_str',
                                'result': return_str
                            }), 200, {
                                'Location': ""
                            }

                    if "noncoding_counts" in sqlite_db and "coding_counts" in sqlite_db:
                        master_filepath_dict[color]["mapped_reads"] += float(
                            sqlite_db["noncoding_counts"])
                        master_filepath_dict[color]["mapped_reads"] += float(
                            sqlite_db["coding_counts"])
                    else:
                        if "normalize" in data:
                            return_str = "One or more selected files is missing values for 'coding_counts' and 'non_coding_counts' so cannot normalize with these files, please report this to tripsvizsite@gmail.com or via the contact page."
                            if app.debug:
                                return return_str, "NO_CELERY", {
                                    'Location': None
                                }
                            else:
                                return jsonify({
                                    'current': 100,
                                    'total': 100,
                                    'status': 'return_str',
                                    'result': return_str
                                }), 200, {
                                    'Location': ""
                                }
                    master_filepath_dict[color]["filepaths"].append(filepath)
                    master_filepath_dict[color]["file_ids"].append(file_id)
                    master_filepath_dict[color]["file_names"].append(file_name)
                    master_filepath_dict[color]["file_descs"].append(result[1])
                    master_filepath_dict[color]["file_type"] = result[2]

    html_args = data["html_args"]
    if html_args["user_short"] == "None" or user_short_passed:
        data["short_code"] = generate_short_code(data, data['organism'],
                                                 html_args["transcriptome"],
                                                 "comparison")
    else:
        data["short_code"] = html_args["user_short"]
        user_short_passed = True

    if current_user.is_authenticated:
        user_id = get_user_id(current_user.name)
        data["user_settings"] = get_table("user_settings").filter(
            pl.col("user_id") == user_id)[0]
    else:
        data["user_settings"] = config.DEFAULT_USER_SETTINGS.copy()

    data['master_filepath_dict'] = master_filepath_dict
    return
    # if data['transcript']:
    # return riboflask_compare.generate_compare_plot(data)

    # return "ERROR! Could not find any transcript corresponding to {tran}"
