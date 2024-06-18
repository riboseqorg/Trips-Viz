from typing import Tuple
from plots import VegaPlot
import altair as alt
from flask import Blueprint, render_template, request, jsonify
from flask import current_app as app
from sqlitedict import SqliteDict
import os
from json import loads
from fetch_shelve_reads2 import get_reads
import polars as pl
import config
from core_functions import (fetch_studies, fetch_files, fetch_study_info,
                            fetch_file_paths, generate_short_code, form_filler,
                            string2other)
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
    data['list'] = 1
    return render_template('index_compare.html', template_dict=data)


def anmol(filepath_list, normalize):
    total_reads = 0
    for filepath in filepath_list:
        try:
            sqlite_db = SqliteDict(f"{filepath}",
                                   autocommit=False,
                                   decode=my_decoder)
        except FileNotFoundError:
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
        total_reads += sqlite_db["noncoding_counts"] + sqlite_db[
            "coding_counts"]
    else:
        if normalize:
            return_str = "One or more selected files is missing values for 'coding_counts' and 'non_coding_counts' so cannot normalize with these files, please report this to tripsvizsite@gmail.com or via the contact page."
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

    pass


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
    data = loads(list(request.form.to_dict().keys())[0])
    print(data)
    data = string2other(data)
    data['primetype'] = "fiveprime"  # NOTE: This is default for here
    data["readscore"] = 1  # NOTE: This is default for here
    print(data)

    if not data["groups"]:
        return "Error: No files in the File list box. To add files to the file list box click on a study in the studies section above. This will populate the Ribo-seq and RNA-Seq sections with a list of files. Click on one of the files and then press the  Add button in the studies section. This will add the file to the File list box. Selecting another file and clicking Add again will add the new file to the same group in the File list. Alternatively to add a new group simply change the selected colour (by clicking on the coloured box in the studies section) and then click the Add file button."
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
    transcripts = sqlquery(
        transhelve,
        "transcripts").filter((pl.col('transcript') == data['transcript'])
                              | (pl.col('gene') == data['transcript'])).unique(
                                  subset=['transcript'])
    print(transcripts, 'yyyy')
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

    files = fetch_file_paths(data)
    plots_list = []
    labels = []
    colors = []

    for group in data["groups"]:
        file_paths = files.filter(pl.col("file_id").is_in(data[group]))
        data['file_paths_dict'] = file_paths
        group_core = group.split("_")[0]
        group_label = group_core + '_label'
        group_color = group_core + '_color'
        labels.append(data[group_label])
        colors.append(data[group_color])
        group_name = data[group_label]
        reads_count = get_reads(data)[0].with_columns(frame=pl.lit(group_name))
        plots_list.append(reads_count)
    plots_list = pl.concat(plots_list)

    # plt = VegaPlot(reads_count)
    # plots_list.append(plt.line("pos", "count"))

    if current_user.is_authenticated:
        user_id = get_user_id(current_user.name)
        data["user_settings"] = get_table("user_settings").filter(
            pl.col("user_id") == user_id)[0]
    else:
        data["user_settings"] = config.DEFAULT_USER_SETTINGS.copy()

    # plot_json = plots_list[0]

    colors = alt.Scale(domain=labels, range=colors)
    # for plot in plots_list[1:]:
    # plot_json = plot_json + plot
    return VegaPlot(plots_list, colors).line("pos", "count").to_json()
