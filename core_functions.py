import string
from math import floor
import sqlite3
from flask import session, request
from flask_login import UserMixin, current_user
from Bio.Seq import Seq
from sqlqueries import sqlquery, table2dict
import uuid
import config
import logging


# User model
class User(UserMixin):

    def __init__(self, id):
        self.id = id
        self.name = str(id)
        self.password = self.name + "_secret"

    def is_authenticated(self):
        return self.is_authenticated

    def get_id(self):
        return str(self.id)

    def __repr__(self):
        return "%d/%s/%s" % (self.id, self.name, self.password)


def fetch_user():
    consent = request.cookies.get("cookieconsent_status")
    #If user rejects cookies then do not track them and delete all other cookies
    if consent == "deny":
        return (None, False)
    logging.debug("fetch_user Connecting to trips.sqlite")
    connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
    connection.text_factory = str
    cursor = connection.cursor()
    session.permanent = True
    if "uid" not in session:
        session["uid"] = uuid.uuid4()
    session_id = session["uid"]
    #Check if this session uid is already in the users table
    cursor.execute(
        "Select username from users where username = '{}';".format(session_id))
    result = cursor.fetchone()
    if result == None:
        #Add session uid to user table
        cursor.execute(
            "INSERT INTO users VALUES (NULL,'{}',NULL,'-1','',0,1);".format(
                session_id))
        cursor.execute(
            "SELECT user_id FROM users WHERE username = '{}'".format(
                session_id))
        user_id = cursor.fetchone()[0]
        #marker_size, axis_label_size, subheading_size, title_size, user_id, background_col, readlength_col, metagene_fiveprime_col, metagene_threeprime_col, nuc_comp_a_col ,
        #nuc_comp_t_col , nuc_comp_g_col , nuc_comp_c_col , uga_col , uag_col , uaa_col , comp_uga_col , comp_uag_col , comp_uaa_col , cds_marker_width ,
        #cds_marker_colour VARCHAR, legend_size , ribo_linewidth
        cursor.execute(
            "INSERT INTO user_settings VALUES ('{0}','{1}','{2}','{3}',{4},'{5}','{6}','{7}','{8}','{9}','{10}','{11}','{12}','{13}','{14}','{15}','{16}','{17}','{18}',{19},'{20}',{21},{22});"
            .format(config.MARKER_SIZE, config.AXIS_LABEL_SIZE,
                    config.SUBHEADING_SIZE, config.TITLE_SIZE, user_id,
                    config.BACKGROUND_COL, config.READLENGTH_COL,
                    config.METAGENE_FIVEPRIME_COL,
                    config.METAGENE_THREEPRIME_COL, config.A_COL, config.T_COL,
                    config.G_COL, config.C_COL, config.UGA_COL, config.UAG_COL,
                    config.UAA_COL, config.UGA_COL, config.UAG_COL,
                    config.UAA_COL, config.CDS_MARKER_SIZE,
                    config.CDS_MARKER_COLOUR, config.LEGEND_SIZE,
                    config.RIBO_LINEWIDTH))
        connection.commit()

    #If user is logged in current_user.name will return a username, otherwise set the user to the session uid
    try:
        user = current_user.name
        logged_in = True
    except:
        user = session_id
        logged_in = False
    logging.debug("fetch_user Closing trips.sqlite connection")
    connection.close()
    return (user, logged_in)


# Given a username and an organism returns a list of relevant studies.
def fetch_studies(organism, transcriptome):
    dbpath = '{}/{}'.format(config.SCRIPT_LOC, config.DATABASE_NAME)
    study_access_list = []
    #get a list of organism id's this user can access
    if current_user.is_authenticated:
        # getting user ID
        result = sqlquery(dbpath, "users")  # users is name of table
        user_id = result[result.username ==
                         current_user.name].user_id.values[0]
        # Study accession list
        result = sqlquery(dbpath, "study_accesion")  # users is name of table
        study_access_list = result.loc[result.user_id == user_id,
                                       "study_id"].values

    # Getting organism id
    result = sqlquery(dbpath, "organisms")  # users is name of table
    organism_id = result.loc[(
        result.organism_name == organism and result.transcriptome ==
        transcriptome), "organism_id"].values[
            0]  # NOTE: Expecting it has to be true as it is encoded in form
    organism_id = int(organism_id)  # TODO: Check if it is int from start

    # Getting studies
    studies = sqlquery(dbpath, "studies")  # users is name of table
    studies = studies.loc[studies.organism_id == organism_id,
                          ["study_id", "study_name", "private"]]
    studies = studies.loc[
        studies.private == 0 or studies.study_id.isin(study_access_list),
        ["study_id", "study_name"]]
    return studies  # Accepted studies


# Create a dictionary of files seperated by type, this allows for file type grouping on the front end.
def fetch_files(accepted_studies):
    dbpath = '{}/{}'.format(config.SCRIPT_LOC, config.DATABASE_NAME)
    files = sqlquery(dbpath, "files")
    files = files.loc[
        files.study_id.isin(accepted_studies['study_id']),
        ["file_id", "study_id", "file_name", "file_description", "file_type"
         ]].sort_values("file_description")  # Files Details
    files = files.merge(accepted_studies, on='study_id')
    files['file_name'] = files['file_name'].apply(
        lambda x: x.replace('.self', ''))
    files['study'] = files.apply(lambda x: (x['study_id'], x['study_name']),
                                 axis=1)  # Paired info (ID, Name)
    files['file'] = files.apply(lambda x: (x['file_id'], x['file_name']),
                                axis=1)  # Paired info (ID, Name)
    files = files.drop(['study_id', 'study_name', 'file_id', 'file_name'],
                       axis=1)

    files = table2dict(files, ['file_type', 'study', 'file'])
    return files


# Gets a list of all studies associated with an organism
def fetch_study_info(organism):
    dbpath = '{}/{}'.format(config.SCRIPT_LOC, config.DATABASE_NAME)
    studies = sqlquery(dbpath, "studies")[[
        "study_id", "paper_authors", "srp_nos,paper_year", "paper_pmid",
        "paper_link", "gse_nos", "adapters", "paper_title", "description",
        "study_name"
        #, "organism_id" # TODO: Consider this for selecting
    ]]
    # "paper_link": row[5].strip('"'),  # generate link using pubmed id
    return table2dict(studies, ["study_id"])


# Given a list of file id's as strings returns a list of filepaths to the sqlite files.
def fetch_file_paths(file_list, organism):
    file_path_dict = {"riboseq": {}, "rnaseq": {}, "proteomics": {}}
    #Convert to a tuple so it works with mysql
    try:
        #Remove empty strings from file list
        file_list[:] = [x for x in file_list if x]
        #Convert each string to an int
        int_file_list = [int(x) for x in file_list]
    except:
        return {}
    dbpath = '{}/{}'.format(config.SCRIPT_LOC, config.DATABASE_NAME)
    studies = sqlquery(dbpath, "studies")[["study_id", "study_name"]]
    files = sqlquery(dbpath, "files")
    files = files[files["file_id"].isin(int_file_list)]
    # ----
    connection_ffp = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                    config.DATABASE_NAME))
    connection_ffp.text_factory = str
    cursor = connection_ffp.cursor()

    study_dict = {}
    cursor.execute("SELECT study_id,study_name from studies;")
    result = cursor.fetchall()
    for row in result:
        study_dict[row[0]] = row[1]
    if int_file_list:
        cursor.execute(
            "SELECT file_type,file_name,study_id,file_id,owner from files WHERE file_id in {};"
            .format(str(int_file_list).replace("[", "(").replace("]", ")")))
        result = cursor.fetchall()
        for row in result:
            study_name = study_dict[row[2]]
            if row[0] not in file_path_dict:
                file_path_dict[row[0]] = {}
            # All files uploaded on our side are owned by user id 1, all others are user uploaded
            # Replace .shelf with .sqlite to deal with legacy files
            if row[4] == 1:
                file_path_dict[row[0]][row[3]] = (
                    "{}/{}/{}/{}/{}/{}.sqlite".format(
                        config.SCRIPT_LOC, config.SQLITES_DIR, row[0],
                        organism, study_name,
                        row[1].replace(".shelf", "").replace(".sqlite", "")))
            else:
                file_path_dict[row[0]][row[3]] = ("{}/{}/{}.sqlite".format(
                    config.UPLOADS_DIR, study_name,
                    row[1].replace(".shelf", "").replace(".sqlite", "")))
    logging.debug("fetch_file_paths closing connection")
    connection_ffp.close()
    return file_path_dict


# Builds a url and inserts it into sqlite database
def generate_short_code(data, organism, transcriptome, plot_type):
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    connection.text_factory = str
    cursor = connection.cursor()
    # build a url so that this plot can be recreated later on
    url = "/{}/{}/{}/?".format(organism, transcriptome, plot_type)
    # To stop the file_list argument in the url being extremely long, pass a riboseq_studies and rnaseq_studies arguments, these are study_ids of any study which
    # has all of its riboseq or rnaseq files checked. This part of the code determines which studies fall into that category.
    riboseq_studies = []
    rnaseq_studies = []
    proteomics_studies = []

    cursor.execute(
        "SELECT organism_id from organisms WHERE organism_name = '{}'".format(
            organism))
    result = cursor.fetchone()
    organism_id = int(result[0])
    cursor.execute(
        "SELECT study_id from studies WHERE organism_id = {}".format(
            organism_id))
    result = cursor.fetchall()

    if plot_type == "interactive_plot" or plot_type == "metainfo_plot" or plot_type == "orf_translation":
        url += "files="
        for row in result:
            study_id = int(row[0])

            #Now get all riboseq files that have this study_id, if all those ids are in file_list, add to riboseq studies and remove those files from file_list, do the same for rnaseq and proteomics
            cursor.execute(
                "SELECT file_id from files WHERE study_id = {} AND file_type = 'riboseq'"
                .format(study_id))
            result = cursor.fetchall()
            all_present = True
            for row in result:
                if str(row[0]) not in data["file_list"]:
                    all_present = False
            if all_present == True:
                riboseq_studies.append(study_id)
                for row in result:
                    data["file_list"].remove(str(row[0]))

            cursor.execute(
                "SELECT file_id from files WHERE study_id = {} AND file_type = 'rnaseq'"
                .format(study_id))
            result = cursor.fetchall()
            all_present = True
            # If there are no files of that type for that study then don't bother adding it to the list
            if not result:
                all_present = False
            for row in result:
                if str(row[0]) not in data["file_list"]:
                    all_present = False
            if all_present == True:
                rnaseq_studies.append(study_id)
                for row in result:
                    data["file_list"].remove(str(row[0]))

            cursor.execute(
                "SELECT file_id from files WHERE study_id = {} AND file_type = 'proteomics'"
                .format(study_id))
            result = cursor.fetchall()
            all_present = True
            # If there are no files of that type for that study then don't bother adding it to the list
            if not result:
                all_present = False
            for row in result:
                if str(row[0]) not in data["file_list"]:
                    all_present = False
            if all_present == True:
                proteomics_studies.append(study_id)
                for row in result:
                    data["file_list"].remove(str(row[0]))

        for filenum in data["file_list"]:
            url += filenum + ","
        if riboseq_studies:
            url += "&ribo_studies="
            for study_id in riboseq_studies:
                url += str(study_id) + ","
        if rnaseq_studies:
            url += "&rna_studies="
            for study_id in rnaseq_studies:
                url += str(study_id) + ","
        if proteomics_studies:
            url += "&proteomics_studies="
            for study_id in proteomics_studies:
                url += str(study_id) + ","
    if plot_type == "interactive_plot":
        url += "&tran={}".format(data['transcript'].upper().strip())
        url += "&minread={}".format(data['minread'])
        url += "&maxread={}".format(data['maxread'])
        url += "&user_dir={}".format(data["primetype"])
        if "ambiguous" in data:
            url += "&ambig=T"
        else:
            url += "&ambig=F"
        if "ribocoverage" in data:
            url += "&cov=T"
        else:
            url += "&cov=F"
        if "varlite" in data:
            url += "&lg=T"
        else:
            url += "&lg=F"
        if "nucseq" in data:
            url += "&nuc=T"
        else:
            url += "&nuc=F"
        if "readscore" in data:
            url += "&rs={}".format(data["readscore"])
        else:
            url += "&rs=1"
        if "color_readlen_dist" in data:
            url += "&crd=T"
        else:
            url += "&crd=F"

    if plot_type == "traninfo_plot":
        url += "&plot={}".format(data['plottype'])
        #nuc_comp_single, nuc_comp_multi, lengths_plot, gene_count, codon_usage, nuc_freq_plot, orfstats
        if data["plottype"] == "nuc_comp_single":
            url += "&metagene_tranlist={}".format(data["metagene_tranlist"])
        if data["plottype"] == "nuc_comp_multi":
            url += "&gc_tranlist={}".format(data["gc_tranlist"])
            url += "&gc_tranlist2={}".format(data["gc_tranlist2"])
            url += "&gc_tranlist3={}".format(data["gc_tranlist3"])
            url += "&gc_tranlist4={}".format(data["gc_tranlist4"])
            url += "&gc_location={}".format(data["gc_location"])
            url += "&nucleotide={}".format(data["nucleotide"])
            url += "&plot_type={}".format(data["plot_type"])
        if data["plottype"] == "lengths_plot":
            url += "&gc_tranlist={}".format(data["gc_tranlist"])
            url += "&gc_tranlist2={}".format(data["gc_tranlist2"])
            url += "&gc_tranlist3={}".format(data["gc_tranlist3"])
            url += "&gc_tranlist4={}".format(data["gc_tranlist4"])
            url += "&gc_location={}".format(data["gc_location"])
            url += "&plot_type={}".format(data["plot_type"])
        if data["plottype"] == "codon_usage":
            url += "&gc_tranlist={}".format(data["gc_tranlist"])
        if data["plottype"] == "nuc_freq_plot":
            url += "&nuc_freq_plot_anchor={}".format(
                data["nuc_freq_plot_anchor"])
            url += "&nuc_freq_plot_window={}".format(
                data["nuc_freq_plot_window"])
            url += "&nuc_freq_plot_tranlist={}".format(
                data["nuc_freq_plot_tranlist"])

    if plot_type == "metainfo_plot":
        url += "&plot={}".format(data['plottype'])

        # Nuc comp plot
        if data["plottype"] == "nuc_comp":
            if "nuc_comp_direction" in data:
                if data["nuc_comp_direction"] != "None":
                    url += "&nc_dir={}".format(data["nuc_comp_direction"])
            if "nuc_comp_type" in data:
                if data["nuc_comp_type"] != "None":
                    url += "&nc_type={}".format(data["nuc_comp_type"])
            if "nuc_minreadlen" in data:
                if data["nuc_minreadlen"] != "None":
                    url += "&nc_minreadlen={}".format(data["nuc_minreadlen"])
            if "nuc_maxreadlen" in data:
                if data["nuc_maxreadlen"] != "None":
                    url += "&nc_maxreadlen={}".format(data["nuc_maxreadlen"])

        # Heatmap plot
        if data["plottype"] == "heatmap":
            if "heatmap_minreadlen" in data:
                if data["heatmap_minreadlen"] != "None":
                    url += "&hm_minreadlen={}".format(
                        data["heatmap_minreadlen"])
            if "heatmap_maxreadlen" in data:
                if data["heatmap_maxreadlen"] != "None":
                    url += "&hm_maxreadlen={}".format(
                        data["heatmap_maxreadlen"])
            if "heatmap_direction" in data:
                if data["heatmap_direction"] != "None":
                    url += "&hm_dir={}".format(data["heatmap_direction"])
            if "log_scale" in data:
                url += "&hm_log=T"
            else:
                url += "&hm_log=F"
            if "reverse_scale" in data:
                url += "&hm_rev=T"
            else:
                url += "&hm_rev=F"
            if "heatmap_metagene_type" in data:
                if data["heatmap_metagene_type"] != "None":
                    url += "&hm_pos={}".format(data["heatmap_metagene_type"])
            if "heatmap_startpos" in data:
                if data["heatmap_startpos"] != "None":
                    url += "&hm_start={}".format(data["heatmap_startpos"])
            if "heatmap_endpos" in data:
                if data["heatmap_endpos"] != "None":
                    url += "&hm_stop={}".format(data["heatmap_endpos"])
            if "color_palette" in data:
                if data["color_palette"] != "None":
                    url += "&hm_col={}".format(data["color_palette"])
            if "maxscaleval" in data:
                url += "&maxscaleval={}".format(data["maxscaleval"])
            if "metagene_tranlist" in data:
                if data["metagene_tranlist"] != "None":
                    url += "&metagene_tranlist={}".format(
                        data["metagene_tranlist"])

        # Triplet periodicity
        if data["plottype"] == "trip_periodicity":
            if "trip_minreadlen" in data:
                if data["trip_minreadlen"] != "None":
                    url += "&tp_minreadlen={}".format(data["trip_minreadlen"])
            if "trip_maxreadlen" in data:
                if data["trip_maxreadlen"] != "None":
                    url += "&tp_maxreadlen={}".format(data["trip_maxreadlen"])

        # mRNA readlen dist
        if data["plottype"] == "mrna_dist_readlen":
            if "mrna_readlen_per" in data:
                if data["mrna_readlen_per"] != "None":
                    url += "&mdr_per={}".format(data["mrna_readlen_per"])
            if "smooth_amount" in data:
                if data["smooth_amount"] != "None":
                    url += "&mdr_smooth={}".format(data["smooth_amount"])

        # metagene
        if data["plottype"] == "metagene_plot":
            if "include_first" in data:
                if data["include_first"] != "None":
                    url += "&include_first=T"
            if "include_last" in data:
                if data["include_last"] != "None":
                    url += "&include_last=T"
            if "exclude_first" in data:
                if data["exclude_first"] != "None":
                    url += "&exclude_first=T"
            if "exclude_last" in data:
                if data["exclude_last"] != "None":
                    url += "&exclude_last=T"
            if "custom_seq_list" in data:
                if data["custom_seq_list"] != "None":
                    url += "&custom_seq_list={}".format(
                        data["custom_seq_list"])
            if "exclude_first_val" in data:
                if data["exclude_first_val"] != "None":
                    url += "&exclude_first_val={}".format(
                        data["exclude_first_val"])
            if "exclude_last_val" in data:
                if data["exclude_last_val"] != "None":
                    url += "&exclude_last_val={}".format(
                        data["exclude_last_val"])
            if "include_first_val" in data:
                if data["include_first_val"] != "None":
                    url += "&include_first_val={}".format(
                        data["include_first_val"])
            if "include_last_val" in data:
                if data["include_last_val"] != "None":
                    url += "&include_last_val={}".format(
                        data["include_last_val"])
            if "metagene_tranlist" in data:
                if data["metagene_tranlist"] != "None":
                    url += "&metagene_tranlist={}".format(
                        data["metagene_tranlist"])
            if "metagene_type" in data:
                if data["metagene_type"] != "None":
                    url += "&mg_pos={}".format(data["metagene_type"])
            if "minreadlen" in data:
                if data["minreadlen"] != "None":
                    url += "&mg_minreadlen={}".format(data["minreadlen"])
            if "maxreadlen" in data:
                if data["maxreadlen"] != "None":
                    url += "&mg_maxreadlen={}".format(data["maxreadlen"])

        # Replicate comparison
        if data["plottype"] == "replicate_comp":
            if "minimum_reads" in data:
                if data["minimum_reads"] != "None":
                    url += "&rp_minreads={}".format(data["minimum_reads"])

    if plot_type == "comparison":
        url += "files="
        file_string = ""

        for color in data["master_file_dict"]:
            for file_id in data["master_file_dict"][color]["file_ids"]:
                file_string += "{},".format(file_id)
            # Can't use # in html args so encode as %23 instead
            file_string += "{}_".format(color.replace("#", "%23"))

        #remove the trailling _ from file_string
        file_string = file_string[:len(file_string) - 1]

        label_string = ""
        for color in data["master_file_dict"]:
            label_string += "{},".format(
                data["master_file_dict"][color]["label"])
            # Can't use # in html args so encode as %23 instead
            label_string += "{}_".format(color.replace("#", "%23"))

        #remove the trailling _ from file_string
        label_string = label_string[:len(label_string) - 1]

        url += file_string
        url += "&labels={}".format(label_string)
        url += "&transcript={}".format(data['transcript'].upper().strip())
        url += "&minread={}".format(data['minread'])
        url += "&maxread={}".format(data['maxread'])

        url += "&hili_start={}".format(data['hili_start'])
        url += "&hili_stop={}".format(data['hili_stop'])

        if "ambiguous" in data:
            url += "&ambig=T"
        else:
            url += "&ambig=F"

        if "coverage" in data:
            url += "&cov=T"
        else:
            url += "&cov=F"

        if "normalize" in data:
            url += "&normalize=T"
        else:
            url += "&normalize=F"

    if plot_type == "differential":
        url += "&minzscore={}".format(data["minzscore"])
        url += "&minread={}".format(data["minreads"])
        url += "&region={}".format(data["region"])
        url += "&riboseq_files_1={}".format(
            str(data["master_file_dict"]["riboseq1"]["file_ids"]).strip(
                "[]").replace("'", "").replace(" ", ""))
        url += "&riboseq_files_2={}".format(
            str(data["master_file_dict"]["riboseq2"]["file_ids"]).strip(
                "[]").replace("'", "").replace(" ", ""))
        url += "&rnaseq_files_1={}".format(
            str(data["master_file_dict"]["rnaseq1"]["file_ids"]).strip(
                "[]").replace("'", "").replace(" ", ""))
        url += "&rnaseq_files_2={}".format(
            str(data["master_file_dict"]["rnaseq2"]["file_ids"]).strip(
                "[]").replace("'", "").replace(" ", ""))
        url += "&riboseq_labels_1={}".format(
            str(data["master_file_dict"]["riboseq1"]["file_names"]).strip(
                "[]").replace("'", "").replace(" ", ""))
        url += "&riboseq_labels_2={}".format(
            str(data["master_file_dict"]["riboseq2"]["file_names"]).strip(
                "[]").replace("'", "").replace(" ", ""))
        url += "&rnaseq_labels_1={}".format(
            str(data["master_file_dict"]["rnaseq1"]["file_names"]).strip(
                "[]").replace("'", "").replace(" ", ""))
        url += "&rnaseq_labels_2={}".format(
            str(data["master_file_dict"]["rnaseq2"]["file_names"]).strip(
                "[]").replace("'", "").replace(" ", ""))
        if "plottype" in data:
            url += "&plottype={}".format(data["plottype"])
        if "min_cov" in data:
            url += "&min_cov={}".format(data["min_cov"])
        if "gene_list" in data:
            url += "&gene_list={}".format(data["gene_list"])
        if "ambiguous" in data:
            url += "&ambig=T"
        else:
            url += "&ambig=F"

    if plot_type == "orf_translation":
        start_codons = []
        if "sc_aug" in data:
            start_codons.append("AUG")
        if "sc_cug" in data:
            start_codons.append("CUG")
        if "sc_gug" in data:
            start_codons.append("GUG")
        if "sc_none" in data:
            start_codons.append("NONE")

        url += "&start_codons={}".format(
            str(start_codons).strip("[]").replace("'", ""))
        url += "&min_cds={}&max_cds={}&min_len={}&max_len={}&min_avg={}&max_avg={}".format(
            data["min_cds"], data["max_cds"], data["min_len"], data["max_len"],
            data["min_avg"], data["max_avg"])
        url += "&tran_list={}".format(data["tran_list"])

        if "start_increase_check" in data:
            url += "&sic=T"
        else:
            url += "&sic=F"

        if "stop_decrease_check" in data:
            url += "&sdc=T"
        else:
            url += "&sdc=F"

        if "lowest_frame_diff_check" in data:
            url += "&lfdc=T"
        else:
            url += "&lfdc=F"

        if "highest_frame_diff_check" in data:
            url += "&hfdc=T"
        else:
            url += "&hfdc=F"

        if "ambig_check" in data:
            url += "&ambig=T"
        else:
            url += "&ambig=F"

        if "saved_check" in data:
            url += "&saved_check=T"
        else:
            url += "&saved_check=F"

    cursor.execute("SELECT MAX(url_id) from urls;")
    result = cursor.fetchone()
    #If the url table is empty result will return none
    if result[0] == None:
        url_id = 0
    else:
        url_id = int(result[0]) + 1
    cursor.execute("INSERT INTO urls VALUES({},'{}')".format(url_id, url))
    connection.commit()
    short_code = integer_to_base62(url_id)
    connection.close()
    return short_code


# Converts an integer to base62, needed to encode short urls
def integer_to_base62(num):
    base = string.digits + string.ascii_lowercase + string.ascii_uppercase
    r = num % 62
    res = base[r]
    q = num // 62
    while q:
        r = q % 62
        q //= 62
        res = base[int(r)] + res
    return res


# Converts a base62 encoded string to an integer, needed to decode short urls
def base62_to_integer(base62_str):
    base = string.digits + string.ascii_lowercase + string.ascii_uppercase
    res = 0
    for i in base62_str:
        res = 62 * res + base.find(base62_str[i])
    return res


#Takes a nucleotide string and returns the amino acide sequence
def nuc_to_aa(nuc_seq):
    return str(Seq(nuc_seq).translate())


# Calculates the coverage of each gene, for 5' leader, cds and 3' trailer for unambiguous and ambigous reads, needed for diff exp
def calculate_coverages(sqlite_db, longest_tran_list, ambig_type, region,
                        traninfo_dict):
    coverage_dict = {
        "unambig_fiveprime_coverage": {},
        "unambig_cds_coverage": {},
        "unambig_threeprime_coverage": {},
        "unambig_all_coverage": {},
        "ambig_fiveprime_coverage": {},
        "ambig_cds_coverage": {},
        "ambig_threeprime_coverage": {},
        "ambig_all_coverage": {}
    }
    total_trans = 0
    for tran in longest_tran_list:
        total_trans += 1
        unambig_dict = {}
        ambig_dict = {}
        coverage_dict["unambig_fiveprime_coverage"][tran] = 0
        coverage_dict["unambig_cds_coverage"][tran] = 0
        coverage_dict["unambig_threeprime_coverage"][tran] = 0
        coverage_dict["unambig_all_coverage"][tran] = 0
        coverage_dict["ambig_fiveprime_coverage"][tran] = 0
        coverage_dict["ambig_cds_coverage"][tran] = 0
        coverage_dict["ambig_threeprime_coverage"][tran] = 0
        coverage_dict["ambig_all_coverage"][tran] = 0
        if tran in sqlite_db and tran in traninfo_dict:
            cds_start = traninfo_dict[tran]["cds_start"]
            cds_stop = traninfo_dict[tran]["cds_stop"]
            tranlen = float(traninfo_dict[tran]["length"])
            for readlen in sqlite_db[tran]["unambig"]:
                for pos in sqlite_db[tran]["unambig"][readlen]:
                    for i in range(pos, pos + readlen):
                        unambig_dict[i] = ""
                        ambig_dict[i] = ""
            for readlen in sqlite_db[tran]["ambig"]:
                for pos in sqlite_db[tran]["ambig"][readlen]:
                    for i in range(pos, pos + readlen):
                        ambig_dict[i] = ""
            coverage_dict["unambig_all_coverage"][tran] = len(
                unambig_dict.keys()) / tranlen
            coverage_dict["ambig_all_coverage"][tran] = len(
                ambig_dict.keys()) / tranlen
            if cds_start != "None":
                cds_start = float(cds_start)
                cds_stop = float(cds_stop)
                cds_len = cds_stop - cds_start
                three_len = tranlen - cds_stop
                #Use list comprehension to count number of entries less than cds_start
                if cds_start > 0:
                    coverage_dict["unambig_fiveprime_coverage"][tran] = sum(
                        i < cds_start for i in unambig_dict.keys()) / cds_start
                    coverage_dict["ambig_fiveprime_coverage"][tran] = sum(
                        i < cds_start for i in ambig_dict.keys()) / cds_start
                coverage_dict["unambig_cds_coverage"][tran] = sum(
                    i > cds_start and i < cds_stop
                    for i in unambig_dict.keys()) / cds_len
                coverage_dict["ambig_cds_coverage"][tran] = sum(
                    i > cds_start and i < cds_stop
                    for i in ambig_dict.keys()) / cds_len
                if three_len > 0:
                    coverage_dict["unambig_threeprime_coverage"][tran] = sum(
                        i > cds_stop for i in unambig_dict.keys()) / three_len
                    coverage_dict["ambig_threeprime_coverage"][tran] = sum(
                        i > cds_stop for i in ambig_dict.keys()) / three_len
        sqlite_db["unambiguous_fiveprime_coverage"] = coverage_dict[
            "unambig_fiveprime_coverage"]
        sqlite_db["unambiguous_cds_coverage"] = coverage_dict[
            "unambig_cds_coverage"]
        sqlite_db["unambiguous_threeprime_coverage"] = coverage_dict[
            "unambig_threeprime_coverage"]
        sqlite_db["unambiguous_all_coverage"] = coverage_dict[
            "unambig_all_coverage"]
        sqlite_db["ambiguous_fiveprime_coverage"] = coverage_dict[
            "ambig_fiveprime_coverage"]
        sqlite_db["ambiguous_cds_coverage"] = coverage_dict[
            "ambig_cds_coverage"]
        sqlite_db["ambiguous_threeprime_coverage"] = coverage_dict[
            "ambig_threeprime_coverage"]
        sqlite_db["ambiguous_all_coverage"] = coverage_dict[
            "ambig_all_coverage"]
        sqlite_db.commit()


# Builds a profile, applying offsets
def build_profile(trancounts, offsets, ambig, minscore=None, scores=None):
    #print ("trancounts", trancounts)
    #print ("minscore", minscore)
    minreadlen = 15
    maxreadlen = 150
    profile = {}
    try:
        unambig_trancounts = trancounts["unambig"]
    except:
        unambig_trancounts = {}
    try:
        ambig_trancounts = trancounts["ambig"]
    except:
        ambig_trancounts = {}
    for readlen in unambig_trancounts:
        #print ("readlen", readlen)
        if minscore != None:
            if readlen in scores:
                if scores[readlen] < minscore:
                    #print ("continuing1",scores[readlen])
                    continue
            else:
                #print ("continuing2")
                continue
        if readlen < minreadlen or readlen > maxreadlen:
            continue
        try:
            offset = offsets[readlen] + 1
        except:
            offset = 15
        for pos in unambig_trancounts[readlen]:
            count = unambig_trancounts[readlen][pos]
            offset_pos = pos + offset
            try:
                profile[offset_pos] += count
            except:
                profile[offset_pos] = count

    if ambig == True:
        for readlen in ambig_trancounts:
            if minscore != None:
                if readlen in scores:
                    if scores[readlen] < minscore:
                        continue
                else:
                    continue
            if readlen < minreadlen or readlen > maxreadlen:
                continue
            try:
                offset = offsets[readlen] + 1
            except:
                offset = 15
            for pos in ambig_trancounts[readlen]:
                count = ambig_trancounts[readlen][pos]
                offset_pos = pos + offset
                try:
                    profile[offset_pos] += 0
                except:
                    profile[offset_pos] = count
    #print ("RETURNING PROFILE", profile)
    return profile


# Builds a profile, applying offsets
def build_proteomics_profile(trancounts, ambig):
    minreadlen = 15
    maxreadlen = 150
    profile = {}
    try:
        unambig_trancounts = trancounts["unambig"]
    except:
        unambig_trancounts = {}
    try:
        ambig_trancounts = trancounts["ambig"]
    except:
        ambig_trancounts = {}
    for readlen in unambig_trancounts:
        if readlen < minreadlen or readlen > maxreadlen:
            continue

        for pos in unambig_trancounts[readlen]:
            # Rather than add the whole count to each position, we divide by the length of the peptide first,
            # That way when we add the reduced count at each position, in total it will add up to the original count
            # and prevent a bias toward longer peptides, this allows us to count a fraction of a peptide that overlaps with an ORF
            # rather than counting an arbitrary position like the 5' end or 3' end which may fall outside the ORF in question.
            count = unambig_trancounts[readlen][pos] / float(readlen / 3)
            for x in range(pos, pos + readlen, 3):
                try:
                    profile[x] += count
                except:
                    profile[x] = count
    return profile


# Creates nucleotide composition counts if they don't already exist.
def get_nuc_comp_reads(sqlite_db, nuccomp_reads, organism, transcriptome):
    transhelve = sqlite3.connect("{0}{1}/{1}.v2.sqlite".format(
        config.ANNOTATION_DIR, organism))
    cursor = transhelve.cursor()
    cursor.execute(
        "SELECT transcript,cds_start,cds_stop,sequence from transcripts WHERE principal = 1"
    )
    result = cursor.fetchall()
    master_dict = {}
    offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
    for row in result:
        tran = row[0]
        cds_start = int(row[1])
        cds_stop = int(row[2])
        seq = row[3].replace("T", "U")
        if tran in sqlite_db:
            counts = sqlite_db[tran]["unambig"]
            for readlen in counts:
                if readlen not in master_dict:
                    master_dict[readlen] = {}
                    for i in range(0, int(readlen)):
                        master_dict[readlen][i] = {
                            "A": 0,
                            "T": 0,
                            "G": 0,
                            "C": 0,
                            "N": 0
                        }
                for pos in counts[readlen]:
                    count = counts[readlen][pos]
                    if readlen in offsets:
                        offset_pos = pos + offsets[readlen]
                    else:
                        offset_pos = pos + 15
                    if offset_pos >= cds_start and offset_pos <= cds_stop:
                        readframe = offset_pos % 3
                        cds_frame = (cds_start + 2) % 3
                        if readframe == cds_frame:
                            inframe = True
                        else:
                            inframe = False
                        if nuccomp_reads == "inframe" and inframe == False:
                            continue
                        if nuccomp_reads == "offrame" and inframe == True:
                            continue
                        readseq = seq[pos:pos + readlen]
                        for i in range(0, len(readseq)):
                            char = readseq[i].replace("U", "T")
                            master_dict[readlen][i][char] += count
    #save results so they won't have to be computed again later
    if "nuc_counts" in sqlite_db:
        new_nuc_counts = sqlite_db["nuc_counts"]
        new_nuc_counts[nuccomp_reads] = master_dict
        sqlite_db["nuc_counts"] = new_nuc_counts
    else:
        sqlite_db["nuc_counts"] = {nuccomp_reads: master_dict}
    sqlite_db.commit()
    transhelve.close()
    return master_dict


# Creates a dictionary of readlength counts
def fetch_rld(sqlite_db, ambig_type):
    rld = {}
    for transcript in sqlite_db:
        try:
            transcript_dict = sqlite_db[transcript]["unambig"]
        except:
            continue
        try:
            for rl in transcript_dict:
                for pos in transcript_dict[rl]:
                    count = transcript_dict[rl][pos]
                    if rl not in rld:
                        rld[rl] = 0
                    rld[rl] += count
                if ambig_type == "ambig":
                    transcript_dict = sqlite_db[transcript]["ambig"]
                    for rl in transcript_dict:
                        for pos in transcript_dict[rl]:
                            count = transcript_dict[rl][pos]
                            if rl not in rld:
                                rld[rl] = 0
                            rld[rl] += count
        except Exception as e:
            continue
    if ambig_type == "unambig":
        sqlite_db["unambig_read_lengths"] = rld
    elif ambig_type == "ambig":
        sqlite_db["read_lengths"] = rld
    sqlite_db.commit()
    return rld


def fetch_filename_file_id(file_id):
    '''
	return the filename from the database given a file id
	'''
    dbpath = '{}/{}'.format(config.SCRIPT_LOC, config.DATABASE_NAME)
    files = sqlquery(dbpath, "files")
    return files.loc[files["file_id"] == file_id, 'file_name'].values[0]
