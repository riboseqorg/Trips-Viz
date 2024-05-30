XCAPTCHA_ENABLED = True
XCAPTCHA_SITE_KEY = ""
XCAPTCHA_SECRET_KEY = ""
PASSWORD_SALT = "password_salt"
FLASK_SECRET_KEY = "flask_secret_key"
SECRET_KEY = "YouAreAFool"
DEFAULT_USER_SETTINGS = {  # This should match the table columns
    "background_col": "#F2F2F2",
    "uag_col": "#000000",
    "uga_col": "#000000",
    "uaa_col": "#000000",
    "readlength_col": "#ff5f5b",
    "metagene_fiveprime_col": "#ff5f5b",
    "metagene_threeprime_col": "#9acaff",
    "nuc_comp_a_col": "#ff5f5b",
    "nuc_comp_t_col": "#90e090",
    "nuc_comp_g_col": "#9acaff",
    "nuc_comp_c_col": "#ffff91",
    "title_size": 28,
    "subheading_size": 20,
    "axis_label_size": 20,
    "marker_size": 20,
    "cds_marker_width": 2,
    "cds_marker_colour": '#000000',
    "legend_size": 17,
    "ribo_linewidth": 2
}
STARTS_CODONS = ['aug', 'cug', 'gug', 'none']
VARIABLE_UPPER = {"tranlist": str.upper}
VARIABLE_CONVERSION = {
    "cds_start": int,
    "cds_stop": int,
    "coverage": int,
    "exclude_first_val": int,
    "exclude_last_val": int,
    "heatmap_endpos": int,
    "heatmap_maxreadlen": int,
    "heatmap_minreadlen": int,
    "heatmap_startpos": int,
    "hili_end": int,
    "hili_start": int,
    "include_first_val": int,
    "include_last_val": int,
    "length": float,
    "mapped_reads": float,
    "max_coverage": float,
    "max_highest_frame_diff": float,
    "max_lowest_frame_diff": float,
    "max_start_increase": float,
    "max_stop_decrease": float,
    "maxread": int,
    "maxreadlen": int,
    "min_coverage": float,
    "min_highest_frame_diff": float,
    "min_lowest_frame_diff": float,
    "min_start_increase": float,
    "min_stop_decrease": float,
    "minimum_reads": int,
    "minread": int,
    "minreadlen": int,
    "minreads": float,
    "nuc_freq_plot_window": int,
    "nuc_maxreadlen": int,
    "nuc_minreadlen": int,
    "readscore": int,
    "rrna_removed": int,
    "start": int,
    "stop": int,
    "sw_diff_min_diff": int,
    "sw_diff_step_size": int,
    "sw_diff_window_size": int,
    "trip_maxreadlen": int,
    "trip_minreadlen": int,
    "unmapped_reads": int,
}
SCRIPT_LOC = "."
CDS_MARKER_COLOUR = "black"
RIBO_LINEWIDTH = 2
SQLITES_DIR = "trips_data_sample"
ANNOTATION_DIR = "trips_annotations_sample"
UPLOADS_DIR = SCRIPT_LOC + "/uploads/"
EMAIL_PASS = ""
LOG_FILE = "./log.txt"
DATABASE_NAME = "trips.sqlite"
