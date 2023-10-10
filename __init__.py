from typing import Union, Any
import pandas as pd
import os
import time
from datetime import date
import sys
import sqlite3
from sqlqueries import sqlquery, get_user_id, get_table, update_table
# import riboflask_datasets
import logging
from flask import (Flask, get_flashed_messages, render_template, request,
                   send_from_directory, flash, redirect, url_for)
from flask import Response
from flask_xcaptcha import XCaptcha

from flask_login import (LoginManager, login_required, login_user, logout_user,
                         current_user)
from flask_mail import Mail
from werkzeug.security import generate_password_hash, check_password_hash
import config
from werkzeug.utils import secure_filename
from sqlitedict import SqliteDict
import smtplib
import json
from core_functions import fetch_file_paths, base62_to_integer, User, fetch_user
from metainfo_routes import metainfo_plotpage_blueprint, metainfoquery_blueprint
from comparison_routes import comparison_plotpage_blueprint, comparisonquery_blueprint
from single_transcript_routes import (single_transcript_plotpage_blueprint,
                                      single_transcript_query_blueprint)
from single_transcript_routes_genomic import (
    single_transcript_plotpage_genomic_blueprint,
    single_transcript_query_genomic_blueprint)
from diff_exp_routes import diff_plotpage_blueprint, diffquery_blueprint
from pause_routes import pause_detection_blueprint, pausequery_blueprint
from traninfo_routes import traninfo_plotpage_blueprint, traninfoquery_blueprint
# , taskstatus_blueprint
from orfquery_routes import translated_orf_blueprint, orfquery_blueprint
import re

from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

root_logger = logging.getLogger()
# lhStdout = root_logger.handlers[0]
root_logger.setLevel(logging.DEBUG)
root_logger.propagate = False

log_format = logging.Formatter(
    "%(asctime)s [%(levelname)s] [%(name)s] %(message)s")

today = date.today()
file_handler = logging.FileHandler(config.LOG_FILE + str(today))
file_handler.setFormatter(log_format)

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(log_format)

root_logger.addHandler(file_handler)
# root_logger.addHandler(stream_handler)
# root_logger.removeHandler(lhStdout)

# logging.basicConfig(filename=config.LOG_FILE,level=logging.DEBUG,format='%(asctime)s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
app = Flask(__name__, static_folder='static')

logging.debug('This will get logged to a file')

user_short_passed = False

app.register_blueprint(metainfo_plotpage_blueprint)
app.register_blueprint(metainfoquery_blueprint)
app.register_blueprint(single_transcript_plotpage_blueprint)
app.register_blueprint(single_transcript_query_blueprint)
app.register_blueprint(single_transcript_plotpage_genomic_blueprint)
app.register_blueprint(single_transcript_query_genomic_blueprint)
app.register_blueprint(comparison_plotpage_blueprint)
app.register_blueprint(comparisonquery_blueprint)
app.register_blueprint(diff_plotpage_blueprint)
app.register_blueprint(diffquery_blueprint)
app.register_blueprint(translated_orf_blueprint)
app.register_blueprint(orfquery_blueprint)
app.register_blueprint(pause_detection_blueprint)
app.register_blueprint(pausequery_blueprint)
# app.register_blueprint(taskstatus_blueprint)
app.register_blueprint(traninfo_plotpage_blueprint)
app.register_blueprint(traninfoquery_blueprint)
app.config.from_pyfile('config.py')

xcaptcha = XCaptcha(app=app)

app.config['UPLOAD_FOLDER'] = '/static/tmp'
app.config['SECURITY_PASSWORD_SALT'] = config.PASSWORD_SALT
app.config['SECRET_KEY'] = config.FLASK_SECRET_KEY

# Modify security messages so people can't tell which users have already signed up.
app.config['SECURITY_MSG_EMAIL_ALREADY_ASSOCIATED'] = ((
    "Thank you. Confirmation instructions have been sent to %(email)s."),
                                                       "error")
app.config['USER_DOES_NOT_EXIST'] = (("Invalid credentials"), "error")
app.config['INVALID_PASSWORD'] = (("Invalid credentials"), "error")

app.config['MAIL_SERVER'] = 'smtp.gmail.com'
app.config['MAIL_PORT'] = 587
app.config['MAIL_USERNAME'] = 'ribopipe@gmail.com'
app.config['MAIL_PASSWORD'] = config.EMAIL_PASS
app.config['MAIL_USE_TLS'] = True
app.config['MAIL_USE_SSL'] = False

mail = Mail(app)


def sanitize_get_request(request: Any | None) -> Any | None:
    '''
    take a get request and remove any XSS attempts
    '''
    if isinstance(request, str):
        request = re.sub("<.*>", "", request)

    return request


# change cookie name and path, this avoids cookies clashing with other flask apps on the same server
try:
    if sys.argv[1] == "true":
        pass
except Exception:
    app.config.update(
        SESSION_REFRESH_EACH_REQUEST=False,
        SESSION_COOKIE_HTTPONLY=False,
    )

# flask-login
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = "login"
login_manager.login_message = None
login_manager.session_protection = 'basic'


# Provides statistics on trips such as number of organisms, number of files, number of studies etc and lists updates.
@app.route('/stats/')
def statisticspage() -> str:
    '''Display statistics'''

    def rename_organism(organism: str) -> str:
        '''
        Rename the organism name to the short name
        '''
        if '_' in organism:
            orgt = organism.split('_')
            organism = f"{orgt[0][0]}.{orgt[1]}"
        return organism.capitalize()

    organisms = get_table('organisms')
    organisms = organisms.loc[organisms.private == 0,
                              ['organism_id', 'organism_name']]
    organisms['organism_name'] = organisms['organism_name'].apply(
        rename_organism)

    no_organisms = len(organisms)

    files = get_table('files')[['organism_id', 'file_type']]
    riboseq_files = files[files.file_type == 'riboseq'].shape[0]
    rnaseq_files = files[files.file_type == 'rnaseq'].shape[0]
    files['file_type'] = files['file_type'].apply(
        lambda x: f"{x.capitalize()} files")
    org_files_count = organisms.merge(files, on='organism_id').groupby(
        ['organism_name',
         'file_type']).size().reset_index().rename(columns={
             0: 'Count',
             'organism_name': 'Organism'
         })

    # NOTE: Till here

    # Create the graph which breaksdown the number of studies per organism

    org_breakdown_graph = stats_plots.org_breakdown_plot(org_files_count)

    # Create the graph which breaks down studies published per year
    public_studies = get_table('studies')
    public_studies = public_studies[((public_studies.private == 0) and
                                     (~pd.isnull(public_studies.paper_year)))]
    no_studies = len(public_studies)
    year_dist = public_studies.groupby(
        ['paper_year']).size().reset_index().rename(columns={0: 'Count'})

    year_plot = stats_plots.year_dist(year_dist)
    updates = get_table('updates')
    updates = updates.sort_values(by='date',
                                  ascending=False).to_html(classes='')

    return render_template('statistics.html',
                           no_organisms=no_organisms,
                           no_studies=no_studies,
                           riboseq_files=riboseq_files,
                           rnaseq_files=rnaseq_files,
                           org_breakdown_graph=org_breakdown_graph,
                           year_plot=year_plot,
                           news_string=updates)


# Contact page
@app.route('/contactus/', methods=["GET", "POST"])
def contactus() -> str | Response:
    if request.method == "POST":
        if xcaptcha.verify():
            fromaddr = "ribopipe@gmail.com"
            toaddr = "tripsvizsite@gmail.com"
            msg = MIMEMultipart()
            msg['From'] = fromaddr
            msg['To'] = toaddr
            msg['Subject'] = request.form['subject']
            msg.attach(
                MIMEText("Name: {}\nEmail: {}\nMessage: {}".format(
                    request.form['name'], request.form['email'],
                    request.form['message'])))
            server = smtplib.SMTP('smtp.gmail.com', 587)
            server.starttls()
            server.login(fromaddr, config.EMAIL_PASS)
            text = msg.as_string()
            server.sendmail(fromaddr, toaddr, text)
            server.quit()
            flash("Message sent successfully")
        else:
            flash("Recaptcha failed")
            return redirect(url_for('contactus'))

    return render_template('contact.html')


# This is the page where users create a new login.
@app.route("/create", methods=["GET", "POST"])
def create() -> str | Response:
    # if user is already logged in then redirect to homepage
    if current_user.is_authenticated:
        return redirect("/")
    error = None
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        password2 = request.form['password2']
        if xcaptcha.verify():
            username_dict = {}
            logging.debug("Connecting to trips.sqlite")

            # Added by anmol
            users = get_table('users')
            max_user_id = users['user_id'].max()
            users = users[['username', 'password']]
            username_dict = users.set_index('username').to_dict()['password']

            logging.debug("Closing trips.sqlite connection")

            if username in username_dict:
                error = "Error: {} is already registered".format(username)
                return render_template('create.html', error=error)
            if password == "":
                error = "Password cannot be empty"
                return render_template('create.html', error=error)
            if password != password2:  # TODO: Add JS to test
                error = "Passwords do not match"
                return render_template('create.html', error=error)
            hashed_pass = generate_password_hash(password)
            user_id = max_user_id + 1
            # Add -1 to study access list, causes problems when adding study id's later if we don't
            # Last value is temp_user, set to 0 because this is not a temporary user (temporary users identified by uuid in session cookie only, no username or pw)
            # TODO: create code for updating tabel
            update_table(
                'users', {
                    'user_id': user_id,
                    'username': username,
                    'password': hashed_pass,
                    'f1': '-1',
                    'f2': '',
                    'f3': 0,
                    'f4': 0
                }, 'insert')  # TODO: Need to check what are f-fields
            update_table('user_settings', config.DEFAULT_USER_SETTINGS,
                         'insert')
            logging.debug("Closing trips.sqlite connection")
            return redirect("/")
        else:
            error = 'Invalid Captcha. Please try again.'
            return render_template('create.html', error=error)
    else:
        return render_template('create.html', error=error)


# Allows users to change some global settings such as plot background colour, title size, tick label size, etc.
@app.route('/settings/')
# @login_required
def settingspage() -> Response:

    user = fetch_user()[0]
    # If user is not logged in and has rejected cookies they cannot use this page, so redirect to the homepage.
    if not user:
        return redirect(
            url_for(
                'homepage',
                message=
                "To use the settings page you either need to be logged in or allow cookies. Click the cookie policy link at the top left of the page to allow cookies."
            ))
    # get user_id
    user_id = get_user_id(user)
    user_settings = get_table('user_settings')
    user_settings = user_settings[user_settings['user_id'] == user_id].to_dict(
        orient='records')[0]

    return render_template('settings.html', user_settings=user_settings)


# Allows users to download fasta files, as well as scripts needed to produce their own sqlite files
@app.route('/downloads/')
def downloadspage() -> str:

    organism_dict = {
        "Scripts": [
            "bam_to_sqlite.py", "tsv_to_sqlite.py",
            "create_annotation_sqlite.py",
            "create_transcriptomic_to_genomic_sqlite.py"
        ]
    }
    try:
        user = current_user.name
    except Exception:
        user = None
    dbpath = '{}/{}'.format(config.SCRIPT_LOC, config.DATABASE_NAME)
    organisms = sqlquery(dbpath, 'organisms')
    organisms.loc[organisms.private == 0, 'organism_name'].values

    for organism in organisms:
        organism_dict[organism] = []
    trips_annotation_dir = "{}/{}/".format(config.SCRIPT_LOC,
                                           config.ANNOTATION_DIR)
    for org in set(os.listdir(trips_annotation_dir)) & set(organisms):
        for filename in os.listdir(trips_annotation_dir + "/" + org):
            ext = os.path.splitext(filename)[-1]
            if (ext in [
                    "fa", "gtf"
            ]) or (ext == "sqlite" and
                   ("transcriptomic" in filename or org in filename)):
                organism_dict[org].append(filename)

    return render_template('downloads.html',
                           user=user,
                           organism_dict=organism_dict)


# Called when user downloads something from the downloads page
@app.route('/downloadquery', methods=['POST'])
def download_file() -> Response:
    organism = request.form["organism"]
    assembly = request.form["assembly"]
    return send_from_directory("{}/{}/{}".format(config.SCRIPT_LOC,
                                                 config.ANNOTATION_DIR,
                                                 organism),
                               assembly,
                               as_attachment=True)


# Allows users to upload their own sqlite files and transcriptomes.
@app.route('/uploads/')
# @login_required
def uploadspage() -> str:

    user, logged_in = fetch_user()
    # If user is not logged in and has rejected cookies they cannot use this page, so redirect to the homepage.
    if not user:
        return redirect(
            url_for(
                'homepage',
                message=
                "To use the uploads page you either need to be logged in or allow cookies. Click the cookie policy link at the top left of the page to allow cookies."
            ))
    if not logged_in:
        flash(
            "You are not logged in, uploaded data will only be kept for a period of one day."
        )

    user_id = get_user_id(user)
    organisms = get_table('organisms')
    organisms_t = organisms.loc[
        ~organisms.private | organisms.owner == user_id,
        ["organism_name", "transcriptome_list", "organism_id"]]

    organism_dict = table_to_dict(
        organisms_t)  # key: organism name, value: transcriptome_list
    org_id_dict = table_to_dict(
        organisms_t)  # key: organism_id, value: [organism, transcriptome]
    # -- Selected one
    organism_access = get_table('organism_access')
    organism_ids = organism_access[organism_access.user_id ==
                                   user_id].organism_id
    organisms_t = organisms[organisms.organism_id.isin(organism_ids)]
    organism_dict = {**organism_dict, **table_to_dict(organisms_t)}
    org_id_dict = {**org_id_dict, **table_to_dict(organisms_t)}

    study_dict = get_table('studies')
    study_dict = study_dict[study_dict.owner == user_id][[
        'study_id', 'study_name', 'organism_id'
    ]]
    organisms_t = organisms[organisms.organism_id.isin(
        set(study_dict.organism_id) - set(org_id_dict))]
    org_id_dict = {**org_id_dict, **table_to_dict(organisms_t)}

    study_dict[int(row[0])] = [  # TODO: Study dict format
        row[1].replace("_{}".format(user_id), "", 1), org_id_dict[row[2]][0],
        org_id_dict[row[2]][1], []
    ]

    transcriptome_dict = organisms[organism_access.user_id == user_id]
    transcriptome_dict = table_to_dict(
        transcriptome_dict
    )  # key: organism_id, value: [organism_name, transcriptome_list]
    study_access = get_table('study_access')
    study_access = study_access[study_access.study_id.isin(study_dict)]
    users = get_table('users')
    users = users[users.user_id.isin(study_access.user_id)]
    users = set(users.username) - set(study_dict.keys())
    # NOTE : Till here

    for study_id in study_dict:
        study_access = get_table('study_access')
        user_ids = study_access[study_access.study_id == study_id].user_id
        users = get_table('users')
        user_names = users[users.user_id.isin(user_ids)].username
        user_names = set(user_names) - set(study_dict[study_id][3])
        study_dict[study_id][3] += list(user_names)
    files = get_table('files')
    # [file_name,study_id,file_id,file_description]
    files = files[files.owner == user_id]
    studies = get_table('studies')
    studies = studies.merge(files, on='study_id')
    studies['study_name'] = studies['study_name'].apply(
        lambda x: x.replace("_{}".format(user_id), "", 1))

    # key: file_name, value: [study_name,file_id,file_description]
    file_dict = {}
    for _, row in studies.iterrows():  # TODO: Table 2 dict
        file_dict[row['file_name']] = [
            row['study_name'], row['file_id'], row['file_description']
        ]
    seq_dict = get_table('seq_rules')
    seq_dict = table_to_dict(seq_dict[seq_dict.user_id == user_id],
                             ['seq_name', 'frame_breakdown'])
    return render_template('uploads.html',
                           user=user,
                           organism_dict=organism_dict,
                           study_dict=study_dict,
                           transcriptome_dict=transcriptome_dict,
                           file_dict=file_dict,
                           seq_dict=seq_dict)


# Called when user uploads something on the uploads page
@app.route('/uploadquery', methods=['POST'])
# @login_required
def upload_file() -> Response:
    # uploaded_files = request.files.getlist("file")
    f = request.files["file"]
    print(request.form)
    print(request.files)
    user, logged_in = fetch_user()
    user_id = get_user_id(user)
    filename = f.filename
    # filename = secure_filename(f.filename)
    print(filename)
    file_ext = filename.split(".")[-1]
    if file_ext != "sqlite":
        flash("Error: File extension should be sqlite not {}".format(file_ext))
        return redirect("/uploads/")
    if user != "public":
        foldername = "{}_{}".format(
            request.form["foldername"].replace(" ", "_"), user_id)
    else:
        foldername = "{}".format(request.form["foldername"].replace(" ", "_"))
    organism = request.form["organism"]
    assembly = request.form["assembly"]
    filetype_radio = request.form["filetype"]
    filetype = None  # TODO: remove it later
    if filetype_radio == "riboseq":
        filetype = "riboseq"
    elif filetype_radio == "rnaseq":
        filetype = "rnaseq"
    elif filetype_radio == "other":  # TODO: check it for correct type
        filetype = (request.form["seq_type"]).lower().strip()

    # if this filetype is new for this user insert a new entry into seq_rules table
    if filetype not in ["riboseq", "rnaseq"]:
        seq_rule = get_table('seq_rules')
        seq_rule = seq_rule[seq_rule.user_id == user_id,
                            seq_rule.seq_name == filetype]
        if seq_rule.empty:
            update_table('seq_rules', {
                'user_id': user_id,
                'seq_name': filetype,
                'f1': 0
            }, 'insert')  # TODO: What is f1

    if not os.path.isdir("{}/uploads/{}".format(config.SCRIPT_LOC,
                                                foldername)):
        os.makedirs("{}/uploads/{}".format(config.SCRIPT_LOC, foldername))
    upload_file_path = "{}/uploads/{}/{}".format(config.SCRIPT_LOC, foldername,
                                                 filename)
    with open(upload_file_path, "wb") as fout:
        while True:
            chunk = f.stream.read(1024)
            if not chunk:
                break
            fout.write(chunk)

    # f.save("{}/uploads/{}/{}".format(config.SCRIPT_LOC, foldername,
    # filename))
    sqlite_db = SqliteDict("{}/uploads/{}/{}".format(config.SCRIPT_LOC,
                                                     foldername, filename))
    try:
        file_description = sqlite_db["description"]
    except Exception:
        file_description = None
    sqlite_db.close()
    # get file id
    max_file_id = get_table('files').file_id.max()
    new_file_id = max_file_id + 1
    organism_id = get_table('organisms')
    organism_id = organism_id.loc[organism_id.organism_name == organism
                                  & organism_id.transcriptome_list == assembly,
                                  'organism_id'].values[0]
    studies = get_table('studies')
    studies = studies.loc[studies.organism_id == organism_id,
                          studies.owner == user_id,
                          studies.study_name == foldername, "study_id"]

    if not studies.empty:
        study_id = studies.values[0]
    else:
        cursor.execute("SELECT MAX(study_id) FROM studies;")
        result = cursor.fetchone()
        study_id = int(result[0]) + 1
        if user != "public":
            cursor.execute(
                "INSERT INTO studies VALUES({},{},'{}','{}','{}','{}','{}','{}','{}','{}','{}','{}',{},{})"
                .format(study_id, organism_id, foldername, 'NULL', 'NULL',
                        'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL',
                        1, user_id))
        else:
            cursor.execute(
                "INSERT INTO studies VALUES({},{},'{}','{}','{}','{}','{}','{}','{}','{}','{}','{}',{},{})"
                .format(study_id, organism_id, foldername, 'NULL', 'NULL',
                        'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL',
                        0, user_id))
        cursor.execute("INSERT INTO study_access VALUES({},{});".format(
            study_id, user_id))
    cursor.execute(
        "INSERT INTO files VALUES({},{},{},'{}','{}','{}',{},{},{},'{}')".
        format(new_file_id, organism_id, study_id, filename, file_description,
               filetype, user_id, 0, 0, ""))
    # If user is not logged in only keep file for a set period of time
    if not logged_in:
        curr_time = time.time()
        # The time to keep the file in seconds, currently set to one day
        keep_time = 60 * 60 * 24
        deletion_time = curr_time + keep_time
        cursor.execute("INSERT INTO deletions VALUES({},'{}',{})".format(
            new_file_id, upload_file_path, deletion_time))
    connection.commit()
    flash("File uploaded successfully")
    connection.close()
    return redirect("/uploads/")


# Called when a user uploads a custom transcriptome
@app.route('/uploadtranscriptome', methods=['GET', 'POST'])
# @login_required
def upload_transcriptome() -> Union[str, Response, None]:
    user, logged_in = fetch_user()
    user_id = get_user_id(user)
    if request.method == 'POST':
        organism = (request.form["organism"]).lower().strip().replace(" ", "_")
        assembly = (request.form["assembly"]).lower().strip().replace(" ", "_")
        default_tran = (request.form["default_tran"]).lower().strip().replace(
            " ", "_")
        uploaded_annotation = request.files.getlist("anno_file")
        fold_path = "{}/uploads/transcriptomes/{}/{}/{}".format(
            config.SCRIPT_LOC, user_id, organism, assembly)
        if not os.path.isdir(fold_path):
            os.makedirs(fold_path, exist_ok=True)
        for f in uploaded_annotation:
            filename = secure_filename(f.filename)
            ext = filename.split(".")[-1]
            if ext != "sqlite":
                return """Error: Expecting extension sqlite but got extension {}. The file generated by the create_annotation_sqlite.py script should be uploaded here.
						This script can be gotten on the downloads page, by selecting the Scripts group.""".format(
                    ext)
            # Instead of using filename of the uploaded file we rename it to organism_assembly.sqlite, to keep things consistent
            filename = "{}_{}.sqlite".format(organism, assembly)
            full_path = "{}/{}".format(fold_path, filename)
            f.save(full_path)
            max_org_id = get_table('organisms').organism_id.max() + 1
            cursor.execute(
                "INSERT INTO organisms VALUES({},'{}','{}','NULL','NULL','NULL','NULL','{}',1,{})"
                .format(max_org_id, organism, assembly, default_tran.upper(),
                        user_id))
            cursor.execute("INSERT INTO organism_access VALUES ({},{})".format(
                max_org_id, user_id))
            if not logged_in:
                # Instead of deleting the file now, add it to deletions table where it will be deleted via cron job, this will give users time to contact in case of accidental deletion
                curr_time = time.time()
                # The time to keep the file in seconds, currently set to 1 day
                keep_time = 60 * 60 * 24
                deletion_time = curr_time + keep_time
                cursor.execute(
                    "INSERT INTO org_deletions VALUES({},'{}',{})".format(
                        max_org_id, full_path, deletion_time))
            connection.commit()
            connection.close()
        flash("File uploaded successfully")
        return redirect("/uploads/")


# Called by flask in case of an error in the code, returns the exception so it can be displayed to user
@app.errorhandler(500)
def handle_bad_request(e: Exception) -> str:
    return_str = 'ERROR: ' + str(
        e
    ) + " please report this to tripsvizsite@gmail.com or via the contact page. "
    return return_str


# This is the page where users login.
@app.route("/user/login", methods=["GET", "POST"])
def login() -> Union[str, Response]:
    # if user is already logged in then redirect to homepage
    if current_user.is_authenticated:
        return redirect("/")
    error = None
    if request.method == 'POST':
        username = request.form['username'].strip(
        )  # TODO: relace this with jquery.Trim on user side
        password = request.form['password'].strip()
        if xcaptcha.verify() or username == "developer":
            logging.debug("login Connecting to trips.sqlite")
            username_dict = table_to_dict(get_table('users'),
                                          ['username', 'password'])
            logging.debug("Closing trips.sqlite connection")
            if username in username_dict:
                if check_password_hash(username_dict[username], password):
                    login_user(User(username))
                    nxt = sanitize_get_request(request.args.get('next'))
                    if nxt:
                        if "<function login" in nxt:
                            nxt = "/"
                    else:
                        nxt = "/"
                    return redirect(nxt)
                else:
                    error = 'Either username or password incorrect. Please try again.'
                    return render_template('login.html', error=error)
            else:
                error = 'Either username or password incorrect. Please try again.'
                return render_template('login.html', error=error)
        else:
            error = 'Invalid Captcha. Please try again.'
            return render_template('login.html', error=error)
    else:
        return render_template('login.html', error=error)


# Allows users to logout
@app.route("/user/logout")
@login_required
def logout() -> Response:
    logout_user()
    return redirect(url_for('homepage'))


# callback to reload the user object
@login_manager.user_loader
def load_user(userid):
    return User(userid)


# Called when user presses the save button on the orf_translation page.
@app.route('/anno_query', methods=['POST'])
def anno_query() -> str:
    data = json.loads(request.data)
    user = fetch_user()[0]
    data["user_id"] = get_user_id(user)
    for key in [
            'START_CODON', 'CDS_OVERLAP', 'START_SCORE', 'STOP_SCORE',
            'ENTROPY', 'TE', 'COVERAGE', 'CDS_RATIO', 'FILE_LIST'
    ]:
        data[key.lower()] = key
    update_table('users_saved_cases', 'insert', data)

    return ""


# This page shows the saved ORFs specific to the signed in user
@app.route('/saved/')
def saved():

    advanced = False
    user, logged_in = fetch_user()
    # If user is not logged in and has rejected cookies they cannot use this page, so redirect to the homepage.
    if not user:
        flash(
            "To use the Saved ORFs page you either need to be logged in or allow cookies."
        )
        return redirect(url_for('/'))

    user_id = -1
    if user and logged_in:
        flash("You are logged in as {}".format(user))
        user_id = get_user_id(user)
    organism_access = get_table('organism_access')
    organism_access_list = organism_access.loc[organism_access.user_id ==
                                               user_id, 'organism_id'].values
    organism_list = get_table('organisms')
    organism_list = organism_list.loc[
        ~organism_list.private
        | organism_list.organism_id.isin(organism_access_list),
        'organism_id'].values
    return render_template('user_saved_cases.html',
                           advanced=advanced,
                           organism_list=organism_list)


# Retrieves saved ORFs
@app.route('/savedquery', methods=['POST'])
def savedquery():
    data = json.loads(request.data)
    user = fetch_user()[0]
    organism = data["organism"]
    label = data["label"]
    # get user_id
    user_id = get_user_id(user)
    # structure of orf dict is transcript[stop][start] = {"length":x,"score":0,"cds_cov":0} each stop can have multiple starts
    user_saved_cases = get_table('users_saved_cases')
    user_saved_cases = user_saved_cases.loc[user_saved_cases.user_id ==
                                            user_id]
    if organism != 'Select an Organism':
        user_saved_cases = user_saved_cases.loc[user_saved_cases.organism ==
                                                organism]
        if label:
            label_list = label.strip().split(",")
            user_saved_cases = user_saved_cases.loc[
                user_saved_cases.label.isin(label_list)]

    user_saved_cases = user_saved_cases.head(1000)
    returnstr = user_saved_cases.apply(
        lambda x:
        f"{x['gene']},{x['transcript']},{x['start']},{x['stop']},{x['length']},{x['score']},{x['label']},{x['start_codon']},{x['trips_link']}",
        axis=1).values
    returnstr = ".,/".join(returnstr)
    return returnstr


# Allows users to delete previously saved cases
@app.route('/del_query', methods=['POST'])
def del_query():
    data = json.loads(request.data)
    try:
        user = current_user.name
    except Exception:
        return "Error user not signed in"
    user_id = get_user_id(user)
    data["user_id"] = user_id
    update_table('users_saved_cases', 'delete', {}, data)
    return ""


# Allows users to logout
# @app.route("/user/logout")
# @login_required
# def logout():
# logout_user()
# return redirect(login)


# Points to robots.txt in static folder
@app.route('/robots.txt')
def static_from_root():
    return send_from_directory(app.static_folder, request.path[1:])


# Points to license.txt in static folder
@app.route('/license.txt')
def static_license_from_root():
    return send_from_directory(app.static_folder, request.path[1:])


# This is the help page, linked from various other pages to explain terms on that page.
@app.route('/help/')
def helppage():
    parent_acc = sanitize_get_request(request.args.get('parent_acc'))
    child_acc = sanitize_get_request(request.args.get('child_acc'))
    logging.debug(type(parent_acc))
    logging.debug(child_acc)
    return render_template('help.html',
                           parent_acc=parent_acc,
                           child_acc=child_acc)


# This is the help page, linked from various other pages to explain terms on that page.
@app.route('/shared/<folder>')
def sharepage(folder):
    filelist = []
    for filename in os.listdir("{}/shared/{}".format(config.SCRIPT_LOC,
                                                     folder)):
        if filename.endswith(".html"):
            filelist.append(filename.replace(".html", ""))
    filelist = sorted(filelist)
    return render_template('shared.html', filelist=filelist, folder=folder)


# Called when user downloads something from the downloads page
@app.route('/shared/<folder>/<filename>', methods=['GET', 'POST'])
def viewfile(folder, filename):
    filename = filename + ".html"
    return send_from_directory("{}/shared/{}/".format(config.SCRIPT_LOC,
                                                      folder),
                               filename,
                               as_attachment=False)


# This is the short url page, user supplies a short code which will be converted to a full url which user will then be redirected to
@app.route('/short/<short_code>/')
def short(short_code):
    # First convert short code to an integer
    integer = base62_to_integer(short_code)
    url = get_table('urls')
    url = url[url.url_id == integer]
    if url.empty:
        return "Short code not recognized."
    url = url.url.values[0]
    # add a keyword to the url to prevent generating another shortcode
    url += "&short={}".format(short_code)
    return redirect(url)


@app.after_request
def after_request_func(response):
    request.cookies.get("cookieconsent_status")
    return response


# This is the home page it show a list of organisms as defined by trips_dict
# TODO: Move pages in order
@app.route('/')
def homepage2() -> str:
    """Home page."""
    sanitize_get_request(request.cookies.get("cookieconsent_status"))

    # For All
    organisms = get_table('organisms')

    # user related details
    user, logged_in = fetch_user()
    private_organisms_id = []
    if logged_in:  # Only if loggned in
        # TODO: Replace login and register part with username
        flash(f"You are logged in as {user}")
        user_id = get_user_id(user)  # TODO: Find a way to pass only str
        organism_access = get_table('organism_access')
        private_organisms_id = organism_access.loc[
            organism_access.user_id == user_id, 'organism_id'].values.tolist()
    organisms = organisms.loc[
        ~organisms.private |
        (organisms.organism_id.isin(private_organisms_id)),
        ["organism_name", "transcriptome_list"]].drop_duplicates()
    organisms = organisms.groupby("organism_name")['transcriptome_list'].apply(
        list).reset_index()

    print("Anmol", organisms)

    # Create species list

    return render_template('landing2.html', organisms=organisms, message="")


def homepage(message=""):
    organism_access_list = []
    logging.debug("homepage Connecting to trips.sqlite")
    sanitize_get_request(request.cookies.get("cookieconsent_status"))
    user, logged_in = fetch_user()

    user_id = -1
    if user:
        if logged_in:
            flash("You are logged in as {}".format(user))
        user_id = get_user_id(user)
        organism_access_list = get_table('organism_access')
        organism_access_list = organism_access_list.loc[
            organism_access_list.user_id == user_id, 'organism_id'].values
        # get a list of organism id's this user can access

    # returns a tuple with each field as a seperate string
    organisms = get_table('organisms')
    orgsnisms = organisms.loc[(
        ~organisms.private | organisms.organism_id.isin(organism_access_list)
        & organisms.organism_name == organism),
                              ['organism_name', 'transcriptome_list']]
    organisms = organisms.groupby("organism_name")['transcriptome_list'].apply(
        list).reset_index()
    # organism_list.sort()

    logging.debug("homepage Closing trips.sqlite connection")
    print(organisms)

    message = sanitize_get_request(request.args.get('message'))
    return render_template('landing.html',
                           organisms=orgsnisms,
                           message=message)


# Updates the settings for a specific user
@app.route('/settingsquery', methods=['POST'])
# @login_required
def settingsquery():
    data = json.loads(request.data)
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    connection.text_factory = str
    cursor = connection.cursor()

    user, logged_in = fetch_user()
    if logged_in:
        new_password = data['new_password']
        new_password2 = data['new_password2']
        curr_password = data['curr_password']
        if new_password and new_password2:
            if new_password != new_password2:
                return "ERROR: New passwords do not match"
            generate_password_hash(curr_password)
            users = get_table('users')
            old_password_hash = users.loc[users.username == user,
                                          'password'].values[0]

            if check_password_hash(old_password_hash, curr_password):
                new_password_hash = generate_password_hash(new_password)
                cursor.execute(
                    "UPDATE users SET password = '{}' WHERE username = '{}'".
                    format(new_password_hash, user))
                connection.commit()
            else:
                return "ERROR: Current password is not correct"
    # get user_id
    user_id = get_user_id(user)
    update_table('users', 'update', {'user_id': user_id},
                 {'advanced': 1 if 'advanced' in data else 0})
    update_table('user_settings', 'update', {'user_id': user_id},
                 data)  # TODO: might need some pruning

    return "Settings have been updated"


# Allows users to delete files
@app.route('/deletequery', methods=['GET', 'POST'])
# @login_required
def deletequery():
    data = json.loads(request.data)

    user = fetch_user()[0]
    user_id = get_user_id(user)
    files_all = get_table("files")
    studies_all = get_table("studies")
    organisms_all = get_table("organisms")
    for key in data:
        file_id = data[key]["file_id"]
        if "filecheck" in data[key]:
            files = files_all[files_all["file_id"] == file_id].iloc[0]
            studies = studies_all[studies_all["study_id"] ==
                                  files.study_id].iloc[0]
            full_path = "{}{}/{}".format(config.UPLOADS_DIR,
                                         studies.study_name, files.filename)
            # Instead of deleting the file now, add it to deletions table where it will be deleted via cron job, this will give users time to contact in case of accidental deletion
            curr_time = time.time()
            # The time to keep the file in seconds, currently set to 14 days
            keep_time = 60 * 60 * 24 * 14
            deletion_time = curr_time + keep_time
            cursor.execute("INSERT INTO deletions VALUES({},'{}',{})".format(
                file_id, full_path, deletion_time))
            if owner == user_id:
                cursor.execute(
                    "DELETE FROM files WHERE file_id = {}".format(file_id))

        cursor.execute(
            "UPDATE files SET file_description = '{}' WHERE file_id = {}".
            format(data[key]["file_desc"], file_id))
        if data[key]["cutadapt_removed"] != '0':
            organism_id = files_all.loc[files_all["file_id"] == file_id,
                                        "organism_id"].values[0]
            organism = organisms_all.loc[organisms_all["organism_id"] ==
                                         organism_id,
                                         "organism_name"].values[0]
            filepath_dict = fetch_file_paths([file_id], organism)
            for seq_type in filepath_dict:
                if file_id in filepath_dict[seq_type]:
                    filepath = filepath_dict[seq_type][file_id]
                    opendict = SqliteDict(filepath, autocommit=True)
                    opendict["cutadapt_removed"] = int(
                        data[key]["cutadapt_removed"])
                    opendict.close()

        if data[key]["rrna_removed"] != '0':
            organism_id = files_all.loc[files.file_id == file_id,
                                        'organism_id'].values[0]
            organism = organism_all.loc[organism_all.organism_id ==
                                        organism_id, 'organism_name'].values[0]
            filepath_dict = fetch_file_paths([file_id], organism)
            for seq_type in filepath_dict:
                if file_id in filepath_dict[seq_type]:
                    filepath = filepath_dict[seq_type][file_id]
                    opendict = SqliteDict(filepath, autocommit=True)
                    opendict["rrna_removed"] = int(data[key]["rrna_removed"])
                    opendict.close()
        if data[key]["unmapped"] != '0':
            organism_id = files_all.loc[files.file_id == file_id,
                                        'organism_id'].values[0]
            organism = organism_all.loc[organism_all.organism_id ==
                                        organism_id, 'organism_name'].values[0]

            filepath_dict = fetch_file_paths([file_id], organism)
            for seq_type in filepath_dict:
                if file_id in filepath_dict[seq_type]:
                    filepath = filepath_dict[seq_type][file_id]
                    opendict = SqliteDict(filepath, autocommit=True)
                    opendict["unmapped_reads"] = int(data[key]["unmapped"])
                    opendict.close()
    flash("File list updated")
    return redirect("https://trips.ucc.ie/uploads")


# Allows users to delete studies,modify access, modify the organism/transcriptome assembly or study name
@app.route('/deletestudyquery', methods=['GET', 'POST'])
# @login_required
def deletestudyquery():
    data = json.loads(request.data)
    user = fetch_user()[0]
    user_id = get_user_id(user)
    for study_id in data:
        studycheck = data[study_id][0]
        # Delete studies where the "delete" checkbox is checked

        if studycheck.split("_")[-1] != "undefined":
            study_id = studycheck.split("_")[-1]
            print("deleting study_id", study_id)
            files = get_table('files')
            files = files[files.study_id == study_id]
            studies = get_table('studies')
            studies_files = studies.merge(files, on='study_id')
            if not studies_files.empty:
                studies_files = studies_files[studies_files.owner == user_id]
                studies_files['full_path'] = stydies_files.apply(
                    lambda x: "{}{}/{}".format(config.UPLOADS_DIR, x.
                                               study_name, x.filename),
                    axis=1)
                studies_files = studies_files[studies_files['full_path'].apply(
                    os.path.isfile)]
                for _, row in studies_files.iterrows():
                    curr_time = time.time()
                    keep_time = 60 * 60 * 24 * 14
                    deletion_time = curr_time + keep_time
                    # TODO: Use cron for deletion
                    update_table('deletions', 'insert', {}, {
                        'file_id': row['file_id'],
                        'full_path': row['full_path'],
                        'deletion_time': deletion_time
                    })

            # Now remove the study and the files associated with it from the db
            study_owner = studies[studies.study_id == study_id].iloc[0].owner
            if study_owner == user_id:
                # os.rename(study_path,study_path+"_REMOVE")
                update_table('studies', {'study_id': study_id})
                update_table('files', {'study_id': study_id})
                continue

        # Modify access list next to studies
        study_access = data[study_id][1].split(",")
        study_access = list(filter(('').__ne__, study_access))
        # check study_access against a list of all users
        all_users = table_to_dict(get_table('users'), ['username', 'user_id'])
        # Check that all users exist
        for username in study_access:
            if username:
                if username not in all_users.keys():
                    flash(
                        "Error: User {} is not registered on Trips-Viz".format(
                            username))
                    return str(get_flashed_messages())
                else:
                    study_access = get_table('study_access')
                    study_access = study_access[
                        study_access.user_id == all_users[username]
                        & study_access.study_id == study_id]
                    if study_access.empty:
                        update_table('study_access', 'insert', {}, {
                            'study_id': study_id,
                            'user_id': all_users[username]
                        })

        # Modify study names if they have changed
        # TODO: Till here
        new_study_name = "{}_{}".format(data[study_id][2], user_id)
        cursor.execute(
            "SELECT study_name FROM studies WHERE study_id = {}".format(
                study_id))
        old_study_name = cursor.fetchone()[0]
        if old_study_name != new_study_name:
            # Update study name in the sqlite
            cursor.execute(
                "UPDATE studies SET study_name = '{}' WHERE study_id = {}".
                format(new_study_name, study_id))
            # If the new_study_name folder does not exist, rename the old study to the new study, else move all files from old folder to new folder
            if not os.path.isdir("{}/uploads/{}".format(
                    config.SCRIPT_LOC, new_study_name)):
                os.rename(
                    "{0}/uploads/{1}".format(config.SCRIPT_LOC,
                                             old_study_name),
                    "{0}/uploads/{1}".format(config.SCRIPT_LOC,
                                             new_study_name))
            else:
                if os.path.isdir("{}/uploads/{}".format(
                        config.SCRIPT_LOC, old_study_name)):
                    for filename in os.listdir("{}/uploads/{}".format(
                            config.SCRIPT_LOC, new_study_name)):
                        if os.path.isfile("{}/uploads/{}/{}".format(
                                config.SCRIPT_LOC, old_study_name, filename)):
                            os.rename(
                                "{}/uploads/{}/{}".format(
                                    config.SCRIPT_LOC, old_study_name,
                                    filename), "{}/uploads/{}/{}".format(
                                        config.SCRIPT_LOC, new_study_name,
                                        filename))
        # Change organism/transcriptome assembly if applicable
        organism_name = data[study_id][3]
        assembly_name = data[study_id][4]
        cursor.execute(
            "SELECT organism_id FROM studies WHERE study_id = {}".format(
                study_id))
        org_id = cursor.fetchone()[0]
        cursor.execute(
            "SELECT organism_name,transcriptome_list FROM organisms WHERE organism_id = {}"
            .format(org_id))
        result = cursor.fetchone()
        if result != None:
            old_organism = result[0]
            old_assembly = result[1]
            if old_organism != organism_name or old_assembly != assembly_name:
                # Check if the new orgnaism and new assembly are a valid combination
                cursor.execute(
                    "SELECT organism_id  FROM organisms WHERE organism_name = '{}' AND transcriptome_list = '{}'"
                    .format(organism_name, assembly_name))
                result = cursor.fetchone()
                if result == None:
                    return "Invalid organism/transcriptome combo for study {}".format(
                        new_study_name)
                else:
                    cursor.execute(
                        "UPDATE studies SET organism_id = {} WHERE study_id = {}"
                        .format(result[0], study_id))
                    pass
                    # update study_id with new org_id

    connection.commit()
    connection.close()
    flash("Update successful")
    return redirect("https://trips.ucc.ie/uploads")


# Allows users to delete transcriptomes
@app.route('/deletetranscriptomequery', methods=['GET', 'POST'])
# @login_required
def deletetranscriptomequery():
    data = json.loads(request.data)

    user = fetch_user()[0]
    organisms = get_table("organisms")
    organism_ids = [
        val[0].split("_")[-1] for _, val in data.items()
        if not val[0].endswith("_undefined")
    ]
    user_id = get_user_id(user)
    organisms = organisms.loc[
        organisms.organism_id.isin(organism_ids) &
        (organisms.owner == user_id),
        ['organism_name', 'transcriptome_list']].drop_duplicates()
    keep_time = 60 * 60 * 24 * 14  # 14 days
    curr_time = time.time()
    deletion_time = curr_time + keep_time
    organisms['sqlite_path'] = organisms.apply(
        lambda x: "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
            config.UPLOADS_DIR, user_id, x['organism_name'], x[
                'transcriptome_list']),
        axis=1)
    organisms = organisms[organisms['sqlite_path'].apply(os.path.isfile)]
    organisms.apply(lambda x: update_table("org_deletions", 'insert', {},
                                           {'file_id': x['file_id']}),
                    axis=1)

    for organism_id in data:
        organism_id = data[organism_id][0].split("_")[-1]
        if organism_id == "undefined":
            continue
        # Delete the annotation sqlite file
        organisms = get_table("organisms")
        organisms = organisms[organisms.organism_id == organism_id
                              & organisms.owner == user_id].iloc[0][
                                  "organism_name", "transcriptome_list"]
        sqlite_path = "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
            config.UPLOADS_DIR, user_id, organisms.organism_name,
            organisms.transcriptome_list)
        if os.path.isfile(sqlite_path):
            # Instead of deleting the file now, add it to deletions table where it will be deleted via cron job, this will give users time to contact in case of accidental deletion
            curr_time = time.time()
            # The time to keep the file in seconds, currently set to 14 days
            keep_time = 60 * 60 * 24 * 14
            deletion_time = curr_time + keep_time
            update_table("org_deletions", 'insert', {}, {
                "organism_id": organism_id,
                "deletion_time": deletion_time,
                "sqlite_path": sqlite_path
            })
        # sqlite_dir = "{0}transcriptomes/{1}/{2}/{3}".format(config.UPLOADS_DIR, user_id, organism_name,transcriptome_list)
        # if os.path.isdir(sqlite_dir):
        # os.rename(sqlite_dir,sqlite_dir+"_REMOVE")
        update_table("organisms", {"organism_id": organism_id})
        files = get_table("files")
        files = files[files.organism_id == organism_id]
        studies = get_table("studies")
        files_studies = files.merge(studies, on="study_id")
        study_ids = files.study_id.unique()
        study_names = files_studies.study_name.unique()

        # delete all files on the server associated with this organism, if there are any
        for _, row in files_studies.iterrows():
            full_path = "{}{}/{}".format(config.UPLOADS_DIR, row.study_name,
                                         row.file_name)
            if os.path.isfile(full_path):
                os.remove(full_path)

        for study_name in study_names:
            full_path = "{}{}".format(config.UPLOADS_DIR, study_name)
            if os.path.isdir(full_path):
                os.rename(full_path, full_path + "_REMOVE")

        # Now remove the study and the files associated with it from the db
        update_table("studies", {"study_id": study_ids})
        update_table("files", {"study_id": study_ids})
    flash("Update successful")
    return redirect("https://trips.ucc.ie/uploads")


# Updates the "sequence rules", for custom sequence types
@app.route('/seqrulesquery', methods=['GET', 'POST'])
# @login_required
def seqrulesquery():
    data = json.loads(request.data)
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    connection.text_factory = str
    cursor = connection.cursor()

    user = fetch_user()[0]

    # get user_id
    cursor.execute(
        "SELECT user_id from users WHERE username = '{}';".format(user))
    result = (cursor.fetchone())
    user_id = result[0]
    for seq_type in data:
        if data[seq_type][0] == 'False':
            cursor.execute(
                "UPDATE seq_rules SET frame_breakdown = 0 WHERE seq_name = '{}' and user_id = {};"
                .format(seq_type, user_id))
        elif data[seq_type][0] == 'True':
            cursor.execute(
                "UPDATE seq_rules SET frame_breakdown = 1 WHERE seq_name = '{}' and user_id = {};"
                .format(seq_type, user_id))
    connection.commit()
    connection.close()
    flash("Update successful")
    return redirect("https://trips.ucc.ie/uploads")


# breaks down the counts from each file for a specific ORF
@app.route('/<organism>/<transcriptome>/dataset_breakdown/')
def dataset_breakdown(organism, transcriptome):
    # ip = request.environ['REMOTE_ADDR']

    organism = str(organism)
    print(organism, transcriptome)
    accepted_orftype = sanitize_get_request(request.args.get("region"))
    transcript = sanitize_get_request(request.args.get("transcript"))
    start = int(sanitize_get_request(request.args.get("start")))
    stop = int(sanitize_get_request(request.args.get("stop")))
    files = sanitize_get_request(request.args.get("files")).strip(",")
    file_list = []
    for file_id in files.split("_"):
        file_list.append(int(file_id))
    file_paths_dict = fetch_file_paths(file_list, organism)
    xlist = []
    ylist = []
    filenames = []
    file_descs = []
    studies = []
    raw_reads = []
    controls = []
    cell_lines = []
    control_colors = []
    study_colors = []
    cell_line_colors = []
    i = 0
    cell_color_dict = {
        "HEK293": "#e00025",
        "HeLa": "#150ed1",
        "BJ fibroblast": "#0cfffa",
        "fibroblast": "#00c6d1",
        "MCF10A-ER-Src": "#00d100",
        "HCT116": "#ffe900",
        "U2OS": "#ffc042"
    }
    orfquery_connection = sqlite3.connect('{}/{}'.format(
        config.SCRIPT_LOC, config.DATABASE_NAME))
    orfquery_cursor = orfquery_connection.cursor()
    for file_id in file_paths_dict["riboseq"]:
        i += 1
        table_name = accepted_orftype
        xlist.append(i)
        sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
        try:
            raw_count = sqlite_db[table_name][transcript][stop][start][
                "profile_count"]
            raw_reads.append(raw_count)
        except Exception:
            raw_count = 0
            raw_reads.append(raw_count)
        orfquery_cursor.execute(
            "SELECT study_id,file_name,file_description from files WHERE file_id = {};"
            .format(file_id))
        result = orfquery_cursor.fetchone()
        study_id = int(result[0])
        filenames.append(result[1].replace(".sqlite", ""))
        file_descs.append(result[2])
        controls.append("True")
        control_colors.append("#4ef46f")
        cell_lines.append("HEK293")
        cell_line_colors.append(cell_color_dict["HEK293"])

        ylist.append(raw_count / (float(10000) / 100000))
        orfquery_cursor.execute(
            "SELECT * from studies WHERE study_id = {};".format(study_id))
        result = orfquery_cursor.fetchone()
        studies.append(result[0])
        study_colors.append('#BABABA')
    orfquery_cursor.close()
    orfquery_connection.close()
    return
    # return riboflask_datasets.generate_plot(xlist, ylist, filenames,
    # file_descs, studies, raw_reads,
    # controls, cell_lines,
    # control_colors, study_colors,
    # cell_line_colors, transcript,
    # start, stop)


if __name__ == '__main__':
    local = False
    try:
        if sys.argv[1] == "true":
            local = True
    except Exception:
        pass
    try:
        port_no = int(sys.argv[2])
    except Exception:
        port_no = 5000
    if local == False:
        app.run(host='0.0.0.0', debug=False)
    else:
        app.run(host='0.0.0.0', port=port_no, debug=True, threaded=True)
