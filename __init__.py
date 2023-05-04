import pandas as pd
import os
import time
from datetime import date
import sys
import sqlite3
from sqlqueries import sqlquery, get_user_id, get_table
import riboflask_datasets
import logging
from flask import (Flask, get_flashed_messages, render_template, request,
                   send_from_directory, flash, redirect, url_for)
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
import stats_plots
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


def sanitize_get_request(request):
    '''
	take a get request and remove any XSS attempts
	'''
    if isinstance(request, str):
        cleaned_request = re.sub("<.*>", "", request)

        return cleaned_request
    else:
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
def statisticspage():
    '''Display statistics'''

    def rename_organism(organism):
        '''
        Rename the organism name to the short name
        '''
        if '_' in organism:
            organism = organism.split('_')
            organism = f"{organism[0][0]}.{organism[1]}"
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
def contactus():
    if request.method == "POST":
        if xcaptcha.verify():
            fromaddr = "ribopipe@gmail.com"
            toaddr = "tripsvizsite@gmail.com"
            msg = MIMEMultipart()
            msg['From'] = fromaddr
            msg['To'] = toaddr
            msg['Subject'] = str(request.form['subject'])
            msg.attach(
                MIMEText("Name: {}\nEmail: {}\nMessage: {}".format(
                    str(request.form['name']), str(request.form['email']),
                    str(request.form['message']))))
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
def create():
    # if user is already logged in then redirect to homepage
    if current_user.is_authenticated:
        return redirect("/")
    error = None
    if request.method == 'POST':
        username = str(request.form['username'])
        password = str(request.form['password'])
        password2 = str(request.form['password2'])
        if xcaptcha.verify() or local:
            username_dict = {}
            logging.debug("Connecting to trips.sqlite")

            # Added by anmol
            dbpath = '{}/trips.sqlite'.format(config.SCRIPT_LOC)
            users = sqlquery(dbpath, 'users')
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
            cursor.execute(
                "INSERT INTO users VALUES ({},'{}','{}','-1','',0,0);".format(
                    user_id, username, hashed_pass))
            cursor.execute(
                "INSERT INTO user_settings VALUES ('{0}','{1}','{2}','{3}',{4},'{5}','{6}','{7}','{8}','{9}','{10}','{11}','{12}','{13}','{14}','{15}','{16}','{17}','{18}',{19},'{20}',{21},{22});"
                .format(config.MARKER_SIZE, config.AXIS_LABEL_SIZE,
                        config.SUBHEADING_SIZE, config.TITLE_SIZE, user_id,
                        config.BACKGROUND_COL, config.READLENGTH_COL,
                        config.METAGENE_FIVEPRIME_COL,
                        config.METAGENE_THREEPRIME_COL, config.A_COL,
                        config.T_COL, config.G_COL, config.C_COL,
                        config.UGA_COL, config.UAG_COL, config.UAA_COL,
                        config.UGA_COL, config.UAG_COL, config.UAA_COL,
                        config.CDS_MARKER_SIZE, config.CDS_MARKER_COLOUR,
                        config.LEGEND_SIZE, config.RIBO_LINEWIDTH))
            connection.commit()
            logging.debug("Closing trips.sqlite connection")
            connection.close()
            return redirect("/")
        else:
            error = 'Invalid Captcha. Please try again.'
            return render_template('create.html', error=error)
    else:
        return render_template('create.html', error=error)


# Allows users to change some global settings such as plot background colour, title size, tick label size, etc.
@app.route('/settings/')
# @login_required
def settingspage():
    global local
    try:
        print(local)
    except Exception:
        local = False

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

    return render_template('settings.html',
                           local=local,
                           user_settings=user_settings)


# Allows users to download fasta files, as well as scripts needed to produce their own sqlite files
@app.route('/downloads/')
def downloadspage():
    global local
    try:
        print(local)
    except Exception:
        local = False
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
                           local=local,
                           user=user,
                           organism_dict=organism_dict)


# Called when user downloads something from the downloads page
@app.route('/downloadquery', methods=['GET', 'POST'])
def download_file():
    if request.method == 'POST':
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
def uploadspage():
    global local
    try:
        print(local)
    except Exception:
        local = False
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
    organisms = organisms.loc[
        not organisms.private | organisms.owner == user_id,
        ["organism_name", "transcriptome_list", "organism_id"]]

    organism_dict = table_to_dict(organisms)
    organism_dict = {}

    cursor.execute(
        "SELECT organism_name,transcriptome_list,organism_id from organisms where private = 0 OR owner = {};"
        .format(user_id))
    result = cursor.fetchall()
    transcriptome = None  # TODO: Remove this part after correcting code
    for row in result:
        organism = row[0]
        transcriptome = row[1]
        if organism not in organism_dict:
            organism_dict[organism] = [transcriptome]
        else:
            organism_dict[organism].append(transcriptome)
        org_id_dict[row[2]] = [organism, transcriptome]
    cursor.execute(
        "SELECT organism_id from organism_access WHERE user_id = '{}';".format(
            user_id))
    result = (cursor.fetchall())
    for row in result:
        organism_id = row[0]
        cursor.execute(
            "SELECT organism_name,transcriptome_list from organisms where organism_id = {};"
            .format(organism_id))
        org_result = cursor.fetchone()
        if org_result:
            organism_name = org_result[0]
            transcriptome = org_result[1]
            if organism_name not in organism_dict:
                organism_dict[organism_name] = [transcriptome]
            else:
                organism_dict[organism_name].append(transcriptome)
            org_id_dict[organism_id] = [organism_name, transcriptome]
    study_dict = {}
    cursor.execute(
        "SELECT study_id,study_name,organism_id from studies where owner = {}".
        format(user_id))
    result = cursor.fetchall()
    for row in result:
        # add to organism dict if not caught earlier (this only happens when database is modified manually)
        if row[2] not in org_id_dict:
            cursor.execute(
                "SELECT organism_name,transcriptome_list,organism_id from organisms where organism_id =  {};"
                .format(row[2]))
            # organism_dict[organism] = [transcriptome]
            org_id_dict[row[2]] = [row[0], transcriptome]
        study_dict[int(row[0])] = [
            row[1].replace("_{}".format(user_id), "", 1),
            org_id_dict[row[2]][0], org_id_dict[row[2]][1], []
        ]

    transcriptome_dict = {}
    cursor.execute(
        "SELECT organism_id,organism_name,transcriptome_list from organisms where owner = {}"
        .format(user_id))
    result = cursor.fetchall()
    for row in result:
        transcriptome_dict[int(row[0])] = [row[1], row[2]]
    for study_id in study_dict:
        cursor.execute(
            "SELECT user_id FROM study_access WHERE study_id = {}".format(
                study_id))
        result = cursor.fetchall()
        for row in result:
            cursor.execute(
                "SELECT username FROM users WHERE user_id = {}".format(row[0]))
            result2 = cursor.fetchone()
            shared_user = result2[0]
            if shared_user not in study_dict[study_id][3]:
                study_dict[study_id][3].append(result2[0])
    file_dict = {}
    cursor.execute(
        "SELECT file_name,study_id,file_id,file_description from files where owner = {}"
        .format(user_id))
    result = cursor.fetchall()
    for row in result:
        cursor.execute(
            "SELECT study_name from studies where study_id = {}".format(
                row[1]))
        study_name = (cursor.fetchone())
        if study_name != None:
            file_dict[row[0]] = [
                study_name[0].replace("_{}".format(user_id), "", 1), row[2],
                row[3]
            ]
    seq_dict = {}
    cursor.execute(
        "SELECT seq_name,frame_breakdown from seq_rules where user_id = {}".
        format(user_id))
    result = cursor.fetchall()
    for row in result:
        seq_dict[row[0]] = [row[1]]
    connection.close()
    return render_template('uploads.html',
                           local=local,
                           user=user,
                           organism_dict=organism_dict,
                           study_dict=study_dict,
                           transcriptome_dict=transcriptome_dict,
                           file_dict=file_dict,
                           seq_dict=seq_dict)


# Called when user uploads something on the uploads page
@app.route('/uploadquery', methods=['GET', 'POST'])
# @login_required
def upload_file():
    if request.method == 'POST':
        # uploaded_files = request.files.getlist("file")
        f = request.files["file"]
        print(request.form)
        print(request.files)
        connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                    config.DATABASE_NAME))
        connection.text_factory = str
        cursor = connection.cursor()
        user, logged_in = fetch_user()

        cursor.execute(
            "SELECT user_id from users WHERE username = '{}';".format(user))
        result = (cursor.fetchone())
        user_id = result[0]
        filename = f.filename
        # filename = secure_filename(f.filename)
        print(filename)
        file_ext = filename.split(".")[-1]
        if file_ext != "sqlite":
            flash("Error: File extension should be sqlite not {}".format(
                file_ext))
            return redirect("https://trips.ucc.ie/uploads")
        if user != "public":
            foldername = "{}_{}".format(
                request.form["foldername"].replace(" ", "_"), user_id)
        else:
            foldername = "{}".format(request.form["foldername"].replace(
                " ", "_"))
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
        if filetype != "riboseq" and filetype != "rnaseq":
            cursor.execute(
                "SELECT * from seq_rules where user_id = {} and seq_name = '{}'"
                .format(user_id, filetype))
            result = cursor.fetchone()
            if result == None:
                cursor.execute(
                    "INSERT INTO seq_rules VALUES ({},'{}',0)".format(
                        user_id, filetype))

        if not os.path.isdir("{}/uploads/{}".format(config.SCRIPT_LOC,
                                                    foldername)):
            os.makedirs("{}/uploads/{}".format(config.SCRIPT_LOC, foldername))
        upload_file_path = "{}/uploads/{}/{}".format(config.SCRIPT_LOC,
                                                     foldername, filename)
        with open(upload_file_path, "wb") as fout:
            while True:
                chunk = f.stream.read(1024)
                if not chunk:
                    break
                fout.write(chunk)

        # f.save("{}/uploads/{}/{}".format(config.SCRIPT_LOC, foldername,
        # filename))
        sqlite_db = SqliteDict("{}/uploads/{}/{}".format(
            config.SCRIPT_LOC, foldername, filename))
        try:
            file_description = sqlite_db["description"]
        except Exception:
            file_description = "NULL"
        sqlite_db.close()
        # get file id
        cursor.execute("SELECT MAX(file_id) FROM files;")
        result = cursor.fetchone()
        new_file_id = int(result[0]) + 1
        cursor.execute(
            "SELECT organism_id FROM organisms WHERE organism_name = '{}' AND transcriptome_list = '{}';"
            .format(organism, assembly))
        result = cursor.fetchone()
        organism_id = int(result[0])
        # get study_id
        cursor.execute(
            "SELECT study_id FROM studies WHERE study_name = '{}' and organism_id = {} and owner = {}"
            .format(foldername, organism_id, user_id))
        result = cursor.fetchone()
        if result != None:
            study_id = int(result[0])
        else:
            cursor.execute("SELECT MAX(study_id) FROM studies;")
            result = cursor.fetchone()
            study_id = int(result[0]) + 1
            if user != "public":
                cursor.execute(
                    "INSERT INTO studies VALUES({},{},'{}','{}','{}','{}','{}','{}','{}','{}','{}','{}',{},{})"
                    .format(study_id, organism_id, foldername, 'NULL', 'NULL',
                            'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL',
                            'NULL', 1, user_id))
            else:
                cursor.execute(
                    "INSERT INTO studies VALUES({},{},'{}','{}','{}','{}','{}','{}','{}','{}','{}','{}',{},{})"
                    .format(study_id, organism_id, foldername, 'NULL', 'NULL',
                            'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL',
                            'NULL', 0, user_id))
            cursor.execute("INSERT INTO study_access VALUES({},{});".format(
                study_id, user_id))
        cursor.execute(
            "INSERT INTO files VALUES({},{},{},'{}','{}','{}',{},{},{},'{}')".
            format(new_file_id, organism_id, study_id, filename,
                   file_description, filetype, user_id, 0, 0, ""))
        # If user is not logged in only keep file for a set period of time
        if logged_in == False:
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
def upload_transcriptome():
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    connection.text_factory = str
    cursor = connection.cursor()
    user, logged_in = fetch_user()
    cursor.execute(
        "SELECT user_id from users WHERE username = '{}';".format(user))
    result = (cursor.fetchone())
    user_id = result[0]
    if request.method == 'POST':
        organism = (request.form["organism"]).lower().strip().replace(" ", "_")
        assembly = (request.form["assembly"]).lower().strip().replace(" ", "_")
        default_tran = (request.form["default_tran"]).lower().strip().replace(
            " ", "_")
        uploaded_annotation = request.files.getlist("anno_file")
        if not os.path.isdir("{}/uploads/transcriptomes/{}".format(
                config.SCRIPT_LOC, user_id)):
            os.makedirs("{}/uploads/transcriptomes/{}".format(
                config.SCRIPT_LOC, user_id))
        if not os.path.isdir("{}/uploads/transcriptomes/{}/{}".format(
                config.SCRIPT_LOC, user_id, organism)):
            os.makedirs("{}/uploads/transcriptomes/{}/{}".format(
                config.SCRIPT_LOC, user_id, organism))
        if not os.path.isdir("{}/uploads/transcriptomes/{}/{}/{}".format(
                config.SCRIPT_LOC, user_id, organism, assembly)):
            os.makedirs("{}/uploads/transcriptomes/{}/{}/{}".format(
                config.SCRIPT_LOC, user_id, organism, assembly))
        for f in uploaded_annotation:
            filename = secure_filename(f.filename)
            ext = filename.split(".")[-1]
            if ext != "sqlite":
                return """Error: Expecting extension sqlite but got extension {}. The file generated by the create_annotation_sqlite.py script should be uploaded here.
						This script can be gotten on the downloads page, by selecting the Scripts group.""".format(
                    ext)
            # Instead of using filename of the uploaded file we rename it to organism_assembly.sqlite, to keep things consistent
            filename = "{}_{}.sqlite".format(organism, assembly)
            full_path = "{}/uploads/transcriptomes/{}/{}/{}/{}".format(
                config.SCRIPT_LOC, user_id, organism, assembly, filename)
            f.save(full_path)
            cursor.execute("SELECT MAX(organism_id) from organisms;")
            max_org_id = (cursor.fetchone()[0]) + 1
            cursor.execute(
                "INSERT INTO organisms VALUES({},'{}','{}','NULL','NULL','NULL','NULL','{}',1,{})"
                .format(max_org_id, organism, assembly, default_tran.upper(),
                        user_id))
            cursor.execute("INSERT INTO organism_access VALUES ({},{})".format(
                max_org_id, user_id))
            if logged_in == False:
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
        return redirect("https://trips.ucc.ie/uploads")


# Called by flask in case of an error in the code, returns the exception so it can be displayed to user
@app.errorhandler(500)
def handle_bad_request(e):
    return_str = 'ERROR: ' + str(
        e
    ) + " please report this to tripsvizsite@gmail.com or via the contact page. "
    return return_str


# This is the page where users login.
@app.route("/user/login", methods=["GET", "POST"])
def login():
    global local
    # if user is already logged in then redirect to homepage
    if current_user.is_authenticated:
        return redirect("/")
    error = None
    if request.method == 'POST':
        username = str(request.form['username']).strip()
        password = str(request.form['password']).strip()
        if xcaptcha.verify() or local == True or username == "developer":
            username_dict = {}
            logging.debug("login Connecting to trips.sqlite")
            connection = sqlite3.connect('{}/trips.sqlite'.format(
                config.SCRIPT_LOC))
            connection.text_factory = str
            cursor = connection.cursor()
            cursor.execute("SELECT username,password from users;")
            result = (cursor.fetchall())
            logging.debug("Closing trips.sqlite connection")
            connection.close()
            for row in result:
                username_dict[row[0]] = row[1]
            if username in username_dict:
                if check_password_hash(username_dict[username],
                                       password) == True or local == True:
                    id = username
                    user = User(id)
                    login_user(user)
                    nxt = sanitize_get_request(request.args.get('next'))
                    if nxt != None:
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
def logout():
    logout_user()
    return redirect(url_for('homepage'))


# callback to reload the user object
@login_manager.user_loader
def load_user(userid):
    return User(userid)


# Forgot password page.
# @app.route("/forgot", methods=["GET", "POST"])
# def forgot():
# if user is already logged in then redirect to homepage
# if current_user.is_authenticated:
# return redirect("/")
# error=None
# if request.method == 'POST':
# emailaddress = str(request.form['emailaddress'])
# print ("address is ", emailaddress)
# return redirect("/")
# else:
# return render_template('forgot_password.html',error=error)


@app.route("/reset", methods=["GET", "POST"])
def reset():
    return None


# Called when user presses the save button on the orf_translation page.
@app.route('/anno_query', methods=['POST'])
def anno_query():
    data = json.loads(request.data)
    user = fetch_user()[0]
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    cursor = connection.cursor()
    cursor.execute(
        "SELECT user_id from users WHERE username = '{}';".format(user))
    result = (cursor.fetchone())
    user_id = result[0]
    cursor.execute(
        "INSERT INTO users_saved_cases VALUES('{}','{}',{},{},{},{},'{}','START_CODON','CDS_OVERLAP','START_SCORE','STOP_SCORE','ENTROPY','TE','COVERAGE','CDS_RATIO','{}','FILE_LIST',{},'{}','{}','{}');"
        .format(data["gene"].strip(), data["transcript"].strip(),
                data["start"].strip(), data["stop"].strip(),
                data["length"].strip(), data["score"].strip(),
                data["type"].strip(), data["trips_link"].strip(), user_id,
                data["label"].strip(), data["organism"].strip(),
                data["transcriptome"].strip()))
    connection.commit()
    connection.close()
    return ""


# This page shows the saved ORFs specific to the signed in user
@app.route('/saved/')
def saved():
    global local
    try:
        print(local)
    except Exception:
        local = False
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    connection.text_factory = str
    cursor = connection.cursor()
    advanced = False
    user, logged_in = fetch_user()
    # If user is not logged in and has rejected cookies they cannot use this page, so redirect to the homepage.
    if user == None:
        flash(
            "To use the Saved ORFs page you either need to be logged in or allow cookies."
        )
        return redirect(url_for('homepage'))

    user_id = -1
    if user != None:
        if logged_in == True:
            flash("You are logged in as {}".format(user))
        cursor.execute(
            "SELECT user_id from users WHERE username = '{}';".format(user))
        result = (cursor.fetchone())
        user_id = result[0]
    cursor.execute(
        "SELECT organism_id from organism_access WHERE user_id = '{}';".format(
            user_id))
    result = (cursor.fetchall())
    organism_access_list = []
    for row in result:
        organism_access_list.append(row[0])
    cursor.execute("SELECT organism_id,organism_name,private from organisms;")
    organism_list = []
    # List of all rows returned
    result = (cursor.fetchall())
    for row in result:
        if row[2] == 0:
            organism_list.append(str(row[1]))
        elif row[2] == 1:
            if row[0] in organism_access_list:
                organism_list.append(str(row[1]))
    connection.close()

    return render_template('user_saved_cases.html',
                           local=local,
                           advanced=advanced,
                           organism_list=organism_list)


# Retrieves saved ORFs
@app.route('/savedquery', methods=['POST'])
def savedquery():
    data = json.loads(request.data)
    user = fetch_user()[0]
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    cursor = connection.cursor()
    organism = data["organism"]
    label = data["label"]
    returnstr = ""
    # get user_id
    cursor.execute(
        "SELECT user_id from users WHERE username = '{}';".format(user))
    result = (cursor.fetchone())
    user_id = result[0]
    start_codons = []
    if "sc_aug" in data:
        start_codons.append("ATG")
    if "sc_cug" in data:
        start_codons.append("CTG")
    if "sc_gug" in data:
        start_codons.append("GTG")
    if "sc_none" in data:
        start_codons.append("None")

    # structure of orf dict is transcript[stop][start] = {"length":x,"score":0,"cds_cov":0} each stop can have multiple starts
    if organism != 'Select an Organism':
        if label == "":
            cursor.execute(
                "SELECT * FROM users_saved_cases WHERE user_id = {} and organism = '{}';"
                .format(user_id, organism))
        else:
            label_list = (label.strip(" ")).split(",")
            str_list = str(label_list).strip("[]")
            cursor.execute(
                "SELECT * FROM users_saved_cases WHERE user_id = {} AND label IN ({}) and organism = '{}';"
                .format(user_id, str_list, organism))
    else:
        if label == "":
            cursor.execute(
                "SELECT * FROM users_saved_cases WHERE user_id = {};".format(
                    user_id))
        else:
            label_list = (label.strip(" ")).split(",")
            str_list = str(label_list).strip("[]")
            cursor.execute(
                "SELECT * FROM users_saved_cases WHERE user_id = {} AND label IN ({});"
                .format(user_id, str_list))
    result = cursor.fetchall()
    total_rows = 0
    for row in result:
        gene = str(row[0])
        tran = str(row[1])
        start = str(row[2])
        stop = str(row[3])
        length = str(row[4])
        score = str(row[5])
        label = str(row[18])
        start_codon = str(row[7])

        trips_link = str(row[15])

        # Limit the number of returned cases to 1000, as datatable is memory intensive
        if total_rows < 1000:
            returnstr += "{},{},{},{},{},{},{},{},{}.,/".format(
                gene, tran, start, stop, length, score, label, start_codon,
                trips_link)
            total_rows += 1
        else:
            break
    connection.close()
    return returnstr


# Allows users to delete previously saved cases
@app.route('/del_query', methods=['POST'])
def del_query():
    data = json.loads(request.data)
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    cursor = connection.cursor()
    try:
        user = current_user.name
    except Exception:
        user = None
        return "Error user not signed in"
    cursor.execute(
        "SELECT user_id from users WHERE username = '{}';".format(user))
    result = (cursor.fetchone())
    user_id = result[0]
    cursor.execute(
        "DELETE from users_saved_cases WHERE tran = '{}' and start = '{}' and stop = '{}' and trips_link = '{}' and label = '{}' and user_id = {};"
        .format(data["transcript"], data["start"], data["stop"],
                data["trips_link"], data["label"], user_id))
    connection.commit()
    connection.close()
    return ""


# Allows users to logout
# @app.route("/user/logout")
# @login_required
# def logout():
#	logout_user()
#	return redirect(login)


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
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    connection.text_factory = str
    cursor = connection.cursor()
    cursor.execute("SELECT url from urls WHERE url_id = '{}';".format(integer))
    result = cursor.fetchone()
    if result == None:
        return "Short code not recognized."
    url = ""
    url += result[0]
    # add a keyword to the url to prevent generating another shortcode
    url += "&short={}".format(short_code)
    connection.close()
    return redirect(url)


@app.after_request
def after_request_func(response):
    consent = request.cookies.get("cookieconsent_status")
    # if consent == "deny":
    #	session.clear()
    #	for cookie_name in request.cookies:
    #		if cookie_name != "cookieconsent_status":
    #			response.delete_cookie(cookie_name)
    #			if cookie_name == "tripsviz_session" or cookie_name == "session":
    #				response.delete_cookie(cookie_name,path='/',domain='trips.ucc.ie')
    #			else:
    #				response.delete_cookie(cookie_name,path='/',domain='.ucc.ie')
    return response


# This is the home page it show a list of organisms as defined by trips_dict
# TODO: Move pages in order
@app.route('/')
def homepage(message=""):
    organism_access_list = []
    organism_list = {}
    logging.debug("homepage Connecting to trips.sqlite")
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    connection.text_factory = str
    cursor = connection.cursor()
    sanitize_get_request(request.cookies.get("cookieconsent_status"))
    user, logged_in = fetch_user()

    user_id = -1
    if user != None:
        if logged_in == True:
            flash("You are logged in as {}".format(user))
        cursor.execute(
            "SELECT user_id from users WHERE username = '{}';".format(user))
        result = (cursor.fetchone())
        user_id = result[0]
        # get a list of organism id's this user can access
        cursor.execute(
            "SELECT organism_id from organism_access WHERE user_id = '{}';".
            format(user_id))
        result = (cursor.fetchall())
        for row in result:
            organism_access_list.append(int(row[0]))

    # returns a tuple with each field as a seperate string
    cursor.execute(
        "SELECT organism_id,organism_name,private,owner from organisms;"
    )  #TODO: Seq uniq
    # List of all rows returned
    result = (cursor.fetchall())
    for row in result:  # TODO: Need to make it more efficient
        if row[2] == 0:
            if row[1] not in organism_list:
                organism_list[row[1]] = []
        elif row[2] == 1:
            if row[0] in organism_access_list or row[3] == user_id:
                if row[1] not in organism_list:
                    organism_list[row[1]] = []
    # organism_list.sort()

    for organism in organism_list:
        cursor.execute(
            "SELECT transcriptome_list from organisms WHERE organism_name = '{}' AND (private = 0 OR organism_id IN ({})) ;"
            .format(organism,
                    str(organism_access_list).strip("[]")))
        result = cursor.fetchall()
        if result:
            for row in result:
                organism_list[organism].append(row[0])

    logging.debug("homepage Closing trips.sqlite connection")
    connection.close()
    print(organism_list)

    message = sanitize_get_request(request.args.get('message'))
    return render_template('landing.html',
                           organisms=organism_list,
                           message=message)


# show a list of plot types
@app.route('/<organism>/<transcriptome>/')
def plogpage(organism, transcriptome):
    try:  # TODO: Integrate login details in main page and remove this function
        user = current_user.name
    except Exception:
        user = None
    return render_template('plot_types.html', current_username=user)


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
    if logged_in == True:
        new_password = data['new_password']
        new_password2 = data['new_password2']
        curr_password = data['curr_password']
        if new_password != "" or new_password2 != "":
            if new_password != new_password2:
                return "ERROR: New passwords do not match"
            generate_password_hash(curr_password)
            cursor.execute(
                "SELECT password FROM users WHERE username = '{}'".format(
                    user))
            result = cursor.fetchone()
            old_password_hash = result[0]
            if check_password_hash(old_password_hash, curr_password) == True:
                new_password_hash = generate_password_hash(new_password)
                cursor.execute(
                    "UPDATE users SET password = '{}' WHERE username = '{}'".
                    format(new_password_hash, user))
                connection.commit()
            else:
                return "ERROR: Current password is not correct"
    # get user_id
    cursor.execute(
        "SELECT user_id from users WHERE username = '{}';".format(user))
    result = (cursor.fetchone())
    user_id = result[0]

    if "advanced" in data:
        cursor.execute(
            "UPDATE users SET advanced = 1 WHERE user_id = '{}';".format(
                user_id))
        connection.commit()
    else:
        cursor.execute(
            "UPDATE users SET advanced = 0 WHERE user_id = '{}';".format(
                user_id))
        connection.commit()
    # TODO: Check the input before updating the settings

    # get a list of organism id's this user can access
    cursor.execute(
        "UPDATE user_settings SET background_col = '{}' WHERE user_id = '{}';".
        format(data["background_colour"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET title_size = '{}' WHERE user_id = '{}';".
        format(data["title_size"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET subheading_size = '{}' WHERE user_id = '{}';"
        .format(data["subheading_size"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET axis_label_size = '{}' WHERE user_id = '{}';"
        .format(data["axis_label_size"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET marker_size = '{}' WHERE user_id = '{}';".
        format(data["marker_size"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET legend_size = '{}' WHERE user_id = '{}';".
        format(data["legend_size"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET ribo_linewidth = '{}' WHERE user_id = '{}';".
        format(data["ribo_linewidth"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET cds_marker_width = '{}' WHERE user_id = '{}';"
        .format(data["cds_marker_width"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET cds_marker_colour = '{}' WHERE user_id = '{}';"
        .format(data["cds_marker_colour"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET readlength_col = '{}' WHERE user_id = '{}';".
        format(data["readlength_colour"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET metagene_fiveprime_col = '{}' WHERE user_id = '{}';"
        .format(data["metagene_fiveprime_colour"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET metagene_threeprime_col = '{}' WHERE user_id = '{}';"
        .format(data["metagene_threeprime_colour"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET nuc_comp_a_col = '{}' WHERE user_id = '{}';".
        format(data["nuc_comp_a_col"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET nuc_comp_t_col = '{}' WHERE user_id = '{}';".
        format(data["nuc_comp_t_col"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET nuc_comp_g_col = '{}' WHERE user_id = '{}';".
        format(data["nuc_comp_g_col"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET nuc_comp_c_col = '{}' WHERE user_id = '{}';".
        format(data["nuc_comp_c_col"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET uag_col = '{}' WHERE user_id = '{}';".format(
            data["uag_col"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET uga_col = '{}' WHERE user_id = '{}';".format(
            data["uga_col"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET uaa_col = '{}' WHERE user_id = '{}';".format(
            data["uaa_col"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET comp_uag_col = '{}' WHERE user_id = '{}';".
        format(data["comp_uag_col"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET comp_uga_col = '{}' WHERE user_id = '{}';".
        format(data["comp_uga_col"].strip(), user_id))
    connection.commit()
    cursor.execute(
        "UPDATE user_settings SET comp_uaa_col = '{}' WHERE user_id = '{}';".
        format(data["comp_uaa_col"].strip(), user_id))
    connection.commit()
    connection.close()
    flash("Settings have been updated")
    return "Settings have been updated"


# Allows users to delete files
@app.route('/deletequery', methods=['GET', 'POST'])
# @login_required
def deletequery():
    data = json.loads(request.data)

    user = fetch_user()[0]
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    connection.text_factory = str
    cursor = connection.cursor()
    cursor.execute(
        "Select user_id from users where username = '{}';".format(user))
    result = cursor.fetchone()
    user_id = result[0]
    for key in data:
        file_id = data[key]["file_id"]
        if "filecheck" in data[key]:
            cursor.execute(
                "SELECT owner, study_id, file_name FROM files WHERE file_id = {}"
                .format(file_id))
            result = cursor.fetchone()
            owner = result[0]
            study_id = result[1]
            filename = result[2]
            cursor.execute(
                "SELECT * FROM studies WHERE study_id = {}".format(study_id))
            result = cursor.fetchone()
            study_name = result[2]
            full_path = "{}{}/{}".format(config.UPLOADS_DIR, study_name,
                                         filename)
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
            cursor.execute(
                "SELECT organism_id FROM files WHERE file_id = {}".format(
                    file_id))
            result = cursor.fetchone()
            cursor.execute(
                "SELECT organism_name FROM organisms WHERE organism_id = {}".
                format(result[0]))
            organism = cursor.fetchone()[0]
            filepath_dict = fetch_file_paths([file_id], organism)
            for seq_type in filepath_dict:
                if file_id in filepath_dict[seq_type]:
                    filepath = filepath_dict[seq_type][file_id]
                    opendict = SqliteDict(filepath, autocommit=True)
                    opendict["cutadapt_removed"] = int(
                        data[key]["cutadapt_removed"])
                    opendict.close()

        if data[key]["rrna_removed"] != '0':
            cursor.execute(
                "SELECT organism_id FROM files WHERE file_id = {}".format(
                    file_id))
            result = cursor.fetchone()
            cursor.execute(
                "SELECT organism_name FROM organisms WHERE organism_id = {}".
                format(result[0]))
            organism = cursor.fetchone()[0]
            filepath_dict = fetch_file_paths([file_id], organism)
            for seq_type in filepath_dict:
                if file_id in filepath_dict[seq_type]:
                    filepath = filepath_dict[seq_type][file_id]
                    opendict = SqliteDict(filepath, autocommit=True)
                    opendict["rrna_removed"] = int(data[key]["rrna_removed"])
                    opendict.close()
        if data[key]["unmapped"] != '0':
            cursor.execute(
                "SELECT organism_id FROM files WHERE file_id = {}".format(
                    file_id))
            result = cursor.fetchone()
            cursor.execute(
                "SELECT organism_name FROM organisms WHERE organism_id = {}".
                format(result[0]))
            organism = cursor.fetchone()[0]
            filepath_dict = fetch_file_paths([file_id], organism)
            for seq_type in filepath_dict:
                if file_id in filepath_dict[seq_type]:
                    filepath = filepath_dict[seq_type][file_id]
                    opendict = SqliteDict(filepath, autocommit=True)
                    opendict["unmapped_reads"] = int(data[key]["unmapped"])
                    opendict.close()
    connection.commit()
    connection.close()
    flash("File list updated")
    return redirect("https://trips.ucc.ie/uploads")


# Allows users to delete studies,modify access, modify the organism/transcriptome assembly or study name
@app.route('/deletestudyquery', methods=['GET', 'POST'])
# @login_required
def deletestudyquery():
    data = json.loads(request.data)
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    connection.text_factory = str
    cursor = connection.cursor()

    user = fetch_user()[0]

    cursor.execute(
        "SELECT user_id from users WHERE username = '{}';".format(user))
    result = (cursor.fetchone())
    user_id = result[0]
    for study_id in data:
        studycheck = data[study_id][0]
        # Delete studies where the "delete" checkbox is checked

        if studycheck.split("_")[-1] != "undefined":
            study_id = studycheck.split("_")[-1]
            print("deleting study_id", study_id)

            # First delete all files on the server associated with this study, if there are any
            cursor.execute(
                "SELECT * FROM files WHERE study_id = {}".format(study_id))
            result = cursor.fetchall()
            if result != None:
                for row in result:
                    file_id = row[0]
                    filename = row[3]
                    cursor.execute(
                        "SELECT owner, study_name FROM studies WHERE study_id = {}"
                        .format(study_id))
                    result = cursor.fetchone()
                    owner = row[0]
                    study_name = row[1]
                    full_path = "{}{}/{}".format(config.UPLOADS_DIR,
                                                 study_name, filename)
                    if owner == user_id:
                        if os.path.isfile(full_path):
                            # Instead of deleting the file now, add it to deletions table where it will be deleted via cron job, this will give users time to contact in case of accidental deletion
                            curr_time = time.time()
                            # The time to keep the file in seconds, currently set to 14 days
                            keep_time = 60 * 60 * 24 * 14
                            deletion_time = curr_time + keep_time
                            cursor.execute(
                                "INSERT INTO deletions VALUES({},'{}',{})".
                                format(file_id, full_path, deletion_time))
            # Now remove the study and the files associated with it from the db
            cursor.execute(
                "SELECT owner FROM studies WHERE study_id = {}".format(
                    study_id))
            result = cursor.fetchone()
            study_owner = result[0]
            if study_owner == user_id:
                # os.rename(study_path,study_path+"_REMOVE")
                cursor.execute(
                    "DELETE FROM studies WHERE study_id = {}".format(study_id))
                cursor.execute(
                    "DELETE FROM files WHERE study_id = {}".format(study_id))
                connection.commit()
                continue

        # Modify access list next to studies
        study_access = data[study_id][1].split(",")
        # check study_access against a list of all users
        all_users = {}
        cursor.execute("SELECT username,user_id FROM users;")
        result = cursor.fetchall()
        for row in result:
            all_users[row[0]] = row[1]
        # Check that all users exist
        for username in study_access:
            if username:
                if username not in all_users.keys():
                    flash(
                        "Error: User {} is not registered on Trips-Viz".format(
                            username))
                    return str(get_flashed_messages())
                else:
                    cursor.execute(
                        "SELECT * FROM study_access WHERE user_id = {} and study_id = {};"
                        .format(all_users[username], study_id))
                    result = cursor.fetchone()
                    if result == None:
                        cursor.execute(
                            "\n\n\nINSERT INTO study_access VALUES({},{});".
                            format(study_id, all_users[username]))

        # Modify study names if they have changed
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
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    connection.text_factory = str
    cursor = connection.cursor()

    user = fetch_user()[0]

    cursor.execute(
        "SELECT user_id from users WHERE username = '{}';".format(user))
    result = (cursor.fetchone())
    user_id = result[0]
    for organism_id in data:
        transcriptomecheck = data[organism_id][0]
        if transcriptomecheck.split("_")[-1] != "undefined":
            organism_id = transcriptomecheck.split("_")[-1]
            # Delete the annotation sqlite file
            cursor.execute(
                "SELECT organism_name, transcriptome_list, owner FROM organisms WHERE organism_id = {}"
                .format(organism_id))
            result = cursor.fetchone()
            organism_name = result[0]
            transcriptome_list = result[1]
            owner = result[2]
            if owner != user_id:
                continue
            sqlite_path = "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
                config.UPLOADS_DIR, user_id, organism_name, transcriptome_list)
            if os.path.isfile(sqlite_path):
                # os.rename(sqlite_path,sqlite_path+"_REMOVE")
                # Instead of deleting the file now, add it to deletions table where it will be deleted via cron job, this will give users time to contact in case of accidental deletion
                curr_time = time.time()
                # The time to keep the file in seconds, currently set to 14 days
                keep_time = 60 * 60 * 24 * 14
                deletion_time = curr_time + keep_time
                cursor.execute(
                    "INSERT INTO org_deletions VALUES({},'{}',{})".format(
                        organism_id, sqlite_path, deletion_time))
            # sqlite_dir = "{0}transcriptomes/{1}/{2}/{3}".format(config.UPLOADS_DIR, user_id, organism_name,transcriptome_list)
            # if os.path.isdir(sqlite_dir):
            #	os.rename(sqlite_dir,sqlite_dir+"_REMOVE")

            cursor.execute(
                "DELETE FROM organisms WHERE organism_id = {}".format(
                    organism_id))
            # delete all files on the server associated with this organism, if there are any
            cursor.execute("SELECT * FROM files WHERE organism_id = {}".format(
                organism_id))
            result = cursor.fetchall()
            study_ids = []
            study_names = []
            if result != None:
                for row in result:
                    filename = row[3]
                    study_id = row[2]
                    study_ids.append(study_id)
                    cursor.execute(
                        "SELECT * FROM studies WHERE study_id = {}".format(
                            study_id))
                    result = cursor.fetchone()
                    study_name = result[2]
                    study_names.append(study_name)
                    full_path = "{}{}/{}".format(config.UPLOADS_DIR,
                                                 study_name, filename)
                    if os.path.isfile(full_path):
                        os.remove(full_path)

            for study_name in study_names:
                full_path = "{}{}".format(config.UPLOADS_DIR, study_name)
                if os.path.isdir(full_path):
                    os.rename(full_path, full_path + "_REMOVE")

            # Now remove the study and the files associated with it from the db
            for study_id in study_ids:
                cursor.execute(
                    "DELETE FROM studies WHERE study_id = {}".format(study_id))
                cursor.execute(
                    "DELETE FROM files WHERE study_id = {}".format(study_id))
            connection.commit()
    connection.close()
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
    global local
    try:
        print(local)
    except Exception:
        local = False

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
    return riboflask_datasets.generate_plot(xlist, ylist, filenames,
                                            file_descs, studies, raw_reads,
                                            controls, cell_lines,
                                            control_colors, study_colors,
                                            cell_line_colors, transcript,
                                            start, stop)


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
