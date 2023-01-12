import os
import time
import sys
import MySQLdb
import sqlite3
from math import floor
import string
from flask import Flask,get_flashed_messages, render_template, json, request, make_response, send_from_directory, flash, Response, abort, redirect, url_for
from flask_recaptcha import ReCaptcha
#from flask_caching import Cache
from threading import Lock
lock = Lock()
import ast
import metainfo_plots
import stats_plots
import numpy as np
from math import log
import collections
import operator
from flask_login import LoginManager, UserMixin, login_required, login_user, logout_user, current_user
from werkzeug.security import generate_password_hash, check_password_hash
from werkzeug.exceptions import BadRequest
import config
import riboflask
import riboflask_compare
import riboflask_diff
from werkzeug import secure_filename
from sqlitedict import SqliteDict
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
import smtplib
from random import shuffle
from random import seed
from bisect import bisect_left 



#connection = MySQLdb.connect (host = "localhost", user = "stephenk", passwd = "ribosome", db = "trips")


trips_shelves_location = "/home/DATA/www/tripsviz/tripsviz/trips_shelves/"
trips_annotation_location = "/home/DATA/www/tripsviz/tripsviz/trips_annotations/"
trips_uploads_location = "/home/DATA/www/tripsviz/tripsviz/uploads/"



# Given a list of file id's as strings returns a list of paths to the shelf files.
def fetch_file_paths(file_list,organism):
	#print "RECEIVED FILE LIST, ORGANISM ", file_list, organism
	file_path_dict = {"riboseq":{},"rnaseq":{}}
	#Convert to a tuple so it works with mysql
	try:
		int_file_list = [int(x) for x in file_list]
	except:
		return {}
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	study_dict ={}
	cursor.execute("SELECT study_id,study_name from studies;")
	result = cursor.fetchall();
	for row in result:
		study_dict[row[0]] = row[1]
	if int_file_list:
		# Mysql expects a tuple but freaks out if you pass a tuple of 1 so we awkwardly replace [ and ] with ( and ) instead
		cursor.execute("SELECT file_type,file_name,study_id,file_id,owner from files WHERE file_id in {};".format(str(int_file_list).replace("[","(").replace("]",")")))
		result = cursor.fetchall()
		for row in result:
			study_name = study_dict[row[2]]
			if row[0] not in file_path_dict:
				file_path_dict[row[0]] = {}

			# All files uploaded on our side are owned by user id 1, all others are user uploaded
			# Used to store extension in sql db but don't anymore, hence we replace .shelf and .sqlite with nothing to deal with old cases, this can be removed eventually
			if row[4] == 1:
				file_path_dict[row[0]][row[3]] = ("{}{}/{}/{}/{}.sqlite".format(trips_shelves_location,row[0],organism,study_name,row[1].replace(".shelf","").replace(".sqlite","")))
			else:
				file_path_dict[row[0]][row[3]] = ("{}/{}/{}.sqlite".format(trips_uploads_location,study_name,row[1].replace(".shelf","").replace(".sqlite","")))
	connection.close()
	#print "returning file path dict", file_path_dict
	return file_path_dict







## Setting up matplotlib sytles using BMH
#s = json.load(open("./static/bmh_matplotlibrc.json"))

app = Flask(__name__, static_folder='static')#, static_url_path='')

#cache = Cache(app, config={'CACHE_TYPE': 'simple'})


app.config.from_pyfile('config.py')
recaptcha = ReCaptcha(app=app)
#recaptcha.init(app, '6LcXTDYUAAAAAMWgED__TFbCOc7MEt7w4z43eLpF', '6LcXTDYUAAAAALM80l8tHmvXqOplXSSB_6L9e5cH', is_enabled=True)

app.config['UPLOAD_FOLDER'] = '/static/tmp'
#change cookie name and path, this avoids cookies clashing with other flask apps on the same server (this doesn't work in local mode)
try:
	if sys.argv[1] == "true":
		pass
except:
	app.config.update(
			SESSION_COOKIE_NAME = 'session_tripsviz',
			SESSION_COOKIE_PATH = '/'
			)






# flask-login
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = "login"
login_manager.login_message = None


# silly user model
class User(UserMixin):

	def __init__(self, id):
		self.id = id
		self.name = str(id)
		self.password = self.name + "_secret"

	def is_authenticated(self):
		return self.is_authenticated

	def get_id(self):
		"""Return the email address to satisfy Flask-Login's requirements."""
		"""Requires use of Python 3"""
		return str(self.id)
	def __repr__(self):
		return "%d/%s/%s" % (self.id, self.name, self.password)





# Provides statistics on trips such as number of organisms, number of files, number of studies etc.
@app.route('/stats/')

def statisticspage():
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	cursor.execute("SELECT organism_id from organisms WHERE private = 0;")
	result = (cursor.fetchall())
	no_organisms = len(result)

	public_studies = []
	cursor.execute("SELECT study_id from studies WHERE private = 0;")
	result = (cursor.fetchall())
	for row in result:
		study_id = row[0]
		public_studies.append(study_id)
	no_studies = len(public_studies)

	riboseq_files = 0
	cursor.execute("SELECT study_id,file_id from files WHERE file_type = 'riboseq';")
	result = (cursor.fetchall())
	for row in result:
		if row[0] in public_studies:
			riboseq_files += 1

	rnaseq_files = 0
	cursor.execute("SELECT study_id,file_id from files WHERE file_type = 'rnaseq';")
	result = (cursor.fetchall())
	for row in result:
		if row[0] in public_studies:
			rnaseq_files += 1


	# Create the graph which breaksdown the number of studies per organism
	org_dict = {}
	cursor.execute("SELECT organism_id,organism_name from organisms WHERE private = 0;")
	result = (cursor.fetchall())
	for row in result:
		if "_" in row[1]:
			org_name = (row[1].split("_")[0][0])+"."+(row[1].split("_")[1])
		else:
			org_name = row[1]
		org_dict[row[0]] = {"organism_name":org_name,
							"riboseq_files":0,
							"rnaseq_files":0}
	cursor.execute("SELECT file_type,organism_id from files;")
	result = (cursor.fetchall())
	for row in result:
		if row[1] in org_dict:
			if row[0] == "riboseq":
				org_dict[row[1]]["riboseq_files"] += 1
			elif row[0] == "rnaseq":
				org_dict[row[1]]["rnaseq_files"] += 1
	read_dict = {"organisms":[],
					 "riboseq_files":[],
					 "rnaseq_files":[]}
	for org in org_dict:
		read_dict["organisms"].append(org_dict[org]["organism_name"])
		read_dict["riboseq_files"].append(org_dict[org]["riboseq_files"])
		read_dict["rnaseq_files"].append(org_dict[org]["rnaseq_files"])
	org_breakdown_graph = stats_plots.org_breakdown_plot(read_dict)



	# Create the graph which breaks down studies published per year
	year_dict = {}
	cursor.execute("SELECT paper_year from studies WHERE private = 0;")
	result = cursor.fetchall()
	for row in result:
		if row[0] == "None" or row[0] == "NULL":
			continue
		if int(row[0]) not in year_dict:
			year_dict[int(row[0])] = 0
		year_dict[int(row[0])] += 1
	final_year_dict = {"year":[],
					   "no_studies":[]}
	min_year = min(year_dict.keys())
	max_year = max(year_dict.keys())
	for i in range(min_year,max_year+1):
		final_year_dict["year"].append(str(i))
		if i in year_dict:
			final_year_dict["no_studies"].append(year_dict[i])
		else:
			final_year_dict["no_studies"].append(0)

	year_plot = stats_plots.year_dist(final_year_dict)
	
	news_string = ""

	cursor.execute("SELECT * from updates ORDER BY date DESC;")
	result = cursor.fetchall()
	
	for row in result:
		news_string += "<tr><td>{}</td><td>{}</td></tr>".format(row[0], row[1])
	
	
	
	connection.close()
	return render_template('statistics.html', no_organisms=no_organisms, no_studies=no_studies, riboseq_files=riboseq_files, rnaseq_files=rnaseq_files,org_breakdown_graph=org_breakdown_graph,year_plot=year_plot,news_string=news_string)




# Provides statistics on trips such as number of organisms, number of files, number of studies etc.
@app.route('/contactus/',methods=["GET", "POST"])

def contactus():
	if request.method == "POST":
		name = str(request.form['name'])
		email = str(request.form['email'])
		fromaddr = "ribopipe@gmail.com"
		toaddr = "tripsvizsite@gmail.com"
		msg = MIMEMultipart()
		msg['From'] = fromaddr
		msg['To'] = toaddr
		msg['Subject'] = str(request.form['subject'])
		msg.attach(MIMEText("Name: {}\nEmail: {}\nMessage: {}".format(str(request.form['name']),str(request.form['email']),str(request.form['message']))))
		server = smtplib.SMTP('smtp.gmail.com', 587)
		server.starttls()
		#TODO, move this to the config file
		server.login(fromaddr, "Ribosome")
		text = msg.as_string()
		server.sendmail(fromaddr, toaddr, text)
		server.quit()
		flash("Message sent successfully")
	return render_template('contact.html')







#Allows users to change some global settings such as plot background colour
@app.route('/settings/')

@login_required
def settingspage():

	global local
	try:
		print local
	except:
		local = False

	user = current_user.name
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	#get user_id
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	cursor.execute("SELECT background_col,readlength_col,metagene_fiveprime_col,metagene_threeprime_col,nuc_comp_a_col,nuc_comp_t_col,nuc_comp_g_col,nuc_comp_c_col,uag_col,uaa_col,uga_col,comp_uag_col,comp_uaa_col,comp_uga_col,title_size,subheading_size,axis_label_size,marker_size,cds_marker_width,cds_marker_colour,legend_size from user_settings WHERE user_id = {};".format(user_id))
	result = (cursor.fetchone())
	background_colour = result[0]
	readlength_colour = result[1]
	metagene_fiveprime_colour = result[2]
	metagene_threeprime_colour = result[3]
	nuc_comp_a_col = result[4]
	nuc_comp_t_col = result[5]
	nuc_comp_g_col = result[6]
	nuc_comp_c_col = result[7]
	uag_col = result[8]
	uaa_col = result[9]
	uga_col = result[10]
	comp_uag_col = result[11]
	comp_uaa_col = result[12]
	comp_uga_col = result[13]
	title_size = result[14]
	subheading_size = result[15]
	axis_label_size = result[16]
	marker_size = result[17]
	cds_marker_size = result[18]
	cds_marker_colour = result[19]
	legend_size = result[20]


	connection.close()
	return render_template('settings.html',
						   local=local,
						   background_colour=background_colour,
						   readlength_colour=readlength_colour,
						   metagene_fiveprime_colour=metagene_fiveprime_colour,
						   metagene_threeprime_colour=metagene_threeprime_colour,
						   nuc_comp_a_col=nuc_comp_a_col,
						   nuc_comp_t_col=nuc_comp_t_col,
						   nuc_comp_g_col=nuc_comp_g_col,
						   nuc_comp_c_col=nuc_comp_c_col,
						   uag_col=uag_col,
						   uaa_col=uaa_col,
						   uga_col=uga_col,
						   comp_uag_col=comp_uag_col,
						   comp_uaa_col=comp_uaa_col,
						   comp_uga_col=comp_uga_col,
						   title_size=title_size,
						   subheading_size=subheading_size,
						   axis_label_size=axis_label_size,
						   marker_size=marker_size,
						   cds_marker_width=cds_marker_size,
						   cds_marker_colour=cds_marker_colour,
						   legend_size=legend_size)





@app.route('/downloads/')
def downloadspage():

	global local
	try:
		print local
	except:
		local = False
	organism_dict = {"Scripts":["bam_to_sqlite.py","tsv_to_sqlite.py","create_annotation_sqlite.py","create_transcriptomic_to_genomic_sqlite.py"]}
	#organism_name            | transcriptome_list
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	try:
		user = current_user.name
	except:
		user = None

	cursor.execute("SELECT organism_name from organisms where private = 0")
	result = cursor.fetchall()
	for row in result:
		organism = row[0]
		organism_dict[organism] = []
	#organism_dict = {"Homo_sapiens":["GEN_V24","RFSQ"],"MOUSE":["gen_vmouse","notRFSQ"]}


	#TODO: Change this to make it relative to the tripsviz root
	trips_annotation_dir = "/home/DATA/www/tripsviz/tripsviz/trips_annotations/"
	for org in os.listdir(trips_annotation_dir):
		if org not in organism_dict:
			continue
		for filename in os.listdir(trips_annotation_dir+"/"+org):
			if "." in filename:
				ext = filename.split(".")[-1]
				if ext == "fa" or ext == "gtf":
					organism_dict[org].append(filename)
				elif ext == "sqlite":
					if "transcriptomic" in filename or org in filename:
						organism_dict[org].append(filename)

	connection.close()
	return render_template('downloads.html',
						   local=local,
						   user=user,
						   organism_dict=organism_dict
						   )


@app.route('/downloadquery', methods = ['GET', 'POST'])
def download_file():
	print "\n\n\n\n\n\n\nUPLOAD FILE CALLED"
	if request.method == 'POST':
		organism = request.form["organism"]
		assembly = request.form["assembly"]
		return send_from_directory("/home/DATA/www/tripsviz/tripsviz/trips_annotations/{}".format(organism), assembly, as_attachment=True)
		#return redirect("https://trips.ucc.ie/downloads")




@app.route('/uploads/')
@login_required
def uploadspage():
	global local
	try:
		print local
	except:
		local = False
	organism_dict = {}
	
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	user = current_user.name
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	org_id_dict = {}
	cursor.execute("SELECT organism_name,transcriptome_list,organism_id from organisms where private = 0 OR owner = {};".format(user_id))
	result = cursor.fetchall()

	for row in result:
		organism = row[0]
		transcriptome = row[1]
		if organism not in organism_dict:
			organism_dict[organism] = [transcriptome]
		else:
			organism_dict[organism].append(transcriptome)
		org_id_dict[row[2]] = [row[0],transcriptome]

	study_dict = {}
	cursor.execute("SELECT study_id,study_name,organism_id from studies where owner = {}".format(user_id))
	result = cursor.fetchall()
	for row in result:
		study_dict[int(row[0])] = [row[1].replace("_{}".format(user_id),"",1),org_id_dict[row[2]][0],org_id_dict[row[2]][1],[]]

	transcriptome_dict = {}
	cursor.execute("SELECT organism_id,organism_name,transcriptome_list from organisms where owner = {}".format(user_id))
	result = cursor.fetchall()
	for row in result:
		transcriptome_dict[int(row[0])] = [row[1],row[2]]

	cursor.execute("SELECT username,study_access from users")
	result = cursor.fetchall()
	for row in result:
		study_access_list = row[1].split(",")

		for study_id in study_access_list:
			if study_id == '':
				continue
			if int(study_id) in study_dict:
				study_dict[int(study_id)][3].append(row[0])

	file_dict = {}

	cursor.execute("SELECT file_name,study_id,file_id from files where owner = {}".format(user_id))
	result = cursor.fetchall()
	for row in result:
		print "study id", row[1]
		cursor.execute("SELECT study_name from studies where study_id = {}".format(row[1]))
		study_name = (cursor.fetchone())
		file_dict[row[0]] = [study_name[0].replace("_{}".format(user_id),"",1),row[2]]

	seq_dict = {}
	cursor.execute("SELECT seq_name,frame_breakdown from seq_rules where user_id = {}".format(user_id))
	result = cursor.fetchall()
	for row in result:
		seq_dict[result[0][0]] = [row[1]]
		
	
		
		
	connection.close()

	return render_template('uploads.html',
						   local=local,
						   user=user,
						   organism_dict=organism_dict,
						   study_dict=study_dict,
						   transcriptome_dict=transcriptome_dict,
						   file_dict=file_dict,
						   seq_dict=seq_dict
						   )



# file_id          | int(6)
 # organism_id      | int(5)
 # study_id         | int(6)
 # file_name        | varchar(300)
 # file_description | varchar(300)
 # file_type        | varchar(300)
#CREATE TABLE user_files (file_id INT(6), organism_id INT(5), study_id INT(6), file_name VARCHAR(300), file_description VARCHAR(300), file_type VARCHAR(300),
#owner INT(6));



@app.route('/uploadquery', methods = ['GET', 'POST'])
@login_required
def upload_file():
	if request.method == 'POST':
		uploaded_files = request.files.getlist("file")
		for f in uploaded_files:
			connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
			connection.text_factory = str
			cursor = connection.cursor()
			user = current_user.name
			cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
			result = (cursor.fetchone())
			user_id = result[0]
			filename = secure_filename(f.filename)
			file_ext = (filename.split("."))[-1]
			if file_ext != "sqlite":
				flash("Error: File extension should be sqlite not {}".format(file_ext))
				return redirect("https://trips.ucc.ie/uploads")
			if user != "public":
				foldername = "{}_{}".format(request.form["foldername"],user_id)
			else:
				foldername = "{}".format(request.form["foldername"])
			organism = request.form["organism"]
			assembly = request.form["assembly"]
			filetype = (request.form["seq_type"]).lower().strip()
			
			#if this filetype is new for this user insert a new entry into seq_rules table
			if filetype != "riboseq" and filetype != "rnaseq":
				cursor.execute("SELECT * from seq_rules where user_id = {} and seq_name = '{}'".format(user_id,filetype))
				result = cursor.fetchone()
				if result == None:
					cursor.execute("INSERT INTO seq_rules VALUES ({},'{}',0)".format(user_id,filetype))


			
			
			if not os.path.isdir("{}/uploads/{}".format(config.SCRIPT_LOC, foldername)):
				os.makedirs("{}/uploads/{}".format(config.SCRIPT_LOC, foldername))
			upload_file_path = "{}/uploads/{}/{}".format(config.SCRIPT_LOC, foldername,filename)
			f.save("{}/uploads/{}/{}".format(config.SCRIPT_LOC, foldername,filename))
			#f.save(secure_filename(f.filename))
			sqlite_db = SqliteDict("{}/uploads/{}/{}".format(config.SCRIPT_LOC, foldername,filename))
			try:
				file_description = sqlite_db["description"]
			except:
				file_description = "NULL"
			sqlite_db.close()
			#get file id
			cursor.execute("SELECT MAX(file_id) FROM files;")
			result = cursor.fetchone();
			new_file_id = int(result[0])+1

			#organism id
			print "ASSEMBLY IS ", assembly
			cursor.execute("SELECT organism_id FROM organisms WHERE organism_name = '{}' AND transcriptome_list = '{}';".format(organism, assembly))
			result = cursor.fetchone();
			organism_id = int(result[0])
			print "result is", result

			#get study_id
			cursor.execute("SELECT study_id FROM studies WHERE study_name = '{}' and organism_id = {} and owner = {}".format(foldername, organism_id, user_id))
			result = cursor.fetchone();
			if result != None:
				study_id = int(result[0])
			else:
				cursor.execute("SELECT MAX(study_id) FROM studies;")
				result = cursor.fetchone();
				study_id = int(result[0])+1
				if user != "public":
					cursor.execute("INSERT INTO studies VALUES({},{},'{}','{}','{}','{}','{}','{}','{}','{}','{}','{}',{},{})".format(study_id,organism_id,foldername,'NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL',1,user_id))
				else:
					cursor.execute("INSERT INTO studies VALUES({},{},'{}','{}','{}','{}','{}','{}','{}','{}','{}','{}',{},{})".format(study_id,organism_id,foldername,'NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL',0,user_id))
				#UPDATE study_access list with this study id
				cursor.execute("SELECT study_access from users WHERE user_id = {}".format(user_id))
				result = cursor.fetchone()
				study_access_list = result[0]
				if str(study_id) not in study_access_list.split(","):
					study_access_list += ",{}".format(study_id)
					cursor.execute("UPDATE users SET study_access = '{}' WHERE user_id = {}".format(study_access_list, user_id))


			cursor.execute("INSERT INTO files VALUES({},{},{},'{}','{}','{}',{})".format(new_file_id,organism_id,study_id,filename,file_description,filetype,user_id))
			connection.commit()
			connection.close()

			flash("File uploaded successfully")
		return redirect("https://trips.ucc.ie/uploads")


'''
CREATE TABLE `organisms` (
  `organism_id` integer NOT NULL PRIMARY KEY AUTOINCREMENT
,  `organism_name` varchar(300) DEFAULT NULL
,  `transcriptome_list` varchar(500) DEFAULT NULL
,  `gwips_databasename` varchar(300) DEFAULT NULL
,  `gwips_clade` varchar(300) DEFAULT NULL
,  `gwips_organism` varchar(300) DEFAULT NULL
,  `gwips_database` varchar(300) DEFAULT NULL
,  `default_transcript` varchar(300) DEFAULT NULL
,  `private` integer DEFAULT NULL
, owner INT(1));
'''




@app.route('/uploadtranscriptome', methods = ['GET', 'POST'])
@login_required
def upload_transcriptome():
	
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	user = current_user.name
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	
	if request.method == 'POST':
		organism = (request.form["organism"]).lower().strip().replace(" ","_")
		assembly = (request.form["assembly"]).lower().strip().replace(" ","_")
		default_tran = (request.form["default_tran"]).lower().strip().replace(" ","_")
		uploaded_annotation = request.files.getlist("anno_file")
		if not os.path.isdir("{}/uploads/transcriptomes/{}".format(config.SCRIPT_LOC, user_id)):
			os.makedirs("{}/uploads/transcriptomes/{}".format(config.SCRIPT_LOC, user_id))
		if not os.path.isdir("{}/uploads/transcriptomes/{}/{}".format(config.SCRIPT_LOC, user_id, organism)):
			os.makedirs("{}/uploads/transcriptomes/{}/{}".format(config.SCRIPT_LOC, user_id, organism))
		if not os.path.isdir("{}/uploads/transcriptomes/{}/{}/{}".format(config.SCRIPT_LOC, user_id, organism, assembly)):
			os.makedirs("{}/uploads/transcriptomes/{}/{}/{}".format(config.SCRIPT_LOC, user_id, organism, assembly))
		
		for f in uploaded_annotation:
			filename = secure_filename(f.filename)
			ext = filename.split(".")[-1]
			if ext != "sqlite":
				return """Error: Expecting extension sqlite but got extension {}. The file generated by the create_annotation_sqlite.py script should be uploaded here.
						This script can be gotten on the downloads page, by selecting the Scripts group.""".format(ext)
			#Instead of using filename of the uploaded file we rename it to organism_assembly.sqlite, to keep things consistent
			filename = "{}_{}.sqlite".format(organism, assembly)
			f.save("{}/uploads/transcriptomes/{}/{}/{}/{}".format(config.SCRIPT_LOC, user_id, organism, assembly,filename))
			
			cursor.execute("SELECT MAX(organism_id) from organisms;")
			max_org_id = (cursor.fetchone()[0])+1
			print "max org id", max_org_id
			
			cursor.execute("INSERT INTO organisms VALUES({},'{}','{}','NULL','NULL','NULL','NULL','{}',1,{})".format(max_org_id, organism, assembly,default_tran,user_id))
			#cursor.execute("INSERT INTO files VALUES({},{},{},'{}','{}','{}',{})".format(new_file_id,organism_id,study_id,filename,file_description,filetype,user_id))
			connection.commit()
			connection.close()

		flash("File uploaded successfully")
		#return "organism {} and assembly {}".format(organism, assembly)
			
		return redirect("https://trips.ucc.ie/uploads")


# Called by flask in case of an error in the code, returns the exception so it can be displayed to user
@app.errorhandler(500)
def handle_bad_request(e):
	return 'ERROR: '+str(e)+" please report this to tripsvizsite@gmail.com or via the contact page. "

# somewhere to login
@app.route("/user/login", methods=["GET", "POST"])

def login():
	global local
	#if user is already logged in then redirect to homepage
	if current_user.is_authenticated:
		return redirect("/")

	error=None
	if request.method == 'POST':
		username = str(request.form['username']).strip()
		password = str(request.form['password']).strip()
		#login_attempts_log = open("/home/DATA/www/FlaskApp/FlaskApp/logins.txt","a")
		#login_attempts_log.write("{},{}\n".format(username, password))
		#login_attempts_log.close()
		if recaptcha.verify() or local == True or username == "developer":
			username_dict = {}
			connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
			connection.text_factory = str
			cursor = connection.cursor()
			cursor.execute("SELECT username,password from users;")
			result = (cursor.fetchall())
			connection.close()
			for row in result:
				username_dict[row[0]] = row[1]
			if username in username_dict:
				if check_password_hash(username_dict[username],password) == True:
					id = username
					user = User(id)
					login_user(user)
					nxt = request.args.get('next')
					if nxt != None:
						if "<function login" in nxt:
							nxt = "/"
					else:
						nxt = "/"
					return redirect(nxt)
				else:
					#flash("No user with that email/password combo")
					#return Response("No user with that email/password combo")
					error = 'Either username or password incorrect. Please try again.'
					return render_template('login.html',error=error)
			else:
				#flash("No user with that email/password combo")
				#return Response("No user with that email/password combo")
				error = 'Either username or password incorrect. Please try again.'
				return render_template('login.html',error=error)
		else:
			error = 'Invalid Captcha. Please try again.'
			return render_template('login.html',error=error)
	else:
		return render_template('login.html',error=error)

# somewhere to login
@app.route("/create", methods=["GET", "POST"])

def create():

	#if user is already logged in then redirect to homepage
	if current_user.is_authenticated:
		return redirect("/")

	error=None
	if request.method == 'POST':
		username = str(request.form['username'])
		password = str(request.form['password'])
		password2 = str(request.form['password2'])
		#
		#login_attempts_log = open("/home/DATA/www/FlaskApp/FlaskApp/logins.txt","a")
		#login_attempts_log.write("{},{}\n".format(username, password))
		#login_attempts_log.close()
		if recaptcha.verify() or local == True:
			username_dict = {}
			connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
			connection.text_factory = str
			cursor = connection.cursor()
			cursor.execute("SELECT username,password from users;")
			result = (cursor.fetchall())
			connection.close()
			for row in result:
				username_dict[row[0]] = row[1]
			if username in username_dict:
				error = "Error: {} is already registered".format(username)
				return render_template('create.html',error=error)
			#if "@" not in username or "." not in username:
			#    error = "Error: {} is not a valid email address".format(username)
			#    return render_template('create.html',error=error)
			if password == "":
				error = "Password cannot be empty"
				return render_template('create.html',error=error)
			if password != password2:
				error = "Passwords do not match"
				return render_template('create.html',error=error)
			hashed_pass = generate_password_hash(password)

			connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
			connection.text_factory = str
			cursor = connection.cursor()

			cursor.execute("SELECT MAX(user_id) from users;")
			result = cursor.fetchone()
			max_user_id = int(result[0])
			user_id = max_user_id+1
			print user_id
			#Add -1 to study access list, causes problems when adding study id's later if we don't
			cursor.execute("INSERT INTO users VALUES ({},'{}','{}','-1','',0);".format(user_id, username, hashed_pass))
			#cursor.execute("INSERT INTO user_settings VALUES ({},'#F2F2F7','#FF5F5B','#FF5F5B','#9ACAFF','#FF5F5B','#90E090','#9ACAFF','#FFFF91','gray','gray','gray','gray','gray','gray');".format(user_id))
			cursor.execute("INSERT INTO user_settings VALUES ('20','20','32','18',{},'#F2F2F7','#FF5F5B','#FF5F5B','#9ACAFF','#FF5F5B','#90E090','#9ACAFF','#FFFF91','gray','gray','gray','gray','gray','gray',2,'gray');".format(user_id))


			connection.commit()
			connection.close()
			return redirect("/")




		else:
			error = 'Invalid Captcha. Please try again.'
			return render_template('create.html',error=error)
	else:
		return render_template('create.html',error=error)






# somewhere to logout
@app.route("/user/logout")
@login_required
def logout():
	logout_user()
	return redirect(login)


# callback to reload the user object
@login_manager.user_loader
def load_user(userid):
	return User(userid)

# Points to robots.txt in static folder to make the site searchable or not
@app.route('/robots.txt')
def static_from_root():
	return send_from_directory(app.static_folder, request.path[1:])

#This is the help page, linked from various other pages to explain terms on that page,
# Do not add a cache to this, it messes up the help links
@app.route('/help/')
def helppage():
	parent_acc = request.args.get('parent_acc')
	child_acc = request.args.get('child_acc')
	return render_template('help.html', parent_acc=parent_acc, child_acc=child_acc)


#This is the help page, linked from various other pages to explain terms on that page
@app.route('/help/single_gene_help/')
def singlegenehelppage():
	return render_template('single_gene_help.html')


#This is the help page, linked from various other pages to explain terms on that page
@app.route('/help/comparison_help/')
def comparisonhelppage():
	return render_template('comparison_help.html')

#This is the help page, linked from various other pages to explain terms on that page
@app.route('/help/metainfo_help/')
def metainfohelppage():
	return render_template('metainfo_help.html')


#This is the help page, linked from various other pages to explain terms on that page
@app.route('/help/differential_help/')
def differentialhelppage():
	return render_template('differential_help.html')


# Converts an integer to base62
def integer_to_base62(num):
	base = string.digits + string.lowercase + string.uppercase
	r = num % 62
	res = base[r];
	q = floor(num / 62)
	while q:
		r = q % 62
		q = floor(q / 62)
		res = base[int(r)] + res
	return res

# Converts a base62 encoded string to an integer
def base62_to_integer(base62_str):
	base = string.digits + string.lowercase + string.uppercase
	limit = len(base62_str)
	res = 0
	for i in xrange(limit):
		res = 62 * res + base.find(base62_str[i])
	return res

#This is the short url page, user supplies a short code which will be converted to a full url which user will then be redirected to
@app.route('/short/<short_code>/')
def short(short_code):
	try:
		user = current_user.name
	except:
		user = None
	#First convert short code to an integer
	integer = base62_to_integer(short_code)

	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT url from urls WHERE url_id = '{}';".format(integer))
	result = cursor.fetchone()
	if result == None:
		return "Short code not recognized."
	url = ""
	url += result[0]

	#add a keyword to the url to prevent generating another shortcode
	url += "&short={}".format(short_code)
	connection.close()
	return redirect(url)

#This is the home page it show a list of organisms as defined by trips_dict
@app.route('/')



def homepage():
	try:
		user = current_user.name
	except:
		user = None
	organism_access_list = []
	organism_list = []
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	user_id = -1
	if user != None:
		flash("You are logged in as {}".format(user))
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]
		#get a list of organism id's this user can access
		cursor.execute("SELECT organism_access from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		if result[0]:
			split_list = (result[0]).split(",")
			for org_id in split_list:
				organism_access_list.append(int(org_id))
				
	#returns a tuple with each field as a seperate string
	cursor.execute("SELECT organism_id,organism_name,private,owner from organisms;")
	# List of all rows returned
	result = (cursor.fetchall())
	for row in result:
		if row[2] == 0:
			organism_list.append(row[1])
		elif row[2] == 1:
			if row[0] in organism_access_list or row[3] == user_id:
				if row[1] not in organism_list:
					organism_list.append(row[1])
	
	organism_list.sort()



	connection.close()
	
	return render_template('landing.html',organisms=organism_list)




#This is the home page it show a list of organisms as defined by trips_dict
@app.route('/predictions/')
def predictionspage():
	studies = ["park","xu","combo","battle"]
	return render_template('index_predictions.html',studies=studies)

#show a list of transcriptomes
@app.route('/<organism>/')

def transcriptomepage(organism):
	transcriptomes = []

	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	cursor.execute("SELECT transcriptome_list from organisms WHERE organism_name = '{}';".format(organism))
	result = (cursor.fetchall())
	if result:
		for row in result:
			transcriptomes.append(row[0])

	#print "transcriptomes", transcriptomes
	transcriptomes = str(transcriptomes).strip("[]").replace("'","")
	return render_template('transcriptomes.html', transcriptomes=transcriptomes)

#show a list of plot types
@app.route('/<organism>/<transcriptome>/')

def plogpage(organism,transcriptome):
	try:
		user = current_user.name
	except:
		user = None
	return render_template('plot_types.html',current_username=user)

# Given a username and an organism returns a list of relevant studies.
def fetch_studies(username, organism, transcriptome):

	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	accepted_studies = {}
	study_access_list = []
	#get a list of organism id's this user can access

	if username != None:
		cursor.execute("SELECT study_access from users WHERE username = '{}';".format(username))
		result = (cursor.fetchone())

		if result[0]:
			split_list = (result[0]).split(",")
			for study_id in split_list:
				study_access_list.append(int(study_id))


	cursor.execute("SELECT organism_id from organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
	result = (cursor.fetchone())
	if result[0]:
		organism_id = int(result[0])

	#keys are converted from int to str as javascript will not accept a dictionary with ints for keys.
	cursor.execute("SELECT study_id,study_name,private from studies WHERE organism_id = '{}';".format(organism_id))
	result = (cursor.fetchall())
	if result != []:
		if result[0]:
			for row in result:
				if row[2] == 0:
					accepted_studies[str(row[0])] = {"filetypes":[],"study_name":row[1]}
				elif row[2] == 1:
					if row[0] in study_access_list:
						accepted_studies[str(row[0])] = {"filetypes":[],"study_name":row[1]}
	connection.close()
	return accepted_studies


# Create a dictionary of files seperated by type, this allows for file type grouping on the front end.
def fetch_files(accepted_studies):
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	accepted_files = {}
	seq_types = []
	file_id_to_name_dict = {}
	cursor.execute("SELECT file_id,study_id,file_name,file_description,file_type from files;")
	result = cursor.fetchall()

	for row in result:
		#print "row", row
		if row[4] not in accepted_files:
			accepted_files[row[4]] = {}
		#accepted file keys are strings due to javascript, so convert to string before checking
		#if the study id is not in accepted_studies skip it
		if str(row[1]) not in accepted_studies.keys():
			continue

		if str(row[1]) not in accepted_files[row[4]]:
			accepted_files[row[4]][str(row[1])] = {}

		#add the seq type to seq_types if it's not already in there
		if row[4] not in seq_types:
			seq_types.append(row[4])
		accepted_files[row[4]][str(row[1])][str(row[0])] = {"file_name":row[2].replace(".shelf",""),
															"file_description":row[3],
															"file_type":row[4]}
		accepted_studies[str(row[1])]["filetypes"].append(row[4])
		file_id_to_name_dict[str(row[0])] = row[2].replace(".shelf","")


	connection.close()
	return file_id_to_name_dict,accepted_studies,accepted_files,seq_types



#This is the interactive plot page
@app.route('/<organism>/<transcriptome>/interactive_plot/')

def interactiveplotpage(organism,transcriptome):
	#ip = request.environ['REMOTE_ADDR']
	global local
	try:
		print local
	except:
		local = False
	try:
		user = current_user.name
	except:
		user = None
	organism = str(organism)
	
	accepted_studies = fetch_studies(user, organism, transcriptome)
	file_id_to_name_dict,accepted_studies,accepted_files,seq_types = fetch_files(accepted_studies)
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()


	cursor.execute("SELECT gwips_clade,gwips_organism,gwips_database,default_transcript from organisms WHERE organism_name = '{}';".format(organism))
	result = (cursor.fetchone())
	gwips_clade = result[0]
	gwips_org = result[1]
	gwips_db = result[2]
	gwips_info = {"organism":gwips_org,
				  "clade": gwips_clade,
				  "database": gwips_db}
	default_tran = result[3]
	studyinfo_dict = fetch_study_info(organism)
	#print "studyinfo dict", studyinfo_dict
	user_transcript = request.args.get('tran')
	user_readscore = request.args.get('rs')
	user_hili = request.args.get('hili')
	user_generate_shorturl = request.args.get('genshort')
	user_files = request.args.get('files')
	user_minread = request.args.get('minread')
	user_maxread = request.args.get('maxread')
	user_dir = request.args.get('dir')
	user_line_graph = request.args.get('lg')
	user_ambig = request.args.get('ambig')
	user_cov = request.args.get('cov')
	user_nuc = request.args.get('nuc')
	user_short = request.args.get('short')
	user_crd = request.args.get('crd')


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

	try:
		user_hili_start = int(user_hili.split("_")[0])
		user_hili_stop = int(user_hili.split("_")[1])
	except:
		user_hili_start = None
		user_hili_stop = None

	try:
		user_minread = int(user_minread)
		user_maxread = int(user_maxread)
	except:
		user_minread = None
		user_maxread = None
	advanced='True'

	connection.close()
	return render_template('index.html', gwips_info=gwips_info, gwips_clade=gwips_clade, gwips_org=gwips_org, gwips_db=gwips_db,organism=organism,transcriptome=transcriptome,default_tran=default_tran,
						   user_transcript=user_transcript, user_readscore=user_readscore, user_hili_start=user_hili_start, user_hili_stop=user_hili_stop,local=local,studies_dict=accepted_studies,
						   accepted_files=accepted_files,user_files=user_files,user_ribo_studies=user_ribo_studies,user_rna_studies=user_rna_studies,user_minread=user_minread,user_maxread=user_maxread,
						   user_dir=user_dir,user_line_graph=user_line_graph,user_ambig=user_ambig,user_cov=user_cov,user_nuc=user_nuc,user_short=user_short, user_crd=user_crd,studyinfo_dict=studyinfo_dict,
						   advanced=advanced,seq_types=seq_types)


def fetch_study_info(organism):
	studyinfo_dict = {}
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT study_id,paper_authors,srp_nos,paper_year,paper_pmid,paper_link,gse_nos,adapters,paper_title,description from studies;".format(organism))
	result = (cursor.fetchall())
	connection.close()

	for row in result:
		# string as javascript cant handle numerical keys
		study_id = "study"+str(row[0])
		studyinfo_dict[study_id] = {"paper_authors":row[1],
									"srp_nos":row[2],
									"paper_year":row[3],
									"paper_pmid":row[4],
									"paper_link":row[5].strip('"'),
									"gse_nos":row[6],
									"adapters":row[7],
									"paper_title":row[8],
									"description":row[9]}

	return studyinfo_dict


#This is metainfo plot page
@app.route('/<organism>/<transcriptome>/metainfo_plot/')
def metainfo_plotpage(organism, transcriptome):
	global local
	try:
		print local
	except:
		local = False
	try:
		user = current_user.name
	except:
		user = None
	organism = str(organism)
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	accepted_studies = fetch_studies(user, organism, transcriptome)
	#print "accepted studies", accepted_studies
	file_id_to_name_dict, accepted_studies,accepted_files,seq_types = fetch_files(accepted_studies)


	cursor.execute("SELECT gwips_clade,gwips_organism,gwips_database,default_transcript from organisms WHERE organism_name = '{}';".format(organism))
	result = (cursor.fetchone())
	studyinfo_dict = fetch_study_info(organism)
	gwips_clade = result[0]
	gwips_org = result[1]
	gwips_db = result[2]
	gwips_info = {"organism":gwips_org,
				  "clade": gwips_clade,
				  "database": gwips_db}
	default_tran = result[3]

	# holds all values the user could possibly pass in the url (keywords are after request.args.get), anything not passed by user will be a string: "None"
	html_args = {"user_short":str(request.args.get('short')),
				 "user_plot_type":str(request.args.get('plot')),
				 "nuc_comp_direction":str(request.args.get('nc_dir')),
				 "nuc_comp_type":str(request.args.get('nc_type')),
				 "nuc_comp_min_readlen":str(request.args.get('nc_minreadlen')),
				 "nuc_comp_max_readlen":str(request.args.get('nc_maxreadlen')),
				 "nuc_comp_mapped":str(request.args.get('nc_mapped')),
				 "nuc_comp_unmapped":str(request.args.get('nc_unmapped')),
				 "trip_periodicity_min_readlen":str(request.args.get('tp_minreadlen')),
				 "trip_periodicity_max_readlen":str(request.args.get('tp_maxreadlen')),
				 "heatmap_dir":str(request.args.get('hm_dir')),
				 "heatmap_log_scale":str(request.args.get('hm_log')),
				 "heatmap_reverse":str(request.args.get('hm_rev')),
				 "heatmap_position":str(request.args.get('hm_pos')),
				 "heatmap_minreadlen":str(request.args.get('hm_minreadlen')),
				 "heatmap_maxreadlen":str(request.args.get('hm_maxreadlen')),
				 "heatmap_start":str(request.args.get('hm_start')),
				 "heatmap_stop":str(request.args.get('hm_stop')),
				 "heatmap_colour":str(request.args.get('hm_col')),
				 "metagene_pos":str(request.args.get('mg_pos')),
				 "metagene_minreadlen":str(request.args.get('mg_minreadlen')),
				 "metagene_maxreadlen":str(request.args.get('mg_maxreadlen')),
				 "replicate_minreads":str(request.args.get('rp_minreads')),
				 "mrna_dist_readlen_per":str(request.args.get('mdr_per')),
				 "transcriptome":str(transcriptome),
				 "maxscaleval":str(request.args.get('maxscaleval')),
				 "mrna_dist_readlen_smooth":str(request.args.get('mdr_smooth')),
				 "te_minreads":str(request.args.get('te_minreads')),
				 "te_tranlist":str(request.args.get('te_tranlist'))}


	user_files = request.args.get('files')
	if user_files != None:
		user_files = user_files.split(",")
		html_args["user_files"] = [str(x) for x in user_files]
	else:
		html_args["user_files"] = []

	user_ribo_studies = request.args.get('ribo_studies')
	if user_ribo_studies != None:
		user_ribo_studies = user_ribo_studies.split(",")
		html_args["user_ribo_studies"] = [str(x) for x in user_ribo_studies]
	else:
		html_args["user_ribo_studies"] = []
	user_rna_studies = request.args.get('rna_studies')
	if user_rna_studies != None:
		user_rna_studies = user_rna_studies.split(",")
		html_args["user_rna_studies"] = [str(x) for x in user_rna_studies]
	else:
		html_args["user_rna_studies"] = []




	connection.close()
	return render_template('metainfo_index.html', gwips_clade=gwips_clade, gwips_org=gwips_org, gwips_db=gwips_db,transcriptome=transcriptome,organism=organism,default_tran=default_tran,current_username=user,local=local,
						   studies_dict=accepted_studies,accepted_files=accepted_files,html_args=html_args,user_files=user_files,user_ribo_studies=user_ribo_studies, user_rna_studies=user_rna_studies,
						   studyinfo_dict=studyinfo_dict,seq_types=seq_types)


@app.route('/<organism>/<transcriptome>/comparison/')
def comparisonpage(organism, transcriptome):
	global local
	try:
		print local
	except:
		local = False
	#ip = request.environ['REMOTE_ADDR']
	try:
		user = current_user.name
	except:
		user = None
	organism = str(organism)
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT gwips_clade,gwips_organism,gwips_database,default_transcript from organisms WHERE organism_name = '{}';".format(organism))
	result = (cursor.fetchone())
	gwips_clade = result[0]
	gwips_org = result[1]
	gwips_db = result[2]
	gwips_info = {"organism":gwips_org,
				  "clade": gwips_clade,
				  "database": gwips_db}
	default_tran = result[3]
	studyinfo_dict = fetch_study_info(organism)
	#?files=227,228,229,#ff1f00_230,#3BFF00
	user_file_dict = {}
	if str(request.args.get('files')) != "None":


		colors = str(request.args.get('files')).split("_")
		for filelist in colors:
			all_items = filelist.split(",")
			files = []
			for item in all_items:
				if "#" in item:
					color = item
				else:
					files.append(item)
			user_file_dict[color] = files

	user_label_dict = {}
	if str(request.args.get('labels')) != "None":
		colors = str(request.args.get('labels')).split("_")
		for label in colors:
			all_items = label.split(",")
			for item in all_items:
				if "#" in item:
					color = item
				else:
					label = item
			user_label_dict[color] = label


	html_args = {"user_short":str(request.args.get('short')),
				 "user_file_dict":user_file_dict,
				 "user_label_dict":user_label_dict,
				 "transcript":str(request.args.get('transcript')),
				 "minread":str(request.args.get('minread')),
				 "maxread":str(request.args.get('maxread')),

				 "hili_start":str(request.args.get('hili_start')),
				 "hili_stop":str(request.args.get('hili_stop')),

				 "ambig":str(request.args.get('ambig')),
				 "cov":str(request.args.get('cov')),
				 "normalize":str(request.args.get('normalize')),
				 "transcriptome":str(transcriptome)}

	accepted_studies = fetch_studies(user, organism, transcriptome)
	file_id_to_name_dict,accepted_studies,accepted_files,seq_types = fetch_files(accepted_studies)

	connection.close()
	return render_template('index_compare.html', studies_dict=accepted_studies, accepted_files=accepted_files, gwips_info=gwips_info,
						gwips_clade=gwips_clade, gwips_org=gwips_org,gwips_db=gwips_db, organism=organism,transcriptome=transcriptome,
						default_tran=default_tran,local=local,html_args=html_args,file_id_to_name_dict=file_id_to_name_dict,studyinfo_dict=studyinfo_dict,
						seq_types=seq_types)






@app.route('/<organism>/<transcriptome>/differential/')
def diffpage(organism,transcriptome):
	global local
	try:
		print local
	except:
		local = False
	#ip = request.environ['REMOTE_ADDR']
	try:
		user = current_user.name
	except:
		user = None
	organism = str(organism)
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	cursor.execute("SELECT gwips_clade,gwips_organism,gwips_database,default_transcript from organisms WHERE organism_name = '{}';".format(organism))
	result = (cursor.fetchone())

	gwips_clade = result[0]
	gwips_org = result[1]
	gwips_db = result[2]
	gwips_info = {"organism":gwips_org,
				  "clade": gwips_clade,
				  "database": gwips_db}
	default_tran = result[3]

	studyinfo_dict = fetch_study_info(organism)
	# holds all values the user could possibly pass in the url (keywords are after request.args.get), anything not passed by user will be a string: "None"
	html_args = {"user_short":str(request.args.get('short')),
				 "minreads":str(request.args.get('minreads')),
				 "minzscore":str(request.args.get('minzscore')),
				 "region":str(request.args.get('region')),
				 "transcriptome":str(transcriptome)}


	html_args["riboseq_files_1"] = "None"
	html_args["riboseq_files_2"] = "None"
	html_args["rnaseq_files_1"] = "None"
	html_args["rnaseq_files_2"] = "None"

	user_files = request.args.get('riboseq_files_1')
	#print "USER FILES 1",len(user_files)
	if user_files != None:
		if len(user_files) != 0:
			user_files = user_files.split(",")
			html_args["riboseq_files_1"] = [str(x) for x in user_files]


	user_files = request.args.get('riboseq_files_2')
	if user_files != None:
		if len(user_files) != 0:
			user_files = user_files.split(",")
			html_args["riboseq_files_2"] = [str(x) for x in user_files]


	user_files = request.args.get('rnaseq_files_1')
	#print "USER files 2", len(user_files)
	if user_files != None:
		if len(user_files) != 0:
			user_files = user_files.split(",")
			html_args["rnaseq_files_1"] = [str(x) for x in user_files]


	user_files = request.args.get('rnaseq_files_2')
	if user_files != None:
		if len(user_files) != 0:
			user_files = user_files.split(",")
			html_args["rnaseq_files_2"] = [str(x) for x in user_files]


	user_labels = request.args.get('riboseq_labels_1')
	if user_labels != None:
		user_labels = user_labels.split(",")
		html_args["riboseq_labels_1"] = [str(x) for x in user_labels]
	else:
		html_args["riboseq_labels_1"] = []

	user_labels = request.args.get('riboseq_labels_2')
	if user_labels != None:
		user_labels = user_labels.split(",")
		html_args["riboseq_labels_2"] = [str(x) for x in user_labels]
	else:
		html_args["riboseq_labels_2"] = []

	user_labels = request.args.get('rnaseq_labels_1')
	if user_labels != None:
		user_labels = user_labels.split(",")
		html_args["rnaseq_labels_1"] = [str(x) for x in user_labels]
	else:
		html_args["rnaseq_labels_1"] = []

	user_labels = request.args.get('rnaseq_labels_2')
	if user_labels != None:
		user_labels = user_labels.split(",")
		html_args["rnaseq_labels_2"] = [str(x) for x in user_labels]
	else:
		html_args["rnaseq_labels_2"] = []
		
	accepted_studies = fetch_studies(user, organism, transcriptome)
	file_id_to_name_dict,accepted_studies,accepted_files,seq_types = fetch_files(accepted_studies)
	connection.close()
	
	return render_template('index_diff.html', studies_dict=accepted_studies, accepted_files=accepted_files,organism=organism, default_tran=default_tran,local=local,transcriptome=transcriptome,html_args=html_args,studyinfo_dict=studyinfo_dict,seq_types=seq_types)




def generate_short_code(data,organism,transcriptome,plot_type):
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	# We now build a url so that this plot can be recreated later on
	# /saccharomyces_cerevisiae/Gencode_v24/interactive_plot/?files=&tran=YPR122W&ribo_studies=9,10,11,12,13,14,15,16,17&rs=0.5&hili=500_1000&minread=15&maxread=100&user_dir=5&lg=T&ambig=F&cov=T&nuc=T

	url = "/{}/{}/{}/?".format(organism,transcriptome,plot_type)

	# To stop the file_list argument in the url being extremely long (and thus saving space), we can pass a riboseq_studies and rnaseq_studies arguments, these are study_ids of any study which
	# has all of its riboseq or rnaseq files checked. This part of the code determines which studies fall into that category.
	riboseq_studies = []
	rnaseq_studies = []

	cursor.execute("SELECT organism_id from organisms WHERE organism_name = '{}'".format(organism))
	result = cursor.fetchone()
	organism_id = int(result[0])
	cursor.execute("SELECT study_id from studies WHERE organism_id = {}".format(organism_id))
	result = cursor.fetchall()



	if plot_type == "interactive_plot" or plot_type == "metainfo_plot":
		url += "files="
		for row in result:
			study_id = int(row[0])

			#Now get all riboseq files that have this study_id, if all those ids are in file_list, add to riboseq studies and remove those files from file_list, do the same for rnaseq
			cursor.execute("SELECT file_id from files WHERE study_id = {} AND file_type = 'riboseq'".format(study_id))
			result = cursor.fetchall()
			all_present = True
			for row in result:
				if str(row[0]) not in data["file_list"]:
					all_present = False
			if all_present == True:
				riboseq_studies.append(study_id)
				for row in result:
					data["file_list"].remove(str(row[0]))

			cursor.execute("SELECT file_id from files WHERE study_id = {} AND file_type = 'rnaseq'".format(study_id))
			result = cursor.fetchall()
			all_present = True

			# If there are no files of that type for that study then don't bother adding it to the list
			if not result:
				all_present=False
			for row in result:
				if str(row[0]) not in data["file_list"]:
					all_present = False
			if all_present == True:
				rnaseq_studies.append(study_id)
				for row in result:
					data["file_list"].remove(str(row[0]))

		for filenum in data["file_list"]:
			url += filenum+","
		if riboseq_studies:
			url += "&ribo_studies="
			for study_id in riboseq_studies:
				url += str(study_id)+","
		if rnaseq_studies:
			url += "&rna_studies="
			for study_id in rnaseq_studies:
				url += str(study_id)+","




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
					url += "&hm_minreadlen={}".format(data["heatmap_minreadlen"])
			if "heatmap_maxreadlen" in data:
				if data["heatmap_maxreadlen"] != "None":
					url += "&hm_maxreadlen={}".format(data["heatmap_maxreadlen"])
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
			file_string += "{}_".format(color.replace("#","%23"))

		#remove the trailling _ from file_string
		file_string = file_string[:len(file_string)-1]

		label_string = ""
		for color in data["master_file_dict"]:
			label_string += "{},".format(data["master_file_dict"][color]["label"])
			# Can't use # in html args so encode as %23 instead
			label_string += "{}_".format(color.replace("#","%23"))

		#remove the trailling _ from file_string
		label_string = label_string[:len(label_string)-1]



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
		print "MASTER FILE DICT",data["master_file_dict"]
		#print "\n\n\n PLOT TYPE IS DIFFERENTIAL"
		url += "&minzscore={}".format(data["minzscore"])
		url += "&minread={}".format(data["minreads"])
		url += "&region={}".format(data["region"])
		url += "&riboseq_files_1={}".format(str(data["master_file_dict"]["riboseq1"]["file_ids"]).strip("[]").replace("'","").replace(" ",""))
		url += "&riboseq_files_2={}".format(str(data["master_file_dict"]["riboseq2"]["file_ids"]).strip("[]").replace("'","").replace(" ",""))
		url += "&rnaseq_files_1={}".format(str(data["master_file_dict"]["rnaseq1"]["file_ids"]).strip("[]").replace("'","").replace(" ",""))
		url += "&rnaseq_files_2={}".format(str(data["master_file_dict"]["rnaseq2"]["file_ids"]).strip("[]").replace("'","").replace(" ",""))
		url += "&riboseq_labels_1={}".format(str(data["master_file_dict"]["riboseq1"]["file_names"]).strip("[]").replace("'","").replace(" ",""))
		url += "&riboseq_labels_2={}".format(str(data["master_file_dict"]["riboseq2"]["file_names"]).strip("[]").replace("'","").replace(" ",""))
		url += "&rnaseq_labels_1={}".format(str(data["master_file_dict"]["rnaseq1"]["file_names"]).strip("[]").replace("'","").replace(" ",""))
		url += "&rnaseq_labels_2={}".format(str(data["master_file_dict"]["rnaseq2"]["file_names"]).strip("[]").replace("'","").replace(" ",""))

	cursor.execute("SELECT MAX(url_id) from urls;")
	result = cursor.fetchone()
	url_id = int(result[0])+1

	cursor.execute("INSERT INTO urls VALUES({},'{}')".format(url_id, url))
	connection.commit()
	short_code = integer_to_base62(url_id)
	connection.close()
	return short_code


@app.route('/query', methods=['POST'])
def query():

	tran_dict = {}
	gene_dict = {}
	ribo_user_files = {}
	data = ast.literal_eval(request.data)

	tran = data['transcript'].upper().strip()
	readscore = data['readscore']
	minread = int(data['minread'])
	maxread = int(data['maxread'])
	minfiles = int(data['minfiles'])
	organism = data['organism']
	seqhili = data['seqhili'].split(",")
	transcriptome = data['transcriptome']
	advanced =  data["advanced"]


	# Send file_list (a list of integers intentionally encoded as strings due to javascript), to be converted to a dictionary with riboseq/rnaseq lists of file paths.
	file_paths_dict = fetch_file_paths(data["file_list"],organism)

	primetype = data["primetype"]
	user_hili_start = data["user_hili_start"]
	user_hili_stop = data["user_hili_stop"]
	user_short = data["user_short"]

	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
	owner = (cursor.fetchone())[0]
	
	#transhelve = shelve.open("{0}{1}/{1}.shelf".format(trips_annotation_location,organism))
	if owner == 1:
		transhelve = SqliteDict("{0}{1}/{1}.sqlite".format(trips_annotation_location,organism), autocommit=False)
	else:
		print "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location,owner,organism,transcriptome)
		transhelve = SqliteDict("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location,owner,organism,transcriptome), autocommit=False)
	inputtran = True

	try:
		newtran = transhelve[tran]
	except:
		#print "setting input tran to false"
		inputtran = False
	if inputtran == False:
		try:
			all_trans = transhelve["genes"][tran]
			if len(all_trans) == 1:
				tran = all_trans[0]
			else:
				return_str = "TRANSCRIPTS"
				tran_list = []
				for transcript in transhelve["genes"][tran]:
					tranlen = transhelve[transcript]["length"]
					cds_start = transhelve[transcript]["cds_start"]
					cds_stop = transhelve[transcript]["cds_stop"]
					if cds_start == "NULL" or cds_stop == "NULL":
						cdslen = 0
						threeutrlen = 0
						cds_start = 0
					else:
						cdslen = cds_stop-cds_start
						threeutrlen = tranlen - cds_stop
					tran_list.append((transcript, tranlen, cds_start, cdslen, threeutrlen))
				sorted_tran_list = sorted(tran_list, key=lambda x: x[3])

				for item in sorted_tran_list[::-1]:
					return_str += (":{},{},{},{},{}".format(item[0], item[1], item[2], item[3], item[4]))
				return return_str
		except Exception as e:
			print "Error",e
			return "ERROR!{} Could not find any transcript corresponding to {}, here is an example of acceptable transcripts for this organism: {}".format(e,tran, transhelve.keys()[:5])

	transhelve.close()
	if 'varlite' in data:
		lite = "y"
	else:
		lite="n"
	if 'preprocess' in data:
		preprocess = True
	else:
		preprocess = False
	if 'uga_diff' in data:
		uga_diff = True
	else:
		uga_diff = False
	if 'color_readlen_dist' in data:
		color_readlen_dist = True
	else:
		color_readlen_dist = False
	if 'ribocoverage' in data:
		ribocoverage = True
	else:
		ribocoverage = False
	if "nucseq" in data:
		nucseq = True
	else:
		nucseq = False
	if "ambiguous" in data:
		ambiguous = "ambig"
	else:
		ambiguous = "unambig"

	if "noisered" in data:
		noisered = True
	else:
		noisered = False

	if "mismatch" in data:
		mismatch = True
	else:
		mismatch = False
	#print "ambig", ambiguous


	if user_short == "None":
		short_code = generate_short_code(data,organism,transcriptome,"interactive_plot")
	else:
		short_code = user_short


	try:
		user = current_user.name
	except:
		user = None
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	background_col = config.BACKGROUND_COL
	uga_col = config.UGA_COL
	uag_col = config.UAG_COL
	uaa_col = config.UAA_COL
	title_size = config.TITLE_SIZE
	subheading_size = config.SUBHEADING_SIZE
	axis_label_size = config.AXIS_LABEL_SIZE
	marker_size = config.MARKER_SIZE
	cds_marker_size = config.CDS_MARKER_SIZE
	cds_marker_colour = config.CDS_MARKER_COLOUR
	legend_size = config.LEGEND_SIZE


	seq_rules = {}

	#get user_id
	if user != None:
		
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]


		#get a list of organism id's this user can access
		cursor.execute("SELECT background_col,uga_col,uag_col,uaa_col,title_size,subheading_size,axis_label_size,marker_size,cds_marker_width,cds_marker_colour,legend_size from user_settings WHERE user_id = '{}';".format(user_id))
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

		print "USER IS NOT NONE results are", result
		#get rules for all custom seq types
		cursor.execute("SELECT * from seq_rules;")
		result = (cursor.fetchall())
		for row in result:
			seq_name = row[1]
			frame_breakdown = row[2]
			seq_rules[seq_name] = {"frame_breakdown":frame_breakdown}
		
		
		
		
		
		connection.close()
	#print "file paths dict", file_paths_dict
	#print "seq rules",seq_rules
	
	if tran != "":
		x = riboflask.plot_profile(tran, ambiguous, minread, maxread, lite, {} , ribocoverage, organism, readscore, noisered,primetype,preprocess,
								   minfiles,nucseq, user_hili_start, user_hili_stop,uga_diff,file_paths_dict,short_code, color_readlen_dist,
								   background_col,uga_col, uag_col, uaa_col,advanced,trips_annotation_location,seqhili,seq_rules,title_size,
								   subheading_size,axis_label_size,marker_size,transcriptome,trips_uploads_location,cds_marker_size,cds_marker_colour,
								   legend_size)
	else:
		x = "ERROR! Could not find any transcript corresponding to whatever you entered"
	return x

@app.route('/settingsquery', methods=['GET','POST'])
@login_required
def settingsquery():
	print "settings query called"
	data = ast.literal_eval(request.data)
	print data

	user = current_user.name
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	new_password = data['new_password']
	new_password2 = data['new_password2']
	curr_password = data['curr_password']
	if new_password != "" or new_password2 != "":
		if new_password != new_password2:
			flash("ERROR: New passwords dont match")
			return redirect("https://trips.ucc.ie/settings")
		curr_password_hash = generate_password_hash(curr_password)
		cursor.execute("SELECT password FROM users WHERE username = '{}'".format(user))
		result = cursor.fetchone()
		old_password_hash = result[0]
		if check_password_hash(old_password_hash,curr_password) == True:
			cursor.execute("UPDATE users SET password = '{}' WHERE username = '{}'".format(curr_password_hash, user))
			cursor.commit()
		else:
			flash("ERROR: Current password does not match")
			return redirect("https://trips.ucc.ie/settings")



	#get user_id
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]

	if "advanced" in data:
		cursor.execute("UPDATE users SET advanced = 1 WHERE user_id = '{}';".format(user_id))
		connection.commit()
	else:
		cursor.execute("UPDATE users SET advanced = 0 WHERE user_id = '{}';".format(user_id))
		connection.commit()

	print "DATA", data

	#get a list of organism id's this user can access
	cursor.execute("UPDATE user_settings SET background_col = '{}' WHERE user_id = '{}';".format(data["background_colour"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET title_size = '{}' WHERE user_id = '{}';".format(data["title_size"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET subheading_size = '{}' WHERE user_id = '{}';".format(data["subheading_size"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET axis_label_size = '{}' WHERE user_id = '{}';".format(data["axis_label_size"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET marker_size = '{}' WHERE user_id = '{}';".format(data["marker_size"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET legend_size = '{}' WHERE user_id = '{}';".format(data["legend_size"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET cds_marker_width = '{}' WHERE user_id = '{}';".format(data["cds_marker_width"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET cds_marker_colour = '{}' WHERE user_id = '{}';".format(data["cds_marker_colour"].strip(),user_id))
	connection.commit()


	cursor.execute("UPDATE user_settings SET readlength_col = '{}' WHERE user_id = '{}';".format(data["readlength_colour"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET metagene_fiveprime_col = '{}' WHERE user_id = '{}';".format(data["metagene_fiveprime_colour"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET metagene_threeprime_col = '{}' WHERE user_id = '{}';".format(data["metagene_threeprime_colour"].strip(),user_id))
	connection.commit()

	cursor.execute("UPDATE user_settings SET nuc_comp_a_col = '{}' WHERE user_id = '{}';".format(data["nuc_comp_a_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET nuc_comp_t_col = '{}' WHERE user_id = '{}';".format(data["nuc_comp_t_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET nuc_comp_g_col = '{}' WHERE user_id = '{}';".format(data["nuc_comp_g_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET nuc_comp_c_col = '{}' WHERE user_id = '{}';".format(data["nuc_comp_c_col"].strip(),user_id))
	connection.commit()

	cursor.execute("UPDATE user_settings SET uag_col = '{}' WHERE user_id = '{}';".format(data["uag_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET uga_col = '{}' WHERE user_id = '{}';".format(data["uga_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET uaa_col = '{}' WHERE user_id = '{}';".format(data["uaa_col"].strip(),user_id))
	connection.commit()

	cursor.execute("UPDATE user_settings SET comp_uag_col = '{}' WHERE user_id = '{}';".format(data["comp_uag_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET comp_uga_col = '{}' WHERE user_id = '{}';".format(data["comp_uga_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET comp_uaa_col = '{}' WHERE user_id = '{}';".format(data["comp_uaa_col"].strip(),user_id))
	connection.commit()


	#result = (cursor.fetchone())
	#background_col = result[0]
	#print "Background col", background_col

	connection.close()
	flash("Settings have been updated")
	return "Settings have been updated"

@app.route('/deletequery', methods=['GET','POST'])
@login_required
def deletequery():
	data = ast.literal_eval(request.data)

	user = current_user.name
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	for key in data:
		file_id = data[key].replace("delete_","")
		cursor.execute("SELECT * FROM files WHERE file_id = {}".format(file_id))
		result = cursor.fetchone()
		study_id = result[2]
		filename = result[3]
		cursor.execute("SELECT * FROM studies WHERE study_id = {}".format(study_id))
		result = cursor.fetchone()
		study_name = result[2]
		full_path = "{}{}/{}".format(trips_uploads_location,study_name, filename)
		os.remove(full_path)
		cursor.execute("DELETE FROM files WHERE file_id = {}".format(file_id))
		connection.commit()

	connection.close()
	flash("Files have been deleted")
	return redirect("https://trips.ucc.ie/uploads")



@app.route('/deletestudyquery', methods=['GET','POST'])
@login_required
def deletestudyquery():
	data = ast.literal_eval(request.data)

	user = current_user.name
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	for study_id in data:
		studycheck = data[study_id][0]
		if studycheck.split("_")[-1] != "undefined":
			study_id = studycheck.split("_")[-1]

			#First delete all files on the server associated with this study, if there are any
			cursor.execute("SELECT * FROM files WHERE study_id = {}".format(study_id))
			result = cursor.fetchall()
			if result != None:
				for row in result:
					file_id = row[0]
					filename = row[3]
					cursor.execute("SELECT * FROM studies WHERE study_id = {}".format(study_id))
					result = cursor.fetchone()
					study_name = row[2]
					full_path = "{}{}/{}".format(trips_uploads_location,study_name, filename)
					if os.path.isfile(full_path):
						os.remove(full_path)

			#Now remove the study and the files associated with it from the db
			cursor.execute("DELETE FROM studies WHERE study_id = {}".format(study_id))
			cursor.execute("DELETE FROM files WHERE study_id = {}".format(study_id))
			connection.commit()


		else:
			study_access = data[study_id][1].split(",")
			if user not in study_access:
				study_access.append(user)
			#check study_access against a list of all users
			all_users = []

			cursor.execute("SELECT username FROM users;")
			result = cursor.fetchall()
			for row in result:
				all_users.append(row[0])
			for user in study_access:
				if user not in all_users:
					flash("Error: User {} is not registered on Trips-Viz".format(user))
					return str(get_flashed_messages())
					#return render_template()
					#return redirect(url_for('uploads'))

	connection.close()
	flash("Update successful")
	return redirect("https://trips.ucc.ie/uploads")


@app.route('/deletetranscriptomequery', methods=['GET','POST'])
@login_required
def deletetranscriptomequery():
	data = ast.literal_eval(request.data)
	user = current_user.name
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	for organism_id in data:
		transcriptomecheck = data[organism_id][0]
		if transcriptomecheck.split("_")[-1] != "undefined":
			organism_id = transcriptomecheck.split("_")[-1]
			#Delete the annotation sqlite file 
			cursor.execute("SELECT organism_name, transcriptome_list FROM organisms WHERE organism_id = {}".format(organism_id))
			result = cursor.fetchone()
			sqlite_path = "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location, user_id, result[0], result[1] )
			#print "deleting sqlite path", sqlite_path
			if os.path.isfile(sqlite_path):
				os.remove(sqlite_path)
				
			sqlite_dir = "{0}transcriptomes/{1}/{2}/{3}".format(trips_uploads_location, user_id, result[0],result[1])
			if os.path.isdir(sqlite_dir):
				os.rmdir(sqlite_dir)
			
			#This removes the organism directory, but needs to handle cases where two transcriptomes are uplaoded under the same organism. 
			#sqlite_dir = "{0}transcriptomes/{1}/{2}".format(trips_uploads_location, user_id, result[0])
			#if os.path.isdir(sqlite_dir):
			#	os.rmdir(sqlite_dir)

				
				
				
			cursor.execute("DELETE FROM organisms WHERE organism_id = {}".format(organism_id))
			
			#delete all files on the server associated with this organism, if there are any
			cursor.execute("SELECT * FROM files WHERE organism_id = {}".format(organism_id))
			result = cursor.fetchall()
			#print "ALL FILES RESULT", result
			study_ids = []
			study_names = []
			if result != None:
				for row in result:
					file_id = row[0]
					filename = row[3]
					study_id = row[2]
					study_ids.append(study_id)
					cursor.execute("SELECT * FROM studies WHERE study_id = {}".format(study_id))
					result = cursor.fetchone()
					study_name = result[2]
					study_names.append(study_name)
					#print row, study_name
					full_path = "{}{}/{}".format(trips_uploads_location,study_name, filename)
					#print "DELETING THE FILE IN THE PATH {}".format(full_path)
					if os.path.isfile(full_path):
						os.remove(full_path)
			
			for study_name in study_names:
					full_path = "{}{}".format(trips_uploads_location,study_name)
					#print "DELETING THE FILE IN THE PATH {}".format(full_path)
					if os.path.isdir(full_path):
						os.rmdir(full_path)
				
			#Now remove the study and the files associated with it from the db
			#print "RUNNING QUERY DELETE FROM studies WHERE study_id = {}".format(study_id)
			#print "DELETE FROM files WHERE study_id = {}".format(study_id)
			for study_id in study_ids:
				cursor.execute("DELETE FROM studies WHERE study_id = {}".format(study_id))
				cursor.execute("DELETE FROM files WHERE study_id = {}".format(study_id))
			connection.commit()
	connection.close()
	flash("Update successful")
	return redirect("https://trips.ucc.ie/uploads")


@app.route('/seqrulesquery', methods=['GET','POST'])
@login_required
def seqrulesquery():
	data = ast.literal_eval(request.data)
	print "DATA", data
	user = current_user.name
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	#get user_id
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]

	for seq_type in data:
		if data[seq_type][0] == 'False':
			cursor.execute("UPDATE seq_rules SET frame_breakdown = 0 WHERE seq_name = '{}' and user_id = {};".format(seq_type, user_id))
		elif data[seq_type][0] == 'True':
			cursor.execute("UPDATE seq_rules SET frame_breakdown = 1 WHERE seq_name = '{}' and user_id = {};".format(seq_type, user_id))

	connection.commit()
	connection.close()
	flash("Update successful")
	return redirect("https://trips.ucc.ie/uploads")



@app.route('/comparequery', methods=['POST'])
def comparequery():
	tran_dict = {}
	data = ast.literal_eval(request.data)

	tran = data['transcript'].upper().strip()
	organism = data['organism']
	transcriptome = data['transcriptome']
	#transhelve = shelve.open("{0}{1}/{1}.shelf".format(trips_annotation_location,organism))
	
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
	
	owner = (cursor.fetchone())[0]
	#transhelve = shelve.open("{0}{1}/{1}.shelf".format(trips_annotation_location,organism))
	if owner == 1:
		transhelve = SqliteDict("{0}{1}/{1}.sqlite".format(trips_annotation_location,organism), autocommit=False)
	else:
		#print "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location,owner,organism,transcriptome)
		transhelve = SqliteDict("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location,owner,organism,transcriptome), autocommit=False)
	inputtran = True
	try:
		newtran = transhelve[tran]
	except:
		inputtran = False
	if inputtran == False:
		try:
			all_trans = transhelve["genes"][tran]
			if len(all_trans) == 1:
				tran = all_trans[0]
			else:
				return_str = "TRANSCRIPTS"
				for transcript in transhelve["genes"][tran]:
					tranlen = transhelve[transcript]["length"]
					cds_start = transhelve[transcript]["cds_start"]
					cds_stop = transhelve[transcript]["cds_stop"]
					if cds_start == "NULL" or cds_stop == "NULL":
						cdslen = "NULL"
						threeutrlen = "NULL"
					else:
						cdslen = cds_stop-cds_start
						threeutrlen = tranlen - cds_stop
					return_str += (":{},{},{},{},{}".format(transcript, tranlen, cds_start, cdslen, threeutrlen))
				return return_str
		except Exception as e:
			print "Error",e
			return "ERROR! Could not find any transcript corresponding to {}".format(tran)
	transhelve.close()
	minread = int(data['minread'])
	maxread = int(data['maxread'])
	hili_start = int(data['hili_start'])
	hili_stop = int(data['hili_stop'])
	master_filepath_dict = {}
	master_file_dict = data['master_file_dict']
	# This section is purely to sort by label alphabetically

	if master_file_dict == {}:
		return "Error: No files in the File list box. To add files to the file list box click on a study in the studies section above. This will populate the Ribo-seq and RNA-Seq sections with a list of files. Click on one of the files and then press the  Add button in the studies section. This will add the file to the File list box. Selecting another file and clicking Add again will add the new file to the same group in the File list. Alternatively to add a new group simply change the selected colour (by clicking on the coloured box in the studies section) and then click the Add file button."

	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	for color in master_file_dict:
		master_filepath_dict[color] = {"filepaths":[],"file_ids":[],"file_names":[],"file_descs":[],"mapped_reads":0}
		for file_id in master_file_dict[color]["file_ids"]:
			cursor.execute("SELECT file_name,file_description from files WHERE file_id = {};".format(file_id))
			result = (cursor.fetchone())
			file_name = master_file_dict[color]["label"] #result[0].replace(".shelf","")
			file_paths = fetch_file_paths([file_id],organism)

			for filetype in file_paths:
				for file_id in file_paths[filetype]:
					filepath = file_paths[filetype][file_id]
					#sqlite_db = SqliteDict(filename, autocommit=False)
					if os.path.isfile(filepath):
						sqlite_db = SqliteDict(filepath, autocommit=False)
					else:
						return "File not found, please report this to tripsvizsite@gmail.com or via the contact page."
					if "noncoding_counts" in sqlite_db and "coding_counts" in sqlite_db:
						master_filepath_dict[color]["mapped_reads"] += float(sqlite_db["noncoding_counts"])
						master_filepath_dict[color]["mapped_reads"] += float(sqlite_db["coding_counts"])
					else:
						if "normalize" in data:
							return "One or more selected files is missing values for 'coding_counts' and 'non_coding_counts' so cannot normalize with these files, please report this to tripsvizsite@gmail.com or via the contact page."
					master_filepath_dict[color]["filepaths"].append(filepath)
					master_filepath_dict[color]["file_ids"].append(file_id)
					master_filepath_dict[color]["file_names"].append(file_name)
					master_filepath_dict[color]["file_descs"].append(result[1])

	if 'ribocoverage' in data:
		ribocoverage = True
	else:
		ribocoverage = False
	if "ambiguous" in data:
		ambiguous = "ambig"
	else:
		ambiguous = "unambig"
	if "normalize" in data:
		normalize = True
	else:
		normalize = False

	html_args = data["html_args"]

	if html_args["user_short"] == "None":
		short_code = generate_short_code(data,organism,html_args["transcriptome"],"comparison")
	else:
		short_code = html_args["user_short"]

	try:
		user = current_user.name
	except:
		user = None
	#set colours to default values, if user logged in these will be overwritten
	background_col = config.BACKGROUND_COL
	comp_uga_col = config.UGA_COL
	comp_uag_col = config.UAG_COL
	comp_uaa_col = config.UAA_COL
	title_size = config.TITLE_SIZE
	subheading_size = config.SUBHEADING_SIZE
	axis_label_size = config.AXIS_LABEL_SIZE
	marker_size = config.MARKER_SIZE
	cds_marker_size = config.CDS_MARKER_SIZE
	cds_marker_colour = config.CDS_MARKER_COLOUR
	legend_size = config.LEGEND_SIZE
	if user != None:
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]
		
		
		
		#get a list of organism id's this user can access
		cursor.execute("SELECT background_col,comp_uga_col,comp_uag_col,comp_uaa_col,title_size,subheading_size,axis_label_size,marker_size,cds_marker_width,cds_marker_colour,legend_size from user_settings WHERE user_id = '{}';".format(user_id))
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
		connection.close()
		

	if tran != "":
		x = riboflask_compare.plot_profile(tran, ambiguous, minread, maxread, master_filepath_dict, "y", {}, ribocoverage, organism,normalize,short_code,background_col,hili_start,
											hili_stop,comp_uag_col,comp_uga_col,comp_uaa_col,trips_annotation_location,title_size, subheading_size,axis_label_size, marker_size,cds_marker_size,cds_marker_colour,
											legend_size)
	else:
		x = "ERROR! Could not find any transcript corresponding to whatever you entered"
	return x





#Given either two or four filepaths, calculates a z-score, places the z-scores in a master dict 
def calculate_zscore(riboseq1_filepath, riboseq2_filepath, rnaseq1_filepath, rnaseq2_filepath,master_dict, longest_tran_list, mapped_reads_norm,label,region,traninfo_dict, minreads, minzscore):
	print "Calculate z-score called with following filepaths {},{},{},{}".format(riboseq1_filepath, riboseq2_filepath, rnaseq1_filepath, rnaseq2_filepath)
	
	riboseq1_tot_reads = 0.001
	riboseq2_tot_reads = 0.001
	rnaseq1_tot_reads = 0.001
	rnaseq2_tot_reads = 0.001
	transcript_dict ={}
	
	groupname = ""
	
	
	
	if riboseq1_filepath:
		groupname += (riboseq1_filepath.split("/")[-1]).replace(".sqlite","")+"_"
		print "riboseq1_filepath", riboseq1_filepath
		if os.path.isfile(riboseq1_filepath):
			sqlite_db = SqliteDict(riboseq1_filepath, autocommit=False)
		else:
			return "File not found"
		#print "shelvename", shelvename, sqlite_db.keys()[0]
		if region == "fiveprime":
			opendict =sqlite_db["unambiguous_fiveprime_totals"]
		elif region == "cds":
			opendict =sqlite_db["unambiguous_cds_totals"]
		elif region == "threeprime":
			opendict =sqlite_db["unambiguous_threeprime_totals"]
		elif region == "all":
			opendict =sqlite_db["unambiguous_all_totals"]

		if mapped_reads_norm == True:
			riboseq1_tot_reads += float(sqlite_db["noncoding_counts"])
			riboseq1_tot_reads += float(sqlite_db["coding_counts"])
		sqlite_db.close()
		for transcript in longest_tran_list:
			#print "transcript",transcript
			if transcript not in transcript_dict:
				transcript_dict[transcript] = {"riboseq1":0.0001,"riboseq2":0.0001, "rnaseq1":0.0001, "rnaseq2":0.0001}
			if transcript in opendict:
				transcript_dict[transcript]["riboseq1"] += opendict[transcript]
	if riboseq2_filepath:
		groupname += (riboseq2_filepath.split("/")[-1]).replace(".sqlite","")+"_"

		if os.path.isfile(riboseq2_filepath):
			sqlite_db = SqliteDict(riboseq2_filepath, autocommit=False)
		else:
			return "File not found, please report this to tripsvizsite@gmail.com or via the contact page."
		if region == "fiveprime":
			opendict =sqlite_db["unambiguous_fiveprime_totals"]
		elif region == "cds":
			opendict =sqlite_db["unambiguous_cds_totals"]
		elif region == "threeprime":
			opendict =sqlite_db["unambiguous_threeprime_totals"]
		elif region == "all":
			opendict =sqlite_db["unambiguous_all_totals"]
		if mapped_reads_norm == True:
			riboseq2_tot_reads += float(sqlite_db["noncoding_counts"])
			riboseq2_tot_reads += float(sqlite_db["coding_counts"])
		sqlite_db.close()
		for transcript in longest_tran_list:
			if transcript not in transcript_dict:
				transcript_dict[transcript] = {"riboseq1":0.0001,"riboseq2":0.0001, "rnaseq1":0.0001, "rnaseq2":0.0001}
			if transcript in opendict:
				transcript_dict[transcript]["riboseq2"] += opendict[transcript]

	if rnaseq1_filepath:
		groupname += (rnaseq1_filepath.split("/")[-1]).replace(".sqlite","")+"_"
		if os.path.isfile(rnaseq1_filepath):
			sqlite_db = SqliteDict(rnaseq1_filepath, autocommit=False)
		else:
			return "File not found, please report this to tripsvizsite@gmail.com or via the contact page."
		if region == "fiveprime":
			opendict =sqlite_db["unambiguous_fiveprime_totals"]
		elif region == "cds":
			opendict =sqlite_db["unambiguous_cds_totals"]
		elif region == "threeprime":
			opendict =sqlite_db["unambiguous_threeprime_totals"]
		elif region == "all":
			opendict =sqlite_db["unambiguous_all_totals"]
		if mapped_reads_norm == True:
			rnaseq1_tot_reads += float(sqlite_db["noncoding_counts"])
			rnaseq1_tot_reads += float(sqlite_db["coding_counts"])
		sqlite_db.close()
		for transcript in longest_tran_list:
			if transcript not in transcript_dict:
				transcript_dict[transcript] = {"riboseq1":0.0001,"riboseq2":0.0001, "rnaseq1":0.0001, "rnaseq2":0.0001}
			if transcript in opendict:
				transcript_dict[transcript]["rnaseq1"] += opendict[transcript]
	if rnaseq2_filepath:
		groupname += (rnaseq2_filepath.split("/")[-1]).replace(".sqlite","")+"_"
		if os.path.isfile(rnaseq2_filepath):
			sqlite_db = SqliteDict(rnaseq2_filepath, autocommit=False)
		else:
			return "File not found, please report this to tripsvizsite@gmail.com or via the contact page."
		if region == "fiveprime":
			opendict =sqlite_db["unambiguous_fiveprime_totals"]
		elif region == "cds":
			opendict =sqlite_db["unambiguous_cds_totals"]
		elif region == "threeprime":
			opendict =sqlite_db["unambiguous_threeprime_totals"]
		elif region == "all":
			opendict =sqlite_db["unambiguous_all_totals"]
		if mapped_reads_norm == True:
			rnaseq2_tot_reads += float(sqlite_db["noncoding_counts"])
			rnaseq2_tot_reads += float(sqlite_db["coding_counts"])
		sqlite_db.close()
		for transcript in longest_tran_list:
			if transcript not in transcript_dict:
				transcript_dict[transcript] = {"riboseq1":0.0001,"riboseq2":0.0001, "rnaseq1":0.0001, "rnaseq2":0.0001}
			if transcript in opendict:
				transcript_dict[transcript]["rnaseq2"] += opendict[transcript]

	current_min_reads_list = []
	diff_expressed = []

	if mapped_reads_norm == True:
		ribo1_modifier = riboseq2_tot_reads/riboseq1_tot_reads
		if ribo1_modifier < 1:
			ribo1_modifier = 1

		ribo2_modifier = riboseq1_tot_reads/riboseq2_tot_reads
		if ribo2_modifier < 1:
			ribo2_modifier = 1

		rna1_modifier = rnaseq2_tot_reads/rnaseq1_tot_reads
		if rna1_modifier < 1:
			rna1_modifier = 1

		rna2_modifier = rnaseq1_tot_reads/rnaseq2_tot_reads
		if rna2_modifier < 1:
			rna2_modifier = 1
	else:
		ribo1_modifier = 1
		ribo2_modifier = 1
		rna1_modifier = 1
		rna2_modifier = 1


	
	
	
	
	#transcript_dict[transcript]["riboseq1"] = log(transcript_dict[transcript]["riboseq1"],2)
	#transcript_dict[transcript]["riboseq2"] = log(transcript_dict[transcript]["riboseq2"],2)
	#transcript_dict[transcript]["rnaseq1"] = log(transcript_dict[transcript]["rnaseq1"],2)
	#transcript_dict[transcript]["rnaseq2"] = log(transcript_dict[transcript]["rnaseq2"],2)


	for transcript in transcript_dict:
		gene = "unknown"
		if transcript in traninfo_dict:
			if "gene" in traninfo_dict[transcript]:
				gene = traninfo_dict[transcript]["gene"]
		#print "LLLLLlabel", label
		if label == "TE":
			current_min_reads = min(transcript_dict[transcript]["riboseq1"],transcript_dict[transcript]["riboseq2"],transcript_dict[transcript]["rnaseq1"],transcript_dict[transcript]["rnaseq2"])
			skip=False
			if current_min_reads < minreads:
				skip=True
			product = abs(transcript_dict[transcript]["riboseq1"])*abs(transcript_dict[transcript]["riboseq2"])*abs(transcript_dict[transcript]["rnaseq1"])*abs(transcript_dict[transcript]["rnaseq2"])
			#print transcript_dict[transcript]["riboseq1"]
			#print transcript_dict[transcript]["riboseq2"]
			#print transcript_dict[transcript]["rnaseq1"]
			#print transcript_dict[transcript]["rnaseq2"]
			#print product
			geometric_mean = log(product**(1/float(4)),2)
			te1 = (float(transcript_dict[transcript]["riboseq1"])*ribo1_modifier)/(float(transcript_dict[transcript]["rnaseq1"])*rna1_modifier)
			te2 = (float(transcript_dict[transcript]["riboseq2"])*ribo2_modifier)/(float(transcript_dict[transcript]["rnaseq2"])*rna2_modifier)
			fold_change = log(te2/te1,2)
			#current_min_reads_list.append((transcript, log(current_min_reads,2), fold_change, gene))
			if skip == False:
				current_min_reads_list.append([transcript, geometric_mean, fold_change, gene,skip])
		
		elif label == "Riboseq":
			current_min_reads = min(transcript_dict[transcript]["riboseq1"],transcript_dict[transcript]["riboseq2"])
			skip=False
			if current_min_reads < minreads:
				skip=True
			product = abs(transcript_dict[transcript]["riboseq1"])*abs(transcript_dict[transcript]["riboseq2"])
			geometric_mean = log(product**(1/float(2)),2)
			ribo1 = float(transcript_dict[transcript]["riboseq1"])*ribo1_modifier
			ribo2 = float(transcript_dict[transcript]["riboseq2"])*ribo2_modifier
			#print "ribo1, ribo2",ribo1, ribo2
			fold_change = log(ribo2/ribo1,2)
			#current_min_reads_list.append((transcript, log(current_min_reads,2), fold_change, gene))
			if skip == False:
				current_min_reads_list.append([transcript, geometric_mean, fold_change, gene, skip])

		elif label == "Rnaseq":
			current_min_reads = min(transcript_dict[transcript]["rnaseq1"],transcript_dict[transcript]["rnaseq2"])
			skip=False
			if current_min_reads < minreads:
				skip=True
			product = abs(transcript_dict[transcript]["rnaseq1"])*abs(transcript_dict[transcript]["rnaseq2"])
			geometric_mean = log(product**(1/float(2)),2)
			rna1 = float(transcript_dict[transcript]["rnaseq1"])*rna1_modifier
			rna2 = float(transcript_dict[transcript]["rnaseq2"])*rna2_modifier
			fold_change = log(rna2/rna1,2)
			#current_min_reads_list.append((transcript, log(current_min_reads,2), fold_change, gene))
			if skip == False:
				current_min_reads_list.append([transcript, geometric_mean, fold_change, gene, skip])




	#print "current_min_reads_list", current_min_reads_list
	positive_fc_final_list = []
	negative_fc_final_list = []
	sorted_current_min_reads_list = sorted(current_min_reads_list, key=lambda x: x[1])
	#print "sorted current min reads list",sorted_current_min_reads_list[:100]
	bin_list = []

	#print "sorted current min reads list,0,10,100,1000", sorted_current_min_reads_list[0],sorted_current_min_reads_list[10],sorted_current_min_reads_list[100],sorted_current_min_reads_list[1000]
	#for entry in sorted_current_min_reads_list:
	#	if entry[0] == "ENST00000400445":
	#		print "SCMRL",entry
	#print "SCMRL",len(sorted_current_min_reads_list)
	for i in range(0, len(sorted_current_min_reads_list), 300):
		# if we are not near the end of the list calculate mean and std dev
		if i < (len(sorted_current_min_reads_list)-300):
			#for every transcript in this bin add the min exp to bin_values list
			bin_values = []
			for x in range(i,i+300):
				bin_values.append(sorted_current_min_reads_list[x][2])
			#work out the mean and standard deviation of bin_values
			mean = np.mean(bin_values)
			#print "mean", mean
			standard_dev = np.std(bin_values)
			# Append mean and std dev to sorted_current_min_reads_list so we can work out z-score for each gene later
			for x in range(i,i+300):
				sorted_current_min_reads_list[x].append(mean)
				sorted_current_min_reads_list[x].append(standard_dev)
			y = minzscore*(standard_dev)
			threshold = y+abs(mean)
			bin_list.append([mean, standard_dev,threshold])
		else:
			bin_list.append(bin_list[-1])
			# Append mean and std dev to sorted_current_min_reads_list so we can work out z-score for each gene later
			for x in range(i,len(sorted_current_min_reads_list)):
				sorted_current_min_reads_list[x].append(bin_list[-1][0])
				sorted_current_min_reads_list[x].append(bin_list[-1][1])

	#['ENST00000304218', 13.862280403639371, -0.15002833793029144, 'HIST1H1E', 0.11898454708383407, 0.4569292631514013]
	for row in sorted_current_min_reads_list:
		tran = row[0]
		if tran not in master_dict:
			master_dict[tran] = {"skip":False}
		if row[4] == True:
			master_dict[tran]["skip"] = True
		else:
			master_dict[tran][groupname] = {"geometric_mean":row[1],"fold_change":row[2],"gene":row[3],"mean":row[5],"standard_dev":row[6]}

	transcript_dict["ribo1_modifier"] = ribo1_modifier
	transcript_dict["ribo2_modifier"] = ribo2_modifier
	transcript_dict["rna1_modifier"] = rna1_modifier
	transcript_dict["rna2_modifier"] = rna2_modifier


	return transcript_dict







@app.route('/diffquery', methods=['POST'])
def diffquery():
	#Elife paper: With the minimum expression level threshold of two reads, transcripts were sorted in ascending order and arranged in bins of size 300.
	#Each bin had transcripts with a similar number of mapped reads and was analysed independently. The mean and standard deviation of change
	#in expression of the transcripts within each bin was used to determine a Z-score for each transcript. For the remaining transcripts of
	#insufficient number to be binned (<300), the mean and standard deviation was obtained from the previous bin. The Z-score determined
	#for each transcript enabled comparison between bins.

	data = ast.literal_eval(request.data)
	#print "diffquery data is ",data
	html_args = data["html_args"]
	organism = data["organism"]
	transcriptome = data["transcriptome"]
	filename = organism+"_differential_translation_"+str(time.time())+".csv"
	script_location = "/home/DATA/www/tripsviz/tripsviz"
	csv_file = open(script_location+"/static/tmp/"+filename,"w")
	master_file_dict = data["master_file_dict"]
	
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
	owner = (cursor.fetchone())[0]
		
		
	if owner == 1:
		traninfo_dict = SqliteDict("{0}{1}/{1}.sqlite".format(trips_annotation_location,organism), autocommit=False)
	else:
		#print "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location,owner,organism,transcriptome)
		traninfo_dict = SqliteDict("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location,owner,organism,transcriptome), autocommit=False)

	
	
	
	#This will hold one or more z-scores for every transcript
	master_dict = {}

	riboseq1_filepaths = {}
	file_paths_dict = fetch_file_paths(master_file_dict["riboseq1"]["file_ids"], organism)
	for seq_type in file_paths_dict:
		for file_id in file_paths_dict[seq_type]:
			riboseq1_filepaths[file_id] = file_paths_dict[seq_type][file_id]
	riboseq2_filepaths = {}
	file_paths_dict = fetch_file_paths(master_file_dict["riboseq2"]["file_ids"], organism)
	for seq_type in file_paths_dict:
		for file_id in file_paths_dict[seq_type]:
			riboseq2_filepaths[file_id] = file_paths_dict[seq_type][file_id]
	try:
		rnaseq1_filepaths = {}
		file_paths_dict = fetch_file_paths(master_file_dict["rnaseq1"]["file_ids"], organism)
		for seq_type in file_paths_dict:
			for file_id in file_paths_dict[seq_type]:
				rnaseq1_filepaths[file_id] = file_paths_dict[seq_type][file_id]
		rnaseq2_filepaths = {}
		file_paths_dict = fetch_file_paths(master_file_dict["rnaseq2"]["file_ids"], organism)
		for seq_type in file_paths_dict:
			for file_id in file_paths_dict[seq_type]:
				rnaseq2_filepaths[file_id] = file_paths_dict[seq_type][file_id]
	except:
		rnaseq1_filepaths = {}
		rnaseq2_filepaths = {}
		
	minreads = float(data["minreads"])
	minzscore = float(data["minzscore"])
	if "mapped_reads_norm" in data:
		mapped_reads_norm = True
	else:
		mapped_reads_norm = False
	#print "riboseq 1 filepaths ", riboseq1_filepaths
	# User can decide to look at just riboseq fold-change, rnaseq fold-change or TE fold-change
	if len(master_file_dict["riboseq1"]["file_ids"]) != len(master_file_dict["riboseq2"]["file_ids"]):
		return "Error: Both Riboseq Condition boxes need to have an equal number of files."
	if len(master_file_dict["rnaseq1"]["file_ids"]) != len(master_file_dict["rnaseq2"]["file_ids"]):
		return "Error: Both RNA-Seq Condition boxes need to have an equal number of files."
	if len(master_file_dict["rnaseq1"]["file_ids"]) >= 1 and len(master_file_dict["riboseq1"]["file_ids"]) >= 1 and (len(master_file_dict["riboseq1"]["file_ids"]) != len(master_file_dict["rnaseq1"]["file_ids"])):
		return "Error: If RNA-Seq boxes are not empty they need to have an equal number of files as the riboseq boxes"
	
	if len(master_file_dict["riboseq1"]["file_ids"]) >= 1 and len(master_file_dict["riboseq2"]["file_ids"]) >= 1 and len(master_file_dict["rnaseq1"]["file_ids"]) >= 1 and len(master_file_dict["rnaseq2"]["file_ids"]) >= 1:
		label = "TE"
	elif len(master_file_dict["riboseq1"]["file_ids"]) >= 1 and len(master_file_dict["riboseq2"]["file_ids"]) >= 1 and len(master_file_dict["rnaseq1"]["file_ids"]) == 0 and len(master_file_dict["rnaseq2"]["file_ids"]) == 0:
		label = "Riboseq"
	elif len(master_file_dict["riboseq1"]["file_ids"]) == 0 and len(master_file_dict["riboseq2"]["file_ids"]) == 0 and len(master_file_dict["rnaseq1"]["file_ids"]) >= 1 and len(master_file_dict["rnaseq2"]["file_ids"]) >= 1:
		label = "Rnaseq"
	else:
		return "ERROR IMBALANCED OR NO FILES: Either all 4 boxes (RIBO-seq files 1, RIBO-seq files 2, mRNA-seq files 1, mRNA-seq files 2) must have a file associated with it OR both riboseq boxes OR both rnaseq boxes, you currently have {} files in RIBO-seq condition 1 files, {} in RIBO-seq condition 2 files, {} in mRNA-seq condition 1 files and {} in mRNA-seq condition 2 files".format(len(riboseq1_filepaths),len(riboseq2_filepaths),len(rnaseq1_filepaths),len(rnaseq2_filepaths))


	if organism == "homo_sapiens" or organism == "homo_sapiens_polio" or organism == "homo_sapiens_hepc":
		principal_db = SqliteDict("{0}/homo_sapiens/principal_isoforms_5ldr3tlr_rnaseq.sqlite".format(trips_annotation_location),autocommit=False)
		longest_tran_list = principal_db["transcripts"]
		principal_db.close()
		#print "longest tran list is", longest_tran_list
		if organism == "homo_sapiens_hepc":
			longest_tran_list.append("hepc")
		if organism == "homo_sapiens_polio":
			longest_tran_list.append("polio")
	elif organism == "mus_musculus":
		principal_db = SqliteDict("{0}/mus_musculus/longest_cds_transcripts.sqlite".format(trips_annotation_location),autocommit=False)
		longest_tran_list = principal_db["transcripts"]
		principal_db.close()

	else:
		longest_tran_list = traninfo_dict.keys()
		
	no_groups = max(len(master_file_dict["riboseq1"]["file_ids"]),len(master_file_dict["rnaseq1"]["file_ids"]))

	region = data["region"] #can be cds,fiveprime, or threeprime



	riboseq1_tot_reads = 0.001
	riboseq2_tot_reads = 0.001
	rnaseq1_tot_reads = 0.001
	rnaseq2_tot_reads = 0.001
	csv_file.write("Transcript,Gene,")
	# Group together files according to the order they were given in, if a mismatch in group number raise an error. 
	for i in range(0,no_groups):
		if label == "TE":
			csv_file.write("Group,Fold change (log2),Geometric Mean (log2), Bin mean, Bin standard deviation, Z-score,RiboSeq_Cond1_count,RiboSeq_Cond2_count,mRNA-Seq_Cond1_count,mRNA-Seq_Cond2_count,")
			if mapped_reads_norm:
				csv_file.write("Normalised_RiboSeq_Cond1_count,Normalised_RiboSeq_Cond2_count,Normalised_mRNA-Seq_Cond1_count,Normalised_mRNA-Seq_Cond2_count,")
		elif label == "Riboseq":
			csv_file.write("Group,Fold change (log2),Geometric Mean (log2), Bin mean, Bin standard deviation, Z-score,RiboSeq_Cond1_count,RiboSeq_Cond2_count,")
			if mapped_reads_norm:
				csv_file.write("Normalised_RiboSeq_Cond1_count,Normalised_RiboSeq_Cond2_count,")
		else:
			csv_file.write("Group,Fold change (log2),Geometric Mean (log2), Bin mean, Bin standard deviation, Z-score,mRNA-Seq_Cond1_count,mRNA-Seq_Cond2_count,")
			if mapped_reads_norm:
				csv_file.write("Normalised_mRNA-Seq_Cond1_count,Normalised_mRNA-Seq_Cond2_count,")

		if len(riboseq1_filepaths) != 0:
			riboseq1_filepath = riboseq1_filepaths[int(master_file_dict["riboseq1"]["file_ids"][i])]
			riboseq2_filepath = riboseq2_filepaths[int(master_file_dict["riboseq2"]["file_ids"][i])]
		else:
			riboseq1_filepath = ""
			riboseq2_filepath = ""
			
		if len(rnaseq1_filepaths) != 0:
			rnaseq1_filepath = rnaseq1_filepaths[int(master_file_dict["rnaseq1"]["file_ids"][i])]
			rnaseq2_filepath = rnaseq2_filepaths[int(master_file_dict["rnaseq2"]["file_ids"][i])]
		else:
			rnaseq1_filepath = ""
			rnaseq2_filepath = ""

		transcript_dict = calculate_zscore(riboseq1_filepath, riboseq2_filepath, rnaseq1_filepath, rnaseq2_filepath, master_dict, longest_tran_list, mapped_reads_norm,label,region,traninfo_dict,minreads, minzscore)
	csv_file.write("Average fold change (log2),Average Z-score")
	if no_groups == 2:
		csv_file.write(",FDR")
	csv_file.write("\n")
	aggregated_values = []
	# ['ENST00000304218', 13.862280403639371, -0.15002833793029144, 'HIST1H1E', 0.11898454708383407, 0.4569292631514013]
	# master_dict[tran][groupname] = {"geometric_mean":row[1],"fold_change":row[2],"gene":row[3],"mean":row[4],"standard_dev":row[5]}

	total_cases = 0





	if no_groups == 2:
		print "checkpoint 1"
		all_genuine_z_scores = []
		false_positives_z_dict = {}
		for i in range(0,no_groups):
			false_positives_z_dict[i] = []
		
		
		print "checkpoint 2"
		for tran in master_dict:
			if master_dict[tran]["skip"] == True:
				continue
			geo_mean_list = []
			fc_list = []
			z_scores = []
			for groupname in master_dict[tran]:
				if groupname == "skip":
					continue
				z_score = (master_dict[tran][groupname]["fold_change"]-master_dict[tran][groupname]["mean"])/master_dict[tran][groupname]["standard_dev"]
				z_scores.append(z_score)
				
				
			average_z_score = sum(z_scores)/len(z_scores)
			master_dict[tran]["z_score"] = average_z_score
			all_genuine_z_scores.append(average_z_score)
			#If z_scores not in same direction for all groups, this counts as a false positive. 
			if not all(item >= 0 for item in z_scores) and not all(item <= 0 for item in z_scores):
				#print z_scores
				for i in range(0,len(z_scores)):
					false_positives_z_dict[i].append(abs(float(z_scores[i])))
		
		
		print "checkpoint 3"
		#return "error"
		p_val_dict = {}
		all_shuffled_z_scores = []
		for i in range(0,10):
			for x in range(0,no_groups):
				shuffle(false_positives_z_dict[x])
			for y in range(0,len(false_positives_z_dict[0])):
				shuffled_z_scores = []
				for t in range(0,no_groups):
					shuffled_z_scores.append(false_positives_z_dict[t][y])
				average_z_score = sum(shuffled_z_scores)/len(shuffled_z_scores)
				all_shuffled_z_scores.append(average_z_score)
		all_shuffled_z_scores.sort()
		genuine_te_set = list(set(all_genuine_z_scores))
		temp_count = 0
		p_val_dict[0] = len(all_shuffled_z_scores)
		if all_shuffled_z_scores != []:
			for gen_te in genuine_te_set:
				if gen_te  > all_shuffled_z_scores[-1]:
					p_val_dict[gen_te] = 0
					continue

				success1s = 0
				success2s = 0
				#for ran_te in enumerate(all_shuffled_z_scores[:-1]):
				#	#print "ran_te", ran_te
				#	if gen_te <= ran_te[1]:
				#		#print "pos, shuffle_z", ran_te[0], all_shuffled_z_scores[ran_te[0]+1]
				#		success1s += 1
				#		if gen_te > all_shuffled_z_scores[ran_te[0]+1]:
				#			#print "success2 at position ",ran_te[0]+1
				#			#print "next possible te", all_shuffled_z_scores[ran_te[0]+2]
				#			success2s += 1
				#			p_val_dict[gen_te] = ran_te[0]+1
				
				pos = bisect_left(all_shuffled_z_scores, abs(gen_te))
				reversed_pos = len(all_shuffled_z_scores)-pos
				p_val_dict[abs(gen_te)] = reversed_pos
			
			print "checkpoint 4"
			#print p_val_dict.keys()
			for tran in master_dict:
				if "z_score" in master_dict[tran]:
					z_score = master_dict[tran]["z_score"]
					vala = float(len([1 for bb2 in all_genuine_z_scores if bb2 >= abs(z_score)]))
					valb = float(p_val_dict[abs(z_score)])/10
					try:
						fdr = round(valb/(vala),2)
					except:
						fdr = 0.0
					#if fdr > 1:
					#	fdr = 1
					master_dict[tran]["fdr"] = fdr
	
	
	
	
	for tran in master_dict:
		#print "Tran", tran
		if master_dict[tran]["skip"] == True:
			#print "skipping"
			continue
		#print "tran", tran
		total_cases += 1
		geo_mean_list = []
		fc_list = []
		z_scores = []
		gene = traninfo_dict[tran]["gene"]
		csv_file.write("{},{},".format(tran, gene))
		for groupname in master_dict[tran]:
			if groupname in ["skip","fdr","z_score"]:
				continue
			#Group,Fold change,Geometric Mean, Bin mean, Bin standard deviation, Z-score\n
			
			geo_mean_list.append(master_dict[tran][groupname]["geometric_mean"])
			fc_list.append(master_dict[tran][groupname]["fold_change"])
			
			fc = master_dict[tran][groupname]["fold_change"]
			
			z_score = (fc-master_dict[tran][groupname]["mean"])/master_dict[tran][groupname]["standard_dev"]
			z_scores.append(z_score)
			#if tran == "ENST00000370089":
			#	print "fold change, mean, standard_dev, z_score", fc, master_dict[tran][groupname]["mean"],master_dict[tran][groupname]["standard_dev"],z_score

			geo_mean = master_dict[tran][groupname]["geometric_mean"]
			#print "LABEL",label
			if label == "TE":
				csv_file.write("{},{},{},{},{},{},{},{},{},{},".format(groupname,fc,geo_mean,master_dict[tran][groupname]["mean"],master_dict[tran][groupname]["standard_dev"],z_score,int(transcript_dict[tran]["riboseq1"]),int(transcript_dict[tran]["riboseq2"]),int(transcript_dict[tran]["rnaseq1"]),int(transcript_dict[tran]["rnaseq2"])))
				if mapped_reads_norm:
					csv_file.write("{},{},{},{},".format(transcript_dict[tran]["riboseq1"]*transcript_dict["ribo1_modifier"]-0.0001,transcript_dict[tran]["riboseq2"]*transcript_dict["ribo2_modifier"]-0.0001,transcript_dict[tran]["rnaseq1"]*transcript_dict["rna1_modifier"]-0.0001,transcript_dict[tran]["rnaseq2"]*transcript_dict["rna2_modifier"]-0.0001))
			elif label == "Riboseq":
				csv_file.write("{},{},{},{},{},{},{},{},".format(groupname,fc,geo_mean,master_dict[tran][groupname]["mean"],master_dict[tran][groupname]["standard_dev"],z_score,int(transcript_dict[tran]["riboseq1"]),int(transcript_dict[tran]["riboseq2"])))
				if mapped_reads_norm:
					csv_file.write("{},{},".format(transcript_dict[tran]["riboseq1"]*transcript_dict["ribo1_modifier"]-0.0001,transcript_dict[tran]["riboseq2"]*transcript_dict["ribo2_modifier"]-0.0001))
			elif label == "Rnaseq":
				csv_file.write("{},{},{},{},{},{},{},{},".format(groupname,fc,geo_mean,master_dict[tran][groupname]["mean"],master_dict[tran][groupname]["standard_dev"],z_score,int(transcript_dict[tran]["rnaseq1"]),int(transcript_dict[tran]["rnaseq2"])))
				if mapped_reads_norm:
					csv_file.write("{},{},".format(transcript_dict[tran]["rnaseq1"]*transcript_dict["rna1_modifier"]-0.0001,transcript_dict[tran]["rnaseq2"]*transcript_dict["rna2_modifier"]-0.0001))
		try:
			average_geo_mean = (sum(geo_mean_list)/len(geo_mean_list))
		except:
			average_geo_mean = 0
		try:
			average_fc_mean = (sum(fc_list)/len(fc_list))
		except:
			average_fc_mean = 0
		average_z_score = sum(z_scores)/len(z_scores)
		master_dict[tran]["z_score"] = average_z_score
		
		fc = average_fc_mean
		aggregated_values.append([tran,average_geo_mean,fc,gene])
		csv_file.write("{},{}".format(fc, average_z_score))
		if no_groups == 2:
			if "fdr" in master_dict[tran]:
				csv_file.write(",{}".format(master_dict[tran]["fdr"]))
		csv_file.write("\n")

					
	sorted_aggregated_values = sorted(aggregated_values, key=lambda x: x[1])
	#print "sorted aggregative values 0,1,100,1000", sorted_aggregated_values[0], sorted_aggregated_values[1],sorted_aggregated_values[100], sorted_aggregated_values[1000]
	#print "SAV leng", len(sorted_aggregated_values)
	#for entry in sorted_aggregated_values:
	#	if entry[0] == "ENST00000400445":
	#		print "SORTED AGGREGATED VALUES",entry
	#print "sorted current min reads list",sorted_aggregated_values[:100]
	bin_list = []
	#print "SAV",sorted_aggregated_values

	for i in range(0, len(sorted_aggregated_values), 300):
		# if we are not near the end of the list calculate mean and std dev
		if i < (len(sorted_aggregated_values)-300):
			#for every transcript in this bin add the min exp to bin_values list
			bin_values = []
			for x in range(i,i+300):
				bin_values.append(sorted_aggregated_values[x][2])
			#work out the mean and standard deviation of bin_values
			mean = np.mean(bin_values)
			standard_dev = np.std(bin_values)
			# Append mean and std dev to sorted_aggregated_values so we can work out z-score for each gene later
			for x in range(i,i+300):
				sorted_aggregated_values[x].append(mean)
				sorted_aggregated_values[x].append(standard_dev)
			y = minzscore*(standard_dev)
			upper_threshold = y+mean
			lower_threshold = (-1*y)+mean
			bin_list.append([mean, standard_dev, upper_threshold,lower_threshold])
		else:
			bin_list.append(bin_list[-1])
			# Append mean and std dev to sorted_aggregated_values so we can work out z-score for each gene later
			for x in range(i,len(sorted_aggregated_values)):
				sorted_aggregated_values[x].append(bin_list[-1][0])
				sorted_aggregated_values[x].append(bin_list[-1][1])


	try:
		user = current_user.name
	except:
		user = None
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	background_col = config.BACKGROUND_COL
	uga_col = config.UGA_COL
	uag_col = config.UAG_COL
	uaa_col = config.UAA_COL
	title_size = config.TITLE_SIZE
	subheading_size = config.SUBHEADING_SIZE
	axis_label_size = config.AXIS_LABEL_SIZE
	marker_size = config.MARKER_SIZE

	#get user_id

	if user != None:
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]

		#get a list of organism id's this user can access
		cursor.execute("SELECT background_col,title_size,subheading_size,axis_label_size,marker_size from user_settings WHERE user_id = '{}';".format(user_id))
		result = (cursor.fetchone())
		background_col = result[0]
		title_size = result[1]
		subheading_size = result[2]
		axis_label_size = result[3]
		marker_size = result[4]
		connection.close()

	if html_args["user_short"] == "None":
		short_code = generate_short_code(data,organism,html_args["transcriptome"],"differential")
	else:
		short_code = html_args["user_short"]

	return riboflask_diff.plot_profile(sorted_aggregated_values,
										bin_list,
										organism,
										label,
										transcriptome,
										master_file_dict["riboseq1"]["file_ids"],
										master_file_dict["riboseq2"]["file_ids"],
										master_file_dict["rnaseq1"]["file_ids"],
										master_file_dict["rnaseq2"]["file_ids"],
										background_col,
										short_code,
										mapped_reads_norm,
										filename,
										no_groups,
										str(title_size)+"pt",
										str(axis_label_size)+"pt", 
										str(subheading_size)+"pt",
										str(marker_size)+"pt")





def aggregate_counts(file_paths_dict, traninfo_dict, longest_tran_list, region, organism, all_seq_types, te_minimum_reads, html_args):
	#print "AGGREGATE COUNTS CALLED"
	file_list = ""
	transcript_dict = {}
	table_str = ""
	filename = organism+"_translation_efficiencies_"+str(time.time())+".csv"
	table_str += filename+"?~"
	#print "file paths dict", file_paths_dict
	for seq_type in file_paths_dict:
		#print "seq type", seq_type
		for file_id in file_paths_dict[seq_type]:
			#print "file_id", file_id
			file_list += "{},".format(file_id)
			filepath = file_paths_dict[seq_type][file_id]
			if os.path.isfile(filepath):
				sqlite_db = SqliteDict(filepath, autocommit=False)
			else:
				return "File not found, please report this to tripsvizsite@gmail.com or via the contact page."
			#print "shelvename", shelvename, sqlite_db.keys()[0]
			if region == "all":
				opendict =sqlite_db["unambiguous_all_totals"]
			elif region == "cds":
				opendict =sqlite_db["unambiguous_cds_totals"]
			elif region == "fiveprime":
				opendict =sqlite_db["unambiguous_fiveprime_totals"]
			elif region == "threeprime":
				opendict =sqlite_db["unambiguous_threeprime_totals"]
			sqlite_db.close()
			for transcript in longest_tran_list:
				#print "transcript",transcript
				if transcript not in transcript_dict:
					transcript_dict[transcript] = {}
				if seq_type not in transcript_dict[transcript]:
					transcript_dict[transcript][seq_type] = 0
				if transcript in opendict:
					transcript_dict[transcript][seq_type] += opendict[transcript]
	#print "transcript_dict", transcript_dict
	total_rows = 0
	# The next line can get the script location if running flask locally, running through apache however sys.argv[0]
	# seems to return empty, so I've hardcoded the path for now
	# This is bad programming but I have deadlines to meet :)
	# script_location =  os.path.dirname(os.path.realpath(sys.argv[0]))
	script_location = "/home/DATA/www/tripsviz/tripsviz"
	tmp_te_file = open(script_location+"/static/tmp/"+filename,"w")
	tmp_te_file.write("Filename, Gene,Transcript,Region,Riboseq count, Rnaseq count, Translation efficiency")
	for seq_type in all_seq_types:
		if seq_type != "riboseq" and seq_type != "rnaseq":
			tmp_te_file.write(",{}".format(seq_type))
	tmp_te_file.write("\n")

	all_rows = []

	for transcript in transcript_dict:
		seq_count_dict = {}
		skip = False
		for seq_type in transcript_dict[transcript]:
			if transcript_dict[transcript][seq_type] < te_minimum_reads:
				skip = True
			try:
				gene = traninfo_dict[transcript]["gene"]
			except:
				gene = "Unknown"
			#total_rows += 1

			if seq_type not in seq_count_dict:
				seq_count_dict[seq_type] = transcript_dict[transcript][seq_type]
		if skip == True:
			continue
		if "riboseq" in seq_count_dict:
			riboseq_count = seq_count_dict["riboseq"]
		else:
			riboseq_count = 0
		if "rnaseq" in seq_count_dict:
			rnaseq_count = seq_count_dict["rnaseq"]
		else:
			rnaseq_count = 0
		if rnaseq_count == 0 or riboseq_count == 0:
			te = 0
		else:
			te = float(transcript_dict[transcript]["riboseq"])/float(transcript_dict[transcript]["rnaseq"])
			te = round(te,2)
		tmp_te_file.write("Aggregate,{},{},{},{},{},{}".format(gene,transcript,region, riboseq_count, rnaseq_count, te))
		input_list = ["Aggregate", gene,transcript,region, riboseq_count, rnaseq_count, te]
		for seq_type in all_seq_types:
			if seq_type != "riboseq" and seq_type != "rnaseq":
				if seq_type in seq_count_dict:
					input_list.append(seq_count_dict[seq_type])
					tmp_te_file.write(",{}".format(seq_count_dict[seq_type]))
				else:
					input_list.append(0)
					tmp_te_file.write(",0")
		input_list.append("<a href='http://trips.ucc.ie/"+organism+"/"+html_args["transcriptome"]+"/interactive_plot/?tran="+transcript+"&files="+file_list+"' target='_blank_' >View plot</a>")
		all_rows.append(input_list)
		tmp_te_file.write("\n")
	tmp_te_file.close()
	os.chmod(script_location+"/static/tmp/"+filename, 0777)
	#if both rnaseq and riboseq files, sort by te, else sort by the relevant count
	anyfile = False
	if len(file_paths_dict["riboseq"]) != 0 and len(file_paths_dict["rnaseq"]) != 0:
		all_sorted_rows = sorted(all_rows, key=lambda x: x[6],reverse=True)
	elif len(file_paths_dict["riboseq"]) != 0:
		all_sorted_rows = sorted(all_rows, key=lambda x: x[4],reverse=True)
	elif len(file_paths_dict["rnaseq"]) != 0:
		all_sorted_rows = sorted(all_rows, key=lambda x: x[5],reverse=True)
	else:
		for seq_type in all_seq_types:
			if len(file_paths_dict[seq_type]) != 0:
				anyfile = True
				all_sorted_rows = sorted(all_rows, key=lambda x: x[6],reverse=True)
		if anyfile == False:
			return "No files selected. Select a file by clicking on a study name in the studies section. Then select one of the files that appear in the files section."
	for row in all_sorted_rows:
		total_rows += 1
		if total_rows <= 1000:
			input_str = ""
			for item in row:
				input_str += "{}.;".format(item)
			input_str += "?~"
			table_str += input_str
	table_str = "TE?~"+str(total_rows)+"?~"+table_str
	#return response
	return table_str





def sample_counts(file_paths_dict, traninfo_dict, longest_tran_list, region, organism, all_seq_types, te_minimum_reads, html_args):
	print "SAMPLE COUNTS CALLED"
	file_list = ""
	transcript_dict = {}
	table_str = ""
	filename = organism+"_translation_efficiencies_"+str(time.time())+".csv"
	table_str += filename+"?~"
	for seq_type in file_paths_dict:
		for file_id in file_paths_dict[seq_type]:
			file_list += "{},".format(file_id)
			filepath = file_paths_dict[seq_type][file_id]
			inputfilename = (filepath.split("/")[-1]).replace(".sqlite","")
			if os.path.isfile(filepath):
				sqlite_db = SqliteDict(filepath, autocommit=False)
			else:
				return "File not found, please report this to tripsvizsite@gmail.com or via the contact page."
			#print "shelvename", shelvename, sqlite_db.keys()[0]
			if region == "all":
				opendict =sqlite_db["unambiguous_all_totals"]
			elif region == "cds":
				opendict =sqlite_db["unambiguous_cds_totals"]
			elif region == "fiveprime":
				opendict =sqlite_db["unambiguous_fiveprime_totals"]
			elif region == "threeprime":
				opendict =sqlite_db["unambiguous_threeprime_totals"]
			sqlite_db.close()
			for transcript in longest_tran_list:
				#print "transcript",transcript
				if transcript not in transcript_dict:
					transcript_dict[transcript] = {}
				if seq_type not in transcript_dict[transcript]:
					transcript_dict[transcript][seq_type] = {}
				if transcript in opendict:
					transcript_dict[transcript][seq_type][inputfilename] = {"count":opendict[transcript],"file_id":str(file_id)}

	total_rows = 0
	# The next line can get the script location if running flask locally, running through apache however sys.argv[0]
	# seems to return empty, so I've hardcoded the path for now
	# This is bad programming but I have deadlines to meet :)
	# script_location =  os.path.dirname(os.path.realpath(sys.argv[0]))
	script_location = "/home/DATA/www/tripsviz/tripsviz"
	tmp_te_file = open(script_location+"/static/tmp/"+filename,"w")
	tmp_te_file.write("Filename,Gene,Transcript,Region,Riboseq count, Rnaseq count, Translation efficiency")
	for seq_type in all_seq_types:
		if seq_type != "riboseq" and seq_type != "rnaseq":
			tmp_te_file.write(",{}".format(seq_type))
	tmp_te_file.write("\n")

	all_rows = []

	for transcript in transcript_dict:
		seq_count_dict = {}
		for seq_type in transcript_dict[transcript]:
			for inputfilename in transcript_dict[transcript][seq_type]:
				count = transcript_dict[transcript][seq_type][inputfilename]["count"]
				file_id = transcript_dict[transcript][seq_type][inputfilename]["file_id"]
				if count < te_minimum_reads:
					continue
				try:
					gene = traninfo_dict[transcript]["gene"]
				except:
					gene = "Unknown"
				riboseq_count = 0.001
				rnaseq_count = 0.001
				if seq_type == "riboseq":
					riboseq_count = count
				if seq_type == "rnaseq":
					rnaseq_count = count
				try:
					te = riboseq_count/rnaseq_count
				except:
					te = 0
				tmp_te_file.write("{},{},{},{},{},{},{}".format(inputfilename,gene,transcript,region, riboseq_count, rnaseq_count, te))
				input_list = [inputfilename, gene,transcript,region, riboseq_count, rnaseq_count, te]
				input_list.append("<a href='http://trips.ucc.ie/"+organism+"/"+html_args["transcriptome"]+"/interactive_plot/?tran="+transcript+"&files="+file_id+"' target='_blank_' >View plot</a>")
				all_rows.append(input_list)
				tmp_te_file.write("\n")
	tmp_te_file.close()
	os.chmod(script_location+"/static/tmp/"+filename, 0777)
	#if both rnaseq and riboseq files, sort by te, else sort by the relevant count
	anyfile = False
	if len(file_paths_dict["riboseq"]) != 0 and len(file_paths_dict["rnaseq"]) != 0:
		all_sorted_rows = sorted(all_rows, key=lambda x: x[6],reverse=True)
	elif len(file_paths_dict["riboseq"]) != 0:
		all_sorted_rows = sorted(all_rows, key=lambda x: x[4],reverse=True)
	elif len(file_paths_dict["rnaseq"]) != 0:
		all_sorted_rows = sorted(all_rows, key=lambda x: x[5],reverse=True)
	else:
		for seq_type in all_seq_types:
			if len(file_paths_dict[seq_type]) != 0:
				anyfile = True
				all_sorted_rows = sorted(all_rows, key=lambda x: x[6],reverse=True)
		if anyfile == False:
			return "No files selected, please report this to tripsvizsite@gmail.com or via the contact page."
	for row in all_sorted_rows:
		total_rows += 1
		if total_rows <= 1000:
			input_str = ""
			for item in row:
				input_str += "{}.;".format(item)
			input_str += "?~"

			table_str += input_str

	table_str = "TE?~"+str(total_rows)+"?~"+table_str
	#return response
	return table_str




def fetch_rld(sqlite_db,ambig_type):
	print sqlite_db
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










@app.route('/metainfoquery', methods=['POST'])
def metainfoquery():
	start_time = time.time()
	print "Metainfo query called at {}".format(time.time()-start_time)
	tran_dict = {}
	gene_dict = {}
	data = ast.literal_eval(request.data)


	plottype = data["plottype"]
	minreadlen = int(data['minreadlen'])
	maxreadlen = int(data['maxreadlen'])
	trip_minreadlen = int(data['trip_minreadlen'])
	trip_maxreadlen = int(data['trip_maxreadlen'])
	mismatch_minreadlen = int(data['mismatch_minreadlen'])
	mismatch_maxreadlen = int(data['mismatch_maxreadlen'])
	nuc_minreadlen = int(data['nuc_minreadlen'])
	nuc_maxreadlen = int(data['nuc_maxreadlen'])
	heatmap_minreadlen = int(data['heatmap_minreadlen'])
	heatmap_maxreadlen = int(data['heatmap_maxreadlen'])
	heatmap_startpos = int(data['heatmap_startpos'])
	heatmap_endpos = int(data['heatmap_endpos'])
	maxscaleval = data['maxscaleval']
	if maxscaleval != "None" and maxscaleval != "":
		try:
			maxscaleval = int(maxscaleval)
		except:
			maxscalval = "None"

	color_palette = (data["color_palette"])
	minimum_reads = int(data['minimum_reads'])
	te_minimum_reads = data['te_minimum_reads']
	te_tranlist = data['te_tranlist'].replace(","," ")
	organism = data['organism']
	transcriptome = data['transcriptome']
	metagene_type = data["metagene_type"]
	heatmap_metagene_type = data["heatmap_metagene_type"]
	contaminant_organism = data["contaminant_organism"]
	nuc_comp_type = data["nuc_comp_type"]
	nuc_comp_frame = data["nuc_comp_frame"]
	nuc_comp_direction = data["nuc_comp_direction"]
	heatmap_direction = data["heatmap_direction"]
	smooth_amount = int(data["smooth_amount"])
	html_args = data["html_args"]
	all_seq_types = data["all_seq_types"]


	try:
		user = current_user.name

	except:
		user = None
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	background_col = config.BACKGROUND_COL
	readlength_col = config.READLENGTH_COL
	metagene_fiveprime_col = config.METAGENE_FIVEPRIME_COL
	metagene_threeprime_col = config.METAGENE_THREEPRIME_COL
	a_col = config.A_COL
	t_col = config.T_COL
	g_col = config.G_COL
	c_col = config.C_COL
	title_size = config.TITLE_SIZE
	subheading_size = config.SUBHEADING_SIZE
	axis_label_size = config.AXIS_LABEL_SIZE
	marker_size = config.MARKER_SIZE


	#get a list of organism id's this user can access
	if user != None:
		#get user_id
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]
		cursor.execute("SELECT background_col,readlength_col,metagene_fiveprime_col,metagene_threeprime_col,nuc_comp_a_col,nuc_comp_t_col,nuc_comp_g_col,nuc_comp_c_col,title_size,subheading_size,axis_label_size,marker_size from user_settings WHERE user_id = '{}';".format(user_id))
		result = (cursor.fetchone())
		background_col = result[0]
		readlength_col = result[1]
		metagene_fiveprime_col = result[2]
		metagene_threeprime_col = result[3]
		a_col = result[4]
		t_col = result[5]
		g_col = result[6]
		c_col = result[7]
		title_size = result[8]
		subheading_size = result[9]
		axis_label_size = result[10]
		marker_size = result[11]
		
		
		
	
	if "nuc_mapped" in data:
		nuc_mapped = True
	else:
		nuc_mapped = False

	if "nuc_unmapped" in data:
		nuc_unmapped = True
	else:
		nuc_unmapped = False
		
	if "readlen_ambig" in data:
		readlen_ambig = True
	else:
		readlen_ambig = False

	if "breakdown_per" in data:
		breakdown_per = True
	else:
		breakdown_per = False

	if "metagene_five" in data:
		metagene_five = True
	else:
		metagene_five = False

	if "metagene_three" in data:
		metagene_three = True
	else:
		metagene_three = False



	total_files = len(data["file_list"])
	# Send file_list (a list of integers intentionally encoded as strings due to javascript), to be converted to a dictionary with riboseq/rnaseq lists of file paths.
	file_paths_dict = fetch_file_paths(data["file_list"],organism)




	if "log_scale" in data:
		log_scale = True
	else:
		log_scale = False
	if "mrna_readlen_per" in data:
		mrna_readlen_per = True
	else:
		mrna_readlen_per = False
	if "count_agg" in data:
		count_agg = True
	else:
		count_agg = False

	if "reverse_scale" in data:
		reverse_scale = True
	else:
		reverse_scale = False

	if html_args["user_short"] == "None":
		short_code = generate_short_code(data,organism,html_args["transcriptome"],"metainfo_plot")
	else:
		short_code = html_args["user_short"]

	print "Print data read in  at {}".format(time.time()-start_time)

	if plottype == "readlen_dist":
		master_dict = {}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page. ".format(filepath)
				if readlen_ambig == True:
					if "read_lengths" not in sqlite_db:
						return "No readlength distribution data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
					else:
						read_lengths = sqlite_db["read_lengths"]
					sqlite_db.close()
					for i in read_lengths:
						if i in master_dict:
							master_dict[i] += read_lengths[i]
						else:
							master_dict[i] = read_lengths[i]
				elif readlen_ambig == False:
					if "unambig_read_lengths" not in sqlite_db:
						return "No unambiguous readlength distribution data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
					else:
						read_lengths = sqlite_db["unambig_read_lengths"]
					sqlite_db.close()
					for i in read_lengths:
						if i in master_dict:
							master_dict[i] += read_lengths[i]
						else:
							master_dict[i] = read_lengths[i]
					
		title = "Readlength distribution"
		connection.close()
		return metainfo_plots.readlen_dist(master_dict,title,short_code, background_col,readlength_col,title_size, axis_label_size, subheading_size,marker_size)
	elif plottype == "mismatches":
		master_dict = {}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page. ".format(filepath)
				if "global_mismatches" not in sqlite_db:
					return "No mismatches data for this file, please report this to tripsvizsite@gmail.com or via the contact page. "
				mismatches = sqlite_db["global_mismatches"]
				sqlite_db.close()
				for readlen in mismatches:
					if readlen < mismatch_minreadlen or readlen > mismatch_maxreadlen:
						continue
					for position in mismatches[readlen]:
						if position not in master_dict:
							master_dict[position] = mismatches[readlen][position]
						else:
							master_dict[position] += mismatches[readlen][position]
		title = "Mismatches"
		return metainfo_plots.mismatches(master_dict, title, short_code, background_col,title_size, axis_label_size, subheading_size,marker_size)
				
				
				
				
			
	elif plottype == "te":
		#traninfo_dict = shelve.open("{0}{1}/{1}.shelf".format(trips_annotation_location,organism))

		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		
		#transhelve = shelve.open("{0}{1}/{1}.shelf".format(trips_annotation_location,organism))
		if owner == 1:
			traninfo_dict = SqliteDict("{0}{1}/{1}.sqlite".format(trips_annotation_location,organism), autocommit=False)
		else:
			#print "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location,owner,organism,transcriptome)
			traninfo_dict = SqliteDict("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location,owner,organism,transcriptome), autocommit=False)

		#print "shelf in"
		#traninfo_dict = dict(traninfo_shelve)
		#traninfo_shelve.close()

		if total_files >= 31 and count_agg == False:
			return "Error: A maximum of 30 files can be selected if not aggregating counts. Reduce number of selected files or click the 'Aggregate counts' checkbox at the top right of the page."
			

		te_minimum_reads = int(te_minimum_reads)
		transcript_list =  te_tranlist.split(" ")
		#print "transcript list", transcript_list

		region = data["region"]


		#For human only return the principal isoforms, for everything else return everything
		if transcript_list == ['']:
			if organism == "homo_sapiens" or organism == "homo_sapiens_polio":
				longest_tran_db = SqliteDict("{0}homo_sapiens/principal_isoforms_5ldr3tlr_rnaseq.sqlite".format(trips_annotation_location),autocommit=True)
				longest_tran_list = longest_tran_db["transcripts"]
				longest_tran_db.close()
				if organism == "homo_sapiens_polio":
					longest_tran_list.append("POLIO")
			else:
				longest_tran_list = traninfo_dict.keys()
		else:
			longest_tran_list = transcript_list
			
		if count_agg == True:
			table_str = aggregate_counts(file_paths_dict, traninfo_dict, longest_tran_list, region, organism,all_seq_types, te_minimum_reads, html_args)
		elif count_agg == False:
			table_str = sample_counts(file_paths_dict, traninfo_dict, longest_tran_list, region, organism,all_seq_types, te_minimum_reads, html_args)
		return table_str




	elif plottype == "mrna_dist":
		#traninfo_shelve = shelve.open("{0}{1}/{1}.shelf".format(trips_annotation_location,organism))
		#traninfo_dict = dict(traninfo_shelve)
		#traninfo_shelve.close()
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		
		if owner == 1:
			traninfo_dict = SqliteDict("{0}{1}/{1}.sqlite".format(trips_annotation_location,organism), autocommit=False)
		else:
			#print "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location,owner,organism,transcriptome)
			traninfo_dict = SqliteDict("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location,owner,organism,transcriptome), autocommit=False)



		st_time = time.time()
		print "start time", st_time
		if organism == "homo_sapiens" or organism == "homo_sapiens_polio":
			longest_tran_db = SqliteDict("{0}homo_sapiens/principal_isoforms_5ldr3tlr_rnaseq.sqlite".format(trips_annotation_location), autocommit=False)
			longest_tran_list = longest_tran_db["transcripts"]
			longest_tran_db.close()
			if organism == "homo_sapiens_polio":
				longest_tran_list.append("POLIO")
		else:
			longest_tran_list = traninfo_dict.keys()
		mrna_dist_dict = {}

		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				filename = (filepath.split("/")[-1]).replace(".sqlite","")
				mrna_dist_dict[filename] = {"5_leader":0,
											"start_codon":0,
											"cds":0,
											"stop_codon":0,
											"3_trailer":0}
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
					#opendict = dict(sqlite_db)
					#sqlite_db.close()
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				if "mrna_dist_dict" in sqlite_db:
					mrna_dist_dict[filename]["5_leader"] = sqlite_db["mrna_dist_dict"]["5_leader"]
					mrna_dist_dict[filename]["start_codon"] = sqlite_db["mrna_dist_dict"]["start_codon"]
					mrna_dist_dict[filename]["cds"] = sqlite_db["mrna_dist_dict"]["cds"]
					mrna_dist_dict[filename]["stop_codon"] = sqlite_db["mrna_dist_dict"]["stop_codon"]
					mrna_dist_dict[filename]["3_trailer"] = sqlite_db["mrna_dist_dict"]["3_trailer"]
				else:
					for transcript in longest_tran_list:
						try:
							transcript_dict = sqlite_db[transcript]["unambig"]
						except:
							continue
							
						try:
							cds_start = int(traninfo_dict[transcript]["cds_start"])
							cds_stop = int(traninfo_dict[transcript]["cds_stop"])
						except:
							continue
						for readlen in transcript_dict:
							for five_pos in transcript_dict[readlen]:
								three_pos = five_pos+readlen
								if three_pos <= cds_start+3:
									mrna_dist_dict[filename]["5_leader"] += transcript_dict[readlen][five_pos]
								elif five_pos <= cds_start-4 and three_pos >= cds_start+4:
									mrna_dist_dict[filename]["start_codon"] += transcript_dict[readlen][five_pos]
								elif five_pos >= cds_start-3 and three_pos <= cds_stop-2:
									mrna_dist_dict[filename]["cds"] += transcript_dict[readlen][five_pos]
								elif five_pos <= cds_stop-9 and three_pos >= cds_stop-1:
									mrna_dist_dict[filename]["stop_codon"] += transcript_dict[readlen][five_pos]
								elif five_pos >= cds_stop-9:
									mrna_dist_dict[filename]["3_trailer"] += transcript_dict[readlen][five_pos]
					sqlite_db["mrna_dist_dict"] = {"5_leader":mrna_dist_dict[filename]["5_leader"],
													"start_codon":mrna_dist_dict[filename]["start_codon"],
													"cds":mrna_dist_dict[filename]["cds"],
													"stop_codon":mrna_dist_dict[filename]["stop_codon"],
													"3_trailer":mrna_dist_dict[filename]["3_trailer"]}
					sqlite_db.commit()
				sqlite_db.close()

		connection.close()
		
		fin_time = time.time() - st_time 
		print "Finish time",fin_time
		return metainfo_plots.mrna_dist(mrna_dist_dict,short_code, background_col,title_size, axis_label_size, subheading_size,marker_size)
	elif plottype == "mrna_dist_readlen":
		minreadlen = 15
		maxreadlen = 100

		#traninfo_shelve = shelve.open("{0}{1}/{1}.shelf".format(trips_annotation_location,organism))
		#traninfo_dict = dict(traninfo_shelve)
		#traninfo_shelve.close()
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		
		
		if owner == 1:
			traninfo_dict = SqliteDict("{0}{1}/{1}.sqlite".format(trips_annotation_location,organism), autocommit=False)
		else:
			#print "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location,owner,organism,transcriptome)
			traninfo_dict = SqliteDict("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(trips_uploads_location,owner,organism,transcriptome), autocommit=False)

		if organism == "homo_sapiens" or organism == "homo_sapiens_polio":
			longest_tran_db = SqliteDict("{0}homo_sapiens/principal_isoforms_5ldr3tlr_rnaseq.sqlite".format(trips_annotation_location),autocommit=True)
			longest_tran_list = longest_tran_db["transcripts"]
			longest_tran_db.close()
			if organism == "homo_sapiens_polio":
				longest_tran_list.append("POLIO")
		else:
			longest_tran_list = traninfo_dict.keys()
		mrna_dist_dict = {"5_leader":collections.OrderedDict(),
						  "start_codon":collections.OrderedDict(),
						  "cds":collections.OrderedDict(),
						  "stop_codon":collections.OrderedDict(),
						  "3_trailer":collections.OrderedDict(),
						  "total":collections.OrderedDict()}

		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				try:
					sqlite_db = SqliteDict(filepath, autocommit=False)
					#opendict = dict(sqlite_db)
					#sqlite_db.close()
				except:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)

				if "mrna_dist_readlen_dict" in sqlite_db:
					print sqlite_db["mrna_dist_readlen_dict"]
					for readlen in range(minreadlen,maxreadlen+1):
						if readlen not in mrna_dist_dict["5_leader"]:
							mrna_dist_dict['5_leader'][readlen]=0
							mrna_dist_dict['start_codon'][readlen]=0
							mrna_dist_dict['cds'][readlen]=0
							mrna_dist_dict['stop_codon'][readlen]=0
							mrna_dist_dict['3_trailer'][readlen]=0
							mrna_dist_dict['total'][readlen]=0
						mrna_dist_dict['5_leader'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["5_leader"][readlen]
						mrna_dist_dict['start_codon'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["start_codon"][readlen]
						mrna_dist_dict['cds'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["cds"][readlen]
						mrna_dist_dict['stop_codon'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["stop_codon"][readlen]
						mrna_dist_dict['3_trailer'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["3_trailer"][readlen]
						mrna_dist_dict['total'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["total"][readlen]
						
					
				else:
					file_specific_mrna_dist_dict = {"5_leader":{},
													"start_codon":{},
													"cds":{},
													"stop_codon":{},
													"3_trailer":{},
													"total":{}}
					
					for readlen in range(minreadlen,maxreadlen+1):
							file_specific_mrna_dist_dict['5_leader'][readlen]=0
							file_specific_mrna_dist_dict['start_codon'][readlen]=0
							file_specific_mrna_dist_dict['cds'][readlen]=0
							file_specific_mrna_dist_dict['stop_codon'][readlen]=0
							file_specific_mrna_dist_dict['3_trailer'][readlen]=0
							file_specific_mrna_dist_dict['total'][readlen]=0

					for transcript in longest_tran_list:
						#if transcript != "YLR044C":
						#    continue
						try:
							transcript_dict = sqlite_db[transcript]["unambig"]
						except:
							continue 
						try:
							cds_start = int(traninfo_dict[transcript]["cds_start"])
							cds_stop = int(traninfo_dict[transcript]["cds_stop"])
						except:
							continue 
						
						

						for readlen in range(minreadlen,maxreadlen+1):
							if readlen not in mrna_dist_dict["5_leader"]:
								mrna_dist_dict['5_leader'][readlen]=0
								mrna_dist_dict['start_codon'][readlen]=0
								mrna_dist_dict['cds'][readlen]=0
								mrna_dist_dict['stop_codon'][readlen]=0
								mrna_dist_dict['3_trailer'][readlen]=0
								mrna_dist_dict['total'][readlen]=0
							if readlen not in transcript_dict:
								continue

							for five_pos in transcript_dict[readlen]:
								'''
								#print "five pos", five_pos
								# As in preiss paper we filter reads that are 20 fold higher than surrounding regions
								cur_count = transcript_dict[readlen][five_pos]
								#print "cur count", cur_count
								# Now we add up all positions in the -1 and +1 positions
								minus_count = 0
								plus_count = 0
								for secondary_readlen in range(readlen,readlen+1):
									if secondary_readlen in transcript_dict:
										if five_pos-1 in transcript_dict[secondary_readlen]:
											minus_count += transcript_dict[secondary_readlen][five_pos-1]
										if five_pos+1 in transcript_dict[secondary_readlen]:
											plus_count += transcript_dict[secondary_readlen][five_pos+1]
								#print "minus_count, plus_count", minus_count, plus_count
								if cur_count > (minus_count*20) or cur_count > (plus_count*20):
									#if transcript == "YLR044C":
									#    print "Skipping readlen {} and five_pos {}, cur_count {}, minus_count {}, plus_count {}".format(readlen,five_pos, cur_count, minus_count, plus_count)
									continue
								else:
									pass
									#if transcript == "YLR044C":
									#    print "ACCEPTING readlen {} and five_pos {}, cur_count {}, minus_count {}, plus_count {}".format(readlen,five_pos, cur_count, minus_count, plus_count)
								#print "five_pos", five_pos
								'''
								three_pos = five_pos+readlen
								#print "three_pos", three_pos
								#mrna_dist_dict["total"][readlen] += transcript_dict[readlen][five_pos]
								if three_pos <= cds_start+3:
									#print "adding to 5' count",readlen
									mrna_dist_dict["5_leader"][readlen] += transcript_dict[readlen][five_pos]
									file_specific_mrna_dist_dict["5_leader"][readlen] += transcript_dict[readlen][five_pos]
								elif five_pos <= cds_start-4 and three_pos >= cds_start+4:
									#print "adding to start codon count",readlen
									mrna_dist_dict["start_codon"][readlen] += transcript_dict[readlen][five_pos]
									file_specific_mrna_dist_dict["start_codon"][readlen] += transcript_dict[readlen][five_pos]
								elif five_pos >= cds_start-3 and three_pos <= cds_stop-2:
									#print "adding to cds count",readlen
									mrna_dist_dict["cds"][readlen] += transcript_dict[readlen][five_pos]
									file_specific_mrna_dist_dict["cds"][readlen] += transcript_dict[readlen][five_pos]
								elif five_pos <= cds_stop-9 and three_pos >= cds_stop-1:
									#print "adding to stop count",readlen
									mrna_dist_dict["stop_codon"][readlen] += transcript_dict[readlen][five_pos]
									file_specific_mrna_dist_dict["stop_codon"][readlen] += transcript_dict[readlen][five_pos]
								elif five_pos >= cds_stop-9:
									#print "adding to 3' trailer count",readlen
									mrna_dist_dict["3_trailer"][readlen] += transcript_dict[readlen][five_pos]
									file_specific_mrna_dist_dict["3_trailer"][readlen] += transcript_dict[readlen][five_pos]
									
									
					sqlite_db["mrna_dist_readlen_dict"] = file_specific_mrna_dist_dict
					sqlite_db.commit()
				sqlite_db.close()



		if smooth_amount != 0:
			smoothed_mrna_dist_dict = {}
			for read_type in ["5_leader","start_codon","cds","stop_codon","3_trailer"]:
				smoothed_mrna_dist_dict[read_type] = {}
				for i in range(minreadlen,maxreadlen+1):
					tot = 0
					for x in range(i-smooth_amount,i+smooth_amount):
						if x in mrna_dist_dict[read_type]:
							tot += (mrna_dist_dict[read_type][x])
					avg = tot/(smooth_amount*2)
					smoothed_mrna_dist_dict[read_type][i] = avg
			mrna_dist_dict = smoothed_mrna_dist_dict

		# If mrna_readlen_per is true normalize everything over the max in it's category
		if mrna_readlen_per == True:
			# For each category find the max value, max value will default to 1 if there is no values in that category.
			max_five = max(1,float(max(mrna_dist_dict["5_leader"].values())))
			max_start = max(1,float(max(mrna_dist_dict["start_codon"].values())))
			max_cds = max(1,float(max(mrna_dist_dict["cds"].values())))
			max_stop = max(1,float(max(mrna_dist_dict["stop_codon"].values())))
			max_three = max(1,float(max(mrna_dist_dict["3_trailer"].values())))
			for readlen in mrna_dist_dict["5_leader"]:
				mrna_dist_dict["5_leader"][readlen] = (float(mrna_dist_dict["5_leader"][readlen])/max_five)*100
			for readlen in mrna_dist_dict["start_codon"]:
				mrna_dist_dict["start_codon"][readlen] = (float(mrna_dist_dict["start_codon"][readlen])/max_start)*100
			for readlen in mrna_dist_dict["cds"]:
				mrna_dist_dict["cds"][readlen] = (float(mrna_dist_dict["cds"][readlen])/max_cds)*100
			for readlen in mrna_dist_dict["stop_codon"]:
				mrna_dist_dict["stop_codon"][readlen] = (float(mrna_dist_dict["stop_codon"][readlen])/max_stop)*100
			for readlen in mrna_dist_dict["3_trailer"]:
				mrna_dist_dict["3_trailer"][readlen] = (float(mrna_dist_dict["3_trailer"][readlen])/max_three)*100




		connection.close()
		return metainfo_plots.mrna_dist_readlen(mrna_dist_dict, mrna_readlen_per,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size)

	elif plottype == "rust_dwell":

		#traninfo_shelve = shelve.open("{0}{1}/{1}.shelf".format(trips_annotation_location,organism))
		#traninfo_dict = dict(traninfo_shelve)
		#traninfo_shelve.close()
		traninfo_dict = SqliteDict("{0}{1}/{1}.sqlite".format(trips_annotation_location,organism), autocommit=False)


		longest_tran_db = SqliteDict("/home/DATA/www/tripsviz/tripsviz/trips_annotations/homo_sapiens/principal_isoforms_5ldr3tlr_rnaseq.sqlite",autocommit=True)
		longest_tran_list = longest_tran_db["transcripts"]

		codon_count_dict = {"TTT":0, "TTC":0, "TTA":0, "TTG":0,
		"TCT":0, "TCC":0, "TCA":0, "TCG":0,
		"TAT":0, "TAC":0, "TAA":0, "TAG":0,
		"TGT":0, "TGC":0, "TGA":0, "TGG":0,
		"CTT":0, "CTC":0, "CTA":0, "CTG":0,
		"CCT":0, "CCC":0, "CCA":0, "CCG":0,
		"CAT":0, "CAC":0, "CAA":0, "CAG":0,
		"CGT":0, "CGC":0, "CGA":0, "CGG":0,
		"ATT":0, "ATC":0, "ATA":0, "ATG":0,
		"ACT":0, "ACC":0, "ACA":0, "ACG":0,
		"AAT":0, "AAC":0, "AAA":0, "AAG":0,
		"AGT":0, "AGC":0, "AGA":0, "AGG":0,
		"GTT":0, "GTC":0, "GTA":0, "GTG":0,
		"GCT":0, "GCC":0, "GCA":0, "GCG":0,
		"GAT":0, "GAC":0, "GAA":0, "GAG":0,
		"GGT":0, "GGC":0, "GGA":0, "GGG":0}


		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
					opendict = dict(sqlite_db)
					sqlite_db.close()
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)


				#TODO CHANGE THIS SO THAT THE CODON COUNT DICT IS OUTSIDE THE FILEPATH FOR LOOP
				offsets = opendict["offsets"]["fiveprime"]["offsets"]

				position_dict = {}

				for transcript in longest_tran_list:
					#print "transcript", transcript
					if transcript in opendict:
						#print "succcess"
						cds_start = traninfo_dict[transcript]["cds_start"]
						cds_stop = traninfo_dict[transcript]["cds_stop"]

						transeq = traninfo_dict[transcript]["seq"]

						#print "transeq", transeq
						for readlen in opendict[transcript]["unambig"]:
							#print "readlen", readlen
							try:
								offset = 15#temp setting offset as 15 to compare with gwips offsets[readlen]
							except:
								offset = 15
							for pos in opendict[transcript]["unambig"][readlen]:
								a_site = (pos+offset)+1
								#print "a_site", a_site
								if a_site > cds_start+120 and a_site < cds_stop-60:
									#print "success"
									codon = transeq[a_site:a_site+3]
									codon_count_dict[codon] += opendict[transcript]["unambig"][readlen][pos]
									'''
									for i in range(a_site-60,a_site+20,3):
										relative_pos = (a_site-i)
										codon = transeq[i:i+3]
										if relative_pos not in position_dict:
											position_dict[relative_pos] = {}
										if codon not in position_dict[relative_pos]:
											position_dict[relative_pos][codon] = 0
										position_dict[relative_pos][codon] += opendict[transcript]["unambig"][readlen][pos]
										codon_count_dict[codon] += opendict[transcript]["unambig"][readlen][pos]
										#print "added to pos dict"
									'''



		return metainfo_plots.rust_dwell(codon_count_dict,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size)

	elif plottype == "unmapped":
		master_dict = {}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath.split("/")[-1])
				if "frequent_unmapped_reads" not in sqlite_db:
					return "No unmapped reads data for {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath.split("/")[-1])

				#unmapped reads list is a list of tuples of length 100, first item in tuple is a sequence second is a count
				unmapped_reads_list = sqlite_db["frequent_unmapped_reads"]
				sqlite_db.close()

				for tup in unmapped_reads_list:
					if tup[0] in master_dict:
						master_dict[tup[0]] += tup[1]
					else:
						master_dict[tup[0]] = tup[1]

		title = "Most frequent unmapped reads ({})".format(short_code)
		studyname = ""

		top_reads = (sorted(master_dict.items(), key=operator.itemgetter(1)))[-50:]
		html_table = "<h1><center>{}</center></h1>".format(title)

		html_table += """<table class="unmapped_table">
		<thead><tr><th>Sequence</th><th>Frequency</th><th>Blast Link</th></tr></thead>"""
		for tup in top_reads[::-1]:
			html_table += ("<tr><td>{0}</td><td>{1}</td>    <td><a href='https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome&QUERY=%3E{2}_unmapped_sequence%0A{0}' target='_blank'>Blast</a></td></tr>".format(tup[0], tup[1],studyname))
		html_table += ("</table>")
		connection.close()
		return html_table
	elif plottype == "contamination":
		print "contaminant organism", contaminant_organism
		count_dict = {}
		master_sequence = ""
		contaminant_file = open("/home/DATA/www/tripsviz/tripsviz/static/contaminants/mycoplasma.fa")
		contaminant_lines = contaminant_file.read()
		contaminant_split = contaminant_lines.split(">")
		for entry in contaminant_split[1:]:
			header = entry.split("\n")[0]
			
			sequence = "".join(entry.split("\n")[1:])
			#sequences[header] = sequence
			master_sequence += sequence
	
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				filename = filepath.split("/")[-1]
				if filename not in count_dict:
					count_dict[filename] = {"count":0,"coverage":[],"unique_reads":0}
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath.split("/")[-1])
				if "frequent_unmapped_reads" not in sqlite_db:
					return "No unmapped reads data for {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath.split("/")[-1])

				#unmapped reads list is a list of tuples of length 100, first item in tuple is a sequence second is a count
				unmapped_reads_list = sqlite_db["frequent_unmapped_reads"]
				sqlite_db.close()

				for tup in unmapped_reads_list:
					read = tup[0]
					count = tup[1]
					readlen = len(read)
					#read = "AGAGCACGTTAAAGTGTGATGGCGTACATCTT"
					#print "read", read
					#for header in sequences:
					#	sequence = sequences[header]
					#position = master_sequence.find(read[1:8])
					#print "position", position
					#if position != -1:
					for x in range(0,len(master_sequence)-readlen):
						mismatches = 0
						for y in range(0,readlen):
							if master_sequence[x+y] != read[y]:
								mismatches += 1
							if mismatches > 2:
								break 
						if mismatches <= 2:
							count_dict[filename]["unique_reads"] += 1
							count_dict[filename]["count"] += count
							for i in range(x,x+readlen):
								if i not in count_dict[filename]["coverage"]:
									count_dict[filename]["coverage"].append(i)
		master_seq_len = len(master_sequence)
		for filename in count_dict:
			coverage = float(len(count_dict[filename]["coverage"]))/float(master_seq_len)
			coverage = round((coverage*100),2)
			count_dict[filename]["coverage"] = coverage
		
		
		title = "Contamination counts ({})".format(short_code)


		top_reads = (sorted(count_dict.items(), key=operator.itemgetter(1)))
		print "top reads", top_reads
		html_table = "<h1><center>{}</center></h1>".format(title)

		html_table += """<table class="unmapped_table">
		<thead><tr><th>Filename</th><th>Counts</th><th>Unique reads</th><th>Percentage coverage</th></tr></thead>"""
		for tup in top_reads[::-1]:
			html_table += ("<tr><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td></tr>".format(tup[0], tup[1]["count"],tup[1]["unique_reads"],tup[1]["coverage"]))
		html_table += ("</table>")
		connection.close()
		
		
		
		
		
		
		return html_table
		#reads = ["CGATGCTAGTCGTAGTCTAGTCGT","GTCGATCGTAGTCGTA","TAGCTAGCT","TAGCTAG","TACGCTAGCTAGCTAGCTG"]
		
		#for read in reads:
		#	m = sequence.find(read)
			
	
	

	elif plottype == "replicate_comp":
		labels = []
		transcript_dict = {}

		#This is used as an index to see which file number we're currently on
		file_count = 0

		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				cursor.execute("SELECT file_description from files where file_id = '{}';".format(file_id))
				result = cursor.fetchone();
				labels.append(result[0])
				filepath = file_paths_dict[filetype][file_id]

				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
					opendict =sqlite_db["unambiguous_all_totals"]
					sqlite_db.close()
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)

				for transcript in opendict.keys():
					if transcript not in transcript_dict:
						transcript_dict[transcript] = [0,0]

					if transcript in opendict:
						try:
							transcript_dict[transcript][file_count] = (log(opendict[transcript],2))
						except:
							transcript_dict[transcript][file_count] = 0
					else:
						transcript_dict[transcript][file_count] = 0


				file_count += 1
				if file_count >= 2:
					break

		minimum_reads = int(minimum_reads)
		if minimum_reads >0:
			min_log_val = log(minimum_reads,2)
		else:
			min_log_val = 0
		connection.close()
		if len(labels) == 1:
			return "Error: Select exactly two files for replicate comparison"
		return metainfo_plots.replicate_comp(labels, transcript_dict, min_log_val,short_code,background_col,str(title_size)+"pt", str(axis_label_size)+"pt", str(subheading_size)+"pt",str(marker_size)+"pt")

	elif plottype == "nuc_comp":
		print "NUC COMP FRAME IS", nuc_comp_frame
		
		#nuc_comp_minusoneframe
		#nuc_comp_cdsframe
		#nuc_comp_plusoneframe
		#nuc_comp_allframe


		
		master_count_dict = {"A":collections.OrderedDict(),"T":collections.OrderedDict(),"G":collections.OrderedDict(),"C":collections.OrderedDict()}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				if "nuc_counts" not in sqlite_db:
					return "No nucleotide counts data for this file, please report this to tripsvizsite@gmail.com or via the contact page."

				mappings = []
				if nuc_mapped == True:
					mappings.append("mapped")
				if nuc_unmapped == True:
					mappings.append("unmapped")

				if nuc_comp_direction == "nuc_comp_five":
					nuc_counts = sqlite_db["nuc_counts"]
					for mapping_type in mappings:
						for readlen in range(nuc_minreadlen, nuc_maxreadlen+1):
							for i in range(0,readlen+1):
								for nuc in ["A","T","G","C"]:
									if i not in master_count_dict[nuc]:
										master_count_dict[nuc][i] = 0
									if mapping_type in nuc_counts:
										if readlen in nuc_counts[mapping_type]:
											if i in nuc_counts[mapping_type][readlen]:
												master_count_dict[nuc][i] += nuc_counts[mapping_type][readlen][i][nuc]
									elif mapping_type == "mapped":
										if readlen in nuc_counts:
											if i in nuc_counts[readlen]:
												master_count_dict[nuc][i] += nuc_counts[readlen][i][nuc]
									else:
										return "No unmapped reads available for this file, please report this to tripsvizsite@gmail.com or via the contact page."

				elif nuc_comp_direction == "nuc_comp_three":
					nuc_counts = sqlite_db["threeprime_nuc_counts"]
					for mapping_type in mappings:
						for readlen in range(nuc_minreadlen, nuc_maxreadlen+1):
							for i in range(-1,-(nuc_maxreadlen),-1):
								for nuc in ["A","T","G","C"]:
									if i not in master_count_dict[nuc]:
										master_count_dict[nuc][i] = 0
									if mapping_type in nuc_counts:
										if readlen in nuc_counts[mapping_type]:
											if i in nuc_counts[mapping_type][readlen]:
												master_count_dict[nuc][i] += nuc_counts[mapping_type][readlen][i][nuc]
									elif mapping_type == "mapped":
										if readlen in nuc_counts:
											if i in nuc_counts[readlen]:
												master_count_dict[nuc][i] += nuc_counts[readlen][i][nuc]
									else:
										return "No unmapped reads available for this file, please report this to tripsvizsite@gmail.com or via the contact page."
				sqlite_db.commit()
				sqlite_db.close()


		#print "master count dict", master_count_dict

		master_dict = {"A":collections.OrderedDict(),
					   "T":collections.OrderedDict(),
					   "G":collections.OrderedDict(),
					   "C":collections.OrderedDict()}

		master_dict["A"][0] = 0
		master_dict["T"][0] = 0
		master_dict["G"][0] = 0
		master_dict["C"][0] = 0


		if nuc_comp_direction == "nuc_comp_five":
			for nuc in ["A","T","G","C"]:
				for i in range(0,nuc_maxreadlen):
					if i in master_count_dict[nuc]:
						thiscount = master_count_dict[nuc][i]
						othercount = 0.01
						for subnuc in ["A","T","G","C"]:
							othercount +=  master_count_dict[subnuc][i]
						if nuc_comp_type == "nuc_comp_per":
							master_dict[nuc][i] = float(thiscount)/float(othercount)
						elif nuc_comp_type == "nuc_comp_count":
							master_dict[nuc][i] = float(thiscount)
		elif nuc_comp_direction == "nuc_comp_three":
			for nuc in ["A","T","G","C"]:
				for i in range(-1,-(nuc_maxreadlen),-1):
					if i in master_count_dict[nuc]:
						thiscount = master_count_dict[nuc][i]
						othercount = 0.01
						for subnuc in ["A","T","G","C"]:
							othercount +=  master_count_dict[subnuc][i]
						if nuc_comp_type == "nuc_comp_per":
							master_dict[nuc][i] = float(thiscount)/float(othercount)
						elif nuc_comp_type == "nuc_comp_count":
							master_dict[nuc][i] = float(thiscount)

		title = "Nucleotide composition"

		connection.close()
		return metainfo_plots.nuc_comp(master_dict, nuc_maxreadlen,title, nuc_comp_type,nuc_comp_direction,short_code,background_col,a_col,t_col,g_col,c_col,title_size, axis_label_size, subheading_size,marker_size)



	elif plottype == "dinuc_bias":
		master_count_dict = collections.OrderedDict([("AA",0), ("AT",0), ("AG",0), ("AC",0),
													 ("TA",0), ("TT",0), ("TG",0), ("TC",0),
													 ("GA",0), ("GT",0), ("GG",0), ("GC",0),
													 ("CA",0), ("CT",0), ("CG",0), ("CC",0)])

		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				dinuc_counts = sqlite_db["dinuc_counts"]
				for readlen in dinuc_counts:
					for dinuc in dinuc_counts[readlen]:
						master_count_dict[dinuc] += dinuc_counts[readlen][dinuc]

		connection.close()
		return metainfo_plots.dinuc_bias(master_count_dict,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size)


	elif plottype == "fastq_screen":
		#print "Plot type is fastqscrren at {}".format(time.time()-start_time)
		html_filepath = ""
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if html_filepath == "":
					html_filepath = filepath.replace(".sqlite","_lessrRNA_screen.html")
				else:
					return "Error: Only one dataset at a time can be selected for fastq screen"
		#print "Have html filepath at {}".format(time.time()-start_time)
		if os.path.isfile(html_filepath):
			openfile = open(html_filepath,"r")
			#print "File opened at {}".format(time.time()-start_time)
			fastq_lines = openfile.readlines()

			#print "Fastq html at {}".format(time.time()-start_time)
			
			#fixed_html = str(fastq_html.replace("padding:0 20px 20px","").replace("max-width:1200px;","padding:0 20px 20px").replace("<html>","").replace("</html>","").replace("<body>","").replace("</body>","")).replace("<!DOCTYPE html>","").replace("<head>","").replace("</head>","").replace("container","container2")
			#print "Fixed html called at {}".format(time.time()-start_time)
			#no_local_plotly =fastq_lines[:150]+fastq_lines[-209:]
			#The base64 encoded png string in the header is too long for firefox, will work for one plot and then crash firefox, this is to prevent that
			fixed_html = ""
			for line in fastq_lines:
				if "iVBORw0KGgoAAAANSUhEUgAAA4wAAAGVCAYAAAHC" in line:
					fixed_html += '<a style="float:left;" href="http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen" target="_blank"><img width="50%" height="50%" alt="FastQ Screen"src="/static/fastq_screen.png"</a>'
					continue
				else:
					fixed_html += (line)
			#print "fixed_html",fixed_html
			#<a style="float:left;" href="http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen" target="_blank"><img alt="FastQ Screen" src=
			
			# Replace serves two functions here, first remove the padding line in the body tag as this affects the header bar and everything else on trips,
			# but removal means fastq screen logo is slightly off screen
			# Second remove the max-width line in the .container class, replace it with the padding line removed from the body tag, as this will now be specific to the container
			# and fix the fastq screen logo.
			fixed_html = str(fixed_html.replace("padding:0 20px 20px","").replace("max-width:1200px;","padding:0 20px 20px").replace("<html>","").replace("</html>","").replace("<body>","").replace("</body>","")).replace("<!DOCTYPE html>","").replace("<head>","").replace("</head>","").replace("container","container2")
			return fixed_html
		else:
			return "No fastq_screen file available for this dataset"
	elif plottype == "explore_offsets":
		readlen_dict = {}
		#traninfo_shelve = shelve.open("{0}{1}/{1}.shelf".format(trips_annotation_location,organism))
		#traninfo_dict = dict(traninfo_shelve)
		#traninfo_shelve.close()
		traninfo_dict = SqliteDict("{0}{1}/{1}.sqlite".format(trips_annotation_location,organism), autocommit=False)
		tranlist = traninfo_dict.keys()[:10000]

		labels = []
		f0_counts = []
		f1_counts = []
		f2_counts = []

		#For the first file in selected files
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
					opendict = dict(sqlite_db)
					sqlite_db.close()
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page".format(filepath)
				tranlist = opendict.keys()

				#For each readlength we display on the final graph
				for readlen in range(25,36):
					try:
						chosen_offset = opendict["offsets"]["fiveprime"]["offsets"][readlen]
					except:
						chosen_offset = 15
					labels.append("")
					labels.append("{}_{}".format(readlen, chosen_offset))
					labels.append("")
					for offset in [chosen_offset-1,chosen_offset,chosen_offset+1]:
						trancount = 0
						inframe_counts = 0
						outframe_counts = 0


						#print "offset", offset
						if readlen not in readlen_dict:
							readlen_dict[readlen] = {chosen_offset-1:[0,0,0],chosen_offset:[0,0,0],chosen_offset+1:[0,0,0]}
							#print "readlen, readlendict[readlen]", readlen, readlen_dict[readlen]


						#frame_counts keeps track of the total counts in each frame for a particular readlength and offset
						#frame_counts = {0:0,1:0,2:0}

						#for each transcript get the frame counts breakdown from the cds given this particular offset
						for tran in tranlist:
							if tran not in traninfo_dict:
								continue
							tempdict = dict(opendict[tran])
							trancount += 1
							if trancount > 5000:
								break

							if "cds_start" not in traninfo_dict[tran]:
								continue
							cds_start = traninfo_dict[tran]["cds_start"]
							cds_stop = traninfo_dict[tran]["cds_stop"]

							if cds_start == "NULL" or cds_stop == "NULL":
								continue
							if cds_start <= 1 or cds_stop <= 1:
								continue
							#to account for 0-based counts ,without this line the frame will be wrong
							cds_start += 1
							cds_frame = cds_start%3
							#first walk through this entry in the opendict for only the readlength in question applying the relevant offset
							count_dict = {}

							if readlen in tempdict["unambig"]:
								#print "readlen", readlen
								for fiveprime_pos in tempdict["unambig"][readlen]:
									count = tempdict["unambig"][readlen][fiveprime_pos]
									new_pos = fiveprime_pos + offset
									count_dict[new_pos] = count

								for i in range(cds_start,cds_stop):
									frame = i%3
									if i in count_dict:
										if frame == cds_frame:
											inframe_counts += count_dict[i]
										else:
											outframe_counts += count_dict[i]

						f0_counts.append(inframe_counts)
						f1_counts.append(outframe_counts)
						f2_counts.append(0)
		connection.close()
		return metainfo_plots.explore_offsets(f0_counts, f1_counts, f2_counts, labels,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size)







	elif plottype == "metagene_plot":
		minpos = -300
		maxpos = 300
		pos_list = []
		for i in range(minpos, maxpos+1):
			pos_list.append(i)
		count_dict = {"fiveprime":{}, "threeprime":{}}
		for primetype in ["fiveprime","threeprime"]:
			for i in range(minpos,maxpos+1):
				count_dict[primetype][i] = 0

		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				if "metagene_counts" not in sqlite_db:
					return "No metagene counts data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
				if metagene_type == "metagene_start":
					mgc = sqlite_db["metagene_counts"]
				elif metagene_type == "metagene_stop":
					mgc = sqlite_db["stop_metagene_counts"]
				elif metagene_type == "metagene_second_aug":
					mgc = sqlite_db["secondary_metagene_counts"]

				sqlite_db.close()
				#print "\n\n\n\n\n\n\n\nMGC", mgc
				for primetype in ["fiveprime", "threeprime"]:
					#print "primetype", primetype
					for i in pos_list:
						#print "checking i", i
						for readlen in range(minreadlen, maxreadlen+1):
							#print "readlen", readlen
							if readlen in mgc[primetype]:
								if i in mgc[primetype][readlen]:
									#print i, " is in mgc", mgc[primetype][readlen]
									count_dict[primetype][i] += (mgc[primetype][readlen][i])

		fiveprime_counts = []
		threeprime_counts = []
		for i in pos_list:
			fiveprime_counts.append(count_dict["fiveprime"][i])
			threeprime_counts.append(count_dict["threeprime"][i])


		title = "Metagene profile"
		#print pos_dist, count_dict["fiveprime"],count_dict["threeprime"]
		#print fiveprime_counts, threeprime_counts
		#print "fetching metagene plot"
		connection.close()
		return metainfo_plots.metagene_plot(pos_list,fiveprime_counts,threeprime_counts,metagene_type,title,minpos, maxpos,short_code,background_col,metagene_fiveprime_col,metagene_threeprime_col,title_size, axis_label_size, subheading_size,marker_size,metagene_five, metagene_three)


	elif plottype == "trip_periodicity":
		read_dict = {"readlengths":[],
					 "frame1":[],
					 "frame2":[],
					 "frame3":[]}
		if trip_maxreadlen < trip_minreadlen:
			return "Error: max read length less than min read length, increase max read length using the input at the top of the page."
		for i in range(trip_minreadlen, trip_maxreadlen+1):
			read_dict["readlengths"].append(i)
			read_dict["frame1"].append(0)
			read_dict["frame2"].append(0)
			read_dict["frame3"].append(0)

		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				if "trip_periodicity" not in sqlite_db:
					return "No triplet periodicity data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
				trip_periodicity_dict = sqlite_db["trip_periodicity"]
				sqlite_db.close()
				readlen_index = 0
				for readlength in read_dict["readlengths"]:
					if readlength in trip_periodicity_dict["fiveprime"]:
						read_dict["frame1"][readlen_index] += trip_periodicity_dict["fiveprime"][readlength]["0"]
						read_dict["frame2"][readlen_index] += trip_periodicity_dict["fiveprime"][readlength]["1"]
						read_dict["frame3"][readlen_index] += trip_periodicity_dict["fiveprime"][readlength]["2"]
					readlen_index += 1

		title = "Triplet periodicity"

		connection.close()
		return metainfo_plots.trip_periodicity_plot(read_dict,title,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size)

	elif plottype == "mapped_reads_plot":
		labels = [""]
		unmapped = [0]
		mapped_coding = [0]
		mapped_noncoding = [0]
		ambiguous = [0]
		cutadapt_removed = [0]
		rrna_removed = [0]

		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]

				cursor.execute("SELECT file_name,file_description from files where file_id = '{}';".format(file_id))
				result = cursor.fetchone();
				file_name = (result[0]).replace(".shelf","")
				labels.append(result[1])

				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)


				if "unmapped_reads" in sqlite_db:
					unmapped.append(sqlite_db["unmapped_reads"])
				else:
					unmapped.append(0)

				if "coding_counts" in sqlite_db:
					mapped_coding.append(sqlite_db["coding_counts"])
				else:
					mapped_coding.append(0)

				if "noncoding_counts" in sqlite_db:
					mapped_noncoding.append(sqlite_db["noncoding_counts"])
				else:
					mapped_noncoding.append(0)

				if "ambiguous_counts" in sqlite_db:
					ambiguous.append(sqlite_db["ambiguous_counts"])
				else:
					ambiguous.append(0)

				if "cutadapt_removed" in sqlite_db:
					cutadapt_removed.append(sqlite_db["cutadapt_removed"])
				else:
					cutadapt_removed.append(0)

				if "rrna_removed" in sqlite_db:
					rrna_removed.append(sqlite_db["rrna_removed"])
				else:
					rrna_removed.append(0)


				sqlite_db.close()


		connection.close()
		labels.append("")


	

		'''
		if breakdown_per == True:
			for i in range(1,len(labels)):
				total = float(unmapped[i]+mapped_coding[i]+mapped_noncoding[i]+ambiguous[i]+cutadapt_removed[i]+rrna_removed[i])
				unmapped[i] = unmapped[i]/total
				mapped_coding[i] = mapped_coding[i]/total
				mapped_noncoding
		'''	
				
				
				
		for listname in [unmapped, mapped_coding, mapped_noncoding, ambiguous, cutadapt_removed, rrna_removed]:
			listname.append(0)
		return metainfo_plots.mapped_reads_plot(unmapped, mapped_coding, mapped_noncoding, labels,ambiguous,cutadapt_removed,rrna_removed,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size,breakdown_per)




	elif plottype == "heatmap":
		min_readlen = heatmap_minreadlen
		max_readlen = heatmap_maxreadlen
		min_pos = heatmap_startpos
		max_pos = heatmap_endpos
		count_list = []
		positions = []
		readlengths = []



		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				if "metagene_counts" not in sqlite_db:
					return "No metagene counts data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
				if heatmap_metagene_type == "metagene_start":
					mgc = sqlite_db["metagene_counts"]
				elif heatmap_metagene_type == "metagene_stop":
					mgc = sqlite_db["stop_metagene_counts"]
				elif heatmap_metagene_type == "metagene_second_aug":
					mgc = sqlite_db["secondary_metagene_counts"]
				sqlite_db.close()



				for primetype in [heatmap_direction]:
					for readlen in range(max_readlen, min_readlen-1, -1):
						for i in range(min_pos, max_pos+1):
							if readlen in mgc[primetype]:
								if i in mgc[primetype][readlen]:
									if mgc[primetype][readlen][i] != 0:
										if log_scale == True:
											count_list.append(log(mgc[primetype][readlen][i],2))
										else:
											count_list.append(mgc[primetype][readlen][i])
									else:
										count_list.append(None)
								else:
									count_list.append(None)
								readlengths.append(readlen)
								positions.append(i)
							else:
								readlengths.append(readlen)
								positions.append(i)
								count_list.append(None)

		title = "Heatmap"


		connection.close()
		return metainfo_plots.heatplot(min_readlen, max_readlen, min_pos, max_pos, positions, readlengths,count_list,heatmap_metagene_type,title,reverse_scale,color_palette,short_code,background_col,maxscaleval,str(title_size)+"pt", str(axis_label_size)+"pt", str(subheading_size)+"pt",str(marker_size)+"pt")






	else:
		if plottype not in ["replicate_comp"]:
			print "ERROR2 plottype is not in list"
		if (plottype.strip(" ").replace("\n")) not in ["replicate_comp"]:
			print "ERROR 3"

	return "Error, unknown plot type selected: {}".format(plottype)




if __name__ == '__main__':
	local=False
	try:
		if sys.argv[1] == "true":
			local = True
	except:
		pass
	if local == False:
		app.run(host='0.0.0.0',debug=False)
	else:
		app.run(host='0.0.0.0',debug=True)
