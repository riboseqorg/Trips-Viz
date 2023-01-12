from flask import Blueprint, render_template, abort, request
import sqlite3
from sqlitedict import SqliteDict
import ast
import os
import time
import logging
import config
from core_functions import fetch_studies, fetch_files,fetch_study_info,fetch_file_paths,generate_short_code,build_profile,build_proteomics_profile,nuc_to_aa,fetch_user
import collections
from flask_login import current_user
from bokeh.plotting import figure, show, output_file
from bokeh.embed import file_html, components
from bokeh.resources import CDN
from bokeh.palettes import YlOrRd9 as palette
from bokeh.palettes import inferno
from bokeh.palettes import all_palettes
from bokeh.models.glyphs import Text
import bokeh.models as bmo
from bokeh.io import show
import random
from scipy.stats import mannwhitneyu,wilcoxon,zscore
from statsmodels.stats.multitest import multipletests
from bokeh.models import (
	TapTool,
	OpenURL,
	Range1d,
	Label,
	FuncTickFormatter,
	LogTicker,
	ColumnDataSource,
	HoverTool,
	LinearColorMapper,
	LogColorMapper,
	BasicTicker,
	PrintfTickFormatter,
	ColorBar
)


from tensorflow.keras.models import Sequential

from tensorflow.keras.layers import Dense, Activation, Dropout, Embedding, LSTM,BatchNormalization
import numpy as np
from tensorflow.keras import optimizers
from tensorflow.keras import backend as K
from tensorflow.keras.optimizers import RMSprop
import email
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText



#This page is used to detect translated open reading frames
translated_orf_blueprint = Blueprint("orf_translationpage", __name__, template_folder="templates")
@translated_orf_blueprint.route('/<organism>/<transcriptome>/orf_translation/')
def orf_translationpage(organism,transcriptome):
	#ip = request.environ['REMOTE_ADDR']
	global local
	try:
		logging.debug(local)
	except:
		local = False

	organism = str(organism)
	user,logged_in = fetch_user()		
	accepted_studies = fetch_studies(user, organism,transcriptome)
	file_id_to_name_dict, accepted_studies,accepted_files,seq_types = fetch_files(accepted_studies)
	advanced = False
			
	# holds all values the user could possibly pass in the url (keywords are after request.args.get), anything not passed by user will be a string: "None"
	html_args = {"user_short":str(request.args.get('short')),
				 "start_codons":str(request.args.get('start_codons')),
				 "min_start_inc":str(request.args.get('min_start_inc')),
				 "max_start_inc":str(request.args.get('max_start_inc')),
				 "min_stop_dec":str(request.args.get('min_stop_dec')),
				 "max_stop_dec":str(request.args.get('max_stop_dec')),
				 #"min_cds_rat":str(request.args.get('min_cds_rat')),
				 #"max_cds_rat":str(request.args.get('max_cds_rat')),
				 "min_lfd":str(request.args.get('min_lfd')),
				 "max_lfd":str(request.args.get('max_lfd')),
				 "min_hfd":str(request.args.get('min_hfd')),
				 "max_hfd":str(request.args.get('max_hfd')),
				 "min_cds":str(request.args.get('min_cds')),
				 "max_cds":str(request.args.get('max_cds')),
				 "min_len":str(request.args.get('min_len')),
				 "max_len":str(request.args.get('max_len')),
				 "min_avg":str(request.args.get('min_avg')),
				 "max_avg":str(request.args.get('max_avg')),
				 "tran_list":str(request.args.get('tran_list')),
				 "ambig":str(request.args.get('ambig')),
				 "saved_check":str(request.args.get('saved_check')),
				 "sic":str(request.args.get('sic')),
				 "sdc":str(request.args.get('sdc')),
				 "crc":str(request.args.get('crc')),
				 "lfdc":str(request.args.get('lfdc')),
				 "hfdc":str(request.args.get('hfdc')),
				 "ambig":str(request.args.get('ambig')),
				 "saved_check":str(request.args.get('saved_check')),
				 }
	
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

	if html_args["start_codons"] != None:
		start_codons = html_args["start_codons"].split(",")
		html_args["start_codons"] = ["#"+str(x).strip(" ") for x in start_codons]
	else:
		html_args["start_codons"] = []
	
	#If ever need to use other seq types than these modify fetch_files to return a proper list
	seq_types = ["riboseq","proteomics"]
	# If seq type not riboseq remove it from accepted files as orf translation is only appicable to riboseq
	del_types = []
	for seq_type in accepted_files:
		if seq_type not in seq_types:
			del_types.append(seq_type)
	for seq_type in del_types:
			del accepted_files[seq_type]
	studyinfo_dict = fetch_study_info(organism)
	return render_template('orf_translation.html', studies_dict=accepted_studies, accepted_files=accepted_files,user=user,
						   organism=organism,default_tran="",local=local,transcriptome=transcriptome,advanced=advanced,
						   seq_types=seq_types,studyinfo_dict=studyinfo_dict,html_args=html_args)







def create_orf_plot(sorted_all_values,organism,transcriptome,file_string,filename):
	file_string = file_string.replace(";",",")
	master_dict = {}
	total_tups = 0
	for tup in sorted_all_values:
		if total_tups > 2000:
			break
		gene = tup[0]
		transcript = tup[1]
		start = tup[2]
		stop = tup[3]
		length = tup[4]
		rank = tup[15]
		orftype = tup[16]
		if orftype not in master_dict:
			master_dict[orftype] = {"x":[],"y":[],"trans":[],"genes":[],"starts":[],"stops":[],"orftype":[]}
		master_dict[orftype]["x"].append(rank)
		master_dict[orftype]["y"].append(length)
		master_dict[orftype]["trans"].append(transcript)
		master_dict[orftype]["genes"].append(gene)
		master_dict[orftype]["starts"].append(start)
		master_dict[orftype]["stops"].append(stop)
		master_dict[orftype]["orftype"].append(orftype)
		#trips_link = '<a href="https://trips.ucc.ie/'+organism+'/'+transcriptome+'/interactive_plot/?tran='+transcript+'&hili='+str(start)+"_"+str(stop)+'&files='+file_string+'" target="_blank_" >View on trips-viz</a>'
		total_tups += 1
	p = figure(plot_width=2200, plot_height=1800,x_axis_label="",  y_axis_label='GC%',title="GC%",toolbar_location="below",
			tools = "reset,pan,box_zoom,save,hover,tap")
	
	color_list = ["red","green","blue","orange","purple"]
	color_index = 0
	for orftype in master_dict:
		#logging.debug("orftype, color_index",orftype, color_index)
		source = ColumnDataSource({'Rank': master_dict[orftype]["x"],'Length':master_dict[orftype]["y"],'transcript':master_dict[orftype]["trans"],'gene':master_dict[orftype]["genes"],'start':master_dict[orftype]["starts"],'stop':master_dict[orftype]["stops"],'orftype':master_dict[orftype]["orftype"]})	
		p.scatter('Rank','Length', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color=color_list[color_index],legend=orftype)
		color_index += 1

	p.legend.location = "top_left"
	p.legend.click_policy="hide"
	p.legend.label_text_font_size = "28px"
	p.title.align="center"
	#p.title.text_font_size = title_size
	#p.xaxis.axis_label_text_font_size = axis_label_size
	#p.xaxis.major_label_text_font_size = marker_size
	#p.yaxis.axis_label_text_font_size = axis_label_size
	#p.yaxis.major_label_text_font_size = marker_size
	#p.background_fill_color = background_col
	p.xgrid.grid_line_color = "#cccccc"
	p.ygrid.grid_line_color = "#cccccc"
	

	hover = p.select(dict(type=HoverTool))
	hover.mode = 'mouse'
	hover.tooltips = [("Length:", "@Length"),("Rank:","@Rank"),("Transcript:","@transcript"),("Gene:","@gene"),("Start:","@start"),("Stop:","@stop"),("Orftype:","@orftype")]
	
	#url = "http://trips.ucc.ie/{}/{}/comparison/?files={}{}&transcript=@trans&normalize={}&cov=T&ambig={}&minread=25&maxread=150".format(organism, transcriptome,file_string,label_string,str(normalized)[0],ambig)
	url ="https://trips.ucc.ie/{}/{}/interactive_plot/?tran=@transcript&hili=@start\_@stop&files={}".format(organism, transcriptome, file_string[:-1])
	#https://trips.ucc.ie/homo_sapiens/Gencode_v25/interactive_plot/?tran=ENST00000376263&hili=???@stop&files=710;
	#https://trips.ucc.ie/homo_sapiens/Gencode_v25/interactive_plot/?tran=ENST00000376263&hili=130_198&files=710,
	taptool = p.select(type=TapTool)
	taptool.callback = OpenURL(url=url)
	
	graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as csv file</b></button></a> </div>".format(filename)
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br> </div>".format(short_code)
	#layout = column(text_input, p)
	graph += file_html(p,CDN)
	return graph






def tran_to_genome(tran, pos, transcriptome_info_dict):
	if tran in transcriptome_info_dict:
		traninfo = transcriptome_info_dict[tran]
	else:
		return ("Null",0)
	chrom = traninfo["chrom"]
	strand = traninfo["strand"]
	exons = traninfo["exons"]
	#logging.debug(exons)
	if strand == "+":
		exon_start = 0
		for tup in exons:
			exon_start = tup[0]
			exonlen = tup[1] - tup[0]
			if pos > exonlen:
				pos = (pos - exonlen)-1
			else:
				break
		genomic_pos = (exon_start+pos)-1
	elif strand == "-":
		exon_start = 0
		for tup in exons[::-1]:
			exon_start = tup[1]
			exonlen = tup[1] - tup[0]
			if pos > exonlen:
				pos = (pos - exonlen)-1
			else:
				break
		genomic_pos = (exon_start-pos)+1
	return "{}_{}".format(chrom, genomic_pos)



def get_highly_expressed_trans(cds_dict, file_paths_dict,no_cases):
	for file_id in file_paths_dict["riboseq"]:
		sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
		#print file_paths_dict["riboseq"][file_id]
		
		if "aggregate_riboseq" not in file_paths_dict["riboseq"][file_id]:
			counts = sqlite_db["unambiguous_cds_totals"]
			#print "aggregate_riboseq not in", file_paths_dict["riboseq"][file_id]
			for transcript in counts:
				if transcript in cds_dict:
					count = counts[transcript]
					cds_dict[transcript]["count"] += count
		else:
			counts = {}
			for key in sqlite_db:
				counts[key] = sqlite_db[key]
			sqlite_db.close()
			for transcript in counts:
				
				if transcript in cds_dict:
					cds_start = cds_dict[transcript]["cds_start"]
					cds_stop = cds_dict[transcript]["cds_stop"]
					count = 0 
					for pos in counts[transcript]["riboseq"]:
						if pos > cds_start and pos < cds_stop:
							count += counts[transcript]["riboseq"][pos]
					cds_dict[transcript]["count"] += count
			
	all_avg_counts = []
	for transcript in cds_dict:
		cds_start = cds_dict[transcript]["cds_start"]
		cds_stop = cds_dict[transcript]["cds_stop"]
		length = cds_dict[transcript]["length"]
		if cds_start > 10 and cds_stop < (length-10):
			cds_len = cds_dict[transcript]["cds_stop"] - cds_dict[transcript]["cds_start"]
			avg_count = cds_dict[transcript]["count"]/cds_len
			all_avg_counts.append((avg_count,transcript))
	top_transcripts = []

	for tup in sorted(all_avg_counts)[(no_cases*-1):]:
		top_transcripts.append(tup[1])
	#print "top transcripts", top_transcripts
	return top_transcripts
		


def create_aggregate(file_paths_dict,study_path, seq_type):
	logging.debug("aggregate filepathsdict", file_paths_dict)
	ambig = "unambig"
	file_count = 0
	offset_dict = {}
	profile_dict = {}
	file_list = []
	for file_id in file_paths_dict[seq_type]:
		file_list.append(file_id)
		file_count += 1
		sqlite_dict = SqliteDict(file_paths_dict[seq_type][file_id])
		sqlite_db = dict(sqlite_dict)
		sqlite_dict.close()
		#offsets = offset_dict[file_id]
		offsets = {}
		if seq_type == "riboseq":
			try:
				offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
			except:
				pass
		for transcript in sqlite_db:
			if transcript not in profile_dict:
				profile_dict[transcript] = {"riboseq":{},"proteomics":{}}
			try:
				counts = sqlite_db[transcript]
			except:
				continue
			if seq_type == "riboseq":
				subprofile = build_profile(counts, offsets,ambig)
			elif seq_type == "proteomics":
				subprofile = build_proteomics_profile(counts,ambig)
			for pos in subprofile:
				try:
					profile_dict[transcript][seq_type][pos] += subprofile[pos]
				except:
					profile_dict[transcript][seq_type][pos] = subprofile[pos]
		logging.debug("{} files read".format(file_count))
	outfile = SqliteDict("{}/aggregate_{}.sqlite".format(study_path, seq_type))	
	for tran in profile_dict:
		outfile[tran] = profile_dict[tran]
	outfile["file_list"] = file_list
	outfile.commit()
	outfile.close()
	


def create_profiles(file_paths_dict,accepted_transcript_list,ambig,total_files,minscore):
	#logging.debug(file_paths_dict)
	#If
	file_count = 0
	# This will be populated with the users chosen file_ids and passed to the table, so that the trips link can use these files aswell.
	file_string = ""
	ribo_study_string = "&ribo_studies="
	proteomics_study_string = "&proteomics_studies="
	offset_dict = {}
	score_dict = {}
	seq_types = ["riboseq","proteomics"]
	#for seq_type in seq_types:
	#	if seq_type in file_paths_dict:
	#		for file_id in file_paths_dict[seq_type]:
	#			file_string += "{};".format(file_id)
	#			sqlite_db = SqliteDict(file_paths_dict[seq_type][file_id])
	#			try:
	#				offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
	#				offset_dict[file_id] = offsets
	#			except:
	#				offset_dict[file_id] = {}
	#			try:
	#				scores  = sqlite_db["offsets"]["fiveprime"]["read_scores"]
	#				score_dict[file_id] = scores
	#			except:
	#				score_dict[file_id] = {}
	#			sqlite_db.close()
	profile_dict = {}
	for seq_type in seq_types:
		if seq_type not in file_paths_dict:
			continue
		for file_id in file_paths_dict[seq_type]:
			#logging.debug("file_id", file_id)
			file_count += 1
			sqlite_db = SqliteDict(file_paths_dict[seq_type][file_id])
			if type(file_id) == str:
				if "STUDY" in file_id:
					subfiles = sqlite_db["file_list"]
					if seq_type == "riboseq":
						ribo_study_string += "{};".format(file_id.replace("STUDY_",""))
					else:
						proteomics_study_string += "{};".format(file_id.replace("STUDY_",""))
					#for sub_file_id in subfiles:
					#	file_string += "{};".format(sub_file_id)
					#logging.debug("success")
					for transcript in accepted_transcript_list:
						if transcript not in profile_dict:
							profile_dict[transcript] = {"riboseq":{},"proteomics":{}}
						if transcript in sqlite_db:
							#logging.debug("tran in sqlite")
							subprofile = sqlite_db[transcript][seq_type]
							#logging.debug("subprofile", subprofile)
							for pos in subprofile:
								try:
									profile_dict[transcript][seq_type][pos] += subprofile[pos]
								except:
									profile_dict[transcript][seq_type][pos] = subprofile[pos]
			else:
				file_string += "{};".format(file_id)
				offsets = {}
				scores = {}
				if seq_type == "riboseq":
					try:
						offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
					except:
						pass
					try:
						scores = sqlite_db["offsets"]["fiveprime"]["read_scores"]
					except:
						pass

				
				for transcript in accepted_transcript_list:
					if transcript not in profile_dict:
						profile_dict[transcript] = {"riboseq":{},"proteomics":{}}
					try:
						counts = sqlite_db[transcript]
					except:
						continue
					if seq_type == "riboseq":
						subprofile = build_profile(counts, offsets,ambig,minscore,scores)
					elif seq_type == "proteomics":
						subprofile = build_proteomics_profile(counts,ambig)
					for pos in subprofile:
						try:
							profile_dict[transcript][seq_type][pos] += subprofile[pos]
						except:
							profile_dict[transcript][seq_type][pos] = subprofile[pos]
			#logging.debug("{}/{} files read".format(file_count, total_files))
			sqlite_db.close()
	#logging.debug("profile_dict", profile_dict["ENST00000233893"])
	file_string+=ribo_study_string
	file_string +=proteomics_study_string
	return (profile_dict,file_string)
	
	
	
def extract_cds_counts(profile_dict, cds_dict):
	cds_count_dict = {}
	for tran in cds_dict:
		if tran in profile_dict:
			ribo_counts =  profile_dict[tran]["riboseq"]
			cds_start = cds_dict[tran]["cds_start"]
			cds_stop = cds_dict[tran]["cds_stop"]
			cds_count = 0
			for i in range(cds_start+1,cds_stop,3):
				if i in ribo_counts:
					cds_count += ribo_counts[i]
			cds_count_dict[tran] = cds_count/float(cds_stop-cds_start)
	return cds_count_dict
	
	
	
def geo_mean(iterable):
    a = np.array(iterable)
    return a.prod()**(1.0/len(a))
	
	
	
def extract_nnet_values(accepted_orf_dict,data,tran_gene_dict,selected_seq_types,profile_dict,cds_count_dict):
	logging.debug("extract nnet values called")
	best_high_frame = 0
	best_low_frame = 0
	best_start_score = 0
	best_stop_score = 0
	best_inframe_cov = 0
	all_values = []
	tot_loc = 0
	for locus in accepted_orf_dict:
		tot_loc += 1
		if tot_loc%100 == 0:
			logging.debug("total transcripts {}".format(tot_loc))
		for stop in accepted_orf_dict[locus]:
			best_values = {"start":-1,"high_frame_count":1,"low_frame_count":1,"start_score":-1,"stop_score":1,"final_score":-10000,"coverage":0,"length":0,
				  "transcript":locus,"stop":0,"proteomics_count":0,"ratio":0,"inframe_count":0,"start_ratio":0,"stop_ratio":0,"high_ratio":0,"low_ratio":0}
			for start in accepted_orf_dict[locus][stop]:
				#print "locus start stop", locus, start, stop
				transcript = accepted_orf_dict[locus][stop][start]["transcript"]
				best_values["transcript"] = transcript
				if best_values["start"] == -1:
					best_values["start"] = start
					best_values["stop"] = stop
				gene = tran_gene_dict[transcript]
				if "riboseq" in selected_seq_types and "proteomics" in selected_seq_types:
					profile = profile_dict[transcript]["riboseq"]
					prot_profile = profile_dict[transcript]["proteomics"]
				elif "proteomics" not in selected_seq_types:
					profile = profile_dict[transcript]["riboseq"]
					prot_profile = {}
				else:
					profile = profile_dict[transcript]["proteomics"]
					prot_profile = {}
				transcriptome_stop = accepted_orf_dict[locus][stop][start]["stop"]
				orftype = accepted_orf_dict[locus][stop][start]["orftype"]
				length = (transcriptome_stop-start)+1
				inframe_count = 0
				proteomics_count = 0
				if_cl = []
				mo_count = 0
				po_count = 0
				if_cov = 0.0
				if_len = (stop-start)/3
				if if_len <= 0:
					print "error, transcript, start, stop, orftype", transcript, start, stop, orftype
					return ""
				
				for p in range(start+1, transcriptome_stop,3):
					if p in prot_profile:
						proteomics_count += prot_profile[p]
				#BREAKPOINT
				'''
				for x in range(start+6, transcriptome_stop-9,3):
					#logging.debug("x", x)
					curr_mo = 0
					curr_po = 0
					if_len += 1
					if x-1 in profile:
						mo_count += profile[x-1]
						curr_mo = profile[x-1]
					if x+1 in profile:
						po_count += profile[x+1]
						curr_po = profile[x+1]
					if x in  profile:
						if_cl.append(profile[x])
						if "highest_frame_diff_check" in data:
							if profile[x] > max(curr_mo, curr_po):
								#logging.debug("profile[x]", profile[x])
								#logging.debug("profile[x] is higher than curr_mo and curr_po", curr_mo, curr_po)
								if_cov += 1
						else:
							if profile[x] > min(curr_mo, curr_po):
								if_cov += 1
				if_cov = if_cov/if_len
				if if_cov != 0:
					if_cov = if_cov*100
				#In frame count discards the highest peak
				if locus in cds_count_dict:
					ratio = (sum(if_cl)/(stop-start))/(cds_count_dict[locus]+1.00)
				else:
					ratio = 0
				inframe_count = sum(sorted(if_cl)[:-1])
				lowest_frame_count = inframe_count - min(mo_count,po_count)
				high_frame_count = inframe_count - max(mo_count, po_count)
				before_start = 0
				after_start = 0
				#logging.debug("inframe count, lowest_frame_count, high_frame_count", inframe_count, lowest_frame_count, high_frame_count)
				for y in range(start-14,start+13,3):
					if y in profile:
						if y < start:
							before_start += profile[y]
						else:
							after_start += profile[y]
				start_score = after_start-before_start
				#logging.debug("before start, after start", before_start, after_start)
				before_stop = 0
				after_stop = 0
				for z in range(transcriptome_stop-13,transcriptome_stop+13,3):
					if z in profile:
						if z <= transcriptome_stop:
							before_stop += profile[z]
						else:
							after_stop += profile[z]
							
				stop_score = before_stop-after_stop
				'''
				#logging.debug("before stop, after stop", before_stop, after_stop)
				#logging.debug("if_cov", if_cov)
				if_cov = 0.0
				inframe_values = []
				minusone_values = []
				plusone_values = []
				for i in range(start-9,stop+12,3):
					if i-1 in profile:
						minusone_values.append(profile[i-1])
					else:
						minusone_values.append(0)
					if i in profile:
						inframe_values.append(profile[i])
					else:
						inframe_values.append(0)
					if i+1 in profile:
						plusone_values.append(profile[i+1])
					else:
						plusone_values.append(0)
				for x in range(4,len(inframe_values)-4):
					#print "x", x
					if "highest_frame_diff_check" in data:
						#print "highest frame check"
						if inframe_values[x] > max(minusone_values[x],plusone_values[x]):
							#print "{} is greater than {} and {}".format()
							if_cov += 1
					else:
						if inframe_values[x] > min(minusone_values[x],plusone_values[x]):
							if_cov += 1
				if locus in cds_count_dict:
					ratio = (sum(inframe_values)/(stop-start))/(cds_count_dict[locus]+1.00)
				else:
					ratio = 0	
				
					
				if_cov = if_cov/if_len
				if_cov = if_cov*100
				start_score = float(sum(inframe_values[4:8]))-float(sum(inframe_values[:4]))
				stop_score = float(sum(inframe_values[-8:-4]))-float(sum(inframe_values[-4:]))
				inframe_sum = float(sum(inframe_values[4:-4]))
				highframe_sum = float(max(sum(minusone_values[4:-4]),sum(plusone_values[4:-4])))
				lowframe_sum = float(min(sum(minusone_values[4:-4]),sum(plusone_values[4:-4])))
				if sum(minusone_values[4:-4]) == highframe_sum:
					highframe_list = minusone_values[4:-4]
				else:
					highframe_list = plusone_values[4:-4]

				if sum(minusone_values[4:-4]) == lowframe_sum:
					lowframe_list = minusone_values[4:-4]
				else:
					lowframe_list = plusone_values[4:-4]
					

				#Ratios 
				start_ratio = float(sum(inframe_values[4:8]))/float((sum(inframe_values[:4])+1))
				stop_ratio = float(sum(inframe_values[-8:-4]))/float((sum(inframe_values[-4:])+1))
				high_ratio = inframe_sum/(highframe_sum+1)
				low_ratio = inframe_sum/(lowframe_sum+1)
				#print start, stop , start_ratio, start_score, inframe_values[:4], inframe_values[4:8]
				#print start_ratio, stop_ratio, high_ratio, low_ratio
				
				high_frame_count = inframe_sum - highframe_sum
				lowest_frame_count = inframe_sum - lowframe_sum
				#BREAKPOINT
				
				final_score_values = []
				if "start_increase_check" in data:
					final_score_values.append(start_score)

				if "stop_decrease_check" in data:
					final_score_values.append(stop_score)

				if "lowest_frame_diff_check" in data:
					final_score_values.append(lowest_frame_count)

				if "highest_frame_diff_check" in data:
					final_score_values.append(high_frame_count)
				if "coverage_check" in data:
					final_score_values.append(if_cov)
				final_score_values.append(float(inframe_sum)/float(if_len))
				

				#TO DO: Instead of just summing these, they should be normalised over the current best value and then compared. The way it works now, coverage
				# has very little affect on the final outcome as it's a number below one. Instead normalise everything to whatever is highest between the current
				# value and the best value, e.g if current start score is 50 and best start score is 100, current start becomes 0.5 and best becomes 1,
				# another e.g if current coverage is 0.8 and best coverage is 0.2, current coverage becomes 1, best coverage becomes 0.25
				final_score = sum(final_score_values)
				#logging.debug("start", start)
				#logging.debug("final score", final_score)
				#print "start score", start_score
				#print "best_score", best_values["start_score"]
				#print "\nHHHHHHHIGHFRAEM VLAUE", high_frame_count
				if np.isnan(high_frame_count):
					high_frame_count = 1.0
				if np.isnan(lowest_frame_count):
					lowest_frame_count = 1.0
				if np.isnan(start_score):
					start_score = 1.0
				if np.isnan(stop_score):
					stop_score = 1.0
				#if transcript == "ENST00000341423":
				#	print "\nLLLLoowf rame count",start, stop, low_frame_count
				if start_score > best_values["start_score"] or (start_score == best_values["start_score"] and start < best_values["start"]):
					#if transcript == "ENST00000358435":
					#	print "updating start score",start_score
					best_values["start"] = start
					best_values["stop"] = transcriptome_stop
					best_values["transcript"] = transcript
					best_values["length"] = length
					best_values["high_frame_count"] = high_frame_count
					best_values["low_frame_count"] = lowest_frame_count
					best_values["start_score"] = start_score
					best_values["stop_score"] = stop_score
					best_values["final_score"] = final_score
					best_values["coverage"] = if_cov
					best_values["proteomics_count"] = proteomics_count
					best_values["ratio"] = ratio
					best_values["inframe_count"] = float(inframe_sum)/float(if_len)
					best_values["start_ratio"] = start_ratio
					best_values["stop_ratio"] = stop_ratio
					best_values["high_ratio"] = high_ratio
					best_values["low_ratio"] = low_ratio
					if high_frame_count > best_high_frame:
						best_high_frame = high_frame_count
					if lowest_frame_count > best_low_frame:
						best_low_frame = lowest_frame_count
					if stop_score > best_stop_score:
						best_stop_score = stop_score
					if start_score > best_start_score:
						best_start_score = start_score
					if if_cov > best_inframe_cov:
						best_inframe_cov = if_cov
			#if best_values["final_score"] > 0:
			if best_values["start"] <= 0:
				continue
			#if transcript == "ENST00000358435":
			#	print "best_values", best_values
			all_values.append([gene,
							best_values["transcript"],
							best_values["start"],
							best_values["stop"],
							best_values["length"],
							best_values["high_frame_count"],
							best_values["low_frame_count"],
							best_values["stop_score"],
							best_values["start_score"],
							best_values["coverage"],
							best_values["inframe_count"],
							0,
							0,
							0,
							0,
							0,
							0,
							0,
							orftype,
							best_values["proteomics_count"],
							best_values["ratio"],
							best_values["start_ratio"],
							best_values["stop_ratio"],
							best_values["high_ratio"],
							best_values["low_ratio"]])
	len_all_rows = float(len(all_values))
	#print "all values", all_values
	logging.debug("LENGTH OF ALL ROWS {}".format(len_all_rows))
	if len_all_rows == 0:
		return None
	'''
	sorted_all_values = sorted(all_values, key=lambda x: x[5],reverse=True)
	#FDR correction
	#logging.debug("FDR correcting")
	#logging.debug(str(sorted_all_values_no_corr))
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	#logging.debug("FDR correcting done")
	rank = 1
	prev_value = sorted_all_values[0][5]
	for row in sorted_all_values:
		curr_value = row[5]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[10] = rank
	sorted_all_values = sorted(all_values, key=lambda x: x[6],reverse=True)
	#FDR correction
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	rank = 1
	prev_value = sorted_all_values[0][6]
	for row in sorted_all_values:
		curr_value = row[6]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[11] = rank
	
	sorted_all_values = sorted(all_values, key=lambda x: x[7],reverse=True)
	#FDR correction
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	rank = 1
	prev_value = sorted_all_values[0][7]
	for row in sorted_all_values:
		curr_value = row[7]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[12] = rank
	
	sorted_all_values = sorted(all_values, key=lambda x: x[8],reverse=True)
	#FDR correction
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	rank = 1
	prev_value = sorted_all_values[0][8]
	for row in sorted_all_values:
		curr_value = row[8]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[13] = rank

	sorted_all_values = sorted(all_values, key=lambda x: x[9],reverse=True)
	rank = 1
	prev_value = sorted_all_values[0][9]
	for row in sorted_all_values:
		curr_value = row[9]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[14] = rank
	
	'''
	
	#FDR correction
	'''
	hfc = []
	lfc = []
	stop_s = []
	start_s = []
	
	for row in all_values:
		hfc.append(row[5])
		lfc.append(row[6])
		stop_s.append(row[7])
		start_s.append(row[8])
	logging.debug(hfc[:10])
	logging.debug(lfc[:10])
	logging.debug(stop_s[:10])
	logging.debug(start_s[:10])
	hfc_corr = multipletests(hfc, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	lfc_corr = multipletests(lfc, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	#print "STOP SCORES INPUT", stop_s
	stop_corr = multipletests(stop_s, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	#print "stop SCORES OUTPUT", stop_corr
	start_corr = multipletests(start_s, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	logging.debug(hfc_corr[1][:10])
	logging.debug(lfc_corr[1][:10])
	logging.debug(stop_corr[1][:10])
	logging.debug(start_corr[1][:10])
	
	logging.debug("first rows")
	#for row in all_values:
	#	if row[1] == "ENST00000594768":
	#		logging.debug(row)
	for i in range(0,len(all_values)):
		all_values[i][5] = hfc_corr[1][i] 
		all_values[i][6] = lfc_corr[1][i] 
		all_values[i][7] = stop_corr[1][i] 
		all_values[i][8] = start_corr[1][i] 
	#for row in all_values:
	#	if row[1] == "ENST00000594768":
	#		logging.debug(row)
	'''
	
	
	
	
	
	
	sorted_all_values = sorted(all_values, key=lambda x: x[5],reverse=False)
	rank = 1
	prev_value = sorted_all_values[0][5]
	for row in sorted_all_values:
		curr_value = row[5]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[11] = rank
	sorted_all_values = sorted(all_values, key=lambda x: x[6],reverse=False)
	#FDR correction
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	rank = 1
	prev_value = sorted_all_values[0][6]
	for row in sorted_all_values:
		curr_value = row[6]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[12] = rank
	
	sorted_all_values = sorted(all_values, key=lambda x: x[7],reverse=False)
	#FDR correction
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	rank = 1
	prev_value = sorted_all_values[0][7]
	for row in sorted_all_values:
		curr_value = row[7]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[13] = rank
	
	sorted_all_values = sorted(all_values, key=lambda x: x[8],reverse=False)
	#FDR correction
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	rank = 1
	prev_value = sorted_all_values[0][8]
	for row in sorted_all_values:
		curr_value = row[8]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[14] = rank

	sorted_all_values = sorted(all_values, key=lambda x: x[9],reverse=True)
	rank = 1
	prev_value = sorted_all_values[0][9]
	for row in sorted_all_values:
		curr_value = row[9]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[15] = rank
		
	
	sorted_all_values = sorted(all_values, key=lambda x: x[10],reverse=True)
	rank = 1
	prev_value = sorted_all_values[0][10]
	for row in sorted_all_values:
		curr_value = row[10]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[16] = rank
		
		
		

	for row in all_values:
		normalised_score_values = []
		#inframe count
		normalised_score_values.append(row[16])
		if "coverage_check" in data:
			normalised_score_values.append(row[15])
		if "start_increase_check" in data:
			#if row[13] >0.05:
			#	continue
			normalised_score_values.append(row[14])
		if "stop_decrease_check" in data:
			#if row[12] >0.05:
			#	continue
			normalised_score_values.append(row[13])
		if "lowest_frame_diff_check" in data:
			#if row[11] >0.05:
			#	continue
			normalised_score_values.append(row[12])
		if "highest_frame_diff_check" in data:
			#if row[10] >0.05:
			#	continue
			normalised_score_values.append(row[11])
		normalised_score = sum(normalised_score_values)
		#logging.debug("normalised_score_values, normalised score", normalised_score_values, normalised_score)
		row.append(round(normalised_score,2))
	logging.debug("sorting all values")
	sorted_all_values = sorted(all_values, key=lambda x: x[-1],reverse=False)
	logging.debug("LENGHT OF ALL SORTED ROWS {}".format(len(sorted_all_values)))
	final_rank = 1
	for tup in sorted_all_values:
		tup[17] = final_rank
		final_rank += 1
	return sorted_all_values	
	
	
	
	
	
	
	
	
	
	
	
	
def extract_values(accepted_orf_dict,data,tran_gene_dict,selected_seq_types,profile_dict,cds_count_dict):
	#logging.debug("extract values called")
	best_high_frame = 0
	best_low_frame = 0
	best_start_score = 0
	best_stop_score = 0
	best_inframe_cov = 0
	all_values = []
	tot_loc = 0
	for locus in accepted_orf_dict:
		tot_loc += 1
		if tot_loc%100 == 0:
			logging.debug("total transcripts {}".format(tot_loc))
		for stop in accepted_orf_dict[locus]:
			best_values = {"start":-1,"high_frame_count":1,"low_frame_count":1,"start_score":-10,"stop_score":1,"final_score":-10000,"coverage":0,"length":0,
				  "transcript":locus,"stop":0,"proteomics_count":0,"ratio":0,"inframe_count":0,"start_ratio":0,"stop_ratio":0,"high_ratio":0,"low_ratio":0,
				  "stop_geo":0,"start_geo":0}
			for start in accepted_orf_dict[locus][stop]:
				#print "start stop", start, stop
				transcript = accepted_orf_dict[locus][stop][start]["transcript"]
				best_values["transcript"] = transcript
				if best_values["start"] == -1:
					best_values["start"] = start
					best_values["stop"] = stop
				gene = tran_gene_dict[transcript]
				if "riboseq" in selected_seq_types and "proteomics" in selected_seq_types:
					profile = profile_dict[transcript]["riboseq"]
					prot_profile = profile_dict[transcript]["proteomics"]
				elif "proteomics" not in selected_seq_types:
					profile = profile_dict[transcript]["riboseq"]
					prot_profile = {}
				else:
					profile = profile_dict[transcript]["proteomics"]
					prot_profile = {}
				transcriptome_stop = accepted_orf_dict[locus][stop][start]["stop"]
				orftype = accepted_orf_dict[locus][stop][start]["orftype"]
				length = (transcriptome_stop-start)+1
				inframe_count = 0
				proteomics_count = 0
				if_cl = []
				mo_count = 0
				po_count = 0
				if_cov = 0.0
				if_len = (stop-start)/3
				
				for p in range(start+1, transcriptome_stop,3):
					if p in prot_profile:
						proteomics_count += prot_profile[p]
				#BREAKPOINT
				'''
				for x in range(start+6, transcriptome_stop-9,3):
					#logging.debug("x", x)
					curr_mo = 0
					curr_po = 0
					if_len += 1
					if x-1 in profile:
						mo_count += profile[x-1]
						curr_mo = profile[x-1]
					if x+1 in profile:
						po_count += profile[x+1]
						curr_po = profile[x+1]
					if x in  profile:
						if_cl.append(profile[x])
						if "highest_frame_diff_check" in data:
							if profile[x] > max(curr_mo, curr_po):
								#logging.debug("profile[x]", profile[x])
								#logging.debug("profile[x] is higher than curr_mo and curr_po", curr_mo, curr_po)
								if_cov += 1
						else:
							if profile[x] > min(curr_mo, curr_po):
								if_cov += 1
				if_cov = if_cov/if_len
				if if_cov != 0:
					if_cov = if_cov*100
				#In frame count discards the highest peak
				if locus in cds_count_dict:
					ratio = (sum(if_cl)/(stop-start))/(cds_count_dict[locus]+1.00)
				else:
					ratio = 0
				inframe_count = sum(sorted(if_cl)[:-1])
				lowest_frame_count = inframe_count - min(mo_count,po_count)
				high_frame_count = inframe_count - max(mo_count, po_count)
				before_start = 0
				after_start = 0
				#logging.debug("inframe count, lowest_frame_count, high_frame_count", inframe_count, lowest_frame_count, high_frame_count)
				for y in range(start-14,start+13,3):
					if y in profile:
						if y < start:
							before_start += profile[y]
						else:
							after_start += profile[y]
				start_score = after_start-before_start
				#logging.debug("before start, after start", before_start, after_start)
				before_stop = 0
				after_stop = 0
				for z in range(transcriptome_stop-13,transcriptome_stop+13,3):
					if z in profile:
						if z <= transcriptome_stop:
							before_stop += profile[z]
						else:
							after_stop += profile[z]
							
				stop_score = before_stop-after_stop
				'''
				#logging.debug("before stop, after stop", before_stop, after_stop)
				#logging.debug("if_cov", if_cov)
				if_cov = 0.0
				inframe_values = []
				minusone_values = []
				plusone_values = []
				for i in range(start-9,stop+12,3):
					if i-1 in profile:
						minusone_values.append(profile[i-1])
					else:
						minusone_values.append(0)
					if i in profile:
						inframe_values.append(profile[i])
					else:
						inframe_values.append(0)
					if i+1 in profile:
						plusone_values.append(profile[i+1])
					else:
						plusone_values.append(0)
				for x in range(4,len(inframe_values)-4):
					#print "x", x
					if "highest_frame_diff_check" in data:
						#print "highest frame check"
						if inframe_values[x] > max(minusone_values[x],plusone_values[x]):
							#print "{} is greater than {} and {}".format()
							if_cov += 1
					else:
						if inframe_values[x] > min(minusone_values[x],plusone_values[x]):
							if_cov += 1
				if locus in cds_count_dict:
					ratio = (sum(inframe_values)/(stop-start))/(cds_count_dict[locus]+1.00)
				else:
					ratio = 0	
				
				
					
					
					
				if_cov = if_cov/if_len
				if_cov = if_cov*100
				#if inframe_values[:4] == inframe_values[4:8]:
				#	start_score = 1
				#else:
				#	stat, start_score = mannwhitneyu(inframe_values[:4], inframe_values[4:8],alternative="less")
				#	#if transcript == "ENST00000358435":
				#	#	print ("start score, inframe_values[:4], inframe_values[4:8]", start, stop, start_score, inframe_values[:4], inframe_values[4:8])
				start_score_raw = (float(sum(inframe_values[4:8]))+1)/(float(sum(inframe_values[:4]))+1)
				start_score = np.log(start_score_raw)
				start_geo = geo_mean(inframe_values[:8])
				#logging.debug(inframe_values[:4], inframe_values[4:8], start_score)
				#start_score = sum(inframe_values[4:8])-sum(inframe_values[:4])
				stop_score_raw = (float(sum(inframe_values[-8:-4]))+1)/(float(sum(inframe_values[-4:]))+1)
				stop_score = np.log(stop_score_raw)
				if stop_score < -1000:
					print stop_score_raw, stop_score, sum(inframe_values[-8:-4])+1, sum(inframe_values[-4:])+1
				stop_geo = geo_mean(inframe_values[-8:])
				
				#stop_score = sum(inframe_values[-8:-4])-sum(inframe_values[-4:])
				inframe_sum = float(sum(inframe_values[4:-4]))
				highframe_sum = float(max(sum(minusone_values[4:-4]),sum(plusone_values[4:-4])))
				lowframe_sum = float(min(sum(minusone_values[4:-4]),sum(plusone_values[4:-4])))
				if sum(minusone_values[4:-4]) == highframe_sum:
					highframe_list = minusone_values[4:-4]
				else:
					highframe_list = plusone_values[4:-4]
					
				
				stat, high_frame_count = wilcoxon(inframe_values[4:-4],highframe_list)
				if inframe_sum > highframe_sum:
					high_frame_count = high_frame_count/2
				else:
					high_frame_count =  1-(high_frame_count/2)
				
				
				
				if sum(minusone_values[4:-4]) == lowframe_sum:
					lowframe_list = minusone_values[4:-4]
				else:
					lowframe_list = plusone_values[4:-4]
					
				
				stat, low_frame_count = wilcoxon(inframe_values[4:-4],lowframe_list)
				if inframe_sum > lowframe_sum:
					lowest_frame_count = low_frame_count/2
				else:
					lowest_frame_count =  1-(low_frame_count/2)
				
				#Ratios 
				start_ratio = float(sum(inframe_values[4:8]))/float((sum(inframe_values[:4])+1))
				stop_ratio = float(sum(inframe_values[-8:-4]))/float((sum(inframe_values[-4:])+1))
				high_ratio = inframe_sum/(highframe_sum+1)
				low_ratio = inframe_sum/(lowframe_sum+1)
				#print start, stop , start_ratio, start_score, inframe_values[:4], inframe_values[4:8]
				#print start_ratio, stop_ratio, high_ratio, low_ratio
				if transcript == "ENST00000309311" and start == 90:
					print "start_score_raw, start_score, start_geo", start_score_raw, start_score, start_geo, inframe_values[4:8], inframe_values[:4]
					for y in range(80,100):
						if y in profile:
							print "y", y, profile[y]
						else:
							print "y", y, 0
				#high_frame_count = sum(inframe_values[4:-4]) - max(sum(minusone_values[4:-4]),sum(plusone_values[4:-4]))
				#lowest_frame_count =  - min(sum(minusone_values[4:-4]),sum(plusone_values[4:-4]))
				#BREAKPOINT
				
				final_score_values = []
				if "start_increase_check" in data:
					final_score_values.append(start_score)

				if "stop_decrease_check" in data:
					final_score_values.append(stop_score)

				if "lowest_frame_diff_check" in data:
					final_score_values.append(lowest_frame_count)

				if "highest_frame_diff_check" in data:
					final_score_values.append(high_frame_count)
				if "coverage_check" in data:
					final_score_values.append(if_cov)
				final_score_values.append(float(inframe_sum)/float(if_len))
				

				#TO DO: Instead of just summing these, they should be normalised over the current best value and then compared. The way it works now, coverage
				# has very little affect on the final outcome as it's a number below one. Instead normalise everything to whatever is highest between the current
				# value and the best value, e.g if current start score is 50 and best start score is 100, current start becomes 0.5 and best becomes 1,
				# another e.g if current coverage is 0.8 and best coverage is 0.2, current coverage becomes 1, best coverage becomes 0.25
				final_score = sum(final_score_values)
				#logging.debug("start", start)
				#logging.debug("final score", final_score)
				#print "start score", start_score
				#print "best_score", best_values["start_score"]
				#print "\nHHHHHHHIGHFRAEM VLAUE", high_frame_count
				if np.isnan(high_frame_count):
					high_frame_count = 1.0
				if np.isnan(lowest_frame_count):
					lowest_frame_count = 1.0
				if np.isnan(start_score):
					start_score = 1.0
				if np.isnan(stop_score):
					stop_score = 1.0
				#if transcript == "ENST00000341423":
				#	print "\nLLLLoowf rame count",start, stop, low_frame_count
				if start_score > best_values["start_score"] or (start_score == best_values["start_score"] and start < best_values["start"]):
					#if transcript == "ENST00000358435":
					#	print "updating start score",start_score
					best_values["start"] = start
					best_values["stop"] = transcriptome_stop
					best_values["transcript"] = transcript
					best_values["length"] = length
					best_values["high_frame_count"] = high_frame_count
					best_values["low_frame_count"] = lowest_frame_count
					best_values["start_score"] = start_score
					best_values["stop_score"] = stop_score
					best_values["final_score"] = final_score
					best_values["coverage"] = if_cov
					best_values["proteomics_count"] = proteomics_count
					best_values["ratio"] = ratio
					best_values["inframe_count"] = float(inframe_sum)/float(if_len)
					best_values["start_ratio"] = start_ratio
					best_values["stop_ratio"] = stop_ratio
					best_values["stop_geo"] = stop_geo
					best_values["start_geo"] = start_geo
					best_values["high_ratio"] = high_ratio
					best_values["low_ratio"] = low_ratio
					if high_frame_count > best_high_frame:
						best_high_frame = high_frame_count
					if lowest_frame_count > best_low_frame:
						best_low_frame = lowest_frame_count
					if stop_score > best_stop_score:
						best_stop_score = stop_score
					if start_score > best_start_score:
						best_start_score = start_score
					if if_cov > best_inframe_cov:
						best_inframe_cov = if_cov
			#if best_values["final_score"] > 0:
		
			#if transcript == "ENST00000358435":
			#	print "best_values", best_values

			all_values.append([gene,
							best_values["transcript"],
							best_values["start"],
							best_values["stop"],
							best_values["length"],
							best_values["high_frame_count"],
							best_values["low_frame_count"],
							best_values["stop_score"],
							best_values["start_score"],
							best_values["coverage"],
							best_values["inframe_count"],
							0,
							0,
							0,
							0,
							0,
							0,
							0,
							orftype,
							best_values["proteomics_count"],
							best_values["ratio"],
							best_values["start_ratio"],
							best_values["stop_ratio"],
							best_values["high_ratio"],
							best_values["low_ratio"],
							best_values["start_geo"],
							best_values["stop_geo"]])
	len_all_rows = float(len(all_values))
	#print "all values", all_values
	logging.debug("LENGTH OF ALL ROWS {}".format(len_all_rows))
	if len_all_rows == 0:
		return None
	'''
	sorted_all_values = sorted(all_values, key=lambda x: x[5],reverse=True)
	#FDR correction
	#logging.debug("FDR correcting")
	#logging.debug(str(sorted_all_values_no_corr))
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	#logging.debug("FDR correcting done")
	rank = 1
	prev_value = sorted_all_values[0][5]
	for row in sorted_all_values:
		curr_value = row[5]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[10] = rank
	sorted_all_values = sorted(all_values, key=lambda x: x[6],reverse=True)
	#FDR correction
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	rank = 1
	prev_value = sorted_all_values[0][6]
	for row in sorted_all_values:
		curr_value = row[6]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[11] = rank
	
	sorted_all_values = sorted(all_values, key=lambda x: x[7],reverse=True)
	#FDR correction
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	rank = 1
	prev_value = sorted_all_values[0][7]
	for row in sorted_all_values:
		curr_value = row[7]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[12] = rank
	
	sorted_all_values = sorted(all_values, key=lambda x: x[8],reverse=True)
	#FDR correction
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	rank = 1
	prev_value = sorted_all_values[0][8]
	for row in sorted_all_values:
		curr_value = row[8]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[13] = rank

	sorted_all_values = sorted(all_values, key=lambda x: x[9],reverse=True)
	rank = 1
	prev_value = sorted_all_values[0][9]
	for row in sorted_all_values:
		curr_value = row[9]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[14] = rank
	
	'''
	
	#FDR correction
	'''
	hfc = []
	lfc = []
	
	for row in all_values:
		hfc.append(row[5])
		lfc.append(row[6])
	hfc_corr = multipletests(hfc, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	lfc_corr = multipletests(lfc, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	
	for i in range(0,len(all_values)):
		all_values[i][5] = hfc_corr[1][i] 
		all_values[i][6] = lfc_corr[1][i] 
	'''

	
	#Calculate z-scores for start score/stop_score
	sorted_all_values = sorted(all_values, key=lambda x: x[25],reverse=True)
	#print len(sorted_all_values)
	for i in range(0,len(sorted_all_values),300):
		#print "i", i
		start_scores = []
		for x in range(i, i+300):
			try:
				start_scores.append(sorted_all_values[x][8])
			except:
				break
			
		#if there is less than 300 scores left in all values, keep adding to start_scores:
		#print len(sorted_all_values)-i
		if len(sorted_all_values)-i < 300:
			#print "adding extras"
			for y in range(i+300, len(sorted_all_values)):
				start_scores.append(sorted_all_values[y][8])
		#print "len start scores", len(start_scores)
		if len(start_scores) < 300:
			#print len(start_scores)
			continue
		
		
		z_scores = zscore(start_scores)
		ind = 0
		for z_score in z_scores:
			#print "i+ind", i+ind
			sorted_all_values[i+ind][8] = z_score
			ind+= 1
			
	sorted_all_values = sorted(all_values, key=lambda x: x[26],reverse=True)
	for i in range(0,len(sorted_all_values),300):
		stop_scores = []
		for x in range(i, i+300):
			try:
				stop_scores.append(sorted_all_values[x][7])
			except:
				break
			
		#if there is less than 300 scores left in all values, keep adding to start_scores:
		if len(sorted_all_values)-i < 300:
			for y in range(i+300, len(sorted_all_values)):
				stop_scores.append(sorted_all_values[y][8])
		if len(stop_scores) < 300:
			print len(stop_scores)
			continue
		
		
		z_scores = zscore(stop_scores)
		#print "sotp scores", stop_scores
		#print "z_scores", z_scores
		ind = 0
		for z_score in z_scores:
			sorted_all_values[i+ind][7] = z_score
			ind+= 1		
	
	
	
	sorted_all_values = sorted(all_values, key=lambda x: x[5],reverse=False)
	rank = 1
	prev_value = sorted_all_values[0][5]
	for row in sorted_all_values:
		curr_value = row[5]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[11] = rank
	sorted_all_values = sorted(all_values, key=lambda x: x[6],reverse=False)
	#FDR correction
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	rank = 1
	prev_value = sorted_all_values[0][6]
	for row in sorted_all_values:
		curr_value = row[6]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[12] = rank
	
	sorted_all_values = sorted(all_values, key=lambda x: x[7],reverse=True)
	#FDR correction
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	rank = 1
	prev_value = sorted_all_values[0][7]
	for row in sorted_all_values:
		curr_value = row[7]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[13] = rank
	
	sorted_all_values = sorted(all_values, key=lambda x: x[8],reverse=True)
	#FDR correction
	#sorted_all_values = multipletests(sorted_all_values_no_corr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	rank = 1
	prev_value = sorted_all_values[0][8]
	for row in sorted_all_values:
		curr_value = row[8]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[14] = rank

	sorted_all_values = sorted(all_values, key=lambda x: x[9],reverse=True)
	rank = 1
	prev_value = sorted_all_values[0][9]
	for row in sorted_all_values:
		curr_value = row[9]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[15] = rank
		
	
	sorted_all_values = sorted(all_values, key=lambda x: x[10],reverse=True)
	rank = 1
	prev_value = sorted_all_values[0][10]
	for row in sorted_all_values:
		curr_value = row[10]
		if curr_value != prev_value:
			rank += 1
			prev_value = curr_value
		row[16] = rank
		
		
		

	for row in all_values:
		normalised_score_values = []
		#inframe count
		normalised_score_values.append(row[16])
		if "coverage_check" in data:
			normalised_score_values.append(row[15])
		if "start_increase_check" in data:
			#if row[13] >0.05:
			#	continue
			normalised_score_values.append(row[14])
		if "stop_decrease_check" in data:
			#if row[12] >0.05:
			#	continue
			normalised_score_values.append(row[13])
		if "lowest_frame_diff_check" in data:
			#if row[11] >0.05:
			#	continue
			normalised_score_values.append(row[12])
		if "highest_frame_diff_check" in data:
			#if row[10] >0.05:
			#	continue
			normalised_score_values.append(row[11])
		normalised_score = sum(normalised_score_values)
		#logging.debug("normalised_score_values, normalised score", normalised_score_values, normalised_score)
		row.append(round(normalised_score,2))
	logging.debug("sorting all values")
	sorted_all_values = sorted(all_values, key=lambda x: x[-1],reverse=False)
	logging.debug("LENGHT OF ALL SORTED ROWS {}".format(len(sorted_all_values)))
	final_rank = 1
	for tup in sorted_all_values:
		tup[17] = final_rank
		final_rank += 1
	return sorted_all_values


def write_to_file(sorted_all_values,filename,sequence_dict,organism,transcriptome,file_string):
	#logging.debug("all sorted all values", sorted_all_values)
	returnstr = "Table|"
	tmp_result_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")
	tmp_result_file.write("Gene,Tran,Start,Stop,Length,Global_Rank,Type,Trips-viz link,Start Codon,Highframe rank,Highframe value,Lowframe rank,Lowframe value,Stop rank,Stop value,Start rank,Start value,Coverage rank,Coverage value,Inframe Count Rank,Inframe Count Value,Amino acid sequence,Proteomics count,CDS ratio,Start_ratio, Stop_ratio, Highframe_ratio, Lowframe_ratio\n")
	tup_count = 0
	#logging.debug("writing to file",len(sorted_all_values))
	for tup in sorted_all_values:
		#logging.debug("tup", tup)
		gene = tup[0]
		transcript = tup[1]
		start = tup[2]
		stop = tup[3]
		length = tup[4]
		global_rank = tup[17]
		orftype = tup[18]
		proteomics_count = tup[19]
		ratio = tup[20]
		#logging.debug("transcript", transcript, "start", start, "stop", stop)
		try:
			seq = sequence_dict[transcript][start-1:stop]
		except:
			seq = ""
		start_codon = seq[:3]
		# Get amino acid sequence of this ORF, excluding the stop codon. 
		while len(seq)%3 != 0:
			seq = seq[:-1]

		aa_seq = nuc_to_aa(seq[:-3])
		trips_link = '<a href="https://trips.ucc.ie/'+organism+'/'+transcriptome+'/interactive_plot/?tran='+transcript+'&hili='+str(start)+"_"+str(stop)+'&files='+file_string+'" target="_blank_" >View on trips-viz</a>'
		tmp_result_file.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(gene,transcript,start,stop,length,global_rank,orftype,trips_link,
																								  start_codon,tup[11],tup[5],tup[12],tup[6],tup[13],tup[7],tup[14],tup[8],tup[15],tup[9],tup[16],tup[10],aa_seq,proteomics_count,ratio,tup[21],tup[22],tup[23],tup[24]))
		if tup_count < 1000:
			returnstr += "{},{},{},{},{},{},{},NULL,NULL,NULL,NULL,{}.,/".format(gene,transcript,start,stop,length,global_rank,orftype,trips_link)
		tup_count += 1	
		#if tup_count%100 == 0:
		#	logging.debug("tup count", tup_count)
	return returnstr


def create_model(columns):
    # create model
    model = Sequential()
    model.add(BatchNormalization(axis=-1,momentum=0.99,epsilon=0.001,center=True,scale=True))
    model.add(Dense(columns, input_dim=(columns-1), kernel_initializer='uniform', activation='elu'))
    #model.add(Dropout(0.1))
    model.add(Dense(columns, kernel_initializer='uniform', activation='elu'))
    model.add(Dense(columns, kernel_initializer='uniform', activation='elu'))
    #model.add(Dense(columns, kernel_initializer='uniform', activation='elu'))
    #model.add(Dense(columns, kernel_initializer='uniform', activation='elu'))
    #model.add(Dense(columns, kernel_initializer='uniform', activation='elu'))
    #model.add(Dropout(0.5))
    model.add(Dense(1, kernel_initializer='uniform', activation='sigmoid'))
    opt = optimizers.Adam(learning_rate=0.001)
    # Compile model
    model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
    return model


def neural_net(sorted_all_values,organism, transcriptome, file_string,training_filename,feature_list,filename):
	returnstr = "Table|"
	columns = len(feature_list)
	model = create_model(columns)
	openfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,training_filename)).readlines()
	headers = openfile[0].split(",")
	index_list = []
	for feature in feature_list:
		if feature in headers:
			ind = headers.index(feature)
			index_list.append(ind)
		else:
			return "Could not find feature {} in training dataset: {}".format(feature, training_filename)

	
	dataset = np.loadtxt("{}/static/tmp/{}".format(config.SCRIPT_LOC,training_filename), delimiter=",",skiprows=1,usecols=index_list)
	
	# split into input (X) and output (Y) variables
	X = dataset[:,1:columns]
	Y = dataset[:,0]
	randomize = np.arange(len(X))
	np.random.shuffle(randomize)
	X = X[randomize]
	Y = Y[randomize]

	# Fit the model
	model.fit(X, Y, nb_epoch=500, batch_size=100, validation_split=0.2)

	# evaluate the model
	scores = model.evaluate(X, Y)
	#print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100))

	testfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename)).readlines()
	outfile = open("{}/static/tmp/{}.nnet.csv".format(config.SCRIPT_LOC,filename),"w")
	outfile.write(testfile[0])
	testset = np.loadtxt("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename), delimiter=",",skiprows=1,usecols=index_list[1:])
	X = testset[:,0:columns]
	predictions = model.predict(X)
	tup_count = 0
	for row in testfile[1:]:
		splitrow = row.split(",")
		prediction = predictions[tup_count][0]
		outfile.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(splitrow[0],splitrow[1],splitrow[2],splitrow[3],splitrow[4],prediction,splitrow[6],splitrow[7],splitrow[8],splitrow[9],splitrow[10],splitrow[11],splitrow[12],splitrow[13],splitrow[14],splitrow[15],splitrow[16],splitrow[17],splitrow[18],splitrow[19],splitrow[20],splitrow[21],splitrow[22]))
		if tup_count < 1000:	
			returnstr += "{},{},{},{},{},{},{},NULL,NULL,NULL,NULL,{}.,/".format(splitrow[0],splitrow[1],splitrow[2],splitrow[3],splitrow[4],prediction,splitrow[6],splitrow[7])
		tup_count += 1

	original_filepath = "{}/static/tmp/{}".format(config.SCRIPT_LOC,filename)
	new_filepath = "{}/static/tmp/{}.nnet.csv".format(config.SCRIPT_LOC,filename)
	os.rename(new_filepath,original_filepath)
	return returnstr



def create_training_dict(top_cases,cds_dict):
	top_cases_len = len(top_cases)
	batch_size = top_cases_len/4
	accepted_orf_dict = {}
	#Positive cases
	for tran in top_cases:
		cds_start = cds_dict[tran]["cds_start"]
		cds_stop = cds_dict[tran]["cds_stop"]
		length = cds_stop-cds_start
		cds_cov = 1
		start_codon = "ATG"
		accepted_orf_dict[tran] = {cds_stop:{}}
		accepted_orf_dict[tran][cds_stop][cds_start] = {"length":length,
												"score":0,
												"cds_cov":cds_cov,
												"start_codon":start_codon,
												"orftype":"1",
												"stop":cds_stop,
												"transcript":tran}

	#Negative cases, wrong start
	for tran in top_cases[:batch_size]:
		cds_start = random.randint(4,cds_dict[tran]["cds_start"]-1)
		cds_stop = cds_dict[tran]["cds_stop"]+3
		length = cds_stop-cds_start
		cds_cov = 1
		start_codon = "ATG"
		accepted_orf_dict[tran][cds_stop] = {}
		accepted_orf_dict[tran][cds_stop][cds_start] = {"length":length,
													"score":0,
													"cds_cov":cds_cov,
													"start_codon":start_codon,
													"orftype":"0",
													"stop":cds_stop,
													"transcript":tran}
	
	#Negative cases, wrong stop
	for tran in top_cases[batch_size:batch_size*2]:
		cds_start = cds_dict[tran]["cds_start"]
		cds_stop = random.randint(cds_dict[tran]["cds_stop"]+4,cds_dict[tran]["length"]-4)
		length = cds_stop-cds_start
		cds_cov = 1
		start_codon = "ATG"
		accepted_orf_dict[tran][cds_stop] = {}
		accepted_orf_dict[tran][cds_stop][cds_start] = {"length":length,
													"score":0,
													"cds_cov":cds_cov,
													"start_codon":start_codon,
													"orftype":"0",
													"stop":cds_stop,
													"transcript":tran}
	#Negative cases wrong frame
	for tran in top_cases[batch_size*2:batch_size*3]:
		cds_start = random.choice([cds_dict[tran]["cds_start"]-1,cds_dict[tran]["cds_start"]+1])
		cds_stop = cds_dict[tran]["cds_stop"]+5
		length = cds_stop-cds_start
		cds_cov = 1
		start_codon = "ATG"
		accepted_orf_dict[tran][cds_stop] = {}
		accepted_orf_dict[tran][cds_stop][cds_start] = {"length":length,
													"score":0,
													"cds_cov":cds_cov,
													"start_codon":start_codon,
													"orftype":"0",
													"stop":cds_stop,
													"transcript":tran}
	#Negative cases, low reads (3' trailer)
	for tran in top_cases[batch_size*3:batch_size*4]:
		cds_start = cds_dict[tran]["cds_stop"]+random.randint(1,100)
		cds_stop = random.randint(cds_start+101,cds_start+300)
		length = cds_stop-cds_start
		cds_cov = 1
		start_codon = "ATG"
		accepted_orf_dict[tran][cds_stop] = {}
		accepted_orf_dict[tran][cds_stop][cds_start] = {"length":length,
													"score":0,
													"cds_cov":cds_cov,
													"start_codon":start_codon,
													"orftype":"0",
													"stop":cds_stop,
													"transcript":tran}

	return accepted_orf_dict

# Returns a table with ranked orf scores
orfquery_blueprint = Blueprint("orfquery", __name__, template_folder="templates")
@orfquery_blueprint.route('/orfquery', methods=['POST'])
def orfquery():
	logging.debug("orfquery called")
	global user_short_passed
	data = ast.literal_eval(request.data)
	connection = sqlite3.connect("{}/trips.sqlite".format(config.SCRIPT_LOC))
	cursor = connection.cursor()
	organism = data["organism"]
	transcriptome = data["transcriptome"]
	cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
	owner = (cursor.fetchone())[0]
	if owner == 1:
		if os.path.isfile("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome)):
			traninfo_connection = sqlite3.connect("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome))
		else:
			return "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
	else:
		traninfo_connection = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
	start_time = time.time()
	acceptable = 0
	
	minscore = float(data["minscore"])
	try:
		user = current_user.name
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]
	except:
		user = None

	file_paths_dict = fetch_file_paths(data["file_list"],organism)
	#Find out which studies have all files of a specific sequence type selected (to create aggregates)
	full_studies = []
	if minscore == 0:
		for seq_type in file_paths_dict:
			all_study_ids = []
			all_file_ids = file_paths_dict[seq_type].keys()
			cursor.execute("SELECT DISTINCT study_id FROM files WHERE file_id IN ({})".format(str(all_file_ids).strip("[").strip("]")))
			result = cursor.fetchall()
			for row in result:
				all_study_ids.append(row[0])
			logging.debug("for seq type {} here are the study ids: {}".format(seq_type, all_study_ids))
			for study_id in all_study_ids:
				cursor.execute("SELECT file_id FROM files WHERE study_id = {} AND file_type = '{}'".format(study_id,seq_type))
				result = cursor.fetchall()
				study_file_ids = []
				for row in result:
					study_file_ids.append(row[0])
				logging.debug("study id {}".format(study_id))
				#Do not create an aggregate if there is only one file in the study:
				if len(study_file_ids) > 1:
					logging.debug("length over 1")
					#Check if all the file_ids in this study are present in all_file_ids, if they are this study is "full" and an aggreagate will be made
					if set(study_file_ids).issubset(all_file_ids):
						#This is a full study, create an aggregate for it if non exits
						#Get path to study directory
						study_path = "/".join((file_paths_dict[seq_type][study_file_ids[0]]).split("/")[:-1])
						logging.debug("study_path {}".format(study_path))
						full_studies.append(study_id)
						#Check if aggregate file exists
						if not os.path.isfile("{}/aggregate_{}.sqlite".format(study_path, seq_type)):
							sub_file_paths_dict = {seq_type:{}}
							for file_id in study_file_ids:
								sub_file_paths_dict[seq_type][file_id] = file_paths_dict[seq_type][file_id]
							create_aggregate(sub_file_paths_dict,study_path, seq_type)
							logging.debug("aggregate created")
						#Remove all file_ids associated with this study from file_paths_dict
						for file_id in study_file_ids:
							del file_paths_dict[seq_type][file_id]
						#Add the aggregate to file_paths_dict
						file_paths_dict[seq_type]["STUDY_"+str(study_id)] = "{}/aggregate_{}.sqlite".format(study_path, seq_type)
						
						
	connection.close()				
	logging.debug("Full studies {}".format(full_studies))
	#return str(full_studies)
	#return str(all_study_ids)
		
	html_args = data["html_args"]
	returnstr = ""
	
	start_codons = []
	
	#Used to extract columns from csv for neural neural_net
	feature_list = ["Type"]
	if "sc_aug" in data:
		start_codons.append("ATG")
	if "sc_cug" in data:
		start_codons.append("CTG")
	if "sc_gug" in data:
		start_codons.append("GTG")
	if "sc_none" in data:
		start_codons.append("any")
	#start_codons.append("any")
	logging.debug("start codons {}".format(start_codons))
	#min_start_increase = float(data["min_start_increase"])
	#max_start_increase = float(data["max_start_increase"])

	#min_stop_decrease = float(data["min_stop_decrease"])
	#max_stop_decrease = float(data["max_stop_decrease"])

	#min_coverage = float(data["min_coverage"])
	#max_coverage = float(data["max_coverage"])

	#min_lowest_frame_diff = float(data["min_lowest_frame_diff"])
	#max_lowest_frame_diff = float(data["max_lowest_frame_diff"])

	#min_highest_frame_diff = float(data["min_highest_frame_diff"])
	#max_highest_frame_diff = float(data["max_highest_frame_diff"])
	output = data["output"]
	
	min_cds = float(data["min_cds"])
	max_cds = float(data["max_cds"])

	min_len = float(data["min_len"])
	if data["max_len"] == "nnet":
		output = "nnet"
		max_len = 10000
	else:
		max_len = float(data["max_len"])

	min_avg = float(data["min_avg"])
	max_avg = float(data["max_avg"])


	accepted_orftypes = []
	
	#Checkbox implementation
	'''
	if "uorf" in data:
		accepted_orftypes.append("uorf")
	if "ouorf" in data:
		accepted_orftypes.append("ouorf")
	if "cds" in data:
		accepted_orftypes.append("cds")
	if "extension" in data:
		accepted_orftypes.append("extension")
	if "nested" in data:
		accepted_orftypes.append("nested")
	if "odorf" in data:
		accepted_orftypes.append("odorf")
	if "plusone" in data:
		accepted_orftypes.append("plusone")
	if "minusone" in data:
		accepted_orftypes.append("minusone")
	if "readthrough" in data:
		accepted_orftypes.append("readthrough")
	if "dorf" in data:
		accepted_orftypes.append("dorf")
	if "noncoding" in data:
		accepted_orftypes.append("noncoding")	
	'''
	
	region = data["region"]
	accepted_orftypes.append(region)
	
	
	logging.debug("accepted orftypes {}".format(accepted_orftypes))
	accepted_orftype = accepted_orftypes[0]
	
	

	try:
		cons_score = data["cons_score"].strip()
	except:
		cons_score = ""
	user_defined_transcripts = []
	#tran_list is a radio button, user can choose between principal, all or a custom list
	tranlist = data["tran_list"]
	#custom_tran_list is the actual comma seperated list of transcripts that user would enter should they choose the custom option in tranlist
	custom_tran_list = data["custom_tran_list"]
	
	ambig = False
	if "ambig_check" in data:
		ambig = True

	filtered_transcripts = {}
	if "saved_check" in data:
		#if filter previously saved cases is turned on, then we query the sqlite database here and remove hits from transcript_list
		cursor.execute("SELECT tran,stop FROM users_saved_cases WHERE user_id = '{}' and organism = '{}';".format(user_id,organism))
		result = cursor.fetchall()
		
		for tran in result:
			if str(tran[0]) not in filtered_transcripts:
				filtered_transcripts[str(tran[0])] = []
			filtered_transcripts[str(tran[0])].append(int(tran[1]))
	settings_string = "{},{},{},{},{},{},{},{},{},{},{},{},{}".format(organism,transcriptome,str(start_codons).strip("[]").replace("'",""),min_cds,max_cds,
						min_len,max_len,min_avg,max_avg,str(accepted_orftypes[0]),str(data["file_list"]).strip("[]").replace("'",""),
					  str(custom_tran_list).strip("[]").replace("'",""),cons_score)


	if "start_increase_check" in data:
		#settings_string += "{},{}".format(min_start_increase, max_start_increase)
		feature_list.append("Start value")

	if "stop_decrease_check" in data:
		#settings_string += "{},{}".format(min_stop_decrease, max_stop_decrease)
		feature_list.append("Stop value")

	if "coverage_check" in data:
		#settings_string += "{},{}".format(min_coverage, max_coverage)
		feature_list.append("Coverage value")

	if "lowest_frame_diff_check" in data:
		#settings_string += "{},{}".format(min_lowest_frame_diff, max_lowest_frame_diff)
		feature_list.append("Lowframe value")

	if "highest_frame_diff_check" in data:
		#settings_string += "{},{}".format(min_highest_frame_diff, max_highest_frame_diff)
		feature_list.append("Highframe value")
	#feature_list.append("Inframe Count Value")
	if html_args["user_short"] == "None":
		short_code = generate_short_code(data,organism,data["transcriptome"],"orf_translation")
	else:
		short_code = html_args["user_short"]
		user_short_passed = True	

	filename = short_code+".csv"
	if os.path.isfile("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename)):
		logging.debug("File exists {}/static/tmp/{}" .format(config.SCRIPT_LOC,filename))
		tmp_result_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"r")
		total_rows = 0
		returnstr = "Table|"
		for row in tmp_result_file:
			splitrow = row.split(",")
			if len(splitrow) <= 1:
				continue
			tran = splitrow[1].strip(" ")
			try:					
				stop = int(splitrow[3].strip(" "))
			except:
				continue
			if tran in filtered_transcripts:
				if stop in filtered_transcripts[tran]:
					continue

			if total_rows <1001:
				returnstr += "{},{},{},{},{},{},{},NULL,NULL,NULL,NULL,{}.,/".format(splitrow[0],splitrow[1],splitrow[2],splitrow[3],splitrow[4],splitrow[5],splitrow[6],splitrow[7])
				total_rows += 1
			else:
				break
		returnstr += "|"
		returnstr += "Min,NULL,NULL,NULL,NULL,NULL,NULL.,/"
		returnstr += "10th_percentile,NULL,NULL,NULL,NULL,NULL,NULL.,/"
		returnstr += "25th_percentile,NULL,NULL,NULL,NULL,NULL,NULL.,/"
		returnstr += "50th_percentile,NULL,NULL,NULL,NULL,NULL,NULL.,/"
		returnstr += "Max,NULL,NULL,NULL,NULL,NULL,NULL.,/"
		returnstr += "|{}".format(filename)
		returnstr += "|{}".format(str(data["file_list"]).strip("[]").replace("'","").replace(" ",""))
		returnstr += "|{}".format(data["html_args"]["user_short"])
		return returnstr
	else:
		logging.debug("File does not exists {}/static/tmp/{}" .format(config.SCRIPT_LOC,filename))
	
	if custom_tran_list != "":
		custom_tran_list  = custom_tran_list.replace(","," ")
		for item in custom_tran_list.split(" "):
			user_defined_transcripts.append(item)

	if accepted_orftypes == []:
		return "Error no ORF type selected"
	tran_gene_dict = {}
	# structure of orf dict is transcript[stop][start] = {"length":x,"score":0,"cds_cov":0} each stop can have multiple starts
	accepted_orf_dict = {}

	if start_codons == []:
		return "Error no start codon types selected:"



	if owner == 1:
		if os.path.isfile("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome)):
			traninfo_connection = sqlite3.connect("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome))
		else:
			return "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
	else:
		traninfo_connection = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))


	#traninfo_connection = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips_annotations/{0}/{0}.{1}.sqlite".format(organism,transcriptome))
	traninfo_cursor = traninfo_connection.cursor()

	principal_transcripts = []
	if tranlist == "prin_trans":
		traninfo_cursor.execute("SELECT transcript,gene FROM transcripts WHERE principal = 1;")
		result = traninfo_cursor.fetchall()
		for row in result:
			principal_transcripts.append(str(row[0]))
			tran_gene_dict[row[0]] = row[1]
	elif tranlist == "all_trans":
		traninfo_cursor.execute("SELECT transcript,gene FROM transcripts;")
		result = traninfo_cursor.fetchall()
		for row in result:
			principal_transcripts.append(str(row[0]))
			tran_gene_dict[row[0]] = row[1]
	elif tranlist == "custom_trans":
		principal_transcripts = user_defined_transcripts
		traninfo_cursor.execute("SELECT transcript,gene FROM transcripts WHERE transcript IN ({});".format(str(principal_transcripts).strip("[]")))
		result = traninfo_cursor.fetchall()
		for row in result:
			tran_gene_dict[row[0]] = row[1]
	
	transcriptome_info_dict = {}
	traninfo_cursor.execute("SELECT transcript,strand,chrom from transcripts WHERE transcript IN ({});".format(str(principal_transcripts).strip("[]")))
	result = traninfo_cursor.fetchall()
	for row in result:
		transcriptome_info_dict[str(row[0])] = {"strand":row[1],"chrom":row[2],"exons":[]}	
	traninfo_cursor.execute("SELECT * from exons WHERE transcript IN ({});".format(str(principal_transcripts).strip("[]")))
	result = traninfo_cursor.fetchall()
	for row in result:
		transcriptome_info_dict[str(row[0])]["exons"].append((row[1],row[2]))	
	logging.debug("building transcriptom info dict")
	sequence_dict = {}
	traninfo_cursor.execute("SELECT transcript,sequence FROM transcripts WHERE transcript IN ({})".format(str(principal_transcripts).strip("[]")))
	result = traninfo_cursor.fetchall()
	for row in result:
		sequence_dict[row[0]] = row[1]
	#logging.debug("sequence dict keys",sequence_dict.keys())
	#Holds a list of all transcripts in accepted_orf_dict
	accepted_transcript_list = []
	for table_name in accepted_orftypes:
		#logging.debug("table_name", table_name)
		#logging.debug("start codons", start_codons)
		if "any" in start_codons:
			logging.debug("selecting any start")
			traninfo_cursor.execute("SELECT transcript,start_codon,length,cds_coverage,start,stop  FROM {} WHERE cds_coverage >= {} AND cds_coverage <= {} AND length >= {} AND length <= {} AND transcript IN ({});".format(table_name,min_cds,max_cds,min_len, max_len, str(principal_transcripts).strip("[]")))
		else:
			logging.debug("selecting aug,cug,gug")
			#logging.debug("SELECT transcript,start_codon,length,cds_coverage,start,stop  FROM {} WHERE start_codon IN ({}) AND cds_coverage >= {} AND cds_coverage <= {} AND length >= {} AND length <= {} AND transcript IN ({});".format(table_name, str(start_codons).strip("[]"),min_cds,max_cds,min_len, max_len, str(principal_transcripts).strip("[]")))
			traninfo_cursor.execute("SELECT transcript,start_codon,length,cds_coverage,start,stop  FROM {} WHERE start_codon IN ({}) AND cds_coverage >= {} AND cds_coverage <= {} AND length >= {} AND length <= {} AND transcript IN ({});".format(table_name, str(start_codons).strip("[]"),min_cds,max_cds,min_len, max_len, str(principal_transcripts).strip("[]")))
		result = traninfo_cursor.fetchall()
		traninfo_connection.close()
		logging.debug("for row in result")
		rows = 0 
		for row in result:
			rows += 1
			if rows%1000 == 0:
				logging.debug("rows {}".format(rows))
			#logging.debug("row", row)
			transcript = str(row[0])
			start_codon = str(row[1])
			length = row[2]
			cds_cov = row[3]
			start = row[4]
			# Orfs with multiple potential starts will be grouped by start codon (so only one potential start is reported)
			# If the user has selected all transcripts then instead group by the genomic stop codon co-ordinates (so only one ORF will be reported,
			# even if it occurs on multiple transcript isoforms)
			if tranlist != "all_trans2":
				stop = row[5]
				locus = transcript
			else:
				stop = tran_to_genome(transcript,row[5],transcriptome_info_dict)
				locus = tran_gene_dict[transcript]
			
			transcriptome_stop = row[5]
			if locus not in accepted_orf_dict:
				accepted_orf_dict[locus] = {}
			if stop not in accepted_orf_dict[locus]:
				accepted_orf_dict[locus][stop] = {}
			if transcript not in accepted_transcript_list:
				accepted_transcript_list.append(transcript)
			accepted_orf_dict[locus][stop][start] = {"length":length,
														"score":0,
														"cds_cov":cds_cov,
														"start_codon":start_codon,
														"orftype":table_name,
														"stop":transcriptome_stop,
														"transcript":transcript}


	#logging.debug("accepted orf dict", accepted_orf_dict)
	logging.debug("accepted orf dict built")

	#Now build a profile for every transcript in accepted_transcripts
	master_profile_dict = {}

	all_scores = []
	all_te = []
	all_start_increases = []
	all_stop_decreases = []
	#all_cds_ratios = []
	all_results = []

	cds_average_dict = {}
	score_dict = {}

	# keeps track of the number of hits per gene
	gene_count_dict = {}
	missing_file_ids = []

	if file_paths_dict["rnaseq"] == {}:
		if "te_check" in data:
			del data["te_check"]
	
	if file_paths_dict["riboseq"] == {} and file_paths_dict["proteomics"] == {}:
		return "Error no files selected"
	



	total_files = 0
	selected_seq_types = []
	if "riboseq" in file_paths_dict:
		total_files += len(file_paths_dict["riboseq"])
		if "riboseq" not in selected_seq_types:
			selected_seq_types.append("riboseq")
	if "proteomics" in file_paths_dict:
		total_files += len(file_paths_dict["proteomics"])
		if "proteomics" not in selected_seq_types:
			selected_seq_types.append("proteomics")


	total_trans = len(accepted_transcript_list)
	logging.debug("total trans".format(total_trans))

	profile_dict,file_string = create_profiles(file_paths_dict,accepted_transcript_list,ambig,total_files,minscore)
	logging.debug("profile dict built")
	if region == "cds" or region == "noncoding":
		cds_count_dict = {}
	else:
		if owner == 1:
			if os.path.isfile("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome)):
				traninfo_connection = sqlite3.connect("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome))
			else:
				return "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
		else:
			traninfo_connection = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
		traninfo_cursor = traninfo_connection.cursor()
		traninfo_cursor.execute("SELECT transcript,cds_start,cds_stop,length  FROM transcripts WHERE principal = 1 AND tran_type = 1;")
		result = traninfo_cursor.fetchall()
		traninfo_connection.close()
		cds_dict = {}
		for row in result:
			transcript = str(row[0])
			cds_start = int(row[1])
			cds_stop = int(row[2])
			length = int(row[3])
			cds_dict[transcript] = {"cds_start":cds_start,"cds_stop":cds_stop,"length":length}
		cds_count_dict = extract_cds_counts(profile_dict,cds_dict)
	logging.debug("extracting values")
	if output != "nnet":
		sorted_all_values = extract_values(accepted_orf_dict,data,tran_gene_dict,selected_seq_types,profile_dict, cds_count_dict)
	else:
		sorted_all_values = extract_values(accepted_orf_dict,data,tran_gene_dict,selected_seq_types,profile_dict, cds_count_dict)
	if sorted_all_values == None:
		return "No results, try making filters less restrictive"
	returnstr = write_to_file(sorted_all_values,filename,sequence_dict,organism,transcriptome,file_string)
	#cursor.execute("INSERT INTO cache VALUES('{}','{}');".format(settings_string, filename))
	#connection.commit()
	#connection.close()
	
	if output == "nnet":
		#Find the top 2000 genes with the highest average riboseq in CDS (used as positive training examples)

		if owner == 1:
			if os.path.isfile("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome)):
				traninfo_connection = sqlite3.connect("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome))
			else:
				return "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
		else:
			traninfo_connection = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
		traninfo_cursor = traninfo_connection.cursor()
		traninfo_cursor.execute("SELECT transcript,cds_start,cds_stop,length  FROM transcripts WHERE principal = 1 AND tran_type = 1;")
		result = traninfo_cursor.fetchall()
		traninfo_connection.close()
		cds_dict = {}
		for row in result:
			transcript = str(row[0])
			cds_start = int(row[1])
			cds_stop = int(row[2])
			length = int(row[3])
			cds_dict[transcript] = {"cds_start":cds_start,"cds_stop":cds_stop,"count":0,"length":length}
		top_cases = get_highly_expressed_trans(cds_dict, file_paths_dict,2000)
		profile_dict = {}
		profile_dict,null_file_string = create_profiles(file_paths_dict,top_cases,ambig,total_files,minscore)
		training_accepted_orf_dict = create_training_dict(top_cases,cds_dict)
		training_sorted_all_values = extract_values(training_accepted_orf_dict,data,tran_gene_dict,selected_seq_types,profile_dict,cds_count_dict)
		training_filename = short_code+".training.csv"
		write_to_file(training_sorted_all_values,training_filename,sequence_dict,organism,transcriptome,file_string)
		returnstr = neural_net(sorted_all_values,organism, transcriptome, file_string,training_filename,feature_list,filename)
		
	
	

	logging.debug("creating returnstr")
	returnstr += "|"
	returnstr += "Min,0,0,0,0,0,0.,/"
	returnstr += "10th_percentile,0,0,0,0,0,0.,/"
	returnstr += "25th_percentile,0,0,0,0,0,0.,/"
	returnstr += "50th_percentile,0,0,0,0,0,0.,/"
	returnstr += "Max,0,0,0,0,0,0.,/"
	returnstr += "|{}".format(filename)
	returnstr += "|{}".format(str(data["file_list"]).strip("[]").replace("'","").replace(" ",""))#file_list is empty after passin through generate short_Code need to make a copy of it beforehand
	returnstr += "|{}".format(short_code)
	total_time = time.time()-start_time
	logging.debug("sending email")
	#If the job was > 5 minutes and the user is using an email address, send an email to say the job is done
	if user != None:
		if total_time > 300 and "@" in user:
			try:
				fromaddr = "ribopipe@gmail.com"
				toaddr = user
				msg = MIMEMultipart()
				msg['From'] = fromaddr
				msg['To'] = toaddr
				msg['Subject'] = "Trips-Viz job completion"
				msg.attach(MIMEText("Your Trips-Viz job is complete: https:trips.ucc.ie/short/{}".format(short_code)))
				server = smtplib.SMTP('smtp.gmail.com', 587)
				server.starttls()
				#TODO, move this to the config file
				server.login(fromaddr, "Ribosome")
				text = msg.as_string()
				logging.debug("sending now")
				server.sendmail(fromaddr, toaddr, text)
				server.quit()
			except:
				pass
	logging.debug("returning result")
	if output == "table" or output == "nnet":
		return returnstr
	elif output == "plot":
		return create_orf_plot(sorted_all_values,organism, transcriptome, file_string,filename)
	else:
		return "Undetermined output {}".format(output)




