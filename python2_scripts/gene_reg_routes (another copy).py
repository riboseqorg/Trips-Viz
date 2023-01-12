from flask import Blueprint, render_template, abort, request
import sqlite3
import ast
import os
import numpy as np
from scipy.stats import  mannwhitneyu,percentileofscore
#from scipy.stats import 

gene_regulation_page = Blueprint("gene_regulation_page", __name__, template_folder="templates2")
@gene_regulation_page.route('/<organism>/<transcriptome>/gene_regulation/')
def gene_regulationpage(organism,transcriptome):
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

	connection = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips.sqlite")
	cursor = connection.cursor()
	advanced = False
	if user != None:
		cursor.execute("SELECT advanced from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		print "\n\n\n\nAdvanced {}".format(result)
		if result[0] == 1:
			advanced = True
		else:
			advanced = False



	return render_template('gene_regulation.html', studies_dict={}, accepted_files={},user=user,
						   organism=organism,default_tran="",local=local,transcriptome=transcriptome,advanced=advanced,
						   seq_types={},studyinfo_dict={})


gene_regulation_query = Blueprint("gene_regulation_query", __name__, template_folder="templates")
@gene_regulation_query.route('/gene_regulation_query', methods=['POST'])
def generegulationquery():
	print "gene regulation query called"
	connection = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips.sqlite")
	cursor = connection.cursor()
	script_location = "/home/DATA/www/tripsviz/tripsviz"
	acceptable = 0
	data = ast.literal_eval(request.data)

	try:
		user = current_user.name
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]
		#cursor.execute("SELECT user_id from users WHERE username = 'master';")
		#result = (cursor.fetchone())
	except:
		user = None

	organism = data["organism"]
	returnstr = ""
	secondary_tablestr = ""
	secondary_tablestr += "{}:{}.,/".format(5,"ATF4")
	secondary_tablestr += "{}:{}.,/".format(6,"ATF4,ATF5")
	secondary_tablestr += "{}:{}.,/".format(7,"ATF4,ATF5,SLC35A4")
	transcriptome = data["transcriptome"]
	gene_list = (data["tran_list"].replace(" ","")).split(",")


			
	gene_list = [gxx.upper() for gxx in gene_list]
	if len(gene_list) == 0:
		#break
		pass
	else :
		#if len(gene_list) == 1:
		if 1:
			infileopen = open("pmid/dip_pair_top")
			start_dict = {}
			for line in infileopen:
				linesplit = line[:-1].split("\t")
				#print linesplit
				if linesplit[0] == linesplit[1]:
					start_dict[linesplit[0]] = float(linesplit[2])

			#infileopen = open("cluster_euclidean_average2b")
			dicta = {}
			dictb = {}
			dictc = {}
			dictd = {}
			used = []
			#infileopen = open("pmid/cluster_euclidean_average2b")
			##show = []
			##for line in infileopen:
				##if line[-1] != "\n": continue
				##linesplit = line[:-1].split(",")
				##list1 = linesplit[2].split("\t")
				##list1 = [nxx for nxx in list1 if nxx != ""]
				##list2 = linesplit[3].split("\t")
				##list2 = [nxx for nxx in list2 if nxx != ""]
				###list3.extend(linesplit[3].split("\t"))

				##list3 = linesplit[2].split("\t")
				##list3.extend(linesplit[3].split("\t"))
				##list3 = [nxx for nxx in list3 if nxx != ""]
				###stb = "".join([for nxx in list3])
				##stb = " ".join(list3)
				###listdd.append((float(linesplit[0]), list3))
				##if set(gene_list) - set(list3)  == set([]):
					##if float(linesplit[0]) > 3:
						###show = [linesplit[0], stb]

			#geen = gene_list[0]
			#dicta = {}
			#dictb = {}
			#dictc = {}
			#for line in infileopen:
				#linesplit = line[:-1].split(",")
				##print linesplit[2], 1, linesplit[3], 2
				##print
				#if line[-1] != "\n": continue
				#list1 = linesplit[2].split("\t")
				#list1 = [nxx for nxx in list1 if nxx != ""]
				#list2 = linesplit[3].split("\t")
				#list2 = [nxx for nxx in list2 if nxx != ""]
				##list3.extend(linesplit[3].split("\t"))

				#list3 = linesplit[2].split("\t")
				#list3.extend(linesplit[3].split("\t"))
				#list3 = [nxx for nxx in list3 if nxx != ""]
				##for key in list3:
					##dictb[key] = float(linesplit[0])
				##print list1, list2
				#if len(list1) == 1 and len(list2):
					#dicta[str(list3)] = "(%s,%s)"%(list1[0],list2[0])
					#dictb[str(list2)] = float(linesplit[0])
				#elif len(list1) ==1:
					#dictb[str(list1)] = float(linesplit[0])
					#dicta[str(list3)] = "(%s,%s)"%(list1[0],dicta[str(list2)])
				#elif len(list2) ==1:
					#dicta[str(list3)] = "(%s,%s)"%(dicta[str(list1)],list2[0])
				#else:
					#dicta[str(list3)] = "(%s,%s)"%(dicta[str(list1)],dicta[str(list2)])

				##dictb[str(list1)] = float(linesplit[0])
				##dictb[str(list2)] = float(linesplit[0])
				##dictb[str(list3)] =
				##print dictb
				#dictc[str(list3)] = float(linesplit[0])
				#if len(list1) != 1:
					#dictd[str(list1)] = dictc[str(list1)]- float(linesplit[0])
				#else:
					#dictd[str(list1)] = start_dict[list1[0]]- float(linesplit[0])

				#if len(list2) != 1:
					#dictd[str(list2)] = dictc[str(list2)]- float(linesplit[0])
				#else:
					#dictd[str(list2)] = start_dict[list2[0]]- float(linesplit[0])

				##print  dictb[str(list1)]- float(linesplit[0])
				##if len(list1) == 1 and len(list2):
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(list1[0],float(linesplit[0]),list2[0],float(linesplit[0]))
				##elif len(list1) ==1:
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(list1[0],float(linesplit[0]),dicta[str(list2),float(linesplit[0])])
				##elif len(list2) ==1:
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(dicta[str(list1)],float(linesplit[0]),list2[0],float(linesplit[0]))
				##else:
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(dicta[str(list1)],float(linesplit[0]),dicta[str(list2),float(linesplit[0])])

				#if geen in list3:
					#if len(list3) > 500 :break

					##print dicta[str(list3)]
				#used.append(list3)
				#try:
					#used.remove(list1)
				#except:
					#continue

				#try:
					#used.remove(list2)
				#except:
					#continue


			#dicta1 = {}
			#infileopen.seek(0)
			#for line in infileopen:
				#linesplit = line[:-1].split(",")
				#if line[-1] != "\n": continue
				##print linesplit[2], 1, linesplit[3], 2
				##print
				#list1 = linesplit[2].split("\t")
				#list1 = [nxx for nxx in list1 if nxx != ""]
				#list2 = linesplit[3].split("\t")
				#list2 = [nxx for nxx in list2 if nxx != ""]
				##list3.extend(linesplit[3].split("\t"))

				#list3 = linesplit[2].split("\t")
				#list3.extend(linesplit[3].split("\t"))
				#list3 = [nxx for nxx in list3 if nxx != ""]
				##for key in list3:
					##dictb[key] = float(linesplit[0])
				##print list1, list2
				#if len(list1) == 1 and len(list2):
					#dicta1[str(list3)] = "(%s:%s,%s:%s)"%(list1[0],dictd[str(list1)],list2[0],dictd[str(list2)])
					##dictb[str(list2)] = float(linesplit[0])
				#elif len(list1) ==1:
					##dictb[str(list1)] = float(linesplit[0])
					#dicta1[str(list3)] = "(%s:%s,%s:%s)"%(list1[0],dictd[str(list1)],dicta1[str(list2)],dictd[str(list2)])
				#elif len(list2) ==1:
					#dicta1[str(list3)] = "(%s:%s,%s:%s)"%(dicta1[str(list1)],dictd[str(list1)],list2[0],dictd[str(list2)])
				#else:
					#dicta1[str(list3)] = "(%s:%s,%s:%s)"%(dicta1[str(list1)],dictd[str(list1)],dicta1[str(list2)],dictd[str(list2)])

				##dictb[str(list1)] = float(linesplit[0])
				##dictb[str(list2)] = float(linesplit[0])
				##dictb[str(list3)] =
				##print dictb
				##dictc[str(list3)] = float(linesplit[0])
				##if len(list1) != 1:
					##dictd[str(list1)] = dictc[str(list1)]- float(linesplit[0])
				##else:
					##dictd[str(list1)] = 0

				##if len(list2) != 1:
					##dictd[str(list2)] = dictc[str(list2)]- float(linesplit[0])
				##else:
					##dictd[str(list2)] = 0

				##print  dictb[str(list1)]- float(linesplit[0])
				##if len(list1) == 1 and len(list2):
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(list1[0],float(linesplit[0]),list2[0],float(linesplit[0]))
				##elif len(list1) ==1:
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(list1[0],float(linesplit[0]),dicta[str(list2),float(linesplit[0])])
				##elif len(list2) ==1:
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(dicta[str(list1)],float(linesplit[0]),list2[0],float(linesplit[0]))
				##else:
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(dicta[str(list1)],float(linesplit[0]),dicta[str(list2),float(linesplit[0])])

				#if geen in list3:
					#if len(list3) > 500 :break

					##print dicta[str(list3)]
				#used.append(list3)
				#try:
					#used.remove(list1)
				#except:
					#continue

				#try:
					#used.remove(list2)
				#except:
					#continue

				##used.remove(list2)
				##else:
			#from ete3 import Tree
			##from ete3 import Tree, TreeStyle
			##circular_style = TreeStyle()
			##circular_style.mode = "c" # draw tree in circular mode
			##circular_style.scale = 20
			#for x in  used:
				#if geen in x:
					#xp = dicta1[str( x)]
					###print xp, "xp"
					#t = Tree("%s;"%xp, format=1)
					###t = Tree("((D:0.723274,F:0.567784)1.000000:0.067192,(B:0.279326,H:0.756049)1.000000:0.807788);",format = 0)
					###print xp
					###print t
					##returnstr += ",{},,,,,,.,/".format(xp)

			#if gene_list[0] in list3:
			#from ete3 import Tree
			#if show != []:
				#returnstr += "{},{},,,,,,.,/".format(show[0], show[1])
					
		infileopen = open("pmid/pmid")
		infileopen.readline()
		dict_pmid = {}
		for line in infileopen:
			linesplit = line[:-1].split("\t")
			dict_pmid[linesplit[0]] = linesplit[1]

		infileopen = open("pmid/description")
		#infileopen.readline()
		dict_des = {}
		for line in infileopen:
			linesplit = line[:-1].split("\t")
			dict_des[linesplit[0]] = linesplit[0]
			if len(linesplit)> 1:
				dict_des[linesplit[0]] = linesplit[1]
		infileopen = open("pmid/new_name")
		#infileopen.readline()
		dict_name = {}
		for line in infileopen:
			linesplit = line[:-1].split("\t")
			dict_name[linesplit[0]] = linesplit[0]
			if len(linesplit)> 1:
				dict_name[linesplit[0]] = linesplit[1]

			#else:

		dict_pathway = {}
		#list_file2 =["mTOR_inhibition"]
		list_file2 = os.listdir("cluster/")
		list_out = []
		for filea in list_file2 :
			infileopen = open("cluster/%s"%filea)
			line = infileopen.readline()
			linesplit = line[:-1].split("\t")
			#dict_pathway[filea] =
			for key in linesplit[1:]:
				try:
					aba = float(key)
				except:
					dict_pathway.setdefault(key,[])
					if filea not in dict_pathway[key]:
							dict_pathway[key].append((0, filea))
							dict_pathway[key].sort()
							dict_pathway[key].reverse()
					#if dict_pathway[key] == "":
						#dict_pathway[key] = "{}".format(filea)
					#else:
						#dict_pathway[key] ="{}_{}".format(dict_pathway[key],filea)
			count = 0
			count1 = ""
			#for line in infileopen:
				#linesplit = line[:-1].split("\t")
				#if linesplit[0] in gene_list:
					##quer2.append(abs(float(linesplit2[1])))
					#count += 1
					##count1.append(linesplit[0])
					#if count1 == "":
						#count1 = linesplit[0]
					#else:
						#count1 = "{}_{}".format(count1, linesplit[0])

			##list_out.append((count, filea,count1))
		#list_out.sort()
		#list_out.reverse()
		#for z_score in list_out:
			#if z_score[0]> 0:
				#returnstr += "{}{},{},None,None.,/".format(z_score[0],z_score[1], z_score[2])
	#returnstr += "None,None,None,None.,/"
		#infileopen = open("patrick_regulated_in")
		#for line in infileopen:
			#linesplit = line.split("\t")
			#if gene_list[0] == linesplit[0]:
				#data = linesplit[1:]
				#data1 = []
				#data2 = []
				#for x_i, x in enumerate(data):
					#if x_i %2 ==1:
						#data1.append(x)
					#else:
						#data2.append(x)
						
				#break
		#for study, z_score in zip(data1,data2):
			##returnstr += "{},{},None,None.,/".format(z_score, study)
			#stra = ""
			#dict_pathway.setdefault(study,[])
			#for tuple1 in dict_pathway[study]:
				#if stra == "":
					#stra = "{}".format(tuple1[1])
				#else:
					#stra = "{}_{}".format(stra, tuple1[1])
			#annotation_popup = "Annotation description here_top"
			#returnstr += "{},{},None,{},{}.,/".format(study, z_score,stra,annotation_popup)
	#else:
	#if 1:
		list_file2 = os.listdir("z_combine/")
		#list_file2 = ["Alvarez-Dominguez_24_v0hr_erythroid_differentiation", "Alvarez-Dominguez_33_v24hr_erythroid_differentiation", "Alvarez-Dominguez_48_v33hr_erythroid_differentiation", "Arango18_NAT10_ko", "Aviner17_hnRNP_c_kd", "Bennett_msi2_shRNA", "bercovich_kinori_influenza_2_4hr", "bercovich_kinori_influenza_2hr", "bercovich_kinori_influenza_4_8hr", "Blair_hesc_npc", "Blair_neu14_hesc", "Blair_neu14_npc", "Blanco16_Nsun1_ko", "castaneda_mael_ko", "Castelo17_zt0", "Castelo17_zt2", "Castelo17_zt4", "Castelo17_zt6", "Castelo17_zt8", "Castelo17_zt10", "Castelo17_zt12", "Castelo17_zt14", "Castelo17_zt16", "Castelo17_zt18", "Castelo17_zt20", "Castelo17_zt22", "dai_2hr_vaccinia_infection", "dai_2v4hr_vaccinia_infection", "dai_4v8hr_vaccinia_infection", "Diaz-Munoz_Hur_KO", "Diaz-Munoz_Hur_KO_LPS", "dima_aicar_f", "dimaelife_arsenite", "dimaGB_ogd20", "dimaGB_ogd20_40", "dimaGB_ogd40_60", "dima_mechano_h1", "dima_mechano_I1", "dima_myxo_anoxia", "dima_myxo_E1", "dima_myxo_g5_hct116sca", "dima_myxo_G1_hct116", "
#dima_myxo_H1", "dima_phendc3h1", "dima_phendc3k", "dima_pieri_f", "dima_pieri_g5_hct116sca", "dima_pieri_G1_hct116", "dima_pieri_H1", "dima_te_myxo_sco2_anoxia", "Eichhorn_mi155ko_2h", "Eichhorn_mi155ko_4h", "Eichhorn_miR-1", "fijalkowska_eif1_kd", "Fradejas_Secisbp2_ko", "Fradejas_Trsp_ko", "Gameiro_EtOH_NoBCAA", "Gameiro_EtOH_NoCys", "Gameiro_EtOH_NoG", "Gameiro_EtOH_NoQ", "Gameiro_TAM_NoBCAA", "Gameiro_TAM_NoCys", "Gameiro_TAM_NoG", "Gameiro_TAM_NoQ", "Gameiro_TAM_torin", "Gao17_liver_fasting", "Gao17_MEF_starvation", "gao_aminoastarvation", "gao_eif2phos_mimetic", "Ginossar12_HCMV_24v5hr", "Gonzalez14_Tumor_v_normal", "goodarzi16_argCCG_oe", "goodarzi16_gluUUC_oe", "Goodarzi_highly_v_lowly_metastatic", "Grabole_AZD-8055_tsc2-cellline", "Grabole_AZD-8055_tsc2+cellline", "Grow_Rec_OE", "Guo14_te_miR-1_v_control", "howard13_A37G_0.1_0_selenocysteine", "howard13_A37G_2_0.1_selenocysteine", "howard13_wt_0.1_0_selenocysteine", "howard13_wt_2_0.1_selenocysteine", "hsieh_pp242", "Irigoyen_coronavirus_1h", "
##Irigoyen_coronavirus_8h", "iwasaki_roca", "Jackson18_LPS_simulation", "Jakobson_KO_EEF1AKMT4", "jan_0", "jan_2", "jan_4", "jan_6", "jan_8", "jan_10", "jan_12", "jan_14", "jan_16", "jan_18", "jan_20", "jan_22", "Janich15_zt0", "Janich15_zt2", "Janich15_zt4", "Janich15_zt6", "Janich15_zt8", "Janich15_zt10", "Janich15_zt12", "Janich15_zt14", "Janich15_zt16", "Janich15_zt18", "Janich15_zt20", "Janich15_zt22", "Ji15_Src_induced_transformation_1h", "Ji15_Src_induced_transformation_4h", "Jiang17_hypoxia1h", "khajuria_rps5", "khajuria_rps19", "Laguesse_Elp3cKO", "Leshkowitz19_eIF1A_sil", "Li_cyscys_deprivation", "Li_glucose_deprivation", "Liu18_FMRP_ko", "Loayza-Puch13_Nutlin-3a_2hv0", "Loayza-Puch13_Nutlin-3a_4hv2h", "Loayza-Puch13_Nutlin-3a_6hv4h", "Loayza-Puch13_Nutlin-3a_6hv19h", "Loayza-Puch13_pre-senescence_senescence", "Loayza-Puch13_profileration_quiescence", "Loayza-Puch14_myc_v_control", "Loayza-Puch16_tumour_v_control", "Loayza-Puch_siSL3A2", "Loayza-Puch_TGFb", "Murat_dhx9_kd", "oh_arsenite_ddx3", "oh_
#arsenite_ddx3_r534h", "oh_arsenite_wt", "Paolini_tunicamycin", "park16_s_v_m", "park_torin24", "Razooky17_RP8_v_mock", "reid_0_6h_denv1", "reid_0_6h_denv2", "reid_6_12h_denv1", "reid_6_12h_denv2", "reid_12_24h_denv1", "reid_12_24h_denv2", "reid_24_48h_denv1", "reid_24_48h_denv2", "reid_48_72h_denv1", "reid_48_72h_denv2", "Reid13_cyto_Thapsigargin_0.5_v_0hr", "Reid13_cyto_Thapsigargin_1_v_0.5hr", "Reid13_cyto_Thapsigargin_2_v_1hr", "Reid13_cyto_Thapsigargin_4_v_2hr", "Reid13_er_Thapsigargin_0.5_v_0hr", "Reid13_er_Thapsigargin_1_v_0.5hr", "Reid13_er_Thapsigargin_2_v_1hr", "Reid13_er_Thapsigargin_4_v_2hr", "reid_cytosol_v_ER", "Ricci13_Stau1_overexpression_v_KnockDown", "rubio_silvestrol_2h", "Rutkowski_hsv1_2h", "Rutkowski_hsv1_4h", "Sendoel_SOX2_induction", "Shalgi_2hr_severe_heatshock", "Shalgi_8hr_mild_heatshock", "shi_siYTHDF3", "sid_tunamycin", "Sims14_tumour_normal", "Simsek17_siPKM", "Slobodin17_Campthotecin_v_control", "Slobodin17_Nutlin", "Stumpf_G1_v_S", "Stumpf_M_v_G1", "Stumpf_S_v_M", "Su15_IFNg", "#tanenbaum_g2_g1", "tanenbaum_g2_m", "tanenbaum_m_g1", "Thoreen_torin1", "Tichon_NORD_si", "Tirosh_HCMV_5h", "Tirosh_HCMV_12h", "Tirosh_HCMV_24h", "Tirosh_irradiated_HCMV", "Tirosh_type1_interferon", "wang_siYTHDF2", "Wein14_mut_wt", "werner_hES_shKBTBD8", "werner_ND1_shKBTBD8", "werner_ND3_shKBTBD8", "werner_ND6_shKBTBD8", "witta_bortezomib1.5_3h", "witta_bortezomib1.5h", "witta_bortezomib3_6h", "witta_bortezomib6_9h", "witta_bortezomib9_12h", "wolfe_silvestrol", "xu_correct", "xu_mutant", "you15_siKrr1", "zang_mi_oe", "Zhou15_heatshock", "Zhou15_heatshock_FTO", "Zhou15_heatshock_YTHDF2", "Zur_M_G1"]
		list_file = []
		for filea in list_file2 :
				infileopen = open("z_combine_rna/%s"%filea)
				#dict2[filea] = 0
				dicta = {}
				quer1, quer2 = [],[]
				for line2 in infileopen:
					linesplit2 = line2.split("\t")
					quer1.append(abs(float(linesplit2[1])))
					if linesplit2[0] in gene_list:
						quer2.append(abs(float(linesplit2[1])))
				rna_score = mannwhitneyu(quer1, quer2)[1]

				infileopen = open("z_combine_ribo/%s"%filea)
				#dict2[filea] = 0
				dicta = {}
				quer1, quer2 = [],[]
				for line2 in infileopen:
					linesplit2 = line2.split("\t")
					quer1.append(abs(float(linesplit2[1])))
					if linesplit2[0] in gene_list:
						quer2.append(abs(float(linesplit2[1])))
				ribo_score = mannwhitneyu(quer1, quer2)[1]
				
				infileopen = open("z_combine/%s"%filea)
				#dict2[filea] = 0
				dicta = {}
				quer1, quer2 = [],[]
				change2b = [0,0, len(gene_list)]
				for line2 in infileopen:
					linesplit2 = line2.split("\t")
					quer1.append(abs(float(linesplit2[1])))
					if linesplit2[0] in gene_list:
						quer2.append(abs(float(linesplit2[1])))
						if float(linesplit2[1]) >0:
							change2b[0] += 1
						else:
							change2b[1] += 1
						#if float(linesplit2[1]) > 0:
							#dict2[filea]  += 1
						#else:
							#dict2[filea]  += -1
				if len(quer2) > 0:
					#if np.mean(quer2) > 1:
					dop= max([ribo_score,mannwhitneyu(quer1, quer2)[1] ])
					#dop = mannwhitneyu(quer1, quer2)[1]
					quer1.sort()
					quer2a = quer1[:len(quer2)]
					quer2b = quer1[-len(quer2):]
					ideal_score = min([mannwhitneyu(quer1, quer2a)[1],mannwhitneyu(quer1, quer2b)[1]])
					#if dop < rna_score:
					if 1:
						
						list_file.append((dop,  filea,round(percentileofscore(quer1, abs(np.mean(quer2))),2),ideal_score,len(quer2),change2b))
				#quer1.sort()
				#list_file_200[filea] = quer1[-1000]
		list_file.sort()

		returnstr += ",{},{},,,,.,/".format( "#Min Pvalue#",list_file[0][3])
		#returnstr += ",{},{},,.,/".format(10**(np.log10(list_file[0][-1])*4/5), "#STRONG SCORE#")

		for study in list_file:
			if study[0] > 0.05: continue
			if study[2] < 70: continue
			#if np.log10(study[0])< np.log10(study[3])*1/2:
			if 1:

				
				dict_pathway.setdefault(study[1],[])
				#dict_pathway2.setdefault(study[1],"")
				stra = ""
				for tuple1 in dict_pathway[study[1]]:
					if stra == "":
						stra = "{}".format(tuple1[1])
					else:
						stra = "{}_{}".format(stra, tuple1[1])
				dict_pmid.setdefault(study[1].split("_")[0], "na")
				annotation_popup = dict_des[study[1]]
				#returnstr += "{},{},{},{},{},.,/".format(study[1], study[0],dict_pmid[study[1].split("_")[0]],stra,annotation_popup)
				missing = study[5][2] - study[5][1]-study[5][0]
				changestr = "{}_{}_{}".format(study[5][0],study[5][1], missing)
				
				#returnstr += "{},{},{},{},{},{},{}.,/".format(dict_name[study[1]], dict_pmid[study[1].split("_")[0]],study[0],study[2],changestr,stra,annotation_popup)
				returnstr += "{},{},{},{},{},{},{},{}.,/".format(dict_name[study[1]], dict_pmid[study[1].split("_")[0]],study[0],study[2],changestr,stra,annotation_popup,"")

				print "return str noW", returnstr
	#print returnstr
	returnstr += "???"+secondary_tablestr
	return returnstr