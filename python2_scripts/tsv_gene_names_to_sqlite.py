import sys
from sqlitedict import SqliteDict


infile = open(sys.argv[1],"r")
inlines = infile.readlines()
traninfo = SqliteDict("/home/DATA/www/tripsviz/tripsviz/trips_annotations/homo_sapiens/homo_sapiens.sqlite",autocommit=False)
gene_dict =  traninfo["genes"]

master_dict = {"unambiguous_fiveprime_totals":{},
				"unambiguous_cds_totals":{},
				"unambiguous_threeprime_totals":{},
				"unambiguous_all_totals":{}}


def nuc_to_amino(nucseq):
	seq_dict = {0:"",1:"",2:""}
	codon_dict = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
	"TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
	"TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
	"TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
	"CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
	"CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
	"CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
	"CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
	"ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
	"ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
	"AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
	"AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
	"GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
	"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
	"GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
	"GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

	for y in range(0,2):
		for i in range(y,len(nucseq),3):
			codon = nucseq[i:i+3]
			seq_dict[y].append(codon_dict[codon])
	return seq_dict






for line in inlines[1:]:
	#print line
	splitline = line.split("\t")
	peptide_seq = splitline[0]
	gene = splitline[28].upper().strip()
	print gene
	if gene not in gene_dict:
		continue
	tran_list = gene_dict[gene]
	print tran_list
	for tran in tran_list:
		nucseq = traninfo[tran]["seq"]
		aminoseq = nuc_to_amino(nucseq)
		print tran, aminoseq
		sys.exit()
		#start_pos = int(splitline[36])
		#end_pos = int(splitline[37])
		#readlength = (end_pos-start_pos)*3
		#if transcript not in master_dict:
		#	master_dict[transcript] = {"unambig":{}, "ambig":{}}
		#if readlength not in master_dict[transcript]["unambig"]:
		#	master_dict[transcript]["unambig"][readlength] = {}
		#if start_pos not in master_dict[transcript]["unambig"][readlength]:
		#	master_dict[transcript]["unambig"][readlength][start_pos] = 0
		#master_dict[transcript]["unambig"][readlength][start_pos] += 1





for transcript in master_dict:
	if transcript in traninfo:
		all_total = 0
		fiveprime_total = 0
		cds_total = 0
		threeprime_total = 0
		cds_start = traninfo[transcript]["cds_start"]
		cds_stop = traninfo[transcript]["cds_stop"]
		for readlength in master_dict[transcript]["unambig"]:
			for pos in master_dict[transcript]["unambig"][readlength]:
				count = master_dict[transcript]["unambig"][readlength][pos]
				all_total += count
				if pos <= cds_start:
					fiveprime_total += count
				elif pos > cds_start and pos <= cds_stop:
					cds_total += count
				elif pos > cds_stop :
					threeprime_total += count
		master_dict["unambiguous_fiveprime_totals"][transcript] = fiveprime_total
		master_dict["unambiguous_cds_totals"][transcript] = cds_total
		master_dict["unambiguous_threeprime_totals"][transcript] = threeprime_total
		master_dict["unambiguous_all_totals"][transcript] = all_total


output = SqliteDict("proteomics.sqlite",autocommit=False)
for key in master_dict:
	output[key] = master_dict[key]
output.commit()
output.close()
