
#gene_data = open("../gene_details/gene_details_human")
#gene_type = {}
#gene_type2 = {}
#cds_set = []
#for line in gene_data :
	#linesplit = line[:-1].split("\t")
	#gene_type[linesplit[0]] = linesplit[5]

	#if linesplit[5] == "CDS" :
		#cds_set.append(linesplit[0])
	##gene_type2.setdefault(linesplit[0],[]).append(linesplit[5])
#gene_data.close()


#gene_data = open("../gene_details/gene_details_mouse")
#gene_type = {}
#gene_type2 = {}
##cds_set = []
#for line in gene_data :
	#linesplit = line[:-1].split("\t")
	#gene_type[linesplit[0]] = linesplit[5]

	#if linesplit[5] == "CDS" :
		#cds_set.append(linesplit[0])
	##gene_type2.setdefault(linesplit[0],[]).append(linesplit[5])
#gene_data.close()
#cds_set = set(cds_set)


import os
files = os.listdir("z_combine_ribo")
for filea in files:
	if filea in ["human","mouse"]: continue
	print filea
	infileopen = open("z_combine_ribo/%s"%filea)
	infileopen.readline()
	outfileopen = open("gene_reg_tables/%s.csv"%filea,"w")
	for line in infileopen:
		linesplit = line.split("\t")
		#print linesplit
		#if linesplit[0]
		outfileopen.write("%s\t%s"%(linesplit[1], linesplit[2]))
	outfileopen.close()
	infileopen.close()
	