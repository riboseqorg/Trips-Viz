import sys

gc_file = open("best_transcripts.csv_GC.csv").readlines()
hipp_file = open("dmso_vs_hipp_TE.csv").readlines()

hipp_trans = []

for line in hipp_file[1:]:
	splitline = (line.replace("\n","")).split(",")
	tran = splitline[1]
	hipp_trans.append(tran)

outfile = open("hipp_gc.csv","w")

for line in gc_file[1:]:
	splitline = (line.replace("\n","")).split(",")
	gene = splitline[0]
	tran = splitline[1]
	fivelen = splitline[3]
	gc = splitline[-1]
	#print tran, gc
	if tran in hipp_trans:
		continue
	else:
		label = "unregulated"	
	outfile.write("{},{},{},{},{}\n".format(gene,tran,label,gc,fivelen))
