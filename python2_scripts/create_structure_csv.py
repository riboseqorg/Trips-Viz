import sys
import sqlite3
from sqlitedict import SqliteDict

transhelve = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips_annotations/mus_musculus/mus_musculus.v2.sqlite")
cursor = transhelve.cursor()
infile = SqliteDict("/home/DATA/www/tripsviz/tripsviz/uploads/structure_surfer_3/structure_surfer_vivo.sqlite")

cursor.execute("SELECT transcript,cds_start,cds_stop,gene from transcripts WHERE principal = 1;")
result = cursor.fetchall()

outfile = open("vivo_gene_avgs.csv","w")


for row in result:
	tran = row[0]
	cds_start = int(row[1])
	cds_stop = int(row[2])
	if cds_start <=0:
		continue
	gene = row[3]
	cds_len = cds_stop - cds_start
	shape_count = 0
	cov = 0.0
	if tran in infile:
		shape_counts = infile[tran]
		for i in range(cds_start,cds_start+300):
			if i in shape_counts["unambig"][1]:
				if shape_counts["unambig"][1][i] > 0:
					if gene == "AGL":
						print "success",i,  shape_counts["unambig"][1][i]
					cov += 1
				shape_count += shape_counts["unambig"][1][i]
	avg_shape_count = shape_count/300
	coverage = cov/300
	if gene == "AGL":
		print "cov, cds_len, coverage", cov, cds_len, coverage
	if avg_shape_count > 0:# and coverage > 0.5:
		outfile.write("{},{}\n".format(gene, avg_shape_count))


