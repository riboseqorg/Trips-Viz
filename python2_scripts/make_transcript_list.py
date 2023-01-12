from sqlitedict import SqliteDict


sqlite_db = SqliteDict("./trips_annotations/homo_sapiens/principal_isoforms_5ldr3tlr_rnaseq.sqlite")
anno_db = SqliteDict("./trips_annotations/homo_sapiens/homo_sapiens.sqlite")



trans = sqlite_db["transcripts"]
outfile = open("principal_transcript_ids.txt","w")

for tran in trans:
	gene = anno_db[tran]["gene"]
	outfile.write("{},{}\n".format(gene,tran))
