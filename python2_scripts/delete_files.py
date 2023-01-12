import sqlite3
import time
import sys

input_folder = sys.argv[1]

del_file = open("{}/deletions.txt".format(input_folder), "a")


connection = sqlite3.connect('{}/trips.sqlite'.format(input_folder))
connection.text_factory = str
cursor = connection.cursor()
cursor.execute("SELECT * from deletions;")
result = cursor.fetchall()



for row in result:
	file_id = row[0]
	file_path = row[1]
	del_time = row[2]
	curr_time = time.time()
	#print "file id, file_path, del_time, curr_time", file_id, file_path, del_time, curr_time
	if curr_time > del_time:
		del_file.write("cursor.execute(DELETE FROM files WHERE file_id = {})\n".format(file_id))			
		del_file.write("os.remove('{}')\n".format(file_path))
		del_file.write("cursor.execute(DELETE from deletions WHERE file_id = '{}')\n".format(file_id))


cursor.execute("SELECT * from org_deletions;")
result = cursor.fetchall()



for row in result:
        org_id = row[0]
        file_path = row[1]
        del_time = row[2]
        curr_time = time.time()
        #print "file id, file_path, del_time, curr_time", file_id, file_path, del_time, curr_time
        if curr_time > del_time:
                del_file.write("cursor.execute(DELETE FROM organisms WHERE organism_id = {})\n".format(org_id))
                del_file.write("os.remove('{}')\n".format(file_path))
                del_file.write("cursor.execute(DELETE from org_deletions WHERE organism_id = '{}')\n".format(org_id))

