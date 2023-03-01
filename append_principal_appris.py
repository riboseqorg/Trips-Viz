import sys
import sqlite3
from sqlitedict import SqliteDict

sqlite_db = SqliteDict(sys.argv[2])
tranlist = sqlite_db["transcripts"]

transhelve = sqlite3.connect(sys.argv[1])
cursor = transhelve.cursor()
cursor.execute("UPDATE transcripts SET principal = 0")
transhelve.commit()

cursor.execute("SELECT transcript,gene,length from transcripts")
result = cursor.fetchall()
print("tranlist done")

#print "UPDATE transcripts SET principal = 1 WHERE transcript IN ({});".format(str(tranlist).strip("[]"))
cursor.execute(
    "UPDATE transcripts SET principal = 1 WHERE transcript IN ({});".format(
        str(tranlist).strip("[]")))

transhelve.commit()
transhelve.close()
