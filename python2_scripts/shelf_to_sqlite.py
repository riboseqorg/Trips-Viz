import shelve
from sqlitedict import SqliteDict
import sys
import os


infile = shelve.open(sys.argv[1])
if infile == {}:
    sys.exit()
if os.path.isfile(sys.argv[1].replace(".shelf",".sqlite").replace(".shelve",".sqlite")):
    print "file exists, exiting"
    sys.exit()



outfile = SqliteDict(sys.argv[1].replace(".shelf",".sqlite").replace(".shelve",".sqlite"),autocommit=False)


for key in infile:
    print key
    if key == "features":
        continue 
    if key == "ENST00000282561":
        print infile["ENST00000282561"]
    outfile[key] = infile[key]


outfile.commit()
outfile.close()
infile.close()


