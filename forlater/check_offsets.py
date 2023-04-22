import sys
from sqlitedict import SqliteDict
import os

indir = (sys.argv[1])

for file in os.listdir(indir):
    if file.endswith(
            ".sqlite") and "aggregate" not in file and "unmapped" not in file:
        filename = (os.path.join(indir, file))
        print("filename", filename)
        infile = SqliteDict(filename)
        offset_dict = infile["offsets"]
        if file in [
                "SRR970490.sqlite", "SRR970538.sqlite", "SRR970561.sqlite",
                "SRR970565.sqlite", "SRR970587.sqlite", "SRR970588.sqlite"
        ]:
            #offset_dict["fiveprime"]["offsets"][32] = 16
            offset_dict["fiveprime"]["offsets"][33] = 16
            offset_dict["fiveprime"]["offsets"][34] = 16
            offset_dict["fiveprime"]["offsets"][35] = 16
            for readlen in offset_dict["fiveprime"]["offsets"]:
                offset = offset_dict["fiveprime"]["offsets"][readlen]
                if offset > 10:
                    print(readlen, offset)
                else:
                    print(readlen, offset, "error")
            infile["offsets"] = offset_dict
            infile.commit()
            infile.close()
