import sys

from Bio import SeqIO
from tqdm import tqdm

present = []

with open(sys.argv[2], "w") as f:
    for record in tqdm(SeqIO.parse(sys.argv[1], "fasta")):
        rec_id = record.id.split(".")[0]
        if rec_id in present:
            continue
        present.append(rec_id)
        f.write(">" + rec_id + "\n")
        f.write(str(record.seq) + "\n")
