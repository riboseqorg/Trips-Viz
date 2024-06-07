#!/usr/bin/env python

# This script is to compress the files to reduce download size of the files

import os
import sqlite3
import config


def compress():
    trips_annotation_dir = "{}/{}/".format(config.SCRIPT_LOC,
                                           config.ANNOTATION_DIR)
    organism_dict = {
        "Scripts": [
            "bam_to_sqlite.py", "tsv_to_sqlite.py",
            "create_annotation_sqlite.py",
            "create_transcriptomic_to_genomic_sqlite.py"
        ]
    }
    connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,
                                                config.DATABASE_NAME))
    connection.text_factory = str
    cursor = connection.cursor()
    cursor.execute("SELECT organism_name from organisms where private = 0")
    result = cursor.fetchall()
    for row in result:
        organism = row[0]
        organism_dict[organism] = []
    for org in os.listdir(trips_annotation_dir):
        if org not in organism_dict:
            continue
        for filename in os.listdir(trips_annotation_dir + "/" + org):
            if "." in filename:
                ext = filename.split(".")[-1]
                if ext == "fa" or ext == "gtf":
                    organism_dict[org].append(filename)
                elif ext == "sqlite":
                    if "transcriptomic" in filename or org in filename:
                        organism_dict[org].append(filename)

    for key, values in organism_dict.items():
        for value in values:
            fl = f"{trips_annotation_dir}{key}/{value}"
            print(f'Compressing {fl}')
            os.system(f"tar -czvf {fl}.tar.gz {fl}")


if __name__ == '__main__':
    compress()
