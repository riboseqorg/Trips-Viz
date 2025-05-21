#! /usr/bin/env python

import sqlite3

import click
import polars as pl
from sqlitedict import SqliteDict
from tqdm import tqdm


@click.command()
@click.argument('annotations', type=click.Path(exists=True))
def run(annotations):



    annotations_helper = SqliteDict(annotations.replace('.sqlite','_helper.sqlite'))

    annotations = pl.read_database(query="select * from transcripts",connection = sqlite3.connect(annotations)) 

    for row in tqdm(annotations.iter_rows(named=True)):
        starts = {0:[],1:[],2:[]}
        orfs = []
        for i in range(len(row["sequence"])-3):
            if row["sequence"][i:i+3] == "ATG": 
                starts[i%3].append(i)
            if row["sequence"][i:i+3] in ["TAG","TGA" ,"TAA"]:
                for start in starts[i%3]:
                    orfs.append([i%3,start,i])
                    break
                starts[i%3] = []
        start_stop_list = []
        # start_stop_list = [[i,-5,"start"] for i in range(3)]
        
        for st in row["start_list"].split(","):
            if '-' in st: continue
            stt = int(st) 
            start_stop_list.append([(stt%3+2)%3,stt,"start"])# TODO: Check why shifts are required
        for st in row["stop_list"].split(","):
            if '-' in st: continue
            stt = int(st) 
            start_stop_list.append([((stt+1)%3+2)%3 ,stt,"stop"]) # TODO: Check why shifts are required
        annotations_helper[row["transcript"]] = {"start_stop_list":start_stop_list,"start_stop_list_cols":["frame","position","type"],"orfs":orfs, 'orf_cols':['frame','start','stop']}

    annotations_helper.commit()
    annotations_helper.close()



if __name__ == "__main__":
    run()
