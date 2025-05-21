import collections
from time import sleep
from typing import Dict, Tuple, Union

import pandas as pd
import polars as pl
from flask import flash
from sqlitedict import SqliteDict

from core_functions import dict2df
from fixed_values import merge_dicts


# Merge two dictionaries
def merge_dicts(
        dict1,
        dict2):  # NOTE: expecting that second dictionary is always smaller
    for readlen in dict2:
        if readlen not in dict1:
            dict1[readlen] = dict2[readlen]
        else:
            for pos in dict2[readlen]:
                if pos in dict1[readlen]:
                    dict1[readlen][pos] += dict2[readlen][pos]
                else:
                    dict1[readlen][pos] = dict2[readlen][pos]
    return dict1


# Create dictionary of read counts at each position in a transcript
def get_reads(data) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """

    Parameters: 

    Returns:

    Example:
    """
    mismatch_dict = []
    master_dict = []

    master_file_dict = {}

    # first make a master dict consisting of all the read dicts from each filename
    offset_dict = {}
    print(data['file_paths_dict']['path'])
    print("====================")
    for fl in data['file_paths_dict']['path']:
        try:
            sqlite_db = SqliteDict(fl)
        except FileNotFoundError:
            flash(f"Sqlite file for {path.split(fl)[1].split('.')[0]} not found.")
            # return pd.DataFrame(), pd.DataFrame()
        try:
            all_offsets_n_scores = sqlite_db["offsets"][ data["offsite"]]
            print("all_offsets_n_scores", all_offsets_n_scores)

        except KeyError:
            read_length = list(range(data["minread"], data["maxread"] + 1))
            range_len = data["maxread"] - data["minread"] + 1
            all_offsets_n_scores = pl.DataFrame({
                'read_len': read_length,
                'offsets': [15] * range_len,
                'read_scores': [1] * range_len
            })
        print(data)
        
        offset_dict[fl] = all_offsets_n_scores.filter(
            pl.col('read_scores') >= data["readscore"])

        if "mismatch" in data:
            try:
                mismatch_dict = pl.DataFrame(sqlite_db[data['transcript']]["seq"], schema=['pos', 'A', 'C', 'G', 'T']).with_columns(pos=pl.col('pos') + 1)

            except Exception:
                mismatch_dict = pl.DataFrame(schema=['pos', 'A', 'C', 'G', 'T'])

        try:
            alltrandict = sqlite_db[data['transcript']]
            tdf = alltrandict["unambig"] 
            # print("alltrandict", alltrandict)
            if ("ambig" in alltrandict):
                # print("hellow ambig")
                tdf[0:0]= alltrandict["ambig"]
            # TODO: Change merge_dicts to take a list of dicts instead of two
            # NOTE: pcr not found in databases, explore for other data
            if "pcr" in data:  # TODO: Convert this value as ambig and unambig
                if "unambig_pcr" in alltrandict:
                    tdf[0:0] = alltrandict["unambig_pcr"]

                if ("ambig" in data) and "ambig_pcr" in alltrandict:
                    tdf[0:0] = alltrandict["ambig_pcr"]

            # print(
            #     "tdfbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
            #     fl, tdf)

            master_file_dict[fl] = pl.DataFrame(tdf, schema=['readlen', 'pos', 'count']).group_by('readlen', 'pos').sum()
        except Exception as e:
            print(e)
            print("hellow")
    range_set = set(range(data["minread"], data["maxread"] + 1))
    if "subcodon" not in data:
        print(master_file_dict, 'aaaaaaaaaaa')
        master_file_dict_values = pl.concat(master_file_dict.values()).filter(
(pl.col("readlen") >= data["minread"]) & (pl.col("readlen") <= data["maxread"])).select(
                'pos', 'count').group_by('pos').sum()
        print(master_file_dict_values, "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK")
        # data["coverage"] = True  # For testing purpose only
        if "coverage" in data:

            # TODO: Simply the coverage value
            master_file_dict_values = master_file_dict_values.with_columns(
                pos=pl.col("pos") + 15)
        else:
            master_file_dict_values = master_file_dict_values.with_columns(
                pos=pl.col("pos") + 16)

    master_dict_sub = []
    if ("subcodon" in data) and ("coverage" not in data):
        for filename in set(offset_dict) & set(master_file_dict):
            for readlen in set(master_file_dict[filename]) & range_set & set(
                    offset_dict[filename]):
                offset = offset_dict[filename][readlen] + 1
                for pos in master_file_dict[filename][readlen]:
                    count = master_file_dict[filename][readlen][pos]
                    if data["primetype"] == "threeprime":
                        pos += readlen
                    if "coverage" in data:  # WARN:this shouldn't be here??
                        for i in range(0, readlen, 3):
                            new_offset_pos = (i + pos) + (offset % 3)
                            master_dict_sub.append([new_offset_pos, count])

                    else:
                        offset_pos = pos + offset
                        master_dict_sub.append([offset_pos, count])

    # master_dict = pl.DataFrame(master_dict, schema=["pos", "count"])
    # master_dict_sub = pl.DataFrame(master_dict_sub,
    #                                schema=["pos", "count"]).filter(
    #                                    pl.col('pos').is_in(master_dict['pos']))
    # master_dict = pl.concat([master_dict,
    #                          master_dict_sub]).group_by('pos').agg(
    #                              pl.col('count').sum()).sort("count")
    # del master_dict_sub
    #
    # mismatch_dict = pl.DataFrame(mismatch_dict, schema=["pos", "nuc", "count"])
    # if 'mismatch' in data:
    #
    #     mismatch_dict = mismatch_dict.filter(
    #         pl.sum_horizontal('A', 'T', 'G', 'C') > 0)
    return master_file_dict_values, mismatch_dict
