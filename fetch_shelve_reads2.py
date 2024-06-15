from typing import Dict, Tuple, Union
from time import sleep
import collections
from fixed_values import merge_dicts
from sqlitedict import SqliteDict
import polars as pl
import pandas as pd


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
def get_reads(data, ) -> Tuple[pl.DataFrame, pl.DataFrame]:
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
    print(data)
    for fl in data['file_paths_dict']['path']:
        try:
            sqlite_db = SqliteDict(fl, autocommit=False)
        except FileNotFoundError:
            return pd.DataFrame(), pd.DataFrame()
        try:
            all_offsets = sqlite_db["offsets"][data["primetype"]]["offsets"]
            all_offsets = pl.DataFrame({
                'read_len': all_offsets.keys(),
                'offset': all_offsets.values()
            })
            scores = sqlite_db["offsets"][data["primetype"]]["read_scores"]
            scores = pl.DataFrame({
                'read_len': scores.keys(),
                'read_score': scores.values()
            })
            all_offsets_n_scores = all_offsets.join(
                scores, on='read_len', how='outer',
                coalesce=True).with_columns([
                    pl.col('offset').fill_null(15),
                    pl.col('read_score').fill_null(1)
                ])

        except KeyError:
            read_length = list(range(data["minread"], data["maxread"] + 1))
            range_len = data["maxread"] - data["minread"] + 1
            all_offsets_n_scores = pl.DataFrame({
                'read_len': read_length,
                'offset': [15] * range_len,
                'read_score': [1] * range_len
            })
        offset_dict[fl] = all_offsets_n_scores.filter(
            pl.col('read_score') >= data["readscore"])

        if "mismatch" in data:
            try:
                sqlite_db_seqvar = sqlite_db[data['transcript']]["seq"]
                print(sqlite_db_seqvar)
                sleep(10)

                for pos in sqlite_db_seqvar:
                    for nuc, count in sqlite_db_seqvar[pos]:
                        # convert to one based
                        mismatch_dict.append([pos + 1, nuc, count])

            except Exception:
                pass

        try:
            alltrandict = sqlite_db[data['transcript']]
            unambig_tran_dict = alltrandict["unambig"]
            ambig_tran_dict = alltrandict["ambig"] if ("ambig"
                                                       in alltrandict) else {}
            # TODO: Change merge_dicts to take a list of dicts instead of two
            trandict = merge_dicts(unambig_tran_dict, ambig_tran_dict)
            if "pcr" in data:  # TODO: Convert this value as ambig and unambig
                if "unambig_pcr" in alltrandict:
                    trandict = merge_dicts(trandict,
                                           alltrandict["unambig_pcr"])
                if ("ambig" in data) and "ambig_pcr" in alltrandict:
                    trandict = merge_dicts(trandict, alltrandict["ambig_pcr"])

            master_file_dict[fl] = trandict
        except Exception:
            pass
    range_set = set(range(data["minread"], data["maxread"] + 1))
    if "subcodon" not in data:
        for filename in master_file_dict:
            for readlen in set(master_file_dict[filename]) & range_set:
                for pos in master_file_dict[filename][readlen]:
                    count = master_file_dict[filename][readlen][pos]
                    if "coverage" in data:
                        for i in range(pos, pos + (readlen + 1)):
                            master_dict.append([i, count])
                        # use this so line graph does not have 'ramps'
                    else:
                        offset = pos + 15
                        master_dict.append([offset + 1, count])

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

    master_dict = pl.DataFrame(master_dict, schema=["pos", "count"])
    master_dict_sub = pl.DataFrame(master_dict_sub,
                                   schema=["pos", "count"]).filter(
                                       pl.col('pos').is_in(master_dict['pos']))
    master_dict = pl.concat([master_dict,
                             master_dict_sub]).group_by('pos').agg(
                                 pl.col('count').sum()).sort("count")
    del master_dict_sub

    mismatch_dict = pl.DataFrame(mismatch_dict, schema=["pos", "nuc", "count"])
    if 'mismatch' in data:

        mismatch_dict = mismatch_dict.filter(
            pl.sum_horizontal('A', 'T', 'G', 'C') > 0)
    return master_dict, mismatch_dict
