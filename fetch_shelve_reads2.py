from typing import Dict, Tuple, Union
import os
import collections
from bokeh.palettes import all_palettes
from fixed_values import merge_dicts
from sqlitedict import SqliteDict
import pandas as pd
from pandas.core.frame import DataFrame


# Create dictionary of read counts at each position in a transcript
def get_reads(
    read_type: str,
    min_read: int,
    max_read: int,
    tran: str,
    user_files: Dict[str, Dict[str, str]],
    tranlen: int,
    coverage,
    # organism,
    subcodon: bool,
    # noisered,
    primetype: str,
    filetype: str,
    readscore: int,
    # secondary_readscore=1,
    pcr: bool = False,
    get_mismatches: bool = False,
    # self_obj=None
) -> Union[None, Tuple[Dict[int, int], Dict[int, int]], Tuple[str, Union[
        str, Dict[str, Dict[int, int]]]], Tuple[str, DataFrame]]:
    mismatch_dict = pd.DataFrame(0,
                                 index=range(0, tranlen + 1),
                                 columns=["A", "T", "G", "C"])

    master_dict = collections.OrderedDict()
    for i in range(0, tranlen + 1):
        master_dict[i] = 0
    master_file_dict = {}

    # first make a master dict consisting of all the read dicts from each filename
    offset_dict = {}
    if read_type == "unambig":
        if filetype in user_files:
            for file_id in user_files[filetype]:
                filename = user_files[filetype][file_id]
                try:
                    sqlite_db = SqliteDict(filename, autocommit=False)
                except FileNotFoundError:
                    return "Could not open file:", filename.split("/")[-1]
                accepted_offsets = {}
                try:
                    all_offsets = sqlite_db["offsets"][primetype]["offsets"]
                    scores = sqlite_db["offsets"][primetype]["read_scores"]
                except KeyError:
                    all_offsets = {}
                    scores = {}
                    for i in range(min_read, max_read):
                        scores[i] = 1
                        all_offsets[i] = 15
                for rl in scores:
                    if scores[rl] >= readscore:
                        if rl in all_offsets:
                            accepted_offsets[rl] = all_offsets[rl]
                        else:
                            accepted_offsets[rl] = 15

                if accepted_offsets:
                    offset_dict[filename] = {}
                    for length in accepted_offsets:
                        offset_dict[filename][length] = accepted_offsets[
                            length]

                if get_mismatches:
                    try:
                        sqlite_db_seqvar = sqlite_db[tran]["seq"]

                        for pos in sqlite_db_seqvar:
                            # convert to one based
                            fixed_pos = pos + 1
                            for char in sqlite_db_seqvar[pos]:
                                mismatch_dict.loc[
                                    fixed_pos,
                                    char] += sqlite_db_seqvar[pos][char]

                    except Exception:
                        pass

                try:
                    alltrandict = sqlite_db[tran]
                    unambig_tran_dict = alltrandict["unambig"]
                    if pcr:
                        if "unambig_pcr" in alltrandict:
                            unambig_tran_dict = merge_dicts(
                                unambig_tran_dict, alltrandict["unambig_pcr"])
                    sqlite_db.close()
                    master_file_dict[filename] = unambig_tran_dict
                except Exception:
                    pass
    else:
        if filetype in user_files:
            for file_id in user_files[filetype]:
                filename = user_files[filetype][file_id]
                if os.path.isfile(filename):
                    sqlite_db = SqliteDict(filename, autocommit=False)
                else:
                    return "File not found", mismatch_dict
                accepted_offsets = {}
                if "offsets" in sqlite_db:
                    all_offsets = sqlite_db["offsets"][primetype]["offsets"]
                    scores = sqlite_db["offsets"][primetype]["read_scores"]
                else:
                    all_offsets = {}
                    scores = {}
                # oyster has no readscores so subcodon profiles not displaying, so i put this in to give every readlength a score of 1
                if scores == {}:
                    for i in range(min_read, max_read):
                        scores[i] = 1
                        all_offsets[i] = 15

                for rl in scores:
                    if float(scores[rl]) >= float(readscore):
                        if rl in all_offsets:
                            accepted_offsets[rl] = all_offsets[rl]
                        else:
                            accepted_offsets[rl] = 15

                if accepted_offsets != {}:
                    offset_dict[filename] = {}
                    for length in accepted_offsets:
                        offset_dict[filename][length] = accepted_offsets[
                            length]
                try:
                    sqlite_db_seqvar = sqlite_db[tran]["seq"]
                    for pos in sqlite_db_seqvar:
                        # convert to one based
                        fixed_pos = pos + 1
                        for char in sqlite_db_seqvar[pos]:
                            mismatch_dict.loc[
                                fixed_pos, char] += sqlite_db_seqvar[pos][char]
                except Exception:
                    pass

                try:
                    alltrandict = sqlite_db[tran]
                    sqlite_db.close()
                    unambig_tran_dict = alltrandict["unambig"]
                    if "ambig" in alltrandict:
                        ambig_tran_dict = alltrandict["ambig"]
                    else:
                        ambig_tran_dict = {}
                    # TODO: Change merge_dicts to take a list of dicts instead of two
                    trandict = merge_dicts(unambig_tran_dict, ambig_tran_dict)
                    if pcr:  # TODO: Convert this value as ambig and unambig
                        if "unambig_pcr" in alltrandict:
                            trandict = merge_dicts(trandict,
                                                   alltrandict["unambig_pcr"])
                        if "ambig_pcr" in alltrandict:
                            trandict = merge_dicts(trandict,
                                                   alltrandict["ambig_pcr"])
                    master_file_dict[filename] = trandict
                except Exception:
                    pass
    # Next check coverage, if that's true then calculate coverage for each rl and return dict
        if not subcodon:
            for filename in master_file_dict:
                for readlen in master_file_dict[filename]:
                    if readlen >= min_read and readlen <= max_read:
                        for pos in master_file_dict[filename][readlen]:
                            count = master_file_dict[filename][readlen][pos]
                            if coverage:
                                if pos != 0 and pos - 1 not in master_dict:
                                    master_dict[pos - 1] = 0
                                i = 0
                                for i in range(pos, pos + (readlen + 1)):
                                    if i in master_dict:
                                        master_dict[i] += count
                                    else:
                                        master_dict[i] = count
                                # use this so line graph does not have 'ramps'
                                if i + 1 not in master_dict:
                                    master_dict[i + 1] = 0
                            else:
                                offset_pos = pos + 15
                                if offset_pos + 1 not in master_dict:
                                    master_dict[offset_pos + 1] = 0
                                master_dict[offset_pos] += count

        # Fetching subcodon reads
        if subcodon:
            if not coverage:
                for filename in master_file_dict:
                    if filename not in offset_dict:
                        continue
                    for readlen in master_file_dict[filename]:
                        if readlen >= min_read and readlen <= max_read:
                            if readlen in offset_dict[filename]:
                                offset = offset_dict[filename][readlen] + 1
                                for pos in master_file_dict[filename][readlen]:
                                    count = master_file_dict[filename][
                                        readlen][pos]
                                    if primetype == "threeprime":
                                        pos += readlen

                                    if coverage:
                                        for i in range(0, readlen, 3):
                                            new_offset_pos = (i + pos) + (
                                                offset % 3)
                                            try:
                                                master_dict[
                                                    new_offset_pos] += count
                                            except KeyError:
                                                pass
                                    else:
                                        offset_pos = pos + offset
                                        try:
                                            master_dict[offset_pos] += count
                                        except KeyError:
                                            print("Error tried adding to "
                                                  f"position {e} but tranlen "
                                                  f"is only {tranlen}")
        if not get_mismatches:
            mismatch_dict = mismatch_dict[mismatch_dict.sum(axis=1) > 0]
        return master_dict, mismatch_dict


# Create dictionary of counts at each position, averged by readlength
def get_readlength_breakdown(
    read_type: str,
    min_read: int,
    max_read: int,
    tran: str,
    user_files: Dict[str, Dict[str, str]],
    # offset_dict,
    tranlen: int,
    # subcodon, noisered, primetype, preprocess,
    coverage: int,  # organism,
    filetype: str,
    colorbar_minread: int,
    colorbar_maxread: int,
) -> Tuple[Dict[int, str], Dict[str, Dict[int, int]]]:
    master_dict = {}
    color_range = float(colorbar_maxread - colorbar_minread)
    color_list = all_palettes["RdYlGn"][10]

    for i in range(0, tranlen + max_read):
        master_dict[i] = {}
        for x in range(min_read, max_read + 1):
            master_dict[i][x] = 0
    # the keys of master readlen dict are readlengths the value is a dictionary
    # of position:count, there is also a colour key
    colored_master_dict = {}
    master_file_dict = {}
    # first make a master dict consisting of all the read dicts from each filename
    if filetype in user_files:
        for file_id in user_files[filetype]:
            filename = user_files[filetype][file_id]
            try:
                openshelf = SqliteDict(filename)
                alltrandict = dict(openshelf[tran])
                trandict = alltrandict["unambig"]
                if read_type == "ambig":
                    trandict = merge_dicts(trandict, alltrandict["ambig"])
                master_file_dict[filename] = trandict
                openshelf.close()
            except KeyError:
                pass

    for filename in master_file_dict:
        for readlen in master_file_dict[filename]:
            if readlen >= min_read and readlen <= max_read:
                if coverage:
                    for pos in master_file_dict[filename][readlen]:
                        count = master_file_dict[filename][readlen][pos]
                        for i in range(pos, pos + (readlen + 1)):
                            master_dict[i][readlen] += count
                else:
                    if os.path.isfile(filename):
                        openshelf = SqliteDict(filename)
                    else:
                        continue
                    offsets = openshelf["offsets"]["fiveprime"]["offsets"]
                    if readlen in offsets:
                        offset = offsets[readlen]
                    else:
                        offset = 15
                    for pos in master_file_dict[filename][readlen]:
                        count = master_file_dict[filename][readlen][pos]
                        if (pos + offset) in master_dict:
                            master_dict[pos + offset][readlen] += count
    sorted_master_dict = collections.OrderedDict()
    for key in sorted(master_dict.keys()):
        sorted_master_dict[key] = master_dict[key]

    for pos in sorted_master_dict:
        count = sum(sorted_master_dict[pos].values())
        tot_readlen = 0.0
        tot_count = 0.0001
        for readlen in sorted_master_dict[pos]:
            tot_count += sorted_master_dict[pos][readlen]
            tot_readlen += (sorted_master_dict[pos][readlen] * readlen)
        avg_readlen = int(tot_readlen / tot_count)
        if avg_readlen < colorbar_minread:
            avg_readlen = colorbar_minread
        if avg_readlen > colorbar_maxread:
            avg_readlen = colorbar_maxread
        # find where this avg readlen lies in the range of min readlen
        # to max readlen and use that to assign a color

        y = avg_readlen - colorbar_minread
        per = y / color_range
        final_per = int(per * 10)
        if final_per > 9:
            final_per = 9
        color = color_list[final_per]
        if color not in colored_master_dict:
            colored_master_dict[color] = collections.OrderedDict()
            for i in range(0, tranlen + max_read + 1):
                colored_master_dict[color][i] = 0
        if pos not in colored_master_dict[color]:
            colored_master_dict[color][pos] = 0
        colored_master_dict[color][pos] += count
    return color_list, colored_master_dict


# Create a dictionary of mismatches at each position
def get_seq_var(
    user_files: Dict[str, Dict[str, str]],
    tranlen: int,
    #organism,
    tran: str,
) -> Union[str, Dict[str, Dict[str, int]]]:
    mismatch_dict = pd.DataFrame(0,
                                 index=range(1, tranlen + 1),
                                 columns=["A", "T", "G", "C"])
    for filetype in ["riboseq", "rnaseq"]:
        if filetype in user_files:
            for file_id in user_files[filetype]:
                filename = user_files[filetype][file_id]
                if os.path.isfile(filename):
                    sqlite_db = SqliteDict(filename, autocommit=False)
                else:
                    return "File not found"
                if tran in sqlite_db:
                    if "seq" in sqlite_db[tran]:
                        sqlite_db_seqvar = dict(sqlite_db[tran]["seq"])

                        for pos in sqlite_db_seqvar:
                            #convert to one based
                            fixed_pos = pos + 1
                            for char in sqlite_db_seqvar[pos]:
                                if char != "N":
                                    mismatch_dict.loc[
                                        fixed_pos,
                                        char] += sqlite_db_seqvar[pos][char]

                    sqlite_db.close()

    return mismatch_dict
