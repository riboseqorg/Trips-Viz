import collections
import os

import numpy as np
import polars as pl
from sqlitedict import SqliteDict

import config
import fixed_values
from core_functions import sequence2rdg
from fetch_shelve_reads2 import get_reads
from fixed_values import get_user_defined_seqs
from sqlqueries_2 import get_table, sqlquery


def generate_plot(data, settings) -> str:

    if ("line" not in data) and ("ribocoverage" in data):
        # TODO: Convert this into notification
        return ("Error: Cannot display Ribo-Seq Coverage when 'Line Graph'" +
                " is turned off")

    # This is a list of booleans that decide if the interactive legends boxes are filled in or not.Needs to be same length as labels
    if data['owner'] == 1:
        sqlpath = "{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC,
                                                      config.ANNOTATION_DIR,
                                                      data["organism"],
                                                      data["transcriptome"])
        if not os.path.isfile(sqlpath):
            return "Cannot find annotation file {}.{}.sqlite".format(
                data["organism"], data["transcriptome"])
    else:
        sqlpath = "{0}/transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
            config.UPLOADS_DIR, data['owner'], data["organism"],
            data["transcriptome"])
    traninfo = sqlquery(sqlpath, "transcripts").filter(
        pl.col("transcript") == data["transcript"])[0].to_dicts()[0]
    data["tranlen"] = traninfo['length']
    for ss in ['start_list', 'stop_list', 'exon_junctions']:
        try:
            traninfo[ss] = np.array([int(x) for x in traninfo[ss].split(",")])
        except Exception:
            traninfo[ss] = []
    # TODO: Replace next with fill na in dataframe
    if not traninfo["cds_start"]:
        traninfo["cds_start"] = 0
    # if not traninfo["cds_stop"]:
        traninfo["cds_stop"] = 0

    try:
        coding_regions = sqlquery(sqlpath, "coding_regions").filter(
            pl.col("transcript") == data["transcript"]).select(
                "coding_start", "coding_stop")
    except Exception:  # pragma: no cover
        coding_regions = pl.DataFrame(schema={
            "coding_start": int,
            "coding_stop": int
        })

    data["coding_regions"] = coding_regions
    data["traninfo"] = traninfo





    # TODO: I need to orf in sqlite file
    orfs = SqliteDict(sqlpath.replace(".sqlite","_helper.sqlite")) # TODO:: Check if it will work 
    start_stop_dataframe = pl.DataFrame(orfs["start_stop_list"], schema=["frame", "pos", "type"])
    orf_dataframe = pl.DataFrame(orfs["orfs"], schema=["frame", "start", "stop"])
    # TODO: add avarible of rnaseq

    all_rna_reads, rna_seqvar_dict = get_reads(
        data
    )  # TODO: keep rna_Seqvar_dict and ribo_seq_var_dict meltated to compine them
    all_rna_reads = all_rna_reads.with_columns(frame=pl.col("pos") % 3 + 1)
    
    rdg = str(sequence2rdg(seq))
    # print(rdg)

    return {"plot": all_rna_reads.sort("pos").write_csv(), "rdg": rdg, 'start_stop': start_stop.to_pandas().to_json(orient="records"),'seq':seq} # all_rna_reads.sort("pos").write_csv()
    # return plt.to_json()
