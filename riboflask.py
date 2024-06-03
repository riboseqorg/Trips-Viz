import config
import os
import numpy as np
import polars as pl
from plots import VegaPlot
from fetch_shelve_reads2 import get_reads
from sqlitedict import SqliteDict
import collections
import matplotlib.pyplot as plt
import matplotlib
import fixed_values
from fixed_values import get_user_defined_seqs
from sqlqueries_2 import get_table, sqlquery

matplotlib.use("agg")


def generate_plot(data, settings) -> str:
    seq_rules = {
        "proteomics": {
            "frame_breakdown": 1
        },
        "conservation": {
            "frame_breakdown": 1
        },
        "tcpseq": {
            "frame_breakdown": 0
        }
    }

    if ("line" not in data) and ("ribocoverage" in data):
        # TODO: Convert this into notification
        return ("Error: Cannot display Ribo-Seq Coverage when 'Line Graph'" +
                " is turned off")

    # Label and visibility for the interactive legends
    labels_visibility = {
        "Frame 1 profiles": True,
        "Frame 2 profiles": True,
        "Frame 3 profiles": True,
        "RNA": True,
        "Exon Junctions": True,
        "CDS markers": True,
    }
    if "mismatches" in data:
        for nuc in "ATGC":
            labels_visibility[f"Mismatches {nuc}"] = False
    # This is a list of booleans that decide if the interactive legends boxes are filled in or not.Needs to be same length as labels
    frame_orfs = {1: [], 2: [], 3: []}
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
    print(traninfo, 'ha ha ha')
    data["tranlen"] = traninfo['length']
    for ss in ['start_list', 'stop_list', 'exon_junctions']:
        try:
            traninfo[ss] = np.array([int(x) for x in traninfo[ss].split(",")])
        except Exception:
            traninfo[ss] = []
    # TODO: Replace next with fill na in dataframe
    if not traninfo["cds_start"]:
        traninfo["cds_start"] = 0
    if not traninfo["cds_stop"]:
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

    seq = traninfo["sequence"].upper(
    )  # NOTE: I guess it is already upper case

    start_dataframe = pl.DataFrame({
        'frame': ((traninfo['start_list'] - 1) % 3) + 1,
        'type': ['start'] * len(traninfo['start_list']),
        'pos':
        traninfo['start_list']
    })
    stop_dataframe = []
    orf_dataframe = []
    for i in range(len(seq)):
        if seq[i:i + 3] in ["TAG", "TAA", "TGA"]:
            frame = (i % 3) + 1
            stop_dataframe.append([frame, seq[i:i + 3], i + 1])
            orf_dataframe.append([
                frame,
                start_dataframe.filter((pl.col('frame') == frame)
                                       & (pl.col('pos') < i + 1)).select(
                                           pl.last('pos'))[0], i + 1
            ])
    stop_dataframe = pl.DataFrame(stop_dataframe,
                                  schema=["frame", "type", "pos"])
    orf_dataframe = pl.DataFrame(orf_dataframe,
                                 schema=["frame", "start", "pos"])
    # Error occurs if one of the frames is empty for any given start/stop, so we initialise with -5 as this won't be seen by user and will prevent the error
    start_stop_dataframe = pl.DataFrame({
        'frame': [1, 2, 3],
        'type': ['start'] * 3,
        'pos': [-5] * 3
    })
    start_stop_dataframe = pl.concat(
        [start_stop_dataframe, start_dataframe, stop_dataframe])
    start_stop = start_stop_dataframe.with_columns(
        pl.col("type").apply(lambda x: x if x == 'start' else 'stop'))

    # TODO: add avarible of rnaseq

    all_rna_reads, rna_seqvar_dict = get_reads(
        data
    )  # TODO: keep rna_Seqvar_dict and ribo_seq_var_dict meltated to compine them
    all_rna_reads = all_rna_reads.with_columns(frame=pl.col("pos") % 3 + 1)
    plt = VegaPlot(all_rna_reads)
    lineplot = plt.lineplot("pos", "count")
    start_stop_plot = []
    for frame in [1, 2, 3]:
        start_stop_plot.append(
            VegaPlot(start_stop.filter(pl.col("frame") == frame)).vline(
                "pos", last=True if frame == 4 else False))
    # plt2 = VegaPlot(all_rna_reads).scatter("pos", "count")
    seq_frames = pl.DataFrame({
        'sequence': list(seq),
        'pos': range(len(seq)),
        'y': [0] * len(seq)
    }).with_columns(frame=pl.col('pos') % 3 + 1)
    print(seq_frames)
    seq_plot = VegaPlot(seq_frames.to_pandas()).seq_plot()
    print(seq_plot.to_json())
    plt = plt.vact_plot_json([lineplot] + start_stop_plot + [seq_plot])

    print(all_rna_reads, rna_seqvar_dict, "Anmol Kiran You are here")
    return plt.to_json()
    # self.update_state(state='PROGRESS',meta={'current': 100, 'total': 100,'status': "Fetching Ribo-Seq Reads"})
    # TODO: Add a variable of RiboSeq
    # all_subcodon_reads, ribo_seqvar_dict = get_reads(data)
    # print(all_subcodon_reads, ribo_seqvar_dict)

    # seq_var_dict = fixed_values.merge_dicts(ribo_seqvar_dict, rna_seqvar_dict)
    # rnamax = all_rna_reads['count'].max()
    # subcodonmax = all_subcodon_reads['count'].max()
    # y_max = max(1, rnamax, subcodonmax) * 1.1
