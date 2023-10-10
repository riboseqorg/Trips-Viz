from typing import List, Dict
import config
import os
import sqlite3
from fetch_shelve_reads2 import get_reads
import pandas as pd
from sqlitedict import SqliteDict
import collections
import mpld3
import matplotlib.pyplot as plt
import matplotlib
import fixed_values
from fixed_values import get_user_defined_seqs

matplotlib.use('agg')


def generate_plot(
        tran: str,
        ambig: str,
        min_read: int,
        max_read: int,
        lite: str,
        ribocoverage: bool,
        organism: str,
        readscore: str,
        noisered: str,
        primetype: str,
        # minfiles,
        nucseq: str,
        user_hili_starts: str,
        user_hili_stops: str,
        # uga_diff,
        file_paths_dict: Dict[str, str],
        short_code: str,
        # color_readlen_dist,
        background_col: str,
        uga_col: str,
        uag_col: str,
        uaa_col: str,
        # advanced, seqhili, seq_rules,
        # title_size,
        subheading_size: int,
        axis_label_size: int,
        marker_size: int,
        transcriptome: str,
        trips_uploads_location: str,
        cds_marker_size: int,
        cds_marker_colour: int,
        legend_size: int,
        ribo_linewidth: int,
        secondary_readscore: str,
        pcr: bool,
        mismatches: bool,
        hili_start: bool,
        hili_stop: bool) -> str:

    if lite == "n" and ribocoverage:
        return ("Error: Cannot display Ribo-Seq Coverage when 'Line Graph'" +
                " is turned off")
    labels = [
        "Frame 1 profiles", "Frame 2 profiles", "Frame 3 profiles", "RNA",
        "Exon Junctions"
    ]
    start_visible = [True] * 5
    if mismatches:
        labels += [f"Mismatches {nuc}" for nuc in "ATGC"]
        start_visible += [False] * 4
    start_visible.append(True)
    labels.append("CDS markers")
    # This is a list of booleans that decide if the interactive legends boxes are filled in or not.Needs to be same length as labels
    stop_codons = ["TAG", "TAA", "TGA"]
    frame_orfs = {1: [], 2: [], 3: []}
    owner = get_table("organisms")
    owner = owner.loc[(owner.organism_name == organism) &
                      (owner.transcriptome_list == transcriptome),
                      "owner"].values[0]
    if owner == 1:
        sqlpath = "{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC,
                                                      config.ANNOTATION_DIR,
                                                      organism, transcriptome)
        if not os.path.isfile(sqlpath):
            return_str = "Cannot find annotation file {}.{}.sqlite".format(
                organism, transcriptome)
            return return_str
    else:
        sqlpath = "{0}/transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
            trips_uploads_location, owner, organism, transcriptome)
    transcripts = get_table("transcripts")
    transcripts = transcripts[transcripts.transcript == tran].iloc[0]
    traninfo = {
        "transcript": result[0],
        "gene": result[1],
        "length": result[2],
        "cds_start": result[3],
        "cds_stop": result[4],
        "seq": result[5],
        "strand": result[6],
        "stop_list": result[7].split(","),
        "start_list": result[8].split(","),
        "exon_junctions": result[9].split(","),
        "tran_type": result[10],
        "principal": result[11]
    }
    try:
        traninfo["stop_list"] = [int(x) for x in traninfo["stop_list"]
                                 ]  #TODO: make it sure that it is list of int
    except Exception:
        traninfo["stop_list"] = []

    try:
        traninfo["start_list"] = [int(x) for x in traninfo["start_list"]]
    except Exception:
        traninfo["start_list"] = []

    if str(traninfo["exon_junctions"][0]) != "":
        traninfo["exon_junctions"] = [
            int(x) for x in traninfo["exon_junctions"]
        ]
    else:
        traninfo["exon_junctions"] = []
    all_cds_regions = []
    # Check if the 'coding_regions' table exists
    cursor.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='coding_regions';"
    )
    result = cursor.fetchone()
    if result:
        cursor.execute(
            "SELECT * from coding_regions WHERE transcript = '{}'".format(
                tran))
        result = cursor.fetchall()
        for row in result:
            all_cds_regions.append((row[1], row[2]))
    transhelve.close()
    gene = traninfo["gene"]
    tranlen = traninfo["length"]
    cds_start = traninfo["cds_start"]
    cds_stop = traninfo["cds_stop"]
    if not cds_start:
        cds_start = 0
    if not cds_stop:
        cds_stop = 0
    all_starts = traninfo["start_list"]
    all_stops = {"TAG": [], "TAA": [], "TGA": []}
    exon_junctions = traninfo["exon_junctions"]
    seq = traninfo["seq"].upper()
    for i in range(0, len(seq)):
        if seq[i:i + 3] in stop_codons:
            all_stops[seq[i:i + 3]].append(i + 1)
    # Error occurs if one of the frames is empty for any given start/stop, so we initialise with -5 as this won't be seen by user and will prevent the error
    start_stop_dict = {
        1: {
            "starts": [-5],
            "stops": {
                "TGA": [-5],
                "TAG": [-5],
                "TAA": [-5]
            }
        },
        2: {
            "starts": [-5],
            "stops": {
                "TGA": [-5],
                "TAG": [-5],
                "TAA": [-5]
            }
        },
        3: {
            "starts": [-5],
            "stops": {
                "TGA": [-5],
                "TAG": [-5],
                "TAA": [-5]
            }
        }
    }
    for start in all_starts:
        rem = ((start - 1) % 3) + 1
        start_stop_dict[rem]["starts"].append(start)
    for stop in all_stops:
        for stop_pos in all_stops[stop]:
            rem = ((stop_pos - 1) % 3) + 1
            start_stop_dict[rem]["stops"][stop].append(stop_pos)
    # find all open reading frames
    for frame in [1, 2, 3]:
        for start in start_stop_dict[frame]["starts"]:
            best_stop_pos = 10000000
            for stop in start_stop_dict[frame]["stops"]:
                for stop_pos in start_stop_dict[frame]["stops"][stop]:
                    if stop_pos > start and stop_pos < best_stop_pos:
                        best_stop_pos = stop_pos
            if best_stop_pos != 10000000:
                frame_orfs[frame].append((start, best_stop_pos))
    # self.update_state(state='PROGRESS',meta={'current': 100, 'total': 100,'status': "Fetching RNA-Seq Reads"})
    all_rna_reads, rna_seqvar_dict = get_reads(ambig,
                                               min_read,
                                               max_read,
                                               tran,
                                               file_paths_dict,
                                               tranlen,
                                               True,
                                               organism,
                                               False,
                                               noisered,
                                               primetype,
                                               "rnaseq",
                                               readscore,
                                               pcr,
                                               get_mismatches=mismatches)
    # self.update_state(state='PROGRESS',meta={'current': 100, 'total': 100,'status': "Fetching Ribo-Seq Reads"})
    all_subcodon_reads, ribo_seqvar_dict = get_reads(ambig,
                                                     min_read,
                                                     max_read,
                                                     tran,
                                                     file_paths_dict,
                                                     tranlen,
                                                     ribocoverage,
                                                     organism,
                                                     True,
                                                     noisered,
                                                     primetype,
                                                     "riboseq",
                                                     readscore,
                                                     secondary_readscore,
                                                     pcr,
                                                     get_mismatches=mismatches)
    seq_var_dict = fixed_values.merge_dicts(ribo_seqvar_dict, rna_seqvar_dict)
    try:
        rnamax = max(all_rna_reads.values())
    except Exception:
        rnamax = 0
    try:
        subcodonmax = max(all_subcodon_reads.values())
    except Exception:
        subcodonmax = 0
    y_max = max(1, rnamax, subcodonmax) * 1.1
    fig = plt.figure(figsize=(13, 8))
    ax_main = plt.subplot2grid((31, 1), (0, 0), rowspan=22)

    ax_main.spines['bottom'].set_visible(False)
    for s in ['bottom', 'left', 'top', 'right']:
        ax_main.spines[s].set_linewidth(15)
        ax_main.spines[s].set_color("red")
    alt_seq_type_vars = []
    # Store the  counts of any alt seq types so they can be written to csv file later
    alt_seq_dict = {}
    # Plot any alternative sequence types if there are any
    for seq_type in file_paths_dict:

        if seq_type != "riboseq" and seq_type != "rnaseq":
            if file_paths_dict[seq_type] == {}:
                continue
            if seq_type in seq_rules:
                if seq_rules[seq_type]["frame_breakdown"] == 1:
                    frame_breakdown = True
                else:
                    frame_breakdown = False
            else:
                frame_breakdown = False
            alt_sequence_reads, _ = get_reads(ambig, min_read, max_read, tran,
                                              file_paths_dict, tranlen, True,
                                              organism, frame_breakdown,
                                              noisered, primetype, seq_type,
                                              readscore)
            alt_seq_dict[seq_type] = alt_sequence_reads
            if not frame_breakdown:
                alt_seq_plot = ax_main.plot(alt_sequence_reads.keys(),
                                            alt_sequence_reads.values(),
                                            alpha=1,
                                            label=seq_type,
                                            zorder=2,
                                            color='#5c5c5c',
                                            linewidth=2)
                labels.append(seq_type)
                start_visible.append(True)
                alt_seq_type_vars.append(alt_seq_plot)
            else:
                alt_frame_counts = {
                    0: collections.OrderedDict(),
                    1: collections.OrderedDict(),
                    2: collections.OrderedDict()
                }
                for key in alt_sequence_reads:
                    frame = (key + 1) % 3  # key is start pos
                    alt_frame_counts[frame][key] = alt_sequence_reads[key]
                frame0_altseqplot = ax_main.plot(alt_frame_counts[0].keys(),
                                                 alt_frame_counts[0].values(),
                                                 alpha=0.75,
                                                 label=seq_type + "frame0",
                                                 zorder=2,
                                                 color="#FF4A45",
                                                 linewidth=2)
                frame1_altseqplot = ax_main.plot(alt_frame_counts[1].keys(),
                                                 alt_frame_counts[1].values(),
                                                 alpha=0.75,
                                                 label=seq_type + "frame1",
                                                 zorder=2,
                                                 color="#64FC44",
                                                 linewidth=2)
                frame2_altseqplot = ax_main.plot(alt_frame_counts[2].keys(),
                                                 alt_frame_counts[2].values(),
                                                 alpha=0.75,
                                                 label=seq_type + "frame2*",
                                                 zorder=2,
                                                 color="#5687F9",
                                                 linewidth=2)
                labels.append(seq_type + "frame 1")
                labels.append(seq_type + "frame 2")
                labels.append(seq_type + "frame 3")
                start_visible.append(True)
                start_visible.append(True)
                start_visible.append(True)
                alt_seq_type_vars.append(frame0_altseqplot)
                alt_seq_type_vars.append(frame1_altseqplot)
                alt_seq_type_vars.append(frame2_altseqplot)
            if max(alt_sequence_reads.values()) > y_max:
                y_max = max(alt_sequence_reads.values())
    label = 'Reads'
    ax_main.set_ylabel(label, fontsize=axis_label_size, labelpad=30)
    label = 'Position (nucleotides)'
    ax_main.set_xlabel(label, fontsize=axis_label_size, labelpad=-10)
    ax_main.set_ylim(0, y_max)

    if lite == "n":
        rna_bars = ax_main.bar(all_rna_reads.keys(),
                               all_rna_reads.values(),
                               alpha=1,
                               label=labels,
                               zorder=1,
                               color='lightgray',
                               linewidth=0,
                               width=1)
    else:
        rna_bars = ax_main.plot(all_rna_reads.keys(),
                                all_rna_reads.values(),
                                alpha=1,
                                label=labels,
                                zorder=1,
                                color='#a7adb7',
                                linewidth=4)

    cds_markers = ax_main.plot((cds_start, cds_start), (0, y_max * 0.97),
                               color=cds_marker_colour,
                               linestyle='solid',
                               linewidth=cds_marker_size)
    ax_main.text(cds_start,
                 y_max * 0.97,
                 "CDS start",
                 fontsize=18,
                 color="black",
                 ha="center")
    cds_markers += ax_main.plot((cds_stop + 1, cds_stop + 1),
                                (0, y_max * 0.97),
                                color=cds_marker_colour,
                                linestyle='solid',
                                linewidth=cds_marker_size)
    ax_main.text(cds_stop,
                 y_max * 0.97,
                 "CDS stop",
                 fontsize=18,
                 color="black",
                 ha="center")
    ax_cds = plt.subplot2grid((31, 1), (23, 0), rowspan=1, sharex=ax_main)
    ax_cds.set_facecolor("white")
    ax_cds.set_ylabel('Merged CDS',
                      labelpad=4,
                      verticalalignment='center',
                      horizontalalignment="right",
                      rotation="horizontal",
                      color="black",
                      fontsize=(axis_label_size / 1.5))

    ax_f1 = plt.subplot2grid((31, 1), (27, 0), rowspan=1, sharex=ax_main)
    ax_f1.set_facecolor(color_dict['frames'][0])
    ax_f2 = plt.subplot2grid((31, 1), (28, 0), rowspan=1, sharex=ax_main)
    ax_f2.set_facecolor(color_dict['frames'][1])
    ax_f3 = plt.subplot2grid((31, 1), (29, 0), rowspan=1, sharex=ax_main)
    ax_f3.set_facecolor(color_dict['frames'][2])
    ax_nucseq = plt.subplot2grid((31, 1), (30, 0), rowspan=1, sharex=ax_main)
    ax_nucseq.set_xlabel('Transcript: {} Length: {} nt'.format(tran, tranlen),
                         fontsize=subheading_size)

    for tup in all_cds_regions:
        ax_cds.fill_between([tup[0], tup[1]], [1, 1],
                            zorder=0,
                            alpha=1,
                            color="#001285")

    # plot a dummy exon junction at postion -1, needed in cases there are no exon junctions, this wont be seen
    allexons = ax_main.plot((-1, -1), (0, 1),
                            alpha=0.01,
                            color='black',
                            linestyle='-.',
                            linewidth=2)
    for exon in exon_junctions:
        allexons += ax_main.plot((exon, exon), (0, y_max),
                                 alpha=0.95,
                                 color='black',
                                 linestyle=':',
                                 linewidth=3)

    # dictionary for each frame in which the keys are the posistions and the values are the counts
    frame_counts = {
        0: collections.OrderedDict(),
        1: collections.OrderedDict(),
        2: collections.OrderedDict()
    }
    for key in all_subcodon_reads:
        rem = key % 3
        if rem == 1:  # frame 1
            frame = 0
        elif rem == 2:  # frame 2
            frame = 1
        elif rem == 0:  # frame 3
            frame = 2
        frame_counts[frame][key] = all_subcodon_reads[key]
        if lite == "n":
            frame_counts[frame][key + 1] = 0
            frame_counts[frame][key + 2] = 0

    if lite == "n":
        frame0subpro = ax_main.bar(frame_counts[0].keys(),
                                   frame_counts[0].values(),
                                   alpha=0.75,
                                   label=labels,
                                   zorder=2,
                                   color="#FF4A45",
                                   edgecolor="#FF4A45",
                                   width=1,
                                   linewidth=4)
        frame1subpro = ax_main.bar(frame_counts[1].keys(),
                                   frame_counts[1].values(),
                                   alpha=0.75,
                                   label=labels,
                                   zorder=2,
                                   color="#64FC44",
                                   edgecolor="#64FC44",
                                   width=1,
                                   linewidth=4)
        frame2subpro = ax_main.bar(frame_counts[2].keys(),
                                   frame_counts[2].values(),
                                   alpha=0.75,
                                   label=labels,
                                   zorder=2,
                                   color="#5687F9",
                                   edgecolor="#5687F9",
                                   width=1,
                                   linewidth=4)
    else:
        frame0subpro = ax_main.plot(frame_counts[0].keys(),
                                    frame_counts[0].values(),
                                    alpha=0.75,
                                    label=labels,
                                    zorder=2,
                                    color="#FF4A45",
                                    linewidth=ribo_linewidth)
        frame1subpro = ax_main.plot(frame_counts[1].keys(),
                                    frame_counts[1].values(),
                                    alpha=0.75,
                                    label=labels,
                                    zorder=2,
                                    color="#64FC44",
                                    linewidth=ribo_linewidth)
        frame2subpro = ax_main.plot(frame_counts[2].keys(),
                                    frame_counts[2].values(),
                                    alpha=0.75,
                                    label=labels,
                                    zorder=2,
                                    color="#5687F9",
                                    linewidth=ribo_linewidth)
    if mismatches:
        a_mismatches = ax_main.plot(seq_var_dict["A"].keys(),
                                    seq_var_dict["A"].values(),
                                    alpha=0.01,
                                    label=labels,
                                    zorder=2,
                                    color="purple",
                                    linewidth=2)
        t_mismatches = ax_main.plot(seq_var_dict["T"].keys(),
                                    seq_var_dict["T"].values(),
                                    alpha=0.01,
                                    label=labels,
                                    zorder=2,
                                    color="yellow",
                                    linewidth=2)
        g_mismatches = ax_main.plot(seq_var_dict["G"].keys(),
                                    seq_var_dict["G"].values(),
                                    alpha=0.01,
                                    label=labels,
                                    zorder=2,
                                    color="orange",
                                    linewidth=2)
        c_mismatches = ax_main.plot(seq_var_dict["C"].keys(),
                                    seq_var_dict["C"].values(),
                                    alpha=0.01,
                                    label=labels,
                                    zorder=2,
                                    color="pink",
                                    linewidth=2)

    xy = 0
    if nucseq:
        ax_nucseq.set_facecolor(background_col)
        mrnaseq = seq.replace("T", "U")
        color_list = ["#FF4A45", "#64FC44", "#5687F9"]
        char_frame = 0
        for _ in mrnaseq:
            ax_nucseq.text((xy + 1) - 0.1,
                           0.2,
                           mrnaseq[xy],
                           fontsize=20,
                           color=color_list[char_frame % 3])
            xy += 1
            char_frame += 1

    # If the user passed a list of sequences to highlight, find and plot them here.
    if seqhili != ['']:
        near_cog_starts, signalhtml = get_user_defined_seqs(seq, seqhili)
        for slip in near_cog_starts[0]:
            try:
                hili_sequences += ax_f1.plot((slip, slip), (0, 0.5),
                                             alpha=1,
                                             label=labels,
                                             zorder=4,
                                             color='black',
                                             linewidth=5)
            except Exception:
                hili_sequences = ax_f1.plot((slip, slip), (0, 0.5),
                                            alpha=1,
                                            label=labels,
                                            zorder=4,
                                            color='black',
                                            linewidth=5)
        for slip in near_cog_starts[1]:
            try:
                hili_sequences += ax_f2.plot((slip, slip), (0, 0.5),
                                             alpha=1,
                                             label=labels,
                                             zorder=4,
                                             color='black',
                                             linewidth=5)
            except Exception:
                hili_sequences = ax_f2.plot((slip, slip), (0, 0.5),
                                            alpha=1,
                                            label=labels,
                                            zorder=4,
                                            color='black',
                                            linewidth=5)
        for slip in near_cog_starts[2]:
            try:
                hili_sequences += ax_f3.plot((slip, slip), (0, 0.5),
                                             alpha=1,
                                             label=labels,
                                             zorder=4,
                                             color='black',
                                             linewidth=5)
            except Exception:
                hili_sequences = ax_f3.plot((slip, slip), (0, 0.5),
                                            alpha=1,
                                            label=labels,
                                            zorder=4,
                                            color='black',
                                            linewidth=5)

        # Plot sequence identifiers which will create a popup telling user what the subsequence is (useful if they have passed multiple subsequences)
        frame1_subsequences = ax_f1.plot(near_cog_starts[0],
                                         [0.25] * len(near_cog_starts[0]),
                                         'o',
                                         color='b',
                                         mec='k',
                                         ms=12,
                                         mew=1,
                                         alpha=0,
                                         zorder=4)
        frame2_subsequences = ax_f2.plot(near_cog_starts[1],
                                         [0.25] * len(near_cog_starts[1]),
                                         'o',
                                         color='b',
                                         mec='k',
                                         ms=12,
                                         mew=1,
                                         alpha=0,
                                         zorder=4)
        frame3_subsequences = ax_f3.plot(near_cog_starts[2],
                                         [0.25] * len(near_cog_starts[2]),
                                         'o',
                                         color='b',
                                         mec='k',
                                         ms=12,
                                         mew=1,
                                         alpha=0,
                                         zorder=4)

        # axisname.set_xticks([])
    for label in ax_main.xaxis.get_majorticklabels():
        label.set_fontsize(36)
    for axis, frame in ((ax_f1, 1), (ax_f2, 2), (ax_f3, 3)):
        axis.set_xlim(1, tranlen)
        starts = [(item, 1) for item in start_stop_dict[frame]['starts']]
        uag_stops = [(item, 1)
                     for item in start_stop_dict[frame]['stops']['TAG']]
        uaa_stops = [(item, 1)
                     for item in start_stop_dict[frame]['stops']['TAA']]
        uga_stops = [(item, 1)
                     for item in start_stop_dict[frame]['stops']['TGA']]
        # Plot start positions
        axis.broken_barh(starts, (0.30, 1),
                         color="white",
                         zorder=2,
                         linewidth=7)
        # Plot stop positions
        axis.broken_barh(uag_stops, (0, 1),
                         color=uag_col,
                         zorder=2,
                         linewidth=4)
        axis.broken_barh(uaa_stops, (0, 1),
                         color=uaa_col,
                         zorder=2,
                         linewidth=4)
        axis.broken_barh(uga_stops, (0, 1),
                         color=uga_col,
                         zorder=2,
                         linewidth=4)
        axis.set_ylim(0, 1)
        axis.set_ylabel('Frame {}'.format(frame),
                        labelpad=4,
                        verticalalignment='center',
                        horizontalalignment="right",
                        rotation="horizontal",
                        color="black",
                        fontsize=(axis_label_size / 1.5))
    title_str = '{} ({})'.format(gene, short_code)
    plt.title(title_str, fontsize=50, y=38)
    line_collections = [
        frame0subpro, frame1subpro, frame2subpro, rna_bars, allexons
    ]
    if mismatches:
        line_collections.append(a_mismatches)
        line_collections.append(t_mismatches)
        line_collections.append(g_mismatches)
        line_collections.append(c_mismatches)
    line_collections.append(cds_markers)

    if not (hili_start == 0 and hili_stop == 0):
        hili_start = int(hili_start)
        hili_stop = int(hili_stop)
        hili = ax_main.fill_between([hili_start, hili_stop], [y_max, y_max],
                                    zorder=0,
                                    alpha=0.75,
                                    color="#fffbaf")
        labels.append("Highligted region")
        start_visible.append(True)
        line_collections.append(hili)

    for alt_plot in alt_seq_type_vars:
        line_collections.append(alt_plot)
    if 'hili_sequences' in locals():
        labels.append("Highligted sequences")
        start_visible.append(True)
        line_collections.append(hili_sequences)
    if user_hili_starts != [] and user_hili_stops != []:
        for i in range(0, len(user_hili_starts)):
            user_hili_start = int(user_hili_starts[i])
            user_hili_stop = int(user_hili_stops[i])
            try:
                hili += ax_main.fill_between([user_hili_start, user_hili_stop],
                                             [y_max, y_max],
                                             alpha=0.75,
                                             color="#fffbaf")
            except Exception:
                hili = ax_main.fill_between([user_hili_start, user_hili_stop],
                                            [y_max, y_max],
                                            alpha=0.75,
                                            color="#fffbaf")
        labels.append("Highligter")
        start_visible.append(True)
        line_collections.append(hili)

    leg_offset = (legend_size - 17) * 5
    if leg_offset < 0:
        leg_offset = 0

    htmllabels = {1: [], 2: [], 3: []}
    all_start_points = {1: [], 2: [], 3: []}
    try:
        con_scores = SqliteDict(
            "{0}/{1}/homo_sapiens/score_dict.sqlite".format(
                config.SCRIPT_LOC, config.ANNOTATION_DIR))
    except Exception:
        con_scores = []
    for frame in [1, 2, 3]:
        orf_list = frame_orfs[frame]
        for tup in orf_list:
            orf_ribo = 0.0
            outframe_ribo = 0.0
            orf_rna = 0.0001
            start = tup[0]
            try:
                context = (seq[start - 7:start + 4].upper()).replace("T", "U")
            except Exception:
                con_score = "?"
            if len(context) != 11 or context[6:9] != "AUG":
                con_score = "?"
            else:
                try:
                    con_score = con_scores[context.upper()]
                except Exception:
                    con_score = "?"
            all_start_points[frame].append(start - 1)
            stop = tup[1]
            for i in range(start + 2, stop, 3):
                for subframe in [0, 1, 2]:
                    if i in frame_counts[subframe]:
                        orf_ribo += frame_counts[subframe][i]

            for i in range(start, stop, 3):
                for subframe in [0, 1, 2]:
                    if i in frame_counts[subframe]:
                        outframe_ribo += frame_counts[subframe][i]

            for i in range(start + 1, stop, 3):
                for subframe in [0, 1, 2]:
                    if i in frame_counts[subframe]:
                        outframe_ribo += frame_counts[subframe][i]

            for i in range(start, stop + 1):
                if i in all_rna_reads:
                    orf_rna += all_rna_reads[i]

            orf_te = float(orf_ribo) / float(orf_rna)
            orf_len = int(stop - start)

            try:
                in_out_ratio = orf_ribo / outframe_ribo
            except Exception:
                in_out_ratio = None

            datadict = {
                'inframe ribo': [orf_ribo],
                'outframe ribo': [outframe_ribo],
                'in/out ratio': [in_out_ratio],
                'rna': [orf_rna],
                'te': [orf_te],
                'len': [orf_len],
                'context_score': [str(con_score) + "/150"]
            }
            df = pd.DataFrame(datadict,
                              columns=([
                                  "inframe ribo", "outframe ribo",
                                  "in/out ratio", "rna", "te", "len",
                                  "context_score"
                              ]))
            label = df.iloc[[0], :].T
            label.columns = ["Start pos: {}".format(start - 1)]
            htmllabels[frame].append(str(label.to_html()))
    points1 = ax_f1.plot(all_start_points[1],
                         [0.75] * len(all_start_points[1]),
                         'o',
                         color='b',
                         mec='k',
                         ms=13,
                         mew=1,
                         alpha=0,
                         zorder=3)
    points2 = ax_f2.plot(all_start_points[2],
                         [0.75] * len(all_start_points[2]),
                         'o',
                         color='b',
                         mec='k',
                         ms=13,
                         mew=1,
                         alpha=0,
                         zorder=3)
    points3 = ax_f3.plot(all_start_points[3],
                         [0.75] * len(all_start_points[3]),
                         'o',
                         color='b',
                         mec='k',
                         ms=13,
                         mew=1,
                         alpha=0,
                         zorder=3)

    # This works but hides axis labels for all axes
    ax_main.tick_params('both', labelsize=marker_size)

    returnstr = "Position,Sequence,Frame 1,Frame 2,Frame 3,RNA-Seq"
    for seq_type in alt_seq_dict:
        returnstr += ",{}".format(seq_type)
    returnstr += "\n"
    for i in range(0, len(seq)):
        f1_count = 0
        f2_count = 0
        f3_count = 0
        rna_count = 0
        if i + 1 in frame_counts[0]:
            f1_count = frame_counts[0][i + 1]
        elif i + 1 in frame_counts[1]:
            f2_count = frame_counts[1][i + 1]
        elif i + 1 in frame_counts[2]:
            f3_count = frame_counts[2][i + 1]
        if i + 1 in all_rna_reads:
            rna_count = all_rna_reads[i + 1]
        returnstr += "{},{},{},{},{},{}".format(i + 1, seq[i], f1_count,
                                                f2_count, f3_count, rna_count)
        for seq_type in alt_seq_dict:
            if i + 1 in alt_seq_dict[seq_type]:
                count = alt_seq_dict[seq_type][i + 1]
            else:
                count = 0
            returnstr += ",{}".format(count)
        returnstr += "\n"

    ax_main.set_facecolor("white")
    # This changes the size of the tick markers, works on both firefox and chrome.

    for axisname in (ax_f1, ax_f2, ax_f3, ax_cds):
        axisname.tick_params(
            axis='both',
            which='both',
            top=False,
            bottom=False,
            left=False,
            right=False,
            labelbottom=False,
        )
        axisname.set_yticklabels([])

    plt.savefig('test.png', dpi=300,
                bbox_inches='tight')  # NOTE: What is it doing here?
    # ax_main.grid(False, color="white", linewidth=30,linestyle="solid")
    # Without this style tag the markers sizes will appear correct on browser but be original size when downloaded via png
    graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(
        marker_size)
    graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(
        short_code)
    graph += mpld3.fig_to_html(fig)
    return graph
