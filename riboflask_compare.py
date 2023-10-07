from fetch_shelve_reads2 import get_reads
import sqlite3
import os
import config

color_dict = {'frames': ['#FF4A45', '#64FC44', '#5687F9']}


def generate_compare_plot(
    tran: str,
    ambig: str,
    master_filepath_dict: dict,
    ribocoverage: bool,
    organism: str,
    normalize: bool,
    short_code: str,
    background_col: str,
    hili_start: str,
    hili_stop: str,
    comp_uag_col: str,
    comp_uga_col: str,
    comp_uaa_col: str,
    title_size: int,
    subheading_size: int,
    axis_label_size: int,
    marker_size: int,
    cds_marker_size: int,
    cds_marker_colour: int,
    legend_size: int,
    transcriptome: str,
) -> str:
    labels = []
    start_visible = []
    line_collections = []
    all_stops = ["TAG", "TAA", "TGA"]
    returnstr = "Position,"
    y_max = 0 if normalize else 50
    organisms = get_table("organisms")
    owner = organisms.loc[organisms.organism_name == organism
                          & organisms.transcriptome_list == transcriptome,
                          "owner"].values[0]
    if owner == 1:
        sqlfile = "{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC,
                                                      config.ANNOTATION_DIR,
                                                      organism, transcriptome)
        if not os.path.isfile(sqlfile):
            return_str = "Cannot find annotation file {}.{}.sqlite".format(
                organism, transcriptome)
            return {
                'current': 400,
                'total': 100,
                'status': 'return_str',
                'result': return_str
            }
    else:
        sqlfile = "{0}/transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
            config.UPLOADS_DIR, owner, organism, transcriptome)
    traninfo = get_table("transcripts", sqlfile)
    traninfo = traninfo[traninfo.transcript == tran].iloc[0]

    # traninfo = {
    # "transcript": result[0],
    # "gene": result[1],
    # "length": result[2],
    # "cds_start": result[3],
    # "cds_stop": result[4],
    # "seq": result[5],
    # "strand": result[6],
    # "stop_list": result[7].split(","),
    # "start_list": result[8].split(","),
    # "exon_junctions": result[9].split(","),
    # "tran_type": result[10],
    # "principal": result[11]
    # }
    traninfo["stop_list"] = [int(x) for x in traninfo["stop_list"]]
    traninfo["start_list"] = [int(x) for x in traninfo["start_list"]]
    if traninfo["exon_junctions"]:  # NOTE: Expecting list here
        traninfo["exon_junctions"] = [
            int(x) for x in traninfo["exon_junctions"]
        ]
    else:
        traninfo["exon_junctions"] = []

    gene = traninfo["gene"]
    tranlen = traninfo["length"]
    cds_start = traninfo["cds_start"]
    cds_stop = traninfo["cds_stop"]

    if not cds_stop:
        cds_start = 0
        cds_stop = 0

    all_starts = traninfo["start_list"]
    all_stops = {"TAG": [], "TAA": [], "TGA": []}
    seq = traninfo["seq"].upper()
    for i in range(0, len(seq)):
        if seq[i:i + 3] in all_stops:
            all_stops[seq[i:i + 3]].append(i + 1)
    start_stop_dict = {}
    for frame in [1, 2, 3]:
        start_stop_dict[frame] = {
            "starts": [0] , # Why do we need zero?
            "stops": {
                "TGA": [0],
                "TAG": [0],
                "TAA": [0]
            }
        }
    for start in all_starts:
        rem = ((start - 1) % 3) + 1
        start_stop_dict[rem]["starts"].append(start - 1)
    for stop in all_stops:
        for stop_pos in all_stops[stop]:
            rem = ((stop_pos - 1) % 3) + 1
            start_stop_dict[rem]["stops"][stop].append(stop_pos - 1)

    fig = plt.figure(figsize=(13, 8))
    ax_main = plt.subplot2grid((30, 1), (0, 0), rowspan=22)
    label = 'Read count' if not normalize else  'Normalized read count'
    ax_main.set_ylabel(label, fontsize=axis_label_size, labelpad=30)
    label = 'Position (nucleotides)'
    ax_main.set_xlabel(label, fontsize=axis_label_size, labelpad=10)

    #if normalize is true work out the factors for each colour
    if normalize:
        all_mapped_reads = []
        for color in master_filepath_dict:
            all_mapped_reads.append(
                    master_filepath_dict[color]["mapped_reads"]) # NOTE: if 'mapped reads would be fixed it would faster'
        for color in master_filepath_dict:
            factor =min(all_mapped_reads) *1. / master_filepath_dict[color]["mapped_reads"])
            master_filepath_dict[color]["factor"] = factor

    # So items can be plotted alphabetically
    unsorted_list = []
    for color in master_filepath_dict:
        input_list = [
            color, master_filepath_dict[color]["file_names"],
            master_filepath_dict[color]["file_descs"],
            master_filepath_dict[color]["file_ids"],
            master_filepath_dict[color]["filepaths"],
            master_filepath_dict[color]["file_type"],
            master_filepath_dict[color]["minread"],
            master_filepath_dict[color]["maxread"]
        ]
        if "factor" in master_filepath_dict[color]:
            input_list.append(master_filepath_dict[color]["factor"])
        unsorted_list.append(input_list)

    sorted_list = sorted(unsorted_list, key=lambda x: x[1][0])
    returndict = {}
    for item in sorted_list:
        # needed to make get_reads accept file_paths
        file_paths = {"riboseq": {}}
        for i in range(0, len(item[3])):
            file_paths["riboseq"][item[3][i]] = item[4][i]
        file_names = item[1][0]
        if item[5] == "riboseq":
            filename_reads, _ = get_reads(ambig, item[6], item[7], tran,
                                          file_paths, tranlen, ribocoverage,
                                          organism, False, False, "fiveprime",
                                          "riboseq", 1)
        else:
            filename_reads, _ = get_reads(ambig, item[6], item[7], tran,
                                          file_paths, tranlen, True, organism,
                                          False, False, "fiveprime", "riboseq",
                                          1)
        if not normalize:
            try:
                max_val = max(filename_reads.values()) * 1.1
                if max_val > y_max:
                    y_max = max_val
            except Exception:
                pass
            labels.append(file_names)
            start_visible.append(True)
            plot_filename = ax_main.plot(filename_reads.keys(),
                                         filename_reads.values(),
                                         alpha=1,
                                         label=labels,
                                         zorder=1,
                                         color=item[0],
                                         linewidth=3)
            line_collections.append(plot_filename)
            returndict[file_names] = {}
            for pos in filename_reads:
                returndict[file_names][pos] = filename_reads[pos]

        else:
            normalized_reads = {}
            for pos in filename_reads:
                normalized_reads[pos] = filename_reads[pos] * item[8]
            try:
                max_val = max(normalized_reads.values()) * 1.1
                if max_val > y_max:
                    y_max = max_val
            except Exception:
                pass
            labels.append(file_names)
            start_visible.append(True)
            plot_filename = ax_main.plot(normalized_reads.keys(),
                                         normalized_reads.values(),
                                         alpha=1,
                                         label=labels,
                                         zorder=1,
                                         color=item[0],
                                         linewidth=3)

            line_collections.append(plot_filename)
            returndict[file_names] = {}
            for pos in filename_reads:
                returndict[file_names][pos] = normalized_reads[pos]

    for plot_filename in returndict:
        returnstr += "{},".format(plot_filename)
    returnstr += "\n"

    for i in range(0, tranlen):
        returnstr += "{},".format(i)
        for plot_filename in returndict:
            returnstr += "{},".format(returndict[plot_filename][i])
        returnstr += "\n"

    ax_main.set_ylim(0, y_max)
    # draw cds start

    # draw cds end

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
    line_collections.append(cds_markers)
    start_visible.append(True)
    labels.append("CDS Markers")
    ax_f1 = plt.subplot2grid((30, 1), (27, 0), rowspan=1, sharex=ax_main)
    ax_f1.set_facecolor('lightgray')
    ax_f2 = plt.subplot2grid((30, 1), (28, 0), rowspan=1, sharex=ax_main)
    ax_f2.set_facecolor('lightgray')
    ax_f6 = plt.subplot2grid((30, 1), (29, 0), rowspan=1, sharex=ax_main)
    ax_f6.set_facecolor('lightgray')
    ax_f6.set_xlabel('Transcript: {}   Length: {} nt'.format(tran, tranlen),
                     fontsize=subheading_size,
                     labelpad=10)

    for axis, frame in ((ax_f1, 1), (ax_f2, 2), (ax_f6, 3)):
        color = color_dict['frames'][frame - 1]
        axis.set_xlim(0, tranlen)
        starts = [(item, 1) for item in start_stop_dict[frame]['starts']]
        axis.broken_barh(starts, (0.5, 1),
                         color='white',
                         zorder=5,
                         linewidth=2)
        uag_stops = [(item, 1)
                     for item in start_stop_dict[frame]['stops']['TAG']]
        uaa_stops = [(item, 1)
                     for item in start_stop_dict[frame]['stops']['TAA']]
        uga_stops = [(item, 1)
                     for item in start_stop_dict[frame]['stops']['TGA']]
        axis.broken_barh(uag_stops, (0, 1),
                         color=comp_uag_col,
                         zorder=2,
                         linewidth=2)
        axis.broken_barh(uaa_stops, (0, 1),
                         color=comp_uaa_col,
                         zorder=2,
                         linewidth=2)
        axis.broken_barh(uga_stops, (0, 1),
                         color=comp_uga_col,
                         zorder=2,
                         linewidth=2)
        axis.set_ylabel('{}'.format(frame),
                        rotation='horizontal',
                        labelpad=10,
                        verticalalignment='center')
        axis.set_ylim(0, 1)
        axis.tick_params(top=False,
                         left=False,
                         right=False,
                         bottom=False,
                         labeltop=False,
                         labelleft=False,
                         labelright=False,
                         labelbottom=False)
    ax_f6.axes.get_yaxis().set_ticks([])
    ax_f2.axes.get_yaxis().set_ticks([])
    ax_f1.axes.get_yaxis().set_ticks([])
    title_str = '{} ({})'.format(gene, short_code)
    plt.title(title_str, fontsize=title_size, y=36)

    if not (hili_start == 0 and hili_stop == 0):
        hili_start = int(hili_start)
        hili_stop = int(hili_stop)
        hili = ax_main.fill_between([hili_start, hili_stop], [y_max, y_max],
                                    zorder=0,
                                    alpha=0.75,
                                    color="#fffbaf")
        labels.append("Highligter")
        start_visible.append(True)
        line_collections.append(hili)

    leg_offset = (legend_size - 17) * 5
    if leg_offset < 0:
        leg_offset = 0
    leg_offset += 230
    ilp = InteractiveLegendPlugin(line_collections,
                                  labels,
                                  alpha_unsel=0,
                                  alpha_sel=0.85,
                                  xoffset=leg_offset,
                                  yoffset=20,
                                  start_visible=start_visible,
                                  fontsize=legend_size)
    plugins.connect(fig, ilp, TopToolbar(yoffset=-50, xoffset=-300),
                    DownloadProfile(returnstr=returnstr),
                    DownloadPNG(returnstr=title_str))
    ax_main.set_facecolor(background_col)
    # This changes the size of the tick markers, works on both firefox and chrome.
    ax_main.tick_params('both', labelsize=marker_size)
    ax_main.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax_main.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax_main.grid(color="white", linewidth=20, linestyle="solid")
    graph = "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(
        short_code)
    graph += mpld3.fig_to_html(fig)
    return graph
