from fetch_shelve_reads2 import get_reads
from sqlqueries_2 import get_table, sqlquery
import polars as pl
import os
import config


def generate_compare_plot(
    data,
) -> str | dict:
    """

    Parameters:

    Example:
    """
    start_visible = []
    line_collections = []
    all_stops = ["TAG", "TAA", "TGA"]
    returnstr = "Position,"
    y_max = 0 if data['normalize'] else 50
    owner = get_table("organisms").filter(
        (pl.col("organism_name") == data['organism'])
        & (pl.col("transcriptome_list") == data['transcriptome']))[0, 'owner']
    if owner == 1:
        sqlfile = "{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC,
                                                      config.ANNOTATION_DIR,
                                                      data['organism'],
                                                      data['transcriptome'])
        if not os.path.isfile(sqlfile):
            return_str = "Cannot find annotation file {}.{}.sqlite".format(
                data['organism'], data['transcriptome'])
            return {
                'current': 400,
                'total': 100,
                'status': 'return_str',
                'result': return_str
            }
    else:
        sqlfile = "{0}/transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
            config.UPLOADS_DIR, owner, data['organism'], data['transcriptome'])
    traninfo = sqlquery(sqlfile, "transcripts").filter(
        pl.col("transcript") == data['transcript'])[0].with_columns(
            seq=pl.col("seq").str.to_uppercase()
    )

    for i in range(0, len(traninfo[0, 'seq'])):
        codon = traninfo[0, 'seq'][i:i + 3]
        if codon in all_stops:
            all_stops[codon].append(i + 1)
    start_stop_dict = {}
    for frame in [1, 2, 3]:
        start_stop_dict[frame] = {
            "starts": [0],  # Why do we need zero?
            "stops": {
                "TGA": [0],
                "TAG": [0],
                "TAA": [0]
            }
        }
    for start in all_starts: 
        rem = start % 3
        rem = rem if rem else 3
        start_stop_dict[rem]["starts"].append(start - 1)
    for stop in all_stops:
        for stop_pos in all_stops[stop]:
            rem = stop_pos % 3
            rem = rem if rem else 3
            start_stop_dict[rem]["stops"][stop].append(stop_pos - 1)

    label = 'Read count' if 'normalize' not in data else 'Normalized read count'

    # if normalize is true work out the factors for each colour
    if data['normalize']:
        all_mapped_reads = []
        for color in data['master_filepath_dict']:
            all_mapped_reads.append(
                data['master_filepath_dict'][color]["mapped_reads"]
            )  # NOTE: if 'mapped reads would be fixed it would faster'
        for color in data['master_filepath_dict']:
            data['master_filepath_dict'][color]["factor"] = (
                min(all_mapped_reads) * 1. /
                data['master_filepath_dict'][color]["mapped_reads"])

    # So items can be plotted alphabetically
    unsorted_list = []
    for color in data['master_filepath_dict']:
        input_list = [
            color, data['master_filepath_dict'][color]["file_names"],
            data['master_filepath_dict'][color]["file_descs"],
            data['master_filepath_dict'][color]["file_ids"],
            data['master_filepath_dict'][color]["filepaths"],
            data['master_filepath_dict'][color]["file_type"],
            data['master_filepath_dict'][color]["minread"],
            data['master_filepath_dict'][color]["maxread"]
        ]
        if "factor" in data['master_filepath_dict'][color]:
            input_list.append(data['master_filepath_dict'][color]["factor"])
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
            filename_reads, _ = get_reads(data['ambigious'], item[6], item[7], data['transcript'],
                                          file_paths, tranlen, data['ribocoverage'],
                                          data['organism'], False, False,
                                          "fiveprime", "riboseq", 1)
        else:
            filename_reads, _ = get_reads(data['ambigious'], item[6], item[7], data['transcript'],
                                          file_paths, tranlen, True,
                                          data['organism'], False, False,
                                          "fiveprime", "riboseq", 1)
        if not data['normalize']:
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
                               color=data['cds_marker_colour'],
                               linestyle='solid',
                               linewidth=data['cds_marker_size'])
    cds_markers += ax_main.plot((cds_stop + 1, cds_stop + 1),
                                (0, y_max * 0.97),
                                color=data['cds_marker_colour'],
                                linestyle='solid',
                                linewidth=data['cds_marker_size'])

    for axis, frame in ((ax_f1, 1), (ax_f2, 2), (ax_f6, 3)):
        color = color_dict['frames'][frame - 1]
        axis.set_xlim(0, tranlen)
        starts = [(item, 1) for item in start_stop_dict[frame]['starts']]
    title_str = '{} ({})'.format(gene, data['short'])

                                alpha = 0.75,
                                color = "#fffbaf")
    line_collections.append(hili)

    leg_offset = (data["legend_size"] - 17) * 5
    if leg_offset < 0:
        leg_offset = 0
    leg_offset += 230
    reurn plot
