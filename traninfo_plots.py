from typing import List, Tuple, Union, Dict
import config
import altair as alt
from plots import VegaPlot
import collections
import polars as pl
import polars.selectors as cs
import numpy as np

#ViennaRNA can be installed from here https://github.com/ViennaRNA/ViennaRNA
try:
    import RNA
    vienna_rna = True
except Exception:
    print(
        "Could not import the RNA module, ViennaRNA needs to be installed (https://github.com/ViennaRNA/ViennaRNA), MFE will not be plotted on traninfo plot"
    )
    vienna_rna = False

redhex = "#FF5F5B"
greenhex = "#90E090"
bluehex = "#9ACAFF"
yellowhex = "#FFFF91"


def nuc_freq_plot(master_dict: Dict[str, Dict[str, str]], title: str,
                  short_code: str, background_col: str, readlength_col: str,
                  title_size: int, axis_label_size: int, subheading_size: int,
                  marker_size: int, filename: str) -> str:
    """

    Parameters:
    - master_dict (Dict[str, Dict[str, str]]): master dictionary
    - title (str): title
    - short_code (str): short code
    - background_col (str): background color
    - readlength_col (str): readlength color
    - title_size (int): title size
    - axis_label_size (int): axis label size
    - subheading_size (int): subheading size
    - marker_size (int): marker size
    - filename (str): filename

    Returns:


    Example:
    """
    labels = ["A", "T", "G", "C"]
    returnstr = "Position,A,T,G,C\n"
    minpos = min(master_dict.keys())
    maxpos = max(master_dict.keys())
    x_pos = []
    a_counts = []
    t_counts = []
    g_counts = []
    c_counts = []
    for i in range(minpos, maxpos):
        returnstr += "{},{:.2f},{:.2f},{:.2f},{:.2f}\n".format(
            i, master_dict[i]["A"], master_dict[i]["T"], master_dict[i]["G"],
            master_dict[i]["C"])
        x_pos.append(i)
        a_counts.append(master_dict[i]["A"])
        t_counts.append(master_dict[i]["T"])
        g_counts.append(master_dict[i]["G"])
        c_counts.append(master_dict[i]["C"])

    ax.set_xlabel('Position (nucleotides)', fontsize=axis_label_size)

    ax = plt.subplot(111)
    title_str = "{} ({})".format(title, short_code)
    ax.set_title(title_str, y=1.05, fontsize=title_size)
    a_line = ax.plot(x_pos, a_counts, label=labels, color="blue", linewidth=4)
    t_line = ax.plot(x_pos, t_counts, label=labels, color="red", linewidth=4)
    g_line = ax.plot(x_pos, g_counts, label=labels, color="green", linewidth=4)
    c_line = ax.plot(x_pos,
                     c_counts,
                     label=labels,
                     color="orange",
                     linewidth=4)
    ax.set_facecolor(background_col)
    ax.tick_params('both', labelsize=marker_size)
    plt.grid(color="white", linewidth=2, linestyle="solid")


def nuc_comp_single(data: dict):
    """

    Parameters:
    - tran (str): transcript
    - master_dict (dict): master dictionary
    - title (str): title
    - short_code (str): short code
    - background_col (str): background color
    - readlength_col (str): read length color
    - title_size (int): title size
    - axis_label_size (int): axis label size
    - subheading_size (int): subheading size
    - marker_size (int): marker size
    - traninfo (dict): transcript information

    Returns:

    Example:

    """
    step_size = 2
    window_size = 60
    nucleotide_content = []
    t_counts = [0, 0, 0, 0]
    seq = data[0, 'sequence']
    for nuc in seq[:window_size]:
        if nuc == "A":
            t_counts[0] += 1
        elif nuc == "C":
            t_counts[1] += 1
        elif nuc == "G":
            t_counts[2] += 1
        elif nuc == "T":
            t_counts[3] += 1
    nucleotide_content.append(t_counts.copy())
    for i in range(step_size, data[0, "length"] - (window_size), step_size):
        # TODO: optimise this
        for nuc in seq[i:i + step_size]:
            if nuc == "A":
                t_counts[0] -= 1
            elif nuc == "C":
                t_counts[1] -= 1
            elif nuc == "G":
                t_counts[2] -= 1
            elif nuc == "T":
                t_counts[3] -= 1
        for nuc in seq[i + window_size - step_size:i + window_size]:
            if nuc == "A":
                t_counts[0] += 1
            elif nuc == "C":
                t_counts[1] += 1
            elif nuc == "G":
                t_counts[2] += 1
            elif nuc == "T":
                t_counts[3] += 1
        nucleotide_content.append(t_counts.copy())
    nucleotide_content = pl.DataFrame(
        nucleotide_content, schema=[
            "A", "C", "G", "T"
        ]).with_columns((pl.col("G") + pl.col("C")).alias("GC")).with_columns(
            pl.all() * 100 / window_size  # 100 in the max based on old code
        ).with_columns(
            pos=pl.arange(0, len(seq) - window_size, step_size, eager=True) +
            (window_size / 2)).melt(id_vars="pos",
                                    value_vars=["A", "C", "G", "T", "GC"],
                                    variable_name="frame",
                                    value_name="count")  #.collect()
    colors = alt.Scale(domain=["A", "C", "G", "T", "GC"],
                       range=config.BOX_COLORS[:4])
    plot = VegaPlot(nucleotide_content, colors)
    return plot.line("pos", "count").to_json()
    print(nucleotide_content)

    plot_mfe = True
    if plot_mfe and vienna_rna:
        step_size = 2
        window_size = 60
        mfe_dict = collections.OrderedDict()
        for i in range(0, len(seq) - (window_size), step_size):
            seq_window = str(seq[i:i + window_size])
            (ss, mfe) = RNA.fold(seq_window)
            mfe_dict[i + (window_size / 2)] = abs(mfe)
    else:
        mfe_dict = {}


def gc_metagene(title: str, short_code: str, background_col: str,
                readlength_col: str, title_size: int, axis_label_size: int,
                subheading_size: int, marker_size: int, traninfo: str) -> str:
    """

    Parameters:
    - title (str): title
    - short_code (str): short code
    - background_col (str): background color
    - readlength_col (str): readlength color
    - title_size (int): title size
    - axis_label_size (int): axis label size
    - subheading_size (int): subheading size
    - marker_size (int): marker size
    - traninfo (str): traninfo

    Returns:

    Example:
    """
    labels = ["CDS markers"]
    start_visible = [True]
    color_dict = {'frames': ['#FF4A45', '#64FC44', '#5687F9']}
    gene = ""
    y_max = 100
    fig = plt.figure(figsize=(13, 8))

    ax_main = plt.subplot2grid((30, 1), (0, 0), rowspan=22)
    ax_main.spines['bottom'].set_visible(False)
    ax_main.set_ylabel("%", fontsize=axis_label_size, labelpad=30)
    ax_main.set_ylim(0, y_max)
    ax_main.set_xlim(0, 1500)
    cds_markers = ax_main.plot((500, 500), (0, y_max - 3),
                               color="black",
                               linestyle='solid',
                               linewidth=2)
    cds_markers += ax_main.plot((1000, 1000), (0, y_max - 3),
                                color="black",
                                linestyle='solid',
                                linewidth=2)
    for label in ax_main.xaxis.get_majorticklabels():
        label.set_fontsize(36)

    title_str = '{} ({})'.format(gene, short_code)
    plt.title(title_str, fontsize=title_size, y=36)
    line_collections = [cds_markers]

    plot_gc = True
    plot_mfe = False
    if plot_mfe and vienna_rna:
        step_size = 10
        window_size = 60
        mfe_dict = collections.OrderedDict()
        for item in traninfo:
            transcript = item[0]
            cds_start = float(item[1])
            cds_stop = float(item[2])
            seq = item[3]
            seqlen = len(seq)
            for i in range(0, len(seq) - (window_size), step_size):
                seq_window = str(seq[i:i + window_size])
                (ss, mfe) = RNA.fold(seq_window)
                if i < cds_start:
                    per = (i + (window_size / 2) / cds_start) * 5
                if i >= cds_start and i <= cds_stop:
                    per = 500 + ((i + (window_size / 2) /
                                  (cds_stop - cds_start)) * 5)
                if i > cds_stop:
                    per = 1000 + ((i + (window_size / 2) /
                                   (seqlen - cds_stop)) * 5)
                if per not in mfe_dict:
                    mfe_dict[per] = [abs(mfe)]
                else:
                    mfe_dict[per].append(abs(mfe))
        for per in mfe_dict:
            mfe_dict[per] = sum(mfe_dict[per]) / len(mfe_dict[per])

    if plot_gc:
        step_size = 10
        window_size = 60
        a_dict = {}
        t_dict = {}
        g_dict = {}
        c_dict = {}
        gc_dict = {}
        sorted_gc_dict = collections.OrderedDict()
        for item in traninfo:
            cds_start = float(item[1])
            cds_stop = float(item[2])
            seq = item[3]
            seqlen = len(seq)
            for i in range(0, len(seq) - (window_size), step_size):
                mid_window = i + (window_size / 2)
                a_count = 0.0
                t_count = 0.0
                g_count = 0.0
                c_count = 0.0
                for x in range(i, i + window_size):
                    if seq[x] == "A":
                        a_count += 1
                    elif seq[x] == "T":
                        t_count += 1
                    elif seq[x] == "G":
                        g_count += 1
                    elif seq[x] == "C":
                        c_count += 1

                gc_count = g_count + c_count
                norm_a = a_count / window_size
                norm_t = t_count / window_size
                norm_g = g_count / window_size
                norm_c = c_count / window_size
                norm_gc = gc_count / window_size

                final_a = norm_a * y_max
                final_t = norm_t * y_max
                final_g = norm_g * y_max
                final_c = norm_c * y_max
                final_gc = norm_gc * y_max
                if mid_window < cds_start:
                    per = mid_window / cds_start
                    per = per * 500
                    per = int(per)
                if mid_window >= cds_start and mid_window <= cds_stop:
                    per = (mid_window - cds_start) / (cds_stop - cds_start)
                    per = per * 500
                    per += 500
                    per = int(per)
                if mid_window > cds_stop:
                    per = (mid_window - cds_stop) / (seqlen - cds_stop)
                    per = per * 500
                    per += 1000
                    per = int(per)
                if per not in a_dict:
                    a_dict[per] = [final_a]
                else:
                    a_dict[per].append(final_a)

                if per not in t_dict:
                    t_dict[per] = [final_t]
                else:
                    t_dict[per].append(final_t)

                if per not in g_dict:
                    g_dict[per] = [final_g]
                else:
                    g_dict[per].append(final_g)

                if per not in c_dict:
                    c_dict[per] = [final_c]
                else:
                    c_dict[per].append(final_c)

                if per not in gc_dict:
                    gc_dict[per] = [final_gc]
                else:
                    gc_dict[per].append(final_gc)

        for per in a_dict:
            a_dict[per] = sum(a_dict[per]) / len(a_dict[per])
        for per in t_dict:
            t_dict[per] = sum(t_dict[per]) / len(t_dict[per])
        for per in g_dict:
            g_dict[per] = sum(g_dict[per]) / len(g_dict[per])
        for per in c_dict:
            c_dict[per] = sum(c_dict[per]) / len(c_dict[per])
        for per in sorted(gc_dict.keys()):
            sorted_gc_dict[per] = sum(gc_dict[per]) / len(gc_dict[per])

        a_plot = ax_main.plot(a_dict.keys(),
                              a_dict.values(),
                              alpha=0.01,
                              label=labels,
                              zorder=1,
                              color=color_dict['frames'][0],
                              linewidth=4)
        t_plot = ax_main.plot(t_dict.keys(),
                              t_dict.values(),
                              alpha=0.01,
                              label=labels,
                              zorder=1,
                              color=color_dict['frames'][1],
                              linewidth=4)
        g_plot = ax_main.plot(g_dict.keys(),
                              g_dict.values(),
                              alpha=0.01,
                              label=labels,
                              zorder=1,
                              color=color_dict['frames'][2],
                              linewidth=4)
        c_plot = ax_main.plot(c_dict.keys(),
                              c_dict.values(),
                              alpha=0.01,
                              label=labels,
                              zorder=1,
                              color='#ffff99',
                              linewidth=4)
        gc_plot = ax_main.plot(sorted_gc_dict.keys(),
                               sorted_gc_dict.values(),
                               alpha=1,
                               label=labels,
                               zorder=1,
                               color='grey',
                               linewidth=4)
        if plot_mfe:
            ax_main.plot(mfe_dict.keys(),
                         mfe_dict.values(),
                         alpha=0.01,
                         label=labels,
                         zorder=1,
                         color='#df8500',
                         linewidth=4)
        for item, lbl, viz in [(a_plot, "A%", False), (t_plot, "T%", False),
                               (g_plot, "G%", False), (c_plot, "C%", False),
                               (gc_plot, "GC%", True)]:
            line_collections.append(item)
            labels.append(lbl)
            start_visible.append(viz)

    leg_offset = (30 - 17) * 5
    if leg_offset < 0:
        leg_offset = 0

    ax_main.set_facecolor(background_col)
    # This changes the size of the tick markers, works on both firefox and chrome.
    ax_main.tick_params('both', labelsize=marker_size)
    ax_main.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax_main.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax_main.grid(True, color="white", linewidth=30, linestyle="solid")
    ax_main.text(500,
                 y_max * 0.97,
                 "CDS start",
                 fontsize=18,
                 color="black",
                 ha="center")
    ax_main.text(1000,
                 y_max * 0.97,
                 "CDS stop",
                 fontsize=18,
                 color="black",
                 ha="center")
    #hide x axis set_ticks


def nuc_comp_scatter(master_dict, filename, title_size, axis_label_size,
                     marker_size, nucleotide, short_code):
    """

    Parameters:
    - master_dict (Dict[str, Dict[str, str]]): master dictionary
    - filename (str): filename
    - title_size (int): title size
    - axis_label_size (int): axis label size
    - marker_size (int): marker size
    - nucleotide (str): nucleotide
    - short_code (str): short code

    Returns:

    Example:
    """
    x_values = []
    gc_list = master_dict[1]["lengths"]
    tran_list = master_dict[1]["trans"]
    for i in range(1, len(gc_list) + 1):
        x_values.append(i)
    source = ColumnDataSource({
        'x': x_values,
        'y': gc_list,
        'trans': tran_list
    })
    x_num = len(gc_list)

    x_values2 = []
    gc_list2 = master_dict[2]["lengths"]
    tran_list2 = master_dict[2]["trans"]
    for i in range(1, len(gc_list2) + 1):
        x_values2.append(x_num + i)
    source2 = ColumnDataSource({
        'x': x_values2,
        'y': gc_list2,
        'trans': tran_list2
    })
    x_num += len(gc_list2)

    x_values3 = []
    gc_list3 = master_dict[3]["lengths"]
    tran_list3 = master_dict[3]["trans"]
    for i in range(1, len(gc_list3) + 1):
        x_values3.append(x_num + i)
    source3 = ColumnDataSource({
        'x': x_values3,
        'y': gc_list3,
        'trans': tran_list3
    })
    x_num += len(gc_list3)

    x_values4 = []
    gc_list4 = master_dict[4]["lengths"]
    tran_list4 = master_dict[4]["trans"]
    for i in range(1, len(gc_list4) + 1):
        x_values4.append(x_num + i)
    source4 = ColumnDataSource({
        'x': x_values4,
        'y': gc_list4,
        'trans': tran_list4
    })
    x_num += len(gc_list4)

    p.scatter('x',
              'y',
              alpha=0.2,
              color="black",
              fill_alpha=1,
              size=12,
              source=source,
              fill_color='green')
    p.scatter('x',
              'y',
              alpha=0.2,
              color="black",
              fill_alpha=1,
              size=12,
              source=source2,
              fill_color='red')
    p.scatter('x',
              'y',
              alpha=0.2,
              color="black",
              fill_alpha=1,
              size=12,
              source=source3,
              fill_color='blue')
    p.scatter('x',
              'y',
              alpha=0.2,
              color="black",
              fill_alpha=1,
              size=12,
              source=source4,
              fill_color='yellow')
    hover = p.select(dict(type=HoverTool))


def lengths_scatter(master_dict, filename, title_size, axis_label_size,
                    marker_size, short_code):
    """

    Parameters:
    - master_dict (Dict[str, Dict[str, str]]): master dictionary
    - filename (str): filename
    - title_size (int): title size
    - axis_label_size (int): axis label size
    - marker_size (int): marker size
    - short_code (str): short code

    Returns:

    Example:
    """
    x_values = []
    gc_list = master_dict[1]["lengths"]
    tran_list = master_dict[1]["trans"]
    for i in range(1, len(gc_list) + 1):
        x_values.append(i)
    source = ColumnDataSource({
        'x': x_values,
        'y': gc_list,
        'trans': tran_list
    })
    x_num = len(gc_list)

    x_values2 = []
    gc_list2 = master_dict[2]["lengths"]
    tran_list2 = master_dict[2]["trans"]
    for i in range(1, len(gc_list2) + 1):
        x_values2.append(x_num + i)
    source2 = ColumnDataSource({
        'x': x_values2,
        'y': gc_list2,
        'trans': tran_list2
    })
    x_num += len(gc_list2)

    x_values3 = []
    gc_list3 = master_dict[3]["lengths"]
    tran_list3 = master_dict[3]["trans"]
    for i in range(1, len(gc_list3) + 1):
        x_values3.append(x_num + i)
    source3 = ColumnDataSource({
        'x': x_values3,
        'y': gc_list3,
        'trans': tran_list3
    })
    x_num += len(gc_list3)

    x_values4 = []
    gc_list4 = master_dict[4]["lengths"]
    tran_list4 = master_dict[4]["trans"]
    for i in range(1, len(gc_list4) + 1):
        x_values4.append(x_num + i)
    source4 = ColumnDataSource({
        'x': x_values4,
        'y': gc_list4,
        'trans': tran_list4
    })
    x_num += len(gc_list4)
    full_title = "Lengths ({})".format(short_code)
    y_lab = "Length"

    p.scatter('x',
              'y',
              alpha=0.2,
              color="black",
              fill_alpha=1,
              size=12,
              source=source,
              fill_color='green')
    p.scatter('x',
              'y',
              alpha=0.2,
              color="black",
              fill_alpha=1,
              size=12,
              source=source2,
              fill_color='red')
    p.scatter('x',
              'y',
              alpha=0.2,
              color="black",
              fill_alpha=1,
              size=12,
              source=source3,
              fill_color='blue')
    p.scatter('x',
              'y',
              alpha=0.2,
              color="black",
              fill_alpha=1,
              size=12,
              source=source4,
              fill_color='yellow')


def nuc_comp_box(master_dict, filename, nucleotide, title_size, box_colour,
                 axis_label_size, marker_size, short_code):
    """

    Parameters:
    - master_dict (Dict[str, Dict[str, str]]): master dictionary
    - filename (str): filename
    - nucleotide (str): nucleotide
    - title_size (int): title size
    - box_colour (str): box colour
    - axis_label_size (int): axis label size
    - marker_size (int): marker size
    - short_code (str): short code

    Returns:

    Example:
    """
    gc_list = master_dict[1]["gc"]
    gc_list2 = master_dict[2]["gc"]
    gc_list3 = master_dict[3]["gc"]
    gc_list4 = master_dict[4]["gc"]

    a_list = master_dict[1]["a"]
    a_list2 = master_dict[2]["a"]
    a_list3 = master_dict[3]["a"]
    a_list4 = master_dict[4]["a"]

    t_list = master_dict[1]["t"]
    t_list2 = master_dict[2]["t"]
    t_list3 = master_dict[3]["t"]
    t_list4 = master_dict[4]["t"]

    g_list = master_dict[1]["g"]
    g_list2 = master_dict[2]["g"]
    g_list3 = master_dict[3]["g"]
    g_list4 = master_dict[4]["g"]

    c_list = master_dict[1]["c"]
    c_list2 = master_dict[2]["c"]
    c_list3 = master_dict[3]["c"]
    c_list4 = master_dict[4]["c"]

    cats = []
    grouplist = []
    gclist = []

    if nucleotide == "A":
        cats.append("Group 1")
        for item in a_list:
            grouplist.append("Group 1")
            gclist.append(item)
        if a_list2 != []:
            cats.append("Group 2")
            for item in a_list2:
                grouplist.append("Group 2")
                gclist.append(item)
        if a_list3 != []:
            cats.append("Group 3")
            for item in a_list3:
                grouplist.append("Group 3")
                gclist.append(item)
        if a_list4 != []:
            cats.append("Group 4")
            for item in a_list4:
                grouplist.append("Group 4")
                gclist.append(item)
    if nucleotide == "C":
        cats.append("Group 1")
        for item in c_list:
            grouplist.append("Group 1")
            gclist.append(item)
        if c_list2 != []:
            cats.append("Group 2")
            for item in c_list2:
                grouplist.append("Group 2")
                gclist.append(item)
        if c_list3 != []:
            cats.append("Group 3")
            for item in c_list3:
                grouplist.append("Group 3")
                gclist.append(item)
        if c_list4 != []:
            cats.append("Group 4")
            for item in c_list4:
                grouplist.append("Group 4")
                gclist.append(item)
    if nucleotide == "GC":
        cats.append("Group 1")
        for item in gc_list:
            grouplist.append("Group 1")
            gclist.append(item)
        if gc_list2 != []:
            cats.append("Group 2")
            for item in gc_list2:
                grouplist.append("Group 2")
                gclist.append(item)
        if gc_list3 != []:
            cats.append("Group 3")
            for item in gc_list3:
                grouplist.append("Group 3")
                gclist.append(item)
        if gc_list4 != []:
            cats.append("Group 4")
            for item in gc_list4:
                grouplist.append("Group 4")
                gclist.append(item)
    if nucleotide == "G":
        cats.append("Group 1")
        for item in g_list:
            grouplist.append("Group 1")
            gclist.append(item)
        if g_list2 != []:
            cats.append("Group 2")
            for item in g_list2:
                grouplist.append("Group 2")
                gclist.append(item)
        if g_list3 != []:
            cats.append("Group 3")
            for item in g_list3:
                grouplist.append("Group 3")
                gclist.append(item)
        if g_list4 != []:
            cats.append("Group 4")
            for item in g_list4:
                grouplist.append("Group 4")
                gclist.append(item)
    if nucleotide == "T":
        cats.append("Group 1")
        for item in t_list:
            grouplist.append("Group 1")
            gclist.append(item)
        if t_list2 != []:
            cats.append("Group 2")
            for item in t_list2:
                grouplist.append("Group 2")
                gclist.append(item)
        if t_list3 != []:
            cats.append("Group 3")
            for item in t_list3:
                grouplist.append("Group 3")
                gclist.append(item)
        if t_list4 != []:
            cats.append("Group 4")
            for item in t_list4:
                grouplist.append("Group 4")
                gclist.append(item)
    df = pd.DataFrame({"group": grouplist, "gc": gclist})
    groups = df.groupby('group')
    q1 = groups.quantile(q=0.25)
    q2 = groups.quantile(q=0.5)
    q3 = groups.quantile(q=0.75)
    iqr = q3 - q1
    upper = q3 + 1.5 * iqr
    lower = q1 - 1.5 * iqr

    def outliers(group):
        """

        Find the outliers for each category
        Parameters:
        - group: dataframe

        Returns:

        Example:
        """
        cat = group.name
        return group[(group.gc > upper.loc[cat]['gc']) |
                     (group.gc < lower.loc[cat]['gc'])]['gc']

    out = groups.apply(outliers).dropna()

    # prepare outlier data for plotting, we need coordinates for every outlier.
    if not out.empty:
        outx = []
        outy = []
        for keys in out.index:
            try:
                outx.append(keys[0])
                outy.append(out.loc[keys[0]].loc[keys[1]])
            except Exception:
                pass
    full_title = "{}% ({})".format(nucleotide, short_code)
    y_lab = '{} %'.format(nucleotide)
    qmin = groups.quantile(q=0.00)
    qmax = groups.quantile(q=1.00)
    upper.gc = [
        min([x, y]) for (x, y) in zip(list(qmax.loc[:, 'gc']), upper.gc)
    ]
    lower.gc = [
        max([x, y]) for (x, y) in zip(list(qmin.loc[:, 'gc']), lower.gc)
    ]

    # stems
    p.segment(cats, upper.gc, cats, q3.gc, line_color="black")
    p.segment(cats, lower.gc, cats, q1.gc, line_color="black")

    # boxes
    p.vbar(cats, 0.7, q2.gc, q3.gc, fill_color=box_colour, line_color="black")
    p.vbar(cats, 0.7, q1.gc, q2.gc, fill_color=box_colour, line_color="black")

    # whiskers (almost-0 height rects simpler than segments)
    p.rect(cats, lower.gc, 0.2, 0.01, line_color="black")
    p.rect(cats, upper.gc, 0.2, 0.01, line_color="black")
    if not out.empty:
        p.circle(outx, outy, size=6, color="#F38630", fill_alpha=0.6)


def lengths_box(master_dict, filename, box_colour, short_code, title_size,
                marker_size, axis_label_size):
    """

    Parameters:
    - master_dict (Dict[str, Dict[str, str]]): master dictionary
    - filename (str): filename
    - box_colour (str): box colour
    - short_code (str): short code
    - title_size (int): title size
    - marker_size (int): marker size
    - axis_label_size (int): axis label size

    Returns:

    Example:
    """
    gc_list = master_dict[1]["lengths"]
    gc_list2 = master_dict[2]["lengths"]
    gc_list3 = master_dict[3]["lengths"]
    gc_list4 = master_dict[4]["lengths"]
    cats = []
    grouplist = []
    gclist = []

    cats.append("Lengths_G1")
    for item in gc_list:
        grouplist.append("Lengths_G1")
        gclist.append(item)
    if gc_list2 != []:
        cats.append("Lengths_G2")
        for item in gc_list2:
            grouplist.append("Lengths_G2")
            gclist.append(item)
    if gc_list3 != []:
        cats.append("Lengths_G3")
        for item in gc_list3:
            grouplist.append("Lengths_G3")
            gclist.append(item)
    if gc_list4 != []:
        cats.append("Lengths_G4")
        for item in gc_list4:
            grouplist.append("Lengths_G4")
            gclist.append(item)

    df = pd.DataFrame({"group": grouplist, "lengths": gclist})
    groups = df.groupby('group')
    q1 = groups.quantile(q=0.25)
    q2 = groups.quantile(q=0.5)
    q3 = groups.quantile(q=0.75)
    iqr = q3 - q1
    upper = q3 + 1.5 * iqr
    lower = q1 - 1.5 * iqr

    # find the outliers for each category
    def outliers(group):
        """
        Find the outliers for each category
        """
        # TODO: Merge with other same function
        cat = group.name
        return group[(group.lengths > upper.loc[cat]['lengths']) |
                     (group.lengths < lower.loc[cat]['lengths'])]['lengths']

    out = groups.apply(outliers).dropna()

    # prepare outlier data for plotting, we need coordinates for every outlier.
    if not out.empty:
        outx = []
        outy = []
        for keys in out.index:
            outx.append(keys[0])
            outy.append(out.loc[keys[0]].loc[keys[1]])

    full_title = "Lengths ({})".format(short_code)
    p = figure(plot_width=1300,
               plot_height=1300,
               tools="reset,pan,box_zoom,save,hover,tap",
               title=full_title,
               background_fill_color="#efefef",
               x_range=cats,
               toolbar_location="below")
    # if no outliers, shrink gcs of stems to be no longer than the minimums or maximums
    qmin = groups.quantile(q=0.00)
    qmax = groups.quantile(q=1.00)
    upper.lengths = [
        min([x, y])
        for (x, y) in zip(list(qmax.loc[:, 'lengths']), upper.lengths)
    ]
    lower.lengths = [
        max([x, y])
        for (x, y) in zip(list(qmin.loc[:, 'lengths']), lower.lengths)
    ]

    # stems
    p.segment(cats, upper.lengths, cats, q3.lengths, line_color="black")
    p.segment(cats, lower.lengths, cats, q1.lengths, line_color="black")

    # boxes
    p.vbar(cats,
           0.7,
           q2.lengths,
           q3.lengths,
           fill_color=box_colour,
           line_color="black")
    p.vbar(cats,
           0.7,
           q1.lengths,
           q2.lengths,
           fill_color=box_colour,
           line_color="black")

    # whiskers (almost-0 height rects simpler than segments)
    p.rect(cats, lower.lengths, 0.2, 0.01, line_color="black")
    p.rect(cats, upper.lengths, 0.2, 0.01, line_color="black")
    if not out.empty:
        p.circle(outx, outy, size=6, color="#F38630", fill_alpha=0.6)


def gene_count(short_code, background_col, title_size, axis_label_size,
               subheading_size, marker_size, coding, noncoding):
    """

    Parameters:
    - short_code (str): short code
    - background_col (str): background color
    - title_size (int): title size
    - axis_label_size (int): axis label size
    - subheading_size (int): subheading size
    - marker_size (int): marker size
    - coding (List[int]): coding
    - noncoding (List[int]): noncoding

    Returns:

    Example:
    """
    title_str = "Reads breakdown ({})".format(short_code)

    if len(labels) > 12:
        marker_size = int(marker_size / (len(labels) / 8))

    totals = []
    for i in range(0, N):
        curr_total = 0
        curr_total += coding[i]
        curr_total += noncoding[i]
        if curr_total > 0:
            all_reads_count += curr_total
            totals.append(float(curr_total))
        else:
            totals.append(1)

    for i in range(0, len(coding)):
        per = (coding[i] / totals[i]) * 100

    for i in range(0, len(noncoding)):
        per = (noncoding[i] / totals[i]) * 100

    p1 = plt.bar(ind, coding, bar_width, color='#80ff80', linewidth=0)
    p2 = plt.bar(ind,
                 noncoding,
                 bar_width,
                 color='#ff7275',
                 bottom=coding,
                 linewidth=0)
    p8 = plt.bar(ind, totals, bar_width, color='#5e0003', linewidth=0, alpha=0)

    #Dummy plot point so we can add total reads to the legend
    plt.plot(0, 0, alpha=0)
    plt.legend((p2[0], p1[0]), ('Non Coding: {:,}'.format(
        sum(noncoding)), 'Coding {:,}'.format(sum(coding))))
