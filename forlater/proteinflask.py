import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.font_manager import FontProperties
import shelve
import mpld4
from mpld4 import plugins
import collections
from fetch_proteomic_reads import get_reads

#import seaborn

#seaborn.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
#seaborn.set_style("whitegrid", {'axes.grid' : False})

#Defaults
#lite = "y"
#min_read = 25
#max_read = 35
#ribo = "a"
#rna = "a"
#subcodon = "d"
#tran = ""
#ambig = "u"
#user_ribo_files  = []
#user_rna_files  = []
#offset_dict = {}

# Define some CSS to control our custom labels
point_tooltip_css = """
table
{
  border-collapse: collapse;
}
th
{
  color: #000000;
  background-color: #d2d4d8;
}
td
{
  background-color: #ffffff;
}
table, th, td
{
  font-family:Arial, Helvetica, sans-serif;
  border: 0px solid black;
  text-align: left;
}
"""


def merge_dict(dict1, dict2):
    master_dict = dict1
    for key in dict2:
        if key in master_dict:
            master_dict[key] += dict2[key]
        else:
            master_dict[key] = dict2[key]
    return master_dict


def set_axis_color(axis, color, alpha=None):
    """Sets the spine color of all sides of an axis (top, right, bottom, left)."""
    for side in ('top', 'right', 'bottom', 'left'):
        spine = axis.spines[side]
        spine.set_color(color)
        if alpha is not None:
            spine.set_alpha(alpha)


def get_color_palette(scheme):
    """Return colors for a given scheme. Default colors are returned for an item
    if undefined in scheme.

    """
    color_schemes = {
        'default': {
            'frames': ['#FF4A45', '#64FC44', '#5687F9'],
            'background': '#ffffff',
            'color': '#616161',
            'ticks': '#757575',
            'start': '#ffffff',
            'stop': '#909090',
            'rna': '#BFBFBF',
            'axis': '#e0e0e0',
            'grey': '#bdbdbd'
        },
        'colorbrewer': {
            'frames': ['#fc8d62', '#66c2a5', '#8da0cb']
        },
        'rgb': {
            'frames': ['red', 'green', 'blue']
        },
        'greyorfs': {}
    }

    colors = {}
    for k, v in color_schemes['default'].items():
        try:
            vals = color_schemes[scheme][k]
        except KeyError:
            vals = v
        colors[k] = vals
    return colors


def get_near_cog_starts(seq):
    near_cog_starts = {0: [], 1: [], 2: []}
    for i in range(0, len(seq)):
        codon = seq[i:i + 3]
        if codon == "CTG":
            near_cog_starts[(i + 1) % 3].append(i)
    return near_cog_starts


def get_fs_signals(seq):
    plusone_signals = [
        "AAATAAA", "CCCT", "CCCTGA", "CTTAGG", "CTTAGGC", "CTTTAAC", "CTTTGAC",
        "GCGA", "TCCTGA", "TTTTGA", "CTTAGTT"
    ]
    minusone_signals = [
        "AAAAAAC", "AAAAAAG", "AAAAAAT", "AAAAG", "AAATTTA", "AAATTTG",
        "AAATTTT", "AAATTTC", "CCCAAAA", "CCCAAAC", "CCCAAAG", "CCCTTTA",
        "CGAAAG", "GGATTTA", "GGATTTT", "GGGAAAA", "GGGAAAC", "GGGAAAG",
        "GGGAAAT", "GGGCCCC", "GGGCCCT", "GGGGAAC", "GGGTTTA", "GGGTTTT",
        "GGTTTTC", "GTTAAAC", "TTTAAAA", "TTTAAAC", "TTTAAAT", "TTTTTTA",
        "TTTTTTC", "TTTTTTG"
    ]
    signal_dict = {
        "plus": {
            1: {
                0: "A"
            },
            2: {
                0: "A"
            },
            3: {
                0: "A"
            }
        },
        "minus": {
            1: {
                0: "A"
            },
            2: {
                0: "A"
            },
            3: {
                0: "A"
            }
        }
    }

    for i in range(0, len(seq)):
        frame = i % 3
        frame += 1

        #Change this to range(4,8) to pick up fs signals of length 4 but they seem so common they are useless
        for x in range(5, 8):
            if seq[i:i + x] in plusone_signals:
                signal_dict["plus"][frame][i] = seq[i:i + x]
            if seq[i:i + x] in minusone_signals:
                signal_dict["minus"][frame][i] = seq[i:i + x]
    return signal_dict


#               (tran, minread, maxread, proteo_user_files, lite,organism,nucseq)
def plot_profile(tran, min_read, max_read, user_files_dict, lite, organism,
                 nucseq, minscore):
    user_proteo_files = []
    for seq_type in user_files_dict:
        for filename in user_files_dict[seq_type]:
            user_proteo_files.append(filename)

    labels = ["Frame 1 profiles", "Frame 2 profiles", "Frame 3 profiles"]
    all_stops = ["TAG", "TAA", "TGA"]
    transhelve = shelve.open(
        "/home/DATA/GWIPS_viz/Transcriptome_annotations/{0}/{0}.shelf".format(
            organism))
    traninfo = dict(transhelve[tran])
    transhelve.close()
    gene = traninfo["gene"]
    tranlen = traninfo["length"]
    cds_start = traninfo["cds_start"]
    cds_stop = traninfo["cds_stop"]
    if cds_start == "NULL":
        cds_start = 0
    if cds_stop == "NULL":
        cds_stop = 0
    all_starts = traninfo["start_list"]
    all_stops = traninfo["stop_list"]
    exon_junctions = traninfo["exon_junctions"]
    seq = traninfo["seq"].upper()

    start_stops = {
        1: {
            "starts": [0],
            "stops": [0]
        },
        2: {
            "starts": [0],
            "stops": [0]
        },
        3: {
            "starts": [0],
            "stops": [0]
        }
    }
    for start in all_starts:
        rem = (start % 3) + 1
        start_stops[rem]["starts"].append(start - 1)
    for stop in all_stops:
        rem = (stop % 3) + 1
        start_stops[rem]["stops"].append(stop - 1)

    all_subcodon_reads = get_reads(min_read, max_read, tran, user_proteo_files,
                                   organism, tranlen, minscore)
    try:
        subcodonmax = max(all_subcodon_reads.values())
    except:
        subcodonmax = 0

    y_max = subcodonmax * 1.25

    fig = plt.figure(figsize=(23, 12))
    colors = get_color_palette(scheme='default')
    gs = gridspec.GridSpec(3, 1, height_ratios=[6, 1.3, 0.5], hspace=0.35)
    font_axis = {
        'family': 'sans-serif',
        'color': 'black',
        'weight': 'bold',
        'size': 20
    }

    # riboseq bar plots
    gs2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0])
    ax2 = plt.subplot(gs2[0])

    label = 'Read count'
    ax2.set_ylabel(label, fontdict=font_axis, labelpad=30)
    label = 'Position (nucleotides)'
    ax2.set_xlabel(label, fontdict=font_axis, labelpad=10)
    ax2.set_ylim(0, y_max)

    #plot a fake exon junction at postion 0, needed in cases there are no exon junctions
    allexons = ax2.plot((-1, -1), (0, 1),
                        alpha=0.01,
                        color='black',
                        linestyle='-.',
                        linewidth=2)
    for exon in exon_junctions:
        allexons += ax2.plot((exon, exon), (0, y_max),
                             alpha=0.01,
                             color='black',
                             linestyle='-.',
                             linewidth=3)

    #dictionary for each frame in which the keys are the posistions and the values are the counts
    frame_counts = {
        0: collections.OrderedDict(),
        1: collections.OrderedDict(),
        2: collections.OrderedDict()
    }
    for key in all_subcodon_reads:
        start = key + 1
        rem = start % 3
        if rem == 1:  # frame 1
            frame = 2
        elif rem == 2:  # frame 2
            frame = 0
        elif rem == 0:  # frame 3
            frame = 1
        frame_counts[frame][key] = all_subcodon_reads[key]
        if lite == "n":
            frame_counts[frame][key + 1] = 0
            frame_counts[frame][key + 2] = 0

    if lite == "n":
        frame0subpro = ax2.fill_between(frame_counts[0].keys(),
                                        frame_counts[0].values(),
                                        alpha=0.75,
                                        label=labels,
                                        zorder=2,
                                        color="#FF4A45",
                                        linewidth=0)
        frame1subpro = ax2.fill_between(frame_counts[1].keys(),
                                        frame_counts[1].values(),
                                        alpha=0.75,
                                        label=labels,
                                        zorder=2,
                                        color="#64FC44",
                                        linewidth=0)
        frame2subpro = ax2.fill_between(frame_counts[2].keys(),
                                        frame_counts[2].values(),
                                        alpha=0.75,
                                        label=labels,
                                        zorder=2,
                                        color="#5687F9",
                                        linewidth=0)
    else:
        frame0subpro = ax2.plot(frame_counts[0].keys(),
                                frame_counts[0].values(),
                                alpha=0.75,
                                label=labels,
                                zorder=2,
                                color="#FF4A45",
                                linewidth=2)
        frame1subpro = ax2.plot(frame_counts[1].keys(),
                                frame_counts[1].values(),
                                alpha=0.75,
                                label=labels,
                                zorder=2,
                                color="#64FC44",
                                linewidth=2)
        frame2subpro = ax2.plot(frame_counts[2].keys(),
                                frame_counts[2].values(),
                                alpha=0.75,
                                label=labels,
                                zorder=2,
                                color="#5687F9",
                                linewidth=2)

    # draw cds start
    plt.plot((cds_start - 2, cds_start - 2), (0, y_max),
             'gray',
             linestyle='solid',
             linewidth=2)

    # draw cds end
    plt.plot((cds_stop + 3, cds_stop + 3), (0, y_max),
             'gray',
             linestyle='solid',
             linewidth=2)

    gs3 = gridspec.GridSpecFromSubplotSpec(5,
                                           1,
                                           subplot_spec=gs[1],
                                           hspace=0.1)
    axisbg = colors['frames']
    ax4 = plt.subplot(gs3[1], sharex=ax2)
    ax4.patch.set_facecolor(colors['frames'][0])
    ax5 = plt.subplot(gs3[2], sharex=ax2)
    ax5.patch.set_facecolor(colors['frames'][1])
    ax6 = plt.subplot(gs3[3], sharex=ax2)
    ax6.patch.set_facecolor(colors['frames'][2])
    ax7 = plt.subplot(gs3[4], sharex=ax2)
    ax7.set_xlabel('Transcript: {} Length: {} nt'.format(tran, tranlen),
                   fontdict=font_axis,
                   labelpad=10)
    ax7.set_ylabel('',
                   fontdict={
                       'family': 'sans-serif',
                       'color': 'black',
                       'weight': 'normal',
                       'size': '11'
                   },
                   rotation='horizontal',
                   labelpad=10,
                   verticalalignment='center')
    #ax8 = ax2.twinx()
    xy = 0
    if nucseq == True:
        ax7.set_axis_bgcolor("#F2F2F7")
        for char in seq:
            ax7.text(xy - 0.1, 0.2, seq[xy], fontsize=20, color="grey")
            xy += 1

    # Legend
    gs4 = gridspec.GridSpecFromSubplotSpec(1,
                                           1,
                                           subplot_spec=gs[2],
                                           hspace=0.1)
    for axisname in (ax4, ax5, ax6):
        axisname.tick_params(top=False,
                             left=False,
                             right=False,
                             labeltop=False,
                             labelleft=False,
                             labelright=False,
                             direction='out',
                             colors=colors['ticks'])
    for axis in (ax4, ax5, ax6, ax7):
        axis.tick_params(top=False,
                         left=False,
                         right=False,
                         bottom=False,
                         labeltop=False,
                         labelleft=False,
                         labelright=False,
                         labelbottom=False)
    axes = [ax2]

    fp = FontProperties(size='5')
    for axis in axes:
        set_axis_color(axis, colors['axis'])
        axis.tick_params(colors=colors['ticks'])
        for item in (axis.get_xticklabels() + axis.get_yticklabels()):
            item.set_fontproperties(fp)
            item.set_color(colors['color'])

    for axis, frame in ((ax4, 1), (ax5, 2), (ax6, 3)):
        color = colors['frames'][frame - 1]
        set_axis_color(axis, color, alpha=1)
        for item in (axis.get_xticklabels()):
            item.set_fontproperties(fp)
            item.set_color(colors['color'])
        axis.set_ylim(0, 0.2)
        axis.set_xlim(0, tranlen)
        starts = [(item, 1) for item in start_stops[frame]['starts']]
        stops = [(item, 1) for item in start_stops[frame]['stops']]
        start_colors = [colors['start'] for item in starts]
        axis.broken_barh(starts, (0.5, 1),
                         facecolors=start_colors,
                         edgecolors=start_colors,
                         label='start',
                         zorder=5)
        stop_colors = [colors['stop'] for item in stops]
        axis.broken_barh(stops, (0, 1),
                         facecolors=stop_colors,
                         edgecolors=stop_colors,
                         label='stop',
                         zorder=5)

        axis.set_ylabel('{}'.format(frame),
                        fontdict={
                            'family': 'sans-serif',
                            'color': colors['color'],
                            'weight': 'normal',
                            'size': '8'
                        },
                        rotation='horizontal',
                        labelpad=10,
                        verticalalignment='center')

        axis.tick_params(top=False,
                         left=False,
                         right=False,
                         labeltop=False,
                         labelleft=False,
                         labelright=False,
                         direction='out',
                         colors=colors['ticks'])

        ax4.set_ylim(0, 1)
        ax5.set_ylim(0, 1)
        ax6.set_ylim(0, 1)
        axis.set_ylabel('{}'.format(frame),
                        fontdict={
                            'family': 'sans-serif',
                            'color': colors['color'],
                            'weight': 'normal',
                            'size': '6'
                        },
                        rotation='horizontal',
                        labelpad=10,
                        verticalalignment='center')

        axis.tick_params(top=False,
                         left=False,
                         right=False,
                         labeltop=False,
                         labelleft=False,
                         labelright=False,
                         direction='out',
                         colors=colors['ticks'])

    #plt.title('{}'.format(gene), fontdict={'family': 'sans-serif', 'color': 'black','weight': 'bold', 'size': 87, 'y': 33})
    plt.title('{}'.format(gene), fontsize=32, y=34)
    line_collections = [frame0subpro, frame1subpro, frame2subpro]

    ilp = plugins.InteractiveLegendPlugin(line_collections,
                                          labels,
                                          alpha_unsel=0,
                                          alpha_sel=0.85)

    ax6.axes.get_yaxis().set_ticks([])
    ax5.axes.get_yaxis().set_ticks([])
    ax4.axes.get_yaxis().set_ticks([])

    plugins.connect(fig, ilp, plugins.TopToolbar())
    ax2.set_axis_bgcolor("#F2F2F7")
    ax2.tick_params('both', labelsize=21)
    ax2.grid(True, color="white", linewidth=20, linestyle="solid")
    return mpld4.fig_to_html(fig)
