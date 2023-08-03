from typing import Dict, List, Tuple
import pickle
import pandas as pd

iupac_dict = {
    "A": "A",
    "U": "U",
    "G": "G",
    "C": "C",
    "R": "AG",
    "Y": "CU",
    "S": "GC",
    "W": "AU",
    "K": "GU",
    "M": "AC",
    "B": "CGU",
    "D": "AGU",
    "H": "ACU",
    "V": "ACG",
    "N": "AUGC"
}
ambig_nucs = "RYKMSWBDHVN"

aa_short_codon_list = {
    "gly": ["GGT", "GGC", "GGA", "GGG"],
    "arg": ["AGA", "AGG", "CGT", "CGC", "CGA", "CGG"],
    "ser": ["AGT", "AGC", "TCT", "TCC", "TCA", "TCG"],
    "trp": ["TGG"],
    "cys": ["TGT", "TGC"],
    "glu": ["GAA", "GAG"],
    "asp": ["GAT", "GAC"],
    "lys": ["AAA", "AAG"],
    "asn": ["AAT", "AAC"],
    "gln": ["CAA", "CAG"],
    "his": ["CAT", "CAC"],
    "tyr": ["TAT", "TAC"],
    "ala": ["GCT", "GCC", "GCA", "GCG"],
    "thr": ["ACT", "ACC", "ACA", "ACG"],
    "pro": ["CCT", "CCC", "CCA", "CCG"],
    "val": ["GTT", "GTC", "GTA", "GTG"],
    "met": ["ATG"],
    "ile": ["ATA", "ATT", "ATC"],
    "leu": ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
    "phe": ["TTT", "TTC"]
}

codon_list = [
    "ATG", "TTT", "TTC", "CTT", "CTC", "CTA", "CTG", "TTA", "TTG", "AGT",
    "AGC", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TGT", "TGC", "TGG",
    "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "AGA", "AGG",
    "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ACT", "ACC", "ACA",
    "ACG", "AAT", "AAC", "AAA", "AAG", "GTT", "GTC", "GTA", "GTG", "GCT",
    "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA",
    "GGG", "TAG", "TAA", "TGA"
]

codon_aa_full = {
    "TTT": "Phenylalanine",
    "TTC": "Phenylalanine",
    "TTA": "Leucine",
    "TTG": "Leucine",
    "TCT": "Serine",
    "TCC": "Serine",
    "TCA": "Serine",
    "TCG": "Serine",
    "TAT": "Tyrosine",
    "TAC": "Tyrosine",
    "TAA": "*",
    "TAG": "*",
    "TGT": "Cysteine",
    "TGC": "Cysteine",
    "TGA": "*",
    "TGG": "Tryptophan",
    "CTT": "Leucine",
    "CTC": "Leucine",
    "CTA": "Leucine",
    "CTG": "Leucine",
    "CCT": "Proline",
    "CCC": "Proline",
    "CCA": "Proline",
    "CCG": "Proline",
    "CAT": "Histidine",
    "CAC": "Histidine",
    "CAA": "Glutamine",
    "CAG": "Glutamine",
    "CGT": "Arginine",
    "CGC": "Arginine",
    "CGA": "Arginine",
    "CGG": "Arginine",
    "ATT": "Isoleucine",
    "ATC": "Isoleucine",
    "ATA": "Isoleucine",
    "ATG": "Methionine",
    "ACT": "Threonine",
    "ACC": "Threonine",
    "ACA": "Threonine",
    "ACG": "Threonine",
    "AAT": "Asparagine",
    "AAC": "Asparagine",
    "AAA": "Lysine",
    "AAG": "Lysine",
    "AGT": "Serine",
    "AGC": "Serine",
    "AGA": "Arginine",
    "AGG": "Arginine",
    "GTT": "Valine",
    "GTC": "Valine",
    "GTA": "Valine",
    "GTG": "Valine",
    "GCT": "Alanine",
    "GCC": "Alanine",
    "GCA": "Alanine",
    "GCG": "Alanine",
    "GAT": "Aspartic Acid",
    "GAC": "Aspartic Acid",
    "GAA": "Glutamic Acid",
    "GAG": "Glutamic Acid",
    "GGT": "Glycine",
    "GGC": "Glycine",
    "GGA": "Glycine",
    "GGG": "Glycine"
}


def merge_dicts(dict1: Dict[str, Dict[int, int]],
                dict2: Dict[str, Dict[int, int]]) -> Dict[str, Dict[int, int]]:
    for nuc in dict2:
        if nuc not in dict1:
            dict1[nuc] = dict2[nuc]
        else:
            for pos in dict2[nuc]:
                if pos not in dict1[nuc]:
                    dict1[nuc][pos] = dict2[nuc][pos]
                else:
                    dict1[nuc][pos] += dict2[nuc][pos]
    return dict1


# ----


def codon_usage(codon_dict, short_code, title_size, axis_label_size,
                marker_size, filename):
    allxvals = []
    allyvals = []
    alllabels = []
    amino_acids = []
    aa_dict = codon_aa_full.copy()

    curr_count = 0
    for codon in codon_list:
        curr_count += 1
        allxvals.append(curr_count)
        if codon in codon_dict:
            allyvals.append(codon_dict[codon])
        else:
            allyvals.append(0)
        alllabels.append(codon)
        amino_acids.append(aa_dict[codon])
    full_title = "Codon usage ({})".format(short_code)
    x_lab = ''
    y_lab = 'Count'
    min_y = min(0, min(allyvals)) - .02
    max_y = max(allyvals) * 1.05
    p = figure(plot_width=1300,
               plot_height=1300,
               x_axis_label=x_lab,
               y_axis_label=y_lab,
               title=full_title,
               toolbar_location="below",
               tools="reset,pan,box_zoom,hover,tap,save",
               y_range=(min_y, max_y))
    p.title.align = "center"
    p.title.text_font_size = title_size
    p.xaxis.axis_label_text_font_size = axis_label_size
    p.xaxis.major_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = axis_label_size
    p.yaxis.major_label_text_font_size = marker_size
    # p.background_fill_color = background_color
    p.xgrid.grid_line_color = "white"
    p.ygrid.grid_line_color = "white"
    color_palette_list = []
    colormap = cm.get_cmap(
        "gist_rainbow")  # choose any matplotlib colormap here
    start_val = 0.75
    for _ in range(0, 21):
        start_val -= 0.0293
        rgb = colormap(start_val)[:3]
        hexval = matplotlib.colors.rgb2hex(rgb)
        color_palette_list.append(hexval)

    color_map = bmo.CategoricalColorMapper(factors=[
        "Phenylalanine", "Leucine", "Serine", "Tyrosine", "*", "Cysteine",
        "Tryptophan", "Proline", "Histidine", "Glutamine", "Arginine",
        "Isoleucine", "Methionine", "Threonine", "Asparagine", "Lysine",
        "Valine", "Alanine", "Aspartic Acid", "Glutamic Acid", "Glycine"
    ],
                                           palette=color_palette_list)
    p.quad(
        top=[
            max_y, max_y, max_y, max_y, max_y, max_y, max_y, max_y, max_y,
            max_y, max_y
        ],
        bottom=[
            min_y, min_y, min_y, min_y, min_y, min_y, min_y, min_y, min_y,
            min_y, min_y
        ],
        left=[0.5, 3.5, 15.5, 19.5, 24.5, 28.5, 37.5, 43.5, 49.5, 55.5, 61.5],
        right=[1.5, 9.5, 17.5, 20.5, 26.5, 34.5, 41.5, 45.5, 53.5, 57.5, 64.5],
        color="#e0e0e0")
    source = ColumnDataSource({
        'x': allxvals,
        'y': allyvals,
        'labels': alllabels,
        'amino_acids': amino_acids
    })
    p.scatter('x',
              'y',
              source=source,
              alpha=1,
              color={
                  'field': 'amino_acids',
                  'transform': color_map
              },
              size=16,
              line_color="black")
    p.xaxis.major_label_overrides = {
        1: "Met",
        2.5: "Phe",
        6.5: "Leu",
        12.5: "Ser",
        16.5: "Tyr",
        18.5: "Cys",
        20: "Trp",
        22.5: "Pro",
        25.5: "His",
        27.5: "Gln",
        31.5: "Arg",
        36: "Ile",
        39.5: "Thr",
        42.5: "Asn",
        44.5: "Lys",
        47.5: "Val",
        51.5: "Ala",
        54.5: "Asp",
        56.5: "Glu",
        59.5: "Gly",
        63: "Stop"
    }
    p.xaxis.ticker = p.xaxis.major_label_overrides.keys()

    hover = p.select(dict(type=HoverTool))
    hover.tooltips = [("Count", "@y"), ("Codon", "@labels"),
                      ("Amino acid", "@amino_acids")]

    output_file("scatter10k.html", title="Codon usage")
    hover = p.select(dict(type=HoverTool))
    hover.mode = 'mouse'
    graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as csv file</b></button></a> </div>".format(
        filename)
    graph += file_html(p, CDN)
    return graph


def my_decoder(obj):
    return pickle.load(open(obj, "rb"))


def get_user_defined_seqs(
        seq: str, seqhili: List[str]
) -> Tuple[Dict[int, List[int]], Dict[int, List[str]]]:
    signalhtml = {0: [], 1: [], 2: []}
    seq = seq.replace("T", "U")
    near_cog_starts = {0: [], 1: [], 2: []}
    for i in range(0, len(seq)):
        for subseq in seqhili:
            subseq = subseq.upper()
            subseq = subseq.replace("T", "U").replace(" ", "")
            partial_seq = list(seq[i:i + len(subseq)])
            if len(partial_seq) != len(subseq):
                continue
            x = 0
            for x in range(0, len(subseq)):
                char = subseq[x]
                if partial_seq[x] in iupac_dict[char]:
                    partial_seq[x] = char
            partial_seq = "".join(partial_seq)
            if partial_seq == subseq:
                near_cog_starts[(i) % 3].append(i + 1)
                datadict = {'sequence': [subseq]}
                df = pd.DataFrame(datadict, columns=(["sequence"]))
                label = df.iloc[[0], :].T
                label.columns = ["Position: {}".format(i)]
                signalhtml[(i) % 3].append(str(label.to_html()))
    return near_cog_starts, signalhtml
