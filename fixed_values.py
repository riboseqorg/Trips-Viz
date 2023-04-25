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


def merge_dicts(dict1, dict2):
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
