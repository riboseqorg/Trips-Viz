# Adding genomic coordinates for start sites
from sqlalchemy import create_engine
import polars as pl
import pandas as pd
from sqlqueries_2 import sqlquery


import sqlite3
from tqdm import tqdm


def tran_to_genome(
    strand: str,
    exons: list,
    pos: int,
):

    if strand == "+":
        exon_start = 0
        exonlen = 0
        for tup in exons:
            exon_start = tup[0]
            exonlen = tup[1] - tup[0]
            if pos > exonlen:
                pos = (pos - exonlen) - 1
            else:
                break
        genomic_pos = (exon_start + pos) - 1
    else:
        exons = exons[::-1]
        exon_start = 0
        for tup in exons:
            exon_start = tup[1]
            exonlen = tup[1] - tup[0]
            if pos > exonlen:
                pos = (pos - exonlen) - 1
            else:
                break
        genomic_pos = (exon_start - pos) + 1
    return genomic_pos


sqlfile = "escherichia_coli.Ensembl_k_12_ASM584v2.sqlite"
transcripts = sqlquery(sqlfile, "transcripts").select("transcript", "strand")
exons = sqlquery(sqlfile, "exons").sort("exon_start")

for tab in [
    "noncoding",
    "uorf",
    "cds",
    "ouorf",
    "extension",
    "nested",
    "odorf",
    "dorf",
]:

    # Transcripts

    # exons
    tabqs = sqlquery(sqlfile, tab).sort("transcript")  # .select("transcript", "start")
    new_list = []
    tx = None
    strand = None
    ex = None
    for tabq in tqdm(tabqs.iter_rows(named=True)):
        if tabq["transcript"] != tx:
            strand = transcripts.filter(pl.col("transcript") == tabq["transcript"])[
                0, "strand"
            ]
            ex = (
                exons.filter(pl.col("transcript") == tabq["transcript"])
                .select("exon_start", "exon_stop")
                .rows()
            )

            tx = tabq["transcript"]

        new_col_value = tran_to_genome(strand, ex, tabq["stop"])
        new_list.append(
            [
                tabq["transcript"],
                tabq["start_codon"],
                tabq["length"],
                tabq["start"],
                tabq["stop"],
                tabq["cds_coverage"],
                new_col_value,
            ]
        )

    df = pd.DataFrame(
        new_list,
        columns=[
            "transcript",
            "start_codon",
            "length",
            "start",
            "stop",
            "cds_coverage",
            "genomic_start",
        ],
    )
    engine = create_engine(
        "sqlite:///escherichia_coli.Ensembl_k_12_ASM584v2.sqlite", echo=False
    )
    df.to_sql(name=tab, con=engine, if_exists="replace", index=False)
