import sqlite3
from collections import defaultdict
from sys import argv

import numpy as np
import pandas as pd
import polars as pl
from sqlitedict import SqliteDict
from tqdm import tqdm


def table2dict(table: pd.DataFrame, keys: list[str]) -> dict:
    '''
    Convert a table to a dictionary of lists. 
    >>> data = {'key1': [1, 2, 3], 'key2': [4, 5, 6], 'key3': [7, 8, 9], 
    ... 'key4': [10, 11, 12], 'key5': [13, 14, 15]}
    >>> table = pd.DataFrame(data)
    >>> table2dict(table, ['key1', 'key2', 'key3'])
    >>> {1:{4:{7:[10,13]}}, 2:{5:{8:[11,14]}}, 3:{6:{9:[12,15]}}}

    '''
    # print(keys)
    if not keys:
        return table.values.tolist()[0][0]
    key = keys[0]

    result = {}
    for k, group in table.groupby(key):
        if k == "":
            return {}
        res= table2dict(group.drop(columns=[key]), keys[1:])
        # print(k,res)
        if res:
            result[k] = res 
    return result




def nuc_count(table):
    data = []
    seq_len = table["seq"].apply(len)
    max_read_len = seq_len.max()
    for row in pl.from_pandas(table).iter_rows(named=True):
        read_len = len(row['seq'])
        data.append(list(row['seq'])+[""]*(max_read_len - read_len)+[row['count']])
    data = pd.DataFrame(data)
    data.columns = data.columns.astype(str)
    data["length"] = seq_len
    nuc_count_df = []
    for k in range(max_read_len):
        data["pos"] = k
        nuc_count_df.append(data.groupby(["length","pos", str(k)])[str(max_read_len)].sum().reset_index().rename(columns={str(k):"nuc"}))
    nuc_count_df = pd.concat(nuc_count_df)
    return table2dict(nuc_count_df,["length","pos","nuc"])


def dinuc_count(row):

    read_len = len(row['seq'])
    di_nuc_dist = []
    for i in range(0, read_len - 1):
        dinuc = row['seq'][i:i + 2]
        di_nuc_dist.append(dinuc)
    di_nuc_dist = pd.DataFrame(di_nuc_dist, columns=["dinuc"]).groupby("dinuc").size().reset_index().rename(columns={0:"count"})
    di_nuc_dist["count"] *= row['count']
    di_nuc_dist["readlen"] = read_len
    return di_nuc_dist




if __name__ == "__main__":
    if len(argv) != 2:
        print("Usage: bam_to_sqlite.py <bam_file> <annotation_sqlite_file>")
        # exit(1)

    sql_dict = SqliteDict("anmol" + ".sqlite", autocommit=True)
    sql_dict["rrna_removed"] = 0
    sql_dict["description"] = 'NULL'
    sql_dict["desc"] = 'NULL'
    sql_dict["cutadapt_removed"] =0

    print("Reading sam file")
    samfile = pd.read_csv(argv[1],sep="\t", comment="@", header=None, names=range(14), low_memory=True).drop([1,4,5,6,
 7,8,10,11], axis='columns').rename(columns = {0:"qname",2:"rname",3:"pos",9:"seq",12:"MD",13:"NM"})
    samfile["transcript"] = samfile.rname.apply(lambda x: x.split("|")[0].split(".")[0])


    # NOTE: read length distribution
    print("Read length distribution")

    sql_dict["read_lengths"] = samfile.drop_duplicates(["qname","seq"])["seq"].apply(len).reset_index().groupby("seq").size().to_dict()

    samfile_unmapped = samfile.loc[samfile["rname"]=="*",["seq"]].groupby("seq").size().reset_index().rename(columns={0:'count'})
    # NOTE: unmapped reads counts
    print("Unmapped reads")

    sql_dict["unmapped_reads"] = samfile_unmapped.shape[0]
    # NOTE: most frequent unmapped reads
    print("Most frequent unmapped reads")
    sql_dict["frequent_unmapped_reads"] = list(samfile_unmapped.sort_values("count",ascending=False).head(2000).itertuples(index=False,name=None))


    samfile_mapped = samfile.loc[samfile["rname"]!="*"]
    del samfile
    # NOTE: mapped reads counts
    print("Mapped reads")

    sql_dict["mapped_reads"] = samfile_mapped.shape[0]
    # NOTE:nuc counts

    print("Nuc counts")
    # sql_dict["nuc_counts"] = {"mapped":nuc_count(samfile_mapped.drop_duplicates(["qname","seq"]).groupby("seq").size().reset_index().rename(columns={0:'count'})), "unmapped":nuc_count(samfile_unmapped)} 



    samfile_unambig = samfile_mapped.drop_duplicates(["qname","seq"], keep=False)
    samfile_ambig = samfile_mapped[~samfile_mapped["qname"].isin(samfile_unambig["qname"])] #.groupby(list(samfile_mapped.columns[1:])).size().reset_index().rename(columns={0:'count'})
    print(samfile_ambig.drop_duplicates("seq").shape)


    # sql_dict['ambiguous_counts'] = samfile_ambig.shape[0]
    sql_dict['ambiguous_counts'] = samfile_ambig.drop_duplicates("qname").shape[0]
    # sql_dict['ambiguous_counts'] = samfile_ambig.drop_duplicates("seq").shape[0]
    samfile_unambig = samfile_unambig.groupby(list(samfile_mapped.columns[1:])).size().reset_index().rename(columns={0:'count'})
    # sql_dict["unambig_read_lengths"] = samfile_unambig.groupby("readlen")["count"].sum().to_dict()["seq"].apply(len).reset_index().groupby("seq").size().to_dict()

    sql_dict["dinuc_count"] = samfile_unambig.groupby("seq").size().reset_index().rename(columns={0:'count'}).apply(dinuc_count, axis=1)
    sql_dict.commit() 
    sql_dict.close()









	##
    conn = sqlite3.connect(argv[2]) 
    transcriptome = pl.read_database(
	    query="SELECT * FROM transcripts", 
	    connection=conn,
	).select(["transcript","cds_start","cds_stop","length","strand","chrom","tran_type"])
	
    exons = pl.read_database(
	    query="SELECT * FROM exons", 
	    connection=conn,
	)
		
	
# def mismat_pos(row):
#
#     def str2tag(string):
#         tags = string.split(" ")
#         tag_dict = "" 
#         for tag in tags:
#             tag_split = tag.split(":")
#             if tag_split[0] == "MD":
#                 tag_dict = tag_split[-1]
#             elif tag_split[0] == "NM" and tag_split[-1] == "0":
#                 return              
#         return tag_dict
#     md_tag = str2tag(row["tag"])
#     mismatches = []
#     if md_tag :
#         nucs = ["A","T","G","C"]
#         total_so_far = 0
#         prev_char = ""
#         for char in md_tag:
#             if char in nucs:
#                 if prev_char != "":
#                     total_so_far += int(prev_char)
#                     prev_char = ""
#                 mismatches.append([total_so_far+len(mismatches) , (readseq[total_so_far+len(mismatches)])])
#             else:
#                 if char != "^" and char != "N":
#                     if prev_char == "":
#                         prev_char = char
#                     else:
#                         total_so_far += int(prev_char+char)
#                         prev_char = ""
#     return pl.DataFrame(mismatches,schema=["pos","nuc"]).with_columns(pl.col("pos") + row["pos"])
#
#
#
#
# def reads_in_gene(gene):
#
#     exons_gene = exons.filter(pl.col("gene") == gene).sort("exon_start").with_columns(
#         shift = (pl.col("exon_stop") - pl.col("exon_start"))
#     )
#     if strand == "+": 
#         exons_gene = exons_gene.with_columns(
#             pl.col("shift").cum_sum().shift(fill_value=0)
#         )
#     else:
#         exons_gene = exons_gene.sort("exon_stop", descending=True).with_columns(
#             pl.col("shift").cum_sum().shift(fill_value=0)
#         )
#
#     samfile_uniq_gene = samfile_uniq.filter(pl.col("rname") == gene)
#     unambig = samfile_uniq_gene.with_columns(readlen = pl.col("seq").len()).select("readlen","pos").group_by(["readlen","pos"]).count()
#
#     samfile_duplicated_gene = samfile_duplicated.filter(pl.col("rname") == gene)
#     ambig = samfile_duplicated_gene.with_columns(readlen = pl.col("seq").len()).select("readlen","pos").group_by(["readlen","pos"]).count()
#     # Mismaches
#     mismatches = pl.concat([
#         samfile_uniq_gene.apply(mismat_pos).collect(),
#         samfile_duplicated_gene.apply(mismat_pos).collect()]).group_by("pos","nuc").count()
#
#     for row in exons_gene.iterrows(name=True):
#         t_samfile_uniq_gene = samfile_uniq_gene.filter(pl.col("pos") <= row["shift"])
#         t_samfile_duplicated_gene = samfile_duplicated_gene.filter(pl.col("pos") <= row["shift"])
#
#         # TODO Think about negative strand
#         pass
#
#     t_transcrips = transcriptome.filter(pl.col("chrom") == gene,)
#     # Select if both start and end are not empty
#     if t_transcrips.shape[0] == 0:
#         pass
#     else:
#         pass
#
#
#
#
#
#
#
#     del samfile_duplicated_gene
#     del samfile_uniq_gene
#     pass
#

# def offset(row):
#     samfile_uniq_gene = samfile_uniq.with_ .join(transcripts, on="transcript").with_columns(p5_read_pos=pl.col("pos")-pl.col("cds_start"), p5_stop_pos = pl.col("pos")-pl.col("cds_stop"), p3_read_pos = pl.pos('pos')+pl.col("seq").len()- pl.col("cds_start"), p3_read_stop_pos = pl.('pos') + pl.col("seq").len() -pl.col("cds_stop")).with_columns(p5_frame = p.col("p5_read_pos")%3, p3_frame=pl.col("p3_read_pos")%3)
#
#
    
