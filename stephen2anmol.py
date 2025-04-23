import pickle
from glob import glob
from os import makedirs, path

import pandas as pd
import polars as pl
from sqlitedict import SqliteDict

for fl in glob("trips_data_sample/*/escherichia_coli/Li14/*.sqlite"):
    print(fl)
    outfile = fl.replace("trips_data_sample", "anmol_trips_data_sample")
    outfold = path.dirname(outfile)
    makedirs(outfold, exist_ok=True)

    sqlite_db = SqliteDict(fl)
    converted = SqliteDict(outfile, autocommit=True)
    fields = {
        "totals",
        "unambiguous_threeprime_totals",
        'unambiguous_fiveprime_totals',
        "unambiguous_all_totals",
        "unambiguous_cds_totals",
        "reads_at_least_one_alignment",
        "threeprime_nuc_counts",
        "stop_metagene_counts",
        "ambiguous_counts",
        "cutadapt_removed",
        "noncoding_counts",
        "total_alignments",
        "trip_periodicity",
        "metagene_counts",
        "removed_minus_m",
        "unmapped_reads",
        "coding_counts",
        "dinuc_counts",
        "mapped_reads",
        "read_lengths",
        "rrna_removed",
        "total_reads",
        "nuc_counts",
        "offsets",
    }



    unambigx = []
    for k in ["fiveprime", "cds","threeprime", "all"]:
        unambigx.append( pl.DataFrame({'gene':sqlite_db[f"unambiguous_{k}_totals"].keys(),'count':sqlite_db[f"unambiguous_{k}_totals"].values(), "type":k}))
    unambigx2 = pl.concat(unambigx).pivot(values='count', index='gene', on='type').fill_null(0)
    unambigx = {}
    for row in unambigx2.rows(named=True):
        unambigx[row['gene']] =[row['fiveprime'],row['cds'],row['threeprime'],row['all']] 

    converted['unambiguous'] = {'data':unambigx,'schema':['fiveprime','cds','threeprime','all']}
    converted["gene_schemas"] = {"ambig":["readlen","pos","count"],"unambig":["readlen","pos","count"],"mismatches":["readlen","pos","count"],"seq":["pos","A","C","G","T"]}









    for gene in set(sqlite_db.keys()) - set(fields):
        dct = {}
        for k in sqlite_db[gene]:
            # print(k)
            dt = []
            if k in ["ambig", "unambig", "mismatches"]:
                for readlen, poses in sqlite_db[gene][k].items():
                    for pos, count in poses.items():
                        dt.append([readlen, pos, count])
                # print(dt)
                dct[k] = dt
            if k == "seq":
                dt = []
                nucs = []

                for pos, bases in sqlite_db[gene][k].items():
                    for base, count in bases.items():
                        if base not in nucs:
                            nucs.append(base)
                        dt.append([pos, base, count])
                for nuc in "ATGC": 
                    if nuc not in nucs:
                        dt.append([0    , nuc, 0])



                dct[k] = pl.DataFrame(dt, schema=["pos", "base", "count"], orient="row").pivot(
                    values="count", index="pos", on="base"
                ).fill_null(0).select(["pos","A", "C", "G", "T"]).to_pandas().values.tolist()

        converted[gene] = dct
        # print(dct)

        pass
    for field in fields:
        if field == "read_lengths":
            dt = [] 
            for k,v in sqlite_db[field].items():
                dt.append([k,v])


            converted[field] = {'data':dt, 'schema':['readlen','count']} 

        
        elif field == "threeprime_nuc_counts":
            all_dfs = []
            for k, v in sqlite_db[field].items():
                for k1, v1 in v.items():
                    all_dfs.append(
                        pl.DataFrame(v1).with_columns(
                            pl.lit(k1).alias("pos"), pl.lit(k).alias("readlen")
                        )
                    )
            converted[field] = pl.concat(all_dfs)
        elif field in [
            "reads_at_least_one_alignment",
            "ambiguous_counts",
            "cutadapt_removed",
            "noncoding_counts",
            "total_alignments",
            "removed_minus_m",
            "unmapped_reads",
            "mapped_reads",
            "rrna_removed",
            "total_reads",
                "totals",
        ]:
            converted[field] = sqlite_db[field]
        elif field in ["stop_metagene_counts", "metagene_counts"]:
            all_dfs_five = []
            for k, v in sqlite_db[field]["fiveprime"].items():
                all_dfs_five.append(
                    pl.DataFrame({"pos": v.keys(), "count": v.values()}).with_columns(
                        pl.lit(k).alias("readlen")
                    )
                )
            all_dfs_three = []
            for k, v in sqlite_db[field]["threeprime"].items():
                all_dfs_three.append(
                    pl.DataFrame({"pos": v.keys(), "count": v.values()}).with_columns(
                        pl.lit(k).alias("readlen")
                    )
                )
            converted[field] = {
                "fiveprime": pl.concat(all_dfs_five),
                "threeprime": pl.concat(all_dfs_three),
            }
        elif field == "trip_periodicity":
            all_dfs_five = []
            for k, v in sqlite_db[field]["fiveprime"].items():
                all_dfs_five.append(
                    pl.DataFrame(v).with_columns(pl.lit(k).alias("readlen"))
                )
            all_dfs_three = []
            for k, v in sqlite_db[field]["threeprime"].items():
                all_dfs_three.append(
                    pl.DataFrame(v).with_columns(pl.lit(k).alias("readlen"))
                )
            converted[field] = {
                "fiveprime": pl.concat(all_dfs_five),
                "threeprime": pl.concat(all_dfs_three),
            }
        elif field == ["dinuc_counts", "nuc_counts"]:
            all_dfs = []
            for k, v in sqlite_db[field].items():
                all_dfs.append(pl.DataFrame(v).with_columns(pl.lit(k).alias("readlen")))
            converted[field] = pl.concat(all_dfs)
        elif field == "offsets":
            converted[field] = {
                "fiveprime": pl.DataFrame(
                    {
                        "readlen": sqlite_db[field]["fiveprime"]["offsets"].keys(),
                        "offset": sqlite_db[field]["fiveprime"]["offsets"].values(),
                    }
                ).join(
                    pl.DataFrame(
                        {
                            "readlen": sqlite_db[field]["fiveprime"][
                                "read_scores"
                            ].keys(),
                            "read_scores": sqlite_db[field]["fiveprime"][
                                "read_scores"
                            ].values(),
                        }
                    ),
                    on="readlen",
                    how="inner",
                ),
                "threeprime": pl.DataFrame(
                    {
                        "readlen": sqlite_db[field]["threeprime"]["offsets"].keys(),
                        "offset": sqlite_db[field]["threeprime"]["offsets"].values(),
                    }
                ).join(
                    pl.DataFrame(
                        {
                            "readlen": sqlite_db[field]["threeprime"][
                                "read_scores"
                            ].keys(),
                            "read_scores": sqlite_db[field]["threeprime"][
                                "read_scores"
                            ].values(),
                        }
                    ),
                    on="readlen",
                    how="inner",
                ),
            }

    for k in sqlite_db:
        if k in converted:
            continue
        converted[k] = sqlite_db[k]
    converted.commit()
    converted.close()
    sqlite_db.close()
