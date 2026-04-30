from sqlitedict import SqliteDict
import polars as pl
import pickle
from glob import glob
from os import path, makedirs

for fl in glob("trips_data_sample/*/escherichia_coli/Li14/*.sqlite"):
    print(fl)
    outfile = fl.replace("trips_data_sample", "anmol_trips_data_sample")
    outfold = path.dirname(outfile)
    makedirs(outfold, exist_ok=True)

    sqlite_db = SqliteDict(fl)
    new_db = {"genes": {}}
    non_gene = [
        "ambiguous_counts",
        "coding_counts",
        "cutadapt_removed",
        "desc",
        "dinuc_counts",
        "mapped_reads",
        "metagene_counts",
        "noncoding_counts",
        "nuc_counts",
        "offsets",
        "read_lengths",
        "reads_at_least_one_alignment",
        "removed_minus_m",
        "rrna_removed",
        "stop_metagene_counts",
        "threeprime_nuc_counts",
        "total_alignments",
        "total_reads",
        "totals",
        "trip_periodicity",
        "unambiguous_all_totals",
        "unambiguous_cds_totals",
        "unambiguous_fiveprime_totals",
        "unambiguous_threeprime_totals",
        "unmapped_reads",
    ]

    for key in non_gene:

        if key in [
            "unambiguous_all_totals",
            "unambiguous_cds_totals",
            "unambiguous_fiveprime_totals",
            "unambiguous_threeprime_totals",
        ]:
            new_db[key] = pl.DataFrame(
                {"gene": sqlite_db[key].keys(), "count": sqlite_db[key].values()}
            )

        elif key == "totals":
            dfs = []
            for key2 in sqlite_db[key]:  # key2 = gene
                # print(sqlite_db[key][key2])
                tdf = pl.DataFrame(
                    [[key2] + sqlite_db[key][key2]],
                    schema=["gene", "fiveprime", "CDS", "threeprime"],
                )  # .with_columns(gene=key2)
                # print(tdf, key2)
                dfs.append(tdf)
            new_db[key] = pl.concat(dfs)

        elif key == "trip_periodicity":
            new_db[key] = {}
            for key2 in sqlite_db[key]:  # key2 = fiveprime, threeprime
                dfs = []
                for k in sqlite_db[key][key2]:  # k = readlen
                    tdf = pl.DataFrame(sqlite_db[key][key2][k]).with_columns(read_len=k)
                    dfs.append(tdf)
                new_db[key][key2] = pl.concat(dfs)

        elif key == "read_lengths":
            new_db[key] = pl.DataFrame(
                {"read_len": sqlite_db[key].keys(), "count": sqlite_db[key].values()}
            )
        elif key == "dinuc_counts":
            dfs = []  # Allow empty dict
            for k in sqlite_db[key]:
                tdf = pl.DataFrame(sqlite_db[key][k]).with_columns(read_len=k)
                dfs.append(tdf)
            new_db[key] = pl.concat(dfs)
        elif key in ["metagene_counts", "stop_metagene_counts"]:
            new_db[key] = {}
            for key2 in sqlite_db[key]:  # key2 = fiveprime, threeprime
                dfs = []
                for k in sqlite_db[key][key2]:  # k = readlen
                    tdf = sqlite_db[key][key2][k]
                    tdf = pl.DataFrame(
                        {"pos": tdf.keys(), "count": tdf.values()}
                    ).with_columns(read_len=k)
                    dfs.append(tdf)
                new_db[key][key2] = pl.concat(dfs)
        elif key == "nuc_counts":
            dfs = []  # Allow empty dict
            for read_len in sqlite_db[key]:
                for k in sqlite_db[key][read_len]:  # k = pos

                    tdf = pl.DataFrame(sqlite_db[key][read_len][k]).with_columns(
                        read_len=read_len, pos=k
                    )
                    dfs.append(tdf)
            new_db[key] = pl.concat(dfs)
        elif key == "threeprime_nuc_counts":
            dfs = []  # Allow empty dict
            for key2 in sqlite_db[key]:  # key2 = readlen
                for k in sqlite_db[key][key2]:  # k = pos
                    tdf = pl.DataFrame(sqlite_db[key][key2][k]).with_columns(
                        read_len=key2, pos=k
                    )
                    dfs.append(tdf)
            new_db[key] = pl.concat(dfs)
        elif key == "offsets":
            new_db[key] = {}
            for key2 in sqlite_db[key]:  # key2 = fiveprime, threeprime
                dfs = []
                for k in sqlite_db[key][key2]:  # k = readlen
                    tdf = sqlite_db[key][key2][k]
                    tdf = pl.DataFrame({"read_len": tdf.keys(), k: tdf.values()})
                    dfs.append(tdf)

                new_db[key][key2] = dfs[0].join(
                    dfs[1], on="read_len", how="outer_coalesce"
                )
                print(new_db[key][key2])

        else:
            new_db[key] = sqlite_db[key]

    for key in sqlite_db:
        if key not in non_gene:
            new_db["genes"][key] = {}
            gene_vals = sqlite_db[key]
            if "mismaches" in gene_vals:
                print(key, gene_vals["mismaches"])
            ## ambig
            dfs = []
            for k in gene_vals["unambig"]:
                tdict = gene_vals["unambig"][k]
                tdf = pl.DataFrame(
                    {"pos": list(tdict.keys()), "count": list(tdict.values())}
                ).with_columns(read_len=k)
                dfs.append(tdf)
            if dfs:
                dfs = pl.concat(dfs)
            else:
                dfs = pl.DataFrame(schema={"pos": int, "count": int, "read_len": int})

            new_db["genes"][key]["unambig"] = dfs

            dfs = []  # Allow empty dict
            for k in gene_vals["ambig"]:
                tdict = gene_vals["ambig"][k]
                tdf = pl.DataFrame(
                    {"pos": list(tdict.keys()), "count": list(tdict.values())}
                ).with_columns(read_len=k)
                dfs.append(tdf)
            if dfs:
                dfs = pl.concat(dfs)
            else:
                dfs = pl.DataFrame(schema={"pos": int, "count": int, "read_len": int})
            new_db["genes"][key]["ambig"] = dfs
            # seq

            dfs = []  # Allow empty dict
            for k in gene_vals["seq"]:
                tdict = gene_vals["seq"][k]
                tdf = pl.DataFrame(
                    {"nuc": list(tdict.keys()), "count": list(tdict.values())}
                ).with_columns(pos=k)
                dfs.append(tdf)
            if dfs:
                dfs = pl.concat(dfs)
            else:
                dfs = pl.DataFrame(schema={"nuc": int, "count": int, "pos": int})
            new_db["genes"][key]["seq"] = dfs

    sqlite_db.close()
    with open(outfile, "wb") as f:
        pickle.dump(new_db, f)

    # print(new_db)
