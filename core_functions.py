import logging
import os
import string
import uuid
from json import dumps
from typing import Any, Dict, List, Tuple

import pandas as pd
import polars as pl
from Bio.Seq import Seq
from flask import flash, request, session
from flask_login import UserMixin, current_user

import config
from sqlqueries_2 import get_table, sqlquery, table2dict, update_table


# User model
class User(UserMixin):

    def __init__(self, id):
        self.id = id
        self.name = str(id)
        self.password = self.name + "_secret"

    def is_authenticated(self) -> bool:
        """
        Checks if user is authenticated.

        Parameters:
        - None

        Returns:
        - bool
        """
        return self.is_authenticated

    def get_id(self) -> str:
        """
        Returns user id.

        Parameters:
        - None

        Returns:
        - str
        """
        return str(self.id)

    def __repr__(self) -> str:
        return "%d/%s/%s" % (self.id, self.name, self.password)


def dict2df(sqldict: Dict, keys: List) -> pl.DataFrame | List:
    dfs = []  # Allow empty dict
    print(sorted(list(sqldict.keys()), reverse=True)[:30], keys)
    if len(keys) == 2:  # For gene and triplet periodicity use two keys
        try:
            tdict = sqldict[keys[0]][keys[1]]
        except KeyError:
            return pl.DataFrame()
        if keys[1] in ["fiveprime", "threeprime"]:  # NOTE: keys[0]=="trip_periodity"
            if keys[0] == "trip_periodicity":
                for read_len, periodicity in tdict.items():
                    dfs.append(
                        pl.DataFrame(periodicity).with_columns(read_len=read_len)
                    )

                if dfs:
                    dfs = pl.concat(dfs)
                else:
                    dfs = pl.DataFrame()
            elif keys[0] == "offsets":
                print("mmmmmmmmmmmmmmmmmmmmm", tdict)
                return tdict
                df = (
                    pl.DataFrame(
                        {
                            "read_lens": tdict["read_scores"].keys(),
                            "read_scores": tdict["read_scores"].values(),
                        }
                    )
                    .join(
                        pl.DataFrame(
                            {
                                "read_lens": tdict["offsets"].keys(),
                                "offsets": tdict["offsets"].values(),
                            }
                        ),
                        on="read_lens",
                        how="outer_coalesce",
                    )
                    .with_columns(
                        [
                            pl.col("offsets").fill_null(15),
                            pl.col("read_scores").fill_null(1),
                        ]
                    )
                )
                print("zzzzzzzzzzzzzzzz", df)
                return df

        elif keys[1] == "seq":
            for pos, nuc_dist in tdict:
                dfs.append(pl.DataFrame(nuc_dist).with_columns(pos=pos))

            if dfs:
                dfs = pl.concat(dfs)
            else:
                dfs = pl.DataFrame()
        elif keys[1] in ["ambig", "unambig"]:
            for read_len, pos_count in tdict.items():
                dfs.append(
                    pl.DataFrame(
                        {"pos": pos_count.keys(), "count": pos_count.values()}
                    ).with_columns(read_len=read_len)
                )

            if dfs:
                dfs = pl.concat(dfs)
            else:
                dfs = pl.DataFrame()

        # TODO: Need to check the configuration of "mismatches"

        else:
            dfs = pl.DataFrame()
    elif keys[0] in [
        "unambiguous_all_totals",
        "unambiguous_cds_totals",
        "unambiguous_fiveprime_totals",
        "unambiguous_threeprime_totals",
    ]:
        pl.DataFrame(
            {"gene": sqldict[keys[0]].keys(), "count": sqldict[keys[0]].values()}
        )
    elif keys[0] == "totals":
        for key2 in sqldict[keys[0]]:  # key2 = gene
            # print(sqlite_db[key][key2])
            tdf = pl.DataFrame(
                [[key2] + sqldict[keys[0]][key2]],
                schema=["gene", "fiveprime", "CDS", "threeprime"],
            )  # .with_columns(gene=key2)
            # print(tdf, key2)
            dfs.append(tdf)
        dfs = pl.concat(dfs)
    elif keys[0] == "read_lengths":
        dfs = pl.DataFrame(
            {"read_len": sqldict[keys[0]].keys(), "count": sqldict[keys[0]].values()}
        )
    elif keys[0] == "dinuc_counts":
        for read_len, dinuc_count in sqldict[keys[0]].items():
            dfs.append(pl.DataFrame(dinuc_count).with_columns(read_len=read_len))
        if dfs:
            dfs = pl.concat(dfs)
        else:
            dfs = pl.DataFrame()
    elif keys[0] in ["nuc_counts", "threeprime_nuc_counts"]:
        for read_len, nuc_counts in sqldict[keys[0]].items():
            for pos, nuc_count in nuc_counts.items():

                dfs.append(
                    pl.DataFrame(nuc_count).with_columns(pos=pos, read_len=read_len)
                )
        if dfs:
            dfs = pl.concat(dfs)
        else:
            dfs = pl.DataFrame()

    else:

        dfs = pl.DataFrame()

    return dfs


def form_filler(organism, transcriptome):
    data: Dict[str, Any] = request.args.to_dict()
    data["organism"] = organism
    data["transcriptome"] = transcriptome
    gwips_info = get_table("organisms").filter(
        (pl.col("organism_name") == organism)
        & (pl.col("transcriptome_list") == transcriptome)
    )[
        0,
        [
            "gwips_clade",
            "gwips_organism",
            "gwips_database",
            "default_transcript",
            "organism_id",
        ],
    ]
    # print(accepted_studies)

    data["transcript"] = gwips_info[0, "default_transcript"]
    data["gwips_info"] = gwips_info
    data["user_hili_starts"] = []
    data["user_hili_stops"] = []
    try:
        for item in data["user_hili"].split(","):
            item_split = item.split("_")
            data["user_hili_starts"].append(int(item_split[0]))
            data["user_hili_stops"].append(int(item_split[1]))
    except Exception:
        pass
    data["user_hili_starts"] = ",".join(data["user_hili_starts"])
    data["user_hili_stops"] = ",".join(data["user_hili_stops"])

    return data


def fetch_user() -> Tuple[str | None, bool]:
    """
    Fetches active user from cookies if present and returns username and login status.

    Parameters:
    - None

    Returns:
    - Tuple
    """
    consent = request.cookies.get("cookieconsent_status")
    # If user rejects cookies then do not track them and delete all other cookies
    if consent == "deny":
        return (None, False)
    session.permanent = True
    print(session.keys())
    if "uid" not in session:
        session["uid"] = uuid.uuid4()
    session_id = str(session["uid"])
    # Check if this session uid is already in the users table
    users = get_table("users")
    if session_id not in users["username"]:
        # Add session uid to user table
        update_table(
            "users",
            {
                "user_id": None,
                "username": session_id,
                "password": None,
                "study_access": "-1",
                "organism_access": "",
                "advanced": 0,
                "temp_user": 1,
            },
        )
        user_id = max(users["user_id"]) + 1
        defaul_user_settings = config.DEFAULT_USER_SETTINGS.copy()
        for stop in ["uaa", "uag", "uga"]:
            defaul_user_settings[f"comp_{stop}_col"] = defaul_user_settings[
                f"{stop}_col"
            ]
        defaul_user_settings["user_id"] = user_id
        update_table("user_settings", defaul_user_settings, "insert")
    # id and login starus
    try:
        return current_user.name, True
    except Exception:
        return session_id, False


# Given a username and an organism returns a list of relevant studies.
def fetch_studies(organism_id) -> pl.DataFrame:
    """
    Fetches studies from database using organism and transcriptome information.

    Parameters:
    - organism (str): name of the organism
    - transcriptome (str): name of the transcript

    Returns:
    - Tuple[organism_id:str, studies:DataFrame[study_id, study_name]]
    """

    print(current_user.is_authenticated, type(current_user), "Anmol")
    # get a list of organism id's this user can access
    study_access_list = (
        get_table("study_access")
        .filter(pl.col("user_id") == current_user.id)["study_id"]
        .to_list()
        if current_user.is_authenticated
        else []
    )
    # TODO: Check if it can work without list
    # users is name of table

    # Getting organism id
    # Getting studies
    studies = (
        get_table("studies")
        .filter(
            (pl.col("organism_id") == organism_id)
            & (pl.col("study_id").is_in(study_access_list) | (pl.col("private") == 0))
        )
        .select("study_id", "study_name")
        .unique(subset=["study_id"])
    )  # users is name of table
    # TODO: Compare with original code and discuss which one too choose
    return studies


# Create a dictionary of files seperated by type, this allows for file type grouping on the front end.
def fetch_files(accepted_studies: pl.DataFrame) -> pl.DataFrame:
    """
    Fetches files from database for give studies.

    Parameters:
    - accepted_studies (DataFrame[study_id, study_name]): list of accepted studies

    Returns: -- Fix this part
    - {'seqtype':{('project_id', 'project_name')
                   : {('file_id', 'file_name'): ['file_description']}}}
    """
    return (
        get_table("files")
        .filter(pl.col("study_id").is_in(accepted_studies["study_id"]))
        .select("file_id", "study_id", "file_name", "file_description", "file_type")
        .with_columns(
            pl.col("file_name").map_elements(lambda x: x.replace(".shelf", ".sqlite"))
        )
        .join(accepted_studies, on="study_id")
    )


def string2other(dct: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert string in dict collected from web page form to right types.

    Parameters:
    - dct (Dict[str, Any]): dictionary of values

    Returns:
    - None
    """
    groups = []
    file_ids = []
    print(dct, "Kiran")
    for key, value in dct.items():
        if key in config.VARIABLE_CONVERSION:
            dct[key] = config.VARIABLE_CONVERSION[key](value)
        if key.endswith("_file_list"):

            file_id = list(map(int, value))
            file_ids.extend(file_id)
            dct[key] = file_id
            groups.append(key)
    dct["groups"] = groups
    dct["file_ids"] = file_ids
    return dct


def fetch_study_info(organism_id: int) -> pl.DataFrame:
    """
    Fetches studies from database for organism.


    Parameters:
    - organism_id (int): id of organism

    Returns:
    - Dict[str, List[str]]

    Example:

    """
    studies = (
        get_table("studies")
        .filter(pl.col("organism_id") == organism_id)
        .select(
            "study_id",
            "paper_authors",
            "srp_nos",
            "paper_year",
            "paper_pmid",
            "paper_link",
            "gse_nos",
            "adapters",
            "paper_title",
            "description",
            "study_name",
        )
    )

    # "paper_link": row[5].strip('"'),  # generate link using pubmed id
    return studies


# Given a list of file id's as strings returns a list of filepaths to the sqlite files.
def fetch_file_paths(data: Dict[str, Any]) -> pl.DataFrame:
    """

    Parameters:
    - data (Dict[str, Any]): dictionary of values

    Returns:
    - DataFrame

    Example:

    """

    studies = get_table("studies").select(
        "study_id", "study_name"
    )  # users is name of table
    files = (
        get_table("files")
        .filter(pl.col("file_id").is_in(data["file_ids"]))
        .join(studies, on="study_id")
        .with_columns(
            pl.col("file_name")
            .map_elements(lambda x: x.replace(".shelf", ".sqlite"))
        )
    )

    files = files.with_columns(
        pl.struct("*")
        .map_elements(
            lambda x: (
                "{}/{}/{}/{}/{}/{}".format(
                    config.SCRIPT_LOC,
                    config.SQLITES_DIR,
                    x["file_type"],
                    data["organism"],
                    x["study_name"],
                    x["file_name"],
                )
                if x["owner"]
                else "{}/{}/{}".format(
                    config.UPLOADS_DIR, data["study"], x["file_name"]
                )
            )
        )
        .alias("path")
    )  # TODO: Fix this
    file_not_found = []
    for fl in files["path"]:
        if not os.path.isfile(fl):
            file_not_found.append(fl.split("/")[-1])

    if file_not_found:
        # TODO: Fix the resturn accorind to original as it mught be used by js
        flash(f"File(s) not found: {','.join(file_not_found)}")

    # logging.debug("fetch_file_paths closing connection")
    return files


# Builds a url and inserts it into sqlite database
def generate_short_code(data) -> str:
    """
    Generates a short code for a plot

    Parameters:
    - data (Dict[str, Any]): dictionary of values
    - organism (str): name of the organism
    - transcriptome (str): name of the transcript
    - plot_type (str): type of plot

    Returns:
    - str : short code

    Example:
    """
    # TODO: Keep the form values to retore it in json format share it here before you plot
    # build a url so that this plot can be recreated later on
    key2remove = []
    for key, value in data.items():
        if type(value) in [pl.DataFrame, pd.DataFrame]:
            key2remove.append(key)
    for key in key2remove:
        del data[key]
    url_id: int = get_table("urls")["url_id"].max() + 1

    # If the url table is empty result will return none
    update_table("urls", {"url_id": url_id, "url": dumps(data)})
    short_code = integer_to_base62(url_id)
    return short_code


# Converts an integer to base62, needed to encode short urls
def integer_to_base62(num: int) -> str:
    """
    Converts an integer to base62, needed to encode short urls.

    Parameters:
    - num (int): integer to be converted

    Returns:
    - str

    Example:
    """
    base = string.digits + string.ascii_lowercase + string.ascii_uppercase
    r = num % 62
    res = base[r]
    q = num // 62
    while q:
        r = q % 62
        q //= 62
        res = base[int(r)] + res
    return res


# Converts a base62 encoded string to an integer, needed to decode short urls
def base62_to_integer(base62_str: str) -> int:
    """
    Converts a base62 encoded string to an integer, needed to decode short urls.

    Parameters:
    - base62_str (str): base62 encoded string

    Returns:
    - int

    Example:
    """
    base = string.digits + string.ascii_lowercase + string.ascii_uppercase
    res = 0
    for i in base62_str:
        res = 62 * res + base.find(i)
    return res


# Takes a nucleotide string and returns the amino acide sequence
def nuc_to_aa(nuc_seq: str) -> str:
    """
    Takes a nucleotide string and returns the amino acide sequence.

    Parameters:
    - nuc_seq (str): nucleotide sequence

    Returns:
    - str: amino acid sequence

    Example:


    """
    return str(Seq(nuc_seq).translate())


# Calculates the coverage of each gene, for 5' leader, cds and 3' trailer for unambiguous and ambigous reads, needed for diff exp
def calculate_coverages(
    sqlite_db: Dict[str, Dict[str, Dict[int, int]]],
    longest_tran_list: List[str],
    traninfo_dict: Dict[str, Dict[str, int]],
) -> None:
    """

    Parameters:
    - sqlite_db (Dict[str, Dict[str, Dict[int, int]]]): sqlite database
    - longest_tran_list (List[str]): list of longest transcripts
    - traninfo_dict (Dict[str, Dict[str, int]]): dictionary of transcript information

    Returns:
    - None

    Example:
    """
    coverage_types = [
        "unambig_fiveprime_coverage",
        "unambig_cds_coverage",
        "unambig_threeprime_coverage",
        "unambig_all_coverage",
        "ambig_fiveprime_coverage",
        "ambig_cds_coverage",
        "ambig_threeprime_coverage",
        "ambig_all_coverage",
    ]
    coverage_dict = {}
    for coverage_type in coverage_types:
        coverage_dict[coverage_type] = {}
        for tran in longest_tran_list:
            coverage_dict[coverage_type][tran] = 0
    for tran in set(sqlite_db) and set(traninfo_dict):
        unambig_dict = {}
        ambig_dict = {}
        ambig_range = [1000000, 0]
        unambig_range = [1000000, 0]
        cds_start = traninfo_dict[tran]["cds_start"]
        cds_stop = traninfo_dict[tran]["cds_stop"]
        tranlen = float(traninfo_dict[tran]["length"])
        for readlen in sqlite_db[tran]["unambig"]:
            for pos in sqlite_db[tran]["unambig"][readlen]:
                if pos < unambig_range[0]:
                    unambig_range[0] = pos
                if pos + readlen > unambig_range[1]:
                    unambig_range[1] = pos + readlen
                # for i in range(pos, pos + readlen):
                # unambig_dict[i] = ""  # TODO: Use list instead of dict
        for readlen in sqlite_db[tran]["ambig"]:
            for pos in sqlite_db[tran]["ambig"][readlen]:
                if pos < ambig_range[0]:
                    ambig_range[0] = pos
                if pos + readlen > ambig_range[1]:
                    ambig_range[1] = pos + readlen
                # for i in range(pos, pos + readlen):
                # ambig_dict[i] = ""  # TODO: Use list instead of dict
        coverage_dict["unambig_all_coverage"][tran] = (
            unambig_range[1] - unambig_range[0] + 1
        ) / tranlen
        coverage_dict["ambig_all_coverage"][tran] = (
            ambig_range[1] - ambig_range[0] + 1
        ) / tranlen
        if cds_start != "None":
            cds_start = float(cds_start)  # Why Float
            cds_stop = float(cds_stop)
            cds_len = cds_stop - cds_start
            three_len = tranlen - cds_stop
            # Use list comprehension to count number of entries less than cds_start
            if cds_start > 0:  # TODO: Fix below as above
                coverage_dict["unambig_fiveprime_coverage"][tran] = (
                    sum(i < cds_start for i in unambig_dict.keys()) / cds_start
                )
                coverage_dict["ambig_fiveprime_coverage"][tran] = (
                    sum(i < cds_start for i in ambig_dict.keys()) / cds_start
                )
            coverage_dict["unambig_cds_coverage"][tran] = (
                sum(i > cds_start and i < cds_stop for i in unambig_dict.keys())
                / cds_len
            )
            coverage_dict["ambig_cds_coverage"][tran] = (
                sum(i > cds_start and i < cds_stop for i in ambig_dict.keys()) / cds_len
            )
            if three_len > 0:
                coverage_dict["unambig_threeprime_coverage"][tran] = (
                    sum(i > cds_stop for i in unambig_dict.keys()) / three_len
                )
                coverage_dict["ambig_threeprime_coverage"][tran] = (
                    sum(i > cds_stop for i in ambig_dict.keys()) / three_len
                )
        for coverage in coverage_types:
            sqlite_db[coverage] = coverage_dict[coverage]
        sqlite_db.commit()


def tran_count(trancounts, typ="unambig"):
    t_trancounts = trancounts[typ] if typ in trancounts else {}
    if t_trancounts:
        lists = []
        for readlen, counts in trancounts.items():
            for pos, count in counts.items():
                lists.append([readlen, pos, count])
        t_trancounts = pl.DataFrame(lists, schema=["readlen", "pos", "count"])
    else:
        t_trancounts = pl.DataFrame({"readlen": [], "pos": [], "count": []})

    return t_trancounts


# Builds a profile, applying offsets
def build_profile(
    trancounts: Dict[str, Dict[int, List[int]]],
    offsets_5p_offsetsNscores: Dict[int, int],
    ambig: bool,
):
    """

    Parameters:
    - trancounts
    - offsets
    - ambig
    - minscore
    - scores

    Returns:

    Example:
    """
    # print ("trancounts", trancounts)
    # print ("minscore", minscore)
    minreadlen = 15
    maxreadlen = 150

    t_trancounts = trancounts["unambig"]  # Unambig

    if ambig:
        t_trancounts = (
            pl.concat([t_trancounts, trancounts["ambig"]])
            .groupby("readlen", "pos")
            .agg(pl.sum("count"))
            .filter(
                (pl.col("readlen") >= minreadlen) & (pl.col("readlen") <= maxreadlen)
            )
        )
    else:
        if not t_trancounts.is_empty():
            t_trancounts = t_trancounts.filter(
                (pl.col("readlen") >= minreadlen) & (pl.col("readlen") <= maxreadlen)
            )
        if t_trancounts.is_empty():
            return pl.DataFrame({"pos": [], "count": []})

    t_trancounts = (
        t_trancounts.join(offsets_5p_offsetsNscores, on="readlen", how="left")
        .fill_null(14)
        .with_columns(pos=pl.col("pos") + pl.col("offset") + 1)
        .select("pos", "count")
        .groupby("pos")
        .agg(pl.sum("count"))
    )

    return t_trancounts


# Builds a profile, applying offsets
def build_proteomics_profile(
    trancounts: Dict[str, Dict[int, List[int]]],
    # , ambig
) -> Dict[int, int]:
    """

    Parameters:
    - trancounts

    Returns:

    Example:
    """
    minreadlen = 15
    maxreadlen = 150
    profile = {}
    t_trancounts = (
        tran_count(trancounts)
        .filter(pl.col("readlen") >= minreadlen and pl.col("readlen") <= maxreadlen)
        .with_columns(count=pl.col("count") / pl.col("readlen") / 3.0)
    )
    profile = []
    for row in t_trancounts.iter_rows(names=True):

        for pos in range(row["pos"], row["pos"] + row["readlen"], 3):
            profile.append([pos, row["count"]])
    profile = (
        pl.DataFrame(profile, schema=["pos", "count"])
        .groupby("pos")
        .agg(pl.sum("count"))
    )
    return profile


def fetch_filename_file_id(file_id: int) -> str:
    """
    Return the filename from the database given a file id.

    Parameters:
    - file_id (int): id of file

    Returns:
    - str

    Example:
    """
    return get_table("files").filter(pl.col("file_id") == file_id)[0, "file_name"]



# Sequence to RDG

def extract_translons(
    sequence: str,
    starts: set[str] = {"ATG", "CTG", "GTG"},
    min_length: int = 30,
) -> list[tuple[int, int]]:
    """
    Extract open reading frames (translons) from a nucleotide sequence.

    Parameters:
    - sequence (str): The input nucleotide sequence.
    - starts (Set[str]): Set of start codons to initiate translon detection.
            Default: {"ATG", "CTG", "GTG"}.
    - min_length (int): Minimum length of translons to be included in
            the result. Default: 10.

    Returns:
    List[Tuple[int, int]]: A list of tuples representing the start and
                        stop positions of detected translons.

    Example:
    ```python
    sequence = "ATGCTAGCATGAATAG"
    translons = extract_translons(sequence)
    print(translons)
    # Output: [(0, 15)]
    ```
    """
    if not sequence:
        return []
    sequence = sequence.upper()
    if not starts:
        starts = {"ATG"}# {"ATG", "CTG", "GTG"}
    aug_index = pd.read_table("AUG.csv", comment="#")[["sequence"]].reset_index().rename(columns={"index": "rank"})
    print(aug_index.head())

    stops = {"TAA", "TAG", "TGA"}
    translons = []
    starts_poses = {0:[],1:[],2:[]}
    seq_len = len(sequence)
    start_positions = []
    max_translone_count = 1000
    translon_count = 0
    df = [] # frame, from, to, type
    rank  = -1
    last_stop_start_indexes = []


    for i in range(seq_len - 3):
        codon = sequence[i: i + 3]
        if codon in starts:
            if i < 6:
                rank = -1
            else:
                rank = aug_index.loc[aug_index["sequence"] == sequence[i-6:i+5], "rank"].values[0] + 1
            starts_poses[i%3].append((i, rank))
            start_positions.append(i)

        elif codon in stops and starts_poses[i%3]:
            for pos, rank in starts_poses[i%3]:
                if (i-pos+3) >= min_length:
                    translons.append((pos, i+2))
                    start_index = start_positions.index(pos)
                    if start_index !=0:
                        df.append([ i%3, pos,[pos, i+2],[ start_positions[start_index-1], pos],[ i+2, seq_len], rank ])
                    else:
                        df.append([ i%3, pos,[pos, i+2],[ 0, pos],[ i+2, seq_len], rank ])

                    translon_count += 1
                    if translon_count == max_translone_count:
                        break
                else:
                    start_index = start_positions.index(pos)
                    start_positions = start_positions[:start_index]
                    break
            starts_poses[i%3] = []
        if translon_count == max_translone_count:
            break
    if translon_count == max_translone_count:
        return translons, df
    for start_pos_list in starts_poses.values():
        for start, rank in start_pos_list:
            if seq_len - start >= min_length:
                translons.append((start, seq_len-1))
                df.append([ start%3, start,[start, seq_len],[ start_positions[start_positions.index(start)-1] , start],[ seq_len, seq_len], rank ])
                translon_count += 1
                if max_translone_count == translon_count:
                    break
            else:
                break


    return translons, df

def sequence2rdg(sequence):
    translons = extract_translons(sequence,starts={"ATG"},min_length=10)
    df = pd.DataFrame(translons[1], columns=['frame',"start","cds", "5utr","3utr","rank"]).sort_values('start', ignore_index=True).reset_index().rename(columns={'index':'order'})
    df = df.melt(id_vars=["order","rank"], value_vars=["cds", "5utr","3utr"]).rename(columns={'variable':'frag'})

    df['x1'] = df['value'].apply(lambda x:x[0])
    df['x2'] = df['value'].apply(lambda x:x[1])
    df["lw"] = 1
    df.loc[df['frag']=='cds', 'lw'] = 5
    df["org_order"] = df["order"]
    # del df['value']

    return df.to_json(orient="records")
