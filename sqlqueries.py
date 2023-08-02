from typing import Dict, List, Tuple, Union, Any
from sqlalchemy.orm import Session
from sqlalchemy.ext.automap import automap_base
from sqlalchemy import create_engine, Table, MetaData  #, insert
import config
import pandas as pd
from pandas.core.frame import DataFrame
from sqlalchemy.orm.query import Query
from sqlalchemy.engine.base import Engine
from sqlalchemy import insert, delete, update
# from sqlalchemy.dialects.sqlite import insert, delete, update


def sqlquery(sqlfilepath: str, tablename: str) -> DataFrame:
    engine = create_engine('sqlite:///' + sqlfilepath).connect()
    table = pd.read_sql_table(table_name=tablename, con=engine)
    return table


def sqldict2table(sqldict: Dict, tablename: str) -> pd.DataFrame:
    return pd.DataFrame(sqldict)


def get_user_id(username: str | None) -> int:
    '''Return the user_id for a given username'''
    users = sqlquery('{}/{}'.format(config.SCRIPT_LOC, config.DATABASE_NAME),
                     'users')

    return users.loc[users['username'] == username, 'user_id'].values[0]


def get_table(table: str) -> pd.DataFrame:
    '''Return a table as a pandas dataframe.'''
    return sqlquery('{}/{}'.format(config.SCRIPT_LOC, config.DATABASE_NAME),
                    table)


def table2dict(table: pd.DataFrame, keys: List[str]) -> Dict[str, List[str]]:
    return {
        n: grp.loc[n].to_dict('index')
        for n, grp in table.set_index(keys).groupby(level=keys[0])
    }


def update_table(table: str,
                 dct: Dict[str, Any],
                 task: str = 'delete') -> None:
    sqlfilepath = '{}/{}'.format(config.SCRIPT_LOC, config.DATABASE_NAME)
    engine = create_engine('sqlite:///' + sqlfilepath)
    metadata = MetaData()
    table_query = Table(table, metadata, autoload=True, autoload_with=engine)
    with engine.connect() as conn:

        # TODO: Take case of list values
        query = ""

        if task == 'insert':
            if table == 'urls':
                query = insert(table_query).values(
                    url_id=dct['url_id'],
                    url=dct['url'],
                )
            elif table == 'users':
                print(dct)
                query = insert(table_query).values(
                    user_id=dct['user_id'],
                    username=dct['username'],
                    password=dct['password'],
                    study_access=dct['study_access'],
                    organism_access=dct['organism_access'],
                    advanced=dct['advanced'],
                    temp_user=dct['temp_user'])
            elif table == 'user_settings':
                query = insert(table_query).values(
                    marker_size=dct["marker_size"],
                    axis_label_size=dct["axis_label_size"],
                    subheading_size=dct["subheading_size"],
                    title_size=dct["title_size"],
                    user_id=dct["user_id"],
                    background_col=dct["background_col"],
                    readlength_col=dct["readlength_col"],
                    metagene_fiveprime_col=dct["metagene_fiveprime_col"],
                    metagene_threeprime_col=dct["metagene_threeprime_col"],
                    nuc_comp_a_col=dct["nuc_comp_a_col"],
                    nuc_comp_t_col=dct["nuc_comp_t_col"],
                    nuc_comp_g_col=dct["nuc_comp_g_col"],
                    nuc_comp_c_col=dct["nuc_comp_c_col"],
                    uga_col=dct["uga_col"],
                    uag_col=dct["uag_col"],
                    uaa_col=dct["uaa_col"],
                    comp_uga_col=dct["comp_uga_col"],
                    comp_uag_col=dct["comp_uag_col"],
                    comp_uaa_col=dct["comp_uaa_col"],
                    cds_marker_width=dct["cds_marker_width"],
                    cds_marker_colour=dct["cds_marker_colour"],
                    legend_size=dct["legend_size"],
                    ribo_linewidth=dct["ribo_linewidth"])
            elif table == 'studies':
                query = insert(table_query).values(
                    study_id=dct["study_id"],
                    organism_id=dct["organism_id"],
                    study_name=dct["study_name"],
                    paper_authors=dct["paper_authors"],
                    srp_nos=dct["srp_nos"],
                    paper_year=dct["paper_year"],
                    paper_pmid=dct["paper_pmid"],
                    paper_link=dct["paper_link"],
                    gse_nos=dct["gse_nos"],
                    adapters=dct["adapters"],
                    paper_title=dct["paper_title"],
                    description=dct["description"],
                    private=dct["private"],
                    owner=dct["owner"])
            elif table == 'study_access':
                query = insert(table_query).values(study_id=dct['study_id'],
                                                   user_id=dct['user_id'])
            elif table == 'deletions':
                query = insert(table_query).values(
                    study_id=dct['study_id'],
                    file_path=dct['file_path'],
                    deletion_date=dct['deletion_date'],
                )
            elif table == 'organism':
                query = insert(table_query).values(
                    organism_id=dct["organism_id"],
                    organism_name=dct["organism_name"],
                    transcriptome_list=dct["transcriptome_list"],
                    gwips_databasename=dct["gwips_databasename"],
                    gwips_clade=dct["gwips_clade"],
                    gwips_organism=dct["gwips_organism"],
                    gwips_database=dct["gwips_database"],
                    default_transcript=dct["default_transcript"],
                    private=dct["private"],
                    owner=dct["owner"])

            elif table == 'organism_access':
                query = insert(table_query).values(
                    organism_id=dct["organism_id"], user_id=dct["user_id"])

            elif table == 'files':
                query = insert(table_query).values(
                    file_id=dct["file_id"],
                    organism_id=dct["organism_id"],
                    study_id=dct["study_id"],
                    file_name=dct["file_name"],
                    file_description=dct["file_description"],
                    file_type=dct["file_type"],
                    owner=dct["owner"],
                    mapped_reads=dct["mapped_reads"],
                    control=dct["control"],
                    cell_line=dct["cell_line"])
            elif table == 'users_saved_cases':
                query = insert(table_query).values()
            else:
                pass
        if task == 'update':
            update(table_query).where(**change).values(**dct)
        conn.execute(query)  # TODO: Need to fix the type
        # delete(table_query).where()


def form2filtered_data(data: pd.DataFrame, form: Dict) -> pd.DataFrame:
    '''Returns a dataframe based on the given form filters'''
    form_keys = set(form.keys()) & set(data.columns)
    # Add filters
    return data
