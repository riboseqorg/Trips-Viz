from typing import Dict, List, Tuple, Union
from sqlalchemy.orm import Session
from sqlalchemy.ext.automap import automap_base
from sqlalchemy import create_engine, insert
import config
import pandas as pd


def _sql(sqlfilepath: str, tablename: str) -> Tuple[Session, Engine]:
    engine = create_engine('sqlite:///' + sqlfilepath)
    Base = automap_base()
    Base.prepare(engine, reflect=True)
    session = Session(engine)
    sql = session.query(Base.classes[tablename])
    return sql, engine


def sqlquery(sqlfilepath: str, tablename: str) -> pd.DataFrame:
    sql, engine = _sql(sqlfilepath, tablename)
    sqlTable = pd.read_sql(sql, engine)
    return sqlTable


def sqldict2table(sqldict: Dict, tablename: str) -> pd.DataFrame:
    return pd.DataFrame(sqldict)


def get_user_id(username: str) -> int:
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


def update_table(table: pd.DataFrame, dct: Dict, task: str = 'delete') -> None:
    sqlfilepath = '{}/{}'.format(config.SCRIPT_LOC, config.DATABASE_NAME)
    sql, engine = _sql(sqlfilepath, table)
    # TODO: Take case of list values
    if task == 'delete':
        sql.delete(dct)
    if task == 'insert':
        sql.insert(dct)
    else:
        sql.update(dct)
    sql.commit()


def form2filtered_data(data: pd.DataFrame, form: Dict) -> pd.DataFrame:
    '''Returns a dataframe based on the given form filters'''
    form_keys = set(form.keys()) & set(data.columns)
    # Add filters
    return data
