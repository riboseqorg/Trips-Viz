from typing import Dict, List, Tuple, Union, Any, Hashable  #, Unknown
from typing_extensions import Literal
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


def table2dict(table: pd.DataFrame, keys: List[str]) -> Dict[Hashable, Any]:
    '''
    Convert a table to a dictionary of lists. 
    >>> data = {'key1': [1, 2, 3], 'key2': [4, 5, 6], 'key3': [7, 8, 9], 
    ... 'key4': [10, 11, 12]}
    >>> table = pd.DataFrame(data)
    >>> table2dict(table, ['key1', 'key2', 'key3'])
    >>> {1:{4:{7:[10]}}, 2:{5:{8:[11]}}, 3:{6:{9:[12]}}}

    '''
    return {
        n: grp.loc[n].to_dict('index')
        for n, grp in table.set_index(keys).groupby(level=keys[0])
    }


def update_table(table: str,
                 task: Literal['insert', 'update', 'delete'] = 'insert',
                 values: Dict[str, Any] = {},
                 where: Dict[str, Any] = {}) -> None:
    '''
    Update a table with the given data. 
    >>> update_table('users', "update",{'user_id': 1},  {'user_id': 2}) 
    >>> update_table('users', "delete",{},  {'user_id': 2}) 
    >>> update_table('users', "insert",  {'user_id': 2}) 


    '''
    sqlfilepath = '{}/{}'.format(config.SCRIPT_LOC, config.DATABASE_NAME)
    engine = create_engine('sqlite:///' + sqlfilepath)
    metadata = MetaData()
    table_query = Table(table, metadata, autoload=True, autoload_with=engine)
    with engine.connect() as conn:

        query: Query | Any = ""

        if task == 'insert' and values:
            query = insert(table_query).values(**values)

        elif task == 'update' and values and where:
            query = update(table_query).where(**where).values(**values)
        elif task == 'delete' and where:
            query = delete(table_query).where(**where)
        if query:
            conn.execute(query)


def form2filtered_data(data: pd.DataFrame, form: Dict) -> pd.DataFrame:
    '''Returns a dataframe based on the given form filters'''
    form_keys = set(form.keys()) & set(data.columns)
    # Add filters
    return data
