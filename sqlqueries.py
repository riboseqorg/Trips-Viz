from sqlalchemy.orm import Session
from sqlalchemy.ext.automap import automap_base
from sqlalchemy import create_engine
import pandas as pd


def sqlquery(sqlfilepath, tablename):
    engine = create_engine('sqlite:///' + sqlfilepath)
    Base = automap_base()
    Base.prepare(engine, reflect=True)
    session = Session(engine)
    sql = session.query(Base.classes[tablename])
    sqlTable = pd.read_sql(sql, engine)
    return sqlTable


def table2dict(table, keys_col, values_col):
    dict = {}
    for index, row in table.iterrows():
        dict[row[keys_col]] = row[values_col]
    return dict
