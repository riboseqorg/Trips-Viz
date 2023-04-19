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


def table2dict(table, keys):
    return {
        n: grp.loc[n].to_dict('index')
        for n, grp in table.set_index(keys).groupby(level=keys[0])
    }
