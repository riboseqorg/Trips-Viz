from sqlalchemy.orm import Session
from sqlalchemy.ext.automap import automap_base
from sqlalchemy import create_engine
import pandas as pd


def sqlquery(sqlfilepath, tablename):
    engine = create_engine('sqlite:///' + sqlfilepath)
    Base = automap_base()
    Base.prepare(engine, reflect=True)
    session = Session(engine)
    sqlTable = pd.read_sql(session.query(Base.classes[tablename]), engine)
    return sqlTable
