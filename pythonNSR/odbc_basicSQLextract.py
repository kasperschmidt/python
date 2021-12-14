import pyodbc
import pandas as pd
from importlib import reload
#-----------------------------------------------------------------------------------------------------------------------
def returndatapull(printdata=False,withpandas=True,verbose=True):
    """
    Parameters
    ----------
    printdata
    verbose

    Returns
    -------

    Example of yse
    -------
    # sys.path.append('C:/Users/kaschm/GitHub/python/pythonNSR/')
    import odbc_basicSQLextract as obe
    datapull = obe.returndatapull(printdata=True)

    """
    if verbose: print(' - Defining server and database names')
    server = 'sv1391'
    database = 'Databank_COMPLIANCE'

    if verbose: print(' - Connecting to server')
    cnxn = pyodbc.connect('DRIVER={ODBC Driver 17 for SQL Server};'
                          'SERVER='+server+';'
                          'DATABASE='+database+';'
                          'Trusted_Connection=yes;')
    #--------------------------------------------------------------------------------
    #                                  SQL Query
    if verbose: print(' - Executing SQL query')
    if withpandas:
        datapull = pd.read_sql_query("""
        select top 100 *
        from DRG_LPR3.drgKontakter_alt
        where SYGEHUS_REGION_Tekst = 'Region Sjælland'
        """, cnxn)

    else:
        datapull = cnxn.cursor()
        datapull.execute("""
        select top 100 *
        from DRG_LPR3.drgKontakter_alt
        where SYGEHUS_REGION_Tekst = 'Region Sjælland'
        """)
    #--------------------------------------------------------------------------------
    if printdata:
        if verbose: print(' - Printing data as "printdata"=True ')
        for ii in datapull:
            print(ii)
            if withpandas:
                print(datapull[['CPR', 'AKTIVITETSAAR']])
                pd.DataFrame(datapull, index=[80], columns=['CPR', 'AKTIVITETSAAR'])

    #--------------------------------------------------------------------------------
    if verbose: print(' - Returning pulled SQL data')
    return datapull
#-----------------------------------------------------------------------------------------------------------------------

