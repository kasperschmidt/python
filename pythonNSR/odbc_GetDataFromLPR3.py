import sys
import os
import pyodbc
import pandas as pd
import odbc_GetDataFromLPR3 as gdf
from importlib import reload

#-----------------------------------------------------------------------------------------------------------------------
def returndatapull(SQLquery, server='sv1391', database='Databank_COMPLIANCE', savefilename=None, overwrite=False, verbose=True):
    """
    Simple function to pull out data from LPR3 provided and SQL code
    The code will be sent to server and database given in parameters

    Parameters
    ----------
    SQLquery        The SQL query to push to the SQL server and retrieve data based on.
                    Provided as a string containing the full SQL query
    server          The server containing the data to pull the data from
    database        The database to pull the data from
    savefilename    To save the data pull as an Excel og CSV file provide a filename
    overwrite       If true, then overwrite savefilename if it already exists
    verbose         Toggle verbosity of function

    Returns
    -------
    Pandas dataframe with the results from the SQL query and potentially a saved file (Excel og CSV).

    Example of yse
    -------
    # sys.path.append('C:/Users/kaschm/GitHub/python/pythonNSR/')
    import odbc_GetDataFromLPR3 as gdf
    SQLquery = "select top 100 * from DRG_LPR3.drgKontakter_alt where SYGEHUS_REGION_Tekst = 'Region Sjælland' "
    dataframe = gdf.runSQL(SQLquery)

    import odbc_basicSQLextract as obe
    dataframe = obe.returndatapull(printdata=True)

    """
    if verbose: print(' - Defining server and database names')

    if verbose: print(' - Connecting to server')
    cnxn = pyodbc.connect('DRIVER={ODBC Driver 17 for SQL Server};'
                          'SERVER='+server+';'
                          'DATABASE='+database+';'
                          'Trusted_Connection=yes;')
    #--------------------------------------------------------------------------------
    if verbose: print(' - Executing SQL query')
    SQLresult = pd.read_sql_query(SQLquery, cnxn)
    dataframe = pd.DataFrame(SQLresult)
    # --------------------------------------------------------------------------------
    if savefilename is not None:
        if verbose: print(' - Saving extracted SQL data in Pandas DataFrame to file')
        gdf.savefile(dataframe, savefilename, format='excel', overwrite=overwrite, verbose=verbose)

    #--------------------------------------------------------------------------------
    if verbose: print(' - Returning pulled SQL data')
    return dataframe
#=======================================================================================================================
def savefile(dataframe,out_file,format='CSV',overwrite=False,verbose=True):
    """
    Saving Pandas DataFrame as Excel or CSV file

    Parameters
    ----------
    dataframe      Pandas DataFrane to store in file
    out_file       Path and name to output Excel/CSV file to generate
    format         chose CSV or Excel file

    """
    if format.lower() == 'csv':
        if not out_file.endswith(".csv"):
            out_file = out_file+'.csv'
    if format.lower() == 'excel':
        if not out_file.endswith(".xlsx"):
            out_file = out_file+'.xlsx'

    if os.path.isfile(out_file) and (overwrite is False):
        sys.exit('The output file '+out_file+' already exists and overwrite=False ')
    elif os.path.isfile(out_file) and (overwrite is True):
        if verbose: print(' - The output file ' + out_file + ' already exists but overwrite=True ')

    if verbose: print(' - Storing DataFrame in '+out_file)
    if format.lower() == 'csv':
        dataframe.to_csv(out_file)
    elif format.lower() == 'excel':
        dataframe.to_excel(out_file, sheet_name="data output")
#=======================================================================================================================
def loadCSV(csvfile,verbose=True):
    """
    Loading CSV file

    Parameters
    ----------
    csvfile      CSV file to load into Pandas DataFrame
    verbose      Toggle verbosity
    """
    if verbose: print(' - Loading CSV data file '+csvfile+' into Pandas DataFrame')
    dataframe = pd.read_csv(csvfile)
    return dataframe
#=======================================================================================================================
def loadExcel(excelfile,verbose=True):
    """
    Loading CSV file

    Parameters
    ----------
    excelfile    Excel file to load into Pandas DataFrame
    verbose      Toggle verbosity
    """
    if verbose: print(' - Loading excel data file '+excelfile+' into Pandas DataFrame')
    dataframe = pd.read_excel(excelfile)
    return dataframe
#=======================================================================================================================

