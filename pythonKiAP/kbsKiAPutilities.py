"""
#
# Collection of functions which might come in-handy in various projects and setups.
#
"""
#-------------------------------------------------------------------------------------------------------------
# importing modules
import kbsKiAPutilities as kku
from importlib import reload
import os
import sys
import pandas as pd
import pyodbc as pySQL

#=======================================================================================================================
def returnSQLdatapull(SQLquery, server='10.19.31.20', database='opslagsdata', savefilename=None, overwrite=False, verbose=True):
    """
    Simple function to pull data from an SQL server provided an SQL code.
    The code will be sent to server and database which are provided at input parameters to function.

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

    Example of use
    -------
    import kbsKiAPutilities as kku
    SQLquery = "select top 100 * from opslagsdata.dbo.icpc_icd10_mapning"
    dataframe = kku.returnSQLdatapull(SQLquery)

    """
    if verbose: print(' - Defining server and database names')

    if verbose: print(' - Connecting to server')
    cnxn = pySQL.connect('DRIVER={ODBC Driver 17 for SQL Server}; SERVER='+server+'; DATABASE='+database+'; Trusted_Connection=yes;')

    if verbose: print(' - Executing SQL query')
    SQLresult = pd.read_sql_query(SQLquery, cnxn)
    dataframe = pd.DataFrame(SQLresult)
    # --------------------------------------------------------------------------------
    if savefilename is not None:
        if verbose: print(' - Saving extracted SQL data in Pandas DataFrame to file')
        kku.savedataframe2file(dataframe, savefilename, format='excel', overwrite=overwrite, verbose=verbose)

    #--------------------------------------------------------------------------------
    if verbose: print(' - Returning pulled SQL data')
    return dataframe
#=======================================================================================================================
def savedataframe2file(dataframe,out_file,format='CSV',overwrite=False,verbose=True):
    """
    Saving Pandas DataFrame as Excel or CSV file

    Parameters
    ----------
    dataframe      Pandas DataFrame to store in file
    out_file       Path and name to output Excel/CSV file to generate
    format         choose CSV or Excel file

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
