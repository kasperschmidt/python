import sys
import os

import numpy as np
import pyodbc
import pandas as pd
import odbc_basicSQLextract as obe
#from importlib import reload

#sys.path.append('C:/Users/kaschm/GitHub/python/pythonNSR/')
import odbc_GetDataFromLPR3 as gdf
import lungemed_beddays_and_visitations as lbv

#-----------------------------------------------------------------------------------------------------------------------
def loadSQL_beddays(filepath='O:\Administration\\02 - Økonomi og PDK\Medarbejdermapper\Kasper\Focus1 - Ad hoc opgaver\Lungemed sengedage og visitationer\Lungemed.sql'):
    """
    Loading SQL query in file "filepath"

    """
    content = open(filepath, 'r').read()
    return content
# -----------------------------------------------------------------------------------------------------------------------
def loadSQL_visitations(filepath='O:\Administration\\02 - Økonomi og PDK\Medarbejdermapper\Kasper\Focus1 - Ad hoc opgaver\Lungemed sengedage og visitationer\Lungemed_visitationsoprindelse_nogroup.sql'):
    """
    Loading SQL query in file "filepath"

    """
    content = open(filepath, 'r').read()
    return content
# -----------------------------------------------------------------------------------------------------------------------
def getdata(verbose=True):
    """
    Function returning the data structures for the beddays and visitations SQL queries loaded in the load functions.

    Example of yse
    -------
    import lungemed_beddays_and_visitations as lbv
    dataframe_days, dataframe_vis = lbv.getdata()

    """

    if verbose: print(' - Getting LPR3 data from parsing SQL query for "bed days" ')
    savefilename_days  = None
    overwrite_days     = False
    dataframe_days     = gdf.returndatapull(lbv.loadSQL_beddays(), verbose=verbose, savefilename=savefilename_days, overwrite=overwrite_days)

    if verbose: print(' - Getting LPR3 data from parsing SQL query for "visitations" ')
    savefilename_vis   = None
    overwrite_vis      = False
    dataframe_vis      = gdf.returndatapull(lbv.loadSQL_visitations(), verbose=verbose, savefilename=savefilename_vis, overwrite=overwrite_vis)

    return dataframe_days, dataframe_vis
#=======================================================================================================================

