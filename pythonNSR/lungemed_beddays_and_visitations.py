import pdb
import sys
import os

import numpy as np
import pyodbc
import pandas as pd
import odbc_basicSQLextract as obe
#from importlib import reload
import datetime

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

# -----------------------------------------------------------------------------------------------------------------------
def count_occurrences_per_day(verbose=True):
    """
    Function to count occurrences per day used for bed occupancy

    Example of use
    -------
    import lungemed_beddays_and_visitations as lbv
    df_results = lbv.count_occurrences_per_day()

    """
    if verbose: print('- Getting the data to look at ')
    dataframe_days, dataframe_vis = lbv.getdata(verbose=verbose)

    measurehour = 23
    start_day   = datetime.datetime.strptime("02-02-2019 "+str(23)+":00:00", "%d-%m-%Y %H:%M:%S")
    end_day     = datetime.datetime.strptime("02-03-2019 "+str(23)+":00:00", "%d-%m-%Y %H:%M:%S")
    #end_day     = datetime.datetime.today()
    date_list   = [start_day + datetime.timedelta(days=x) for x in range(0, (end_day - start_day).days)]
    count_list  = [0] * len(date_list)
    occupancy_available  = [0] * len(date_list)
    occupancy_actual     = [0] * len(date_list)


    if verbose: print(' - Counting how many patients are in beds at any given day between '+
                      start_day.strftime("%d-%m-%Y")+' and '+start_day.strftime("%d-%m-%Y")+' at '+str(23)+" o'clock")
    for pp, patient in enumerate(dataframe_days['CPR']):
        intime  = dataframe_days['INDTIDSPUNKT_DRGKONTAKT'][pp]
        outtime = dataframe_days['UDTIDSPUNKT_DRGKONTAKT'][pp]

        for dd, datecheck in enumerate(np.asarray(date_list)):
            if verbose:
                infostr = '   Checking the date '+date_list[0].strftime("%d-%m-%Y")+' for patient number '+str(pp+1)
                sys.stdout.write("%s\r" % infostr)
                sys.stdout.flush()

            if (intime <= datecheck) and (datecheck <= outtime):
                count_list[dd] = count_list[dd] + 1

    if verbose: print('\n - Estimating the occupancy in the available and actual beds ')
    for dd, datecheck in enumerate(np.asarray(date_list)):
        if datecheck < datetime.datetime.strptime("10-06-2021 00:00:00", "%d-%m-%Y %H:%M:%S"):
            occupancy_available[dd] = count_list[dd] / 24.*100
        else:
            occupancy_available[dd] = count_list[dd] / 16. * 100

        if datecheck < datetime.datetime.strptime("01-03-2021 00:00:00", "%d-%m-%Y %H:%M:%S"):
            occupancy_actual[dd] = count_list[dd] / 24.*100
        else:
            occupancy_actual[dd] = count_list[dd] / 16. * 100


    if verbose: print('\n - returning count of patients')
    dict = {'date':date_list, 'count':count_list, 'occupancy_available':occupancy_available, 'occupancy_actual':occupancy_actual}
    df_results = pd.DataFrame(dict)

    return df_results
#=======================================================================================================================

