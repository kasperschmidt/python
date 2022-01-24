#from importlib import reload
import pdb
import sys
import os

import numpy as np
import datetime
import odbc_GetDataFromLPR3 as gdf

#=======================================================================================================================
#-----------------------------------------------------------------------------------------------------------------------
def return_ydernummer(cprlist,verbose=True):
    """
    Function returning the "ydernummer" for a list of CPR numbers.

    cprlist  = list of CPR numbers as strings to return ydernummer for

    cprlist = ['112233-8888','112233-9999']
    pull_ydernummer.return_ydernummer(cprlist)

    """
    cprlist_string = '('
    for cprnumber in cprlist:
        cprlist_string = cprlist_string + " '" + str(cprnumber).replace('-', '')+"',"
    cprlist_string = cprlist_string[:-1] + ')'

    querystring = """
    select *
    from Databank_COMPLIANCE.CPR.DTSIKREDE
    where
    Personnummer in
    """ + cprlist_string
    """
    order
    by
    Personnummer
    """

    if verbose: print(' - Getting LPR3 data from parsing SQL query ')
    dataframe = gdf.returndatapull(querystring, verbose=verbose)

    print(' cprnummer     ydernummer')
    for ii, cprnumber in enumerate(cprlist):
        cprent = np.where(dataframe['Personnummer'] == cprnumber.replace('-', ''))[0]
        if len(cprent) == 0:
            print(str(cprnumber) + '   MISSINGinDB ')
            #print('There is no entry with CPR number '+str(cprnumber)+' in returned dataframe... there should be')
            #pdb.set_trace()
        else:
            print(str(cprnumber)+'   '+str(int(dataframe['Ydernummer'][cprent[0]])))
#=======================================================================================================================

