#from importlib import reload
import pdb
import sys
import os

import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
from matplotlib.ticker import PercentFormatter
import seaborn
import datetime
import loadMDCgroups
import scipy

import Emento_datakombination as ed

# -----------------------------------------------------------------------------------------------------------------------
def generate_Excel_output(ExcelOutputName,verbose=True):
    """
    Main function loading datafiles, addeing and sorting data, combining files and writing Excel output.

    ---EXMAPLE OF USE---
    import Emento_datakombination as ed
    import datetime
    todaystring = datetime.datetime.strftime(datetime.date.today(),"%y%m%d")
    ExcelOutputName = 'Emento_datagrundlag_'+todaystring+'.xlsx'
    outputfile = generate_Excel_output(ExcelOutputName,verbose=True)

    """


    data_base, data_contacts, data_emento, data_ementoKey = ed.load_datafiles(verbose=verbose)

    data_SP = ed.combine_SP_dataarrays(verbpse=verbose)

    data_Emento_wCPR = ed.addCPR2Emento(data_emento, data_ementoKey, verbose=verbose)

    data_wEmento = ed.add_emento_columns(data_SP, data_Emento, verbose=verbose)

    outputfile = ed.save_Excel_file(dataarray, verbose=verbose)

    return outputfile
# -----------------------------------------------------------------------------------------------------------------------
def load_datafiles(verbose=True):
    """
    Wrapper to load data from SP and Emento into arrays

    ---EXMAPLE OF USE---
    import Emento_datakombination as ed
    data_base, data_contacts, data_emento, data_ementoKey = ed.load_datafiles()

    """
    file_path = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Emento/python input and output/'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading SP data base on diagnoses and procedures (data_base)')
    file_base = file_path+'kontakter baseret på dia og pro.csv'
    data_base = pd.read_csv(file_base, delimiter=";", dtype=None, decimal=',')
    datearr = np.asarray([datetime.datetime.strptime(dstr, '%d-%m-%Y %H:%M') for dstr in data_base['Kontakt startdato Dato-tid']])
    date_min = np.min(datearr)
    date_max = np.max(datearr)
    if verbose: print('   The contact start dates fall between '+str(date_min)+' and '+str(date_max))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading SP data matched to the base data (data_contacts)')
    file_contacts = file_path+'kontakter cpr match.csv'
    data_contacts = pd.read_csv(file_contacts, delimiter=";", dtype=None, decimal=',')
    datearr = np.asarray([datetime.datetime.strptime(dstr, '%d-%m-%Y %H:%M') for dstr in data_contacts['Kontakt startdato Dato-tid']])
    date_min = np.min(datearr)
    date_max = np.max(datearr)
    if verbose: print('   The contact start dates fall between '+str(date_min)+' and '+str(date_max))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading Emento tableau csv file (data_emento)')
    file_emento = file_path+'tableau_data230213.csv'
    data_emento = pd.read_csv(file_emento, delimiter=",", dtype=None, decimal='.')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading Emento ID-CPR key for matching (data_ementoKey)')
    file_ementoKey = file_path+'uidMapCSV_230214.csv'
    data_ementoKey = pd.read_csv(file_ementoKey, delimiter=",", dtype=None, decimal='.')

    return data_base, data_contacts, data_emento, data_ementoKey
# -----------------------------------------------------------------------------------------------------------------------
def add_emento_columns(data_SP,data_Emento,verbose=True):
    """
    Function to add the columns from the Emento data to entries with dates
    after the first apperance of a given CPR in the Emento data.
    """

    data_SP = np.concatenate((data_base,data_contacts))

    return data_SP


# -----------------------------------------------------------------------------------------------------------------------
def addCPR2Emento(data_emento, data_ementoKey, verbose=verbose):
    """
    Add CPR column to Emento data

    ---EXAMPLE OF USE---
    import Emento_datakombination as ed
    data_Emento_wCPR = ed.addCPR2Emento(data_emento, data_ementoKey, verbose=True)

    """
    if verbose: print(' - Matching CPR from data_ementoKey file to info from data_emento')
    data_emento_wCPR        = data_emento
    data_emento_wCPR['CPR'] = np.zeros(len(data_emento_wCPR))

    for dd, Eid in enumerate(data_emento_wCPR['id']):
        ent_match = np.where(data_ementoKey['courseId'] == Eid)[0]
        Nmatches = len(ent_match)
        if Nmatches == 0:
            if verbose: print('   The following Emento ID did not have a match in data_ementoKey: '+str(Eid))
        elif len(Nmatches) == 1:
            cprstr_last4  = str(data_ementoKey['uniqueIdentifier'][ent_match[0]])[-4:]
            cprstr_first6 = str("%.d6" % str(data_ementoKey['uniqueIdentifier'][ent_match[0]])[:-4])
            cprstr = cprstr_first6+'-'+cprstr_last4
            data_emento_wCPR['CPR'] = cprstr
        else:
            if verbose: print('   The following Emento ID had --'+str(Nmatches)+'-- mastches in data_ementoKey: ' + str(Eid))

    return data_emento_wCPR
# -----------------------------------------------------------------------------------------------------------------------
def combine_SP_dataarrays(data_base,data_contacts,verbose=True):
    """
    Function to combine the two dataarrays from SP

    ---EXAMPLE OF USE---
    import Emento_datakombination as ed
    data_SP = ed.combine_SP_dataarrays(data_base,data_contacts)

    """

    data_SP = np.concatenate((data_base,data_contacts))

    return data_SP
# -----------------------------------------------------------------------------------------------------------------------
def save_Excel_file(dataarray, verbose=True):
    """
    Function to store the final output Excel file, that can be used for further analysis
    """
    outputname = "PathToFile..."

    return outputname


#=======================================================================================================================

