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

from sys import stdout as sysstdout
import Emento_datakombination as ed

# -----------------------------------------------------------------------------------------------------------------------
def generate_Excel_output(ExcelOutputName,testing=False,verbose=True):
    """
    Main function loading datafiles, addeing and sorting data, combining files and writing Excel output.

    ---EXMAPLE OF USE---
    import Emento_datakombination as ed
    import datetime
    todaystring = datetime.datetime.strftime(datetime.date.today(),"%y%m%d")
    ExcelOutputName = 'Emento_datagrundlag_'+todaystring+'.xlsx'
    outputfile = ed.generate_Excel_output(ExcelOutputName,testing=True,verbose=True)

    """
    todaystring  = datetime.datetime.strftime(datetime.date.today(), "%y%m%d")
    outpath      = 'O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2022 - Emento forløbsapp/python input and output/'

    if verbose: print('-->Loading data files to combine (names hardcoded in Emento_datakombination.load_datafiles')
    data_base, data_contacts, data_emento, data_ementoKey = ed.load_datafiles(testing=testing,verbose=verbose)


    #if verbose: print('-->Combining SP data arrays')
    #data_SP = ed.combine_SP_dataarrays(data_base, data_contacts, verbose=verbose)
    data_SP = data_contacts  # 240314 Use all contacts matched to selected CPRs as SP reference data

    if verbose: print('-->Adding unique identifier to Emento data')
    data_Emento_wCPR = ed.addCPR2Emento(data_emento, data_ementoKey, verbose=verbose)
    if verbose: print('-->Saving Emento with uid to Excel')
    outputfilename_Emento = 'EmentoWuid_'+todaystring+'.xlsx'
    data_Emento_wCPR.to_excel(outpath+'/'+outputfilename_Emento, sheet_name="SP og Emento data")
    if verbose: print('   saved to '+outpath+'/'+outputfilename_Emento)

    if verbose: print('-->Adding Emento data to SP data')
    data_SPandEmento = ed.add_emento_columns(data_SP, data_Emento_wCPR, verbose=verbose)

    if verbose: print('-->Saving final dataframe to Excel')
    outputfilename = 'EmentoOgSPDatakombination_'+todaystring+'.xlsx'
    data_SPandEmento.to_excel(outpath+'/'+outputfilename, sheet_name="SP og Emento data")
    if verbose: print('   saved to '+outpath+'/'+outputfilename)

    return outputfilename
# -----------------------------------------------------------------------------------------------------------------------
def load_datafiles(testing=False,verbose=True):
    """
    Wrapper to load data from SP and Emento into arrays

    ---EXMAPLE OF USE---
    import Emento_datakombination as ed
    data_base, data_contacts, data_emento, data_ementoKey = ed.load_datafiles()

    """
    file_path = 'O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2022 - Emento forløbsapp/python input and output/'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading SP data base on diagnoses and procedures (data_base)')
    if testing:
        file_base = file_path + 'kontakter baseret på dia og pro 230608_testing.csv'
    else:
        file_base = file_path+'kontakter baseret på dia og pro 230608.csv' #<------------------------- Edit filename for updates
    data_base = pd.read_csv(file_base, delimiter=";", dtype=None, decimal=',')
    datearr = np.asarray([datetime.datetime.strptime(dstr, '%d-%m-%Y %H:%M') for dstr in data_base['Kontakt startdato Dato-tid']])
    date_min = np.min(datearr)
    date_max = np.max(datearr)
    if verbose: print('   The contact start dates fall between '+str(date_min)+' and '+str(date_max))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading SP data matched to the base data (data_contacts)')
    if testing:
        file_contacts = file_path + 'kontakter cpr match 230608_testing.csv'
    else:
        file_contacts = file_path+'kontakter cpr match 230608.csv' #<--------------------------------- Edit filename for updates
    data_contacts = pd.read_csv(file_contacts, delimiter=";", dtype=None, decimal=',')
    datearr = np.asarray([datetime.datetime.strptime(dstr, '%d-%m-%Y %H:%M') for dstr in data_contacts['Kontakt startdato Dato-tid']])
    date_min = np.min(datearr)
    date_max = np.max(datearr)
    if verbose: print('   The contact start dates fall between '+str(date_min)+' and '+str(date_max))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading Emento tableau csv file (data_emento)')
    file_emento = file_path+'tableau_data230608.csv' #<------------------------------------------- Edit filename for updates
    data_emento = pd.read_csv(file_emento, delimiter=",", dtype=None, decimal='.')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading Emento ID-CPR key for matching (data_ementoKey)')
    file_ementoKey = file_path+'uidMapCSV_230608.csv' #<------------------------------------------ Edit filename for updates
    data_ementoKey = pd.read_csv(file_ementoKey, delimiter=",", dtype='O', decimal='.')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Creating log with info on loaded data')
    todaystring = datetime.datetime.strftime(datetime.date.today(), "%y%m%d")
    fout = open(file_path+'/'+'Emento_datakombination_load_datafiles_log'+todaystring+'.txt','w')
    fout.write('Log over datafiler der blev brugt i Emento_datakombination.load_datafiles() da kombineret fil '+todaystring+' blev genereret.\n')
    fout.write('\nFølgende filer blev loadet og kombineret:\n')
    fout.write('    SP data af patientkontakter med diagnose og procedure restriktioner    : ' + file_base +' \n')
    fout.write('    SP data af alle kontakter for CPR numre i file_base                    : ' + file_contacts + ' \n')
    fout.write('    Emento Tableau datafile fra https://storage-bi.emento.dk bucket        : ' + file_emento + ' \n')
    fout.write('    Emento UniqueIdentifier-CPR nøgle fra https://sla-mtk.emento.dk/ server: ' + file_ementoKey + ' \n')
    fout.close()

    if verbose: print(' - Done loading SP data base on diagnoses and procedures (data_base)')
    return data_base, data_contacts, data_emento, data_ementoKey
# -----------------------------------------------------------------------------------------------------------------------
def add_emento_columns(data_SP,data_emento,verbose=True):
    """
    Function to add the columns from the Emento data to entries with dates
    after the first appearance of a given CPR in the Emento data.
    """
    SP_CPR = data_SP['Patient CPR-nr.']

    #SP_time = np.asarray([datetime.datetime.strptime(dstr, '%d-%m-%Y %H:%M') for dstr in data_SP['Behandlingskontakt starttidspunkt Dato-tid']])
    #SP_time = pd.to_datetime(data_SP['Behandlingskontakt starttidspunkt Dato-tid'])
    SP_time = pd.to_datetime(data_SP['Kontakt startdato Dato-tid'])
    #Em_time = np.asarray([datetime.datetime.strptime(dstr[:-5], '%Y-%m-%dT%H:%M:%S') for dstr in data_emento['createdAt']])
    Em_time_cph = pd.to_datetime(data_emento['createdAt'],utc=True).map(lambda x: x.tz_convert('Europe/Copenhagen'))
    Em_time = pd.to_datetime(Em_time_cph.map(lambda x: x.tz_localize(None)))

    NlinesSP = len(SP_CPR)
    multiEmentoCPRlist = [] # list to contain CPRs for cases with more than 1 Emento ID associated with them

    for cc, col_Emento in enumerate(data_emento.columns):
        if verbose: print(' - Adding Emento data column : '+col_Emento)
        if data_emento[col_Emento].dtype == 'O':
            col_arr = np.asarray([str('%100s' % '')]*NlinesSP)
            stringcolumn = True
        else:
            col_arr = np.zeros(NlinesSP)*np.nan
            stringcolumn = False

        for ss, CPR_SP in enumerate(SP_CPR):
            if verbose:
                infostr = '   checking SP CPR number ' + str(ss + 1) + ' / ' + str(NlinesSP)+' for column '+str(cc+1)+'/'+str(len(data_emento.columns))
                sysstdout.write("%s\r" % infostr)
                sysstdout.flush()

            ent_cpr_match = np.where(data_emento['uniqueIdentifier'] == CPR_SP)[0]
            Nmatches = len(ent_cpr_match)

            if Nmatches == 1:
                if SP_time[ss] > Em_time[ent_cpr_match[0]]:
                    if stringcolumn:
                        col_arr[ss] = str(data_emento[col_Emento][ent_cpr_match[0]]).strip()
                    else:
                        col_arr[ss] = data_emento[col_Emento][ent_cpr_match[0]]
            elif Nmatches > 1:
                #if verbose: print('    Found '+str(Nmatches)+' matches to cpr '+str(CPR_SP)+' for ids '+idlist+' in Emento file')
                if cc == 0: # only count CPRs for first column
                    multiEmentoCPRlist.append(CPR_SP)
                idlist = str([data_emento[col_Emento][matchent] for matchent in ent_cpr_match])
                min_Em_time = np.min(Em_time[ent_cpr_match])
                if SP_time[ss] > min_Em_time: # comparing SP date with oldest Emento date
                    if stringcolumn:
                        col_arr[ss] = str(data_emento[col_Emento][ent_cpr_match[0]]).strip() # assigning data from first appearance of Emento match
                    else:
                        col_arr[ss] = data_emento[col_Emento][ent_cpr_match[0]]  # assigning data from first appearance of Emento match

        NmultiEmentoCPRlist = len(np.unique(np.asarray(multiEmentoCPRlist)))
        data_SP[col_Emento] = col_arr
    if verbose: print(' done ... found that '+str(NmultiEmentoCPRlist)+' CPRs had more than 1 emento ID associated with them')
    return data_SP
# -----------------------------------------------------------------------------------------------------------------------
def addCPR2Emento(data_emento, data_ementoKey, verbose=True):
    """
    Add CPR column to Emento data

    ---EXAMPLE OF USE---
    import Emento_datakombination as ed
    data_Emento_wCPR = ed.addCPR2Emento(data_emento, data_ementoKey, verbose=True)

    """
    if verbose: print(' - Matching CPR from data_ementoKey file to info from data_emento')
    data_emento_wCPR        = data_emento.copy()
    CPRcol = np.asarray(['']*len(data_emento_wCPR))

    Nnomatch = 0 # counter to keep track of the numer of Emento IDs without matches.
    for dd, Eid in enumerate(data_emento_wCPR['id']):
        ent_match = np.where(data_ementoKey['courseId'].values == Eid)[0]
        Nmatches = len(ent_match)

        if Nmatches == 0:
            #if verbose: print('   The following Emento ID did not have a match in data_ementoKey: '+str(Eid))
            Nnomatch = Nnomatch + 1
        elif Nmatches == 1:
            cprstr_last4  = str(data_ementoKey['uniqueIdentifier'][ent_match[0]])[-4:]
            cprstr_first6 = str("%06.d" % int(str(data_ementoKey['uniqueIdentifier'][ent_match[0]])[:-4]))
            cprstr = cprstr_first6+'-'+cprstr_last4
            CPRcol[dd] = cprstr

        else:
            if verbose: print('   The following Emento ID had --'+str(Nmatches)+'-- mastches in data_ementoKey: ' + str(Eid))

    data_emento_wCPR['uniqueIdentifier'] = CPRcol
    if verbose: print(' - There were '+str(Nnomatch)+' Emento IDs which did not have a match in the ID2CPR list')

    return data_emento_wCPR
# -----------------------------------------------------------------------------------------------------------------------
def combine_SP_dataarrays(data_base,data_contacts,verbose=True):
    """
    Function to combine the two dataarrays from SP

    ---EXAMPLE OF USE---
    import Emento_datakombination as ed
    data_SP = ed.combine_SP_dataarrays(data_base,data_contacts)

    """
    ktk_suffix = ' alle ktk'
    # manually adding column with Pt. record IDs to prevent it disappearing in merge
    data_contacts['Patientkontakt record ID'+' alle ktk'] = data_contacts['Patientkontakt record ID']

    if verbose: print(' Starting pandas merger of the two tables')
    #data_SP = data_base.merge(data_contacts, left_on='Patientkontakt record ID', right_on='Patientkontakt record ID',suffixes=('', ktk_suffix), how='left', copy=False)
    data_SP = pd.merge(data_base, data_contacts, on='Patientkontakt record ID', how='left')
    if verbose: print(' Done merging and returning the combined table')

    return data_SP
#=======================================================================================================================

