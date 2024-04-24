#from importlib import reload
import pdb

import glob
import numpy as np
import pandas as pd
import datetime

from sys import stdout as sysstdout
import Emento_SPdataberigelse as esd

# -----------------------------------------------------------------------------------------------------------------------
def generate_Excel_output(contactstamp='230608',ementostamp='230608',testing=False,verbose=True):
    """
    Main function loading datafiles, adding and sorting data, combining files and writing Excel output.

    ---EXMAPLE OF USE---
    import Emento_SPdataberigelse as esd
    outputfile = esd.generate_Excel_output(contactstamp='230608',ementostamp='230608',testing=False,verbose=True)

    """
    todaystring  = datetime.datetime.strftime(datetime.date.today(), "%y%m%d")
    outpath      = 'O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2022 - Emento forløbsapp/python input and output/'

    if verbose: print('-->Loading data files to combine')
    data_SPcontacts, data_SPbase, data_emento, data_ementoKey = esd.load_datafiles(contactstamp=contactstamp, ementostamp=ementostamp,ementoonly=False, testing=testing, verbose=verbose)

    SPbase_cpr = np.unique(data_SPbase['Patient CPR-nr.'])
    SPcontacts_cpr = np.unique(data_SPcontacts['Patient CPR-nr.'])
    E_cpr  = np.unique(data_ementoKey['uniqueIdentifier'])
    if verbose: print(' - Found the following numer of CPRs in data:')
    if verbose: print('   SP basedata    : ' + str(len(SPbase_cpr)))
    if verbose: print('   SP contactdata : '+str(len(SPcontacts_cpr)))
    if verbose: print('   Emento data    : '+str(len(E_cpr)))

    if verbose: print('-->Adding unique identifier to Emento data')
    data_Emento_wCPR = esd.addCPR2Emento(data_emento, data_ementoKey, verbose=verbose)

    if verbose: print('-->Saving Emento with unique id (CPR) to Excel')
    if testing:
        outputfilename_Emento = 'EmentoWuid_'+todaystring+'_testing.xlsx'
    else:
        outputfilename_Emento = 'EmentoWuid_' + todaystring + '.xlsx'
    data_Emento_wCPR.to_excel(outpath+'/'+outputfilename_Emento, sheet_name="Emento data med CPR", index=False)
    if verbose: print('   saved to '+outpath+'/'+outputfilename_Emento)

    #
    # === Note on further development 240322: ===
    # To add more columns of Emento data loop over esd.load_datafiles with ementoonly=True
    # globbing the Emento tableau datafiles and pulling out Emento_column_key timestamps
    # (make sure there are matching Emento CPR-ID key files i datadirectory)
    # and then run esd.add_emento_columns() in the loop
    #
    if verbose: print('-->Adding Emento data to SP data')
    outname_extensions = ['SPallektk','SPbasis']
    for dd, data_SP in enumerate([data_SPcontacts, data_SPbase]):
        data_SPandEmento, multiEmentoCPRlist = esd.add_emento_columns(data_SP, data_Emento_wCPR, Emento_column_key='', verbose=verbose)
        if verbose: print('-->Saving final dataframe ('+outname_extensions[dd]+') to Excel')
        if testing:
            outputfilename = 'EmentoOgSPDatakombination_'+outname_extensions[dd]+'_'+todaystring+'_testing.xlsx'
        else:
            outputfilename = 'EmentoOgSPDatakombination_'+outname_extensions[dd]+'_' + todaystring + '.xlsx'
        if verbose: print(' - Sort output dataframe by CPR and contact date and time')
        data_SPandEmento['Kontakt startdato Dato-tid'] = pd.to_datetime(data_SPandEmento['Kontakt startdato Dato-tid'], format='%d-%m-%Y %H:%M') # ensure proper sorting
        data_SPandEmento = data_SPandEmento.sort_values(by=["Patient CPR-nr.", "Kontakt startdato Dato-tid"], ascending=[True, True]).copy()
        data_SPandEmento.to_excel(outpath+'/'+outputfilename, sheet_name="SP koblet med Emento data", index=False)
        if verbose: print('   saved to '+outpath+'/'+outputfilename)

    return outputfilename

# -----------------------------------------------------------------------------------------------------------------------
def load_datafiles(contactstamp='230608',ementostamp='230608',ementoonly=False,testing=False,verbose=True):
    """
    Wrapper to load data from SP and Emento into arrays

    ---EXMAPLE OF USE---
    import Emento_SPdataberigelse as esd
    data_contacts, data_emento, data_ementoKey = esd.load_datafiles(contactstamp='230608',ementostamp='230608',ementoonly=False,testing=False,verbose=True)

    """
    file_path = 'O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2022 - Emento forløbsapp/python input and output/'

    if not ementoonly:
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print(' - Loading SP contacts (data_contacts)')
        if testing:
            file_contacts = file_path + 'Kontakter for 10 ementopatiente.csv'
        else:
            file_contacts = file_path + 'kontakter cpr match ' + contactstamp+'.csv'
        data_contacts = pd.read_csv(file_contacts, delimiter=";", dtype=None, decimal=',')
        datearr = np.asarray([datetime.datetime.strptime(dstr, '%d-%m-%Y %H:%M') for dstr in data_contacts['Kontakt startdato Dato-tid']])
        date_min = np.min(datearr)
        date_max = np.max(datearr)
        if verbose: print('   The contact start dates fall between '+str(date_min)+' and '+str(date_max))
    else:
        file_contacts = '"Not loaded as ementoonly=True"'
        data_contacts = []

    if not ementoonly:
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print(' - Loading SP base contacts (data_base)')
        if testing:
            file_base = file_path + 'Kontakter for 10 ementopatiente.csv'
        else:
            file_base = file_path + 'kontakter baseret på dia og pro ' + contactstamp + '.csv'
        data_base = pd.read_csv(file_base, delimiter=";", dtype=None, decimal=',')
        datearr = np.asarray([datetime.datetime.strptime(dstr, '%d-%m-%Y %H:%M') for dstr in data_base['Kontakt startdato Dato-tid']])
        date_min = np.min(datearr)
        date_max = np.max(datearr)
        if verbose: print('   The base contact start dates fall between '+str(date_min)+' and '+str(date_max))
    else:
        file_base = '"Not loaded as ementoonly=True"'
        data_base = []

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading Emento tableau csv file (data_emento)')
    file_emento = file_path+'tableau_data_'+ementostamp+'.csv'
    data_emento = pd.read_csv(file_emento, delimiter=",", dtype=None, decimal='.')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading Emento ID-CPR key for matching (data_ementoKey)')
    file_ementoKey = file_path+'uidMapCSV_'+ementostamp+'.csv'
    data_ementoKey = pd.read_csv(file_ementoKey, delimiter=",", dtype='O', decimal='.')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Creating log with info on loaded data')
    todaystring = datetime.datetime.strftime(datetime.date.today(), "%y%m%d")
    if testing:
        fout = open(file_path+'/'+'Emento_datakombination_load_datafiles_log'+todaystring+'_testing.txt','w')
    else:
        fout = open(file_path + '/' + 'Emento_datakombination_load_datafiles_log' + todaystring + '.txt', 'w')
    fout.write('Log over datafiler der blev brugt i Emento_datakombination.load_datafiles() da kombineret fil '+todaystring+' blev genereret.\n')
    fout.write('\nFølgende filer blev loadet og kombineret:\n')
    fout.write('    SP data af alle kontakter for CPR numre i file_base                    : ' + file_contacts + ' \n')
    fout.write('    SP data med basis ud fra dia og proc (file_base)                       : ' + file_base + ' \n')
    fout.write('    Emento Tableau datafile fra https://storage-bi.emento.dk bucket        : ' + file_emento + ' \n')
    fout.write('    Emento UniqueIdentifier-CPR nøgle fra https://sla-mtk.emento.dk/ server: ' + file_ementoKey + ' \n')
    fout.close()

    if verbose: print(' - Done loading SP data base on diagnoses and procedures (data_base)')
    return data_contacts, data_base, data_emento, data_ementoKey

# -----------------------------------------------------------------------------------------------------------------------
def addCPR2Emento(data_emento, data_ementoKey, verbose=True):
    """
    Add CPR column to Emento data

    ---EXAMPLE OF USE---
    import Emento_SPdataberigelse as esd
    data_Emento_wCPR = esd.addCPR2Emento(data_emento, data_ementoKey, verbose=True)

    """
    if verbose: print(' - Matching CPR from data_ementoKey file to info from data_emento')
    data_emento_wCPR        = data_emento.copy()
    CPRcol = []

    Nnomatch = 0 # counter to keep track of the numer of Emento IDs without matches.
    for dd, Eid in enumerate(data_emento_wCPR['id']):
        ent_match = np.where(data_ementoKey['courseId'].values == Eid)[0]
        Nmatches = len(ent_match)

        if Nmatches == 0:
            #if verbose: print('   The following Emento ID did not have a match in data_ementoKey: '+str(Eid))
            Nnomatch = Nnomatch + 1
            CPRcol.append('')
        elif Nmatches == 1:
            cprstr_last4  = str(data_ementoKey['uniqueIdentifier'][ent_match[0]])[-4:]
            cprstr_first6 = str("%06.d" % int(str(data_ementoKey['uniqueIdentifier'][ent_match[0]])[:-4]))
            cprstr = cprstr_first6+'-'+cprstr_last4
            CPRcol.append(cprstr)
        else:
            if verbose: print('   The following Emento ID had --'+str(Nmatches)+'-- mastches in data_ementoKey: ' + str(Eid))

    data_emento_wCPR['uniqueIdentifier'] = CPRcol
    if verbose: print(' - There were '+str(Nnomatch)+' Emento IDs which did not have a match in the ID2CPR list')

    return data_emento_wCPR
# -----------------------------------------------------------------------------------------------------------------------
def add_emento_columns(data_SP,data_emento,Emento_column_key='_2403',verbose=True):
    """
    Function to add the columns from the Emento data to entries with dates
    after the first appearance of a given CPR in the Emento data.
    """
    SP_CPR = data_SP['Patient CPR-nr.']

    #SP_time = np.asarray([datetime.datetime.strptime(dstr, '%d-%m-%Y %H:%M') for dstr in data_SP['Behandlingskontakt starttidspunkt Dato-tid']])
    #SP_time = pd.to_datetime(data_SP['Behandlingskontakt starttidspunkt Dato-tid'])
    SP_time = pd.to_datetime(data_SP['Kontakt startdato Dato-tid'], format='%d-%m-%Y %H:%M')
    #Em_time = np.asarray([datetime.datetime.strptime(dstr[:-5], '%Y-%m-%dT%H:%M:%S') for dstr in data_emento['createdAt']])
    Em_time_cph = pd.to_datetime(data_emento['createdAt'],utc=True).map(lambda x: x.tz_convert('Europe/Copenhagen'))
    Em_time = pd.to_datetime(Em_time_cph.map(lambda x: x.tz_localize(None)))

    NlinesSP = len(SP_CPR)
    NuniqueSPCPR = len(np.unique(SP_CPR))

    multiEmentoCPRlist = [] # list to contain CPRs for cases with more than 1 Emento ID associated with them

    if verbose: print(' - defining output columns to append SP data')
    data_out_Ektkmatch = np.empty(NlinesSP)
    data_out_Ektkmatch[:] = np.nan

    data_out_Ecompscore = np.empty(NlinesSP)
    data_out_Ecompscore[:] = np.nan

    data_out_Eownscore = np.empty(NlinesSP)
    data_out_Eownscore[:] = np.nan

    data_out_Eid            = np.asarray([str('%100s' % '')]*NlinesSP)
    data_out_EcreatedAt     = np.asarray([str('%100s' % '')]*NlinesSP)
    data_out_Etitle         = np.asarray([str('%100s' % '')]*NlinesSP)
    data_out_EcontextTitle  = np.asarray([str('%100s' % '')]*NlinesSP)
    data_out_Eorganization  = np.asarray([str('%100s' % '')]*NlinesSP)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for ss, CPR_SP in enumerate(np.unique(SP_CPR)):
        if verbose:
            infostr = '   checking SP CPR number ' + str(ss + 1) + ' / ' + str(NuniqueSPCPR)
            sysstdout.write("%s\r" % infostr)
            sysstdout.flush()

        #finding entries for given CPR in data sets
        ent_SP  = np.where(SP_CPR == CPR_SP)[0]
        ent_E_all   = np.where(data_emento['uniqueIdentifier'] == CPR_SP)[0]

        Nmatches_E = len(ent_E_all)
        if Nmatches_E > 0:
            if Nmatches_E > 1: # handling case with multiple Emento matches
                multiEmentoCPRlist.append(CPR_SP) # count CPRs with more than one match

            for ee, ent_E in enumerate(ent_E_all):
                sel_time_SP = SP_time[ent_SP]
                sel_time_E = Em_time[ent_E]
                dt_arr = (sel_time_SP - sel_time_E).astype('timedelta64[m]') # timedifference in minutes
                ent_FirstAfterEmento = (dt_arr == np.min(dt_arr[dt_arr >= 0]))

                data_out_Ektkmatch[ent_SP[sel_time_SP > sel_time_E]] =+ 0
                data_out_Ektkmatch[ent_SP[ent_FirstAfterEmento]] =+ 1

                data_out_Ecompscore[ent_SP[sel_time_SP > sel_time_E]] = data_emento['complianceScore'][ent_E]
                data_out_Ecompscore[ent_SP[sel_time_SP > sel_time_E]] = data_emento['complianceScore'][ent_E]

                data_out_Eid[ent_SP[sel_time_SP > sel_time_E]] = data_emento['id'][ent_E]
                data_out_EcreatedAt[ent_SP[sel_time_SP > sel_time_E]] = data_emento['createdAt'][ent_E]
                data_out_Etitle[ent_SP[sel_time_SP > sel_time_E]] = data_emento['title'][ent_E]
                data_out_EcontextTitle[ent_SP[sel_time_SP > sel_time_E]] = data_emento['contextTitle'][ent_E]
                data_out_Eorganization[ent_SP[sel_time_SP > sel_time_E]] = data_emento['organization'][ent_E]

                inputlist = [data_emento['numActivitiesSeen'][ent_E],data_emento['numActivities'][ent_E]]
                val_ownscore = esd.cal_Ementoscore(inputlist, scoredef=1,verbose=False)
                data_out_Eownscore[ent_SP[sel_time_SP > sel_time_E]] = val_ownscore
            ''
            #if CPR_SP == '190660-2845':
            #if CPR_SP == '100101-7450':
            #    pdb.set_trace()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    NmultiEmentoCPRlist = len(multiEmentoCPRlist)
    if verbose: print('\n done ... found that '+str(NmultiEmentoCPRlist)+' CPRs had more than 1 emento ID associated with them')

    if verbose: print(' - Attaching data to output dataframe ')
    data_SP['Em_EmentoKontaktBesteMatch' + Emento_column_key] = data_out_Ektkmatch
    data_SP['Em_complianceScore' + Emento_column_key] = data_out_Ecompscore
    data_SP['Em_id' + Emento_column_key] = data_out_Eid
    data_SP['Em_createdAt' + Emento_column_key] = data_out_EcreatedAt
    data_SP['Em_title' + Emento_column_key] = data_out_Etitle
    data_SP['Em_contextTitle' + Emento_column_key] =data_out_EcontextTitle
    data_SP['Em_organization' + Emento_column_key] =data_out_Eorganization
    data_SP['Em_ownscore' + Emento_column_key] = data_out_Eownscore

    return data_SP,multiEmentoCPRlist
# -----------------------------------------------------------------------------------------------------------------------
def cal_Ementoscore(inputlist,scoredef=1,verbose=True):
    """
    Function to calculate a compliance score, but define from other criteria than the official score from Emento

    "scoredef" decides what version to return
    """

    if scoredef == 1:
        if verbose: print(' - Calculating Emento score assmunig inputlist = [numActivisiesSeen,numActivitets] as :')
        if verbose: print(' - score = numActivisiesSeen / numActivitets')
        val_score = inputlist[0]/inputlist[1]
    elif scoredef == 2:
        val_score = inputlist[0]*inputlist[1]
    else:
        val_score = -99

    return val_score
#=======================================================================================================================

