#=======================================================================================================================
# Funktion der udvælger genindlæggelser der er resultat af primærindlæggelser med en given aktionsdiagnose og hvor
# udskrivning på genindlæggelsen er på enten R0 (Akutten) eller R8 (Lungemed)
#=======================================================================================================================
#from importlib import reload
import easygui
import pandas as pd
import datetime
from time import sleep
import pdb
import os
import sys
import numpy as np
import gc
import BGI_count_readmissions_based_on_PDIA as bcr
#=======================================================================================================================
def count_bgi(excelfile=None):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Example of use:
    # import BGI_count_readmissions_based_on_PDIA as bcr
    # bcr.count_bgi(excelfile="O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed KOL patienter/NSR LIS Behandlingskontaktgenindlæggelser med personoplysninger.xlsx")
    #
    # Here the excelfil contains the results from a BGI calculation based on scripts similar to behandlingskontaktgeindlæggelser.py
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if excelfile is None:
        excelfile = "O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed KOL patienter/NSR LIS Behandlingskontaktgenindlæggelser med personoplysninger - kort.xlsx"
        print('\n Advarsel: Ingen Excelfil med BGI beregning givet så input sat til en kort NSR LIS fil:\n   '+excelfile)
    # indlæs Excel data i pandas dataframe
    print('\n - Indlæser Excel datafil...')
    df_bgi_in = pd.read_excel(excelfile)
    print('   Indlæste dataframe med ' + str(len(df_bgi_in)) + ' rækker')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    colnames_GI = []
    for cc, colname in enumerate(df_bgi_in.columns):
        if str(colname).startswith('Genindlæggelse fra afsn.'):
            colnames_GI.append(colname)

    adialist = ['DJ44', 'DJ96', 'DJ13', 'DJ14', 'DJ15', 'DJ16', 'DJ17', 'DJ18']

    ADIA_on_PI  = [''] * len(df_bgi_in)
    wards_on_PI = [''] * len(df_bgi_in)
    for ii, ptkid in enumerate(df_bgi_in['PTK ID']):
        infostr = '   Tjekker patientkontakt ' + str(ii + 1) + ' / ' + str(len(df_bgi_in['PTK ID']))
        sys.stdout.write("%s\r" % infostr)
        sys.stdout.flush()

        # first check if the row satisfies the following conditions to be counted
        if (df_bgi_in['Behandlingskontakt udskrivning'][ii] > datetime.datetime.strptime("01-01-2020 00:00:00", "%d-%m-%Y %H:%M:%S")) and \
            (df_bgi_in['Behandlingskontakt udskrivning'][ii] < datetime.datetime.strptime("01-01-2025 00:00:00", "%d-%m-%Y %H:%M:%S")) and \
            (df_bgi_in['Hændelsestype'][ii] == 'UDSKRIVNING') and \
            ( (df_bgi_in['Hændelsesansvarlig overafdeling'][ii] == 'SLA R8, SLA AKUTAFDELING- OVERAFDELING') or \
              (df_bgi_in['Hændelsesansvarlig overafdeling'][ii] == 'SLA R0, SLA LUNGEMEDICIN - OVERAFDELING') ):

            PDIAs = []  # reset PDIAs list
            wards = []  # reset ward list
            # then check if the conditions for readmission of interest are satisfied for any of the readmission columns in the ii'th row
            for gg, GIcol in enumerate(colnames_GI):
                afsn   = GIcol.split('afsn. ')[-1]
                GIflag = False  # reset flag for readmission
                if  ~np.isnan(df_bgi_in[GIcol][ii]) and (str(df_bgi_in['PIDIA for GI fra ' + afsn][ii]) != 'nan'): # check that row is a readmission
                    if (df_bgi_in['PIDIA for GI fra ' + afsn][ii][:4] in adialist): # check that the diagnose on the primary admission is in desired list
                        #print(df_bgi_in['PIDIA for GI fra ' + afsn][ii]) # printing the diagnoses that fullfille the criteria
                        GIflag = True # mark that there was a readmission fulfilling criteria
                        PDIAs.append(df_bgi_in['PIDIA for GI fra '+afsn][ii])
                        wards.append(afsn)
            ADIA_on_PI[ii]  = '; '.join(PDIAs)
            wards_on_PI[ii] = '; '.join(wards)
    print('\n - Færdig med patientkontakttjek ')
    ADIA_on_PI = np.asarray(ADIA_on_PI)
    wards_on_PI = np.asarray(wards_on_PI)
    goodent = np.where(ADIA_on_PI != '')[0]
    print(' - Fandt '+str(len(goodent))+' genindlæggelser der opfylder afsnit og diagnose kriterierne')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' - Forbereder output')
    outputdata = {}
    outputdata['CPR'] = list(df_bgi_in['CPR'][goodent])
    outputdata['Genindlæggelse PTK ID'] = list(df_bgi_in['PTK ID'][goodent])
    outputdata['Genindlæggelse BHK ID'] = list(df_bgi_in['BHK ID'][goodent])
    outputdata['Genindlæggelse aktionsdiagnose'] = list(df_bgi_in['Aktionsdiagnose'][goodent])
    outputdata['Genindlæggelse udskrivende afsnit '] = list(df_bgi_in['Ansvarligt afsnit'][goodent])
    outputdata['Genindlæggelse BHK udskrivelsesdatotid'] = list(df_bgi_in['Behandlingskontakt udskrivning'][goodent])
    outputdata['Genindlæggelse PTK hændelsestype'] = list(df_bgi_in['Hændelsestype'][goodent])
    outputdata['Primærindlæggelse ADIA'] = list(ADIA_on_PI[goodent])
    outputdata['Primærindlæggelse afsnit'] = list(wards_on_PI[goodent])

    df_output = pd.DataFrame(outputdata)
    outfilename = excelfile.split('.xls')[0]+' BGIcountReadmissions.xlsx'
    df_output.to_excel(outfilename, sheet_name="data output")
    print('\n - Output gemt i filen "'+outfilename+'"')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#=======================================================================================================================