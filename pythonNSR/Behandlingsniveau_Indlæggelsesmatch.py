#=======================================================================================================================
#from importlib import reload
import pandas as pd
import datetime
import pdb
from sys import stdout as sysstdout
import numpy as np
pd.options.mode.chained_assignment = None # surpress 'SettingWithCopyError' warnings
import Behandlingsniveau_Indlæggelsesmatch as bim
#=======================================================================================================================
def load_data(pathKMT,verbose=True):
    """

    -- EXAMPLE OF USE --
    import Behandlingsniveau_Indlæggelsesmatch as bim
    pathKMT = "O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2023 - Kvalitetsmonitoreringstavel KMT/"
    df_bo, df_ad = bim.load_data(pathKMT)
    """
    fileBestOrd    = pathKMT+"Behandlingsniveau_BestOrd_stillingtagen.xlsx"
    fileAdmissions = pathKMT+"Antal_indlæggelser_via_hændelse.xlsx"

    print('\n - Indlæser BestOrd Excel datafil \n   ('+fileBestOrd+')')
    df_bo = pd.read_excel(fileBestOrd)

    print('\n - Indlæser Indlæggelses Excel datafil \n   ('+fileAdmissions+')')
    df_ad = pd.read_excel(fileAdmissions)

    print('\n - Færdig med at indlæse filer')
    return df_bo, df_ad
#=======================================================================================================================
def generate_output(verbose=True):
    """
    function to combine BestOrd and Admission data, tjecken for best match to BestOrd timestamps.

    -- EXAMPLE OF USE --
    import Behandlingsniveau_Indlæggelsesmatch as bim
    bim.generate_output(verbose=True)

    """
    pathKMT = "O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2023 - Kvalitetsmonitoreringstavel KMT/"

    df_bo, df_ad = bim.load_data(pathKMT,verbose=verbose)

    outputfilename = 'StillingtagenTilBehandlingsniveau.xlsx'
    #----------------------------------------------------------------------------
    if verbose: print(' - Forbereder output struktur og array')
    outputdata = {}
    NaNlist  = [np.nan] * len(df_ad['Hændelsesansvarlig Afsnit navn'])
    Strlist = [''] * len(df_ad['Hændelsesansvarlig Afsnit navn'])
    zerolist = [0] * len(df_ad['Hændelsesansvarlig Afsnit navn'])
    outputdata['CPR'] = list(df_ad['Patient CPR-nr.'])
    outputdata['Afsnit BestOrd'] = Strlist
    outputdata['Afsnit indlæggende'] = list(df_ad['Patientkontakt på hændelsestidspunkt Kontaktansvarlig afsnit'])
    outputdata['Dato-tid BestOrd'] = Strlist
    outputdata['Dato-tid indlæggelse'] = list(df_ad['Hændelsestidspunkt Dato-tid'])
    outputdata['Tidsforskel [timer]'] = NaNlist
    outputdata['Patientalder'] = list(df_ad['Patient alder ved Behandlingskontaktens start'])

    df_output = pd.DataFrame(outputdata)
    # ----------------------------------------------------------------------------
    if verbose: print(' - Gennemgår BestOrd og matcher til indlæggelsesliste')
    for bb, efs_bestord in enumerate(df_bo['Afsnit']):
        infostr = '   Matcher BestOrd ' + str(bb + 1) + ' / ' + str(len(df_bo['Afsnit']))
        sysstdout.write("%s\r" % infostr)
        sysstdout.flush()

        ment = np.where(df_output['CPR'] == df_bo['CPR'][bb])[0]
        if len(ment) == 1:
            dt = df_output['Dato-tid indlæggelse'][ment]-df_bo['BestOrd datotid'][bb]
            if dt < 0:
                print(" WARNING: somethings fishy - admission for BestOrd missing as dt<0; "+str(df_output['Dato-tid indlæggelse'][ment])+';'+str(df_bo['BestOrd datotid'][bb]))
            else:
                outputdata['Tidsforskel'][ment] = dt.hours
                outputdata['Afsnit BestOrd'][ment] = df_bo['Afsnit'][bb]

        elif len(ment) > 1:
            dts = df_output['Dato-tid indlæggelse'][ment] - df_bo['BestOrd datotid'][bb]
            dt = np.min(dts[dts>0])
            tent = np.where(dts == dt)[0]

            outputdata['Tidsforskel'][ment[tent]] = dt.hours
            outputdata['Afsnit BestOrd'][ment[tent]] = df_bo['Afsnit'][bb]

    if verbose: print('\n   ... done matching')



    # ----------------------------------------------------------------------------
    if verbose: print(' - Genererer output')
    df_output.to_excel(pathKMT+outputfilename, sheet_name="behandlingsniveau")

    print('\n - Output gemt i mappen '+pathKMT+'/')
    print('   med filnavnet '+outputfilename+'\n')
#=======================================================================================================================

