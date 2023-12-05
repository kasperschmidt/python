#=======================================================================================================================
#from importlib import reload
import pandas as pd
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
    fileBestOrd    = pathKMT+"NSR KMT behandlingsniveau BestOrd stillingtagen.xlsx"
    fileAdmissions = pathKMT+"NSR KMT antal indlæggelser via hændelse.xlsx"
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
    pathKMT = "O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2023 - Kvalitetsmonitoreringstavle KMT/"

    df_bo, df_ad = bim.load_data(pathKMT,verbose=verbose)

    outputfilename = 'NSR KMT Stillingtagen til behandlingsniveau med personoplysninger.xlsx'
    #----------------------------------------------------------------------------
    if verbose: print(' - Forbereder output struktur og array')
    outputdata = {}
    NaNlist  = [np.nan] * len(df_ad['Hændelsesansvarlig Afsnit navn'])
    Strlist = [''] * len(df_ad['Hændelsesansvarlig Afsnit navn'])
    zerolist = [0] * len(df_ad['Hændelsesansvarlig Afsnit navn'])
    outputdata['CPR'] = list(df_ad['Patient CPR-nr.'])
    outputdata['CPR BestOrd'] = Strlist
    outputdata['Afsnit indlæggende'] = list(df_ad['Patientkontakt på hændelsestidspunkt Kontaktansvarlig afsnit'])
    for aa, hafs in enumerate(df_ad['Hændelsesansvarlig Afsnit navn']):
        if hafs.startswith('TRANSIT'):
            outputdata['Afsnit indlæggende'][aa] = hafs

    outputdata['Afsnit BestOrd'] = Strlist
    outputdata['Dato-tid indlæggelse'] = list(df_ad['Hændelsestidspunkt Dato-tid'])
    outputdata['Dato-tid BestOrd'] = Strlist
    outputdata['Dato-tid udskrivning'] = list(df_ad['Behandlingskontakt udskrivningsdato Dato-tid'])
    outputdata['Tidsforskel [timer]'] = NaNlist
    outputdata['Hændelsestype'] = list(df_ad['Hændelsestype navn'])
    outputdata['Patientkontakttype'] = list(df_ad['Patientkontakttype navn'])
    outputdata['Patientalder'] = list(df_ad['Patient alder ved Behandlingskontaktens start'])
    outputdata['Patientalder gruppering'] = Strlist.copy()
    for pa, alder in enumerate(df_ad['Patient alder ved Behandlingskontaktens start']):
        if (alder < 18):
            outputdata['Patientalder gruppering'][pa] = '0-17'
        elif (alder > 17) & (alder < 30):
            outputdata['Patientalder gruppering'][pa] = '18-29'
        elif (alder > 29) & (alder < 50):
            outputdata['Patientalder gruppering'][pa] = '30-49'
        elif (alder > 49) & (alder < 70):
            outputdata['Patientalder gruppering'][pa] = '40-69'
        elif (alder > 69) & (alder < 90):
            outputdata['Patientalder gruppering'][pa] = '70-89'
        else:
            outputdata['Patientalder gruppering'][pa] = '90+'

    df_output = pd.DataFrame(outputdata)
    # ----------------------------------------------------------------------------
    if verbose: print(' - Gennemgår BestOrd og matcher til indlæggelsesliste '
                      '(match på CPR og best.ord. tid mellem indlæggelse og udskrivning)')
    Nnomatch = 0
    Nmissingadmission = 0
    for bb, afs_bestord in enumerate(df_bo['Afsnit']):
        infostr = '   Matcher BestOrd ' + str(bb + 1) + ' / ' + str(len(df_bo['Afsnit']))
        sysstdout.write("%s\r" % infostr)
        sysstdout.flush()

        # søg efter best.ord der er udført for samme CPR mellem indlæggelse og udskrivning

        ment = np.where((df_bo['CPR'][bb].strip() == df_output['CPR']) &
                        (df_bo['BestOrd datotid'][bb] > df_output['Dato-tid indlæggelse']) &
                        (df_bo['BestOrd datotid'][bb] < df_output['Dato-tid udskrivning']) )[0]
        if len(ment) == 1:
            dtime = df_bo['BestOrd datotid'][bb] - df_output['Dato-tid indlæggelse'][ment]
            dtime_hours = dtime.dt.total_seconds().values[0] / 3600.0

            if dtime_hours < 0:
                print("   INFO: Indlæggelse af " + str(
                    df_bo["CPR"][bb]) + " mangler for BestOrd da dt<0:\n             Indlæggelse: " + str(
                    df_output['Dato-tid indlæggelse'][ment].values) + ' og BestOrd: ' + str(
                    df_bo['BestOrd datotid'][bb]))
                Nmissingadmission = Nmissingadmission + 1
            else:
                df_output['Tidsforskel [timer]'][ment[0]]    = dtime_hours
                df_output['Afsnit BestOrd'][ment[0]] = afs_bestord
                df_output['Dato-tid BestOrd'][ment[0]] = df_bo['BestOrd datotid'][bb]
                df_output['CPR BestOrd'][ment[0]] = df_bo['CPR'][bb].strip()

        elif len(ment) > 1:
            dtimes = df_bo['BestOrd datotid'][bb] - df_output['Dato-tid indlæggelse'][ment]
            dtimes_hours = dtimes.dt.total_seconds().values / 3600.0
            if (dtimes_hours > 0).any():
                dtime_hours = np.min(dtimes_hours[dtimes_hours>0])
                tent = np.where(dtimes_hours == dtime_hours)[0]

                df_output['Tidsforskel [timer]'][ment[tent[0]]] = dtime_hours
                df_output['Afsnit BestOrd'][ment[tent[0]]] = afs_bestord
                df_output['Dato-tid BestOrd'][ment[tent[0]]] = df_bo['BestOrd datotid'][bb]
                df_output['CPR BestOrd'][ment[tent[0]]] = df_bo['CPR'][bb].strip()
            else:
                print("   INFO: Indlæggelse af " + str(
                    df_bo["CPR"][bb]) + " mangler for BestOrd da alle dt<0:\n             Indlæggelser: " + str(
                    df_output['Dato-tid indlæggelse'][ment].values) + ' og BestOrd: ' + str(
                    df_bo['BestOrd datotid'][bb]))
                Nmissingadmission = Nmissingadmission + 1
        else:
            print("   INFO: BestOrd fra " + str(df_bo['BestOrd datotid'][bb]) + " for " + str(
                df_bo["CPR"][bb]) + " på " + afs_bestord + " havde intet match jf. kriterier på indlæggelsesfil")
            Nnomatch = Nnomatch + 1

    if verbose: print('\n   ... Færdig med match')

    if verbose: print('\n - Fandt ' + str(Nnomatch) + ' BestOrd uden CPR match i indlæggelsesfil')
    if verbose: print('\n - Fandt ' + str(Nmissingadmission) + ' BestOrd der lå før indlæggelser i indlæggelsesfil')

    # ----------------------------------------------------------------------------
    if verbose: print(' - Genererer output')
    df_output.to_excel(pathKMT+outputfilename, sheet_name="behandlingsniveau")

    print('\n - Output gemt i mappen '+pathKMT+'/')
    print('   med filnavnet "'+outputfilename+'"\n')

#=======================================================================================================================

