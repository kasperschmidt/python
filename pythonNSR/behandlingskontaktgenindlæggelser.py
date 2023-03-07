#=======================================================================================================================
#from importlib import reload
import easygui
import pandas as pd
import datetime
from time import sleep
from os import path as ospath
from sys import stdout as sysstdout
import numpy as np
from collections import OrderedDict
pd.options.mode.chained_assignment = None # surpress 'SettingWithCopyError' warnings
import pdb
from gc import collect as gccollect
#=======================================================================================================================
from memory_profiler import profile as memprofile # use @memprofile placed before "def" to memory profile of a function
#=======================================================================================================================
#Switches til kontrol af kode
GUIinput   = True # Aktiver GUI som beder om at indlæse Excel fil?
#=======================================================================================================================
def tjek_for_PI(dia_udsk,dia_alle,afsnit_alle,maade_udsk):
    """
    Tjek indhold af patientkotakt for at bestæmme om der er tale om en primærindlæggelse (PI)
    """
    if (str(dia_alle.values[-1])[:2] != 'DF') & \
            all(['hospice' not in str(afsn).lower() for afsn in afsnit_alle.values]) & \
            all([str(dia)[:5] != 'DZ763' for dia in dia_alle.values]) & \
            (str(dia_udsk.values)[0][:3] != 'DO4') & \
            (str(maade_udsk.values[0]) != 'DØD') & \
            all([str(dia)[:2] != 'DC' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD00' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD01' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD02' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD03' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD04' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD05' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD06' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD07' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD08' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD09' for dia in dia_alle.values]):
        result = 1
    else:
        result = 0
    return result
#=======================================================================================================================
def tjek_for_GI(ptktktype_start,dia_alle,afsnit_alle):
    """
    Tjek indhold af patientkontakt for at bestæmme om der er tale om en genindlæggelse (GI)
    """
    if (ptktktype_start.lower() == 'akut') & \
            all([str(dia)[:5] != 'DZ763' for dia in dia_alle.values]) & \
            all([str(dia)[:2] != 'DF' for dia in [dia_alle.values[0],dia_alle.values[-1]]]) & \
            all(['hospice' not in str(afsn).lower() for afsn in afsnit_alle.values]) & \
            all([str(dia)[:2] != 'DC' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD00' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD01' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD02' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD03' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD04' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD05' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD06' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD07' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD08' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DD09' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DO80' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DO81' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DO82' for dia in dia_alle.values]) & \
            all([str(dia)[:4] != 'DO84' for dia in dia_alle.values]) & \
            (str(dia_alle.values[0])[:2] != 'DS'  ) & \
            (str(dia_alle.values[0])[:3] != 'DT1' ) & \
            (str(dia_alle.values[0])[:3] != 'DT2' ) & \
            (str(dia_alle.values[0])[:4] != 'DT30') & \
            (str(dia_alle.values[0])[:4] != 'DT31') & \
            (str(dia_alle.values[0])[:4] != 'DT32') & \
            (str(dia_alle.values[0])[:4] != 'DT33') & \
            (str(dia_alle.values[0])[:4] != 'DT34') & \
            (str(dia_alle.values[0])[:4] != 'DT35') & \
            (str(dia_alle.values[0])[:4] != 'DT51') & \
            (str(dia_alle.values[0])[:4] != 'DT52') & \
            (str(dia_alle.values[0])[:4] != 'DT53') & \
            (str(dia_alle.values[0])[:4] != 'DT54') & \
            (str(dia_alle.values[0])[:4] != 'DT55') & \
            (str(dia_alle.values[0])[:4] != 'DT56') & \
            (str(dia_alle.values[0])[:4] != 'DT57') & \
            (str(dia_alle.values[0])[:4] != 'DT58') & \
            (str(dia_alle.values[0])[:4] != 'DT59') & \
            (str(dia_alle.values[0])[:3] != 'DT6' ) & \
            (str(dia_alle.values[0])[:3] != 'DT7' ) & \
            (str(dia_alle.values[0])[:3] != 'DT9' ) & \
            (str(dia_alle.values[0])[:2] != 'DX'  ) & \
            (str(dia_alle.values[0])[:2] != 'DY'  ):
        result = 1
    else:
        result = 0
    return result
#=======================================================================================================================
#@memprofile
def beregn_bgi():
    nowstring   = datetime.datetime.strftime(datetime.datetime.now(),"%d-%m-%Y %H:%M:%S")
    todaystring = datetime.datetime.strftime(datetime.date.today(),"%y%m%d")
    print("\n\n - Program til indentificering af behandlingskontaktgenindlæggelser startet "+nowstring)
    if GUIinput:
        title   = "Personhenførbare data?"
        message = "Skal personhenførbare data inkluderes i outputtet?"
        choices = ["Ja tak", "Nej tak"]
        output  = easygui.ynbox(message, title, choices)
        if output: # Hvis der trykkes Ja
            inkluderpersonoplysninger = True # Tilføj CPR og ID info i output?
        else: # Hvis der trykkes Nej
            inkluderpersonoplysninger = False
    else:
        inkluderpersonoplysninger = True
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if GUIinput:
        title="Behandlingskontaktgenindlæggelsesberegner"
        msg = "Identificer Excel datafilen som beregningerne skal baseres på."
        choices = ["Vælg Excel datafil"]
        reply = easygui.buttonbox(msg, title , choices=choices)
        if reply == choices[0]:
            excelfile = easygui.fileopenbox()
        else:
            print(' - fejl i angivelse af Excel fil')
    else:
        excelfile = "O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/genindlæggelser og genbesøg/" \
                    "Behandlingskontaktgenindlæggelser/Behandlingskontaktgenindlæggelser_datatræk_220401-220601_minimalt_kortversion.xlsx"

    outpath = ospath.dirname(excelfile)
    print(" - Excel datafil angivet: \n   "+ excelfile)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if GUIinput:
        wintitle   = "Genindlæggelsestidsinterval efter primær indlæggelse"
        wintext    = "Angiv tidsintervallet efter primær indlæggelse, hvori der skal tælles genindlæggelser. Værdier skal angives i timer. \nNationale kriterier: tmin = 0 og tmaks = 720 (30 dage)."
        input_list = ["tmin [timer]", "tmax [timer]"]

        # creating a integer box
        output     = easygui.multenterbox(wintext, wintitle, input_list, values=['0','720'])
        GItimeMin  = float(output[0])    # mindste tid i timer efter primær indlæggelse en genindlæggelse kan registreres
        GItimeMaks = float(output[1])  # maksimale tid i timer efter primær indlæggelse en genindlæggelse kan registreres
    else:
        GItimeMin  = 0.0    # mindste tid i timer efter primær indlæggelse en genindlæggelse kan registreres
        GItimeMaks = 720.0  # maksimale tid i timer efter primær indlæggelse en genindlæggelse kan registreres
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # indlæs Excel data i pandas dataframe
    print('\n - Indlæser Excel datafil...')
    df_kontakter_in = pd.read_excel(excelfile)
    print('   Sikrer at sorteringen er CPR > BHK udskrivning > Kontakt start > PTK ID ')
    df_kontakter = df_kontakter_in.sort_values(by=["Patient CPR-nr.", "Behandlingskontakt udskrivningsdato Dato-tid", "Kontakt startdato Dato-tid", "Patientkontakt record ID"], ascending=[True, True, True, True]).copy()
    df_kontakter = df_kontakter.reset_index(drop=True) # Sikrer at indeks følger sortering
    print('   Indlæste dataframe med ' + str(len(df_kontakter)) + ' rækker')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if GUIinput:
        wintitle   = "Kolonner der skal tilføjes output"
        wintext    = "Angiv kolonner fra input datafil der skal tilføjes outputtet"
        col_list   = df_kontakter.columns

        # creating a multi choice box
        include_col = easygui.multchoicebox(wintext, wintitle, col_list, preselect=None)
    else:
        include_col = None
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # definer en række vektorer med indhold til senere brug
    uniqueCPR = np.unique(df_kontakter['Patient CPR-nr.'])
    uniqueBHKID = np.unique(df_kontakter['Behandlingskontakt record ID'])
    udskrivningBHK = df_kontakter['Behandlingskontakt udskrivningsdato Dato-tid']
    indlaeggelseBHK = df_kontakter['Behandlingskontakt indlæggelsesdato Dato-tid']
    kontaktStart = df_kontakter['Kontakt startdato Dato-tid']

    #len(np.where(df_output['Kontakt startdato Dato-tid'] > datetime.datetime.strptime("25-05-2022 00:00:00", "%d-%m-%Y %H:%M:%S"))[0])
    #len(np.where(df_kontakter['Behandlingskontakt udskrivningsdato Dato-tid'] > datetime.datetime.strptime("25-05-2022 00:00:00", "%d-%m-%Y %H:%M:%S"))[0])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # definer output dataframe
    outputfilename = 'BGIberegninger'+todaystring+'.xlsx'
    #outputfilename = excelfile.replace('.xls','_BGIberegninger'+todaystring+'.xls')
    outputdata = {}
    NaNlist  = [np.nan] * len(df_kontakter['Hændelsesansvarlig Afsnit navn'])
    zerolist = [0] * len(df_kontakter['Hændelsesansvarlig Afsnit navn'])
    if inkluderpersonoplysninger:
        outputdata['CPR'] = list(df_kontakter['Patient CPR-nr.'])
        outputdata['PTK ID'] = list(df_kontakter['Patientkontakt record ID'])
        outputdata['BHK ID'] = list(df_kontakter['Behandlingskontakt record ID'])
        outputdata['Indlæggelsesmåde'] = list(df_kontakter['Indlæggelsesmåde navn'])
    if include_col is not None:
        for cc, colname in enumerate(include_col):
            outputdata[colname] = list(df_kontakter[colname])
    outputdata['Aktionsdiagnose'] = list(df_kontakter['Aktionsdiagnosekode'])
    outputdata['Aktionsdiagnose tekst'] = list(df_kontakter['Aktionsdiagnose kodetekst'])
    outputdata['Kontakt start'] = list(kontaktStart)
    outputdata['Kontakt slut'] = list(df_kontakter['Kontakt slutdato Dato-tid'])
    outputdata['Behandlingskontakt indlæggelse'] = list(indlaeggelseBHK)
    outputdata['Behandlingskontakt udskrivning'] = list(udskrivningBHK)
    outputdata['Hændelsestype'] = list(df_kontakter['Hændelsestype navn'])
    outputdata['Hændelsesansvarlig overafdeling'] = list(df_kontakter['Hændelsesansvarlig Overafdeling navn'])
    outputdata['Hændelsesansvarligt afsnit'] = list(df_kontakter['Hændelsesansvarlig Afsnit navn'])

    # lav ny kolonne der erstatter "SJ SLAGELSE MODTAGELSE" med SLAAKI1 og SLAAKI2 navne fra kontakt afsnit
    outputdata['Ansvarligt afsnit'] = [np.nan] * len(df_kontakter['Hændelsesansvarlig Afsnit navn'])

    for hh, hafsn in enumerate(df_kontakter['Hændelsesansvarlig Afsnit navn']):
        if ('SJ SLAGELSE MOD' in hafsn):
            outputdata['Ansvarligt afsnit'][hh] = df_kontakter['Kontaktansvar Afsnit navn'][hh]
        else:
            outputdata['Ansvarligt afsnit'][hh] = hafsn

    outputdata['Ansvarligt afsnit'] = np.asarray(outputdata['Ansvarligt afsnit'])

    outputdata['Primærindlæggelse'] = zerolist

    # Defniner kort-navne og lave de tomme kolonner der skal fyldes i outputtabellen
    unique_overafd     = np.unique(df_kontakter['Hændelsesansvarlig Overafdeling navn'])
    unique_overafd_NSR = []
    shortnames_overafd = {}
    for oo, overafd in enumerate(unique_overafd):
        shortnames_overafd[overafd] = overafd.split(' - ')[0].split(', ')[-1]
        if ('SLA ' in overafd) or ('NAE ' in overafd):
            outputdata['Genindlæggelse fra overafd. '+shortnames_overafd[overafd]] = NaNlist
            outputdata['PI for GI fra ' + shortnames_overafd[overafd]] = NaNlist
            unique_overafd_NSR.append(overafd)

    unique_afsn     = np.unique(outputdata['Ansvarligt afsnit'])

    unique_afsn_NSR = []
    shortnames_afsn = {}
    for aa, afsn in enumerate(unique_afsn):
        shortnames_afsn[afsn] = afsn.split(', ')[0]
        if ('SJ SLA' in afsn) or ('SJ NAE' in afsn):
            outputdata['Genindlæggelse fra afsn. '+shortnames_afsn[afsn]] = NaNlist
            outputdata['PI for GI fra ' + shortnames_afsn[afsn]] = NaNlist
            outputdata['PIDIA for GI fra ' + shortnames_afsn[afsn]] = NaNlist
            unique_afsn_NSR.append(afsn)

    df_output   = pd.DataFrame(outputdata)
    #print('--->'+str(sys.getsizeof(df_output)))
    print(' - Justerer dtypes i pandas dataframe ')
    for kk, key in enumerate(df_output.keys()):
        if 'Genindlæg' in key:
            df_output[key] = df_output[key].astype('float16')
    #print('--->' + str(sys.getsizeof(df_output)))

    del outputdata # reset output data array to save on memory use.
    gccollect()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Løkke over behandlingskontakter for at bestæmme primærindlæggelser
    print('\n - Identificerer primærindlæggelser')
    for bb, bid in enumerate(uniqueBHKID):
        infostr = '   Evaluerer udskrivende afsnit af behandlingskontakt ' + str(bb + 1) + ' / ' + str(len(uniqueBHKID))
        sysstdout.write("%s\r" % infostr)
        sysstdout.flush()

        # Vælger sidste (udskrivende) patientkontakt i behandlingskontakt
        ent_ptk_alle = np.where((df_kontakter['Behandlingskontakt record ID'] == bid))[0]
        ent_ptk_udsk = np.where((df_kontakter['Behandlingskontakt record ID'] == bid) & (df_kontakter['Hændelsestype navn'] == 'UDSKRIVNING'))[0]

        if len(ent_ptk_udsk) > 1: # Hvis der er mere end en udskrivning per behandlingskontakt bruges den seneste
            ent_seneste = np.where(df_kontakter['Kontakt slutdato Dato-tid'][ent_ptk_udsk] == np.max(df_kontakter['Kontakt slutdato Dato-tid'][ent_ptk_udsk]))[0]
            ent_ptk_udsk = ent_ptk_udsk[ent_seneste]

        # Kun tjek for primærindlæggelser for udskrivende overafdelinger på NSR
        if df_output['Ansvarligt afsnit'][ent_ptk_udsk].values in unique_afsn_NSR:
            isPI = tjek_for_PI(df_kontakter['Aktionsdiagnosekode'][ent_ptk_udsk],
                               df_kontakter['Aktionsdiagnosekode'][ent_ptk_alle],
                               df_kontakter['Hændelsesansvarlig Overafdeling navn'][ent_ptk_alle],
                               df_kontakter['Afslutningsmåde navn'][ent_ptk_udsk])
            df_output['Primærindlæggelse'][ent_ptk_udsk] = isPI

    print('\n - Færdig med evaluering af alle '+str(len(uniqueBHKID))+' behandlingskontakter')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' ´ Sletter arrays fra løkke')
    del ent_ptk_alle
    del ent_ptk_udsk
    del ent_seneste
    gccollect()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Løkke over patienter for at evaluere dem med henblik på genindlæggelser
    print('\n - Identificerer behandlingskontaktgenindlæggelser')
    for cc, cpr in enumerate(uniqueCPR):
        infostr = '   Evaluerer forløbene for patient ' + str(cc + 1) + ' / ' + str(len(uniqueCPR))
        sysstdout.write("%s\r" % infostr)
        sysstdout.flush()

        ent_cpr = np.where(df_kontakter['Patient CPR-nr.'] == cpr)[0]
        u_BHKID = np.unique(df_kontakter['Behandlingskontakt record ID'][ent_cpr])

        if len(u_BHKID) > 1:  # kun tjek genindlæggelser for patienter med mere end een behandlingskontakt
            ent_PI = np.where((df_kontakter['Patient CPR-nr.'] == cpr) & (df_output['Primærindlæggelse'] == 1))[0]
            #u_PIafsnit, ent_PIafsnit = np.unique(df_output['Ansvarligt afsnit'][ent_PI], return_index=True)
            u_PIafsnit = np.unique(df_output['Ansvarligt afsnit'][ent_PI])

            for aa, PIafs in enumerate(u_PIafsnit):
                #PIoverafd = df_kontakter['Hændelsesansvarlig Overafdeling navn'][ent_PI].values[ent_PIafsnit][0]

                for ee, ent_ktk in enumerate(ent_cpr): # løkke over alle patientkontakter for givent CPR nummer
                    # Find afsnit for alle kontakter markeret som primærindlæggelser for givent CPR nummer
                    ent_forPIogAfsnit = np.where((df_kontakter['Patient CPR-nr.'] == cpr) & (df_output['Primærindlæggelse'] == 1) & (df_output['Ansvarligt afsnit'] == PIafs))[0]

                    if len(ent_forPIogAfsnit) > 0: # Kun check for genindlæggelser hvis CPR har primærindlæggelse på afsnittet PIafs
                        PIoverafd         = np.unique(df_kontakter['Hændelsesansvarlig Overafdeling navn'][ent_forPIogAfsnit])  # Overafdeling for PI kontakter
                        PIdia             = np.asarray(df_kontakter['Aktionsdiagnosekode'][ent_forPIogAfsnit])  # Aktionsdiagnoser for PI kontakter

                        if len(PIoverafd) > 1: # Der burde kun være en unik overafdeling på dette tidspunkt, så advar hvis det ikke er tilfældet
                            print(' - ADVARSEL: For patientkontakt '+str(ent_ktk)+' var der '+str(len(PIoverafd))+' forskellige overafdelinger; bruger den første')
                        PIoverafd = PIoverafd[0]

                        BHKforPIogAfsnit  = df_kontakter['Behandlingskontakt record ID'][ent_forPIogAfsnit]
                        diffPItider       = (indlaeggelseBHK[ent_ktk] - udskrivningBHK[ent_forPIogAfsnit]) / np.timedelta64(1,'h')

                        if any((GItimeMin < diffPItider) & (diffPItider < GItimeMaks ) & (df_kontakter['Behandlingskontakt record ID'][ent_ktk] != BHKforPIogAfsnit)):
                            #ent_BHKID_alle = df_kontakter['Behandlingskontakt record ID'][ent_cpr]
                            ent_ptk_alle   = np.where(df_kontakter['Behandlingskontakt record ID'] == df_kontakter['Behandlingskontakt record ID'][ent_ktk])[0]

                            # df_kontakter['Aktionsdiagnosekode'][ent_ktk]
                            isGI = tjek_for_GI(df_kontakter['Indlæggelsesmåde navn'][ent_ptk_alle].values[0],
                                               df_kontakter['Aktionsdiagnosekode'][ent_ptk_alle],
                                               df_kontakter['Hændelsesansvarlig Overafdeling navn'][ent_ptk_alle])

                            df_output['Genindlæggelse fra afsn. ' + shortnames_afsn[PIafs]][ent_ktk] = isGI
                            df_output['Genindlæggelse fra overafd. ' + shortnames_overafd[PIoverafd]][ent_ktk] = isGI
                            if isGI == 1:
                                df_output['PI for GI fra ' + shortnames_afsn[PIafs]][ent_ktk] = PIafs
                                df_output['PI for GI fra ' + shortnames_overafd[PIoverafd]][ent_ktk] = PIoverafd

                                ent_senestePIdia = np.where((diffPItider[diffPItider > 0]) == np.min(diffPItider[diffPItider > 0])) # index for seneste PIs diagnose før geninglæggelsen
                                df_output['PIDIA for GI fra ' + shortnames_afsn[PIafs]][ent_ktk] = PIdia[ent_senestePIdia][0]

    print('\n - Færdig med evaluering af all '+str(len(uniqueCPR))+' patienter')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' ´ Sletter arrays fra løkke')
    del ent_PI
    del u_PIafsnit
    del PIoverafd
    del PIdia
    del ent_forPIogAfsnit
    del BHKforPIogAfsnit
    del diffPItider
    del ent_ptk_alle
    del ent_senestePIdia
    gccollect()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print('\n - Sikrer at der for hvert behandlingskontakt ID kun er 1 genindlæggelsesflag per overafdeling og afsnit')
    for bb, BHKID in enumerate(uniqueBHKID):
        infostr = '   Korrigerer genindlæggelsesflag for ID ' + str(bb + 1) + ' / ' + str(len(uniqueBHKID))
        sysstdout.write("%s\r" % infostr)
        sysstdout.flush()
        ent_BHKID_alle = np.where(df_kontakter['Behandlingskontakt record ID'] == BHKID)[0]

        for oo, overafd in enumerate(unique_overafd_NSR):
            if any(df_output['Genindlæggelse fra overafd. ' + shortnames_overafd[overafd]][ent_BHKID_alle] == 1):
                BHKI_overafd = np.unique(df_output['Hændelsesansvarlig overafdeling'][ent_BHKID_alle])
                for bo, bhk_overafd in enumerate(BHKI_overafd):
                    ent_BHK_overafd = np.where((df_output['Hændelsesansvarlig overafdeling'] == bhk_overafd) &
                                               (df_kontakter['Behandlingskontakt record ID'] == BHKID))[0]
                    selection = (df_output['Genindlæggelse fra overafd. ' + shortnames_overafd[overafd]][ent_BHK_overafd] != 1)

                    # Erstat værdier der ikke matcher 'selection' ovenfor; så find GI flag og erstat dem.
                    df_output['Genindlæggelse fra overafd. ' + shortnames_overafd[overafd]][ent_BHK_overafd] = \
                        df_output['Genindlæggelse fra overafd. ' + shortnames_overafd[overafd]][ent_BHK_overafd].where( selection , 0)  # replace values not matching selection
                    df_output['PI for GI fra ' + shortnames_overafd[overafd]][ent_BHK_overafd] = \
                        df_output['PI for GI fra ' + shortnames_overafd[overafd]][ent_BHK_overafd].where( selection, 'Udeladt') # replace values not matching selection

                    # Erstat første GI flag så det kun er registreret en gang per BHKID
                    df_output['Genindlæggelse fra overafd. ' + shortnames_overafd[overafd]][ent_BHK_overafd[0]] = 1
                    df_output['PI for GI fra ' + shortnames_overafd[overafd]][ent_BHK_overafd[0]] = overafd

        for aa, afsn in enumerate(unique_afsn_NSR):
            if any(df_output['Genindlæggelse fra afsn. ' + shortnames_afsn[afsn]][ent_BHKID_alle] == 1):
                BHKI_afsn = np.unique(df_output['Hændelsesansvarligt afsnit'][ent_BHKID_alle])
                for bo, bhk_afsn in enumerate(BHKI_afsn):
                    ent_BHK_afsn = np.where((df_output['Hændelsesansvarligt afsnit'] == bhk_afsn) &
                                               (df_kontakter['Behandlingskontakt record ID'] == BHKID))[0]
                    selection = (df_output['Genindlæggelse fra afsn. ' + shortnames_afsn[afsn]][ent_BHK_afsn] != 1)

                    # Erstat værdier der ikke matcher 'selection' ovenfor; så find GI flag og erstat dem.
                    df_output['Genindlæggelse fra afsn. ' + shortnames_afsn[afsn]][ent_BHK_afsn] = \
                        df_output['Genindlæggelse fra afsn. ' + shortnames_afsn[afsn]][ent_BHK_afsn].where( selection , 0)  # replace values not matching selection
                    df_output['PI for GI fra ' + shortnames_afsn[afsn]][ent_BHK_afsn] = \
                        df_output['PI for GI fra ' + shortnames_afsn[afsn]][ent_BHK_afsn].where( selection, 'Udeladt') # replace values not matching selection

                    # Erstat første GI flag så det kun er registreret en gang per BHKID
                    df_output['Genindlæggelse fra afsn. ' + shortnames_afsn[afsn]][ent_BHK_afsn[0]] = 1
                    df_output['PI for GI fra ' + shortnames_afsn[afsn]][ent_BHK_afsn[0]] = afsn

    print('\n - Færdig med korrigering af flag')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' ´ Sletter arrays fra løkke')
    del ent_BHKID_alle
    del BHKI_overafd
    del ent_BHK_overafd
    del selection
    del BHKI_afsn
    del ent_BHK_afsn
    gccollect()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print('\n - Tilføjer kolonner med samlet sum af genindlæggelser')
    df_output['Genindlæggelse afsnit total (sum)'] = df_output.loc[:, [cc for cc in df_output.columns if cc.startswith('Genindlæggelse fra afsn.')]].sum(axis=1)
    ent_lt0 = (df_output['Genindlæggelse afsnit total (sum)'] < 1)
    df_output['Genindlæggelse afsnit optælling'] = df_output['Genindlæggelse afsnit total (sum)'].where( ent_lt0 , 1)  # replace values not matching ent_lt0 with 1

    df_output['Genindlæggelse overafdeling total (sum)'] = df_output.loc[:, [cc for cc in df_output.columns if cc.startswith('Genindlæggelse fra overafd.')]].sum(axis=1)
    ent_lt0 = (df_output['Genindlæggelse overafdeling total (sum)'] < 1)
    df_output['Genindlæggelse overafdeling optælling'] = df_output['Genindlæggelse overafdeling total (sum)'].where( ent_lt0 , 1)  # replace values not matching ent_lt0 with 1

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print('\n - Opretter kolonne der markerer behandlingskontakt IDer med genindlæggelse')
    df_output['Behandlingskontakter med mindst en genindindlæggelse'] = zerolist
    for bb, BHKID in enumerate(uniqueBHKID):
        infostr = '   Checker ID ' + str(bb + 1) + ' / ' + str(len(uniqueBHKID))
        sysstdout.write("%s\r" % infostr)
        sysstdout.flush()

        ent_BHKID_alle = np.where(df_kontakter['Behandlingskontakt record ID'] == BHKID)[0]

        ent_BHKID_udskrivning = np.where((df_kontakter['Behandlingskontakt record ID'] == BHKID) & (df_kontakter['Hændelsestype navn'] == 'UDSKRIVNING'))[0]

        if any(df_output['Genindlæggelse overafdeling optælling'][ent_BHKID_alle] == 1):
            df_output['Behandlingskontakter med mindst en genindindlæggelse'][ent_BHKID_udskrivning[-1]] = 1

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print('\n - Gemmer output dataframe til Excelfil')
    if inkluderpersonoplysninger:
        print('   Sikrer at sorteringen er CPR > BHK udskrivning > Kontakt start')
        df_output = df_output.sort_values(by=['CPR','Behandlingskontakt udskrivning', 'Kontakt start'], ascending=[True, True, True])
    else:
        print('   Sikrer at sorteringen er BHK udskrivning > Kontakt start')
        df_output = df_output.sort_values(by=['Behandlingskontakt udskrivning','Kontakt start'], ascending=[True, True])

    df_output.to_excel(outpath+'/'+outputfilename, sheet_name="data output")

    print('\n - Output gemt i mappen '+outpath+'/')
    print('   med filnavnet '+outputfilename+'\n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    loggingfile = outpath+'/'+outputfilename.replace('.xlsx','_log.txt')
    print(" - Skriver stats til filen "+loggingfile)
    fout = open(loggingfile, 'w')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fout.write("\n - Program til indentificering af behandlingskontaktgenindlæggelser startet "+nowstring)
    fout.write('\n\n - Output gemt i mappen '+outpath+'/')
    fout.write('\n   med filnavnet '+outputfilename+'\n')

    inputfile = ospath.basename(excelfile)
    fout.write('\n - Indholdet i det beregnede output er baseret på filen '+inputfile)
    fout.write('\n   som ligger i mappen '+outpath+'/\n')
    fout.write('\n - Den gemte output ful indeholder indeholdere: ')
    fout.write('\n    o Antal primærindlæggelser (total)                                 = '+str("%12.2f" % df_output['Primærindlæggelse'].sum())+ ' (ud af '+str(len(uniqueBHKID))+' BHKIDer)')
    fout.write('\n    o Antal behandlingskontakter med mindst en genindlæggelse (total)  = '+str("%12.2f" % df_output['Behandlingskontakter med mindst en genindindlæggelse'].sum())+ ' (ud af '+str(len(uniqueBHKID))+' BHKIDer)')
    fout.write('\n    o Behandlingskontagenindlæggelser på afsnitsniveau (total)         = '+str("%12.2f" % df_output['Genindlæggelse afsnit optælling'].sum()))
    fout.write('\n    o Behandlingskontagenindlæggelser på overafdelingsniveu (total)    = '+str("%12.2f" % df_output['Genindlæggelse overafdeling optælling'].sum())+'\n')

    for oo, overafd in enumerate(unique_overafd_NSR):
        fout.write(('\n    o Behandlingskontagenindlæggelser fra overafd. ' + shortnames_overafd[overafd]).ljust(100)+' = ' +
              str("%12.2f" % df_output['Genindlæggelse fra overafd. ' + shortnames_overafd[overafd]].sum()))

    for aa, afsn in enumerate(unique_afsn_NSR):
        fout.write(('\n    o Behandlingskontagenindlæggelser fra afsn. ' + shortnames_afsn[afsn]).ljust(100)+' = ' +
              str("%12.2f" % df_output['Genindlæggelse fra afsn. ' + shortnames_afsn[afsn]].sum()))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nowstring   = datetime.datetime.strftime(datetime.datetime.now(),"%d-%m-%Y %H:%M:%S")
    todaystring = datetime.datetime.strftime(datetime.date.today(),"%y%m%d")
    print("\n\n - Program til indentificering af behandlingskontaktgenindlæggelser sluttede "+nowstring)
    fout.write("\n\n - Program til indentificering af behandlingskontaktgenindlæggelser sluttede "+nowstring)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fout.close() # closing log file
    print("\n - Log er gemt. Afslutter om ", end='')
    for ii in range(0,10):
        sleep(0.25)
        print(str(10-ii)+' ', end='')
    print(' --->')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#=======================================================================================================================
def append_new_bgi_calc(dateminadd = '31-12-2022', GUIinput=True):
    """
    Function appending a new set of BGI calculations (matching column names) to an existing BGI output.
    Data input is used to specify the earliest data (included) in the addition, that should be included.

    --- INPUT ---
    dateminadd          The minimum date to include in appended data (usually it would be the max date of mainBGIfile)
                        Use the format dd-mm-yyyy for the input string.
    GUIinput            To provide files interactively via GUIs set this to True; otherwise hardocede test files will
                        be used

    --- EXAMPLE OF USE ---
    import behandlingskontaktgenindlæggelser as bgi
    bgi.append_new_bgi_calc(dateminadd = '01-02-2023', GUIinput=True)

    bgi.append_new_bgi_calc(dateminadd = '15-01-2022', GUIinput=False) # for testing

    """
    nowstring   = datetime.datetime.strftime(datetime.datetime.now(),"%d-%m-%Y %H:%M:%S")
    todaystring = datetime.datetime.strftime(datetime.date.today(),"%y%m%d")
    print("\n - Tilføjer BGI beregninger til eksisterende BGI output; startede "+nowstring)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if GUIinput:
        title="BGImain Excel fil"
        msg = "Vælg Excel BGI datafil der skal tilføjes nye rækker til."
        choices = ["Vælg Excel datafil"]
        reply = easygui.buttonbox(msg, title , choices=choices)
        if reply == choices[0]:
            excelfile_main = easygui.fileopenbox()
        else:
            print(' - fejl i angivelse af Excel fil')
        excelfile_main = excelfile_main.replace('\\', '/')
    else:
        excelfile_main = "C:/Users/kaschm/Desktop/NSRLISmar2023/BGI230301/NSR LIS Behandlingskontaktgenindlæggelser med personoplysninger 230206 - kort.xlsx"

    print(" - BGI datafil angivet: \n   "+ excelfile_main)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if GUIinput:
        title="BGIadd Excel fil"
        msg = "Vælg Excel BGI datafil der skal tilføjes data fra."
        choices = ["Vælg Excel datafil"]
        reply = easygui.buttonbox(msg, title , choices=choices)
        if reply == choices[0]:
            excelfile_add = easygui.fileopenbox()
        else:
            print(' - fejl i angivelse af Excel fil')
        excelfile_add = excelfile_add.replace('\\','/')
    else:
        excelfile_add = "C:/Users/kaschm/Desktop/NSRLISmar2023/BGI230301/BGIberegninger230301kort.xlsx"

    print(" - BGI datafil der skal tilføjes data fra: \n   "+ excelfile_add)

    outpath    = ospath.dirname(excelfile_add)
    outputfile = outpath+'/'+excelfile_main.split('/')[-1].replace('.xls','_dataadd'+todaystring+'.xls')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datemin     = datetime.datetime.strptime(dateminadd, "%d-%m-%Y")
    print("   Tilføjer data fra denne fil der ligger efter "+ str(datetime.datetime.strftime(datemin,"%d-%m-%Y %H:%M:%S")))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # indlæs Excel data i pandas dataframe
    print('\n - Indlæser Excel datafiler...')
    print(' - BGImain: fil der skal tilføjes data til')
    df_BGImain = pd.read_excel(excelfile_main)
    df_BGImain = df_BGImain.drop('Unnamed: 0', axis='columns')
    print('   Indlæste dataframe med ' + str(len(df_BGImain)) + ' rækker')

    print(' - BGIadd:  fil der skal tilføjes data fra')
    df_BGIadd  = pd.read_excel(excelfile_add)
    df_BGIadd  = df_BGIadd.drop('Unnamed: 0', axis='columns')
    print('   Indlæste dataframe med ' + str(len(df_BGIadd)) + ' rækker')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datecol_main = df_BGImain["Behandlingskontakt udskrivning"]
    datecol_add  = df_BGIadd["Behandlingskontakt udskrivning"]

    dropval      = np.where(datecol_add < datemin)[0]
    if len(dropval) == 0:
        print('\n - ADVARSEL: 0 rækker fra BGI fil, der ikke skal tilføjes baseret på dato input')
    else:
        print('\n - Fjerner '+str(len(dropval))+' rækker fra BGI fil, der ikke skal tilføjes baseret op dato input')
        df_BGIadd = df_BGIadd.drop(df_BGIadd.index[dropval], axis='index') # removing the rows dropindex
    print('   Dette resulterer i dataframe med ' + str(len(df_BGIadd)) + ' rækker')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    print('\n   Gennemgår kolonner i BGImain og bygger dataframe der skal tilføjes ud fra BGIadd ')
    colnames_main = df_BGImain.columns
    colnames_add  = df_BGIadd.columns

    outputdata = OrderedDict()
    NaNlist = [np.nan] * len(df_BGIadd['CPR'])

    for cc, colname in enumerate(colnames_main):
        print('   Tilføjer data til kolonne '+str(cc+1)+'/'+str(len(colnames_main))+': '+colname)
        if colname in colnames_add:

            outputdata[colname] = list(df_BGIadd[colname])
        else:
            print('   ADVARSEL: kolonne ikke til stede i BGIadd; tilføjer NaNs')
            outputdata[colname] = NaNlist

    df_append   = pd.DataFrame(outputdata)

    print('\n   Tilføjer data ti BGImain ')
    df_BGImain  = df_BGImain.append(df_append)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print('\n - Gemmer output dataframe til Excelfil')
    print(' - Outputtet vil indeholde ' + str(len(df_BGImain)) + ' rækker')
    print('   Sikrer at sorteringen er CPR > BHK udskrivning > Kontakt start')
    df_BGImain = df_BGImain.sort_values(by=['CPR','Behandlingskontakt udskrivning', 'Kontakt start'], ascending=[True, True, True])
    df_BGImain.to_excel(outputfile, sheet_name="data output")
    print(' - Output gemt i filen '+outputfile)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nowstring   = datetime.datetime.strftime(datetime.datetime.now(),"%d-%m-%Y %H:%M:%S")
    print("\n - Tilføjelse af BGI beregninger til eksisterende BGI output sluttede " + nowstring)
#=======================================================================================================================
if __name__ == '__main__':
    beregn_bgi()
#=======================================================================================================================