#=======================================================================================================================
#from importlib import reload
import easygui
import pandas as pd
import datetime
import pdb
# import os
import sys
import numpy as np
pd.options.mode.chained_assignment = None # surpress 'SettingWithCopyError' warnings
#=======================================================================================================================
#Switches til kontrol af kode
GUIinput = False # Aktiver GUI som beder om at indlæse Excel fil?
inkluderpersonoplysninger = True # Tilføj CPR og ID info i output?
#=======================================================================================================================
def tjek_for_PI(dia_udsk,dia_alle,afsnit_alle):
    """
    Tjek indhold af patientkotakt for at bestæmme om der er tale om en primærindlæggelse (PI)
    """
    if (str(dia_alle.values[-1])[:2] != 'DF') & \
            all(['hospice' not in str(afsn).lower() for afsn in afsnit_alle.values]) & \
            all([str(dia)[:5] != 'DZ763' for dia in dia_alle.values]) & \
            (str(dia_udsk.values)[0][:3] != 'DO4') & \
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
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
nowstring   = datetime.datetime.strftime(datetime.datetime.now(),"%d-%m-%Y %H:%M:%S")
todaystring = datetime.datetime.strftime(datetime.date.today(),"%y%m%d")
print("\n\n - Program til indentificering af behandlingskontaktgenindlæggelser startet "+nowstring)
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
                "Behandlingskontaktgenindlæggelser/Behandlingskontaktgenindlæggelser_datatræk_220401-220601_ingenBiDiag_kortversion.xlsx"

print(" - Excel datafil angivet: \n   "+ excelfile)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# indlæs Excel data i pandas dataframe
print('\n - Indlæser Excel datafil...')
df_kontakter_in = pd.read_excel(excelfile)
print('   Sikrer at sorteringen er CPR > BHK udskrivning > Kontakt start > PTK ID ')
df_kontakter = df_kontakter_in.sort_values(by=["Patient CPR-nr.", "Behandlingskontakt udskrivningsdato Dato-tid", "Kontakt startdato Dato-tid", "Patientkontakt record ID"], ascending=[True, True, True, True]).copy()
df_kontakter = df_kontakter.reset_index(drop=True) # Sikrer at indeks følger sortering
print('   Indlæste dataframe med ' + str(len(df_kontakter)) + ' rækker')
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
outputdata['Aktionsdiagnose'] = list(df_kontakter['Aktionsdiagnosekode'])
outputdata['Aktionsdiagnose tekst'] = list(df_kontakter['Aktionsdiagnose kodetekst'])
outputdata['Kontakt start'] = list(kontaktStart)
outputdata['Kontakt slut'] = list(df_kontakter['Kontakt slutdato Dato-tid'])
outputdata['Behandlingskontakt indlæggelse'] = list(indlaeggelseBHK)
outputdata['Behandlingskontakt udskrivning'] = list(udskrivningBHK)
outputdata['Hændelsestype'] = list(df_kontakter['Hændelsestype navn'])
outputdata['Hændelsesansvarlig overafdeling'] = list(df_kontakter['Hændelsesansvarlig Overafdeling navn'])
outputdata['Hændelsesansvarligt afsnit'] = list(df_kontakter['Hændelsesansvarlig Afsnit navn'])
outputdata['Primærindlæggelse'] = zerolist

unique_overafd     = np.unique(df_kontakter['Hændelsesansvarlig Overafdeling navn'])
unique_overafd_NSR = []
shortnames_overafd = {}
for oo, overafd in enumerate(unique_overafd):
    shortnames_overafd[overafd] = overafd.split(' - ')[0].split(', ')[-1]
    if ('SLA ' in overafd) or ('NAE ' in overafd):
        outputdata['Genindlæggelse fra '+shortnames_overafd[overafd]] = NaNlist
        outputdata['PI for GI fra ' + shortnames_overafd[overafd]] = NaNlist
        unique_overafd_NSR.append(overafd)

unique_afsn     = np.unique(df_kontakter['Hændelsesansvarlig Afsnit navn'])
unique_afsn_NSR = []
shortnames_afsn = {}
for aa, afsn in enumerate(unique_afsn):
    shortnames_afsn[afsn] = afsn.split(', ')[0]
    if ('SJ SLA' in afsn) or ('SJ NAE' in afsn):
    #if ('SJ SLAGELSE MOD' in afsn):
        outputdata['Genindlæggelse fra '+shortnames_afsn[afsn]] = NaNlist
        outputdata['PI for GI fra ' + shortnames_afsn[afsn]] = NaNlist
        unique_afsn_NSR.append(afsn)

df_output = pd.DataFrame(outputdata)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Løkke over behandlingskontakter for at bestæmme primærindlæggelser
print('\n - Identificerer primærindlæggelser')
for bb, bid in enumerate(uniqueBHKID):
    infostr = '   Evaluerer udskrivende afsnit af behandlingskontakt ' + str(bb + 1) + ' / ' + str(len(uniqueBHKID))
    sys.stdout.write("%s\r" % infostr)
    sys.stdout.flush()

    # Vælger sidste (udskrivende) patientkontakt i behandlingskontakt
    ent_ptk_alle = np.where((df_kontakter['Behandlingskontakt record ID'] == bid))[0]
    ent_ptk_udsk = np.where((df_kontakter['Behandlingskontakt record ID'] == bid) & (df_kontakter['Hændelsestype navn'] == 'UDSKRIVNING'))[0]

    if len(ent_ptk_udsk) > 1: # Hvis der er mere end en udskrivning per behandlingskontakt bruges den seneste
        ent_seneste = np.where(df_kontakter['Kontakt slutdato Dato-tid'][ent_ptk_udsk] == np.max(df_kontakter['Kontakt slutdato Dato-tid'][ent_ptk_udsk]))[0]
        ent_ptk_udsk = ent_ptk_udsk[ent_seneste]

    # Kun tjek for primærindlæggelser for udskrivende overafdelinger på NSR
    if df_kontakter['Hændelsesansvarlig Afsnit navn'][ent_ptk_udsk].values in unique_afsn_NSR:
        isPI = tjek_for_PI(df_kontakter['Aktionsdiagnosekode'][ent_ptk_udsk],
                           df_kontakter['Aktionsdiagnosekode'][ent_ptk_alle],
                           df_kontakter['Hændelsesansvarlig Overafdeling navn'][ent_ptk_alle])
        df_output['Primærindlæggelse'][ent_ptk_udsk] = isPI

print('\n - Færdig med evaluering af alle '+str(len(uniqueBHKID))+' behandlingskontakter')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Løkke over patienter for at evaluere dem med henblik på genindlæggelser
print('\n - Identificerer behandlingskontaktgenindlæggelser')
for cc, cpr in enumerate(uniqueCPR):
    infostr = '   Evaluerer forløbene for patient ' + str(cc + 1) + ' / ' + str(len(uniqueCPR))
    sys.stdout.write("%s\r" % infostr)
    sys.stdout.flush()

    ent_cpr = np.where(df_kontakter['Patient CPR-nr.'] == cpr)[0]
    ent_PI  = np.where((df_kontakter['Patient CPR-nr.'] == cpr) & (df_output['Primærindlæggelse'] == 1))[0]
    u_BHKID = np.unique(df_kontakter['Behandlingskontakt record ID'][ent_cpr])
    u_PIafsnit, ent_PIafsnit = np.unique(df_kontakter['Hændelsesansvarlig Afsnit navn'][ent_PI], return_index=True)

    if len(u_BHKID) > 1:  # kun tjek genindlæggelser for patienter med mere end een behandlingskontakt
        for aa, PIafs in enumerate(u_PIafsnit):
            PIoverafd = df_kontakter['Hændelsesansvarlig Overafdeling navn'][ent_PI].values[ent_PIafsnit][0]
            for ee, ent_ktk in enumerate(ent_cpr):
                ent_forPIogAfsnit = np.where((df_kontakter['Patient CPR-nr.'] == cpr) & (df_output['Primærindlæggelse'] == 1) & (df_kontakter['Hændelsesansvarlig Afsnit navn'] == PIafs))[0]
                if len(ent_forPIogAfsnit) > 0: # Kun check for genindlæggelser hvis CPR har primærindlæggelse på afsnittet PIafs
                    BHKforPIogAfsnit    = df_kontakter['Behandlingskontakt record ID'][ent_forPIogAfsnit]
                    diffPItider = (indlaeggelseBHK[ent_ktk] - udskrivningBHK[ent_forPIogAfsnit]) / np.timedelta64(1,'h')

                    if any((0 < diffPItider) & (diffPItider < 720 ) & (df_kontakter['Behandlingskontakt record ID'][ent_ktk] != BHKforPIogAfsnit)):
                        ent_BHKID_alle = df_kontakter['Behandlingskontakt record ID'][ent_cpr]
                        ent_ptk_alle   = np.where(df_kontakter['Behandlingskontakt record ID'] == df_kontakter['Behandlingskontakt record ID'][ent_ktk])[0]

                        # df_kontakter['Aktionsdiagnosekode'][ent_ktk]
                        isGI = tjek_for_GI(df_kontakter['Indlæggelsesmåde navn'][ent_ptk_alle].values[0],
                                           df_kontakter['Aktionsdiagnosekode'][ent_ptk_alle],
                                           df_kontakter['Hændelsesansvarlig Overafdeling navn'][ent_ptk_alle])

                        df_output['Genindlæggelse fra ' + shortnames_afsn[PIafs]][ent_ktk] = isGI
                        df_output['Genindlæggelse fra ' + shortnames_overafd[PIoverafd]][ent_ktk] = isGI
                        if isGI == 1:
                            df_output['PI for GI fra ' + shortnames_afsn[PIafs]][ent_ktk] = PIafs
                            df_output['PI for GI fra ' + shortnames_overafd[PIoverafd]][ent_ktk] = PIoverafd

print('\n - Færdig med evaluering af all '+str(len(uniqueCPR))+' patienter')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print('\n - Sikrer at der for hvert behandlingskontakt ID kun er 1 genindlæggelsesflag per overafdeling og afsnit')
#print('\n\n\n\nIKKE KODET\n\n\n\n\n\n')
for bb, BHKID in enumerate(uniqueBHKID):
    infostr = '   Korrigerer genindlæggelsesflag for ID ' + str(bb + 1) + ' / ' + str(len(uniqueBHKID))
    sys.stdout.write("%s\r" % infostr)
    sys.stdout.flush()
    ent_BHKID_alle = np.where(df_kontakter['Behandlingskontakt record ID'] == BHKID)[0]

    for oo, overafd in enumerate(unique_overafd_NSR):
        if any(df_output['Genindlæggelse fra ' + shortnames_overafd[overafd]][ent_BHKID_alle] == 1):
            BHKI_overafd = np.unique(df_output['Hændelsesansvarlig overafdeling'][ent_BHKID_alle])
            for bo, bhk_overafd in enumerate(BHKI_overafd):
                ent_BHK_overafd = np.where((df_output['Hændelsesansvarlig overafdeling'] == bhk_overafd) &
                                           (df_kontakter['Behandlingskontakt record ID'] == BHKID))[0]
                selection = (df_output['Genindlæggelse fra ' + shortnames_overafd[overafd]][ent_BHK_overafd] != 1)

                # Erstat værdier der ikke matcher 'selection' ovenfor; så find GI flag og erstat dem.
                df_output['Genindlæggelse fra ' + shortnames_overafd[overafd]][ent_BHK_overafd] = \
                    df_output['Genindlæggelse fra ' + shortnames_overafd[overafd]][ent_BHK_overafd].where( selection , 99)  # replace values not matching selection
                df_output['PI for GI fra ' + shortnames_overafd[overafd]][ent_BHK_overafd] = \
                    df_output['PI for GI fra ' + shortnames_overafd[overafd]][ent_BHK_overafd].where( selection, 'Udeladt') # replace values not matching selection

                # Erstat første GI flag så det kun er registreret en gang per BHKID
                df_output['Genindlæggelse fra ' + shortnames_overafd[overafd]][ent_BHK_overafd[0]] = 1
                df_output['PI for GI fra ' + shortnames_overafd[overafd]][ent_BHK_overafd[0]] = overafd

    for aa, afsn in enumerate(unique_afsn_NSR):
        if any(df_output['Genindlæggelse fra ' + shortnames_afsn[afsn]][ent_BHKID_alle] == 1):
            BHKI_afsn = np.unique(df_output['Hændelsesansvarligt afsnit'][ent_BHKID_alle])
            for bo, bhk_afsn in enumerate(BHKI_afsn):
                ent_BHK_afsn = np.where((df_output['Hændelsesansvarligt afsnit'] == bhk_afsn) &
                                           (df_kontakter['Behandlingskontakt record ID'] == BHKID))[0]
                selection = (df_output['Genindlæggelse fra ' + shortnames_afsn[afsn]][ent_BHK_afsn] != 1)

                # Erstat værdier der ikke matcher 'selection' ovenfor; så find GI flag og erstat dem.
                df_output['Genindlæggelse fra ' + shortnames_afsn[afsn]][ent_BHK_afsn] = \
                    df_output['Genindlæggelse fra ' + shortnames_afsn[afsn]][ent_BHK_afsn].where( selection , 99)  # replace values not matching selection
                df_output['PI for GI fra ' + shortnames_afsn[afsn]][ent_BHK_afsn] = \
                    df_output['PI for GI fra ' + shortnames_afsn[afsn]][ent_BHK_afsn].where( selection, 'Udeladt') # replace values not matching selection

                # Erstat første GI flag så det kun er registreret en gang per BHKID
                df_output['Genindlæggelse fra ' + shortnames_afsn[afsn]][ent_BHK_afsn[0]] = 1
                df_output['PI for GI fra ' + shortnames_afsn[afsn]][ent_BHK_afsn[0]] = afsn

print('\n - Færdig med korrigering af flag')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print('\n - Gemmer output dataframe til Excelfil')

if inkluderpersonoplysninger:
    print('   Sikrer at sorteringen er CPR > BHK udskrivning > Kontakt start')
    df_output = df_output.sort_values(by=['CPR','Behandlingskontakt udskrivning', 'Kontakt start'], ascending=[True, True, True])
else:
    print('   Sikrer at sorteringen er BHK udskrivning > Kontakt start')
    df_output = df_output.sort_values(by=['Behandlingskontakt udskrivning','Kontakt start'], ascending=[True, True])

path = "O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/genindlæggelser og genbesøg/Behandlingskontaktgenindlæggelser/"
df_output.to_excel(path+outputfilename, sheet_name="data output")
print('   Output skrevet til '+path+outputfilename)
#=======================================================================================================================

