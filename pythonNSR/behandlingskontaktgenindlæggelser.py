#=======================================================================================================================
#from importlib import reload
import easygui
import pandas as pd
import datetime
import pdb
# import os
import sys
import numpy as np

#=======================================================================================================================
def tjek_for_PI(dia_udsk,dia_alle,afsnit_alle):
    """
    Tjek indhold af patientkotakt for at bestæmme om der er tale om en primærindlæggelse (PI)
    """
    result = 0
    if (str(dia_udsk.values)[0][:3] != 'DO4') & \
            any([str(dia)[:5] != 'DZ763' for dia in dia_alle.values]) & \
            any([str(dia)[:2] != 'DC' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD00' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD02' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD03' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD04' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD05' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD06' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD07' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD08' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD09' for dia in dia_alle.values]) & \
            any(['hospice' not in str(afsn).lower() for afsn in afsnit_alle.values]):
        result = 1
    return result
#=======================================================================================================================
def tjek_for_GI(dia_ptktk,ptktktype_start,dia_alle,afsnit_alle):
    """
    Tjek indhold af patientkontakt for at bestæmme om der er tale om en genindlæggelse (GI)
    """
    result = 0
    if (startkontaktype == 'Akut') & \
            any([str(dia)[:5] != 'DZ763' for dia in dia_all.values]) & \
            any([str(dia)[:2] != 'DF' for dia in [dia_alle[0],dia_alle[-1]]]) & \
            any([str(dia)[:2] != 'DC' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD00' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD02' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD03' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD04' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD05' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD06' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD07' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD08' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DD09' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DO80' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DO81' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DO82' for dia in dia_alle.values]) & \
            any([str(dia)[:4] != 'DO84' for dia in dia_alle.values]) & \
            (str(dia_alle[0])[:2] != 'DS'  ) & \
            (str(dia_alle[0])[:3] != 'DT1' ) & \
            (str(dia_alle[0])[:3] != 'DT2' ) & \
            (str(dia_alle[0])[:4] != 'DT30') & \
            (str(dia_alle[0])[:4] != 'DT31') & \
            (str(dia_alle[0])[:4] != 'DT32') & \
            (str(dia_alle[0])[:4] != 'DT33') & \
            (str(dia_alle[0])[:4] != 'DT51') & \
            (str(dia_alle[0])[:4] != 'DT52') & \
            (str(dia_alle[0])[:4] != 'DT53') & \
            (str(dia_alle[0])[:4] != 'DT54') & \
            (str(dia_alle[0])[:4] != 'DT55') & \
            (str(dia_alle[0])[:4] != 'DT56') & \
            (str(dia_alle[0])[:4] != 'DT57') & \
            (str(dia_alle[0])[:4] != 'DT58') & \
            (str(dia_alle[0])[:4] != 'DT59') & \
            (str(dia_alle[0])[:3] != 'DT6' ) & \
            (str(dia_alle[0])[:3] != 'DT7' ) & \
            (str(dia_alle[0])[:4] != 'DT81') & \
            (str(dia_alle[0])[:3] != 'DT9' ) & \
            (str(dia_alle[0])[:2] != 'DX'  ) & \
            (str(dia_alle[0])[:2] != 'DY'  ) & \
            any(['hospice' not in str(afsn).lower() for afsn in afsnit_alle.values]):
        result = 1
    return result
#=======================================================================================================================
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
nowstring   = datetime.datetime.strftime(datetime.datetime.now(),"%d-%m-%Y %H:%M:%S")
todaystring = datetime.datetime.strftime(datetime.date.today(),"%y%m%d")
print("\n\n - Program til indentificering af behandlingskontaktgenindlæggelser startet "+nowstring)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
GUIinput = False
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
    excelfile = "O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/genindlæggelser og genbesøg/Behandlingskontaktgenindlæggelser/Behandlingskontaktgenindlæggelser_datatræk_220401-220601_small.xlsx"

print(" - Excel datafil angivet: \n   "+ excelfile)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# indlæs Excel data i pandas dataframe
print('\n - Indlæser Excel datafil...')
df_kontakter = pd.read_excel(excelfile)
print('   Indlæste dataframe med ' + str(len(df_kontakter)) + ' rækker')
print('   Sikrer at sorteringen er CPR > BHK udskrivning > PTK ID > Kontakt start')
df_kontakter.sort_values(by=["Patient CPR-nr.", "Behandlingskontakt udskrivningsdato Dato-tid", "Patientkontakt record ID", "Kontakt startdato Dato-tid"], ascending=[True, True, True, True])

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# definer output dataframe
outputfilename = excelfile.replace('.xls','_BGIberegninger'+todaystring+'.xls')
outputdata = {}
NaNlist  = [np.nan] * len(df_kontakter['Hændelsesansvarlig Afsnit navn'])
zerolist = [0] * len(df_kontakter['Hændelsesansvarlig Afsnit navn'])
outputdata['Kontakt startdato Dato-tid'] = list(df_kontakter['Kontakt startdato Dato-tid'])
outputdata['Ansvarligt afsnit'] = list(df_kontakter['Hændelsesansvarlig Afsnit navn'])
outputdata['Primærindlæggelse'] = zerolist

unique_wards = np.unique(df_kontakter['Hændelsesansvarlig Overafdeling navn'])
shortnames = {}
for ww, ward in enumerate(unique_wards):
    shortnames[ward] = ward.split(' - ')[0].split(', ')[-1]
    if ('SLA ' in ward) or ('NAE ' in ward):
        outputdata['Genindlæggelse fra '+shortnames[ward]] = NaNlist
        outputdata['PI for GI fra ' + shortnames[ward]] = NaNlist

df_output = pd.DataFrame(outputdata)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# definer en række vektorer med indhold til senere brug
uniqueCPR = np.unique(df_kontakter['Patient CPR-nr.'])
uniqueBHKID = np.unique(df_kontakter['Behandlingskontakt record ID'])
udskrivningBHK = df_kontakter['Behandlingskontakt udskrivningsdato Dato-tid']
kontaktStart = df_output['Kontakt startdato Dato-tid']

len(np.where(df_output['Kontakt startdato Dato-tid'] > datetime.datetime.strptime("25-05-2022 00:00:00", "%d-%m-%Y %H:%M:%S"))[0])
len(np.where(df_kontakter['Behandlingskontakt udskrivningsdato Dato-tid'] > datetime.datetime.strptime("25-05-2022 00:00:00", "%d-%m-%Y %H:%M:%S"))[0])

unique_wards = np.unique(df_kontakter['Behandlingskontaktansvarlig Overafdeling navn'])

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Løkke over behandlingskontakter for at bestæmme primærindlæggelser
print('\n - Identificerer primærindlæggelser')
for bb, bid in enumerate(uniqueBHKID):
    infostr = '   Evaluerer udskrivende afsnit af behandlingskontakt ' + str(bb + 1) + ' / ' + str(len(uniqueBHKID))
    sys.stdout.write("%s\r" % infostr)
    sys.stdout.flush()

    # Vælger sidste (udskrivende) patientkontakt i behandlingskontak
    ent_ptk_alle = np.where((df_kontakter['Behandlingskontakt record ID'] == bid))[0]
    ent_ptk_udsk = np.where((df_kontakter['Behandlingskontakt record ID'] == bid) & (df_kontakter['Hændelsestype navn'] == 'UDSKRIVNING'))[0]

    isPI = tjek_for_PI(df_kontakter['Aktionsdiagnosekode'][ent_ptk_udsk],
                       df_kontakter['Aktionsdiagnosekode'][ent_ptk_alle],
                       df_kontakter['Hændelsesansvarlig Afdeling navn'][ent_ptk_alle])
    df_output['Primærindlæggelse'][ent_ptk_udsk] = isPI

print('\n - Færdig med evaluering af alle '+str(len(uniqueBHKID))+' behandlingskontakter')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Løkke over patienter for at evaluere dem med henblik på genindlæggelser
print('\n - Identificerer primærindlæggelser')
for cc, cpr in enumerate(uniqueCPR):
    infostr = '   Evaluerer forløb for patient ' + str(cc + 1) + ' / ' + str(len(uniqueCPR))
    sys.stdout.write("%s\r" % infostr)
    sys.stdout.flush()

    ent_cpr = np.where(df_kontakter['Patient CPR-nr.'] == cpr)[0]
    ent_PI  = np.where((df_kontakter['Patient CPR-nr.'] == cpr) & (df_output['Primærindlæggelse'] == 1))[0]
    u_BHKID = np.unique(df_kontakter['Behandlingskontakt record ID'][ent_cpr])
    u_PIafsnit = np.unique(df_kontakter['Hændelsesansvarlig Overafdeling navn'][ent_PI])

    if len(u_BHKID) > 1:  # kun tjek genindlæggelser for patienter med mere end en behandlingskontakt
        for aa, PIafs in enumerate(u_PIafsnit):
            for ee, ent_ktk in enumerate(ent_cpr):
                ent_BHKID = df_kontakter['Behandlingskontakt record ID'][ent_cpr]

                diffPItider = kontaktstart[ent_ktk] - udskrivningBHK[np.where((df_kontakter['Patient CPR-nr.'] == cpr) &
                                                                              (df_output['Primærindlæggelse'] == 1) &
                                                                              (df_kontakter['Primærindlæggelse'] == df_kontakter['Hændelsesansvarlig Overafdeling navn'][ent_cpr]))]

                if any(diffPItider < 720 hours... HERE221006)
                result = tjek_for_GI()

                    df_output['Genindlæggelse fra ' + shortnames[PIafs] = result

print('\n - Færdig med evaluering af all '+str(len(uniqueCPR))+' patienter')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print('\n - Gemmer output dataframe til Excelfil')
path = "O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/genindlæggelser og genbesøg/Behandlingskontaktgenindlæggelser/"
df_output.to_excel(outputfilename, sheet_name="data output")
print('   Output skrevet til '+outputfilename)
#=======================================================================================================================

