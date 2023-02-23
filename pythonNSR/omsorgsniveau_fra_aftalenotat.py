#=======================================================================================================================
import pdb

import numpy as np
import pandas as pd
import sys
#=======================================================================================================================
def getON(excelfile=None):
    """
    Getting "omsorgsniveau" from input excel file with CPR, aftaledato, aftaleafsnit, aftale notat

    - Example of use:
    import omsorgsniveau_fra_aftalenotat as ofa
    ofa.getON()

    """
    if excelfile is None:
        excelfile = "O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/GynObs omsorgsniveau/Aftalanotater_GynObs_til_omsorgsniveau_230101-230223.xlsx"
        print('\n Advarsel: Ingen Excelfil givet så input sat til:\n   '+excelfile)
    # indlæs Excel data i pandas dataframe
    print('\n - Indlæser Excel datafil...')
    df_in = pd.read_excel(excelfile)
    print('   Indlæste dataframe med ' + str(len(df_in)) + ' rækker')

    out_cpr  = []
    out_date = []
    out_ward = []
    out_note = []
    out_on   = []

    for aa, aftalenotat in enumerate(df_in['Aftale notat']):
        infostr = '   Tjekker aftalenotat ' + str(aa + 1) + ' / ' + str(len(df_in['Aftale notat']))
        sys.stdout.write("%s\r" % infostr)
        sys.stdout.flush()
        aftalenotat = str(aftalenotat)
        if ("niveau" in aftalenotat.lower()) or ("niv" in aftalenotat.lower()):
            if "niveau" in aftalenotat.lower():
                on_val = aftalenotat.lower().split('niveau')[-1].strip().split(' ')[0]
            elif "niv." in aftalenotat.lower():
                on_val = aftalenotat.lower().split('niv.')[-1].strip().split(' ')[0]
            elif "niveua" in aftalenotat.lower():
                on_val = aftalenotat.lower().split('niveua')[-1].strip().split(' ')[0]
            elif "nivaeu" in aftalenotat.lower():
                on_val = aftalenotat.lower().split('nivaeu')[-1].strip().split(' ')[0]
            elif "niv" in aftalenotat.lower():
                on_val = aftalenotat.lower().split('niv')[-1].strip().split(' ')[0]

            if type(on_val) is not float:
                on_val_trimmed = on_val.split('-')[0].split('(')[0]
                stringstoremove = [')','.',',',':','?',';','+','fra']
                for ss, remstr in enumerate(stringstoremove):
                    on_val_trimmed = on_val_trimmed.replace(remstr,'')

                if on_val_trimmed == '':
                    on_out = np.nan
                else:
                    on_out = float(on_val_trimmed)

            out_cpr.append(df_in['CPR-nummer'][aa])
            out_date.append(df_in['Aftale dato og tid'][aa])
            out_ward.append(df_in['Aftale Afsnit'][aa])
            out_note.append(df_in['Aftale notat'][aa])
            out_on.append(on_out)

    print('\n - Færdig med tjek af aftalenotater ')

    print(' - Fandt '+str(len(out_on))+' omsorgsniveauer blandt aftalenotaterne')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' - Forbereder output')
    outputdata = {}
    outputdata['cpr'] = out_cpr
    outputdata['aftale dato og tid'] = out_date
    outputdata['aftale afsnit'] = out_ward
    outputdata['aftale notat'] = out_note
    outputdata['omsorgsniveau'] = out_on

    df_output = pd.DataFrame(outputdata)
    outfilename = excelfile.split('.xls')[0]+' - omsorgsniveau.xlsx'
    df_output.to_excel(outfilename, sheet_name="data output")
    print('\n - Output gemt i filen "'+outfilename+'"')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#=======================================================================================================================