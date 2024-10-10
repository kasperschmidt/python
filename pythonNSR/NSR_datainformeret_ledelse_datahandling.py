#=======================================================================================================================
#from importlib import reload
import os.path

import pandas
import pandas as pd
import glob
import datetime
from os import path as ospath
from sys import stdout as sysstdout
import numpy as np
import pdb

import NSR_datainformeret_ledelse_datahandling as ndld
#=======================================================================================================================
def combine_output_to_TARGITupload(filepath,datestamp_infile,outdatafil_version,verbose=True):
    """
    filepath                stinavn til bibliotek med filer til kørsel (både in- og output)
    datestamp_infile        Dato-stempel på formen 'yymmdd' for input indikatorkataloget
    outdatafil_version      Versionsstempel for indiaktorer i '*outdata.csv' filer der skal bruges til opdatering

    --- EXAMPLE OF USE ---
    import NSR_datainformeret_ledelse_datahandling as ndld
    filepath='O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2023 - Datainformeret ledelse/Data til DIL/'
    outBIfile = ndld.combine_output_to_TARGITupload(filepath,'240201','2024-03')

    """
    nowstr  = datetime.datetime.strftime(datetime.datetime.now(), "%d-%m-%Y %H:%M:%S")
    datestr = datetime.datetime.strftime(datetime.date.today(), "%y%m%d")
    BIfile  = 'NSR datainformeret ledelse - indikatoroversigt '+datestamp_infile+'.xlsx'

    if verbose: print(" - combine: Indlæser input TARGIT fil der skal opdateres: "+filepath+BIfile)
    df_BI = pd.read_excel(filepath+'/../indikatoroversigt/'+BIfile)
    #df_BI.to_excel(filepath+BIfile.replace('.xlsx', '_BCKUP'+datestr+'.xlsx'), sheet_name="Indikatoroversigt")

    if verbose: print(' - combine: Indlæser *.outdata.csv filer; begrænser til indholdet for version '+outdatafil_version)
    filenamelist, dataframelist = ndld.load_data(filepath,outdatafil_version,verbose=verbose)

    df_outdataload = pandas.concat(dataframelist)
    if verbose: print(' - combine: De '+str(len(df_outdataload.columns))+
                      ' kolonner i det samlede dataset er:'+str(df_outdataload.columns))
    if verbose: print(' - combine: Som et tjek bør der være 11 kolonner (240313)')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' - Forbereder output')
    outputdata = {}

    for colname in dataframelist[0].columns:
        outputdata[colname.replace('_',' ')] = df_outdataload[colname]

    outputdata["Målopfyldelse seneste måned"] = [999]*len(outputdata['Enhed'])
    outputdata["Målopfyldelse"] = [999]*len(outputdata['Enhed'])
    outputdata["Udvikling ifht. Seneste måned"] = [999]*len(outputdata['Enhed'])
    outputdata["Seneste måned"] = [999]*len(outputdata['Enhed'])
    #outputdata["outdataCSVversion"] = [outdatafil_version]*len(outputdata['Enhed'])
    outputdata["senestemånedDatostreng"] = [datestamp_infile] * len(outputdata['Enhed'])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' - combine: Sætter inidkatornumre manuelt (overskriver fra input)')
    afdlist = ndld.afdelingsliste()
    for aa, afdeling in enumerate(afdlist): # looping through departments
        ent_afd = np.where(outputdata['Enhed'] == afdeling)[0]
        list_nummer = outputdata['Indikatornummer'].values[ent_afd]
        list_navn = outputdata['Indikatornavn'].values[ent_afd]
        runningnumber = 4
        for ll, lnavn in enumerate(list_navn):
            # if statement manually setting mandatory indicators
            if lnavn == 'Budget':
                list_nummer[ll] = 1
            elif lnavn == 'Sygefravær':
                list_nummer[ll] = 2
            elif lnavn == 'Vikar, FEA, OA, MA':
                list_nummer[ll] = 3
            else:
                list_nummer[ll] = runningnumber
                runningnumber = runningnumber + 1 # increment running number for non-mandatory indicators

        #outputdata['Indikatornummer'].values[ent_afd] = list_nummer # overskriv Indikatornummer med ny nummerliste
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    print(' - combine: Beregninger målopfyldelse for alle indikatorer')
    for cc, lavtgodt in enumerate(outputdata['Lavt er godt']):
        #beregn målopfyldelse
        minpoint = outputdata["Minpunkt"].values[cc]
        goal_redyellow = outputdata["Mål rød gul"].values[cc]
        goal_yellowgreen = outputdata["Mål gul grøn"].values[cc]
        maxpoint = outputdata["Maxpunkt"].values[cc]
        current_val = outputdata["Aktuel værdi"].values[cc]
        mo_val = ndld.maalopfyldelse(lavtgodt, minpoint, goal_redyellow, goal_yellowgreen, maxpoint, current_val, verbose=False)

        outputdata["Målopfyldelse"][cc] = mo_val

        # led efter indikator i eksisterende BI fil
        indicatorkey     = outputdata['Enhed'].values[cc]+outputdata['Indikatornavn'].values[cc]+outputdata['Version'].values[cc]
        indicatorkeys_BI = df_BI['Enhed'] + df_BI['Indikatornavn']
        ind_entBI = np.where(indicatorkeys_BI == indicatorkey[:-7])[0]
        if len(ind_entBI) == 1:
            outputdata["Målopfyldelse seneste måned"][cc] = df_BI["Målopfyldelse"].values[ind_entBI][0]
            outputdata["Udvikling ifht. Seneste måned"][cc] = np.round(outputdata["Målopfyldelse"][cc] - df_BI["Målopfyldelse"].values[ind_entBI][0], 2)
            outputdata["Seneste måned"][cc] = df_BI["Aktuel værdi"].values[ind_entBI][0]
        elif len(ind_entBI) > 1:
            if verbose: print(' - combine WARNING: indikatornøgle "' + indicatorkey + '" optræder '+str(len(ind_entBI))+
                              ' gange i inputfil')
        else:
            if verbose: print(' - combine WARNING: indikatornøgle "'+indicatorkey+'" findes ikke inputfil')

    df_output = pd.DataFrame(outputdata).sort_values(by=['Enhed','Indikatornummer'])
    outfilename = filepath+'/../indikatoroversigt/'+BIfile.replace(datestamp_infile,datestr)
    df_output.to_excel(outfilename, sheet_name="Indikatoroversigt",index=False, encoding='Windows-1252')
    print(' - combine: Opdaterer indikatoroversigt til BI oversigt gemt i filen "'+outfilename+'"')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    return outfilename
#=======================================================================================================================
def load_data(filepath,outdatafil_version,verbose=True):
    """
    -- EXAMPLE OF USE --
    import NSR_datainformeret_ledelse_datahandling as ndld
    filepath = 'O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2023 - Datainformeret ledelse/Data til DIL/'
    filenamelist, dataframelist = ndld.load_data(filepath,'2024-03')
    """
    dataframelist  = []
    filenamelist   = []
    filelist       = glob.glob(filepath+'*outdata.csv')
    filelist       = [fl.replace('\\','/') for fl in filelist]
    Nfiles = len(filelist)

    if verbose: print(' - load_data: Løkke over de '+str(Nfiles)+' fundne filer')
    if verbose: print(' - load_data: Begrænser indlæst data til version ' + outdatafil_version)
    for ff, filepathname in enumerate(filelist):
        if verbose: print(' -            Indlæser '+filepathname.split('/')[-1])
        df_ad = pd.read_csv(filepathname,sep=';',encoding='Windows-1252')

        for colname in df_ad.columns:
            if 'Mål_rød' in colname:
                df_ad.rename(columns={colname: 'Mål rød gul'}, inplace=True)
            if 'rød-gul' in colname:
                df_ad.rename(columns={colname: 'Mål rød gul'}, inplace=True)
            if 'gul-grøn' in colname:
                df_ad.rename(columns={colname: 'Mål gul grøn'}, inplace=True)
            if 'Mål_gul' in colname:
                df_ad.rename(columns={colname: 'Mål gul grøn'}, inplace=True)
            if 'godt ' in colname:
                df_ad.rename(columns={colname: 'Lavt er godt'}, inplace=True)
            if 'Lavt_er' in colname:
                df_ad.rename(columns={colname: 'Lavt er godt'}, inplace=True)
            if 'Aktuel_værdi' in colname:
                df_ad.rename(columns={colname: 'Aktuel værdi'}, inplace=True)

        filenamelist.append(filepathname.split('/')[-1])
        dataframelist.append(df_ad[df_ad['Version'].values == outdatafil_version])

    if verbose: print(' - load_data: Færdig med at indlæse filer')
    return filenamelist, dataframelist

#=======================================================================================================================
def farveintervaller(verbose=True):
    """
    Intervaller brugt til at bestemme ændring i farver for indiaktoroversigter og beregning af målopfyldelse
    """
    farve_dic = {}
    farve_dic['minpoint']  = 0.00
    farve_dic['redyellow']   = 20.00
    farve_dic['yellowgreen']  = 80.00
    farve_dic['maxpoint'] = 100.00

    return farve_dic

#=======================================================================================================================
def tjektal(calcparam,verbose=True):
    """
    Tjek om parameter er tekstren og så konverter til tal
    """
    if isinstance(calcparam, str):
        calcparam = float(calcparam.replace(',', '.'))

    return calcparam
#=======================================================================================================================
def maalopfyldelse(low_is_good, minpoint, goal_redyellow, goal_yellowgreen, maxpoint, current_val,verbose=True):
    """
    Beregning af målopfyldelse for givent input

    D = low_is_good
    H = minpoint
    I = goal_redyellow
    J = goal_yellowgreen
    K = maxpoint
    M = current_val


    """
    if verbose: print(' - maalopfyldelse: Beregner maalopfyldelse for input variable')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - maalopfyldelse: Sikrer at input parametre er floats')
    low_is_good = ndld.tjektal(low_is_good)
    minpoint = ndld.tjektal(minpoint)
    goal_redyellow = ndld.tjektal(goal_redyellow)
    goal_yellowgreen = ndld.tjektal(goal_yellowgreen)
    maxpoint = ndld.tjektal(maxpoint)
    current_val = ndld.tjektal(current_val)
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    colintervals = ndld.farveintervaller(verbose=verbose)

    if low_is_good == 1:
        if current_val >= goal_redyellow:
            mo_val = colintervals['redyellow'] - colintervals['redyellow'] / (minpoint-goal_redyellow) * (current_val-goal_redyellow)
        elif (current_val < goal_redyellow) and (current_val > goal_yellowgreen):
            mo_val = colintervals['yellowgreen'] - (colintervals['yellowgreen']-colintervals['redyellow']) / (goal_yellowgreen-goal_redyellow) * (goal_yellowgreen-current_val)
        else:
            mo_val = colintervals['yellowgreen'] + (colintervals['maxpoint']-colintervals['yellowgreen']) / (goal_yellowgreen-maxpoint) * (goal_yellowgreen-current_val)
    else:
        if current_val <= goal_redyellow:
            mo_val = colintervals['redyellow'] - colintervals['redyellow'] / (minpoint-goal_redyellow) * (current_val-goal_redyellow)
        elif (current_val > goal_redyellow) and (current_val < goal_yellowgreen):
            mo_val = colintervals['yellowgreen'] - (colintervals['yellowgreen']-colintervals['redyellow']) / (goal_yellowgreen-goal_redyellow) * (goal_yellowgreen-current_val)
        else:
            mo_val = colintervals['yellowgreen'] + (colintervals['maxpoint']-colintervals['yellowgreen']) / (goal_yellowgreen-maxpoint) * (goal_yellowgreen-current_val)

    if verbose: print(' - maalopfyldelse: Resultat = '+str(mo_val)+' for følgende input parametre:')
    if verbose: print([low_is_good, minpoint, goal_redyellow, goal_yellowgreen, maxpoint, current_val])

    return mo_val
#=======================================================================================================================
def afdelingsliste(version=None):
    """

    """
    if version == 'budget':
        afdlist = ['Administrationen','Akutafdelingen','Anæstesiologisk afdeling','Børneafdelingen',
                   'Driftsafdelingen','Gynækologisk/Obstetrisk afdeling',
                   'Klinisk Biokemisk afdeling','Kirurgisk afdeling','Medicin 1','Medicin 2','Medicin 3 - Fysergo',
                   'CMKS', 'Total NSR', 'Ortopædkirurgisk afdeling']
    else:
        afdlist = ['Administrationen', 'AKA', 'Anæstesi', 'BogU', 'Driftsafdelingen', 'Garantiklinikken',
                   'GynObs', 'KBA', 'Kirurgi', 'M1', 'M2', 'M3', 'Multisygdom', 'NSR', 'Ort.Kir.']

    return afdlist

#=======================================================================================================================
def forbered_budgettal(sheetname,filepathname,dataversion,outpath,verbose=True):
    """
    Loader budegettal og forbereder outdata.csv fil

    -- EXAMPLE OF USE --
    import NSR_datainformeret_ledelse_datahandling as ndld

    sheetname='Afdelinger samlet'
    filename = 'Økonomistatus 2024 - marts.xlsx'
    filepath = 'O:/Administration/02 - Økonomi og Planlægning/02 Økonomi/01 Budget og Regnskab/2024/Prognoser/'
    filepathname=filepath+filename
    outpath = 'O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2023 - Datainformeret ledelse/Data til DIL/'
    dataversion = '2024-04'
    csvfilename = ndld.forbered_budgettal(sheetname,filepathname,dataversion,outpath,verbose=True)


    sheetname='ArkMedOversigt'
    filename = 'NSRbudgetoversigt.xlsx'
    filepath = 'O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2023 - Datainformeret ledelse/Data til DIL/'
    filepathname=filepath+filename
    outpath = filepath
    dataversion = '2024-03'
    csvfilename = ndld.forbered_budgettal(sheetname,filepathname,dataversion,outpath,verbose=True)

    """
    csvfilename = filepathname.split('/')[-1].split(' - ')[0]+'_outdata.csv'
    if verbose: print(' - Budget: Indlæser budget data i fanen '+sheetname+' fra filen: \n            '+filepathname)
    dic_in = pd.read_excel(filepathname,sheet_name=[sheetname],usecols='C:F',header=2,decimal=",")
    df_in  = dic_in[sheetname][:25]

    Nrows_in = len(df_in['Difference'])
    procentKRoverholdt = ['ikke_udfyldt']*Nrows_in
    for ii, res in enumerate(df_in['Difference']):
        if res > 500000:
            procentKRoverholdt[ii] = 'kr'
        elif res < 0:
            procentKRoverholdt[ii] = 'overholdt'
        else:
            procentKRoverholdt[ii] = 'procent'

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Budget: Forbereder output')
    outputdata = {}

    outputdata['Enhed'] = [999]*Nrows_in
    outputdata['Indikatornummer'] = [1]*Nrows_in
    outputdata['Indikatornavn'] = ['Budget']*Nrows_in
    outputdata['Lavt er godt'] = [1]*Nrows_in
    outputdata['Minpunkt'] = [100]*Nrows_in
    outputdata['Mål rød gul'] = [0.5]*Nrows_in
    outputdata['Mål gul grøn'] = [0.001]*Nrows_in
    outputdata['Maxpunkt'] = [0.0]*Nrows_in
    outputdata['Aktuel værdi'] = [999]*Nrows_in
    outputdata['Kommentar'] = ['']*Nrows_in
    outputdata['Version'] = [dataversion]*Nrows_in

    if verbose: print(' - Budget: Løkke over budget resultater ')
    Outname = ndld.afdelingsliste()
    Inname  = ndld.afdelingsliste(version='budget')

    for rr, res in enumerate(df_in['Difference']):
        unitindex = np.where(np.asarray(Inname) == df_in['Afdelinger'].values[rr])[0]
        if len(unitindex) == 1:
            outputdata['Enhed'][rr] = Outname[unitindex[0]]
        else:
            outputdata['Enhed'][rr] = df_in['Afdelinger'].values[rr]

        if (res > 500000) and (df_in['Afdelinger'].values[rr] in Inname): # use MegaKr for evaluation
            procentKRoverholdt[rr] = 'kr'
            current_val = res / 1e6
        elif (res < 0) and (df_in['Afdelinger'].values[rr] in Inname): # use 0.0 for evaluation as budget is kept
            procentKRoverholdt[rr] = 'overholdt'
            current_val = 0.0
        elif df_in['Afdelinger'].values[rr] in Inname: # use percentage for evaluation
            current_val = (df_in['Forventet Forbrug'].values[rr] / df_in['Forventet Budget'].values[rr] - 1) *100
        else:
            if verbose: print(' - Budget WARNING: Afdeling '+str(df_in['Afdelinger'].values[rr])+' optræder ikke i afdelingsliste og bliver derfor ikke gemt i output')
            current_val = "ukendt_afdeling"

        if current_val == "ukendt_afdeling":
            outputdata['Aktuel værdi'][rr] = current_val
        else:
            outputdata['Aktuel værdi'][rr] = float("{:.4f}".format(current_val))

    outputfilename = outpath+csvfilename
    if verbose: print(' - Budget: output bliver skrevet til CSV filen '+outputfilename)
    df_output = pd.DataFrame(outputdata).sort_values(by=['Enhed', 'Indikatornummer'])
    df_output = df_output[df_output['Aktuel værdi'] != 'ukendt_afdeling']

    #- - - - - - - - - - - - - - - - - - - - - - - - - - -
    # check if csv file exists and correct output according to what already exists in output
    if os.path.isfile(outputfilename):
        if verbose: print('\n - Budget: CSV fil eksisterer allerede; tilføjer output til eksisterende fil')
        df_outdata = pd.read_csv(outputfilename, sep=';', encoding='Windows-1252',decimal=".")

        ukeys_exists = df_outdata['Enhed']+df_outdata['Indikatornummer'].values.astype(int).astype(str)+df_outdata['Indikatornavn'].values+df_outdata['Version']
        ukeys_thisrun = df_output['Enhed'] + df_output['Indikatornummer'].values.astype(int).astype(str) + df_output['Indikatornavn'].values + df_output['Version']
        key_remove_list = []
        for kk, ukey in enumerate(ukeys_thisrun):
            if ukey in ukeys_exists.values:
                if verbose: print(' - Budget WARNING: Kombinationen "'+ukey+'" findes allerede i CSV fil og vil derfor ikke blive overskrevet. Brug alternativ "Version" eller ret i eksisterende CSV fil.')
                key_remove_list.append(df_output.index[kk])
        df_output = df_output.drop(key_remove_list)

        df_output  = pandas.concat([df_outdata,df_output])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    df_output.to_csv(outputfilename, sep=';',encoding='Windows-1252', decimal=".", index=False)
    print(' - Budget: Opdaterer indikatoroversigt til BI oversigt gemt i filen "'+outputfilename+'"')

    return csvfilename

#=======================================================================================================================
def generer_tekstbidder(indikatorkatalog,outpath,verbose=True):
    """
    Funktion til at generere teksbidder til dokumentaion og sider i datainformeret ldelse baseret på indikatorkatalog

    --- EXAMPLE OF USE ---
    import NSR_datainformeret_ledelse_datahandling as ndld
    filename = 'NSR datainformeret ledelse - indikatoroversigt 240313.xlsx'
    filepath = 'O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2023 - Datainformeret ledelse/Data til DIL/'
    tekstfil = ndld.generer_tekstbidder(filepath+filename,filepath,verbose=True)

    """
    if verbose: print(' - gentekst: Indlæser indikatorkatalog')
    df_indi = pandas.read_excel(indikatorkatalog, sheet_name="Indikatoroversigt")

    if verbose: print(' - gentekst: Forbereder output ')
    outfile = outpath+'tekst_til_doc.txt'
    fout = open(outfile,'w')

    if verbose: print(' - gentekst: Indlæser indikatorkatalog')
    for ee, afdeling in enumerate(ndld.afdelingsliste()):
        fout.write('\n\n===================== '+afdeling+' ===================== \n\n')
        afd_ent = np.where(df_indi['Enhed'] == afdeling)
        afd_inidi = df_indi['Indikatornavn'][afd_ent[0]]
        goal_redyellow = df_indi['Mål rød gul'][afd_ent[0]].values
        goal_yellowgreen = df_indi['Mål gul grøn'][afd_ent[0]].values
        lavtgodt = df_indi['Lavt er godt'][afd_ent[0]].values
        for ii, indikator in enumerate(afd_inidi):
            if indikator.lower() == 'budget':
                fout.write("""
- Budget:
    Grøn: Budget overholdt
    Gul: Budget overskredet med op til 0,5% eller 500.000kr.
    Rød: Budget overskredet med mere end 0,5% eller mere end 500.000kr.
""")
            elif indikator.lower() == 'sygefravær':
                fout.write("""
- Sygefravær:
    Grøn: Reduktion > 1,0%-point ifht. sidste år
    Gul: Reduktion mellem 1,0 og 0,0%-point ifht. sidste år.
    Rød: Forøgelse af fravær ifht. sidste år
""")
            elif indikator.lower() == 'vikar, fea, oa, ma':
                fout.write("""
- Vikar, FEA, OA, MA:
    Grøn: Reduktion >46%-point ifht. 2022 niveau
    Gul: Reduktion over 2023 niveau men under 46% af 2022 niveau
    Rød: Reduktion mellem 2022 og 2023 niveau
""")
            else:
                if lavtgodt[ii] == 1:
                    fout.write("""
- {}:
    Grøn: Under {}
    Gul: Mellem {} og {}
    Rød: Over {}
""".format(indikator,goal_yellowgreen[ii],goal_yellowgreen[ii],goal_redyellow[ii],goal_redyellow[ii]))
                else:
                    fout.write("""
- {}:
    Grøn: Over {}
    Gul: Mellem {} og {}
    Rød: Under {}
""".format(indikator,goal_yellowgreen[ii],goal_redyellow[ii],goal_yellowgreen[ii],goal_redyellow[ii]))


    if verbose: print(' - gentekst: Output skrevet til '+outfile)
    fout.close()
#=======================================================================================================================

