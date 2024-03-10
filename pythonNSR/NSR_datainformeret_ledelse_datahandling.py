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
    datestamp_infile        Dato-stempel på formen 'yymmdd' for input inidkatorkataloget
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
    df_BI = pd.read_excel(filepath+BIfile)
    #df_BI.to_excel(filepath+BIfile.replace('.xlsx', '_BCKUP'+datestr+'.xlsx'), sheet_name="Indikatoroversigt")

    if verbose: print(' - combine: Indlæser *.outdata.csv filer; begrænser til indholdet for version '+outdatafil_version)
    filenamelist, dataframelist = ndld.load_data(filepath,outdatafil_version,verbose=verbose)

    df_outdataload = pandas.concat(dataframelist)
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
            outputdata["Udvikling ifht. Seneste måned"][cc] = outputdata["Målopfyldelse"][cc] - df_BI["Målopfyldelse"].values[ind_entBI][0]
            outputdata["Seneste måned"][cc] = df_BI["Aktuel værdi"].values[ind_entBI][0]
        elif len(ind_entBI) > 1:
            if verbose: print(' - combine WARNING: indikatornøgle "' + indicatorkey + '" optræder '+str(len(ind_entBI))+
                              ' gange i inputfil')
        else:
            if verbose: print(' - combine WARNING: indikatornøgle "'+indicatorkey+'" findes ikke inputfil')

    df_output = pd.DataFrame(outputdata).sort_values(by=['Enhed','Indikatornummer'])
    outfilename = filepath+BIfile.replace(datestamp_infile,datestr)
    df_output.to_excel(outfilename, sheet_name="Indikatoroversigt",index=False)
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
        df_ad = pd.read_csv(filepathname,sep=';',encoding='Windows-1252')
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
def afdelingsliste():
    """

    """
    afdlist = ['Administrationen', 'AKA', 'Anæstesi', 'BogU', 'Driftsafdelingen', 'Garantiklinikken',
               'GynObs', 'KBA', 'Kirurgi', 'M1', 'M2', 'M3', 'Multisygdom', 'NSR', 'Ort.Kir.']
    return afdlist

#=======================================================================================================================
def forbered_budgettal(sheetname,filepathname,dataversion,outpath,verbose=True):
    """
    Loader budegettal og forbereder outdata.csv fil

    -- EXAMPLE OF USE --
    import NSR_datainformeret_ledelse_datahandling as ndld

    sheetname='ArkMedOversigt'
    filename = 'NSRbudgetoversigt.xlsx'
    filepath = 'O:/Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2023 - Datainformeret ledelse/Data til DIL/'
    filepathname=filepath+filename
    outpath = filepath
    dataversion = '2024-03'
    csvfilename = ndld.forbered_budgettal(sheetname,filepathname,dataversion,outpath,verbose=True)

    """
    csvfilename = filepathname.split('/')[-1].replace('.xlsx','_outdata.csv')

    if verbose: print(' - Budget: Indlæser budget data i fanen '+sheetname+' fra filen: \n            '+filepathname)
    dic_in = pd.read_excel(filepathname,sheet_name=[sheetname])
    df_in  = dic_in[sheetname]

    Nrows_in = len(df_in['Prognosticeret resultat'])
    procentKRoverholdt = ['ikke_udfyldt']*Nrows_in
    for ii, res in enumerate(df_in['Prognosticeret resultat']):
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
    outputdata['Lavt_er_godt'] = [1]*Nrows_in
    outputdata['Minpunkt'] = [100]*Nrows_in
    outputdata['Mål rød gul'] = [0.5]*Nrows_in
    outputdata['Mål gul grøn'] = [0.001]*Nrows_in
    outputdata['Maxpunkt'] = [0.0]*Nrows_in
    outputdata['Aktuel værdi'] = [999]*Nrows_in
    outputdata['Kommentar'] = ['']*Nrows_in
    outputdata['Version'] = [dataversion]*Nrows_in

    if verbose: print(' - Budget: Løkke over budget resultater ')
    for rr, res in enumerate(df_in['Prognosticeret resultat']):
        outputdata['Enhed'][rr] = df_in['Navn afdeling'].values[rr]
        if (res > 500000) and (df_in['Navn afdeling'].values[rr] in ndld.afdelingsliste()): # use MegaKr for evaluation
            procentKRoverholdt[rr] = 'kr'
            current_val = res / 1e6
        elif (res < 0) and (df_in['Navn afdeling'].values[rr] in ndld.afdelingsliste()): # use 0.0 for evaluation as budget is kept
            procentKRoverholdt[rr] = 'overholdt'
            current_val = 0.0
        elif df_in['Navn afdeling'].values[rr] in ndld.afdelingsliste(): # use percentage for evaluation
            current_val = df_in['Forhold budget'].values[rr]
        else:
            if verbose: print(' - Budget WARNING: Afdeling '+str(df_in['Navn afdeling'].values[rr])+' optræder ikke i afdelingsliste og bliver derfor ikke gemt i output')
            current_val = "ukendt_afdeling"

        outputdata['Aktuel værdi'][rr] = current_val

    outputfilename = outpath+csvfilename
    if verbose: print(' - Budget: output bliver skrevet til CSV filen '+outputfilename)
    df_output = pd.DataFrame(outputdata).sort_values(by=['Enhed', 'Indikatornummer'])
    df_output = df_output[df_output['Aktuel værdi'] != 'ukendt_afdeling']

    #- - - - - - - - - - - - - - - - - - - - - - - - - - -
    # check if csv file exists and correct output according to what already exists in output
    if os.path.isfile(outputfilename):
        if verbose: print('\n - Budget: CSV fil eksisterer allerede; tilføjer output til eksisterende fil')
        df_outdata = pd.read_csv(outputfilename, sep=';', encoding='Windows-1252')

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

    df_output.to_csv(outputfilename, sep=';',encoding='Windows-1252', index=False)
    print(' - Budget: Opdaterer indikatoroversigt til BI oversigt gemt i filen "'+outputfilename+'"')

    return csvfilename

#=======================================================================================================================