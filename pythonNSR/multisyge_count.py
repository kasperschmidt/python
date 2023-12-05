#=======================================================================================================================
#from importlib import reload
import pandas as pd
import datetime
from os import path as ospath
import numpy as np
import pdb
import sys
import multisyge_count as mc
#=======================================================================================================================
def count_multimorbid_patients(verbose=True, testrun=True):
    """

    -- EXAMPLE OF USE --
    import multisyge_count as mc
    mc.count_multimorbid_patients(testrun=True)
    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    excelpath = "O:\Administration/02 - Økonomi og Planlægning/01 Fælles/05 Arbejdsgrupper og projekter/2023 - Multisyge patienter og kronikere/"
    if testrun:
        excelfile = excelpath + 'Kronikere - lille til test.xlsx'
    else:
        excelfile = excelpath + 'Kronikere.xlsx'

    print(" - Excel datafil angivet: \n   " + excelfile)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # læs Excel data af kronikere in i dataframe
    print('\n - Indlæser Excel datafil...')
    sheet2read = 'Datatræk LPR3'
    df_input = pd.read_excel(excelfile,sheet_name='Datatræk LPR3', skiprows=2)

    if verbose: print(' - Fandt '+str(len(df_input))+' rækker i Excelfil i arket '+sheet2read)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    todaystring = datetime.datetime.strftime(datetime.date.today(), "%y%m%d")
    outpath = ospath.dirname(excelfile)
    outputfile = outpath+'/multisygepatienter_'+todaystring+'.xlsx'
    if verbose: print(' - Forbereder output. Multisyge patienter vil blive skrevet til Excelfil:\   '+outputfile)
    out_ktkstart = []
    out_cprlist = []
    out_dialist = []
    out_diagrouplist = []
    out_groupnamelist = []
    out_SOR1 = []  #"Org_Ansv_SOROrgEnhedNiv1Tekst"
    out_SOR2 = []  #"Org_Ansv_SOROrgEnhedNiv2Tekst"
    out_SOR3 = []  # "Org_Ansv_SOROrgEnhedNiv3Tekst"
    out_SOR4 = []  # "Org_Ansv_SOROrgEnhedNiv4Tekst"


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    uniqueCPR = np.unique(df_input['CPRnummer'])
    NuniqueCPR = len(uniqueCPR)
    Nmultimorbid = 0 # counter to keep track of multimorbid patients

    if verbose: print('\n - Tjekker CPR numre for multisygdom')
    for cc, cprno in enumerate(uniqueCPR):
        infostr = '   Tjekker CPR nummer ' + str(cc + 1) + ' / ' + str(NuniqueCPR)
        sys.stdout.write("%s\r" % infostr)
        sys.stdout.flush()

        ent = np.where(df_input['CPRnummer'] == cprno)[0]
        if len(ent) >= 2:
            ktk_start = np.asarray([datetime.datetime.strptime(datestr, "%Y-%m-%d") for datestr in df_input['Konstaktstart dato'][ent]])
            diagnoses = df_input['aDiag_Diagnose_Kode'][ent]
            diaGroups = np.asarray([mc.diaggroup(dia)[0] for dia in diagnoses])
            diaNames  = np.asarray([mc.diaggroup(dia)[1] for dia in diagnoses])
            uniqueDiaGroups = np.unique(diaGroups[np.where(diaGroups != -99)])
            Ngroups = len(uniqueDiaGroups)
            #SOR sections
            pt_SOR1 = df_input['Org_Ansv_SOROrgEnhedNiv1Tekst'][ent]
            pt_SOR2 = df_input['Org_Ansv_SOROrgEnhedNiv2Tekst'][ent]
            pt_SOR3 = df_input['Org_Ansv_SOROrgEnhedNiv3Tekst'][ent]
            pt_SOR4 = df_input['Org_Ansv_SOROrgEnhedNiv4Tekst'][ent]

            if Ngroups >=2:
                Nmultimorbid = Nmultimorbid + 1

                for dd, diagroup in enumerate(diaGroups):
                    out_cprlist.append(cprno)
                    out_ktkstart.append(datetime.datetime.strftime(ktk_start[dd],"%Y-%m-%d"))
                    out_dialist.append(diagnoses.values[dd])
                    out_diagrouplist.append(diagroup)
                    out_groupnamelist.append(diaNames[dd])
                    out_SOR1.append(pt_SOR1.values[dd])
                    out_SOR2.append(pt_SOR2.values[dd])
                    out_SOR3.append(pt_SOR3.values[dd])
                    out_SOR4.append(pt_SOR4.values[dd])

            elif Ngroups == 0:
                print(' - ADVARSEL: Diagnoser for '+str(cprno)+
                      ' der ikke falder i en af de definerede grupper: '+str(diaGroups))

    if verbose: print('\n')
    if verbose: print('\n - Færdig med patientkontakttjek ')
    print(' - Fandt ' + str(Nmultimorbid) + ' multisyge patienter i data')
    multimorbidratio = float(Nmultimorbid)/float(NuniqueCPR)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    logfile = outputfile.replace('.xlsx', '_log.txt')
    fout = open(logfile,'w')
    print('\n - Skriver log med resultat til '+logfile)
    fout.write('Resultat fra optælling af multisyge patienter vha. multisyge_count.py den '+todaystring+'\n')
    fout.write('\n\nInput data læst fra filen: ' + excelfile + '\n')
    fout.write('Antal multisyge og deres diagnosegrupper skrevet til: ' + outputfile + '\n\n')
    fout.write('Antal unikke CPR i input : ' + str(NuniqueCPR) + '\n')
    fout.write('Antal multisyge heraf    : ' + str(Nmultimorbid) + '\n')
    fout.write('Svarende til en andel på : ' + str("%.2f" % multimorbidratio) + '\n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print('\n - Gemmer output dataframe til Excelfil')
    outputdata = {}
    outputdata['Kontakt startdato']     = np.asarray(out_ktkstart)
    outputdata['CPR nummer']            = np.asarray(out_cprlist)
    outputdata['Diagnose']              = np.asarray(out_dialist)
    outputdata['Diagnosegruppe']        = np.asarray(out_diagrouplist)
    outputdata['Diagnosegruppe navn']   = np.asarray(out_groupnamelist)
    outputdata['Org_Ansv_SOROrgEnhedNiv1Tekst'] = np.asarray(out_SOR1)
    outputdata['Org_Ansv_SOROrgEnhedNiv2Tekst'] = np.asarray(out_SOR2)
    outputdata['Org_Ansv_SOROrgEnhedNiv3Tekst'] = np.asarray(out_SOR3)
    outputdata['Org_Ansv_SOROrgEnhedNiv4Tekst'] = np.asarray(out_SOR4)

    df_output   = pd.DataFrame(outputdata)

    df_output.to_excel(outputfile, sheet_name="data output")
    print(' - Output gemt i filen '+outputfile)
#=======================================================================================================================
def diaggroup(diagnose):
    """
    Cronical diagnosis group. Following Schiøtz et al. BMC Public Health (2017) 17:422
    """
    # Diabetes
    if diagnose[:4] in ('DE10','DE11','DE12','DE13','DE14'):
        diagroup = [1,'Diabetes']

    # Cancer
    elif (diagnose[:3] in ('DC5','DC6','DC7','DC8','DC0','DC1','DC2','DC3')) or \
         (diagnose[:4] in ('DC40','DC41','DC42','DC43','DC45','DC46','DC47','DC48','DC49',
                           'DC90','DC91','DC92','DC93','DC94','DC95','DC96','DC97')):
        diagroup = [2,'Cancer']

    # Back pain
    elif (diagnose[:3] in ('DM4')) or (diagnose[:4] in ('DM50','DM51','DM52','DM53','DM54')):
        diagroup = [3,'Back pain']

    # Osteoarthritis
    elif diagnose[:4] in ('DM15','DM16','DM17','DM18','DM19'):
        diagroup = [4,'Osteoarthritis']

    # Osteoporosis
    elif diagnose[:4] in ('DM80','DM81','DM82'):
        diagroup = [5,'Osteoporosis']

    # Joint disease
    elif (diagnose[:4] in ('DM05')) or (diagnose[:5] in ('DM060','DM068','DM070','DM071','DM073','DM100','DM109')):
        diagroup = [6,'Joint disease']

    # Allergies
    elif (diagnose[:5] in ('DJ301','DJ302','DJ303','DJ304')):
        diagroup = [7,'Allergies']

    # Chronic obstructive pulmonary disease (COPD)
    elif (diagnose[:4] in ('DJ40','DJ41','DJ42','DJ43','DJ44','DJ47','DJ96')):
        diagroup = [8,'Chronic obstructive pulmonary disease']

    elif (diagnose[:4] in ('DJ40','DJ41','DJ42','DJ43','DJ44','DJ47','DJ96')):
        diagroup = [9,'Chronic obstructive pulmonary disease']

    # Dementia
    elif (diagnose[:4] in ('DF00', 'DG30', 'DF01')) or \
            (diagnose[:5] in ('DF020','DF039','DG319')) or \
            (diagnose[:6] in ('DG318B','DG318E','DG310B')):
        diagroup = [10,'Dementia']

    # Schizophrenia
    elif (diagnose[:4] in ('DF20','DF21','DF22','DF25','DF28','DF29','DF31')):
        diagroup = [11,'Schizophrenia']

    # Anxiety
    elif (diagnose[:5] in ('DF401','DF411')):
        diagroup = [12,'Anxiety']

    # High cholesterol
    elif (diagnose[:5] in ('DE780','DE782','DE784','DE785')):
        diagroup = [13,'High cholesterol']

    # Hypertension
    elif (diagnose[:4] in ('DI10','DI11','DI12','DI13','DI15')):
        diagroup = [14,'Hypertension']

    # Stroke
    elif (diagnose[:4] in ('DG45','DG46','DI60','DI61','DI62','DI63','DI64','DI65','DI66','DI67','DI68','DI69')):
        diagroup = [15,'Stroke']

    # Heart disease
    elif (diagnose[:4] in ('DI20','DI21','DI23','DI24','DI25','DI50','DI11','DI13')):
        diagroup = [16,'Heart disease']

    # Default "no group"
    else:
        diagroup = [-99,'No group']

    return diagroup
#=======================================================================================================================
