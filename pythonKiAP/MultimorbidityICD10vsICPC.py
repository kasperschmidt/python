"""
Loading a mapping between ICPC2 and ICD10 from an SQL database and generating af classification into
chronical deseases and multimorbidity. The classification and grouping follows Schiøtz et al. BMC Public Health (2017) 17:422.
This work is based on ICD10 codes, which is translated into a grouping in ICPC2 in the current code.
"""
import numpy as np
import MultimorbidityICD10vsICPC as mmii
import kbsKiAPutilities as kku

#--------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    verbose = True
    mmii.main(verbose=verbose)

#--------------------------------------------------------------------------------------------------------------------
def main(verbose=True):
    """
    The main wrapper to run MultimorbidityICD10vsICPC
    """
        
    print(" - Defining multimorbidity classes for ICPC based on ICD10 calssification in Schiøtz et al. BMC Public Health (2017) 17:422")

    map_ICPCICD10 = mmii.load_ICPC_ICD10_mapping(verbose=verbose)
    map_ICPCICD10 = mmii.add_multimorbitity_column(map_ICPCICD10,verbose=verbose)

    ICD10_list = ['DM15','de10']
    for ii, ICD10 in enumerate(ICD10_list):
        print(' - Looking op map and multimorbidity category for : '+ICD10)
        diagnosis_mapping, multimorbidity_category = mmii.map_ICD10toICPC(ICD10,map_ICPCICD10,ICD10in=True,verbose=verbose)
        print('   The returned mapping:     \n'+str(diagnosis_mapping))
        print('   The returned mm category: \n'+str(multimorbidity_category))

    ICPC_list = ['A01','b80','T89']
    for jj, ICPC in enumerate(ICPC_list):
        print(' - Looking op map and multimorbidity category for : '+ICPC)
        diagnosis_mapping, multimorbidity_category = mmii.map_ICD10toICPC(ICPC,map_ICPCICD10,ICD10in=False,verbose=verbose)
        print('   The returned mapping:     \n'+str(diagnosis_mapping))
        print('   The returned mm category: \n'+str(multimorbidity_category))
    

#--------------------------------------------------------------------------------------------------------------------
def map_ICD10toICPC(diagnosis,map_ICPCICD10,ICD10in=True,verbose=True):
    """
    Function returning (by default) ICPC code for ICD10 diagnoses.
    If ICD10in=False the input should be an ICPC code, and the returned value, will then be the ICD10 code.

    --- INPUT ---
    diagnosis           String containing diagnosis to map
    map_ICPCICD10       Pandas dataframe with mapping between ICPC and ICD10 from load_ICPC_ICD10_mapping()
    ICD10in             Indicating whether the content of 'diagnosis' is an ICD10 coder or an ICPC code
    verbose             Toggle verbosity of function

    """
    if ICD10in:
        diagnosis = diagnosis.lower()
        if len(diagnosis) == 3:
            in_col    = 'icd10_level3'
        elif len(diagnosis) == 4:
            in_col    = 'icd10_level4'
        elif len(diagnosis) == 5:
            in_col    = 'icd10_level5'
        else:
            in_col    = 'icd10'
        out_col   = 'icpc2'
    else:
        in_col    = 'icpc2'
        out_col   = 'icd10'
        diagnosis = diagnosis.upper()
    
    goodent = np.where(map_ICPCICD10[in_col] == diagnosis)[0]
    if len(goodent) > 0:
        diagnosis_mapping       = map_ICPCICD10.loc[goodent,out_col]
        multimorbidity_category = map_ICPCICD10.loc[goodent,'multimobidity_category']
    else:
        if verbose: print(' - WARNING: There were not matches to the input diagnosis (ICD10in='+str(ICD10in)+'):'+in_col)
        diagnosis_mapping       = np.asarray([])
        multimorbidity_category = np.asarray([])
    
    return diagnosis_mapping, multimorbidity_category
#--------------------------------------------------------------------------------------------------------------------
def load_ICPC_ICD10_mapping(verbose=True):
    """
    Function loading overall ICPC vs ICD10 code mapping
    """
    if verbose: print(' - Loading ICPC-ICD10 mapping from SQL database')
    SQLquery = "SELECT * FROM opslagsdata.dbo.icpc_icd10_mapning"
    server='10.19.31.20' 
    database='opslagsdata'
    savefilename=None
    overwrite=False
    df_map = kku.returnSQLdatapull(SQLquery,server=server,database=database,savefilename=savefilename,overwrite=overwrite,verbose=verbose)

    if verbose: print(' - Aadding convenience columns to dataframe (ICD10 level 3, 4 and 5)')
    df_map['icd10_level3'] = df_map['icd10'].str.slice(0,3)
    df_map['icd10_level4'] = df_map['icd10'].str.slice(0,4)
    df_map['icd10_level5'] = df_map['icd10'].str.slice(0,5)

    return df_map
#--------------------------------------------------------------------------------------------------------------------
def add_multimorbitity_column(map_ICPCICD10,verbose=True):
    """
    Function adding column with multimorbidity category to ICPC-ICD10 mapping dataframe
    """
    if verbose: print(' - Looping over ICD10 codes to add multimorbidity classification ')
    unique_codes = np.unique(map_ICPCICD10['icd10'].dropna())

    map_ICPCICD10['multimobidity_category']        = map_ICPCICD10['id']*0.0 # dummy column to be filled containing 0's and Nones from ICD10 column
    map_ICPCICD10['multimobidity_category_text']   = map_ICPCICD10['icd10'].str.slice(0,1) # dummy column to be filled containing d's and Nones from ICD10 column
    
    for ii, ucode in enumerate(unique_codes):
        goodent = np.where(map_ICPCICD10['icd10'] == ucode)[0]
        if ucode is None:
            map_ICPCICD10.loc[goodent,'multimobidity_category']         = None
            map_ICPCICD10.loc[goodent,'multimobidity_category_text']    = None

        else:
            diagroup = mmii.multimorbidity_diagnosis_grouping_ICD10(ucode)
            map_ICPCICD10.loc[goodent,'multimobidity_category']         = diagroup[0]
            map_ICPCICD10.loc[goodent,'multimobidity_category_text']    = diagroup[1]

    return map_ICPCICD10
#--------------------------------------------------------------------------------------------------------------------
def multimorbidity_diagnosis_grouping_ICD10(ICD10code):
    """
    Cronical diagnosis groups following Tabel 1 of Schiøtz et al. BMC Public Health (2017) 17:422.
    As the diabetes group is not given ICD10 codes in Schiøtz et al. they are manually added.
    """
    ICD10code = ICD10code.upper()
    # Diabetes
    if ICD10code[:4] in ('DE10','DE11','DE12','DE13','DE14'):
        diagroup = [1,'Diabetes']

    # Cancer
    elif (ICD10code[:3] in ('DC5','DC6','DC7','DC8','DC0','DC1','DC2','DC3')) or \
         (ICD10code[:4] in ('DC40','DC41','DC42','DC43','DC45','DC46','DC47','DC48','DC49',
                           'DC90','DC91','DC92','DC93','DC94','DC95','DC96','DC97')):
        diagroup = [2,'Cancer']

    # Back pain
    elif (ICD10code[:3] in ('DM4')) or (ICD10code[:4] in ('DM50','DM51','DM52','DM53','DM54')):
        diagroup = [3,'Back pain']

    # Osteoarthritis
    elif ICD10code[:4] in ('DM15','DM16','DM17','DM18','DM19'):
        diagroup = [4,'Osteoarthritis']

    # Osteoporosis
    elif ICD10code[:4] in ('DM80','DM81','DM82'):
        diagroup = [5,'Osteoporosis']

    # Joint disease
    elif (ICD10code[:4] in ('DM05')) or (ICD10code[:5] in ('DM060','DM068','DM070','DM071','DM073','DM100','DM109')):
        diagroup = [6,'Joint disease']

    # Allergies
    elif (ICD10code[:5] in ('DJ301','DJ302','DJ303','DJ304')):
        diagroup = [7,'Allergies']

    # Chronic obstructive pulmonary disease (COPD)
    elif (ICD10code[:4] in ('DJ40','DJ41','DJ42','DJ43','DJ44','DJ47','DJ96')):
        diagroup = [8,'Chronic obstructive pulmonary disease']

    elif (ICD10code[:4] in ('DJ40','DJ41','DJ42','DJ43','DJ44','DJ47','DJ96')):
        diagroup = [9,'Chronic obstructive pulmonary disease']

    # Dementia
    elif (ICD10code[:4] in ('DF00', 'DG30', 'DF01')) or \
            (ICD10code[:5] in ('DF020','DF039','DG319')) or \
            (ICD10code[:6] in ('DG318B','DG318E','DG310B')):
        diagroup = [10,'Dementia']

    # Schizophrenia
    elif (ICD10code[:4] in ('DF20','DF21','DF22','DF25','DF28','DF29','DF31')):
        diagroup = [11,'Schizophrenia']

    # Anxiety
    elif (ICD10code[:5] in ('DF401','DF411')):
        diagroup = [12,'Anxiety']

    # High cholesterol
    elif (ICD10code[:5] in ('DE780','DE782','DE784','DE785')):
        diagroup = [13,'High cholesterol']

    # Hypertension
    elif (ICD10code[:4] in ('DI10','DI11','DI12','DI13','DI15')):
        diagroup = [14,'Hypertension']

    # Stroke
    elif (ICD10code[:4] in ('DG45','DG46','DI60','DI61','DI62','DI63','DI64','DI65','DI66','DI67','DI68','DI69')):
        diagroup = [15,'Stroke']

    # Heart disease
    elif (ICD10code[:4] in ('DI20','DI21','DI23','DI24','DI25','DI50','DI11','DI13')):
        diagroup = [16,'Heart disease']

    # Default "no group"
    else:
        diagroup = [-99,'No group']

    return diagroup
#--------------------------------------------------------------------------------------------------------------------