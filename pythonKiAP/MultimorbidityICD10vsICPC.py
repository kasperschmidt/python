import numpy as np


#--------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print(" - Defining multimorbidity classes for ICPC based on ICD10 calssification in Schiøtz et al. BMC Public Health (2017) 17:422")


#--------------------------------------------------------------------------------------------------------------------
def map_ICD10toICPC(diagnosis,ICD10in=True):
    """
    Function returning (by default) ICPC code for ICD10 diagnoses.
    If ICD10in=False the input should be an ICPC code, and the returned value, will then be the ICD10 code.

    """


#--------------------------------------------------------------------------------------------------------------------
def load_ICPC_ICD10_mapping():
    """
    Function loading overall ICPC vs ICD10 code mapping
    """

#--------------------------------------------------------------------------------------------------------------------
def multimorbidity_diagnosis_grouping_ICD10(ICD10code):
    """
    Cronical diagnosis groups following Schiøtz et al. BMC Public Health (2017) 17:422
    """
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
#=======================================================================================================================