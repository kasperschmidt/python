import pdb

import numpy as np
import pandas as pd
import loadMDCgroups
#-----------------------------------------------------------------------------------------------------------------------
def load_into_dataframe(year='2022', verbose=True):
    """
    Loading MDC groups into pandas dataframe

    import loadMDCgroups
    mdcgroups = loadMDCgroups.load_into_dataframe()
    """
    pathname = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/datateam/'
    filename = 'MDCgrupperinger'+year+'.xlsx'
    if verbose: print(' - loading dataframe from '+filename)
    df_mdc   = pd.read_excel(pathname + filename, sheet_name='MDCgrupper')
    return df_mdc
#-----------------------------------------------------------------------------------------------------------------------
def get_group_indices(mdcgroups, verbose=True):
    """
    Splitting group list up into dictionary and returning some stats and group indices.

    import loadMDCgroups
    mdcgroups    = loadMDCgroups.load_into_dataframe()
    groupindices = loadMDCgroups.get_group_indices(mdcgroups)

    """
    groupindices = {}

    print('Found the following MDC groups in dataframe')
    for groupno in np.arange(1, 27, 1):
        groupindex = np.where(mdcgroups['MDCgruppe'] == groupno)[0]
        groupindices['group'+str(groupno)] = groupindex
        if verbose:
            print('   Gruppe '+str("%3.f" % groupno)+': '+
                  str("%50s" % mdcgroups['MDCbeskrivelse'][groupindex].values[0])+
                  '    inkluderer '+str("%5.f" % len(groupindex))+' diagnoser')

    return groupindices
#-----------------------------------------------------------------------------------------------------------------------
