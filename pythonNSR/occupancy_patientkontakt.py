import pdb
import numpy as np
import pandas as pd
import datetime as dt

import occupancy_patientkontakt as op

#------------------------------------------------------------------------------------
def load_patientkontakter(filepath,verbose=True):
    """
    filepath                stinavn til fil indeholdende afsnit, datotid_start og datotid_slut for patientkontakter

    --- EXAMPLE OF USE ---
    import occupancy_patientkontakt as op
    filepath = "O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Helle/Lungemed/Opgørelse ILS og sarkoidose august 2024.xlsx"
    df_ptktk = load_patientkontakter(filepath,verbose=True)

    """
    if verbose: print(" - Indlæser excel fil med patientkontakter: " + filepath)
    df_patientkontakter = pd.read_excel(filepath)
    return df_patientkontakter
#------------------------------------------------------------------------------------
def save_file(filepath,verbose=True):
    """

    --- EXAMPLE OF USE ---
    import occupancy_patientkontakt as op
    filepath = "O:\xx\yy\filnavn.xlsx"
    df_ptktk = load_patientkontakter(filepath,verbose=True)

    """
    if verbose: print(" - Gemmer output til fil: " + filepath)
    df_patientkontakter = pd.to_excel(filepath)
    return df_patientkontakter
#------------------------------------------------------------------------------------
def calcfunction(filepath,verbose=True):

    """
    filepath                stinavn til fil indeholdende afsnit, datotid_start og datotid_slut for patientkontakter

    --- EXAMPLE OF USE ---
    import occupancy_patientkontakt as op
    filepath = "O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Helle/Lungemed/Opgørelse ILS og sarkoidose august 2024.xlsx"
    op.calcfunction(filepath,verbose=true)
    """
    df_ptktk = load_patientkontakter(filepath, verbose=True)

    Afsnitsliste = np.unique(df_ptktk[''])
    Nafsnit = len(Afsnitsliste)

    # Define the start and end dates
    start_date = '2024-01-01 00:00:00'
    end_date = '2024-03-31 23:59:00'

    # Generate a date range
    date_vector = pd.date_range(start=start_date, end=end_date, freq='H')

    # Extract day, month, and year into separate columns
    day = date_vector.day
    month = date_vector.month
    year = date_vector.year
    hour = date_vector.hour
    len_datevector = len(date_vector)

    outarray = np.zeros((len_datevector*Nafsnit,5))

    if verbose: print(' - Fylder år, måned, dag, time og afsnit ind i output array')
    for aa, afsnit in enumerate(Afsnitsliste):
        if verbose: print(' - Kigger på afsnit: '+str(aa)+' '+afsnit)
        outarray[aa*len_datevector:(aa+1)*len_datevector, 0] = year
        outarray[aa*len_datevector:(aa+1)*len_datevector, 1] = month
        outarray[aa*len_datevector:(aa+1)*len_datevector, 2] = day
        outarray[aa*len_datevector:(aa+1)*len_datevector, 3] = hour
        outarray[aa * len_datevector:(aa + 1) * len_datevector, 4] = aa

    if verbose: print(' - Fylder outarray med optælling af patienter på afsnit for givne år, måned, dag, timer og afsnit')


    for aa, afsnit in enumerate(Afsnitsliste):
        ent_sel = np.where(patientkontakt_afsnit == afsnit)[0]
        for pp, patientkontakt in enumerate(patientkontaktlise[ent_sel]):
            outvec = np.zeros(len_datevector)
            ent_indlagt = np.where(patientkontakt_datostart >= date_vector & patientkontakt_datoslut <= date_vector)[0]
            outvec[ent_indlagt] = 1

            outarray[aa * len_datevector:(aa + 1) * len_datevector, 5] = outarray[aa * len_datevector:(aa + 1) * len_datevector, 5] + outvec


    #outname = filepath.split('/')[]
    #op.save_file(outname)