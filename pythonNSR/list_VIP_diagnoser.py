import numpy as np
#-----------------------------------------------------------------------------------------------------------------------
def listdiagnoser():
    """

    import list_VIP_diagnoser
    diagnoser = list_VIP_diagnoser.listdiagnoser()

    """
    #filename = "O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/datateam/NSR VIP patienter/diagnoser_220210.txt"
    filename = "O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/datateam/NSR VIP patienter/diagnoser_ringsted_220222.txt"
    diagstrings = np.genfromtxt(filename, dtype=None, comments='#', delimiter='zzz')

    diagarr_all = []
    for dstr in diagstrings:
        dstrings = dstr.split(';')
        dlist = [dd.replace("Ã…", "Å").replace("Ã¥", "å").replace("Ã˜", "Ø").replace("Ã¸", "ø").replace("Ã†", "Æ").replace("Ã¦", "æ") for dd in dstrings]
        diagarr_all = diagarr_all + dlist

    diagarr = np.unique(np.sort(np.asarray(diagarr_all)))

    print('########################DIAGNOSER##########################')
    for diag in diagarr:
        Ndiag = len(np.where(np.asarray(diagarr_all) == diag)[0])
        print(diag+'   ;   '+str(Ndiag))
    print('###########################################################')
    print(' Fandt '+str(len(diagarr))+' unikke diagnoser i \n '+filename)
    return diagarr
#-----------------------------------------------------------------------------------------------------------------------

