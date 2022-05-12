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

def aldersfordeling():
    """

    import list_VIP_diagnoser
    diagnoser = list_VIP_diagnoser.listdiagnoser()

    """
    vals = np.array([84.0, 74.5, 58.4, 49.13, 39.7, 58.6, 48.9, 83.9, 70.4, 4.0, 82.4, 65.5, 68.0, 67.0, 91.5, 71.0, 89.3, 88.4,
            75.8, 79.5, 59.0, 58.1, 58.1, 44.2, 61.3, 86.7, 58.7, 49.0, 42.5, 19.0, 75.0, 63.9, 59.5, 56.9, 0.3, 83.8,
            78.1, 88.5, 75.9, 73.3, 54.3, 76.3, 89.0, 82.3, 67.9, 47.6, 79.7, 75.4, 70.0, 74.3, 68.0, 73.0, 72.5, 46.0,
            74.2, 70.1, 68.5, 63.0, 48.8, 80.2, 80.1, 75.4, 72.9, 23.2, 54.8, 64.0, 94.4, 85.4, 74.0, 78.9, 70.5, 93.8,
            69.0, 66.5, 83.1, 73.6, 67.1, 77.1, 76.9, 75.1, 1.3, 78.4, 37.2, 73.0, 70.0, 55.0, 71.8, 84.8, 46.5, 43.6,
            25.9, 89.1, 77.1, 78.5, 69.0, 78.5, 74.9, 78.6, 77.8, 64.0, 83.3, 71.7, 51.2, 85.2, 76.8, 75.2, 75.0, 77.4,
            78.7, 76.0, 68.9, 75.7, 61.0, 70.7, 2.1, 68.0, 78.8, 62.8, 55.0, 77.7, 70.1, 69.5, 72.4, 78.3, 44.9, 91.4,
            89.5, 70.6, 71.0, 28.9, 92.8, 87.7, 82.0, 76.0, 63.8, 76.3, 63.7, 26.9, 53.4, 79.0, 84.6, 78.0, 76.4, 55.7,
            86.5, 87.2, 75.8, 72.1, 58.2, 47.2, 84.5, 79.2, 77.0, 50.4, 87.7, 79.8, 78.3, 68.7, 87.6, 49.8, 76.0, 89.0,
            50.8, 91.6, 54.3, 75.0, 82.8, 66.5, 49.4, 91.8, 84.9, 74.9, 45.3, 29.5, 75.1, 24.3, 88.1, 23.1, 75.2, 68.8,
            65.6, 58.0, 48.7, 77.0, 64.9, 76.7, 67.9, 57.5, 85.0, 60.5, 45.4, 55.5, 53.0, 0.0, 75.0, 55.7, 33.0, 76.6,
            60.4, 52.0, 80.7, 76.3, 67.9, 68.0, 86.1, 24.5, 97.7, 84.6, 83.6, 88.9, 78.4, 77.5, 64.0, 30.2, 0.3, 83.9,
            73.7, 60.6, 53.1, 87.0, 66.0, 45.0, 66.5, 56.1, 51.2, 72.0, 91.0, 50.0, 80.0, 72.7, 67.8, 67.8, 57.0, 29.4,
            85.3, 43.1, 80.3, 78.2, 79.7, 73.2, 26.4, 85.1, 83.4, 78.3, 77.5, 38.1, 82.4, 71.6, 17.2, 83.3, 72.0, 59.0,
            77.7, 4.0, 94.8, 78.1, 78.9, 81.7, 76.4, 47.1, 82.4, 64.3, 44.0, 81.0, 75.2, 74.1, 57.5, 67.3, 88.3, 82.5,
            67.5, 79.8, 68.9, 79.7, 53.5, 50.6, 80.0, 78.0, 72.1, 73.3, 66.5, 73.5, 47.6, 69.7, 85.0, 75.6, 72.7, 61.1,
            86.5, 59.1, 59.3, 57.5, 78.8, 73.0, 67.7, 65.9, 56.9, 77.5, 29.7, 76.3, 66.0, 41.2, 87.2, 77.6, 62.9, 22.4,
            91.5, 86.4, 81.0, 65.4, 46.9, 44.3, 82.9, 70.0, 50.0, 23.0, 81.2, 65.3, 46.5, 83.7, 77.3, 72.0, 93.1, 73.6,
            54.6, 77.0, 67.8, 68.8, 21.7, 81.4, 74.0, 73.2, 63.7, 91.8, 88.0, 79.8, 73.0, 71.8, 67.3, 79.5, 75.0, 80.7,
            55.3, 51.4, 81.4, 75.0, 36.3, 69.8, 52.6, 83.5, 71.9, 28.3, 82.2, 81.3, 79.9, 79.0, 78.3, 73.0, 72.7, 60.4,
            0.5, 95.2, 73.2, 54.0, 81.6, 78.0, 17.9, 92.5, 80.4, 84.4, 74.2, 38.9, 75.0, 75.5, 77.5, 74.5, 80.3, 69.9,
            73.4, 54.6, 52.5, 62.6, 60.0, 53.0, 84.4, 69.0, 73.4, 74.5, 77.5, 86.9, 83.7, 79.4, 53.9, 68.1, 67.9, 0.4,
            46.0, 75.8, 63.0, 52.8, 85.6, 78.1, 83.0, 69.4, 58.6, 87.5, 77.3, 0.2, 79.5, 63.1, 27.2, 23.3, 85.1, 51.2,
            32.1, 84.7, 83.0, 74.5, 62.9, 73.4, 77.8, 64.1, 54.0, 89.2, 76.7, 74.9, 56.8, 84.8, 80.9, 69.8, 46.5, 82.8,
            82.2, 57.8, 88.5, 75.0, 73.3, 67.9, 64.2, 77.8, 78.5, 80.6, 54.0, 29.8, 84.3, 83.4, 77.0, 21.8, 83.2, 78.0,
            75.1, 72.2, 78.4, 78.2, 48.6, 2.4, 86.1, 86.6, 66.3, 54.6, 62.0, 60.6, 79.0, 78.5, 38.8, 19.1, 67.4, 64.4,
            63.0, 86.9, 52.5, 67.1, 80.5, 86.6, 78.0, 77.4, 73.3, 65.1, 1.6, 81.4, 71.0, 64.0, 78.0, 60.8, 60.6, 78.6,
            70.9, 67.0, 30.1, 81.0, 76.0, 69.3, 25.1, 87.9, 56.0, 72.0, 69.8, 88.7, 42.0, 75.2, 69.7, 49.0, 56.0, 80.4,
            69.1, 69.0, 66.2, 58.0, 76.0, 66.0, 80.7, 78.7, 91.0, 45.5, 26.2, 8.1, 51.1, 64.5, 79.8, 75.7, 75.2, 71.4,
            81.3, 79.4, 83.2, 72.3, 4.5, 2.2, 1.5, 84.3, 60.5, 52.1, 85.7, 76.9, 85.6, 80.4, 41.0, 29.4, 75.7, 61.2,
            80.6, 66.7, 68.2, 89.2, 65.5, 84.3, 82.7, 68.9, 48.2, 26.2, 73.1, 80.0, 47.5, 26.9, 97.1, 76.0, 72.3, 45.6,
            37.0, 83.0, 70.0, 55.2, 21.2, 55.3, 77.1, 92.5, 69.8, 54.3, 88.5, 83.3, 65.1, 53.6, 88.1, 81.0, 71.5, 82.3,
            78.8, 70.0, 90.5, 69.0, 79.5, 72.6, 21.0, 81.7, 70.8, 61.7, 90.1, 77.6, 89.6, 67.1, 75.0, 66.6, 80.8, 62.0,
            50.0, 81.0, 79.4, 86.8, 79.8, 79.8, 76.0, 89.3, 81.8, 76.6, 65.0, 54.0, 69.0, 83.0, 71.6, 57.3, 47.9, 86.0,
            78.7, 56.0, 0.8, 65.4, 65.0, 92.0, 85.0, 92.0, 71.8, 95.9, 74.2, 73.3, 66.6, 61.8, 88.8, 76.1, 78.6, 49.8,
            79.9, 92.0, 73.6, 69.2, 72.6, 65.8, 36.9, 80.9, 75.9, 70.5, 63.5, 49.0, 63.0, 1.2, 63.9, 77.5, 75.6, 10.3,
            72.6, 65.0, 86.5, 70.6, 80.6, 90.6, 84.9, 56.1, 80.3, 75.4, 75.6, 65.3, 60.9, 56.8, 84.0, 64.3, 46.5, 88.9,
            83.8, 60.2, 85.7, 81.2, 60.5, 55.5, 39.0, 96.6, 84.6, 78.5, 67.3, 37.8, 84.2, 77.8, 61.8, 54.8, 76.5, 63.4,
            56.4, 15.9, 79.8, 0.8, 75.5, 74.8, 73.1, 80.0, 67.6, 49.8, 58.5, 54.1, 74.5, 54.0, 89.0, 85.0, 60.6, 56.4,
            91.4, 75.4, 73.4, 73.5, 29.3, 93.5, 88.4, 86.8, 79.1, 69.7, 25.8, 19.3, 75.0, 52.9, 49.6, 91.8, 81.5, 75.7,
            72.0, 67.5, 46.7, 86.3, 83.6, 79.5, 54.1, 88.0, 87.8, 86.8, 77.0, 39.7, 82.9, 74.0, 84.4, 60.0, 54.3, 78.0,
            0.9, 76.3, 56.7, 23.0, 55.1, 88.0, 88.0, 58.5, 53.0, 80.6, 69.4, 66.0, 87.0, 34.1, 81.5, 59.4, 48.8, 83.4,
            82.4, 80.6, 74.7, 73.7, 48.1, 7.6, 76.2, 59.3, 56.8, 83.0, 77.1, 62.6, 83.0, 74.8, 43.2, 40.1, 66.6, 64.3,
            30.4, 27.7, 75.8, 72.4, 72.5, 82.9, 80.9, 78.1, 65.9, 76.9, 68.0, 84.4, 23.0, 86.6, 53.0, 44.1, 80.0, 64.3,
            38.1, 76.6, 75.8, 88.4, 75.5, 90.5, 86.8, 77.7, 78.4, 71.6, 53.0, 3.0, 81.0, 77.3, 46.5, 21.3, 78.0, 89.0,
            84.6, 33.3, 87.0, 81.0, 53.1, 82.6, 47.4, 80.7, 0.0, 86.0, 82.9, 72.7, 49.1, 33.2, 73.5, 73.0, 91.9, 86.0,
            79.0, 76.6, 74.4, 72.7, 60.4, 76.7])


= If [Akut PI]        = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [Anæ PI]       = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [Anæ NAE PI]   = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [Endo PI]      = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [Ger PI]       = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [GynObs PI]    = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [Kar PI]       = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [Kir PI]       = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [Lung PI]      = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [MedGas PI]    = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [Multi PI]     = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [Neu PI]       = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [Neu NAE PI]   = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [O-kir PI]     = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [O-kir NAE PI] = 1 Then [Kontaktansvar Overafdeling navn]
ElseIf [Pæd PI]       = 1 Then [Kontaktansvar Overafdeling navn]

#-----------------------------------------------------------------------------------------------------------------------


