# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Functions and scripts to generate photoionization "PDFs" in comparison to a set of observational emission
# line constraints. Based on initial "PDF" functions in
# /Users/kschmidt/work/catalogs/BPASSbasedNebularEmission /Users/kschmidt/work/GitHub/python/NEOGALmodels.py
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import numpy as np
import glob
import astropy.io.fits as afits
import sys
import pdb
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import matplotlib
import itertools
import literaturecollection_emissionlinestrengths as lce
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter
import NEOGALmodels as nm
import BPASSmodels as bm
import photoionizationPDFs as pp
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def estimate_object_PDFs(fluxratiodictionarylist,generatePDFplots=False,maxPDFyscale=False,basename='photoionizationmodelPDFs',
                         col_NEOGAL_AGN='blue',col_NEOGAL_SF='red',col_BPASS_bin='green',col_BPASS_sin='orange',verbose=True):
    """
    Function to estimate "PDFs" from the SF and AGN NEOGAL models given a set of flux ratio measurements.

    --- INPUT ---
    fluxratiodictionarylist   List of flux ratio dictionaris to get "PDFs" for.
                              Length of list corresponds to objects with measurements to get "PDFs" for.
    AGNcol                    Color of AGN model related data in plots
    SFcol                     Color of SF model related data in plots
    verbose                   Toggle verbosity


    --- EXAMPLE OF USE ---
    import photoionizationPDFs as pp
    basename= '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/photoionizationPDFs/testing/photoionizationmodelPDFs'

    FRdic = [{'id':111111111111, 'HeII1640/OIII1663':[0.04,0.45],'CIII1908/CIV1550':[1.0,10.0]}, {'id':222222222222, 'OIII1663/HeII1640':[1e-1,1.0],'CIII1908/CIV1550':[0.1,10.0]}, {'id':333333, 'OIII1663/HeII1640':[1e2,1e3],'CIII1908/CIV1550':[1e-3,1e-2]}, {'id':444444, 'OIII1663/HeII1640':[1e-2,1e-1],'CIII1908/CIV1550':[5e-1,1e-0]}, {'id':555555, 'OIII1663/HeII1640':[1e-2,1e35],'CIII1908/CIV1550':[5e-1,1e-0], 'OIII1663/CIII1908':[0.,10.0], 'OIII1663/CIV1550':[1e-2,1e35], 'OIII1663/SiIII1888':[1e-2,1e1], 'CIII1908/SiIII1888':[1e-2,1e35], 'CIV1550/SiIII1888':[1e-2,1e35]}]
    paramcollections, collectionstats = pp.estimate_object_PDFs(FRdic, basename=basename, generatePDFplots=True, maxPDFyscale=True)

    FRdicNC = [{'id':99}] # run for a single objects with no constraints to get instrinsic distribution
    paramcollections, collectionstats = pp.estimate_object_PDFs(FRdicNC, basename=basename+'NOobsCONSTRAINTS', generatePDFplots=True, maxPDFyscale=True)

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading NEOGAL models ')
    NEOGAL_SF           = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/nebular_emission/')
    NEOGAL_AGN          = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/AGN_NLR_nebular_feltre16/')
    BPASS_bin, filename = bm.load_model('combined',filepath='dummy',binaries=True)
    BPASS_sin, filename = bm.load_model('combined',filepath='dummy',binaries=False)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Define all possible line ratios from the lines:\n    '
                      'NV1240, CIV1550, CIII1908, HeII1640, OIII1663, and SiIII1888')
    fluxratiodic = {}
    fluxratiodic['NV1240/CIV1550']     =   [0,1e35]
    fluxratiodic['NV1240/CIII1908']    =   [0,1e35]
    fluxratiodic['NV1240/HeII1640']    =   [0,1e35]
    fluxratiodic['NV1240/OIII1663']    =   [0,1e35]
    fluxratiodic['NV1240/SiIII1888']   =   [0,1e35]

    fluxratiodic['CIV1550/NV1240']     =   [0,1e35]
    fluxratiodic['CIV1550/CIII1908']   =   [0,1e35]
    fluxratiodic['CIV1550/HeII1640']   =   [0,1e35]
    fluxratiodic['CIV1550/OIII1663']   =   [0,1e35]
    fluxratiodic['CIV1550/SiIII1888']  =   [0,1e35]

    fluxratiodic['CIII1908/NV1240']     =  [0,1e35]
    fluxratiodic['CIII1908/CIV1550']    =  [0,1e35]
    fluxratiodic['CIII1908/HeII1640']   =  [0,1e35]
    fluxratiodic['CIII1908/OIII1663']   =  [0,1e35]
    fluxratiodic['CIII1908/SiIII1888']  =  [0,1e35]

    fluxratiodic['HeII1640/NV1240']     =  [0,1e35]
    fluxratiodic['HeII1640/CIV1550']    =  [0,1e35]
    fluxratiodic['HeII1640/CIII1908']   =  [0,1e35]
    fluxratiodic['HeII1640/OIII1663']   =  [0,1e35]
    fluxratiodic['HeII1640/SiIII1888']  =  [0,1e35]

    fluxratiodic['OIII1663/NV1240']     =  [0,1e35]
    fluxratiodic['OIII1663/CIV1550']    =  [0,1e35]
    fluxratiodic['OIII1663/CIII1908']   =  [0,1e35]
    fluxratiodic['OIII1663/HeII1640']   =  [0,1e35]
    fluxratiodic['OIII1663/SiIII1888']  =  [0,1e35]

    fluxratiodic['SiIII1888/NV1240']    =  [0,1e35]
    fluxratiodic['SiIII1888/CIV1550']   =  [0,1e35]
    fluxratiodic['SiIII1888/CIII1908']  =  [0,1e35]
    fluxratiodic['SiIII1888/HeII1640']  =  [0,1e35]
    fluxratiodic['SiIII1888/OIII1663']  =  [0,1e35]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Set up mode line flux vectors')
    BPASSconversion = 1e-45 # Factor to convert BPASS to 1e45 erg/s like the NEOGAL models fluxes
    fluxdic = {}
    fluxdic['NV1240']    = [NEOGAL_SF['NV1240'],
                            NEOGAL_AGN['NV1240'],
                            None,
                            None]

    fluxdic['CIV1550']   = [NEOGAL_SF['CIV1548']+NEOGAL_SF['CIV1551'],
                            NEOGAL_AGN['CIV1548']+NEOGAL_AGN['CIV1551'],
                            BPASS_bin['CIV1548']*BPASSconversion+BPASS_bin['CIV1551']*BPASSconversion,
                            BPASS_sin['CIV1548']*BPASSconversion+BPASS_sin['CIV1551']*BPASSconversion]

    fluxdic['CIII1908']  = [NEOGAL_SF['CIII1908'],
                            NEOGAL_AGN['CIII1907']+NEOGAL_AGN['CIII1910'],
                            BPASS_bin['CIII1907']*BPASSconversion+BPASS_bin['CIII1910']*BPASSconversion,
                            BPASS_sin['CIII1907']*BPASSconversion+BPASS_sin['CIII1910']*BPASSconversion]

    fluxdic['HeII1640']  = [NEOGAL_SF['HeII1640'],
                            NEOGAL_AGN['HeII1640'],
                            BPASS_bin['HeII1640']*BPASSconversion,
                            BPASS_sin['HeII1640']*BPASSconversion]

    fluxdic['OIII1663']  = [NEOGAL_SF['OIII1661']+NEOGAL_SF['OIII1666'],
                            NEOGAL_AGN['OIII1661']+NEOGAL_AGN['OIII1666'],
                            BPASS_bin['OIII1661']*BPASSconversion+BPASS_bin['OIII1666']*BPASSconversion,
                            BPASS_sin['OIII1661']*BPASSconversion+BPASS_sin['OIII1666']*BPASSconversion]

    fluxdic['SiIII1888'] = [NEOGAL_SF['SiIII1888'],
                            NEOGAL_AGN['SiIII1888'],
                            None,
                            None]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Get ranges of model flux ratios')
    Nobj = len(fluxratiodictionarylist)
    if verbose: print(' - Get model selection given flux ratio ranges according to '+
                      str(Nobj)+" object's data provided ")
    if verbose: print('   Selection based on the total number of photoionization models: \n '
                      '   NEOGAL SF      = '+str(len(NEOGAL_SF))+'\n '+
                      '   NEOGAL AGN     = '+str(len(NEOGAL_AGN))+'\n '+
                      '   BPASS binaries = '+str(len(BPASS_bin))+'\n '+
                      '   BPASS singular = '+str(len(BPASS_sin)))

    parametercollection_SF   = [{'id':0, 'Zgas':[],'logUs':[],'xid':[],'nh':[],'COCOsol':[],'mup':[]}]*Nobj
    parametercollection_AGN  = [{'id':0, 'Zgas':[],'logUs':[],'xid':[],'nh':[],'alpha':[]}]*Nobj
    parametercollection_bin  = [{'id':0, 'Zgas':[],'logUs':[],'nh':[],'logAge':[]}]*Nobj
    parametercollection_sin  = [{'id':0, 'Zgas':[],'logUs':[],'nh':[],'logAge':[]}]*Nobj

    for oo, FRdic_input in enumerate(fluxratiodictionarylist):
        objid = FRdic_input['id']

        # ------ resetting flux ratio dictionary for object ------
        fluxratioconstrains_obj = {}
        for key in fluxratiodic.keys():
            fluxratioconstrains_obj[key] = fluxratiodic[key]
        # --------------------------------------------------------

        for FR in FRdic_input.keys():
            if FR in fluxratiodic.keys():
                fluxratioconstrains_obj[FR] = FRdic_input[FR]
                # print(str(objid)+':'+FR+'  -->'+str(fluxratiodic_obj[FR]))
            elif FR == 'id':
                pass
            else:
                print(' WARNING pp.estimate_object_PDFs(): The flux ratio entry '+FR+' is not availble in the \n'
                      '                                    dictionary from the models. Define that flux \n'
                      '                                    ratio or correct input data.')

        goodent_SF   = np.arange(len(NEOGAL_SF))
        goodent_AGN  = np.arange(len(NEOGAL_AGN))
        goodent_bin  = np.arange(len(BPASS_bin))
        goodent_sin  = np.arange(len(BPASS_sin))
        for FR in fluxratiodic.keys():
            numerator   = FR.split('/')[0]
            denominator = FR.split('/')[1]

            goodent_FR_SF  = np.where( (fluxdic[numerator][0]/fluxdic[denominator][0] >= fluxratioconstrains_obj[FR][0]) &
                                       (fluxdic[numerator][0]/fluxdic[denominator][0] <= fluxratioconstrains_obj[FR][1]))[0]
            goodent_SF = np.intersect1d(goodent_SF,goodent_FR_SF)

            goodent_FR_AGN = np.where( (fluxdic[numerator][1]/fluxdic[denominator][1] >= fluxratioconstrains_obj[FR][0]) &
                                       (fluxdic[numerator][1]/fluxdic[denominator][1] <= fluxratioconstrains_obj[FR][1]))[0]
            goodent_AGN    = np.intersect1d(goodent_AGN,goodent_FR_AGN)


            if (fluxdic[numerator][2] is not None) & (fluxdic[denominator][2] is not None):
                zerovals = np.where(fluxdic[denominator][2] == 0.0)[0]
                fluxdic[denominator][2][zerovals] = 1e-25

                goodent_FR_bin = np.where( (fluxdic[numerator][2]/fluxdic[denominator][2] >= fluxratioconstrains_obj[FR][0]) &
                                           (fluxdic[numerator][2]/fluxdic[denominator][2] <= fluxratioconstrains_obj[FR][1]))[0]
                goodent_bin    = np.intersect1d(goodent_bin,goodent_FR_bin)

            else:
                if verbose & (Nobj < 10): print(' WARNING: No constraints on flux ratio '+FR+' for BPASS binary models')

            if (fluxdic[numerator][2] is not None) & (fluxdic[denominator][2] is not None):
                zerovals = np.where(fluxdic[denominator][3] == 0.0)[0]
                fluxdic[denominator][3][zerovals] = 1e-25

                goodent_FR_sin = np.where( (fluxdic[numerator][3]/fluxdic[denominator][3] >= fluxratioconstrains_obj[FR][0]) &
                                           (fluxdic[numerator][3]/fluxdic[denominator][3] <= fluxratioconstrains_obj[FR][1]))[0]

                # if len(goodent_FR_sin) < 40131:
                #     setdiff         = np.setdiff1d(goodent_sin,goodent_FR_sin)
                #     print(setdiff)
                #     setdiff_ratios  = fluxdic[numerator][3][setdiff]/fluxdic[denominator][3][setdiff]
                #     print(setdiff_ratios)
                #     pdb.set_trace()

                goodent_sin    = np.intersect1d(goodent_sin,goodent_FR_sin)
            else:
                if verbose & (Nobj < 10): print(' WARNING: No constraints on flux ratio '+FR+' for BPASS singular models')

        parametercollection_SF[oo]  = {'id'     : FRdic_input['id'],
                                       'Zgas'   : NEOGAL_SF['Zgas'][goodent_SF],
                                       'logUs'  : NEOGAL_SF['logUs'][goodent_SF],
                                       'xid'    : NEOGAL_SF['xid'][goodent_SF],
                                       'nh'     : NEOGAL_SF['nh'][goodent_SF],
                                       'COCOsol': NEOGAL_SF['COCOsol'][goodent_SF],
                                       'mup'    : NEOGAL_SF['mup'][goodent_SF],
                                       'alpha'  : None,
                                       'logAge' : None}

        parametercollection_AGN[oo] = {'id'     : FRdic_input['id'],
                                       'Zgas'   : NEOGAL_AGN['Zgas'][goodent_AGN],
                                       'logUs'  : NEOGAL_AGN['logUs'][goodent_AGN],
                                       'xid'    : NEOGAL_AGN['xid'][goodent_AGN],
                                       'nh'     : NEOGAL_AGN['nh'][goodent_AGN],
                                       'COCOsol': None,
                                       'mup'    : None,
                                       'alpha'  : NEOGAL_AGN['alpha'][goodent_AGN],
                                       'logAge' : None}

        parametercollection_bin[oo] = {'id'     : FRdic_input['id'],
                                       'Zgas'   : BPASS_bin['Zgas'][goodent_bin],
                                       'logUs'  : BPASS_bin['logU'][goodent_bin],
                                       'xid'    : None,
                                       'nh'     : 10**(BPASS_bin['lognH'][goodent_bin]),
                                       'COCOsol': None,
                                       'mup'    : None,
                                       'alpha'  : None,
                                       'logAge' : BPASS_bin['logAge'][goodent_bin]}

        parametercollection_sin[oo] = {'id'     : FRdic_input['id'],
                                       'Zgas'   : BPASS_sin['Zgas'][goodent_sin],
                                       'logUs'  : BPASS_sin['logU'][goodent_sin],
                                       'xid'    : None,
                                       'nh'     : 10**(BPASS_sin['lognH'][goodent_sin]),
                                       'COCOsol': None,
                                       'mup'    : None,
                                       'alpha'  : None,
                                       'logAge' : BPASS_sin['logAge'][goodent_sin]}

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Getting distribution ranges (percentiles) for parameter collections ')
    stat_SF   = pp.get_collectionranges(parametercollection_SF,parametercollection_SF[0].keys())
    stat_AGN  = pp.get_collectionranges(parametercollection_AGN,parametercollection_AGN[0].keys())
    stat_bin  = pp.get_collectionranges(parametercollection_bin,parametercollection_bin[0].keys())
    stat_sin  = pp.get_collectionranges(parametercollection_sin,parametercollection_sin[0].keys())

    paramcollections = [parametercollection_SF,parametercollection_AGN,parametercollection_bin,parametercollection_sin]
    statdics         = [stat_SF,stat_AGN,stat_bin,stat_sin]
    for ss, statdic in enumerate(statdics):
        paramcoll                  = paramcollections[ss]
        stat_idlist                =  [statdic[oo]['id'] for oo in np.arange(len(statdic))]
        parametercollection_idlist =  [paramcoll[oo]['id'] for oo in np.arange(len(paramcoll))]

        if stat_idlist != parametercollection_idlist:
            sys.exit(' photoionizationPDFs.estimate_object_PDFs(): Wait a minute... the ID list is not identical between \n'
                     '                                             the parameter collection ('+str(parametercollection_idlist)+') and \n'
                     '                                             the stats ('+str(stat_idlist)+')')

    paramcollections  = [parametercollection_SF, parametercollection_AGN, parametercollection_bin, parametercollection_sin]
    collectionstats   = [stat_SF, stat_AGN, stat_bin, stat_sin]
    collectionscolors = [col_NEOGAL_SF,col_NEOGAL_AGN,col_BPASS_bin,col_BPASS_sin]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = basename+'_Stats.pdf'
    pp.plot_stat(plotname,collectionstats,collectionscolors,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if generatePDFplots:
        if verbose: print(' - Plotting the extracted model parameter collections')
        plotname = basename+'_PDFs.pdf'
        pp.plot_modelparametercollections(plotname, paramcollections, collectionstats, collectionscolors,
                                          fluxratiodictionarylist=fluxratiodictionarylist, maxPDFyscale=maxPDFyscale,
                                          verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return paramcollections, collectionstats
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_collectionranges(parametercollection,parameterlist):
    """

    """
    statparam = []
    Nobj = len(parametercollection)
    for oo in np.arange(Nobj):
        outputdic = {}
        for paramkey in parameterlist:
            outputdic[paramkey] = []
        outputdic['id'] = parametercollection[oo]['id']
        statparam.append(outputdic)

        for key in statparam[oo].keys():
            if key == 'id': continue

            if (parametercollection[oo][key] is not None):
                if (len(parametercollection[oo][key]) > 0):
                    meanval   = np.mean(parametercollection[oo][key])
                    stdval    = np.std(parametercollection[oo][key])
                    medianval = np.median(parametercollection[oo][key])
                    perc2p5   = np.sort(parametercollection[oo][key])[int(len(parametercollection[oo][key])*0.025)]
                    perc16    = np.sort(parametercollection[oo][key])[int(len(parametercollection[oo][key])*0.16)]
                    perc25    = np.sort(parametercollection[oo][key])[int(len(parametercollection[oo][key])*0.25)]
                    perc50    = np.sort(parametercollection[oo][key])[int(len(parametercollection[oo][key])*0.50)]
                    perc75    = np.sort(parametercollection[oo][key])[int(len(parametercollection[oo][key])*0.75)]
                    perc84    = np.sort(parametercollection[oo][key])[int(len(parametercollection[oo][key])*0.84)]
                    perc97p5  = np.sort(parametercollection[oo][key])[int(len(parametercollection[oo][key])*0.975)]

                    statparam[oo][key] = [meanval,stdval,medianval,perc2p5,perc16,perc25,
                                          perc50,perc75,perc84,perc97p5]
                else:
                    statparam[oo][key] = [np.nan]*10
            else:
                statparam[oo][key] = [np.nan]*10

    return statparam
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_stat(plotname, collectionstats, collectionscolors, verbose=True):
    """

    """
    if verbose: print(' - Plot stats of parameter selection')
    Nobj     = len(collectionstats[0])
    objidall = np.array([statobj['id'] for statobj in collectionstats[0]])

    xranges = {}
    xranges['Zgas']    =  [1e-6,0.2]
    xranges['logUs']   =  [-5.5,-0.5]
    xranges['xid']     =  [0.02,0.58]
    xranges['nh']      =  [1,1e5]
    xranges['COCOsol'] =  [0,1.5]
    xranges['mup']     =  [50,350]
    xranges['alpha']   =  [-2.2,-0.9]
    xranges['logAge']  =  [5,10]

    keylists = []
    for colstat in collectionstats:
        keylists.append([key for key in colstat[0].keys()])

    uniquekeys = np.unique(keylists[0]+keylists[1]+keylists[2]+keylists[3])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for key in uniquekeys:
        if key == 'id': continue
        showmedian = False
        if showmedian:
            xval_ent    = 2 # median
            xstr        = 'median'
        else:
            xval_ent    = 0 # mean
            xstr        = 'mean'

        errval_ent  = 1 # std

        xval_SF     = np.array([])
        errval_SF   = np.array([])
        xval_AGN = np.array([])
        errval_AGN  = np.array([])
        xval_bin = np.array([])
        errval_bin  = np.array([])
        xval_sin = np.array([])
        errval_sin  = np.array([])

        range50_SF  = np.array([])
        range50_AGN = np.array([])
        range50_bin = np.array([])
        range50_sin = np.array([])

        for oo in np.arange(Nobj):
            listent = 0 # NEOGAL SF
            if key in keylists[listent]:
                xval_SF = np.append(xval_SF,collectionstats[listent][oo][key][xval_ent])
                errval_SF = np.append(errval_SF,collectionstats[listent][oo][key][1])
                range50_SF = np.append(range50_SF,collectionstats[listent][oo][key][7]-collectionstats[listent][oo][key][7])

            listent = 1 # NEOGAL AGN
            if key in keylists[listent]:
                xval_AGN = np.append(xval_AGN,collectionstats[listent][oo][key][xval_ent])
                errval_AGN = np.append(errval_AGN,collectionstats[listent][oo][key][1])
                range50_AGN = np.append(range50_AGN,collectionstats[listent][oo][key][7]-collectionstats[listent][oo][key][7])

            listent = 2 # BPASS binaries
            if key in keylists[listent]:
                xval_bin = np.append(xval_bin,collectionstats[listent][oo][key][xval_ent])
                errval_bin = np.append(errval_bin,collectionstats[listent][oo][key][1])
                range50_bin = np.append(range50_bin,collectionstats[listent][oo][key][7]-collectionstats[listent][oo][key][7])

            listent = 3 # BPASS singles
            if key in keylists[listent]:
                xval_sin = np.append(xval_sin,collectionstats[listent][oo][key][xval_ent])
                errval_sin = np.append(errval_sin,collectionstats[listent][oo][key][1])
                range50_sin = np.append(range50_sin,collectionstats[listent][oo][key][7]-collectionstats[listent][oo][key][7])

        xvalues   = [xval_SF,xval_AGN,xval_bin,xval_sin]
        errvalues = [errval_SF,errval_AGN,errval_bin,errval_sin]
        sampstr   = ['NEOGAL SF', 'NEOGAL AGN', 'BPASS binary', 'BPASS single']
        #-------- PLOT 'EM --------
        plotname_stat = plotname.replace('.pdf','_'+xstr+'VSerr_'+key+'.pdf')

        fig = plt.figure(1, figsize=(6, 6))
        # definitions for the axes
        left, width = 0.15, 0.60
        bottom, height = 0.15, 0.60
        bottom_h = left_h = left + width + 0.01

        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=left, right=left+width, bottom=bottom, top=bottom+height)
        rect_histx = [left, bottom_h, width, 0.2]
        rect_histy = [left_h, bottom, 0.2, height]

        Fsize    = 14
        lthick   = 2
        marksize = 6
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()

        for oo in np.arange(Nobj): # loop necessary for coloring and markers
            objid = objidall[oo]
            mfc   = True
            if (objid < 6e8): # CDFS and COSMOS
                markersym   = 'o'
            elif (objid < 7e8) & (objid > 6e8): # UDF
                markersym   = 'D'
            elif (objid < 9e8) & (objid > 7e8): # UDF10
                markersym   = 'X'
            elif (objid > 1e9): # Literature objects
                markersym   = lce.get_reference_fromID(objid,verbose=False)[4]
                mfc         = False
            else:
                print(' WARNING - stopped as could not assing a marker symbol to the id '+str(objid))
                pdb.set_trace()

            ms          = marksize
            limsizefrac = 0.05

            # change color of limits
            markerzorder = 20
            #--------------------------- Plot scatter plot ---------------------------
            for xx, xvalue in enumerate(xvalues):
                if mfc:
                    markerfacecolor = collectionscolors[xx]
                else:
                    markerfacecolor = 'None'

                if len(xvalue) > 0:
                    plt.errorbar(xvalue[oo],errvalues[xx][oo]/np.abs(xvalue[oo]),xerr=None,yerr=None,capthick=0.5,
                                 marker=markersym,lw=lthick/2., markersize=ms, alpha=0.5,
                                 markerfacecolor=markerfacecolor,ecolor=collectionscolors[xx],
                                 markeredgecolor=collectionscolors[xx],zorder=markerzorder)

        plt.xlim(xranges[key])
        xminsys, xmaxsys = plt.xlim()

        # plt.ylim(yranges[key])
        yminsys, ymaxsys = plt.ylim()

        plt.xlabel(xstr+'( '+pp.keylabels(key)+' )')
        plt.ylabel('Normalized standard deviation of parameter PDF')

        if (key == 'Zgas') or (key == 'nh'):
            xlog = True
            plt.xscale('log')
        else:
            xlog = False

        #--------------------------- Estimate wheighted averages ---------------------------

        for xx, xvalue in enumerate(xvalues):
            if any(np.isfinite(errvalues[xx])):
                goodents       = np.where((objidall < 1e9) & np.isfinite(xvalue) & np.isfinite(errvalues[xx]) &
                                          (np.nan_to_num(errvalues[xx]) > 0.0))[0]
                xvalues_uves   = xvalue[goodents]
                errvalues_uves = errvalues[xx][goodents]
                whtav_uves     = np.sum( xvalues_uves/errvalues_uves**2 ) / np.sum(1.0/errvalues_uves**2)
                whtav_uves_err = 1/np.sqrt(np.sum(errvalues_uves))
                if verbose: print('\n - wheighted average UVES('+sampstr[xx]+','+key+') = '+str(whtav_uves)+'+/-'+str(whtav_uves_err)+
                                  ' (Nmeasurements='+str(len(goodents))+')')

                goodents       = np.where((objidall > 1e9) & np.isfinite(xvalue) & np.isfinite(errvalues[xx]) &
                                          (np.nan_to_num(errvalues[xx]) > 0.0))[0]
                xvalues_lit    = xvalue[goodents]
                errvalues_lit  = errvalues[xx][goodents]
                whtav_lit      = np.sum( xvalues_lit/errvalues_lit**2 ) / np.sum(1.0/errvalues_lit**2)
                whtav_lit_err  = 1/np.sqrt(np.sum(errvalues_lit))
                if verbose: print(' - wheighted average literature('+sampstr[xx]+','+key+') = '+str(whtav_lit)+'+/-'+str(whtav_lit_err)+
                                  ' (Nmeasurements='+str(len(goodents))+')')

                goodents       = np.where((objidall > 1.0) & np.isfinite(xvalue) & np.isfinite(errvalues[xx]) &
                                          (np.nan_to_num(errvalues[xx]) > 0.0))[0]
                xvalues_all    = xvalue[goodents]
                errvalues_all  = errvalues[xx][goodents]
                whtav_all      = np.sum( xvalues_all/errvalues_all**2 ) / np.sum(1.0/errvalues_all**2)
                whtav_all_err  = 1/np.sqrt(np.sum(errvalues_all))
                if verbose: print(' - wheighted average all('+sampstr[xx]+','+key+') = '+str(whtav_all)+'+/-'+str(whtav_all_err)+
                                  ' (Nmeasurements='+str(len(goodents))+')'+'\n')

            else:
                if verbose: print(' - wheighted averages not estimated as no finite (error) values for ('+sampstr[xx]+','+key+') \n')

        #--------------------------- Plot X-axis Histograms ---------------------------
        Nbins   = 50
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)

        axHistx.xaxis.set_major_formatter(NullFormatter())
        axHisty.yaxis.set_major_formatter(NullFormatter())

        binwidth_x = np.diff([xminsys,xmaxsys])/Nbins
        bindefs    = np.arange(xminsys, xmaxsys+binwidth_x, binwidth_x)
        if xlog:
            bindefs = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))
            axHistx.set_xscale('log')

        for xx, xvalue in enumerate(xvalues):
            if len(xvalue) > 0:
                axHistx.hist(xvalue[objidall < 1e9][np.isfinite(xvalue[objidall < 1e9])],linestyle='-',
                             bins=bindefs,histtype='step',color=collectionscolors[xx])
                axHistx.hist(xvalue[np.isfinite(xvalue)],linestyle=':',
                             bins=bindefs,histtype='step',color=collectionscolors[xx])


        #--------------------------- Plot Y-axis Histograms ---------------------------
        axHistx.set_xticks([])
        axHistx.set_xlim([xminsys,xmaxsys])

        binwidth_y = np.diff([yminsys,ymaxsys])/Nbins
        bindefs    = np.arange(yminsys, ymaxsys+binwidth_y, binwidth_y)

        for xx, xvalue in enumerate(xvalues):
            if len(xvalue) > 0:
                histvalues  = errvalues[xx][objidall < 1e9][np.isfinite(errvalues[xx][objidall < 1e9])]/ \
                              np.abs(xvalue[objidall < 1e9][np.isfinite(errvalues[xx][objidall < 1e9])])
                axHisty.hist(histvalues,linestyle='-',bins=bindefs,histtype='step',color=collectionscolors[xx],
                             orientation='horizontal')
                histvalues  = errvalues[xx][np.isfinite(errvalues[xx])]/ \
                              np.abs(xvalue[np.isfinite(errvalues[xx])])
                axHisty.hist(histvalues,linestyle=':',bins=bindefs,histtype='step',color=collectionscolors[xx],
                             orientation='horizontal')

        axHisty.set_yticks([])
        axHisty.set_ylim([yminsys,ymaxsys])

        if verbose: print('   Saving plot to '+plotname_stat)
        plt.savefig(plotname_stat)
        plt.clf()
        plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    xkey   = 'Zgas'
    xlabel = 'Z'
    ykey   = 'logUs'
    ylabel = 'log(U)'

    xrange = xranges[xkey]
    yrange = xranges[ykey]

    plotname_obj = plotname.replace('.pdf','_meanANDstd_ZgasVSlogUs.pdf')
    pp.plot_stat_diagram(plotname_obj,collectionstats,xkey,ykey,xlabel,ylabel,collectionscolors,xrange,yrange,
                         xlog=True,ylog=False,datsetup='meanstd',verbose=verbose)

    plotname_obj = plotname.replace('.pdf','_median68percent_ZgasVSlogUs.pdf')
    pp.plot_stat_diagram(plotname_obj,collectionstats,xkey,ykey,xlabel,ylabel,collectionscolors,xrange,yrange,
                         xlog=True,ylog=False,datsetup='med68',verbose=verbose)

    plotname_obj = plotname.replace('.pdf','_medianANDstd_ZgasVSlogUs.pdf')
    pp.plot_stat_diagram(plotname_obj,collectionstats,xkey,ykey,xlabel,ylabel,collectionscolors,xrange,yrange,
                         xlog=True,ylog=False,datsetup='medstd',verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_stat_diagram(plotname,collectionstats,xkey,ykey,xlabel,ylabel,collectionscolors,xrange,yrange,
                      xlog=False,ylog=False,datsetup='med68',verbose=True):

    objidall = np.array([statobj['id'] for statobj in collectionstats[0]])

    if datsetup == 'med68':
        cen_ent     = 2
        errlow_ent  = 4
        errhigh_ent = 8
    elif datsetup == 'medqart':
        cen_ent     = 2
        errlow_ent  = 5
        errhigh_ent = 7
    elif datsetup == 'med95':
        cen_ent     = 2
        errlow_ent  = 3
        errhigh_ent = 9
    elif datsetup == 'medstd':
        cen_ent     = 2
        errlow_ent  = 1
        errhigh_ent = 1
    elif datsetup == 'meanstd':
        cen_ent     = 0
        errlow_ent  = 1
        errhigh_ent = 1

    xvalues  = [[],[],[],[]]
    xerrlow  = [[],[],[],[]]
    xerrhigh = [[],[],[],[]]
    yvalues  = [[],[],[],[]]
    yerrlow  = [[],[],[],[]]
    yerrhigh = [[],[],[],[]]

    for oo in np.arange(len(collectionstats[0])):
        for cc, colstat in enumerate(collectionstats):
            xvalues[cc].append(colstat[oo][xkey][cen_ent])
            xerrlow[cc].append(colstat[oo][xkey][errlow_ent])
            xerrhigh[cc].append(colstat[oo][xkey][errhigh_ent])

            yvalues[cc].append(colstat[oo][ykey][cen_ent])
            yerrlow[cc].append(colstat[oo][ykey][errlow_ent])
            yerrhigh[cc].append(colstat[oo][ykey][errhigh_ent])

    fig = plt.figure(1, figsize=(6, 6))
    # definitions for the axes
    left, width = 0.15, 0.60
    bottom, height = 0.15, 0.60
    bottom_h = left_h = left + width + 0.01

    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=left, right=left+width, bottom=bottom, top=bottom+height)
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    Fsize    = 14
    lthick   = 2
    marksize = 6
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    plt.xlim(xrange)
    xminsys, xmaxsys = plt.xlim()

    plt.ylim(yrange)
    yminsys, ymaxsys = plt.ylim()

    for oo in np.arange(len(collectionstats[0])): # loop necessary for coloring and markers
        objid = collectionstats[0][oo]['id']
        mfc   = True
        if (objid < 6e8): # CDFS and COSMOS
            markersym   = 'o'
        elif (objid < 7e8) & (objid > 6e8): # UDF
            markersym   = 'D'
        elif (objid < 9e8) & (objid > 7e8): # UDF10
            markersym   = 'X'
        elif (objid > 1e9): # Literature objects
            markersym   = lce.get_reference_fromID(objid,verbose=False)[4]
            mfc         = False
        else:
            print(' WARNING - stopped as could not assign a marker symbol to the id '+str(objid))
            pdb.set_trace()

        ms          = marksize
        limsizefrac = 0.05

        # change color of limits
        markerzorder = 20

        for cc, colcol in enumerate(collectionscolors):
            if mfc:
                markerfacecolor = colcol
            else:
                markerfacecolor = 'None'

            if datsetup == 'meanstd':
                plt.errorbar(np.asarray(xvalues[cc])[oo],np.asarray(yvalues[cc])[oo],
                             xerr=np.asarray(xerrlow[cc])[oo],yerr=np.asarray(yerrlow[cc])[oo],capthick=0.5,
                             marker=markersym,lw=lthick/2., markersize=ms, alpha=0.5,
                             markerfacecolor=markerfacecolor,ecolor=colcol,
                             markeredgecolor=colcol,zorder=markerzorder)
            else:
                plt.errorbar(np.asarray(xvalues[cc])[oo],np.asarray(yvalues[cc])[oo],
                             xerr=[[np.asarray(xvalues[cc])[oo]-np.asarray(xerrlow[cc])[oo],
                                    np.asarray(xerrhigh[cc])[oo]-np.asarray(xvalues[cc])[oo]]],
                             yerr=[[np.asarray(yvalues[cc])[oo]-np.asarray(yerrlow[cc])[oo],
                                    np.asarray(yerrhigh[cc])[oo]-np.asarray(yvalues[cc])[oo]]],capthick=0.5,
                             marker=markersym,lw=lthick/2., markersize=ms, alpha=0.5,
                             markerfacecolor=markerfacecolor,ecolor=colcol,
                             markeredgecolor=colcol,zorder=markerzorder)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if ylog:
        plt.yscale('log')
    if xlog:
        plt.xscale('log')

    # --------- HISTOGRAMS ---------
    Nbins   = 50
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    axHistx.xaxis.set_major_formatter(NullFormatter())
    axHisty.yaxis.set_major_formatter(NullFormatter())

    binwidth_x = np.diff([xminsys,xmaxsys])/Nbins
    bindefs    = np.arange(xminsys, xmaxsys+binwidth_x, binwidth_x)
    if xlog:
        bindefs = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))
        axHistx.set_xscale('log')

    for cc, colcol in enumerate(collectionscolors):
        if len(xvalues[cc]) > 0:
            axHistx.hist(np.asarray(xvalues[cc])[objidall < 1e9][np.isfinite(np.asarray(xvalues[cc])[objidall < 1e9])],linestyle='-',
                         bins=bindefs,histtype='step',color=colcol)
            axHistx.hist(np.asarray(xvalues[cc])[np.isfinite(np.asarray(xvalues[cc]))],linestyle=':',
                         bins=bindefs,histtype='step',color=colcol)

    axHistx.set_xticks([])
    axHistx.set_xlim([xminsys,xmaxsys])

    binwidth_y = np.diff([yminsys,ymaxsys])/Nbins
    bindefs    = np.arange(yminsys, ymaxsys+binwidth_y, binwidth_y)
    if ylog:
        bindefs = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))
        axHisty.set_yscale('log')

    for cc, colcol in enumerate(collectionscolors):

        if len(yvalues[cc]) > 0:
            axHisty.hist(np.asarray(yvalues[cc])[objidall < 1e9][np.isfinite(np.asarray(yvalues[cc])[objidall < 1e9])],linestyle='-',
                         bins=bindefs,histtype='step',color=colcol, orientation='horizontal')
            axHisty.hist(np.asarray(yvalues[cc])[np.isfinite(np.asarray(yvalues[cc]))],linestyle=':',
                         bins=bindefs,histtype='step',color=colcol, orientation='horizontal')

    axHisty.set_yticks([])
    axHisty.set_ylim([yminsys,ymaxsys])

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_modelparametercollections(plotname, paramcollections, collectionstats, collectionscolors, maxPDFyscale=False,
                                   constraintsstr=None, fluxratiodictionarylist=None, verbose=True):
    """

    --- INPUT ---
    plotname

    --- EXAMPLE OF USE ---

    """
    yrange = None # defualt until changed with maxPDFyscale.
    Nobj = len(paramcollections[0])
    if verbose: print(' - Will generate plots of photoionization "PDFs" for all '+str(Nobj)+' objects in parameter collections')
    for oo in np.arange(Nobj):
        objid        = paramcollections[0][oo]['id']
        if verbose:
            infostr = '   plotting info for '+str(objid)+' ('+str("%.5d" % (oo+1))+' / '+str("%.5d" % Nobj)+')        '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()
        plotname_obj = plotname.replace('.pdf','_id'+str(objid)+'.pdf')
        # if verbose: print(' - Generating the figure '+plotname_obj)
        figuresize_x = 8
        figuresize_y = 5
        fig          = plt.figure(figsize=(figuresize_x,figuresize_y))
        Fsize        = 8
        LW           = 2
        plt.rc('text', usetex=True)                         # enabling LaTex rendering of text
        plt.rc('font', family='serif',size=Fsize)           # setting text font
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()

        left   = 0.10   # the left side of the subplots of the figure
        right  = 0.95   # the right side of the subplots of the figure
        bottom = 0.10   # the bottom of the subplots of the figure
        top    = 0.85   # the top of the subplots of the figure
        wspace = 1.50   # the amount of width reserved for blank space between subplots
        hspace = 0.50   # the amount of height reserved for white space between subplots
        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        Nrows, Ncols      = 3, 8

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Nmodels_list      = []
        for paramcol in paramcollections:
            Nmodels_list.append(float(len(paramcol[oo]['Zgas'])))

        ylabel            = 'Number of models satisfying observational constraints'
        # \colorbox{BurntOrange}

        N_SFtot  = 10621.
        N_AGNtot = 5184.
        N_bintot = 40131.
        N_sintot = 40131.

        titlestr = 'Observational contraints for ID='+str(objid)+': '

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # if (Nmodels_list[0] > 0) & (Nmodels_list[1] > 0):
        #     Nmodels_ratio     = Nmodels_list[0]/(Nmodels_list[1]*2.0)
        #     titlestr_addition = '; SF/(2*AGN)='+str("%.4f" % Nmodels_ratio)
        #     titlestr          = titlestr+titlestr_addition

        if fluxratiodictionarylist is not None:
            constraints     = fluxratiodictionarylist[oo]

            constraintslist = []
            for key in constraints.keys():
                if key not in ['id']:
                    lowstr = str("%.2f" % constraints[key][0])

                    if constraints[key][1] == 1e35:
                        highstr = '1e35'
                    else:
                        highstr = str("%.2f" % constraints[key][1])

                    constraintslist.append(key+':['+lowstr+','+highstr+']')
            # constraintslist = [key+':['+str("%.2e" % constraints[key][0])+','+str("%.2f" % constraints[key][1])+']'
            #                    for key in constraints.keys() if key not in ['id']]
            #
            # if 'OIII1663/HeII1640' in constraints.keys():
            #     if constraints['OIII1663/HeII1640'][1] > 1e30: pdb.set_trace()

        Ncons = 0
        for constraint in constraintslist:
            if '[0.00,1e35]' not in constraint:
                if (Ncons == 3) & ~titlestr.endswith('\n '):
                    titlestr = titlestr+'\n '

                if (Ncons == 7) & ~titlestr.endswith('\n '):
                    titlestr = titlestr+'\n '

                if (Ncons == 11) & ~titlestr.endswith('\n '):
                    titlestr = titlestr+'\n '

                if constraint.endswith(',1e35]'):
                    titlestr = titlestr+constraint.replace(',1e35]','').replace('[',' $>$ ').replace(':','')+'; '
                elif constraint.split(':')[-1].startswith('[0.00,'):
                    titlestr = titlestr+constraint.replace('[0.00,',' $<$ ').replace(']','').replace(':','')+'; '
                else:
                    titlestr = titlestr+constraint+'; '
                Ncons += 1

        if Ncons == 0:
            titlestr = titlestr+' None'

        # titlestr = titlestr.replace('10000000000.00','1e10')
        if '10000000000.00' in titlestr:
            pdb.set_trace()
        fig.suptitle(titlestr,fontsize=Fsize)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Zgas
        plt.subplot(Nrows, Ncols, (1,4))

        # bindefs = np.array([5e-6,5e-5,1e-4, 0.0002, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.010, 0.014,
        #                     0.017, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08])
        bindefs = 10**np.arange(-6,-0.8,0.1)

        if maxPDFyscale:
            yrange = [0,3500]

        plotkey = 'Zgas'
        pp.plot_modelparametercollections_addhist(paramcollections,collectionstats,collectionscolors,oo,plotkey,
                                                  LW,bindefs=bindefs,Nbins=None,yrange=yrange)

        plt.xscale('log')
        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([5e-6,0.1])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # logUs
        plt.subplot(Nrows, Ncols, (5,8))

        bindefs = np.arange(-5.5, -0.0, 0.05)-0.025

        if maxPDFyscale:
            yrange = [0,2200]

        plotkey = 'logUs'
        pp.plot_modelparametercollections_addhist(paramcollections,collectionstats,collectionscolors,oo,plotkey,
                                                  LW,bindefs=bindefs,Nbins=None,yrange=yrange)

        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([np.min(bindefs),np.max(bindefs)])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # xid
        plt.subplot(Nrows, Ncols, (9,12))

        bindefs = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])-0.05

        if maxPDFyscale:
            yrange = [0,4000]

        plotkey = 'xid'
        pp.plot_modelparametercollections_addhist(paramcollections,collectionstats,collectionscolors,oo,plotkey,
                                                  LW,bindefs=bindefs,Nbins=None,yrange=yrange)

        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([-0.05,0.65])
        plt.ylabel(ylabel)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # nh
        plt.subplot(Nrows, Ncols, (13,16))

        bindefs = 10**(np.arange(0.0, 5.0, 0.25)-0.125)

        if maxPDFyscale:
            yrange = [0,6300]

        plotkey = 'nh'
        pp.plot_modelparametercollections_addhist(paramcollections,collectionstats,collectionscolors,oo,plotkey,
                                                  LW,bindefs=bindefs,Nbins=None,yrange=yrange)
        plt.xscale('log')
        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([0.5,1e5])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # COCOsol
        plt.subplot(Nrows, Ncols, (17,18))

        #bindefs = np.array([0.10, 0.14, 0.20, 0.27, 0.38, 0.52, 0.72, 1.00, 1.40])
        bindefs = 10**np.arange(-1.075,1.0,0.05) #np.arange(0.05,1.5,0.06)

        if maxPDFyscale:
            yrange = [0,1500]

        plotkey = 'COCOsol'
        pp.plot_modelparametercollections_addhist(paramcollections,collectionstats,collectionscolors,oo,plotkey,
                                                  LW,bindefs=bindefs,Nbins=None,yrange=yrange)
        plt.xscale('log')
        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([0.07,2.0])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # mup
        plt.subplot(Nrows, Ncols, (19,20))

        bindefs = np.arange(0,400,50)-25

        if maxPDFyscale:
            yrange = [0,6500]

        plotkey = 'mup'
        pp.plot_modelparametercollections_addhist(paramcollections,collectionstats,collectionscolors,oo,plotkey,
                                                  LW,bindefs=bindefs,Nbins=None,yrange=yrange)

        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([-10,410])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # alpha
        plt.subplot(Nrows, Ncols, (21,22))

        bindefs = np.arange(-2.05,-1.0,0.1)

        if maxPDFyscale:
            yrange = [0,1500]

        plotkey = 'alpha'
        pp.plot_modelparametercollections_addhist(paramcollections,collectionstats,collectionscolors,oo,plotkey,
                                                  LW,bindefs=bindefs,Nbins=None,yrange=yrange)

        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([-2.2,-0.9])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # logAge
        plt.subplot(Nrows, Ncols, (23,24))

        bindefs = np.arange(5.8,8.2,0.05)-0.025

        if maxPDFyscale:
            yrange = [0,2200]

        plotkey = 'logAge'
        pp.plot_modelparametercollections_addhist(paramcollections,collectionstats,collectionscolors,oo,plotkey,
                                                  LW,bindefs=bindefs,Nbins=None,yrange=yrange)

        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([5.8,8.2])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #--------- LEGEND ---------
        plt.plot([-1,-1],[0,0],color=collectionscolors[0],linestyle='-',lw=LW,
                 label='NEOGAL SF (Gutkin+16)\n'
                       'N(model)='+ str(int(Nmodels_list[0]))+' ('+str("%.2f" % (Nmodels_list[0]/N_SFtot *100.0))+'\%)')
        plt.plot([-1,-1],[0,0],color=collectionscolors[1],linestyle='-',lw=LW,
                 label='NEOGAL AGN (Feltre+16)\n'
                       'N(model)='+ str(int(Nmodels_list[1]))+' ('+str("%.2f" % (Nmodels_list[1]/N_AGNtot*100.0))+'\%)')
        plt.plot([-1,-1],[0,0],color=collectionscolors[2],linestyle='-',lw=LW,
                 label='BPASS binary (Xiao+18)\n'
                       'N(model)='+ str(int(Nmodels_list[2]))+' ('+str("%.2f" % (Nmodels_list[2]/N_bintot*100.0))+'\%)')
        plt.plot([-1,-1],[0,0],color=collectionscolors[3],linestyle='-',lw=LW,
                 label='BPASS single (Xiao+18)\n'
                       'N(model)='+ str(int(Nmodels_list[3]))+' ('+str("%.2f" % (Nmodels_list[3]/N_sintot*100.0))+'\%)')

        leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize},ncol=4,numpoints=1,
                         bbox_to_anchor=(-1.75, 4.45),)  # add the legend
        leg.get_frame().set_alpha(0.7)
        #--------------------------

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.savefig(plotname_obj)
        plt.clf()
        plt.close('all')
        # if verbose: print(' - Successfully saved figure to file')
    if verbose: print('\n   done...')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_modelparametercollections_addhist(paramcollections,collectionstats,collectionscolors,objindex,plotkey,LW,
                                           bindefs=None,Nbins=None,yrange=None):
    """

    --- INPUT ---

    --- EXAMPLE OF USE ---
    see pp.plot_modelparametercollections() above

    """
    paramcol = []
    for cc, pc in enumerate(paramcollections):
        if pc[objindex][plotkey] is not None:
            paramcol.append(np.asarray(pc[objindex][plotkey]))
        else:
            paramcol.append(None)
    colstat  = [np.asarray(cs[objindex][plotkey]) for cs in collectionstats]

    if (paramcol[0] is None) & (paramcol[1] is None):#<------------------------------------------ BPASS logAge histogram
        if (len(paramcol[2]) == 0) & (len(paramcol[3]) == 0):
            pass
        else:
            colindices = [2,3]
            if bindefs is None:
                allval = []
                for colindex in colindices:
                    allval.append(allval,paramcol[colindex])
                xmin       = np.min(allval)*0.95
                xmax       = np.max(allval)*1.05
                binwidth_x = np.diff([xmin,xmax])/Nbins
                bindefs    = np.arange(xmin, xmax+binwidth_x, binwidth_x)

            for colindex in colindices:
                plt.hist(paramcol[colindex], bins=bindefs,histtype='step',color=collectionscolors[colindex])

            if yrange is None:
                yminsys, ymaxsys = plt.ylim()
            else:
                yminsys, ymaxsys = yrange[0], yrange[1]

            for colindex in colindices:
                plt.plot([colstat[colindex][0]]*2, [yminsys, ymaxsys], '--',lw=LW,alpha=0.5,
                         color=collectionscolors[colindex])  # mean
                plt.plot([colstat[colindex][2]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,
                         color=collectionscolors[colindex])  # median
                plt.plot([colstat[colindex][4]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,
                         color=collectionscolors[colindex])  # 68% confidence interval lower
                plt.plot([colstat[colindex][8]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,
                         color=collectionscolors[colindex])  # 68% confidence interval lower

    elif (paramcol[0] is None) or (paramcol[1] is None):
        if (paramcol[0] is None) & (paramcol[1] is not None):#<------------------------------ NEOGAL AGN only parameters
            colindex = 1
        elif (paramcol[0] is not None) & (paramcol[1] is None):#<----------------------------- NEOGAL SF only parameters
            colindex = 0
        else:
            sys.exit('no setup for parameter collection combination with both NEOGAL none collection')

        if len(paramcol[colindex]) > 0:
            if bindefs is None:
                xmin       = np.min(paramcol[colindex])-np.abs(np.min(paramcol[colindex]))*0.05
                xmax       = np.max(paramcol[colindex])+np.abs(np.min(paramcol[colindex]))*0.05
                binwidth_x = np.diff([xmin,xmax])/Nbins
                bindefs    = np.arange(xmin, xmax+binwidth_x, binwidth_x)

            plt.hist(paramcol[colindex], bins=bindefs,histtype='step',color=collectionscolors[colindex])

            if yrange is None:
                yminsys, ymaxsys = plt.ylim()
            else:
                yminsys, ymaxsys = yrange[0], yrange[1]

            plt.plot([colstat[colindex][0]]*2, [yminsys, ymaxsys], '--',lw=LW,alpha=0.5,
                     color=collectionscolors[colindex])  # mean
            plt.plot([colstat[colindex][2]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,
                     color=collectionscolors[colindex])  # median
            plt.plot([colstat[colindex][4]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,
                     color=collectionscolors[colindex])  # 68% confidence interval lower
            plt.plot([colstat[colindex][8]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,
                     color=collectionscolors[colindex])  # 68% confidence interval lower

    else:
        lenghtlist = []
        for pc in paramcol:
            if pc is not None:
                lenghtlist.append(len(pc))
            else:
                lenghtlist.append(0)

        if lenghtlist == [0]*len(lenghtlist):
            pass
        else:
            colindices = [0,1,2,3]
            if bindefs is None:
                allval = []
                for colindex in colindices:
                    allval.append(allval,paramcol[colindex])
                xmin       = np.min(allval)*0.95
                xmax       = np.max(allval)*1.05
                binwidth_x = np.diff([xmin,xmax])/Nbins
                bindefs    = np.arange(xmin, xmax+binwidth_x, binwidth_x)

            for colindex in colindices:
                if paramcol[colindex] is not None:
                    plt.hist(paramcol[colindex], bins=bindefs,histtype='step',color=collectionscolors[colindex])

            if yrange is None:
                yminsys, ymaxsys = plt.ylim()
            else:
                yminsys, ymaxsys = yrange[0], yrange[1]

            for colindex in colindices:
                if paramcol[colindex] is not None:
                    plt.plot([colstat[colindex][0]]*2, [yminsys, ymaxsys], '--',lw=LW,alpha=0.5,
                             color=collectionscolors[colindex])  # mean
                    plt.plot([colstat[colindex][2]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,
                             color=collectionscolors[colindex])  # median
                    plt.plot([colstat[colindex][4]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,
                             color=collectionscolors[colindex])  # 68% confidence interval lower
                    plt.plot([colstat[colindex][8]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,
                             color=collectionscolors[colindex])  # 68% confidence interval lower

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def keylabels(keyinput):
    """
    Function returning LaTeX label for photionisation keyword

    """
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def keylabels(keyinput):
    """
    Function returning LaTeX label for photoionization keywords

    """
    labeldic = {}
    labeldic['Zgas']    =  'Z'
    labeldic['logUs']   =  'log(U)'
    labeldic['xid']     =  '$\\xi_\\textrm{d}$'
    labeldic['nh']      =  '$n_\\textrm{H}$ [cm$^{-3}$]'
    labeldic['COCOsol'] =  'C/O [(C/O)$_\\textrm{sun}$]'
    labeldic['mup']     =  'm$_\\textrm{up}$ [M$_\\textrm{sun}$]'
    labeldic['alpha']   =  '$\\alpha$'
    labeldic['logAge']  =  'log(Age) [yr]'
    return labeldic[keyinput]
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =