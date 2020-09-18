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
def estimate_object_PDFs(fluxratiodictionarylist,generatePDFplots=False,basename='photoionizationmodelPDFs',
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
    FRdic = [{'id':111111111111, 'HeII1640/OIII1663':[0.04,0.45],'CIII1908/CIV1550':[1.0,10.0]}, {'id':222222222222, 'OIII1663/HeII1640':[1e-1,1.0],'CIII1908/CIV1550':[0.1,10.0]}, {'id':333333, 'OIII1663/HeII1640':[1e2,1e3],'CIII1908/CIV1550':[1e-3,1e-2]}, {'id':444444, 'OIII1663/HeII1640':[1e-2,1e-1],'CIII1908/CIV1550':[5e-1,1e-0]}, {'id':555555, 'OIII1663/HeII1640':[1e-2,1e10],'CIII1908/CIV1550':[5e-1,1e-0], 'OIII1663/CIII1908':[1e-2,1e10], 'OIII1663/CIV1550':[1e-2,1e10], 'OIII1663/SiIII1888':[1e-2,1e1], 'CIII1908/SiIII1888':[1e-2,1e10], 'CIV1550/SiIII1888':[1e-2,1e10]}]

    import photoionizationPDFs as pp
    basename= '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/photoionizationPDFs/photoionizationmodelPDFs'
    paramcollections, collectionstats = pp.estimate_object_PDFs(FRdic, basename=basename, generatePDFplots=True)

    FRdic = [{'id':111111111111, 'HeII1640/OIII1663':[0.0,1e10]}] # run for a single objects with no constraints to get instrinsic distribution

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
    fluxratiodic['NV1240/CIV1550']     =   [0,1e10]
    fluxratiodic['NV1240/CIII1908']    =   [0,1e10]
    fluxratiodic['NV1240/HeII1640']    =   [0,1e10]
    fluxratiodic['NV1240/OIII1663']    =   [0,1e10]
    fluxratiodic['NV1240/SiIII1888']   =   [0,1e10]

    fluxratiodic['CIV1550/NV1240']     =   [0,1e10]
    fluxratiodic['CIV1550/CIII1908']   =   [0,1e10]
    fluxratiodic['CIV1550/HeII1640']   =   [0,1e10]
    fluxratiodic['CIV1550/OIII1663']   =   [0,1e10]
    fluxratiodic['CIV1550/SiIII1888']  =   [0,1e10]

    fluxratiodic['CIII1908/NV1240']     =  [0,1e10]
    fluxratiodic['CIII1908/CIV1550']    =  [0,1e10]
    fluxratiodic['CIII1908/HeII1640']   =  [0,1e10]
    fluxratiodic['CIII1908/OIII1663']   =  [0,1e10]
    fluxratiodic['CIII1908/SiIII1888']  =  [0,1e10]

    fluxratiodic['HeII1640/NV1240']     =  [0,1e10]
    fluxratiodic['HeII1640/CIV1550']    =  [0,1e10]
    fluxratiodic['HeII1640/CIII1908']   =  [0,1e10]
    fluxratiodic['HeII1640/OIII1663']   =  [0,1e10]
    fluxratiodic['HeII1640/SiIII1888']  =  [0,1e10]

    fluxratiodic['OIII1663/NV1240']     =  [0,1e10]
    fluxratiodic['OIII1663/CIV1550']    =  [0,1e10]
    fluxratiodic['OIII1663/CIII1908']   =  [0,1e10]
    fluxratiodic['OIII1663/HeII1640']   =  [0,1e10]
    fluxratiodic['OIII1663/SiIII1888']  =  [0,1e10]

    fluxratiodic['SiIII1888/NV1240']    =  [0,1e10]
    fluxratiodic['SiIII1888/CIV1550']   =  [0,1e10]
    fluxratiodic['SiIII1888/CIII1908']  =  [0,1e10]
    fluxratiodic['SiIII1888/HeII1640']  =  [0,1e10]
    fluxratiodic['SiIII1888/OIII1663']  =  [0,1e10]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Set up mode line flux vectors')
    fluxdic = {}
    fluxdic['NV1240']    = [NEOGAL_SF['NV1240'],
                            NEOGAL_AGN['NV1240'],
                            None,
                            None]

    fluxdic['CIV1550']   = [NEOGAL_SF['CIV1548']+NEOGAL_SF['CIV1551'],
                            NEOGAL_AGN['CIV1548']+NEOGAL_AGN['CIV1551'],
                            BPASS_bin['CIV1548']+BPASS_bin['CIV1551'],
                            BPASS_sin['CIV1548']+BPASS_sin['CIV1551']]

    fluxdic['CIII1908']  = [NEOGAL_SF['CIII1908'],
                            NEOGAL_AGN['CIII1907']+NEOGAL_AGN['CIII1910'],
                            BPASS_bin['CIII1907']+BPASS_bin['CIII1910'],
                            BPASS_sin['CIII1907']+BPASS_sin['CIII1910']]

    fluxdic['HeII1640']  = [NEOGAL_SF['HeII1640'],
                            NEOGAL_AGN['HeII1640'],
                            BPASS_bin['HeII1640'],
                            BPASS_sin['HeII1640']]

    fluxdic['OIII1663']  = [NEOGAL_SF['OIII1661']+NEOGAL_SF['OIII1666'],
                            NEOGAL_AGN['OIII1661']+NEOGAL_AGN['OIII1666'],
                            BPASS_bin['OIII1661']+BPASS_bin['OIII1666'],
                            BPASS_sin['OIII1661']+BPASS_sin['OIII1666']]

    fluxdic['SiIII1888'] = [NEOGAL_SF['SiIII1888'],
                            NEOGAL_AGN['SiIII1888'],
                            None,
                            None]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Get ranges of model flux ratios')
    Nobj = len(fluxratiodictionarylist)
    if verbose: print(' - Get model selection given flux ratio ranges according to '+
                      str(Nobj)+" object's data provided ")
    if verbose: print('   Selection based on the total number of photoionization models: \n'
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
                goodent_FR_bin = np.where( (fluxdic[numerator][2]/fluxdic[denominator][2] >= fluxratioconstrains_obj[FR][0]) &
                                           (fluxdic[numerator][2]/fluxdic[denominator][2] <= fluxratioconstrains_obj[FR][1]))[0]
                goodent_bin    = np.intersect1d(goodent_bin,goodent_FR_bin)
            else:
                if verbose: print(' WARNING: No constraints on flux ratio '+FR+' for BPASS binary models')

            if (fluxdic[numerator][2] is not None) & (fluxdic[denominator][2] is not None):
                goodent_FR_sin = np.where( (fluxdic[numerator][3]/fluxdic[denominator][3] >= fluxratioconstrains_obj[FR][0]) &
                                           (fluxdic[numerator][3]/fluxdic[denominator][3] <= fluxratioconstrains_obj[FR][1]))[0]
                goodent_sin    = np.intersect1d(goodent_sin,goodent_FR_sin)
            else:
                if verbose: print(' WARNING: No constraints on flux ratio '+FR+' for BPASS singular models')


        parametercollection_SF[oo]  = {'id'     : FRdic_input['id'],
                                       'Zgas'   : NEOGAL_SF['Zgas'][goodent_SF],
                                       'logUs'  : NEOGAL_SF['logUs'][goodent_SF],
                                       'xid'    : NEOGAL_SF['xid'][goodent_SF],
                                       'nh'     : NEOGAL_SF['nh'][goodent_SF],
                                       'COCOsol': NEOGAL_SF['COCOsol'][goodent_SF],
                                       'mup'    : NEOGAL_SF['mup'][goodent_SF]}

        parametercollection_AGN[oo] = {'id'     : FRdic_input['id'],
                                       'Zgas'   : NEOGAL_AGN['Zgas'][goodent_AGN],
                                       'logUs'  : NEOGAL_AGN['logUs'][goodent_AGN],
                                       'xid'    : NEOGAL_AGN['xid'][goodent_AGN],
                                       'nh'     : NEOGAL_AGN['nh'][goodent_AGN],
                                       'alpha'  : NEOGAL_AGN['alpha'][goodent_AGN]}

        parametercollection_bin[oo] = {'id'     : FRdic_input['id'],
                                       'Zgas'   : BPASS_bin['Zgas'][goodent_bin],
                                       'logUs'  : BPASS_bin['logU'][goodent_bin],
                                       'nh'     : 10**(BPASS_bin['lognH'][goodent_bin]),
                                       'logAge' : BPASS_bin['logAge'][goodent_bin]}

        parametercollection_sin[oo] = {'id'     : FRdic_input['id'],
                                       'Zgas'   : BPASS_sin['Zgas'][goodent_sin],
                                       'logUs'  : BPASS_sin['logU'][goodent_sin],
                                       'nh'     : 10**(BPASS_sin['lognH'][goodent_sin]),
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
        pp.plot_modelparametercollections(plotname, parametercollection_SF, parametercollection_AGN,
                                          stat_SF, stat_AGN, AGNcol=col_NEOGAL_AGN,SFcol=col_NEOGAL_SF,
                                          fluxratiodictionarylist=fluxratiodictionarylist,
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

            if len(parametercollection[oo][key]) > 0:

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

    return statparam
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_stat(plotname, collectionstats, collectionscolors, verbose=True):
    """

    """
    if verbose: print(' - Plot stats of parameter selection')
    Nobj     = len(collectionstats[0])
    objidall = np.array([statobj['id'] for statobj in collectionstats[0]])

    xranges = {}
    xranges['Zgas']    =  [0.00001,0.2]
    xranges['logUs']   =  [-5.2,-0.4]
    xranges['xid']     =  [0.02,0.58]
    xranges['nh']      =  [1e1,1e5]
    xranges['COCOsol'] =  [0,1.5]
    xranges['mup']     =  [50,350]
    xranges['alpha']   =  [-2.2,-0.9]
    xranges['logAge']  =  [1,15]

    keylists = []
    for colstat in collectionstats:
        keylists.append([key for key in colstat[0].keys()])

    uniquekeys = np.unique(keylists[0]+keylists[1]+keylists[2]+keylists[3])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for key in uniquekeys:
        if key == 'id': continue

        meanval_SF  = np.array([])
        stdval_SF   = np.array([])
        range50_SF  = np.array([])
        meanval_AGN = np.array([])
        stdval_AGN  = np.array([])
        range50_AGN = np.array([])
        meanval_bin = np.array([])
        stdval_bin  = np.array([])
        range50_bin = np.array([])
        meanval_sin = np.array([])
        stdval_sin  = np.array([])
        range50_sin = np.array([])

        for oo in np.arange(Nobj):
            listent = 0 # NEOGAL SF
            if key in keylists[listent]:
                meanval_SF = np.append(meanval_SF,collectionstats[listent][oo][key][0])
                stdval_SF = np.append(stdval_SF,collectionstats[listent][oo][key][1])
                range50_SF = np.append(range50_SF,collectionstats[listent][oo][key][7]-collectionstats[listent][oo][key][7])

            listent = 1 # NEOGAL AGN
            if key in keylists[listent]:
                meanval_AGN = np.append(meanval_AGN,collectionstats[listent][oo][key][0])
                stdval_AGN = np.append(stdval_AGN,collectionstats[listent][oo][key][1])
                range50_AGN = np.append(range50_AGN,collectionstats[listent][oo][key][7]-collectionstats[listent][oo][key][7])

            listent = 2 # BPASS binaries
            if key in keylists[listent]:
                meanval_bin = np.append(meanval_bin,collectionstats[listent][oo][key][0])
                stdval_bin = np.append(stdval_bin,collectionstats[listent][oo][key][1])
                range50_bin = np.append(range50_bin,collectionstats[listent][oo][key][7]-collectionstats[listent][oo][key][7])

            listent = 3 # BPASS singles
            if key in keylists[listent]:
                meanval_sin = np.append(meanval_sin,collectionstats[listent][oo][key][0])
                stdval_sin = np.append(stdval_sin,collectionstats[listent][oo][key][1])
                range50_sin = np.append(range50_sin,collectionstats[listent][oo][key][7]-collectionstats[listent][oo][key][7])

        #-------- PLOT 'EM --------
        plotname_stat = plotname.replace('.pdf','_meanVSerr_'+key+'.pdf')

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


            SFcol = collectionscolors[0]
            if mfc:
                markerfacecolor = SFcol
            else:
                markerfacecolor = 'None'

            if len(meanval_SF) > 0:
                plt.errorbar(meanval_SF[oo],stdval_SF[oo],xerr=None,yerr=None,capthick=0.5,
                             marker=markersym,lw=lthick/2., markersize=ms, alpha=0.5,
                             markerfacecolor=markerfacecolor,ecolor=SFcol,
                             markeredgecolor=SFcol,zorder=markerzorder)

            AGNcol = collectionscolors[1]
            if mfc:
                markerfacecolor = AGNcol
            else:
                markerfacecolor = 'None'

            if len(meanval_AGN) > 0:
                plt.errorbar(meanval_AGN[oo],stdval_AGN[oo],xerr=None,yerr=None,capthick=0.5,
                             marker=markersym,lw=lthick/2., markersize=ms, alpha=0.5,
                             markerfacecolor=markerfacecolor,ecolor=AGNcol,
                             markeredgecolor=AGNcol,zorder=markerzorder)

            bincol = collectionscolors[2]
            if mfc:
                markerfacecolor = bincol
            else:
                markerfacecolor = 'None'

            if len(meanval_bin) > 0:
                plt.errorbar(meanval_bin[oo],stdval_bin[oo],xerr=None,yerr=None,capthick=0.5,
                             marker=markersym,lw=lthick/2., markersize=ms, alpha=0.5,
                             markerfacecolor=markerfacecolor,ecolor=bincol,
                             markeredgecolor=bincol,zorder=markerzorder)

            sincol = collectionscolors[3]
            if mfc:
                markerfacecolor = sincol
            else:
                markerfacecolor = 'None'

            if len(meanval_sin) > 0:
                plt.errorbar(meanval_sin[oo],stdval_sin[oo],xerr=None,yerr=None,capthick=0.5,
                             marker=markersym,lw=lthick/2., markersize=ms, alpha=0.5,
                             markerfacecolor=markerfacecolor,ecolor=sincol,
                             markeredgecolor=sincol,zorder=markerzorder)

        plt.xlim(xranges[key])
        xminsys, xmaxsys = plt.xlim()

        # plt.ylim(yranges[key])
        yminsys, ymaxsys = plt.ylim()

        plt.xlabel(pp.keylabels(key))
        plt.ylabel('Standard deviation of parameter distribution')

        if (key == 'Zgas') or (key == 'nh'):
            xlog = True
            plt.xscale('log')
        else:
            xlog = False

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

        if len(meanval_SF) > 0:
            axHistx.hist(meanval_SF[objidall < 1e9][np.isfinite(meanval_SF[objidall < 1e9])],linestyle='-',
                         bins=bindefs,histtype='step',color=SFcol)
            axHistx.hist(meanval_SF[np.isfinite(meanval_SF)],linestyle=':',
                         bins=bindefs,histtype='step',color=SFcol)

        if len(meanval_AGN) > 0:
            axHistx.hist(meanval_AGN[objidall < 1e9][np.isfinite(meanval_AGN[objidall < 1e9])],linestyle='-',
                         bins=bindefs,histtype='step',color=AGNcol)
            axHistx.hist(meanval_AGN[np.isfinite(meanval_AGN)],linestyle=':',
                         bins=bindefs,histtype='step',color=AGNcol)

        if len(meanval_bin) > 0:
            axHistx.hist(meanval_bin[objidall < 1e9][np.isfinite(meanval_bin[objidall < 1e9])],linestyle='-',
                         bins=bindefs,histtype='step',color=bincol)
            axHistx.hist(meanval_bin[np.isfinite(meanval_bin)],linestyle=':',
                         bins=bindefs,histtype='step',color=bincol)

        if len(meanval_sin) > 0:
            axHistx.hist(meanval_sin[objidall < 1e9][np.isfinite(meanval_sin[objidall < 1e9])],linestyle='-',
                         bins=bindefs,histtype='step',color=sincol)
            axHistx.hist(meanval_sin[np.isfinite(meanval_sin)],linestyle=':',
                         bins=bindefs,histtype='step',color=sincol)


        axHistx.set_xticks([])
        axHistx.set_xlim([xminsys,xmaxsys])

        binwidth_y = np.diff([yminsys,ymaxsys])/Nbins
        bindefs    = np.arange(yminsys, ymaxsys+binwidth_y, binwidth_y)
        # if ylog:
        #     bindefs = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))
        #     axHisty.set_yscale('log')
        if len(meanval_SF) > 0:
            axHisty.hist(stdval_SF[objidall < 1e9][np.isfinite(stdval_SF[objidall < 1e9])],linestyle='-',
                         bins=bindefs,histtype='step',color=SFcol, orientation='horizontal')
            axHisty.hist(stdval_SF[np.isfinite(stdval_SF)],linestyle=':',
                         bins=bindefs,histtype='step',color=SFcol, orientation='horizontal')
        if len(meanval_AGN) > 0:
            axHisty.hist(stdval_AGN[objidall < 1e9][np.isfinite(stdval_AGN[objidall < 1e9])],linestyle='-',
                         bins=bindefs,histtype='step',color=AGNcol, orientation='horizontal')
            axHisty.hist(stdval_AGN[np.isfinite(stdval_AGN)],linestyle=':',
                         bins=bindefs,histtype='step',color=AGNcol, orientation='horizontal')
        if len(meanval_bin) > 0:
            axHisty.hist(stdval_bin[objidall < 1e9][np.isfinite(stdval_bin[objidall < 1e9])],linestyle='-',
                         bins=bindefs,histtype='step',color=bincol, orientation='horizontal')
            axHisty.hist(stdval_bin[np.isfinite(stdval_bin)],linestyle=':',
                         bins=bindefs,histtype='step',color=bincol, orientation='horizontal')
        if len(meanval_sin) > 0:
            axHisty.hist(stdval_sin[objidall < 1e9][np.isfinite(stdval_sin[objidall < 1e9])],linestyle='-',
                         bins=bindefs,histtype='step',color=sincol, orientation='horizontal')
            axHisty.hist(stdval_sin[np.isfinite(stdval_sin)],linestyle=':',
                         bins=bindefs,histtype='step',color=sincol, orientation='horizontal')


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
            print(' WARNING - stopped as could not assing a marker symbol to the id '+str(objid))
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
def plot_modelparametercollections(plotname, parametercollection_SF, parametercollection_AGN,
                                   stat_SF, stat_AGN, AGNcol='blue',SFcol='red', constraintsstr=None,
                                   fluxratiodictionarylist=None, verbose=True):
    """

    --- INPUT ---
    plotname
    parametercollection_SF
    parametercollection_AGN
    stat_SF
    stat_AGN
    AGNcol
    SFcol
    verbose

    --- EXAMPLE OF USE ---

    """

    Nobj = len(parametercollection_SF)
    if verbose: print(' - Will generate plots of NEOGAL "PDFs" for all '+str(Nobj)+' objects in parameter collections')
    for oo in np.arange(Nobj):
        objid        = parametercollection_SF[oo]['id']
        if verbose:
            infostr = '   plotting info for '+str(objid)+' ('+str("%.5d" % (oo+1))+' / '+str("%.5d" % Nobj)+')        '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()
        plotname_obj = plotname.replace('.pdf','_id'+str(objid)+'.pdf')
        # if verbose: print(' - Generating the figure '+plotname_obj)
        figuresize_x = 6
        figuresize_y = 5
        fig          = plt.figure(figsize=(figuresize_x,figuresize_y))
        Fsize        = 9
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
        top    = 0.90   # the top of the subplots of the figure
        wspace = 1.50   # the amount of width reserved for blank space between subplots
        hspace = 0.50   # the amount of height reserved for white space between subplots
        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        Nrows, Ncols      = 3, 6
        ylabel            = 'Number of NEOGAL SF ('+str(SFcol)+') and AGN ('+str(AGNcol)+') models'

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Nmodels_SF  = float(len(parametercollection_SF[oo]['Zgas']))
        Nmodels_AGN = float(len(parametercollection_AGN[oo]['Zgas']))

        titlestr = 'Models satisfying ID='+str(objid)+' cuts: SF='+str(Nmodels_SF)+'; AGN='+str(Nmodels_AGN)
        if (Nmodels_AGN > 0) & (Nmodels_SF > 0):
            Nmodels_ratio     = Nmodels_SF/Nmodels_AGN
            titlestr_addition = '; SF/AGN='+str("%.4f" % Nmodels_ratio)
            titlestr          = titlestr+titlestr_addition

        if fluxratiodictionarylist is not None:
            constraints     = fluxratiodictionarylist[oo]
            constraintslist = [key+':['+str("%.2f" % constraints[key][0])+','+str("%.2f" % constraints[key][1])+']'
                               for key in constraints.keys() if key not in ['id']]

            if len(constraintslist) < 4:
                constraintsstr = '; '.join(constraintslist)
            elif (len(constraintslist) > 3) & (len(constraintslist) < 7):
                constraintsstr = '; '.join(constraintslist[:3])+'\n'+'; '.join(constraintslist[3:6])
            elif (len(constraintslist) > 6) & (len(constraintslist) < 10):
                constraintsstr = '; '.join(constraintslist[:3])+'\n'+'; '.join(constraintslist[3:6])+\
                                 '\n'+'; '.join(constraintslist[6:])
            else:
                constraintsstr = '; '.join(constraintslist[:3])+'\n'+'; '.join(constraintslist[3:6])+\
                                 '\n'+'; '.join(constraintslist[6:9])+'\n'+'; '.join(constraintslist[9:])

            constraintsstr = constraintsstr.replace('10000000000.00','1e10')

            titlestr = titlestr+'\n'+constraintsstr
            # titlestr = r'{\fontsize{'+str(Fsize)+'pt}{3em}\selectfont{}{'+titlestr+'\r}{\fontsize{'+str((Fsize-2.))+'pt}{3em}\selectfont{}('+constraintsstr+'}'

        # plt.text(x=0.5, y=0.94, s=titlestr, fontsize=Fsize, ha="center", transform=fig.transFigure)
        # plt.text(x=0.5, y=0.88, s=constraintsstr, fontsize=Fsize-2, ha="center", transform=fig.transFigure)
        # fig.title(titlestr,fontsize=Fsize)
        fig.suptitle(titlestr,fontsize=Fsize-2)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Zgas
        plt.subplot(Nrows, Ncols, (1,3))

        bindefs = np.array([0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.014,
                            0.017, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07])-0.00001

        plotkey = 'Zgas'
        pp.plot_modelparametercollections_addhist(parametercollection_SF[oo][plotkey],parametercollection_AGN[oo][plotkey],
                                                  stat_SF[oo][plotkey],stat_AGN[oo][plotkey],
                                                  SFcol,AGNcol,LW,bindefs=bindefs,Nbins=None)

        plt.xscale('log')
        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([0.00001,0.1])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # logUs
        plt.subplot(Nrows, Ncols, (4,6))

        bindefs = np.arange(-4.75, -0.25, 0.5)

        plotkey = 'logUs'
        pp.plot_modelparametercollections_addhist(parametercollection_SF[oo][plotkey],parametercollection_AGN[oo][plotkey],
                                                  stat_SF[oo][plotkey],stat_AGN[oo][plotkey],
                                                  SFcol,AGNcol,LW,bindefs=bindefs,Nbins=None)

        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([-5,-0.5])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # xid
        plt.subplot(Nrows, Ncols, (7,9))

        bindefs = np.array([0.0, 0.2, 0.4, 0.6])

        plotkey = 'xid'
        pp.plot_modelparametercollections_addhist(parametercollection_SF[oo][plotkey],parametercollection_AGN[oo][plotkey],
                                                  stat_SF[oo][plotkey],stat_AGN[oo][plotkey],
                                                  SFcol,AGNcol,LW,bindefs=bindefs,Nbins=None)

        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([-0.05,0.65])
        plt.ylabel(ylabel)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # nh
        plt.subplot(Nrows, Ncols, (10,12))

        bindefs = 10**np.array([1.5, 2.5, 3.5, 4.5])

        plotkey = 'nh'
        pp.plot_modelparametercollections_addhist(parametercollection_SF[oo][plotkey],parametercollection_AGN[oo][plotkey],
                                                  stat_SF[oo][plotkey],stat_AGN[oo][plotkey],
                                                  SFcol,AGNcol,LW,bindefs=bindefs,Nbins=None)
        plt.xscale('log')
        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([10,1e5])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # COCOsol
        plt.subplot(Nrows, Ncols, (13,14))

        #bindefs = np.array([0.10, 0.14, 0.20, 0.27, 0.38, 0.52, 0.72, 1.00, 1.40])
        bindefs = np.arange(0.05,1.5,0.06)

        plotkey = 'COCOsol'
        pp.plot_modelparametercollections_addhist(parametercollection_SF[oo][plotkey],None,
                                                  stat_SF[oo][plotkey],None,
                                                  SFcol,AGNcol,LW,bindefs=bindefs,Nbins=None)

        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([0.00,1.55])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # mup
        plt.subplot(Nrows, Ncols, (15,16))

        bindefs = np.array([0,200,400])

        plotkey = 'mup'
        pp.plot_modelparametercollections_addhist(parametercollection_SF[oo][plotkey],None,
                                                  stat_SF[oo][plotkey],None,
                                                  SFcol,AGNcol,LW,bindefs=bindefs,Nbins=None)

        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([-10,410])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # alpha
        plt.subplot(Nrows, Ncols, (17,18))

        bindefs = np.array([-2.15,-1.85,-1.55,-1.25,-0.95])

        plotkey = 'alpha'
        Nbins   = 10
        pp.plot_modelparametercollections_addhist(None,parametercollection_AGN[oo][plotkey],
                                                  None,stat_AGN[oo][plotkey],
                                                  SFcol,AGNcol,LW,bindefs=None,Nbins=Nbins)

        plt.xlabel(pp.keylabels(plotkey))
        plt.xlim([-2.2,-0.9])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.savefig(plotname_obj)
        plt.clf()
        plt.close('all')
        # if verbose: print(' - Successfully saved figure to file')
    if verbose: print('\n   done...')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_modelparametercollections_addhist(paramcol_SF,paramcol_AGN,stat_SF,stat_AGN,
                                           SFcol,AGNcol,LW,bindefs=None,Nbins=10):
    """

    --- INPUT ---
    parametercollection_SF
    parametercollection_AGN
    stat_SF
    stat_AGN
    SFcol
    AGNcol
    LW
    bindefs
    Nbins

    --- EXAMPLE OF USE ---
    see pp.plot_modelparametercollections() above

    """
    if (paramcol_SF is None) & (paramcol_AGN is not None):
        if len(paramcol_AGN) > 0:
            if bindefs is None:
                xmin       = np.min(paramcol_AGN)-np.abs(np.min(paramcol_AGN))*0.05
                xmax       = np.max(paramcol_AGN)+np.abs(np.min(paramcol_AGN))*0.05
                binwidth_x = np.diff([xmin,xmax])/Nbins
                bindefs    = np.arange(xmin, xmax+binwidth_x, binwidth_x)

            plt.hist(paramcol_AGN, bins=bindefs,histtype='step',color=AGNcol)

            yminsys, ymaxsys = plt.ylim()

            plt.plot([stat_AGN[0]]*2, [yminsys, ymaxsys], '--',lw=LW,alpha=0.5,color=AGNcol)  # mean
            plt.plot([stat_AGN[2]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=AGNcol)  # median
            plt.plot([stat_AGN[4]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=AGNcol)  # 68% confidence interval lower
            plt.plot([stat_AGN[8]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=AGNcol)  # 68% confidence interval lower

    elif (paramcol_SF is not None) & (paramcol_AGN is None):
        if len(paramcol_SF) > 0:

            if bindefs is None:
                xmin       = np.min(paramcol_SF)-np.abs(np.min(paramcol_SF))*0.05
                xmax       = np.max(paramcol_SF)+np.abs(np.min(paramcol_SF))*0.05
                binwidth_x = np.diff([xmin,xmax])/Nbins
                bindefs    = np.arange(xmin, xmax+binwidth_x, binwidth_x)

            plt.hist(paramcol_SF,  bins=bindefs,histtype='step',color=SFcol)

            yminsys, ymaxsys = plt.ylim()

            plt.plot([stat_SF[0]]*2,  [yminsys, ymaxsys], '--',lw=LW,alpha=0.5,color=SFcol)  # mean
            plt.plot([stat_SF[2]]*2,  [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=SFcol)  # median
            plt.plot([stat_SF[4]]*2,  [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=SFcol)  # 68% confidence interval lower
            plt.plot([stat_SF[8]]*2,  [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=SFcol)  # 68% confidence interval lower

    else:

        if (len(paramcol_SF) == 0) & (len(paramcol_AGN) == 0):
            pass
        else:
            if bindefs is None:
                xmin       = np.min(np.append(paramcol_SF,paramcol_AGN))*0.95
                xmax       = np.max(np.append(paramcol_SF,paramcol_AGN))*1.05
                binwidth_x = np.diff([xmin,xmax])/Nbins
                bindefs    = np.arange(xmin, xmax+binwidth_x, binwidth_x)
                # if xlog:
                #     bindefs = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))
                #     axHistax.set_xscale('log')

            plt.hist(paramcol_SF,  bins=bindefs,histtype='step',color=SFcol)
            plt.hist(paramcol_AGN, bins=bindefs,histtype='step',color=AGNcol)

            yminsys, ymaxsys = plt.ylim()

            if len(paramcol_SF) > 0:
                plt.plot([stat_SF[0]]*2,  [yminsys, ymaxsys], '--',lw=LW,alpha=0.5,color=SFcol)  # mean
                plt.plot([stat_SF[2]]*2,  [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=SFcol)  # median
                plt.plot([stat_SF[4]]*2,  [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=SFcol)  # 68% confidence interval lower
                plt.plot([stat_SF[8]]*2,  [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=SFcol)  # 68% confidence interval lower

            if len(paramcol_AGN) > 0:
                plt.plot([stat_AGN[0]]*2, [yminsys, ymaxsys], '--',lw=LW,alpha=0.5,color=AGNcol)  # mean
                plt.plot([stat_AGN[2]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=AGNcol)  # median
                plt.plot([stat_AGN[4]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=AGNcol)  # 68% confidence interval lower
                plt.plot([stat_AGN[8]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=AGNcol)  # 68% confidence interval lower
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