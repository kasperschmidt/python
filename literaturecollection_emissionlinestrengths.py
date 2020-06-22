# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Collection of literature emission line measurements combined into single table
# Main catalog is generated with: literaturecollection_emissionlinestrengths.generate_literature_fitscatalog()
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import sys
import os
from fits2ascii import ascii2fits
import MiGs
from astropy import units as u
import astropy
import astropy.io.fits as afits
import astropy.coordinates as acoord
import astropy.cosmology as acosmo
import numpy as np
import collections
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import literaturecollection_emissionlinestrengths as lce
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def generate_literature_fitscatalog(verbose=True):
    """
    Collecting literature values and storing them in a single file

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce
    lce.generate_literature_fitscatalog()

    # updating plots
    lce.plot_literature_fitscatalog(showphotoionizationmodels=False,secondarydat_fits=None,logaxes=False,shownames=False)
    lce.plot_literature_fitscatalog_legend(legendshape=(15.0, 3.0),ncol=3,extra_textlist=['MUSE-Wide (Schmidt et al. 2020)','MUSE UDF mosaic (Schmidt et al. 2020)', 'MUSE UDF10 (Schmidt et al. 2020)'],extra_symlist=['o','D','X'],showkeynames=False)

    """
    outdir      = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/'
    outputfile  = outdir+'literaturecollection_emissionlinestrengths.txt'

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Initializing the output file: \n   '+outputfile)
    fout = open(outputfile,'w')
    fout.write('# Flux, EW and line ratios collected from the literature:\n')
    fout.write('# References in the column "reference" can be expanded with lce.referencedictionary()["reference"][1] \n')
    fout.write('# Upper and lower limits are given as values with uncertainty of +99 or -99, respectively. \n')
    fout.write('# \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Determine columns to fill in output')
    fluxratiodic = lce.build_master_datadictionary(verbose=False)
    Ncols        = len(fluxratiodic.keys())
    if verbose: print('   The output file will contain '+str(Ncols)+' columns ')
    fout.write('# This file contains the following '+str(Ncols)+' columns:\n')
    fout.write('# id                                  name   reference'+
               ' '.join([str("%20s" % colname) for colname in list(fluxratiodic.keys())[3:]])+'  \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Collecting literature data arrays and appending them')
    refdic                 = lce.referencedictionary()

    dataref, dataarray     = lce.data_nan19(verbose=True)
    outputdataarray        = dataarray
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_sch17(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_sen17(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_rig14(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_rig15(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_erb10(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_sta14(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_sta15(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_mai18(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_sha03(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_bay14(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_sch16(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_lef19(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_amo17(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    dataref, dataarray     = lce.data_ber19(verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0].decode('UTF-8')][1])

    if verbose: print('\n - Hence the total number of objects in the literture collection is '+str(len(outputdataarray['id'])))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Writing output data array to ascii file \n   '+outputfile)
    for oo, id in enumerate(outputdataarray['id']):
        outstr = str(int(id))+' '+str("%30s" % outputdataarray['name'][oo].decode('UTF-8'))+'   '+\
                 outputdataarray['reference'][oo].decode('UTF-8')+' '+\
                 str("%20.10f" % outputdataarray['ra'][oo])+' '+str("%20.10f" % outputdataarray['dec'][oo])+\
                 ' '.join([str("%20.4f" % ff) for ff in outputdataarray[oo].tolist()[5:]])
        fout.write(outstr+' \n')
    fout.close()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Creating fits version of output: \n   '+outputfile.replace('.txt','.fits'))
    fitsformat = ['K','30A','20A'] + ['D']*(Ncols-3)
    fitsoutput = ascii2fits(outputfile,asciinames=True,skip_header=5,fitsformat=fitsformat,verbose=False)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def emissionlinelist(verbose=True):
    """

    """
    linelistdic = MiGs.linelistdic(listversion='full')
    if verbose: print(' - Loaded the line dictionary from MiGs containing:\n   '+linelistdic.keys())
    return linelistdic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_reference_fromID(idlist,verbose=True):
    """
    Function to return the reference(s) given a list of literature IDs

    --- INPUT ---
    idlist       List of "literature ids". They should be 11 digits long where the first 3 digits refers
                 to the reference the IDs belong to as provided by lce.referencedictionary()
    verbose      Toggle verbosity of function

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce

    idlist     = [10000001024,10000001036,20022480005]
    references = lce.get_reference_fromID(idlist)
    """
    refdic     = lce.referencedictionary()
    returnrefs = []

    if type(idlist) is not list:
        idlist = [int(idlist)]

    if verbose: print(' - Looping over the '+str(len(idlist))+' IDs provide in list\n')
    for id in idlist:
        idstr  = str(id)
        if (len(idstr) < 10) or idstr.startswith('0'):
            if verbose: print(' - WARNING: literaturecollection_emissionlinestrenths.get_reference_fromID():\n'
                              '            The id '+idstr+' has less than 10 digitsis which \n            '
                                                          'does not match any object from the literature collection')
        else:
            baseid = int(idstr[:-9].ljust(len(idstr),'0'))

            foundref = False
            for key in refdic.keys():
                if int(refdic[key][0]) == baseid:
                    returnrefs.append([idstr,key,refdic[key][0],refdic[key][1],refdic[key][2]])
                    foundref = True

            if not foundref:
                if verbose: print(' - WARNING: literaturecollection_emissionlinestrenths.get_reference_fromID():\n'
                                  '            The id '+idstr+' has no match in literature collection. '
                                                              'Setting marker for plots to "$??$"')

                pdb.set_trace()
                returnrefs.append([idstr,'NoRef',baseid,'NoRef','$??$'])

    if len(idlist) != len(returnrefs):
        sys.exit('The provided ID list and the found references have different lengths:\n   len(idlist) = '+str(len(idlist))+':\n  '+str(idlist)+'\n   len(references) = '+str(len(returnrefs))+':\n   '+str(returnrefs)+'\n')
    if len(idlist) == 1:
        return returnrefs[0]
    else:
        return returnrefs
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def referencedictionary():
    """

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce
    refval = 'sch17'
    lce.referencedictionary()[refval][1]

    """
    refdic = collections.OrderedDict()
    #                  baseid   reference                                             plotsymbol
    refdic['nan19'] = [01e9,    'Nanaykkara et al. (2019)',                               '^']
    refdic['sch17'] = [02e9,    'Schmidt et al. (2017) & Mainali et al. (2017)',          '<']
    refdic['sen17'] = [03e9,    'Senchyna et al. (2017)',                                 'v']
    refdic['rig14'] = [04e9,    'Rigby et al. (2014)',                                    '>']
    refdic['rig15'] = [05e9,    'Rigby et al. (2015)',                                    '8']
    refdic['erb10'] = [06e9,    'Erb et al. (2010)',                                      's']
    refdic['sta14'] = [07e9,    'Stark et al. (2014)',                                    'p']
    refdic['sta15'] = [08e9,    'Stark et al. (2015a,b, 2017)',                           'P']
    refdic['mai18'] = [09e9,    'Mainali et al. (2018) & Stark et al. (2017)',            '*']
    refdic['sha03'] = [10e9,    'Shapley et al. (2003)',                                  'h']
    refdic['bay14'] = [11e9,    'Bayliss et al. (2014)',                                  'H']
    refdic['dummy'] = [12e9,    'dummy',                                                  '+']
    refdic['lef19'] = [13e9,    'Le Fevre et al. (2019)',                                 'x']
    refdic['ber19'] = [14e9,    'Berg et al. (2016, 2019a,b)',                            'D']
    refdic['dummy'] = [15e9,    'dummy',                                                  'd']
    refdic['dummy'] = [16e9,    'dummy',                                                  (4, 1, 0)]   # 4-point star
    refdic['dummy'] = [17e9,    'dummy',                                                  (4, 1, 45)]  # 4-point star
    refdic['sch16'] = [18e9,    'Schmidt et al. (2016)',                                  (7, 1, 0)]   # 7-point star
    refdic['dummy'] = [19e9,    'dummy',                                                  (7, 1, 90)]  # 7-point star
    refdic['dummy'] = [20e9,    'dummy',                                                  (7, 1, 180)] # 7-point star
    refdic['amo17'] = [51e9,    'Amorin et al. (2017)',                                   (5, 0, 180)] # pentagon rotated 180deg
    refdic['dummy'] = [52e9,    'dummy',                                                  (5, 0, 90)]  # pentagon rotated 90deg
    refdic['dummy'] = [53e9,    'dummy',                                                  (5, 0, 270)] # pentagon rotated 270deg

    # --- MUSE-Wide def: ---
    # CDFS and COSMOS:  'o'
    # UDF:              'D'
    # UDF10:            'X'

    # --- potentially outer symbols def: ---
    #  D    GLASS
    #  s    stack/composite
    #  *    AGN
    #  o    MUSE data

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Ravindranath et al. (2020) Green pea CIII + OIII + Lya emitters
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Senchyna et al. (2019); Rigby+18, Steidel+16, Christensen+12   <--- Byler+20
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Plat et al. (2019) compilation... check
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Huang et al. (2016b)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Malkan et al. (1996)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Smit et al. (2017)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Vanzella et al. (2016) - X-shooter + MUSE
    # Vanzella et al. (2017) - X-shooter + MUSE
    # Vanzella et al. (2020) - X-shooter Ion2 (SunburstArc comparison)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # SILVERRUSH - Shibuya+2017b collection of CIV, NV, HeII limits - no CIII limits. Table 6 in appendix
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Matthee+2017 compilation
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Nakajima+2017/18 Composite spectra
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Harikane+2017 have a predicted model for EWCIII vs EWLYA...
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Who+ (Stark talk at SakuraCLAW - Malhotra on it?)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Berg, Brammer 18: Lensed source
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Verhamme+18 collection of sys + Lya objects.
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Schaerer+18 single obj from Izotov with CIII
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Marchi+19 Vandels CIII vs Lya EWs
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Yang+19 Two z~7 LAGER confirmations with UV line upper limits
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Saxena+19 VANDELS HeII emitters
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # R. Marques-Chaves+19 UV of LAEs from the BELLS survey
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Alendroff+13 type 2 QSO candidates
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Hutchison+19 CIII at redshift 7.5 (Finkelstein's pet object)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Herenz+20 Lya blob with Lya, HeII and CIV(limit) detections at three locations
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Guo+20 MUSE AGN with CIV, CIII and HeII but only Ly fluxes provided - they show L(CIV) so maybe ask... or check when published
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


    return refdic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def lineinfo(verbose=True):
    """

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce
    linetypedic, linelistdic, componentlist = lineinfo()

    """
    linetypedic = collections.OrderedDict()
    linetypedic['Lya']   = 'single'
    #linetypedic['NV']    = 'doublet'
    linetypedic['CIV']   = 'doublet'
    linetypedic['HeII']  = 'single'
    linetypedic['OIII']  = 'doublet'
    linetypedic['SiIII'] = 'doublet'
    linetypedic['CIII']  = 'doublet'
    linetypedic['MgII']  = 'doublet'

    componentlist = []
    for line in linetypedic.keys():
        if linetypedic[line] == 'single':
            componentlist.append(line)
        elif linetypedic[line] == 'doublet':
            componentlist.append(line)
            componentlist.append(line+'1')
            componentlist.append(line+'2')

    linelistdic   = MiGs.linelistdic(listversion='full')

    return linetypedic, linelistdic, componentlist
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def set_ratios(numstr,denomstr,numerator,numeratorerr,denominator,denominatorerr):
    """

    Function to set the flux ratios of emission lines.
    Upper limits (denominator is Nsigma level) are returned with errors of +99
    Lower limits (numerator is Nsigma level) are returned with errors of -99

    --- INPUT ---
    numstr            String to indicate whether a reliable numerator exists,
                      i.e. if numstr='None' provide the Nsigma limit for the numerator.
                      If any other string is provided it is assumed that the numerator contains a reliable flux value.
    denomstr          Indicate whether a reliable denominator exists,
                      i.e. if tempstr='None' provide the Nsigma limit for the denominator
                      If any other string is provided it is assumed that the denominator contains a reliable flux value.
    numerator         Flux value of the numerator in the flux ratio
    numeratorerr      Uncertainty on the numerator
    denominator       Flux value of the denominator in the flux ratio
    denominatorerr    Uncertainty on the denominator

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce
    fluxratio, fluxratio_err = lce.set_ratios('None','good',22.,'dummy',150.,19.) # flux ratio upper limit
    fluxratio, fluxratio_err = lce.set_ratios('good','good',22.,3.,150.,5.)

    """
    ratio      = np.nan
    ratioerr   = np.nan

    if (numstr.lower() != 'none') & (denomstr.lower() != 'none' ):
        ratio      = numerator/denominator
        ratioerr   = np.sqrt( (numeratorerr/numerator)**2+(denominatorerr/denominator)**2) * np.abs(ratio)
    elif (numstr.lower() == 'none') & (denomstr.lower() != 'none' ):
        ratio      = numerator/denominator
        ratioerr   = +99
    elif (numstr.lower() != 'none') & (denomstr.lower() == 'none' ):
        ratio      = numerator/denominator
        ratioerr   = -99
    else:
        raise Exception(' Something went wrong in lce.set_ratios(); '
                        'inputs were: (numstr,denomstr,numerator,numeratorerr,denominator,denominatorerr)='+
                        str([numstr,denomstr,numerator,numeratorerr,denominator,denominatorerr]))

    return ratio, ratioerr
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def calc_doubletValuesFromSingleComponents(f1,f1err,f2,f2err,
                                           EW1=None,EW1err=None,EW2=None,EW2err=None,
                                           verbose=True):
    """
    Function to calculate the combined doublet flux from single measurements.

    --- EXAMPLE OF USE ---
    f_doublet, ferr_doublet, fratio, fratio_err = calc_doubletValuesFromSingleComponents(f1,f1err,f2,f2err)

    """
    f_doublet    = np.array([])
    ferr_doublet = np.array([])
    fratio       = np.array([])
    fratio_err   = np.array([])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for ff, f1_obj in enumerate(f1):
        if (f1err[ff] != 99) & (f2err[ff] != 99):
            f_obj      = f1[ff]+f2[ff]
            ferr_obj   = np.sqrt( f1err[ff]**2.0 + f2err[ff]**2.0 )
            FR, FR_err = lce.set_ratios('good','good',f1[ff],f1err[ff],f2[ff],f2err[ff])
        elif (f1err[ff] == 99) & (f2err[ff] != 99):
            f_obj      = f1[ff]+f2[ff]
            ferr_obj   = 99
            FR, FR_err = lce.set_ratios('None','good',f1[ff],f1err[ff],f2[ff],f2err[ff])
        elif (f1err[ff] != 99) & (f2err[ff] == 99):
            f_obj      = f1[ff]+f2[ff]
            ferr_obj   = 99
            FR, FR_err = lce.set_ratios('good','None',f1[ff],f1err[ff],f2[ff],f2err[ff])
        else:
            f_obj      = f1[ff]+f2[ff]
            ferr_obj   = 99
            FR, FR_err = np.nan, np.nan

        f_doublet    = np.append(f_doublet,f_obj)
        ferr_doublet = np.append(ferr_doublet,ferr_obj)
        fratio       = np.append(fratio,FR)
        fratio_err   = np.append(fratio_err,FR_err)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (EW1 is not None) & (EW2 is not None):
        EW     = np.array([])
        EWerr  = np.array([])
        for ee, EW1_obj in enumerate(EW1):

            if (np.abs(EW1err[ee]) != 99) & (np.abs(EW2err[ee]) != 99):
                cont1     = f1[ee]/EW1[ee]
                cont1err  = np.sqrt( (f1err[ee]/f1[ee])**2 + (EW1err[ee]/EW1[ee])**2 ) * np.abs(cont1)

                cont2     = f2[ee]/EW2[ee]
                cont2err  = np.sqrt( (f2err[ee]/f2[ee])**2 + (EW2err[ee]/EW2[ee])**2 ) * np.abs(cont2)

                cont      = (cont1+cont2) / 2.
                conterr   = np.sqrt(cont1err**2.0 + cont2err**2.0)

                if (ferr_doublet[ee] == 99):
                    EW_obj    = f_doublet[ee] / cont
                    EWerr_obj = 99
                else:
                    EW_obj    = f_doublet[ee] / cont
                    EWerr_obj = np.sqrt( (ferr_doublet[ee]/f_doublet[ee])**2 + (conterr/cont)**2 ) * np.abs(EW_obj)

            elif (np.abs(EW1err[ee]) == 99) & (np.abs(EW2err[ee]) != 99):
                cont      = f2[ee]/EW2[ee]

                if (ferr_doublet[ee] == 99) & (EW2err[ee] == -99):
                    EW_obj    = f_doublet[ee] / cont
                    EWerr_obj = 99
                elif (ferr_doublet[ee] == 99) & (EW2err[ee] == 99):
                    EW_obj    = np.nan
                    EWerr_obj = np.nan
                else:
                    EW_obj    = f_doublet[ee] / cont
                    EWerr_obj = EW1err[ee]

            elif (np.abs(EW1err[ee]) != 99) & (np.abs(EW2err[ee]) == 99):
                cont      = f1[ee]/EW1[ee]

                if (ferr_doublet[ee] == 99) & (EW1err[ee] == -99):
                    EW_obj    = f_doublet[ee] / cont
                    EWerr_obj = 99
                elif (ferr_doublet[ee] == 99) & (EW1err[ee] == 99):
                    EW_obj    = np.nan
                    EWerr_obj = np.nan
                else:
                    EW_obj    = f_doublet[ee] / cont
                    EWerr_obj = EW2err[ee]

            elif (EW1err[ee] == -99) & (EW2err[ee] == -99):
                cont      = (f1[ee]/EW1[ee] + f2[ee]/EW2[ee]) / 2.
                EW_obj    = f_doublet[ee] / cont
                EWerr_obj = -99
            else:
                EW_obj     = np.nan
                EWerr_obj  = np.nan

            EW        = np.append(EW,EW_obj)
            EWerr     = np.append(EWerr,EWerr_obj)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (EW1 is not None) & (EW2 is not None):
        return f_doublet, ferr_doublet, fratio, fratio_err, EW, EWerr
    else:
        return f_doublet, ferr_doublet, fratio, fratio_err
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_dataarray(catreference, datadic, S2Nlim=np.nan, verbose=True):
    """
    Function building data array from literature data

    --- INPUT ---
    ids             The ids of the objects
    redshifts       The redshifts of the objects
    catreference    Catalog reference (should exist in lce.referencedictionary())
    datadic         Dictionary containing data to fill into array
    S2Nlim          Sigma of limits in datadic.
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---
    dataarray = lce.build_dataarray(catreference, linetypedic, datadic) # see example in lce.data_sen17()

    """
    linetypedic, linelistdic, componentlist = lineinfo(verbose=verbose)
    Ndataobj      = len(datadic['id'])

    masterdic     = lce.build_master_datadictionary(verbose=verbose)
    columns       = masterdic.keys()

    dtypebuild    = [(colname, 'd') for colname in columns]
    dtypebuild[1] = ('name', '30a')
    dtypebuild[2] = ('reference', '20a')
    dataarray     = np.array(np.zeros(Ndataobj)*np.nan,dtype=dtypebuild)

    if verbose: print('   Inserting data from literature into master dictionary ')
    for ll, datacol in enumerate(datadic.keys()):
        if datacol in masterdic.keys():
            dataarray[datacol] = datadic[datacol]
        else:
            print('   WARNING: The data dictionary column "'+datacol+'" not found in master dic. so not stored in outut.')

    if verbose: print('   Estimating further entries based on data from literature')
    for ii, id in enumerate(dataarray['id']):
        for ll, numerator_line in enumerate(componentlist):
            fnum    = dataarray['f_'+numerator_line][ii]
            ferrnum = dataarray['ferr_'+numerator_line][ii]

            if ferrnum == 99:
                dataarray['s2n_'+numerator_line][ii]     = S2Nlim
            else:
                dataarray['s2n_'+numerator_line][ii]     = fnum / ferrnum

            for kk, denominator_line in enumerate(componentlist):
                fdenom    = dataarray['f_'+denominator_line][ii]
                ferrdenom = dataarray['ferr_'+denominator_line][ii]

                if numerator_line == denominator_line:
                    continue
                elif numerator_line in denominator_line:
                    continue
                elif denominator_line in numerator_line:
                    continue
                elif np.isnan(fnum) or np.isnan(fdenom):
                    continue
                else:
                    if (ferrnum == 99) & (ferrdenom == 99):
                        FR, FRerr = np.nan, np.nan
                    elif (ferrnum == 99) & (ferrdenom != 99):
                        FR, FRerr = lce.set_ratios('None','denominator_exists',
                                                    fnum,  'dummy', fdenom, ferrdenom)
                    elif (ferrnum != 99) & (ferrdenom == 99):
                        FR, FRerr = lce.set_ratios('numerator_exists','none',
                                                    fnum,  ferrnum, fdenom, 'dummy')
                    else:
                        FR, FRerr = lce.set_ratios('numerator_exists','denominator_exists',
                                                    fnum,  ferrnum, fdenom, ferrdenom)

                    rationame         = numerator_line+denominator_line
                    rationame_inverse = denominator_line+numerator_line
                    if 'FR_'+rationame_inverse not in dataarray.dtype.names: # and ratios withput inverse versions

                        # only attempt to fill empty columns to avoid overwriting values from input
                        if np.isnan(dataarray['FR_'+numerator_line+denominator_line][ii]) & \
                                np.isnan(dataarray['FRerr_'+numerator_line+denominator_line][ii]):

                            dataarray['FR_'+rationame][ii]     = FR
                            dataarray['FRerr_'+rationame][ii]  = FRerr
                            if np.abs(FRerr) == 99:
                                dataarray['FRs2n_'+rationame][ii]  = S2Nlim
                            else:
                                dataarray['FRs2n_'+rationame][ii]  = FR/FRerr

    return dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_master_datadictionary(verbose=True):
    """
    Defining the master dictionary to fill with data from the literature
    """
    datadictionary = collections.OrderedDict()
    datadictionary['id']        = np.array([])
    datadictionary['name']      = np.array([])
    datadictionary['reference'] = np.array([])
    datadictionary['ra']        = np.array([])
    datadictionary['dec']       = np.array([])
    datadictionary['redshift']  = np.array([])

    linetypedic, linelistdic, componentlist = lineinfo(verbose=verbose)

    for ll, numerator_line in enumerate(componentlist):
        datadictionary['f_'+numerator_line]       = np.array([])
        datadictionary['ferr_'+numerator_line]    = np.array([])
        datadictionary['s2n_'+numerator_line]     = np.array([])
        datadictionary['EW0_'+numerator_line]      = np.array([])
        datadictionary['EW0err_'+numerator_line]   = np.array([])
        datadictionary['sigma_'+numerator_line]   = np.array([])     # FWHM     / 2.355
        datadictionary['sigmaerr_'+numerator_line]   = np.array([])  # FWHM_err / 2.355
        datadictionary['vshift_'+numerator_line]  = np.array([])
        datadictionary['vshifterr_'+numerator_line]  = np.array([])

        for kk, denominator_line in enumerate(componentlist):
            if numerator_line == denominator_line:
                continue
            elif numerator_line in denominator_line:
                continue
            elif denominator_line in numerator_line:
                continue
            else:
                rationame         = numerator_line+denominator_line
                rationame_inverse = denominator_line+numerator_line

                if 'FR_'+rationame_inverse not in datadictionary.keys():
                    datadictionary['FR_'+rationame]     = np.array([])
                    datadictionary['FRerr_'+rationame]  = np.array([])
                    datadictionary['FRs2n_'+rationame]  = np.array([])

    if verbose: print(' - Assembled and build a master dictionary with '+str(len(datadictionary.keys()))+' keys')
    return datadictionary
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_literature_fitscatalog(secondarydat_fits=None,logaxes=True, shownames=False,
                                showphotoionizationmodels=True,overwrite=True,verbose=True):
    """
    Wrapper for the plot_literature_fitscatalog_cmd() defining the samples

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce

    secondarydat_fits = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/results_master_catalog_version191128.fits'

    lce.plot_literature_fitscatalog(showphotoionizationmodels=False,secondarydat_fits=secondarydat_fits)
    lce.plot_literature_fitscatalog(showphotoionizationmodels=True)

    """
    maindir = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/'
    litdat  = afits.open(maindir+'/literaturecollection_emissionlinestrengths.fits')[1].data

    if secondarydat_fits is not None:
        dat_2nd = afits.open(secondarydat_fits)[1].data

        file_info = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
        dat_info  = afits.open(file_info)[1].data

        z_2nd   = []
        for id in dat_2nd['id']:
            objent = np.where(dat_info['id'] == id)[0]
            if len(objent) == 1:
                z_2nd.append(dat_info['redshift'][objent[0]])
            else:
                z_2nd.append(0.0)
        z_2nd = np.asarray(z_2nd)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # generate vector of plotting symbols
    refdic = lce.referencedictionary()
    psyms  = []
    for objref in litdat['reference']:
        psyms.append(refdic[objref][2])
    # psyms = np.asarray(psyms)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fluxes_range  = None #[1,1e8]
    EW0_range     = None
    ratios_range  = None #[1e-4,1e3]

    input_lists = []

    # ------------------------------------ Flux vs Flux plots ------------------------------------
    input_lists.append(['f_Lya','f_CIII','f(Ly$\\alpha$) [1e-20 erg/s/cm$^2$/\AA]', 'f(CIII) [1e-20 erg/s/cm$^2$/\AA]',
                        fluxes_range, fluxes_range, 'redshift',  ''])
    input_lists.append(['f_Lya','f_CIV', 'f(Ly$\\alpha$) [1e-20 erg/s/cm$^2$/\AA]', 'f(CIV) [1e-20 erg/s/cm$^2$/\AA]',
                        fluxes_range, fluxes_range, 'redshift',  ''])
    input_lists.append(['f_CIII','f_CIV','f(CIII) [1e-20 erg/s/cm$^2$/\AA]', 'f(CIV) [1e-20 erg/s/cm$^2$/\AA]',
                        fluxes_range, fluxes_range, 'redshift',  ''])
    input_lists.append(['f_OIII','f_HeII','f(OIII) [1e-20 erg/s/cm$^2$/\AA]', 'f(HeII) [1e-20 erg/s/cm$^2$/\AA]',
                        fluxes_range, fluxes_range, 'redshift',  ''])

    # ------------------------------------ EW vs EW plots ------------------------------------
    input_lists.append(['EW0_Lya','EW0_CIII','EW$_0$(Ly$\\alpha$) [\AA]', 'EW$_0$(CIII) [\AA]',
                        EW0_range, EW0_range, 'redshift',  ''])
    input_lists.append(['EW0_Lya','EW0_CIV', 'EW$_0$(Ly$\\alpha$) [\AA]', 'EW$_0$(CIV) [\AA]',
                        EW0_range, EW0_range, 'redshift',  ''])
    input_lists.append(['EW0_CIII','EW0_CIV','EW$_0$(CIII) [\AA]', 'EW$_0$(CIV) [\AA]',
                        EW0_range, EW0_range, 'redshift',  ''])
    input_lists.append(['EW0_OIII','EW0_HeII','EW$_0$(OIII) [\AA]', 'EW$_0$(HeII) [\AA]',
                        EW0_range, EW0_range, 'redshift',  ''])

    # ------------------------------------ Flux artio vs Flux ratio plots ------------------------------------
    input_lists.append(['FR_CIVCIII','FR_CIVHeII','CIV/CIII','CIV/HeII',
                        ratios_range,ratios_range  , 'redshift', 'Schmidt+17 fig. 7 top,   Feltre+16 fig A2a'])
    input_lists.append(['FR_OIIICIII','FR_HeIICIII','OIII/CIII','CIII/HeII',
                        ratios_range,ratios_range  , 'redshift', ' No Title '])

    # ------------------------------------ Handling data and plotting ------------------------------------
    Nhistbins   = 30
    histaxes    = False
    for il in input_lists:
        xval   = litdat[il[0]]
        yval   = litdat[il[1]]
        xerr   = litdat[il[0].replace('_','err_')]
        yerr   = litdat[il[1].replace('_','err_')]
        colval = litdat[il[6]]

        plotname = maindir+'plots/literaturecollection_emissionlinestrengths_plot_'+il[0]+'_VS_'+il[1]+'.pdf'

        if secondarydat_fits is None:
            xval_2nd   = None
            yval_2nd   = None
            xerr_2nd   = None
            yerr_2nd   = None
            colval_2nd = None
        else:
            xval_2nd   = dat_2nd[il[0]]
            yval_2nd   = dat_2nd[il[1]]
            xerr_2nd   = dat_2nd[il[0].replace('_','err_')]
            yerr_2nd   = dat_2nd[il[1].replace('_','err_')]
            colval_2nd = z_2nd

        if showphotoionizationmodels:
            if ('FR_' in il[0]) & ('FR_' in il[1]):
                #x2plot, y2plot, varyparam, cutSFmodels, markersize, SFmarker, AGNmarker = piplotparam
                photoionizationplotparam = lce.colname2NEOGAL(il[0]),lce.colname2NEOGAL(il[1]),\
                                           'Zgas', False, 1.5, 's', 'D'
            else:
                photoionizationplotparam = None
                if verbose: print(' - WARNING Not plotting photoionization models as flux units of these not well-defined. '
                                  'Focus on ratios.')
        else:
            photoionizationplotparam = None

        if photoionizationplotparam is None:
            cmapselected = 'viridis_r' #'autumn_r'
        else:
            cmapselected = 'summer_r'  # 'Greens',#'autumn_r'

        if shownames:
            point_text = litdat['name']
        else:
            point_text = None

        if ('EW' in il[2]) & ('EW' in il[3]):
            drawcurves = ['onetoone','onetothree','threetoone','onetoten','tentoone']
        else:
            drawcurves = ['onetoone']

        lce.plot_literature_fitscatalog_cmd(plotname,
                                            # - - - - - - - Literature Data Setup - - - - - - -
                                            xval=xval,yval=yval,xerr=xerr,yerr=yerr,
                                            colval=colval,psym=psyms,fixcolor='gray',
                                            point_text=point_text,ids=None,
                                            # - - - - - - - Secondary Data Setup - - - - - - -
                                            xval_2nd=xval_2nd,yval_2nd=yval_2nd,xerr_2nd=xerr_2nd,yerr_2nd=yerr_2nd,
                                            colval_2nd=colval_2nd,psym_2nd='o',fixcolor_2nd='red',
                                            point_text_2nd=None,ids_2nd=None,
                                            # - - - - - - - Coloring Setup - - - - - - -
                                            colortype='redshift',colmap=cmapselected,
                                            # - - - - - - - Plot Setup - - - - - - -
                                            drawcurves=drawcurves,histaxes=histaxes,Nbins=Nhistbins,
                                            photoionizationplotparam=photoionizationplotparam,
                                            xlabel=il[2],ylabel=il[3],yrange=il[4],xrange=il[5],ylog=logaxes,xlog=logaxes,
                                            title=il[7],overwrite=overwrite,verbose=verbose)

    # lce.plot_literature_fitscatalog_cmd(maindir+'plots/test.pdf',overwrite=overwrite) #empty

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_literature_fitscatalog_cmd(plotname,
                                    # - - - - - - - Literature Data Setup - - - - - - -
                                    xval=None,yval=None,xerr=None,yerr=None,
                                    colval=None,psym='s',fixcolor='gray',point_text=None,ids=None,
                                    # - - - - - - - Secondary Data Setup - - - - - - -
                                    xval_2nd=None,yval_2nd=None,xerr_2nd=None,yerr_2nd=None,
                                    colval_2nd=None,psym_2nd='o',fixcolor_2nd='red',point_text_2nd=None,ids_2nd=None,
                                    # - - - - - - - Coloring Setup - - - - - - -
                                    colortype=None,colmap='viridis_r', # 'autumn_r'
                                    # - - - - - - - Plot Setup - - - - - - -
                                    drawcurves=['onetoone'],histaxes=False,Nbins=50,
                                    photoionizationplotparam=None,
                                    xlabel='x-axis',ylabel='y-axis',yrange=None,xrange=None,ylog=False,xlog=False,
                                    title=None,overwrite=False,verbose=True):
    """
    Collecting literature values and storing them in a single file
    (initially based on uvEmissionlineSearch.plot_mocspecFELISresults_summary_plotcmds() )

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce
    lce.plot_literature_fitscatalog_cmd()

    """
    if verbose: print(' - Setting up and generating plot ')
    if xval is None:
        xval   = np.array([0,0])
        yval   = np.array([0,0])
        xerr   = np.array([np.nan,np.nan])
        yerr   = np.array([np.nan,np.nan])

    # - - - - - - LITERATURE DATA - - - - -
    if (type(psym) == str):
        psym = np.asarray([psym]*len(xval))

    if ylog & xlog:
        goodent = np.where(np.isfinite(xval) & np.isfinite(yval) & (yval > 0) & (xval > 0))[0]
    elif xlog & (not xlog):
        goodent = np.where(np.isfinite(xval) & np.isfinite(yval) & (xval > 0))[0]
    elif (not xlog) & ylog:
        goodent = np.where(np.isfinite(xval) & np.isfinite(yval) & (yval > 0))[0]
    else:
        goodent = np.where(np.isfinite(xval) & np.isfinite(yval))[0]

    if len(goodent) == 0:
        if verbose: print('WARNING: Did not find any good entries in main data, i.e. where values are finite and >0 if log axes chosen')
        psym = ['.','.']
        xval = np.array([0,0])
        yval = np.array([0,0])
        xerr = np.array([np.nan,np.nan])
        yerr = np.array([np.nan,np.nan])
    else:
        psym   = [psym[ii] for ii in goodent]
        xval   = xval[goodent]
        yval   = yval[goodent]
        xerr   = xerr[goodent]
        yerr   = yerr[goodent]
        if colval is not None:
            colval = colval[goodent]
        if point_text is not None:
            point_text = point_text[goodent]
    # - - - - - - SECONDARY DATA - - - - -
    if xval_2nd is not None:
        if (type(psym_2nd) == str) & (xval_2nd is not None):
            psym_2nd = np.asarray([psym_2nd]*len(xval_2nd))

        if ylog & xlog:
            goodent_2nd = np.where(np.isfinite(xval_2nd) & np.isfinite(yval_2nd) & (yval_2nd > 0) & (xval_2nd > 0))[0]
        elif xlog & (not xlog):
            goodent_2nd = np.where(np.isfinite(xval_2nd) & np.isfinite(yval_2nd) & (xval_2nd > 0))[0]
        elif (not xlog) & ylog:
            goodent_2nd = np.where(np.isfinite(xval_2nd) & np.isfinite(yval_2nd) & (yval_2nd > 0))[0]
        else:
            goodent_2nd = np.where(np.isfinite(xval_2nd) & np.isfinite(yval_2nd))[0]

        if len(goodent_2nd) == 0:
            if verbose: print('WARNING: Did not find any good entries in secondary data, i.e. where values are finite and >0 if log axes chosen')
            psym = ['.','.']
            xval_2nd   = np.array([0,0])
            yval_2nd   = np.array([0,0])
            xerr_2nd   = np.array([np.nan,np.nan])
            yerr_2nd   = np.array([np.nan,np.nan])
        else:
            psym_2nd   = psym_2nd[goodent_2nd]
            xval_2nd   = xval_2nd[goodent_2nd]
            yval_2nd   = yval_2nd[goodent_2nd]
            xerr_2nd   = xerr_2nd[goodent_2nd]
            yerr_2nd   = yerr_2nd[goodent_2nd]

            if colval_2nd is not None:
                colval_2nd = colval_2nd[goodent_2nd]

            if point_text_2nd is not None:
                point_text_2nd = point_text_2nd[goodent_2nd]

    # - - - - - - SETTING UP PLOT - - - - -
    if os.path.isfile(plotname) & (not overwrite):
        if verbose: print('\n - WARNING: the plot '+plotname+' exists and overwrite=False so moving on \n')
    else:
        if histaxes:
            fig = plt.figure(1, figsize=(6, 6))
            left, width = 0.15, 0.60
            bottom, height = 0.15, 0.60
            bottom_h = left_h = left + width + 0.01

            fig.subplots_adjust(wspace=0.1, hspace=0.1,left=left, right=left+width, bottom=bottom, top=bottom+height)
            rect_histx = [left, bottom_h, width, 0.2]
            rect_histy = [left_h, bottom, 0.2, height]
        else:
            fig = plt.figure(2, figsize=(6, 5))
            fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.97, bottom=0.15, top=0.95)
            if photoionizationplotparam:
                fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.10, right=1.2, bottom=0.10, top=0.95)

        Fsize         = 14.0
        lthick        = 2.0
        marksize      = 6.0
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()

        #---------  setting up colors for litearture data ---------
        if (colval is not None) or (colval_2nd is not None):
            cmap    = plt.cm.get_cmap(colmap) # 'autumn_r','viridis_r'

            if (colval is not None) & (colval_2nd is None):
                cdatvec = colval
            elif (colval is None) & (colval_2nd is not None):
                cdatvec = colval_2nd
            else:
                cdatvec = np.append(colval,colval_2nd)

            if colortype.lower() == 'redshift':
                clabel  = '$z$'
                cmin    = 0.0 # 0.0, 1.4
                cmax    = 9.0 # 10.2, 6.2
                cextend = 'neither'
            elif colortype.lower() == 's2n':
                clabel  = 'S/N'
                cmin    = 3.0
                cmax    = 10.0
                cextend = 'both'
            elif colortype.lower() == 'vshift':
                clabel  = 'Velocity shift (spec vs. template match) [km/s]'
                cmin    = 0.0
                cmax    = 300.0
                cextend = 'both'
            elif colortype.lower() == 'ew0lya':
                clabel  = 'EW$_0$(Ly$\\alpha$) [\AA]'
                cmin    = np.min(cdatvec[np.isfinite(cdatvec)])
                cmax    = np.max(cdatvec[np.isfinite(cdatvec)])
                cextend = 'neither'
            elif colortype.lower() == 'sigma':
                clabel  = '$\sigma$ [\AA]'
                cmin    = np.min(cdatvec[np.isfinite(cdatvec)])
                cmax    = np.max(cdatvec[np.isfinite(cdatvec)])
                cextend = 'neither'
            else:
                cmin    = np.min(cdatvec[np.isfinite(cdatvec)])
                cmax    = np.max(cdatvec[np.isfinite(cdatvec)])
                cextend = 'neither'

            colnorm = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
            cmaparr = np.linspace(cmin, cmax, num=50)
            m       = plt.cm.ScalarMappable(cmap=cmap)
            m.set_array(cmaparr)

            if not histaxes:
                colshrink = 1.0
                colaspect = 30
                if photoionizationplotparam is None:
                    colanchor = (0.0,0.5)
                else:
                    colbarscale = 2.1
                    colanchor   = (-1.8,0.0)
                    colshrink   = colshrink/colbarscale
                    colaspect   = colaspect/colbarscale

                cb      = plt.colorbar(m,extend=cextend,orientation='vertical',
                                       pad=0.01,aspect=colaspect,shrink=colshrink,anchor=colanchor,use_gridspec=False)
                cb.set_label(clabel)

            if (colval is not None):
                colvec   = []
                for ii,xv in enumerate(xval):
                    colvec.append(cmap(colnorm(colval[ii])))
                facecol  = colvec
            if (colval_2nd is not None):
                colvec_2nd  = []
                for ii,xv in enumerate(xval_2nd):
                    colvec_2nd.append(cmap(colnorm(colval_2nd[ii])))
                facecol_2nd  = colvec_2nd

        else:
            facecol  = [fixcolor]*len(xval)

        #--------- RANGES ---------
        if not xrange:
            if xlog:
                xmin   = np.min(xval[np.isfinite(xval) & (xval > 0)])
                xmax   = np.max(xval[np.isfinite(xval) & (xval > 0)])
                xrange = [float(str("%.e" % xmin)),float(str("%.e" % xmax))*10.]
            else:
                xmin   = np.min(xval[np.isfinite(xval)])
                xmax   = np.max(xval[np.isfinite(xval)])
                dx     = xmax-xmin
                xrange = [xmin-dx*0.05,xmax+dx*0.05]

        plt.xlim(xrange)
        xminsys, xmaxsys = plt.xlim() # use to get automatically expanded axes if xmin = xmax

        if not yrange:
            if ylog:
                ymin   = np.min(yval[np.isfinite(yval)])
                ymax   = np.max(yval[np.isfinite(yval)])
                yrange = [float(str("%.e" % ymin)),float(str("%.e" % ymax))*10.]
            else:
                ymin   = np.min(yval[np.isfinite(yval)])
                ymax   = np.max(yval[np.isfinite(yval)])
                dy     = ymax-ymin
                yrange = [ymin-dy*0.05,ymax+dy*0.05]

        plt.ylim(yrange)
        yminsys, ymaxsys = plt.ylim() # use to get automatically expanded axes if xmin = xmax

        #--------- Plot Literature data with X and Y limits ---------
        ms          = marksize
        limsizefrac = 0.05
        y_uplimarr  = (np.asarray(yerr).astype(int) == +99)
        y_lolimarr  = (np.asarray(yerr).astype(int) == -99)
        x_uplimarr  = (np.asarray(xerr).astype(int) == +99)
        x_lolimarr  = (np.asarray(xerr).astype(int) == -99)
        for ii,xv in enumerate(xval): # loop necessary for coloring and upper/lower limits markers

            if y_uplimarr[ii]:
                if ylog:
                    dlog     = np.abs(np.diff(np.log10(plt.ylim()))) * limsizefrac
                    yerr[ii] = np.abs(yval[ii] - 10.**(np.log10(yval[ii])-dlog))
                else:
                    yerr[ii] = np.abs(np.diff(plt.ylim())) * limsizefrac
            if y_lolimarr[ii]:
                if ylog:
                    dlog     = np.abs(np.diff(np.log10(plt.ylim()))) * limsizefrac
                    yerr[ii] = np.abs(yval[ii] - 10.**(np.log10(yval[ii])+dlog))
                else:
                    yerr[ii] = np.abs(np.diff(plt.ylim())) * limsizefrac
            if x_uplimarr[ii]:
                if xlog:
                    dlog     = np.abs(np.diff(np.log10(plt.xlim()))) * limsizefrac
                    xerr[ii] = np.abs(xval[ii] - 10.**(np.log10(xval[ii])-dlog))
                else:
                    xerr[ii] = np.abs(np.diff(plt.xlim())) * limsizefrac
            if x_lolimarr[ii]:
                if xlog:
                    dlog     = np.abs(np.diff(np.log10(plt.xlim()))) * limsizefrac
                    xerr[ii] = np.abs(xval[ii] - 10.**(np.log10(xval[ii])+dlog))
                else:
                    xerr[ii] = np.abs(np.diff(plt.xlim())) * limsizefrac

            if y_uplimarr[ii] or y_lolimarr[ii] or x_uplimarr[ii] or x_lolimarr[ii]:
                mfc      = facecol[ii]
                mec      = 'gray'
                ecolor   = 'gray'
                zorder   = 9
            else:
                mfc      = facecol[ii]
                mec      = facecol[ii]
                ecolor   = facecol[ii]
                zorder   = 10

            plt.errorbar(xval[ii],yval[ii],xerr=xerr[ii],yerr=yerr[ii],capthick=0.5,
                         uplims=y_uplimarr[ii],lolims=y_lolimarr[ii],xuplims=x_uplimarr[ii],xlolims=x_lolimarr[ii],
                         marker=psym[ii],lw=lthick/2., markersize=ms,alpha=1.0,
                         markerfacecolor=mfc,ecolor=ecolor,
                         markeredgecolor=mec,zorder=zorder)

            if point_text is not None:
                plt.text(xval[ii]*1.03,yval[ii]*1.03,
                         point_text[ii].replace('_','\_'),color='white',fontsize=Fsize*0.5,zorder=30,
                         bbox=dict(boxstyle="round",edgecolor='k',facecolor=facecol[ii],linewidth=lthick*0.2))

        #--------- Plot secondary data with X and Y limits ---------
        if xval_2nd is not None:
            y_uplimarr  = (np.asarray(yerr_2nd).astype(int) == +99)
            y_lolimarr  = (np.asarray(yerr_2nd).astype(int) == -99)
            x_uplimarr  = (np.asarray(xerr_2nd).astype(int) == +99)
            x_lolimarr  = (np.asarray(xerr_2nd).astype(int) == -99)
            for ii,xv in enumerate(xval_2nd): # loop necessary for coloring and upper/lower limits markers

                if y_uplimarr[ii]:
                    if ylog:
                        dlog     = np.abs(np.diff(np.log10(plt.ylim()))) * limsizefrac
                        yerr_2nd[ii] = np.abs(yval_2nd[ii] - 10.**(np.log10(yval_2nd[ii])-dlog))
                    else:
                        yerr_2nd[ii] = np.abs(np.diff(plt.ylim())) * limsizefrac
                if y_lolimarr[ii]:
                    if ylog:
                        dlog     = np.abs(np.diff(np.log10(plt.ylim()))) * limsizefrac
                        yerr_2nd[ii] = np.abs(yval_2nd[ii] - 10.**(np.log10(yval_2nd[ii])+dlog))
                    else:
                        yerr_2nd[ii] = np.abs(np.diff(plt.ylim())) * limsizefrac
                if x_uplimarr[ii]:
                    if xlog:
                        dlog     = np.abs(np.diff(np.log10(plt.xlim()))) * limsizefrac
                        xerr_2nd[ii] = np.abs(xval_2nd[ii] - 10.**(np.log10(xval_2nd[ii])-dlog))
                    else:
                        xerr_2nd[ii] = np.abs(np.diff(plt.xlim())) * limsizefrac
                if x_lolimarr[ii]:
                    if xlog:
                        dlog     = np.abs(np.diff(np.log10(plt.xlim()))) * limsizefrac
                        xerr_2nd[ii] = np.abs(xval_2nd[ii] - 10.**(np.log10(xval_2nd[ii])+dlog))
                    else:
                        xerr_2nd[ii] = np.abs(np.diff(plt.xlim())) * limsizefrac

                plt.errorbar(xval_2nd[ii],yval_2nd[ii],xerr=xerr_2nd[ii],yerr=yerr_2nd[ii],capthick=0.5,
                             uplims=y_uplimarr[ii],lolims=y_lolimarr[ii],xuplims=x_uplimarr[ii],xlolims=x_lolimarr[ii],
                             marker=psym_2nd[ii],lw=lthick/2., markersize=ms,alpha=1.0,
                             markerfacecolor=facecol_2nd[ii],ecolor=facecol_2nd[ii],
                             markeredgecolor=facecol_2nd[ii],zorder=20)

                if point_text_2nd is not None:
                    plt.text(xval_2nd[ii]*1.03,yval_2nd[ii]*1.03,
                             point_text_2nd[ii].replace('_','\_'),color='white',fontsize=Fsize*0.5,zorder=30,
                             bbox=dict(boxstyle="round",edgecolor='k',facecolor=facecol_2nd[ii],linewidth=lthick*0.2))


        #--------- Draw reference curves on plot ---------
        if drawcurves is not None:
            tot_maxval = np.max([xmaxsys, ymaxsys])
            tot_minval = np.min([xminsys, yminsys])
            for dc in drawcurves:
                if dc == 'zero_horizontal':
                    plt.plot([xminsys,xmaxsys],[0,0],'--',color='black',lw=lthick,zorder=10)
                elif dc == 'zero_vertical':
                    plt.plot([0,0],[yminsys,ymaxsys],'--',color='black',lw=lthick,zorder=10)
                elif dc == 'onetoone':
                    plt.plot([tot_minval,tot_maxval],[tot_minval,tot_maxval],'-',color='black',lw=lthick,zorder=10)
                elif dc == 'tentoone':
                    plt.plot([tot_minval,tot_maxval],[tot_minval/10.,tot_maxval/10.],':',color='black',lw=lthick,zorder=10)
                elif dc == 'onetoten':
                    plt.plot([tot_minval,tot_maxval],[tot_minval*10,tot_maxval*10],':',color='black',lw=lthick,zorder=10)
                elif dc == 'threetoone':
                    plt.plot([tot_minval,tot_maxval],[tot_minval/3.,tot_maxval/3.],'--',color='black',lw=lthick,zorder=10)
                elif dc == 'onetothree':
                    plt.plot([tot_minval,tot_maxval],[tot_minval*3,tot_maxval*3],'--',color='black',lw=lthick,zorder=10)
                elif dc == 'plus':
                    plt.plot([xminsys,xmaxsys],[0,0],'--',color='black',lw=lthick,zorder=10)
                    plt.plot([0,0],[yminsys,ymaxsys],'--',color='black',lw=lthick,zorder=10)
                else:
                    sys.exit(' No setup for the string "'+str(dc)+'" in the list "drawcurves" so nothing drawn')

        #--------- Add Photoionization Grids ---------
        if photoionizationplotparam is not None:
            plotname      = plotname.replace('.pdf','_wPhotIoModels.pdf')
            titleaddition = lce.add_photoionization_models_to_plot(photoionizationplotparam)
        else:
            titleaddition = ''

        if (title is not None) & (histaxes == False):
            plt.title(title+titleaddition,fontsize=Fsize-4)

        #--------- Plot setup ---------
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if ylog:
            plt.yscale('log')
        if xlog:
            plt.xscale('log')
        #--------- Add hist axes ---------
        if histaxes:
            axHistx = plt.axes(rect_histx)
            if (title is not None):
                plt.title(title+titleaddition,fontsize=Fsize-4)

            axHisty = plt.axes(rect_histy)

            axHistx.xaxis.set_major_formatter(NullFormatter())
            axHisty.yaxis.set_major_formatter(NullFormatter())

            binwidth_x = np.diff([xminsys,xmaxsys])/Nbins
            bindefs    = np.arange(xminsys, xmaxsys+binwidth_x, binwidth_x)
            if xlog:
                bindefs = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))
                axHistx.set_xscale('log')

            if xval is not None:
                axHistx.hist(xval, bins=bindefs,histtype='step',color=fixcolor)
            if xval_2nd is not None:
                axHistx.hist(xval_2nd, bins=bindefs,histtype='step',color=fixcolor_2nd)

            axHistx.set_xticks([])
            axHistx.set_xlim([xminsys,xmaxsys])

            binwidth_y = np.diff([yminsys,ymaxsys])/Nbins
            bindefs    = np.arange(yminsys, ymaxsys+binwidth_y, binwidth_y)
            if ylog:
                bindefs = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))
                axHisty.set_yscale('log')

            if yval is not None:
                axHisty.hist(yval, bins=bindefs,histtype='step',color=fixcolor, orientation='horizontal')
            if yval_2nd is not None:
                axHisty.hist(yval_2nd, bins=bindefs,histtype='step',color=fixcolor_2nd, orientation='horizontal')
            axHisty.set_yticks([])
            axHisty.set_ylim([yminsys,ymaxsys])

            if colval is not None:
                cb      = plt.colorbar(m,extend=cextend,orientation='vertical',
                                       pad=0.01,aspect=10,shrink=0.35,anchor=(-15.0,1.58),use_gridspec=False)
                cb.set_label(clabel)
        else:
            pass

    if verbose: print('   Saving plot to \n   '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_literature_fitscatalog_legend(legendshape=(10, 3),ncol=5,extra_textlist=[],extra_symlist=[],showkeynames=True,verbose=True):
    """
    Generating the legend for a set of plotting symboles

    --- EXAMPLE OF USE ---
    lce.plot_literature_fitscatalog_legend(legendshape=(14.0, 2.0),ncol=4,extra_textlist=['UVES restults'],extra_symlist=['o'])

    """
    refdic = lce.referencedictionary()

    maindir = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/'
    plotname = maindir+'plots/literaturecollection_emissionlinestrengths_plot_legend.pdf'


    fig = plt.figure(figsize=legendshape)
    #fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.97, bottom=0.15, top=0.95)
    Fsize         = 14.0
    lthick        = 2.0
    marksize      = 6.0
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    ax  = fig.add_subplot(111)
    ax.spines['bottom'].set_color('None')
    ax.spines['top'].set_color('None')
    ax.spines['right'].set_color('None')
    ax.spines['left'].set_color('None')
    ax.xaxis.label.set_color('None')
    ax.yaxis.label.set_color('None')
    ax.tick_params(axis='x', colors='White')
    ax.tick_params(axis='y', colors='White')

    for key in refdic.keys():
        if key == 'dummy':
            continue

        label = refdic[key][1].replace('&','\&').replace('_','\_')
        if showkeynames:
            label = key+': '+label

        ax.errorbar(-5000,-5000,xerr=None,yerr=1,marker=refdic[key][2],lw=0, markersize=marksize,alpha=1.0,
                    markerfacecolor='None',ecolor='k',markeredgecolor='black',zorder=1,label=label)

    for ee, etl in enumerate(extra_textlist):
        ax.errorbar(-5000,-5000,xerr=None,yerr=1,marker=extra_symlist[ee],lw=0, markersize=marksize,alpha=1.0,
                    markerfacecolor='k',ecolor='k',markeredgecolor='black',zorder=1,label=etl)

    leg = ax.legend(fancybox=True, loc='upper center',prop={'size':Fsize/1.0},ncol=ncol,numpoints=1,
                    bbox_to_anchor=(0.48, 1.1),)  # add the legend
    leg.get_frame().set_alpha(1.0)


    if verbose: print('   Saving plot to \n   '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def colname2NEOGAL(colname):
    """
    Translate between the column names of emission line table and the NEOGAL data column names

    """
    translatedic = {}

    translatedic['f_NV']               = 'NV1240'
    translatedic['f_CIV1']             = 'CIV1549'
    translatedic['f_CIV']              = 'CIV1550'# Not a NEOGAL column but calculated in lce.add_photoionization_models_to_plot()
    translatedic['f_CIV2']             = 'CIV1551'
    translatedic['f_CIII']             = 'CIII1908'
    translatedic['f_CIII1']            = 'CIII1907'
    translatedic['f_CIII2']            = 'CIII1910'
    translatedic['f_HeII']             = 'HeII1640'
    translatedic['f_OIII1']            = 'OIII1661'
    translatedic['f_OIII']             = 'OIII1663'
    translatedic['f_OIII2']            = 'OIII1666'
    translatedic['f_SiIII1']           = 'SiIII1883'
    translatedic['f_SiIII']            = 'SiIII1888'
    translatedic['f_SiIII2']           = 'SiIII1892'# Not a NEOGAL column but calculated in lce.add_photoionization_models_to_plot()


    translatedic['FR_CIIICIV']         = 'CIII1908/CIV1550'
    translatedic['FR_CIVCIII']         = 'CIV1550/CIII1908'
    translatedic['FR_CIVHeII']         = 'CIV1550/HeII1640'

    if colname not in translatedic.keys():
        sys.exit('No available translation between NEOGAL data and colname='+colname+
                 ' in lce.colname2NEOGAL(); add one to proceed.')

    return translatedic[colname]

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def add_photoionization_models_to_plot(piplotparam,verbose=True):
    """
    Wrapper to add NEOGAL photoionization model grids to flux ratio plots

    --- INPUT ---
    piplotparam            The photoionization plot parameters.

    --- EXAMPLE OF RUN ---

    """
    import NEOGALmodels as nm
    modeldataSF  = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/nebular_emission/')
    modeldataAGN = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/AGN_NLR_nebular_feltre16/')

    x2plot, y2plot, varyparam, cutSFmodels, markersize, SFmarker, AGNmarker = piplotparam
    # varyparam options: 'Zgas','logUs','xid','nh','COratio','Mcutoff'
    logcolors = ['Zgas']

    CIIIdoubletratio = 1.5
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # defining columns to look for in model data; also determines whether to show ratios or not
    xcol = []
    if '/' in x2plot:
        xcol.append(x2plot.split('/')[0])
        xcol.append(x2plot.split('/')[1])
    else:
        xcol.append(x2plot)

    ycol = []
    if '/' in y2plot:
        ycol.append(y2plot.split('/')[0])
        ycol.append(y2plot.split('/')[1])
    else:
        ycol.append(y2plot)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if cutSFmodels:
        if verbose: print(' - Performing cut on model SF model grid')
        xid     = 'dummy'
        nh      = 'dummy'
        COratio = 1.00
        Mcutoff = 100
    else:
        if verbose: print(' - Showing all SF model grids, i.e., setting xid, nh, COratio and Mcutoff to dummy values')
        xid     = 'dummy'
        nh      = 'dummy'
        COratio = 'dummy'
        Mcutoff = 'dummy'

    # - - - - - - - - - - - - - - - - - - - - - - - -
    legenddic = {}
    legenddic['Zgas']     = r'Z$_\textrm{gas}$'
    legenddic['logUs']    = r'log U'
    legenddic['xid']      = r'$\xi_\textrm{d}$'
    legenddic['nh']       = r'n$_\textrm{H}$  [cm$^3$]'
    legenddic['COCOsol']  = r'C/O / [C/O]$_\odot$'
    legenddic['mup']      = r'M$_\textrm{cut IMF}$ / [M$_\odot]$'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if cutSFmodels:
        goodentSF  = np.where( (modeldataSF['mup'] == Mcutoff) &
                                #(modeldata['xid'] == xid) &
                                #(modeldata['nh'] == nh) &
                                (modeldataSF['COCOsol'] == COratio) )[0]
        modeldataSF  = modeldataSF[goodentSF]
        infostrSFcut = '(Mcutoff(SF)='+str(Mcutoff)+', COratio(SF)='+str(COratio)+') '#+\
                       #' Showing Zgas=all, zid=all, nh=all, logU=all '
    else:
        infostrSFcut = ''
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # collect data to plot fot star formation models
    NgoodentSF  = len(modeldataSF)

    if NgoodentSF > 1:
        if verbose: print(' - Getting data for '+str(NgoodentSF)+' data points satisfying (SFR)model selection ')
        varydatSF  = modeldataSF[varyparam]
        if varyparam in logcolors:
            varydatSF  = np.log10(modeldataSF[varyparam])

        xdat_SF = []
        for ll, linestr in enumerate(xcol):
            if linestr == 'CIV1550':
                xdat_SF.append(modeldataSF['CIV1548']+modeldataSF['CIV1551'])
            elif linestr == 'CIII1907':
                xdat_SF.append(modeldataSF['CIII1908']/(1.+1./CIIIdoubletratio))
            elif linestr == 'CIII1910':
                ciii1907 = modeldataSF['CIII1908']/(1.+1./CIIIdoubletratio)
                xdat_SF.append(ciii1907/CIIIdoubletratio)
            elif linestr == 'OIII1663':
                xdat_SF.append(modeldataSF['OIII1661']+modeldataSF['OIII1666'])
            elif linestr == 'SiIII1892':
                xdat_SF.append(modeldataSF['SiIII1888']-modeldataSF['SiIII1883'])
            else:
                xdat_SF.append(modeldataSF[linestr])
        if '/' in x2plot:
            xval_SF   = xdat_SF[0]/xdat_SF[1]
        else:
            xval_SF   = xdat_SF[0] # FLUX UNITS!!

        ydat_SF = []
        for ll, linestr in enumerate(ycol):
            if linestr == 'CIV1550':
                ydat_SF.append(modeldataSF['CIV1548']+modeldataSF['CIV1551'])
            elif linestr == 'CIII1907':
                ydat_SF.append(modeldataSF['CIII1908']/(1.+1./CIIIdoubletratio))
            elif linestr == 'CIII1910':
                ciii1907 = modeldataSF['CIII1908']/(1.+1./CIIIdoubletratio)
                ydat_SF.append(ciii1907/CIIIdoubletratio)
            elif linestr == 'OIII1663':
                ydat_SF.append(modeldataSF['OIII1661']+modeldataSF['OIII1666'])
            elif linestr == 'SiIII1892':
                ydat_SF.append(modeldataSF['SiIII1888']-modeldataSF['SiIII1883'])
            else:
                ydat_SF.append(modeldataSF[linestr])

        if '/' in y2plot:
            yval_SF   = ydat_SF[0]/ydat_SF[1]
        else:
            yval_SF   = ydat_SF[0] # FLUX UNITS!!

    else:
        print(' WARNING: lce.add_photoionization_models_to_plot >>>'
              ' Less than 2 (SFR)model grid points to plot; no data added to plot')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # collect data to plot fot AGN models
    NgoodentAGN  = len(modeldataAGN)

    if NgoodentAGN > 1:
        if verbose: print(' - Getting data for '+str(NgoodentAGN)+' data points satisfying (AGN)model selection ')
        varydatAGN  = modeldataAGN[varyparam]
        if varyparam in logcolors:
            varydatAGN = np.log10(modeldataAGN[varyparam])


        xdat_AGN = []
        for ll, linestr in enumerate(xcol):
            if linestr == 'CIV1550':
                xdat_AGN.append(modeldataAGN['CIV1548']+modeldataAGN['CIV1551'])
            elif linestr == 'OIII1663':
                xdat_AGN.append(modeldataAGN['OIII1661']+modeldataAGN['OIII1666'])
            elif linestr == 'CIII1908':
                xdat_AGN.append(modeldataAGN['CIII1907']+modeldataAGN['CIII1910'])
            elif linestr == 'SiIII1892':
                xdat_AGN.append(modeldataAGN['SiIII1888']-modeldataAGN['SiIII1883'])
            else:
                xdat_AGN.append(modeldataAGN[linestr])
        if '/' in x2plot:
            xval_AGN   = xdat_AGN[0]/xdat_AGN[1]
        else:
            xval_AGN   = xdat_AGN[0] # FLUX UNITS!!

        ydat_AGN = []
        for ll, linestr in enumerate(ycol):
            if linestr == 'CIV1550':
                ydat_AGN.append(modeldataAGN['CIV1548']+modeldataAGN['CIV1551'])
            elif linestr == 'OIII1663':
                ydat_AGN.append(modeldataAGN['OIII1661']+modeldataAGN['OIII1666'])
            elif linestr == 'CIII1908':
                ydat_AGN.append(modeldataAGN['CIII1907']+modeldataAGN['CIII1910'])
            elif linestr == 'SiIII1892':
                ydat_AGN.append(modeldataAGN['SiIII1888']-modeldataAGN['SiIII1883'])
            else:
                ydat_AGN.append(modeldataAGN[linestr])
        if '/' in y2plot:
            yval_AGN   = ydat_AGN[0]/ydat_AGN[1]
        else:
            yval_AGN   = ydat_AGN[0] # FLUX UNITS!!

    else:
        print(' WARNING: lce.add_photoionization_models_to_plot >>>'
              ' Less than 2 (AGN)model grid points to plot; no data added to plot')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    cmap    = plt.cm.get_cmap('Reds') #'plasma', 'Greys' 'copper_r'
    cmin    = np.min(np.append(varydatSF,varydatAGN))
    cmax    = np.max(np.append(varydatSF,varydatAGN))
    colnorm = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, 30) #cmax-cmin)
    mm      = plt.cm.ScalarMappable(cmap=cmap)
    mm.set_array(cmaparr)

    cmapAGN    = plt.cm.get_cmap('Blues') # 'copper_r' 'plasma', 'Reds_r' 'spring'
    cminAGN    = np.min(np.append(varydatSF,varydatAGN))
    cmaxAGN    = np.max(np.append(varydatSF,varydatAGN))
    colnormAGN = matplotlib.colors.Normalize(vmin=cminAGN,vmax=cmaxAGN)
    cmaparrAGN = np.linspace(cminAGN, cmaxAGN, 30) #cmax-cmin)
    mmAGN      = plt.cm.ScalarMappable(cmap=cmapAGN)
    mmAGN.set_array(cmaparrAGN)

    if varyparam == 'Zgas':
        # colortickvals   = [1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 0.003048, 4e-3, 6e-3, 8e-3, 0.01, 0.01524, 0.02, 0.03, 0.04, 0.07] # 0.014, 0.017,
        Zsolar          = 0.01524
        colortickvals   = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0])*Zsolar
        #colortickvals   = [1e-4, 3e-4, 7e-4, 1e-4, 0.003048, 0.007, 0.01524, 0.03, 0.07]
        # colorlabels     = [ str(ct) for ct in colortickvals]
        # colorlabels[4]  =  '0.2Z$_\odot$' # = 0.003048
        # colorlabels[6] =  'Z$_\odot$'     # = 0.01524
        colorlabels     = ['0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1.0', '2.0', '4.0']
        colortickvals   = np.log10(np.asarray(colortickvals))
        cbarlegend      = r'Z$_\textrm{gas}$/Z$_\odot$' # legenddic[varyparam]
    elif varyparam == 'logUs':
        colortickvals = [-5.0,-4.5,-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0]
        colorlabels   = [ str(ct) for ct in colortickvals]
        cbarlegend    = legenddic[varyparam]
    else:
        ncolticks     = 10.
        colortickvals = np.arange(cmin,cmax,np.abs(cmax-cmin)/ncolticks)
        colorlabels   = [ str(ct) for ct in colortickvals]

        if varyparam in logcolors:
            cbarlegend = 'log10('+legenddic[varyparam]+')'
        else:
            cbarlegend = legenddic[varyparam]

    colshrink   = 1.0
    colaspect   = 30
    colbarscale = 2.1
    colanchor   = (-0.9,1.0)
    colshrink   = colshrink/colbarscale
    colaspect   = colaspect/colbarscale
    cextend     = 'neither'

    cbSF      = plt.colorbar(mm,extend=cextend,orientation='vertical',ticks=colortickvals,
                             pad=0.01,aspect=colaspect,shrink=colshrink,anchor=colanchor,use_gridspec=False)
    cbSF.ax.set_yticklabels(len(colorlabels)*[''],rotation=0)
    cbSF.set_label('')     #Zgas is unitless as it is a mass ratio (see Gutkin et al. 2016 Eq. 10)

    cbAGN      = plt.colorbar(mmAGN,extend=cextend,orientation='vertical',ticks=colortickvals,
                              pad=0.01,aspect=colaspect,shrink=colshrink,anchor=(0.7,1.0),use_gridspec=False)
    cbAGN.ax.set_yticklabels(colorlabels,rotation=0)
    cbAGN.set_label(cbarlegend)

    for vdSF in np.unique(varydatSF):
        SFcol    = cmap(colnorm(vdSF))
        SFcolent = np.where(varydatSF == vdSF)

        plt.scatter(xval_SF[SFcolent],yval_SF[SFcolent],s=markersize,
                    marker=SFmarker,lw=0.2, facecolor='None',edgecolor=SFcol, zorder=5)

    for vdAGN in np.unique(varydatAGN):
        AGNcol    = cmapAGN(colnormAGN(vdAGN)) # 'gray'
        AGNcolent = np.where(varydatAGN == vdAGN)

        plt.scatter(xval_AGN[AGNcolent],yval_AGN[AGNcolent],s=markersize,
                    marker=AGNmarker,lw=0.2, facecolor='None',edgecolor=AGNcol, zorder=5)

    titleaddition = infostrSFcut
    return titleaddition

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_TEMPLATE(fluxscale=1.0,verbose=True):
    """
    Data collected from REFERENCE+YEAR

    NOTE ON DATA

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'TEMP'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']      = np.array(['Q2343-BX418'])
    datadic['id']        = np.array([9999]) + baseid
    rasex                = np.array(['04:22:00.81'])
    decsex               = np.array(['-38:37:03.59'])
    datadic['ra']        = acoord.Angle(rasex, u.hour).degree
    datadic['dec']       = acoord.Angle(decsex, u.degree).degree
    datadic['redshift']  = np.array([5.55])
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    datadic['f_Lya']          = np.array([])
    datadic['ferr_Lya']       = np.array([])
    datadic['sigma_Lya']      = np.array(['FWHM']) / 2.355 # in Angstrom otherwise: * 1215.6701 / (astropy.constants.c.value/1000.)
    datadic['sigmaerr_Lya']   = np.array(['FWHM']) / 2.355 # in Angstrom otherwise: * 1215.6701 / (astropy.constants.c.value/1000.)
    datadic['EW0_Lya']        = np.array([])
    datadic['EW0err_Lya']     = np.array([])

    datadic['f_CIII1']        = np.array([])
    datadic['ferr_CIII1']     = np.array([])
    datadic['sigma_CIII1']    = np.array(['FWHM']) / 2.355 # in Angstrom
    datadic['sigmaerr_CIII1'] = np.array(['FWHM']) / 2.355 # in Angstrom
    datadic['EW0_CIII1']      = np.array([])
    datadic['EW0err_CIII1']   = np.array([])

    datadic['f_CIII2']        = np.array([])
    datadic['ferr_CIII2']     = np.array([])
    datadic['sigma_CIII2']    = np.array(['FWHM']) / 2.355 # in Angstrom
    datadic['sigmaerr_CIII2'] = np.array(['FWHM']) / 2.355 # in Angstrom
    datadic['EW0_CIII2']      = np.array([])
    datadic['EW0err_CIII2']   = np.array([])

    linename = 'CIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    # ---------------------------------------------------------------------------------
    if verbose: print('   Converting fluxes to 1e-20 erg/s/cm2 using fluxscale = '+str(fluxscale))
    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    if verbose: print('   Making sure EWs are rest-frame EWs, i.e., EW0')
    for key in datadic.keys():
        if key.startswith('EW0'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99] / \
                                                       (1 + datadic['redshift'][np.abs(datadic[key]) != 99])

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_sen17(fluxscale=1e5,verbose=True):
    """
    Data collected from Senchyna+2017

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'sen17'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['id']        = np.array([2, 36, 80, 82, 110, 111, 179, 182, 191, 198]) + baseid
    datadic['name']      = (datadic['id']-baseid).astype(int).astype(str)

    rasex                = np.array(['9:44:01.87','10:24:29.25 ','9:42:56.74 ','11:55:28.34','9:42:52.78',
                                     '12:30:48.60','11:29:14.15 ','11:48:27.34 ','12:15:18.60 ','12:22:25.79'])
    decsex               = np.array(['-0:38:32.2','5:24:51.0','9:28:16.2','57:39:52.0','35:47:26.0',
                                     '12:02:42.8','20:34:52.0','25:46:11.8','20:38:26.7','4:34:04.8'])
    datadic['ra']        = acoord.Angle(rasex, u.hour).degree
    datadic['dec']       = acoord.Angle(decsex, u.degree).degree

    distances            = np.array([19. ,141. ,46. ,76. ,63. ,16. ,25. ,191. ,10. ,16.])
    datadic['redshift']  = np.array([acoord.Distance(objdist, u.Mpc).compute_z(cosmology=acosmo.Planck15) for objdist in distances])

    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    datadic['f_CIV1']       = np.array([0.20,3.3,1.6,1.92,0.20,0.79,2.0,1.54,1.0,3.1])
    datadic['ferr_CIV1']    = np.array([0.02,99,99,0.05,0.05,0.04,99,0.04,99,99])
    datadic['EW0_CIV1']     = np.array([0.25,0.4,0.5,0.67,0.08,0.41,0.5,0.96,0.5,0.7])
    datadic['EW0err_CIV1']  = np.array([0.02,99,99,0.02,0.02,0.02,99,0.03,99,99])

    datadic['f_CIV2']       = np.array([0.19,3.3,1.6,1.32,1.14,0.49,2.0,0.94,1.0,3.1])
    datadic['ferr_CIV2']    = np.array([0.01,99,99,0.05,0.09,0.03,99,0.05,99,99])
    datadic['EW0_CIV2']     = np.array([0.22,0.4,0.4,0.41,0.37,0.24,0.4,0.55,0.4,0.6])
    datadic['EW0err_CIV2']  = np.array([0.02,99,99,0.02,0.03,0.01,99,0.03,99,99])

    linename = 'CIV'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    datadic['f_HeII']    = np.array([1.45, 1.1, 0.8, 1.22, 1.2,2.37, 0.8, 1.42, 1.0, 0.8])
    datadic['ferr_HeII'] = np.array([0.06, +99.0, +99.0, 0.05, +99.0, 0.05, +99.0, 0.05, +99.0, +99.0])
    datadic['EW0_HeII']    = np.array([1.7, 0.4, 0.7, 0.44, 0.4, 1.42, 0.6, 0.94, 0.5, 0.5])
    datadic['EW0err_HeII'] = np.array([0.06, +99.0, +99.0, 0.02, +99.0, 0.03, +99.0, 0.04, +99.0, +99.0])

    datadic['f_OIII1']       = np.array([1.40,0.70,0.60,2.33,0.9,0.56,1.1,0.99,1.19,1.1])
    datadic['ferr_OIII1']    = np.array([0.04,0.06,0.04,0.05,99,0.06,99,0.04,0.07,00])
    datadic['EW0_OIII1']     = np.array([1.94,0.26,0.55,0.85,0.3,0.35,0.9,0.7,0.64,0.8])
    datadic['EW0err_OIII1']  = np.array([0.06,0.02,0.04,0.02,99,0.04,99,0.03,0.04,99])

    datadic['f_OIII2']       = np.array([3.35,2.16,1.69,5.05,1.28,2.68,1.36,2.27,np.nan,1.1])
    datadic['ferr_OIII2']    = np.array([0.04,0.08,0.07,0.05,0.05,0.07,0.06,0.05,np.nan,99])
    datadic['EW0_OIII2']     = np.array([5.05,0.84,1.68,1.89,0.48,1.73,1.17,1.73,np.nan,0.9])
    datadic['EW0err_OIII2']  = np.array([0.1,0.03,0.08,0.02,0.02,0.05,0.06,0.04,np.nan,99])

    linename = 'OIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    datadic['f_CIII']       = np.array([11.06,7.60,4.1,22.08,5.2,3.9,6.8,12.30,15.91,3.79])
    datadic['ferr_CIII']    = np.array([0.63,0.88,99,0.40,99,99,0.23,0.33,0.39,0.30])
    datadic['EW0_CIII']     = np.array([14.68,4.98,4.0,12.09,2.8,3.3,8.71,13.35,11.33,3.38])
    datadic['EW0err_CIII']  = np.array([1.07,0.59,99,0.30,99,99,0.42,0.52,0.34,0.31])

    # ---------------------------------------------------------------------------------
    if verbose: print('   Converting fluxes to 1e-20 erg/s/cm2 using fluxscale = '+str(fluxscale))
    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_nan19(fluxscale=1.0,verbose=True):
    """
    Data collected from Nanaykkara et al. (2019)

    EW values are obtained via privat communication as the published table 3 is faulty.

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'nan19'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['id']        = np.array([1024, 1036, 1045, 1079, 1273, 3621, 87,   109, 144,    97,   39,  84, 161]) + baseid
    datadic['name']      = datadic['name']      = (datadic['id']-baseid).astype(int).astype(str)
    rasex                = np.array(['03:32:31.45','03:32:43.39','03:32:33.78','03:32:37.88','03:32:35.48',
                                     '03:32:39.52','22:32:54.86','22:32:56.16','22:32:58.93','10:00:34.01',
                                     '04:22:00.81','04:22:01.45','04:22:01.52'])
    decsex               = np.array(['-27:47:25.12','-27:47:10.54','-27:48:14.35','-27:47:56.75','-27:46:16.91',
                                     '-27:48:53.53','-60:33:42.12','-60:34:11.74','-60:34:00.07','+02:03:57.99',
                                     '-38:37:03.59','-38:37:20.79','-38:37:19.68'])
    datadic['ra']        = acoord.Angle(rasex, u.hour).degree
    datadic['dec']       = acoord.Angle(decsex, u.degree).degree
    datadic['redshift']  = np.array([2.87, 2.69, 2.61, 2.68, 2.17, 3.07, 2.67, 2.2, 4.02, 2.11, 3.96, 3.1, 3.1])
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    datadic['f_HeII']       = np.array([177., 142., 156., 290., 217., 213., 59., 54., 48., 306., 153., 161., 318.])
    datadic['ferr_HeII']    = np.array([52., 53., 58., 91., 79., 45., 12., 13., 12., 55., 33., 29., 39.])
    datadic['sigma_HeII']   = np.array([5.,4.,4.,11.,5.,6.,4.,3.,4.,5.,7.,4.,9.]) / 2.355
    # datadic['EW0_HeII']     = np.array([18.9,12.9,12.6,35.6,13.7,6.8,11.4,11.0,5.2,5.5,3.6,8.5,28.3])
    # datadic['EW0err_HeII']  = np.array([3.5,2.9,2.2,11.5,3.7,-99,2.1,1.0,1.5,2.3,-99,2.9,-99])
    datadic['EW0_HeII']     = np.array([3.1 ,3.0 ,3.6 ,9.1 ,8.2 ,10.3,4.8 ,3.2 ,10.7,17.9,24.7,26.9,50.3])
    datadic['EW0err_HeII']  = np.array([0.6  ,1.1  ,0.9  ,2.6  ,2.0  ,-99.0,0.8  ,0.4  ,2.2  ,2.4  ,-99.0,3.5  ,-99.0])

    datadic['f_CIII1']       = np.array([308.,436.,388.,111.,402.,252.,78.,71.,np.nan,369.,np.nan,174.,143.])
    datadic['ferr_CIII1']    = np.array([49.,36.,49.,99,48.,99,11.,12.,np.nan,50.,np.nan,47.,99])
    datadic['sigma_CIII1']   = np.array([6.,4.,4.,4.,3.,4.,3.,3.,np.nan,3.,np.nan,3.,3.]) / 2.355
    # datadic['EW0_CIII1']     = np.array([18.5,4.4,6.6,16.6,6.5,9.3,5.5,8.7,np.nan,14.0,np.nan,8.9,11.7])
    # datadic['EW0err_CIII1']  = np.array([3.1,1.0,1.2,0.7,2.9,-99,0.7,0.7,np.nan,2.7,np.nan,-99,-99])
    datadic['EW0_CIII1']     = np.array([ 6.0, 12.4, 11.5, 0.7, 21.3, 7.2, 7.1, 5.7, np.nan, 28.9, np.nan, 18.6, 10.4])
    datadic['EW0err_CIII1']  = np.array([0.9,0.8,1.2,0.6,1.8,-99.0 ,0.8,0.5,np.nan,3.1,np.nan,-99.0 ,-99.0])

    datadic['f_CIII2']       = np.array([206.,300.,211.,111.,271.,212.,33.,63.,np.nan,258.,np.nan,72.,62.])
    datadic['ferr_CIII2']    = np.array([42.,41.,54.,99,47.,99,11.,12.,np.nan,62.,np.nan,23.,99])
    datadic['sigma_CIII2']   = np.array([np.nan]*len(datadic['id']))
    # datadic['EW0_CIII2']     = np.array([20.5,8.3,11.8,15.8,0.4,15.3,9.6,6.4,np.nan,5.3,np.nan,3.4,1.5])
    # datadic['EW0err_CIII2']  = np.array([3.5,1.4,1.7,0.7,1.0,-99,0.7,0.7,np.nan,4.4,np.nan,-99,-99])
    datadic['EW0_CIII2']     = np.array([4.0,8.5,6.3,1.5,14.4,12.4,3.0,5.0,np.nan,20.2,np.nan,12.0,0.3])
    datadic['EW0err_CIII2']  = np.array([0.5,1.0,1.1,0.7,1.5,-99.0 ,0.8,0.5,np.nan,3.8,np.nan,-99.0 ,-99.0])

    linename = 'CIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    datadic['f_OIII1']       = np.array([151.,155.,186.,162.,165.,106.,33.,38.,37.,148.,71.,66.,64.])
    datadic['ferr_OIII1']    = np.array([99,99,99,99,99,99,99,99,11.,43.,99,99,99])
    datadic['sigma_OIII1']   = np.array([np.nan]*len(datadic['id']))
    # datadic['EW0_OIII1']     = np.array([23.2,14.7,14.9,17.1,11.2,11.5,10.5,12.9,2.9,6.1,4.7,4.1,5.7])
    # datadic['EW0err_OIII1']  = np.array([3.6,1.1,1.3,0.6,1.6,-99,0.7,0.8,1.9,1.5,-99,2.2,-99])
    datadic['EW0_OIII1']     = np.array([1.3 ,2.1 ,3.2 ,0.2 ,3.6 ,10.2,2.1 ,1.5 ,9.7 ,8.8 ,5.1 ,8.5 ,1.7])
    datadic['EW0err_OIII1']  = np.array([0.6  ,0.7  ,0.8  ,0.7  ,1.1  ,-99.0,0.5  ,0.5  ,1.8  ,1.4  ,-99.0,2.4  ,-99.0])

    datadic['f_OIII2']       = np.array([161.,230.,200.,81.,195.,105.,49.,72.,160.,284.,66.,54.,58.])
    datadic['ferr_OIII2']    = np.array([48.,54.,56.,99,56.,99,11.,13.,21.,44.,99,99,19.])
    datadic['sigma_OIII2']   = np.array([np.nan]*len(datadic['id']))
    # datadic['EW0_OIII2']     = np.array([21.7,12.0,13.5,17.5,7.6,11.4,8.8,10.1,24.6,2.0,9.6,6.3,2.7])
    # datadic['EW0err_OIII2']  = np.array([3.3,1.1,1.4,0.7,1.1,-99,0.6,0.6,3.0,2.0,-99,1.6,-99])
    datadic['EW0_OIII2']     = np.array([2.8 ,4.8 ,4.6 ,0.2 ,7.2 ,10.1,3.8 ,4.3 ,37.2,16.9,1.3 ,6.3 ,6.9])
    datadic['EW0err_OIII2']  = np.array([0.6  ,0.7  ,0.8  ,0.4  ,1.4  ,-99.0,0.5  ,0.5  ,4.3  ,1.9  ,-99.0,2.0  ,-99.0])

    linename = 'OIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    datadic['f_SiIII1']       = np.array([238.,167.,141.,64.,141.,122.,32.,57.,np.nan,131.,np.nan,393.,181.])
    datadic['ferr_SiIII1']    = np.array([99,55.,38.,99,55.,99,10.,12.,np.nan,99,np.nan,59.,99])
    datadic['sigma_SiIII1']   = np.array([np.nan]*len(datadic['id']))
    # datadic['EW0_SiIII1']     = np.array([20.9,12.2,14.0,17.3,7.5,16.6,9.7,9.9,np.nan,7.2,np.nan,32.0,6.8])
    # datadic['EW0err_SiIII1']  = np.array([3.2,0.9,1.3,0.6,1.6,-99,0.6,1.0,np.nan,2.2,np.nan,-99,-99])
    datadic['EW0_SiIII1']     = np.array([3.6 ,4.6 ,4.1 ,0.0 ,7.3 ,14.5,2.9 ,4.5 ,np.nan,7.7 ,np.nan,39.4,6.6 ])
    datadic['EW0err_SiIII1']  = np.array([0.7,0.6,0.7,0.4,1.9,-99.0 ,0.5,0.7,np.nan,1.9,np.nan,-99.0 ,-99.0])

    datadic['f_SiIII2']       = np.array([317.,184.,129.,118.,153.,122.,47.,36.,np.nan,134.,np.nan,156.,130.])
    datadic['ferr_SiIII2']    = np.array([99,99,44.,99,99,99,13.,99,np.nan,99,np.nan,10.,99])
    datadic['sigma_SiIII2']   = np.array([np.nan]*len(datadic['id']))
    # datadic['EW0_SiIII2']     = np.array([21.3,14.3,14.3,16.9,11.8,10.2,8.5,12.2,np.nan,11.9,np.nan,14.2,12.6])
    # datadic['EW0err_SiIII2']  = np.array([3.0,1.2,1.7,0.7,1.7,-99,0.9,0.9,np.nan,2.0,np.nan,-99,-99])
    datadic['EW0_SiIII2']     = np.array([3.2,2.5,3.8,0.4,3.0,7.0 ,4.1,2.2,np.nan,3.0,np.nan,23.3,12.0])
    datadic['EW0err_SiIII2']  = np.array([1.1,0.8,0.8,0.6,1.4,-99.0,0.7,0.7,np.nan,1.8,np.nan,-99.0,-99.0])

    linename = 'SiIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    # ---------------------------------------------------------------------------------
    if verbose: print('   Converting fluxes to 1e-20 erg/s/cm2 using fluxscale = '+str(fluxscale))
    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    if verbose: print('   Making sure EWs are rest-frame EWs, i.e., EW0')
    for key in datadic.keys():
        if key.startswith('EW0'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99] / \
                                                       (1 + datadic['redshift'][np.abs(datadic[key]) != 99])

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)

    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_sch17(fluxscale=1e3,verbose=True):
    """
    Data collected from Schmidt+2017 and Mainali+2017 on lensed quintet behind RXJ2248 (AS1063)

    Combined data from Schmidt+17 table 3 and Mainali+17 OIII detection are provided

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'sch17'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    # datadic['id']        = np.array([224801, 224802, 224803, 224804, 224805]) + baseid
    # datadic['ra']        = np.array([342.18408,342.1890479,342.171296,342.18104,342.190889,342.18407])
    # datadic['dec']       = np.array([-44.5316378,-44.5300249,-44.519812,-44.53461,-44.537461,-44.53532])
    # datadic['redshift']  = np.array([6.11,6.11,6.11,6.11,6.11])
    datadic['name']      = np.array(['RXJ2248_6.11_quintet'])
    datadic['id']        = np.array([22480005]) + baseid
    datadic['ra']        = np.array([342.18407])
    datadic['dec']       = np.array([-44.53532])
    datadic['redshift']  = np.array([6.11])

    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    datadic['f_Lya']       = np.array([0.89])# Schmidt+17
    datadic['ferr_Lya']    = np.array([0.24])# Schmidt+17
    datadic['EW0_Lya']     = np.array([68.0])# Schmidt+17
    datadic['EW0err_Lya']  = np.array([6.0])# Schmidt+17

    datadic['f_NV']       = np.array([0.27]) # Mainali+17
    datadic['ferr_NV']    = np.array([99])
    datadic['EW0_NV']     = np.array([3.45]) # Mainali+17
    datadic['EW0err_NV']  = np.array([99])

    datadic['f_CIV']       = np.array([0.21])# Schmidt+17
    datadic['ferr_CIV']    = np.array([0.08])# Schmidt+17
    datadic['EW0_CIV']     = np.array([24.0])# Schmidt+17
    datadic['EW0err_CIV']  = np.array([4.0])# Schmidt+17

    datadic['f_HeII']       = np.array([0.63])# Schmidt+17
    datadic['ferr_HeII']    = np.array([99])# Schmidt+17
    datadic['EW0_HeII']     = np.array([78.0])# Schmidt+17
    datadic['EW0err_HeII']  = np.array([99])# Schmidt+17

    # datadic['f_OIII']       = np.array([0.75])# Schmidt+17
    # datadic['ferr_OIII']    = np.array([99])# Schmidt+17
    # datadic['EW0_OIII']     = np.array([81.0])# Schmidt+17
    # datadic['EW0err_OIII']  = np.array([99])# Schmidt+17

    datadic['f_OIII1']       = np.array([0.17])# Mainali+17
    datadic['ferr_OIII1']    = np.array([0.06])# Mainali+17
    datadic['EW0_OIII1']     = np.array([2.9])# Mainali+17
    datadic['EW0err_OIII1']  = np.array([1.4])# Mainali+17

    datadic['f_OIII2']       = np.array([0.27])# Mainali+17
    datadic['ferr_OIII2']    = np.array([0.06])# Mainali+17
    datadic['EW0_OIII2']     = np.array([4.6])# Mainali+17
    datadic['EW0err_OIII2']  = np.array([1.6])# Mainali+17

    linename = 'OIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    datadic['f_CIII']       = np.array([0.45])# Schmidt+17
    datadic['ferr_CIII']    = np.array([99])# Schmidt+17
    datadic['EW0_CIII']     = np.array([60.0])# Schmidt+17
    datadic['EW0err_CIII']  = np.array([99])# Schmidt+17

    # ---------------------------------------------------------------------------------
    if verbose: print('   Converting fluxes to 1e-20 erg/s/cm2 using fluxscale = '+str(fluxscale))
    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)

    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_rig14(fluxscale=1.0,verbose=True):
    """
    Data collected from Rigby+2014

    Only EW provided

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'rig14'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']      = np.array(['S0004','S0108','R0327_Knot_E','R0327_Knot_B','R0327_Knot_G','R0327_Knot_U','S0957','S1441'])
    datadic['id']        = np.array([4,108,3270001,3270002,3270003,3270004,957,1441]) + baseid
    rasex                = np.array(['00:04:51.7','01:08:42.2','03:27:27','03:27:27','03:27:27',
                                     '03:27:27','09:57:38.7','14:41:33.2'])
    decsex               = np.array(['-01:03:21','+06:24:44','-13:26:09','-13:26:09','-13:26:09',
                                     '-13:26:09','+05:09:29','-00:54:01'])
    datadic['ra']        = acoord.Angle(rasex, u.hour).degree
    datadic['dec']       = acoord.Angle(decsex, u.degree).degree
    datadic['redshift']  = np.array([1.6811,1.91921,1.7034,1.7034,1.7034,1.7034,1.82042,1.666])
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    datadic['EW0_Lya']     = np.array([3.3,8.9,1.8,3.75,5.7,1.65,8.1,5.1])
    datadic['EW0err_Lya']  = np.array([0.4,0.6,99,99,99,99,1.2,99])

    datadic['EW0_CIII1']     = np.array([0.2,1.2,0.8,1.6,1.65,1.5,2.0,0.93])
    datadic['EW0err_CIII1']  = np.array([0.1,0.1,0.1,0.4,99,0.1,0.3,99])

    datadic['EW0_CIII2']     = np.array([0.25,0.7,1.2,1.4,1.575,0.9,1.2,0.93])
    datadic['EW0err_CIII2']  = np.array([0.1,0.1,0.1,0.35,99,0.1,0.3,99])

    datadic['EW0_CIII']      = datadic['EW0_CIII1'] + datadic['EW0_CIII2']
    datadic['EW0err_CIII']   = np.sqrt(datadic['EW0err_CIII1']**2 + datadic['EW0err_CIII2']**2)
    datadic['EW0err_CIII'][4]= 99
    datadic['EW0err_CIII'][7]= 99

    datadic['EW0_MgII1']     = np.array([0.45,0.9,1.3,1.1,1.6,2.8,2.8,0.95])
    datadic['EW0err_MgII1']  = np.array([0.1,0.6,0.16,0.5,0.6,0.3,0.9,0.4])

    datadic['EW0_MgII2']     = np.array([0.81,0.87,1.00,0.8,1.4,2.35,1.81,1.2])
    datadic['EW0err_MgII2']  = np.array([0.17,0.25,0.19,0.3,0.6,0.2,0.77,99])

    datadic['EW0_MgII']      = datadic['EW0_MgII1'] + datadic['EW0_MgII2']
    datadic['EW0err_MgII']   = np.sqrt(datadic['EW0err_MgII1']**2 + datadic['EW0err_MgII2']**2)
    datadic['EW0err_MgII'][7]= 99

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_rig15(fluxscale=1.0,verbose=True):
    """
    Data collected from Rigby+2015

    Only EW provided. Coordinates "converted" from names (kept in names list) in Rigby table; not precise.
    A few updated via NED and literature.

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'rig15'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']     = np.array(['RCSGA_032727-13260_Knot_E','RCSGA_032727-13260_Knot_U','RCSGA_032727-13260_Knot_B','RCSGA_032727-13260_Knot_G','SGAS_J000451.7-010321','SGAS_J010842.2+062444','SGAS_J095738.7+050929','SGAS_J090003.3+223408','SGAS_J105039.6+001730','SGAS_J122651.3+215220','SGAS_J142954.9+120239','SGAS_J152745.1+065219','SGAS_J211118.9-011431','Cosmic_Eye','Haro_15','IC_0214','Mrk_26','Mrk_347','Mrk_496','Mrk_499','Mrk_66','Pox_120','Pox_124','Tol_1924-416','Tol_41','NGC_1741','NGC_1741','Mrk_960','SBS_0218+003','Mrk_1087','Mrk_5','Mrk_1199','IRAS_08208+2816','IRAS_08339+6517','SBS_0926+606A','Arp_252','SBS_0948+532','Tol_9','SBS_1054+365','Pox_4','SBS_1319+579','SBS_1415+437','Tol_1457-262','III_Zw_107','IZw18','IZw18-NW_HIIR','IZw18-SE_HIIR','NGC_1569','Mrk_71','NGC_2403-vs38','NGC_2403-vs44','NGC_2403-vs9','NGC_3690','NGC_4214','NGC_4569','NGC_4670','NGC_4861','NGC_5055','NGC_5253-HIIR-1','NGC_5253-HIIR-2','NGC_5253-UV1','NGC_5253-UV2','NGC_5253-UV3','NGC_5457-NGC_5455','NGC_5457-NGC_5471','NGC_5457-Searle5','NGC_7552','SBS_1415+437','Tol_1214-277','Tol_1345-420','UM_469'])[4:]
    datadic['id']        = np.array([3270001,3270004,3270002,3270003,452,10842,95738,90003,105040,122651,142955,152745,211119,999,15,214,26,347,496,499,66,120,124,1924,41,17410001,17410002,960,218,1087,5,1199,8208,8339,926,252,948,9,1054,4,1319,1415,1457,107,180001,180002,180003,1569,71,240338,240344,240309,3690,4214,4569,4670,4861,5055,52530001,52530002,52530003,52530004,52530005,54575455,54575471,54575,7552,14150437,1214,1345,469])[4:] + baseid
    rasex                = ['03:27:27','03:27:27','03:27:27','03:27:27','00:04:51.7','01:08:42.2','09:57:38.7','09:00:03.3','10:50:39.6','12:26:51.3','14:29:54.9','15:27:45.1','21:11:18.9',np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,'19:24:00',np.nan,np.nan,np.nan,np.nan,'02:18:00',np.nan,np.nan,np.nan,'08:20:08','08:38:23','09:30:06.5',np.nan,'09:48:00',np.nan,'10:54:00',np.nan,'13:21:10.0','14:17:01.7','14:57:00',np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,'14:17:01.7','12:17:17.093','13:45:00',np.nan]
    decsex               = ['-13:26:0','-13:26:0','-13:26:0','-13:26:0','-01:03:21','+06:24:44','+05:09:29','+22:34:08','+00:17:30','+21:52:20','+12:02:39','+06:52:19','-01:14:31',np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,'-41:60:00',np.nan,np.nan,np.nan,np.nan,'+00:30:00',np.nan,np.nan,np.nan,'+28:16:00','+65:07:15','+60:26:52',np.nan,'+53:20:00',np.nan,'+36:50:00',np.nan,'+57:39:41','+43:30:13','-26:20:00',np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,'+43:30:13','-28:02:32.67','-42:00:00',np.nan]
    datadic['ra']        = acoord.Angle(rasex, u.hour).degree[4:]
    datadic['dec']       = acoord.Angle(decsex, u.degree).degree[4:]
    datadic['redshift']  = np.array([1.703745,1.703884,1.70360 ,1.70385 ,1.6811  ,1.91021 ,1.82042 ,2.0315  ,2.38115 ,2.9227  ,2.8219  ,2.7604  ,2.8577  ,3.07331 ,0.021371,0.030224,0.030428,0.019417,0.029304,0.025978,0.020972,0.020748,0.024217,0.009453,0.023083,0.013473,0.013473,0.021371,0.058420,0.027809,0.002642,0.013540,0.046776,0.019113,0.013642,0.032989,0.046240,0.010641,0.002010,0.011970,0.006870,0.002031,0.016812,0.019127,0.002505,0.002505,0.002505,0.000347,0.000233,0.000445,0.000445,0.000445,0.010411,0.000970,0.000784,0.003566,0.002785,0.001614,0.001358,0.001358,0.001358,0.001358,0.001358,0.000804,0.000804,0.000804,0.005365,0.002031,0.026001,0.008000,0.058146])[4:]
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    datadic['EW0_Lya']     = np.array([-1.2,-1.1,-2.5,-3.8,-3.3,-8.9,-8.1,-7.0,-19.5,-8.3,-30.0,-13.0,-1.5,np.nan,-0.69,5.97,5.46,6.19,-3.95,8.19,-4.70,-70.2,2.64,-22.3,-38.2,-0.45,1.94,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,-4.45,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,-82.6,np.nan,np.nan])[4:]
    datadic['EW0err_Lya']  = np.array([-99.,-99.,-99.,-99.,0.4,0.6,1.2,0.5,2.0,0.5,3.0,1.0,0.15,np.nan,0.50,5.70,4.83,2.88,1.13,2.88,0.79,4.11,29.2,5.10,10.38,0.13,0.53,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,0.5,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,99.,np.nan,np.nan])[4:]

    datadic['EW0_CIII']     = np.array([-2.0,-2.4,-3.0,-2.2,-0.45,-1.9,-3.2,-0.45,-2.2,-0.2,-0.1,-0.7,-2.0,-0.3,-0.11,-1.9,-3.1,0.04,-0.54,-0.31,-1.1,-14.,-3.2,-4.0,-1.7,-3.5,-3.5,-1.9,-2.6,-1.3,-2.5,-0.3,-0.4,-0.7,-3.1,-5.2,-1.2,-1.0,-0.5,-1.5,-8.3,-0.2,-1.8,-8.0,-4.2,-1.3,-4.4,-0.46,-18.7,-0.87,-1.40,-0.875,0.13,-0.23,-0.11,-0.68,-8.1,0.71,-7.5,-8.5,-2.3,0.10,0.17,-1.4,-6.4,-0.033,-0.069,-2.9,-27.,-6.2,-1.5])[4:]
    datadic['EW0err_CIII']  = np.array([0.14,0.14,0.53,-99.,0.14,0.14,0.42,0.1,0.3,0.2,0.2,0.1,-99.,-99.,0.3,1.5,1.8,0.3,0.4,0.3,0.3,1.8,6.0,0.7,2.85,0.2,0.2,0.1,0.1,0.4,0.1,0.1,0.1,0.1,0.2,0.3,0.1,0.1,0.1,0.3,0.5,0.1,0.1,0.5,0.5,0.1,0.5,0.06,1.6,0.08,0.07,0.06,0.03,0.02,0.015,0.07,1.0,0.14,0.8,0.9,0.2,0.013,0.03,0.1,1.2,0.011,0.009,0.3,13.,0.9,0.1])[4:]

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_erb10(fluxscale=1e3,verbose=True):
    """
    Data collected from Erb et al. (2010) (see also Erb+2006)

    NOTE ON DATA

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'erb10'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']      = np.array(['Q2343-BX418'])
    datadic['id']        = np.array([23430418]) + baseid
    rasex                = np.array(['23:46:18.57'])
    decsex               = np.array(['12:47:47.38'])
    datadic['ra']        = acoord.Angle(rasex, u.hour).degree
    datadic['dec']       = acoord.Angle(decsex, u.degree).degree
    datadic['redshift']  = np.array([2.3052])
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------
    # Erb+2010 table 3
    datadic['EW0_Lya']       = np.array([54.0])
    datadic['EW0err_Lya']    = np.array([1.2])
    datadic['f_Lya']         = np.array([29.3])
    datadic['ferr_Lya']      = np.array([0.4])
    datadic['sigma_Lya']     = np.array([840.0]) / 2.355   * 1215.6701 / (astropy.constants.c.value/1000.)
    datadic['sigmaerr_Lya']  = np.array([17.0]) / 2.355    * 1215.6701 / (astropy.constants.c.value/1000.)
    datadic['vshift_Lya']    = np.array([307.0])
    datadic['vshifterr_Lya'] = np.array([3.0])

    datadic['EW0_CIII']      = np.array([7.1])
    datadic['EW0err_CIII']   = np.array([0.4])
    datadic['f_CIII']        = np.array([1.4])
    datadic['ferr_CIII']     = np.array([0.1])
    datadic['sigma_CIII']    = np.array([225.0]) / 2.355   * 1907.709  / (astropy.constants.c.value/1000.)
    datadic['sigmaerr_CIII'] = np.array([36.0]) / 2.355    * 1907.709  / (astropy.constants.c.value/1000.)
    datadic['vshift_CIII']   = np.array([81.0])
    datadic['vshifterr_CIII']= np.array([21.0])

    datadic['EW0_HeII']      = np.array([2.7])
    datadic['EW0err_HeII']   = np.array([0.2])
    datadic['f_HeII']        = np.array([0.8])
    datadic['ferr_HeII']     = np.array([0.1])
    datadic['sigma_HeII']    = np.array([612]) / 2.355     * 1640.42  / (astropy.constants.c.value/1000.)
    datadic['sigmaerr_HeII'] = np.array([64]) / 2.355      * 1640.42  / (astropy.constants.c.value/1000.)
    datadic['vshift_HeII']   = np.array([-17.0])
    datadic['vshifterr_HeII']= np.array([55.0])

    datadic['EW0_OIII1']      = np.array([1.0])
    datadic['EW0err_OIII1']   = np.array([0.2])
    datadic['f_OIII1']        = np.array([0.3])
    datadic['ferr_OIII1']     = np.array([0.07])
    datadic['sigma_OIII1']    = np.array([276]) / 2.355    * 1660.809  / (astropy.constants.c.value/1000.)
    datadic['sigmaerr_OIII1'] = np.array([99])
    datadic['vshift_OIII1']   = np.array([193.0])
    datadic['vshifterr_OIII1']= np.array([57.0])

    datadic['EW0_OIII2']      = np.array([1.3])
    datadic['EW0err_OIII2']   = np.array([0.2])
    datadic['f_OIII2']        = np.array([0.4])
    datadic['ferr_OIII2']     = np.array([0.08])
    datadic['sigma_OIII2']    = np.array([235.0]) / 2.355   * 1666.150  / (astropy.constants.c.value/1000.)
    datadic['sigmaerr_OIII2'] = np.array([93.0]) / 2.355    * 1666.150  / (astropy.constants.c.value/1000.)
    datadic['vshift_OIII2']   = np.array([-2.0])
    datadic['vshifterr_OIII2']= np.array([41.0])

    linename = 'OIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    # ---------------------------------------------------------------------------------
    if verbose: print('   Converting fluxes to 1e-20 erg/s/cm2 using fluxscale = '+str(fluxscale))
    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')

    return catreference, dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_sta14(fluxscale=1.0,verbose=True):
    """
    Data collected from Stark+2014

    NOTE ON DATA

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'sta14'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']      = np.array(['1.1','6.2','4.1','3.1','C4','C20b','881_329','899_340','883_357','860_359','885_354','863_348','876_330','869_328','854_344','854_362','846_340'])
    datadic['id']        = np.array([11,62,41,31,4,20,881329,899340,883357,860359,885354,863348,876330,869328,854344,854362,846340]) + baseid
    rasex                = np.array(['04:51:53.399','04:51:53.592','04:51:54.488','04:51:55.438','00:37:07.657','00:37:05.405','13:11:31.543','13:11:35.705','13:11:31.882','13:11:26.426','13:11:32.405','13:11:27.350','13:11:30.320','13:11:28.690','13:11:24.982','13:11:24.982','13:11:23.141'])
    decsex               = np.array(['+00:06:40.31','+00:06:24.96','+00:06:49.01','+00:06:41.16','+09:09:05.90','+09:09:59.14','-01:19:45.88','-01:20:25.22','-01:21:26.10','-01:21:31.22','-01:21:15.98','-01:20:54.82','-01:19:51.13','-01:19:42.69','-01:20:41.57','-01:21:43.53','-01:20:23.08'])
    datadic['ra']        = acoord.Angle(rasex, u.hour).degree
    datadic['dec']       = acoord.Angle(decsex, u.degree).degree
    datadic['redshift']  = np.array([2.060,1.405,1.810,1.904,2.622,2.689,1.559,1.599,1.702,1.702,1.705,1.834,1.834,2.543,2.663,2.731,2.976])
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    # Stark et al. (2014) table 2
    datadic['EW0_Lya']        = np.array([np.nan,np.nan,np.nan,np.nan,36.6,19.9*3./2.,np.nan,np.nan,76.0,163.8,35.7,73.1,50.0*3./2.,4.3,45.3,129.6,86.4])
    datadic['EW0err_Lya']     = np.array([np.nan,np.nan,np.nan,np.nan,5.7,-99,np.nan,np.nan,18.8,25.5,4.6,8.6,-99,0.4,2.7,21.8,24.8])

    datadic['EW0_CIII']       = np.array([6.7,1.2*3./2.,10.0,2.0,6.7,10.4*3./2.,7.1,5.1,6.5,12.4,3.9,13.5,10.0,1.8,4.0*3./2.,12.0,10.3*3./2.])
    datadic['EW0err_CIII']    = np.array([0.6,+99,2.4,0.7,2.1,-99,3.1,1.4,0.7,1.5,0.6,1.6,2.7,0.3,-99,3.2,-99])


    # Stark et al. (2014) table 4
    # flux ratios between CIII and NIV1487, NIII1750 also exists
    datadic['FR_NVCIII']      = np.array([ np.nan   , np.nan, np.nan, np.nan, np.nan, np.nan,  2.9    , np.nan  , np.nan  ,    0.6  , np.nan  ,   0.2   ,   1.1   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FRerr_NVCIII']   = np.array([ np.nan   , np.nan, np.nan, np.nan, np.nan, np.nan,  +99    , np.nan  , np.nan  ,    +99  , np.nan  ,   +99   ,   +99   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FR_CIVCIII']     = np.array([ np.nan   , np.nan, np.nan, np.nan, np.nan, np.nan,  0.5    , np.nan  , np.nan  ,    0.4  , np.nan  ,   0.8   ,   0.6   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FRerr_CIVCIII']  = np.array([ np.nan   , np.nan, np.nan, np.nan, np.nan, np.nan,  0.3    , np.nan  , np.nan  ,    0.1  , np.nan  ,   0.1   ,   +99   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FR_HeIICIII']    = np.array([ 0.5      , np.nan, np.nan, np.nan, np.nan, np.nan,  0.5    , np.nan  , np.nan  ,    0.2  , np.nan  ,   0.1   ,   0.5   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FRerr_HeIICIII'] = np.array([ 0.1      , np.nan, np.nan, np.nan, np.nan, np.nan,  +99    , np.nan  , np.nan  ,    +99  , np.nan  ,   +99   ,   +99   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FR_OIII1CIII']   = np.array([ 0.2      , np.nan, np.nan, np.nan, np.nan, np.nan,  0.5    , np.nan  , np.nan  ,    0.2  , np.nan  ,   0.2   ,   0.4   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FRerr_OIII1CIII']= np.array([ 0.1      , np.nan, np.nan, np.nan, np.nan, np.nan,  +99    , np.nan  , np.nan  ,    +99  , np.nan  ,   0.1   ,   +99   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FR_OIII2CIII']   = np.array([ 0.3      , np.nan, np.nan, np.nan, np.nan, np.nan,  0.5    , np.nan  , np.nan  ,    0.3  , np.nan  ,   0.6   ,   0.5   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FRerr_OIII2CIII']= np.array([ 0.1      , np.nan, np.nan, np.nan, np.nan, np.nan,  +99    , np.nan  , np.nan  ,    0.1  , np.nan  ,   0.1   ,   0.2   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FR_SiIII1CIII']= np.array([ 0.2      , np.nan, np.nan, np.nan, np.nan, np.nan,  0.5    , np.nan  , np.nan  ,    0.1  , np.nan  ,   0.2   ,   0.3   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FRerr_SiIII1CIII']= np.array([ 0.1      , np.nan, np.nan, np.nan, np.nan, np.nan,  +99    , np.nan  , np.nan  ,    0.1  , np.nan  ,   0.1   ,   +99   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FR_SiIII2CIII']= np.array([ 0.2      , np.nan, np.nan, np.nan, np.nan, np.nan,  0.5    , np.nan  , np.nan  ,    0.1  , np.nan  ,   0.1   ,   0.3   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FRerr_SiIII2CIII']= np.array([ 0.1      , np.nan, np.nan, np.nan, np.nan, np.nan,  +99    , np.nan  , np.nan  ,    0.1  , np.nan  ,   +99   ,   +99   , np.nan  , np.nan  , np.nan  , np.nan])

    # etimated manually from the entries in table 4
    datadic['FR_OIIICIII']    = np.array([ 0.5      , np.nan, np.nan, np.nan, np.nan, np.nan,  1.0    , np.nan  , np.nan  ,    0.5  , np.nan  ,   0.8   ,   0.9   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FRerr_OIIICIII'] = np.array([ 0.14     , np.nan, np.nan, np.nan, np.nan, np.nan,  +99    , np.nan  , np.nan  ,    +99  , np.nan  ,   0.14  ,   +99   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FR_SiIIICIII']    = np.array([ 0.4      , np.nan, np.nan, np.nan, np.nan, np.nan,  1.0    , np.nan  , np.nan  ,    0.2  , np.nan  ,   0.3   ,   0.6   , np.nan  , np.nan  , np.nan  , np.nan])
    datadic['FRerr_SiIIICIII'] = np.array([ 0.14     , np.nan, np.nan, np.nan, np.nan, np.nan,  +99    , np.nan  , np.nan  ,    0.14 , np.nan  ,   +99   ,   +99   , np.nan  , np.nan  , np.nan  , np.nan])

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_sta15(fluxscale=1.0,verbose=True):
    """
    Data collected from Stark et al. (2015a,b, 2017)

    GN-108036 Lya measurements from Ono+12.
    A1703-zd4 and A1703-zD1 limits on UV lines taken from Mainali+18.
    CIII measurements from A1703-zd6 also taken from Mainali+18, who has it from Stark et al. (2017).

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'sta15'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']      = np.array(['A383-5.2','GN-108036','A1703-zd1a','A1703-zd4','A1703-zd6','A1703-23'])
    datadic['id']        = np.array([38352,108036,17031,17034,17036,170323]) + baseid
    rasex                = np.array(['02:48:04.600','12:36:22.68','13:14:59.418','13:15:07.189','13:15:01.007','13:15:01.469'])
    decsex               = np.array([['-03:31:58.47','+62:08:07.5','+51:50:00.8','+51:50:23.6','+51:50:04.4','+51:48:26.5']])
    datadic['ra']        = acoord.Angle(rasex, u.hour).degree
    datadic['dec']       = acoord.Angle(decsex, u.degree).degree
    datadic['redshift']  = np.array([6.0294,7.213  ,6.75,8.15,7.043,5.828])
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_Lya']          = np.array([11000    , 2500.0  , np.nan   , np.nan  , 2840.0  , 2500.0])
    datadic['ferr_Lya']       = np.array([0.0      , 0.0     , np.nan   , np.nan  , 530.0   , 0.0])
    datadic['EW0_Lya']        = np.array([138.0    , 33.0    , np.nan   , np.nan  , 65.0    , 10.3])
    datadic['EW0err_Lya']     = np.array([0.0      , 0.0     , np.nan   , np.nan  , 12.0    , 0.0])
    datadic['sigma_Lya']      = np.array([np.nan   , 15.0    , np.nan   , np.nan  , np.nan  , np.nan])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_CIII1']        = np.array([518.0    , 130.0   , np.nan   , np.nan  , 330.0   , 210.0 ])
    datadic['ferr_CIII1']     = np.array([150.0    , 40.0    , np.nan   , np.nan  , 99      , 99    ])
    datadic['EW0_CIII1']      = np.array([np.nan   , np.nan  , np.nan   , np.nan  , 19.8    , 2.7   ])
    datadic['EW0err_CIII1']   = np.array([np.nan   , np.nan  , np.nan   , np.nan  , 99      , 99    ])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_CIII2']        = np.array([370.     , 90.0    , np.nan   , np.nan  , 330.0   , 336.0 ])
    datadic['ferr_CIII2']     = np.array([110.     , 30.0    , np.nan   , np.nan  , 99      , 99    ])
    datadic['EW0_CIII2']      = np.array([np.nan   , np.nan  , np.nan   , np.nan  , 19.8    , 4.3   ])
    datadic['EW0err_CIII2']   = np.array([np.nan   , np.nan  , np.nan   , np.nan  , 99      , 99    ])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_CIII']         = np.array([890.     , 220.0   , 530.0    , 570.0   , np.nan  , np.nan])
    datadic['ferr_CIII']      = np.array([270      , 50.0    , 99.0     , 99.0    , np.nan  , np.nan])
    datadic['EW0_CIII']       = np.array([22.5     , 7.6     , 6.8      , 20.0    , np.nan  , np.nan])
    datadic['EW0err_CIII']    = np.array([7.1      , 2.8     , 99.0     , 99.0    , np.nan  , np.nan])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_CIV1']         = np.array([np.nan   , np.nan  , 180.0    , 570.0   , 410.0   , 213.0 ])
    datadic['ferr_CIV1']      = np.array([np.nan   , np.nan  , 99       , 99.0    , 60.0    , 99    ])
    datadic['EW0_CIV1']       = np.array([np.nan   , np.nan  , 1.38     , 20.0    , 19.9    , np.nan])
    datadic['EW0err_CIV1']    = np.array([np.nan   , np.nan  , 99       , 99.0    , 3.6     , np.nan])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_CIV2']         = np.array([np.nan   , np.nan  , 180.0    , 570.0   , 380.0   , 213.0 ])
    datadic['ferr_CIV2']      = np.array([np.nan   , np.nan  , 99       , 99.0    , 90.0    , 99    ])
    datadic['EW0_CIV2']       = np.array([np.nan   , np.nan  , 1.38     , 20.0    , 18.1    , np.nan])
    datadic['EW0err_CIV2']    = np.array([np.nan   , np.nan  , 99       , 99.0    , 4.6     , np.nan])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    linename = 'CIV'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_HeII']         = np.array([np.nan   , np.nan  , np.nan   , 570.0   , 210.0   , np.nan])
    datadic['ferr_HeII']      = np.array([np.nan   , np.nan  , np.nan   , 99.0    , 99      , np.nan])
    datadic['EW0_HeII']       = np.array([np.nan   , np.nan  , np.nan   , 20.0    , 11.4    , np.nan])
    datadic['EW0err_HeII']    = np.array([np.nan   , np.nan  , np.nan   , 99.0    , 99      , np.nan])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_OIII1']        = np.array([np.nan   , np.nan  , np.nan   , 570.0   , 180.0   , np.nan])
    datadic['ferr_OIII1']     = np.array([np.nan   , np.nan  , np.nan   , 99.0    , 70.0    , np.nan])
    datadic['EW0_OIII1']      = np.array([np.nan   , np.nan  , np.nan   , 20.0    , 9.8     , np.nan])
    datadic['EW0err_OIII1']   = np.array([np.nan   , np.nan  , np.nan   , 99.0    , 3.9     , np.nan])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_OIII2']        = np.array([np.nan   , np.nan  , np.nan   , 570.0   , np.nan  , np.nan])
    datadic['ferr_OIII2']     = np.array([np.nan   , np.nan  , np.nan   , 99.0    , np.nan  , np.nan])
    datadic['EW0_OIII2']      = np.array([np.nan   , np.nan  , np.nan   , 20.0    , np.nan  , np.nan])
    datadic['EW0err_OIII2']   = np.array([np.nan   , np.nan  , np.nan   , 99.0    , np.nan  , np.nan])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    linename = 'OIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # ---------------------------------------------------------------------------------
    if verbose: print('   Converting fluxes to 1e-20 erg/s/cm2 using fluxscale = '+str(fluxscale))
    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_mai18(fluxscale=1e2,verbose=True):
    """
    Data collected from Mainali et al (2018) and Stark et al. (2017)

    Not including A1703-zd4, A1703-zD1 and A1703-zd6 as those are included in the sta15 lists.
    Lya  info for EGS-zs8-1 from Oesch+15
    Lya  info for EGS-zs8-2 from Roberts-Borsani+16 and Stark+17
    Lya  info for EGSY8p7   from Zitrin+15
    CIII info for EGS-zs8-1 from Stark+17
    CIII info for EGS-zs8-2 from Stark+17

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'mai18'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']      = np.array(['EGSY8p7','EGS-zs8-1','EGS-zs8-2','A383-2211','Abell2218-S3.a','Abell2218-S3.b','J14144+5446','Abell2218-C1.b','Abell2218-C1.c'])
    datadic['id']        = np.array([87,81,82,383211,221831,221832,5446,221812,221813]) + baseid
    rasex                = np.array(['14:20:08.50','14:20:34.89','14:20:12.09','02:48:01.39','16:35:51.96','16:35:52.08','14:14:46.82','16:35:54.40','16:35:48.92'])
    decsex               = np.array(['+52:53:26.6','+53:00:15.4','+53:00:27.0','-03:32:58.4','+66:12:45.9','+66:12:51.8','+54:46:31.9','+66:12:32.8','+66:12:02.4'])
    datadic['ra']        = acoord.Angle(rasex, u.hour).degree
    datadic['dec']       = acoord.Angle(decsex, u.degree).degree
    datadic['redshift']  = np.array([8.683,7.730,7.477,6.031,5.576,5.576,5.426,6.7,6.7])
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    #                               ['EGSY8p7','EGS-zs8-1','EGS-zs8-2','A383-2211','Abell2218-S3.a','Abell2218-S3.b','J14144+5446','Abell2218-C1.b','Abell2218-C1.c'])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_Lya']         = np.array([17.0     ,  17.0     ,  7.4      ,  16.4     ,  np.nan        ,  np.nan        ,  np.nan     ,  np.nan        ,  np.nan        ])
    datadic['ferr_Lya']      = np.array([0.0      ,  3.0      ,  1.0      ,  1.4      ,  np.nan        ,  np.nan        ,  np.nan     ,  np.nan        ,  np.nan        ])
    datadic['EW0_Lya']       = np.array([28.0     ,  21.0     ,  8.3      ,  21.1     ,  239.0         ,  239.0         ,  260.0      ,  np.nan        ,  np.nan        ])
    datadic['EW0err_Lya']    = np.array([0.0      ,  4.0      ,  1.4      ,  2.8      ,  25.0          ,  25.0          ,  0.0        ,  np.nan        ,  np.nan        ])
    datadic['sigma_Lya']     = np.array([np.nan   ,  13.0     ,  np.nan   ,  np.nan   ,  np.nan        ,  np.nan        ,  np.nan     ,  np.nan        ,  np.nan        ]) / 2.355
    datadic['sigmaerr_Lya']  = np.array([np.nan   ,  3.0      ,  np.nan   ,  np.nan   ,  np.nan        ,  np.nan        ,  np.nan     ,  np.nan        ,  np.nan        ]) / 2.355
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_NV1']         = np.array([np.nan   ,  6.0      ,  1.7      ,  3.8      ,  np.nan        ,  np.nan        ,  np.nan     ,  np.nan        ,  np.nan        ])
    datadic['ferr_NV1']      = np.array([np.nan   ,  +99      ,  +99      ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  np.nan        ,  np.nan        ])
    datadic['EW0_NV1']       = np.array([np.nan   ,  7.4      ,  2.2      ,  4.9      ,  np.nan        ,  np.nan        ,  np.nan     ,  np.nan        ,  np.nan        ])
    datadic['EW0err_NV1']    = np.array([np.nan   ,  +99      ,  +99      ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  np.nan        ,  np.nan        ])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_NV2']         = np.array([2.8      ,  24.3     ,  3.6      ,  3.8      ,  np.nan        ,  np.nan        ,  np.nan     ,  np.nan        ,  np.nan        ])
    datadic['ferr_NV2']      = np.array([0.6      ,  +99      ,  +99      ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  np.nan        ,  np.nan        ])
    datadic['EW0_NV2']       = np.array([4.2      ,  29.6     ,  4.5      ,  4.9      ,  np.nan        ,  np.nan        ,  np.nan     ,  np.nan        ,  np.nan        ])
    datadic['EW0err_NV2']    = np.array([0.9      ,  +99      ,  +99      ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  np.nan        ,  np.nan        ])

    linename = 'NV'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_CIV1']        = np.array([1.6      ,  11.0     ,  1.2      ,  2.5      ,  np.nan        ,  np.nan        ,  np.nan     ,  7.1           ,  7.1           ])
    datadic['ferr_CIV1']     = np.array([+99      ,  +99      ,  +99      ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  +99           ,  +99           ])
    datadic['EW0_CIV1']      = np.array([4.6      ,  17.0     ,  2.4      ,  4.6      ,  np.nan        ,  np.nan        ,  np.nan     ,  5.5           ,  5.5           ])
    datadic['EW0err_CIV1']   = np.array([+99      ,  +99      ,  +99      ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  +99           ,  +99           ])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - -
    datadic['f_CIV2']        = np.array([1.6      ,  np.nan   ,  1.2      ,  1.3      ,  np.nan        ,  np.nan        ,  np.nan     ,  7.1           ,  7.1           ])
    datadic['ferr_CIV2']     = np.array([+99      ,  np.nan   ,  +99      ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  +99           ,  +99           ])
    datadic['EW0_CIV2']      = np.array([4.6      ,  np.nan   ,  2.4      ,  2.3      ,  np.nan        ,  np.nan        ,  np.nan     ,  5.5           ,  5.5           ])
    datadic['EW0err_CIV2']   = np.array([+99      ,  np.nan   ,  +99      ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  +99           ,  +99           ])

    linename = 'CIV'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - -
    datadic['f_HeII']        = np.array([5.3      ,  np.nan   ,  np.nan   ,  5.3      ,  np.nan        ,  np.nan        ,  np.nan     ,  7.1           ,  7.1           ])
    datadic['ferr_HeII']     = np.array([+99      ,  np.nan   ,  np.nan   ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  +99           ,  +99           ])
    datadic['EW0_HeII']      = np.array([14.8     ,  np.nan   ,  np.nan   ,  12.2     ,  np.nan        ,  np.nan        ,  np.nan     ,  5.5           ,  5.5           ])
    datadic['EW0err_HeII']   = np.array([+99      ,  np.nan   ,  np.nan   ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  +99           ,  +99           ])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - -
    datadic['f_OIII1']       = np.array([1.5      ,  np.nan   ,  np.nan   ,  2.5      ,  np.nan        ,  np.nan        ,  np.nan     ,  7.1           ,  7.1           ])
    datadic['ferr_OIII1']    = np.array([+99      ,  np.nan   ,  np.nan   ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  +99           ,  +99           ])
    datadic['EW0_OIII1']     = np.array([4.6      ,  np.nan   ,  np.nan   ,  5.7      ,  np.nan        ,  np.nan        ,  np.nan     ,  5.5           ,  5.5           ])
    datadic['EW0err_OIII1']  = np.array([+99      ,  np.nan   ,  np.nan   ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  +99           ,  +99           ])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - -
    datadic['f_OIII2']       = np.array([1.5      ,  np.nan   ,  np.nan   ,  4.9      ,  np.nan        ,  np.nan        ,  np.nan     ,  7.1           ,  7.1           ])
    datadic['ferr_OIII2']    = np.array([+99      ,  np.nan   ,  np.nan   ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  +99           ,  +99           ])
    datadic['EW0_OIII2']     = np.array([4.6      ,  np.nan   ,  np.nan   ,  11.2     ,  np.nan        ,  np.nan        ,  np.nan     ,  5.5           ,  5.5           ])
    datadic['EW0err_OIII2']  = np.array([+99      ,  np.nan   ,  np.nan   ,  +99      ,  np.nan        ,  np.nan        ,  np.nan     ,  +99           ,  +99           ])

    linename = 'OIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_CIII1']       = np.array([np.nan   ,  4.5      ,  2.3      ,  1.7      ,  3.2           ,  3.2           ,  12.0       ,  np.nan        ,  np.nan        ])
    datadic['ferr_CIII1']    = np.array([np.nan   ,  0.5      ,  +99      ,  +99      ,  +99           ,  +99           ,  +99        ,  np.nan        ,  np.nan        ])
    datadic['EW0_CIII1']     = np.array([np.nan   ,  12.0     ,  7.1      ,  5.8      ,  15.6          ,  15.6          ,  7.3        ,  np.nan        ,  np.nan        ])
    datadic['EW0err_CIII1']  = np.array([np.nan   ,  2.0      ,  +99      ,  +99      ,  +99           ,  +99           ,  +99        ,  np.nan        ,  np.nan        ])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - -  - - - - - - - - - - - - - - - - - - - - - - - -
    datadic['f_CIII2']       = np.array([np.nan   ,  3.6      ,  2.3      ,  1.2      ,  3.2           ,  3.2           ,  8.4        ,  np.nan        ,  np.nan        ])
    datadic['ferr_CIII2']    = np.array([np.nan   ,  0.5      ,  +99      ,  +99      ,  +99           ,  +99           ,  +99        ,  np.nan        ,  np.nan        ])
    datadic['EW0_CIII2']     = np.array([np.nan   ,  10.0     ,  7.1      ,  4.1      ,  15.6          ,  15.6          ,  5.2        ,  np.nan        ,  np.nan        ])
    datadic['EW0err_CIII2']  = np.array([np.nan   ,  1.0      ,  +99      ,  +99      ,  +99           ,  +99           ,  +99        ,  np.nan        ,  np.nan        ])

    linename = 'CIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    # ---------------------------------------------------------------------------------
    if verbose: print('   Converting fluxes to 1e-20 erg/s/cm2 using fluxscale = '+str(fluxscale))
    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_sha03(fluxscale=1.0,verbose=True):
    """
    Data collected from Shapley+2003

    NOTE ON DATA

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'sha03'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']      = np.array(['LAEgr1','LAEgr2','LAEgr3','LAEgr4'])
    datadic['id']        = np.array([1,2,3,4]) + baseid
    datadic['ra']        = np.array([np.nan,np.nan,np.nan,np.nan])
    datadic['dec']       = np.array([np.nan,np.nan,np.nan,np.nan])
    datadic['redshift']  = np.array([0.0,0.0,0.0,0.0])
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    datadic['EW0_Lya']        = np.array([-14.92  ,  -1.10      ,  11.0          ,  52.63        ])
    datadic['EW0err_Lya']     = np.array([0.56    ,  0.38       ,  0.71          ,  2.74        ])

    datadic['EW0_OIII']       = np.array([0.23    ,  0.01       ,  0.10          ,  1.16        ])
    datadic['EW0err_OIII']    = np.array([0.19    ,  0.16       ,  0.12          ,  0.56        ])

    datadic['EW0_CIII2']      = np.array([0.41    ,  2.89       ,  1.90          ,  5.37        ])
    datadic['EW0err_CIII2']   = np.array([0.39    ,  1.04       ,  1.07          ,  2.99        ])

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_bay14(fluxscale=1e2,verbose=True):
    """
    Data collected from Bayliss+2014

    It is not clear if the OIII4363 is a 3sigma limit.
    Including only the measrued fluxes (corrected for lensing magnification). De-reddened fluxes also available.

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'bay14'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']      = np.array(['SGAS_J105039.6+001730'])
    datadic['id']        = np.array([1]) + baseid
    rasex                = np.array(['10:50:41.336'])
    decsex               = np.array(['+00:17:23.35'])
    datadic['ra']        = acoord.Angle(rasex, u.hour).degree
    datadic['dec']       = acoord.Angle(decsex, u.degree).degree
    datadic['redshift']  = np.array([3.6253])
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    datadic['f_HeII']      = np.array([0.29   ])
    datadic['f_OIII1']     = np.array([0.10   ])
    datadic['f_OIII2']     = np.array([0.29   ])
    datadic['f_NIII']      = np.array([0.08   ])
    datadic['f_SiIII1']    = np.array([0.26   ])
    datadic['f_SiIII2']    = np.array([0.15   ])
    datadic['f_CIII1']     = np.array([1.14   ])
    datadic['f_CIII2']     = np.array([0.69   ])
    datadic['f_OII1']      = np.array([2.26   ])
    datadic['f_OII2']      = np.array([2.22   ])
    datadic['f_NeIII']     = np.array([1.83   ])
    datadic['f_Hgamma']    = np.array([2.96   ])
    datadic['f_OIII4363']  = np.array([0.81   ])
    datadic['f_Hbeta']     = np.array([5.71   ])
    datadic['f_OIII4960']  = np.array([13.44  ])
    datadic['f_OIII5008']  = np.array([46.12  ])

    datadic['ferr_HeII']      = np.array([   0.04  ])
    datadic['ferr_OIII1']     = np.array([   0.04  ])
    datadic['ferr_OIII2']     = np.array([   0.04  ])
    datadic['ferr_NIII']      = np.array([   0.03  ])
    datadic['ferr_SiIII1']    = np.array([   0.04  ])
    datadic['ferr_SiIII2']    = np.array([   0.04  ])
    datadic['ferr_CIII1']     = np.array([   0.07  ])
    datadic['ferr_CIII2']     = np.array([   0.08  ])
    datadic['ferr_OII1']      = np.array([   0.06  ])
    datadic['ferr_OII2']      = np.array([   0.06  ])
    datadic['ferr_NeIII']     = np.array([   0.04  ])
    datadic['ferr_Hgamma']    = np.array([   0.12  ])
    datadic['ferr_OIII4363']  = np.array([   +99   ])
    datadic['ferr_Hbeta']     = np.array([   0.06  ])
    datadic['ferr_OIII4960']  = np.array([   0.10  ])
    datadic['ferr_OIII5008']  = np.array([   0.12  ])

    linename = 'SiIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'])

    linename = 'OIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'])

    linename = 'CIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'])

    # ---------------------------------------------------------------------------------
    if verbose: print('   Converting fluxes to 1e-20 erg/s/cm2 using fluxscale = '+str(fluxscale))
    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_sch16(fluxscale=1e3,verbose=True):
    """
    Data collected from Schmidt+2016

    Stacked values published in Matthee+17 with f_lya calculated from Schmidt+16 fluxes.
    RXJ1347_01037 has DEIMOS Lya from Huang+2016a

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'sch16'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']      = np.array(['GLASSstack_z7p2','RXJ1347_01037','MACS2129_00899'])
    datadic['id']        = np.array([1,2,3]) + baseid
    datadic['ra']        = np.array([np.nan, 206.900859670, 322.343220360])
    datadic['dec']       = np.array([np.nan, -11.754209621, -7.693382243])
    datadic['redshift']  = np.array([7.2, 6.76, 8.10])
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    # Fluxes used for estimating GLASS stack f_lya
    # f_lya  ferr_lya
    # 1.91   0.7   A2744-135-150223_00463-G102.2D.fits
    # 2.76   0.94  A2744-135-150223_00844-G102.2D.fits
    # 2.0    0.56  MACS1423.8+2404-008-150223_00648-G102.2D.fits
    # 1.87   0.63  MACS1423.8+2404-008-150223_01102-G102.2D.fits
    # 3.45   0.87  MACS2129.4-0741-050-150223_00677-G102.2D.fits
    # 3.42   0.70  MACS2129.4-0741-328-150223_00677-G102.2D.fits
    # 0.74   0.52  MACS2129.4-0741-050-150223_00899-G102.2D.fits
    # 1.26   0.47  MACS2129.4-0741-050-150223_00899-G141.2D.fits
    # 2.7    0.85  MACS2129.4-0741-050-150223_01516-G102.2D.fits
    # 2.55   1.07  RXJ2248-053-150223_00207-G141.2D.fits

    f_lya_stack    = np.sum(np.array([1.91,2.76,2.0 ,1.87,3.45,3.42,0.74,1.26,2.7,2.55]))/10
    ferr_lya_stack = np.sqrt(np.sum(np.array([0.7 ,0.94,0.56,0.63,0.87,0.70,0.52,0.47,0.85,1.07])**2))/10


    datadic['f_Lya']          = np.array([f_lya_stack,     2.6,    1.00])
    datadic['ferr_Lya']       = np.array([ferr_lya_stack,  0.4,    0.25])
    datadic['EW0_Lya']        = np.array([np.nan,          74.33, 59.0])
    datadic['EW0err_Lya']     = np.array([np.nan,          12.09, 21.0])

    datadic['f_CIV']          = np.array([0.3*datadic['f_Lya'][0], 0.36*datadic['f_Lya'][1], 0.64*datadic['f_Lya'][2]])*3./2.
    datadic['ferr_CIV']       = np.array([+99, +99, +99 ])
    datadic['EW0_CIV']        = np.array([0.3*datadic['EW0_Lya'][0], 0.36*datadic['EW0_Lya'][1], 0.64*datadic['EW0_Lya'][2]])*3./2.
    datadic['EW0err_CIV']     = np.array([+99, +99, +99 ])

    datadic['f_CIII']         = np.array([.2*datadic['f_Lya'][0]*3./2., 0.25 *datadic['f_Lya'][1]*3./2. , np.nan ])
    datadic['ferr_CIII']      = np.array([+99, +99, np.nan ])
    datadic['EW0_CIII']       = np.array([.2*datadic['EW0_Lya'][0]*3./2., 0.25 *datadic['EW0_Lya'][1]*3./2. , np.nan ])
    datadic['EW0err_CIII']    = np.array([+99, +99, np.nan ])

    datadic['f_HeII']         = np.array([0.2*datadic['f_Lya'][0] * 3./2. , np.nan, np.nan ])
    datadic['ferr_HeII']      = np.array([+99, np.nan, np.nan ])
    datadic['EW0_HeII']       = np.array([0.2*datadic['EW0_Lya'][0] * 3./2. , np.nan, np.nan ])
    datadic['EW0err_HeII']    = np.array([+99, np.nan, np.nan ])

    datadic['f_OIII']         = np.array([0.2*datadic['f_Lya'][0] * 3./2. , np.nan, np.nan ])
    datadic['ferr_OIII']      = np.array([+99, np.nan, np.nan ])
    datadic['EW0_OIII']       = np.array([0.2*datadic['EW0_Lya'][0] * 3./2. , np.nan, np.nan ])
    datadic['EW0err_OIII']    = np.array([+99, np.nan, np.nan ])

    datadic['f_NV']           = np.array([0.4*datadic['f_Lya'][0] * 3./2. , np.nan, np.nan ])
    datadic['ferr_NV']        = np.array([+99, np.nan, np.nan ])
    datadic['EW0_NV']         = np.array([0.4*datadic['EW0_Lya'][0] * 3./2. , np.nan, np.nan ])
    datadic['EW0err_NV']      = np.array([+99, np.nan, np.nan ])

    # ---------------------------------------------------------------------------------
    if verbose: print('   Converting fluxes to 1e-20 erg/s/cm2 using fluxscale = '+str(fluxscale))
    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_lef19(fluxscale=1e2,verbose=True):
    """
    Data collected from Le Fevre+2019 VUDS CIII emitter stacks

    Guessing fluxes to be 1e-18 as that's what they are in Amorin+2017
    ----> Fluxes are relative to CIII which is set to one  <----
    Individual objects not provided in paper...

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'lef19'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']      = np.array(['ewGT20','ewLT20GT20','ewLT10GT5','AGNtypeII','SFGs_z2.5','SFGs_z3.5'])
    datadic['id']        = np.array([1,2,3,4,5,6]) + baseid
    # rasex                = np.array(['04:22:00.81'])
    # decsex               = np.array(['-38:37:03.59'])
    # datadic['ra']        = acoord.Angle(rasex, u.hour).degree
    # datadic['dec']       = acoord.Angle(decsex, u.degree).degree
    datadic['ra']        = np.array([np.nan,np.nan,np.nan,np.nan,np.nan,np.nan])
    datadic['dec']       = np.array([np.nan,np.nan,np.nan,np.nan,np.nan,np.nan])
    datadic['redshift']  = np.array([3.0,3.0,3.0,3.0,2.5,3.5])
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    # ----> Fluxes are relative to CIII which is set to one; guessing fluxes to be 1e-18 as that's what they are in Amorin+2017 <----
    # datadic['f_CIII']         =   np.array([  1.00  , 1.00    , 1.00    , 1.00  , 1.0     , 1.0     ])
    # datadic['f_Lya']          =   np.array([  4.6  *datadic['f_CIII'][0] , 6.40 *datadic['f_CIII'][1]     , 5.2    *datadic['f_CIII'][2]   , 6.25 *datadic['f_CIII'][3]   , 3.5    *datadic['f_CIII'][4]   , 14.8   *datadic['f_CIII'][5]   ])
    # datadic['f_NV']           =   np.array([  0.05 *datadic['f_CIII'][0] , 0.02 *datadic['f_CIII'][1]     , np.nan *datadic['f_CIII'][2]   , 1.31 *datadic['f_CIII'][3]   , np.nan *datadic['f_CIII'][4]   , np.nan *datadic['f_CIII'][5]   ])
    # datadic['f_SiII']         =   np.array([  0.07 *datadic['f_CIII'][0] , 0.04 *datadic['f_CIII'][1]     , np.nan *datadic['f_CIII'][2]   , 0.12 *datadic['f_CIII'][3]   , np.nan *datadic['f_CIII'][4]   , np.nan *datadic['f_CIII'][5]   ])
    # datadic['f_SiIV']         =   np.array([  0.13 *datadic['f_CIII'][0] , 0.02 *datadic['f_CIII'][1]     , np.nan *datadic['f_CIII'][2]   , 0.14 *datadic['f_CIII'][3]   , np.nan *datadic['f_CIII'][4]   , np.nan *datadic['f_CIII'][5]   ])
    # datadic['f_NIV']          =   np.array([  0.19 *datadic['f_CIII'][0] , 0.08 *datadic['f_CIII'][1]     , np.nan *datadic['f_CIII'][2]   , 0.19 *datadic['f_CIII'][3]   , np.nan *datadic['f_CIII'][4]   , np.nan *datadic['f_CIII'][5]   ])
    # datadic['f_CIV']          =   np.array([  0.31 *datadic['f_CIII'][0] , 0.36 *datadic['f_CIII'][1]     , 0.01   *datadic['f_CIII'][2]   , 3.00 *datadic['f_CIII'][3]   , -2.1   *datadic['f_CIII'][4]   , -1.2   *datadic['f_CIII'][5]   ])
    # datadic['f_HeII']         =   np.array([  0.27 *datadic['f_CIII'][0] , 0.49 *datadic['f_CIII'][1]     , 0.30   *datadic['f_CIII'][2]   , 0.92 *datadic['f_CIII'][3]   , np.nan *datadic['f_CIII'][4]   , np.nan *datadic['f_CIII'][5]   ])
    # datadic['f_OIII']         =   np.array([  0.27 *datadic['f_CIII'][0] , 0.14 *datadic['f_CIII'][1]     , np.nan *datadic['f_CIII'][2]   , 0.13 *datadic['f_CIII'][3]   , np.nan *datadic['f_CIII'][4]   , np.nan *datadic['f_CIII'][5]   ])
    # datadic['f_NIII']         =   np.array([  0.07 *datadic['f_CIII'][0] , 0.04 *datadic['f_CIII'][1]     , np.nan *datadic['f_CIII'][2]   , 0.13 *datadic['f_CIII'][3]   , np.nan *datadic['f_CIII'][4]   , np.nan *datadic['f_CIII'][5]   ])
    # datadic['f_SiIII']        =   np.array([  0.18 *datadic['f_CIII'][0] , 0.31 *datadic['f_CIII'][1]     , np.nan *datadic['f_CIII'][2]   , 0.19 *datadic['f_CIII'][3]   , np.nan *datadic['f_CIII'][4]   , np.nan *datadic['f_CIII'][5]   ])
    #
    #
    # datadic['ferr_CIII']      =   np.array([  0.05  ,  0.04   , 0.15    , 0.09  ,  0.05   , 0.05    ])
    # datadic['ferr_Lya']       =   np.array([  0.07    ,  0.04    , 0.2      , 0.06   ,  0.1     , 0.12      ])
    # datadic['ferr_NV']        =   np.array([  0.02    ,  0.01    , np.nan   , 0.10   ,  np.nan  , np.nan    ])
    # datadic['ferr_SiII']      =   np.array([  0.02    ,  0.02    , np.nan   , 0.04   ,  np.nan  , np.nan    ])
    # datadic['ferr_SiIV']      =   np.array([  0.03    ,  0.02    , np.nan   , 0.04   ,  np.nan  , np.nan    ])
    # datadic['ferr_NIV']       =   np.array([  0.02    ,  0.02    , np.nan   , 0.05   ,  np.nan  , np.nan    ])
    # datadic['ferr_CIV']       =   np.array([  0.03    ,  0.02    , 0.1      , 0.08   ,  0.05    , 0.05      ])
    # datadic['ferr_HeII']      =   np.array([  0.03    ,  0.03    , 0.06     , 0.06   ,  np.nan  , np.nan    ])
    # datadic['ferr_OIII']      =   np.array([  0.03    ,  0.02    , np.nan   , 0.05   ,  np.nan  , np.nan    ])
    # datadic['ferr_NIII']      =   np.array([  0.03    ,  0.02    , np.nan   , 0.06   ,  np.nan  , np.nan    ])
    # datadic['ferr_SiIII']     =   np.array([  0.12    ,  0.10    , np.nan   , 0.02   ,  np.nan  , np.nan    ])

    datadic['EW0_CIII']       =   np.array([  23.5  , 11.1   , 7.1     , 13.2  , 2.0     , 2.24    ])
    datadic['EW0_Lya']        =   np.array([  71.7  , 58.8   , 30.3    , 98.9  , 7.9     , 14.8    ])
    datadic['EW0_NV']         =   np.array([  0.9   , 0.3    , np.nan  , 19.1  , np.nan  , np.nan  ])
    datadic['EW0_SiII']       =   np.array([  0.6   , 0.2    , np.nan  , 1.7   , np.nan  , np.nan  ])
    datadic['EW0_SiIV']       =   np.array([  1.4   , 0.1    , np.nan  , 1.9   , np.nan  , np.nan  ])
    datadic['EW0_NIV']        =   np.array([  2.2   , 0.6    , np.nan  , 2.5   , np.nan  , np.nan  ])
    datadic['EW0_CIV']        =   np.array([  4.4   , 2.4    , 0.1     , 38.8  , 2.78    , 2.84    ])
    datadic['EW0_HeII']       =   np.array([  4.3   , 3.2    , 1.7     , 11.7  , np.nan  , np.nan  ])
    datadic['EW0_OIII']       =   np.array([  5.5   , 1.4    , 0.7     , 1.7   , np.nan  , np.nan  ])
    datadic['EW0_NIII']       =   np.array([  1.1   , 0.4    , np.nan  , 1.5   , np.nan  , np.nan  ])
    datadic['EW0_SiIII']      =   np.array([  5.0   , 4.0    , np.nan  , 2.5   , np.nan  , np.nan  ])

    datadic['EW0err_CIII']    =   np.array([  1.8   , 1.2    , 1.3     , 1.2   , 0.2     , 0.3    ])
    datadic['EW0err_Lya']     =   np.array([  5.0   , 7.0    , 2.30    , 3.0   , 1.3     , 1.5     ])
    datadic['EW0err_NV']      =   np.array([  0.3   , 0.15   , np.nan  , 2.0   , np.nan  , np.nan ])
    datadic['EW0err_SiII']    =   np.array([  0.2   , 0.15   , np.nan  , 0.4   , np.nan  , np.nan ])
    datadic['EW0err_SiIV']    =   np.array([  0.2   , 0.15   , np.nan  , 0.5   , np.nan  , np.nan ])
    datadic['EW0err_NIV']     =   np.array([  0.3   , 0.10   , np.nan  , 0.5   , np.nan  , np.nan ])
    datadic['EW0err_CIV']     =   np.array([  0.6   , 0.5    , 0.3     , 2.1   , 0.11    , 0.15   ])
    datadic['EW0err_HeII']    =   np.array([  0.7   , 0.5    , 0.2     , 1.0   , np.nan  , np.nan ])
    datadic['EW0err_OIII']    =   np.array([  0.6   , 0.4    , 0.2     , 0.6   , np.nan  , np.nan ])
    datadic['EW0err_NIII']    =   np.array([  0.2   , 0.2    , np.nan  , 0.7   , np.nan  , np.nan ])
    datadic['EW0err_SiIII']   =   np.array([  0.4   , 0.5    , np.nan  , 0.9   , np.nan  , np.nan ])

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_amo17(fluxscale=1e2,verbose=True):
    """
    Data collected from Amorin+2017

    NOTE ON DATA

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'amo17'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']      = np.array(['5100534435',   '5100565880',   '5100750978',   '5100994378',   '5100998761',   '5101421970',   '5101444192',   '510583858' ,   '511267982' ,   '510838687', 'composite' ])
    datadic['id']        = np.array([   534435,      565880,     750978,     994378,     998761,     1421970,     1444192,     583858 ,     1267982 ,     838687,      1]) + baseid
    datadic['ra']        = np.array([np.nan]*11) # coordinates not public but see email 170207
    datadic['dec']       = np.array([np.nan]*11) # coordinates not public but see email 170207
    datadic['redshift']  = np.array([2.9635,3.0505,2.963 ,2.797,2.446,2.465 ,3.424,2.4141,2.8256,2.5539,np.nan])
    datadic['reference'] = [catreference]*len(datadic['id'])
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    datadic['EW0_Lya']        =   np.array([  111.0,   35.0,   56.0,  47.0,  168.0,  82.0,  89.0,  185.0,  261.0,  267.0,  131.0,   ])
    datadic['EW0err_Lya']     =   np.array([  16.0,   12.0,   12.0,  11.0,  10.0,  12.0,  10.0,  35.0,  30.0,  35.0,  11.0,   ])

    datadic['f_Lya']          =   np.array([  126.3,   19.5,   67.7,  44.0,  567.1,  120.9,  377.0,  361.5,  738.4,  182.3,  260.0,   ])
    datadic['ferr_Lya']       =   np.array([  6.8,   3.2,   4.2,  4.0,  8.1,  5.9,  11.7,  21.4,  31.5,  3.5,  2.0,   ])

    datadic['f_CIV']        =   np.array([  4.1,5.7,4.6,1.4,41.3,1.6,6.1,8.3,27.7,3.5,9.1   ])
    datadic['ferr_CIV']     =   np.array([1.4,1.1,1.1,+99,5.3,0.8,1.1,1.4,4.9,0.9,1.1   ])

    datadic['f_HeII']         =   np.array([10.2,1.4,3.8,2.0,9.3,2.0,10.6,6.9,16.6,2.7,7.0   ])
    datadic['ferr_HeII']      =   np.array([5.1,+99,+99,+99,1.2,1.2,1.9,3.1,3.5,0.4,1.2   ])

    datadic['f_OIII']         =   np.array([7.2,4.6,6.9,4.0,13.0,4.0,7.6,8.6,22.5,4.9,8.2   ])
    datadic['ferr_OIII']      =   np.array([1.4,1.1,1.5,0.8,2.0,1.2,1.1,1.7,3.5,0.9,1.0   ])

    datadic['f_CIII']         =   np.array([8.2,12.4,16.6,4.6,36.9,7.9,22.0,21.4,26.0,11.5,15.2   ])
    datadic['ferr_CIII']      =   np.array([1.4,1.8,3.8,1.9,6.5,1.2,1.5,2.4,5.2,1.3,1.2   ])


    # ---------------------------------------------------------------------------------
    if verbose: print('   Converting fluxes to 1e-20 erg/s/cm2 using fluxscale = '+str(fluxscale))
    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_ber19(fluxscale=1e4,verbose=True):
    """
    Data collected from Berg et al. (2016,2019a,b)

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'ber19'
    # ---------------------------- GENERAL SETUP --------------------------------------
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
    baseid              = lce.referencedictionary()[catreference][0]
    datadic = {}
    datadic['name']      = np.array(['J223831','J141851','J120202','J121402','J084236','J171236','J113116','J133126','J132853','J095430','J132347','J094718','J150934','J100348','J025346','J015809','J104654','J093006','J092055','J082555','J104457','J120122','J124159','J122622','J122436','J124827'])
    datadic['id']        = np.array([223831,141851,120202,121402,84236,171236,113116,133126,132853,95430,132347,94718,150934,100348,25346,15809,104654,93006,92055,82555,104457,120122,124159,122622,122436,124827]) + baseid
    rasex                = np.array(['22:38:31.11','14:18:51.12','12:02:02.49','12:14:02.40','08:42:36.48','17:12:36.72','11:31:16.32','13:31:26.88','13:28:53.04','09:54:30.48','13:23:47.52','09:47:18.24','15:09:34.08','10:03:48.72','02:53:46.70','01:58:09.38','10:46:54.00','09:30:06.48','09:20:55.92','08:25:55.52','10:44:57.79','12:01:22.31','12:41:59.34','12:26:22.71','12:24:36.71','12:48:27.79'])
    decsex               = np.array(['+14:00:28.29','+21:02:39.84','+54:15:51.05','+53:45:17.28','+10:33:14.04','+32:16:33.60','+57:03:58.68','+41:51:48.24','+15:59:34.44','+09:52:12.11','-01:32:51.94','+41:38:16.44','+37:31:46.20','+45:04:57.72','-07:23:43.98','-00:06:37.23','+13:46:45.84','+60:26:53.52','+52:34:07.32','+35:32:31.9','+03:53:13.1','+02:11:08.3','-03:40:02.4','-01:15:12.2','+37:24:36.5','+48:23:03.3'])
    datadic['ra']        = acoord.Angle(rasex, u.hour).degree
    datadic['dec']       = acoord.Angle(decsex, u.degree).degree
    datadic['redshift']  = np.array(['0.021','0.009','0.012','0.003','0.010','0.012','0.006','0.012','0.023','0.005','0.022','0.005','0.033','0.009','0.004','0.012','0.011','0.014','0.008','0.003','0.013','0.003','0.009','0.007','0.040','0.030'])
    datadic['reference'] = [catreference]*len(datadic['id'])
    # datadic['reference'] = ['ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber19','ber16','ber16','ber16','ber16','ber16','ber16','ber16']
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')
    # ---------------------------------------------------------------------------------

    datadic['EW0_CIV1']     = np.array([ 2.00,1.11,np.nan,1.89,3.68,1.93,3.03,np.nan,np.nan,0.76,3.37,np.nan,1.20,0.93,np.nan,np.nan,np.nan,np.nan,np.nan,1.06,6.17,np.nan,np.nan,np.nan,np.nan,1.51])
    datadic['EW0_CIV2']     = np.array([ 2.13,1.41,0.56,0.98,2.84,1.99,0.63,np.nan,np.nan,2.17,1.75,np.nan,0.84,0.86,np.nan,np.nan,1.44,np.nan,0.41,0.76,4.52,np.nan,np.nan,np.nan,np.nan,0.74])
    datadic['EW0_HeII']     = np.array([ np.nan,3.36,0.56,np.nan,np.nan,1.49,1.05,1.24,np.nan,1.28,2.02,np.nan,0.92,np.nan,np.nan,np.nan,np.nan,np.nan,0.78,0.78,2.36,np.nan,np.nan,np.nan,np.nan,1.97])
    datadic['EW0_OIII1']    = np.array([ 1.46,2.26,0.45,1.00,0.73,1.78,np.nan,1.30,0.64,1.50,2.39,0.33,np.nan,np.nan,np.nan,np.nan,np.nan,0.36,np.nan,1.16,2.98,np.nan,np.nan,np.nan,0.92,np.nan])
    datadic['EW0_OIII2']    = np.array([ 2.65,5.31,2.80,1.87,3.10,2.91,1.79,2.55,1.79,2.29,5.71,1.49,2.41,1.54,0.62,1.54,2.07,np.nan,1.19,1.74,5.25,1.72,3.54,2.55,1.72,1.60])
    datadic['EW0_SiIII1']   = np.array([ 3.58,3.09,2.12,3.27,np.nan,6.02,np.nan,2.94,np.nan,3.61,4.13,3.84,2.01,1.59,1.94,3.03,1.97,np.nan,3.14,2.98,3.12,np.nan,2.53,np.nan,np.nan,np.nan])
    datadic['EW0_SiIII2']   = np.array([ 1.10,2.22,1.06,np.nan,np.nan,2.73,np.nan,np.nan,np.nan,3.52,2.39,2.79,1.06,np.nan,np.nan,2.10,np.nan,np.nan,0.67,3.26,2.73,np.nan,3.13,np.nan,1.81,np.nan])
    datadic['EW0_CIII1']    = np.array([ 9.31,10.95,6.37,6.48,7.51,8.65,3.10,6.10,3.94,10.31,5.72,7.16,3.64,4.70,3.66,6.27,7.10,1.73,4.33,7.15,11.70,7.82,6.56,5.52,4.72,5.02])
    datadic['EW0_CIII2']    = np.array([ 5.86,7.46,5.63,10.13,2.28,7.43,3.30,7.02,2.55,5.82,7.00,13.27,6.74,6.24,3.16,7.69,5.04,2.85,5.51,9.34,4.65,4.10,4.00,2.60,4.12,2.49])

    datadic['EW0err_CIV1']     = np.array([np.nan]*26)
    datadic['EW0err_CIV2']     = np.array([np.nan]*26)
    datadic['EW0err_HeII']     = np.array([np.nan]*26)
    datadic['EW0err_OIII1']    = np.array([np.nan]*26)
    datadic['EW0err_OIII2']    = np.array([np.nan]*26)
    datadic['EW0err_SiIII1']   = np.array([np.nan]*26)
    datadic['EW0err_SiIII2']   = np.array([np.nan]*26)
    datadic['EW0err_CIII1']    = np.array([np.nan]*26)
    datadic['EW0err_CIII2']    = np.array([np.nan]*26)

    ### Line fluxes from Berg+19 (scaled by CIII flux) ###
    # Ion  ,  J223831  ,  J141851  ,  J120202  ,  J121402  ,  J084236  ,  J171236  ,  J113116  ,  J133126  ,  J132853  ,  J095430  ,  J132347  ,  J094718  ,  J150934  ,  J100348  ,  J025346  ,  J015809  ,  J104654  ,  J093006  ,  J092055  ,  J084956  ,
    f_CIII = np.array([50.9,80.3,47.2,74.4,23.2,42.6,20.4,100.9,23.9,39.9,59.9,55.0,63.8,45.5,46.5,26.8,54.5,58.9,44.5]) # ,44.3 # 1e-16 erg/s/cm2

    f_CIV1_b19        = np.array([  21.6   ,  19.6   , np.nan  ,  18.5   ,  53.6  ,  21.5    ,  62.2   , np.nan , np.nan ,  7.3    ,  49.0   ,  20.1    ,  17.3   ,  10.8  ,  10.7   ,  20.0   , np.nan  , np.nan ,  7.3   ]) *f_CIII / 100.0 # , np.nan ])
    f_CIV2_b19        = np.array([  22.8   ,  13.0   ,  8.1    ,  9.5    ,  41.4  ,  22.1    ,  13.0   , np.nan , np.nan ,  21.1   ,  25.2   ,  10.2    ,  12.1   ,  9.9   ,  11.0   ,  31.2   ,  18.9   ,  21.7  ,  5.3   ]) *f_CIII / 100.0 # , np.nan ])
    f_HeII_b19        = np.array([ np.nan  ,  27.9   ,  7.3    , np.nan  , np.nan ,  15.0    ,  19.2   ,  13.4  , np.nan ,  11.2   ,  24.7   , np.nan   ,  12.9   , np.nan , np.nan  , np.nan  , np.nan  , np.nan ,  12.2  ]) *f_CIII / 100.0 # , np.nan ])
    f_OIII1_b19       = np.array([  13.7   ,  18.6   ,  5.6    ,  8.5    ,  11.2  ,  17.7    , np.nan  ,  13.7  ,  20.0  ,  13.0   ,  28.4   ,  3.0     , np.nan  , np.nan , np.nan  , np.nan  , np.nan  ,  10.6  , np.nan ]) *f_CIII / 100.0 # , np.nan ])
    f_OIII2_b19       = np.array([  25.3   ,  43.5   ,  34.2   ,  15.8   ,  47.7  ,  29.0    ,  33.1   ,  26.7  ,  52.8  ,  19.7   ,  67.7   ,  11.7    ,  34.0   ,  17.9  ,  16.1   ,  17.5   ,  24.8   , np.nan ,  18.2  ]) *f_CIII / 100.0 # ,  24.6  ])
    f_SiIII1_b19      = np.array([  23.8   ,  17.7   ,  18.1   ,  16.6   , np.nan ,  36.7    , np.nan  ,  21.6  , np.nan ,  22.5   ,  30.4   ,  16.8    ,  19.1   ,  12.9  ,  31.2   ,  20.2   ,  16.1   ,  13.0  ,  33.9  ]) *f_CIII / 100.0 # , np.nan ])
    f_SiIII2_b19      = np.array([  7.3    ,  12.5   ,  8.9    , np.nan  , np.nan ,  16.7    , np.nan  , np.nan , np.nan ,  21.8   ,  18.1   ,  12.8    ,  10.2   , np.nan , np.nan  ,  14.5   , np.nan  ,  18.7  ,  7.4   ]) *f_CIII / 100.0 # ,  21.9  ])
    f_CIII1_b19       = np.array([  61.4   ,  59.6   ,  53.1   ,  38.6   ,  76.9  ,  53.5    ,  48.3   ,  46.4  ,  56.9  ,  63.9   ,  44.8   ,  34.8    ,  35.0   ,  42.8  ,  53.8   ,  44.6   ,  58.4   ,  38.4  ,  50.0  ]) *f_CIII / 100.0 # ,  63.9  ])
    f_CIII2_b19       = np.array([  38.6   ,  40.4   ,  46.9   ,  61.4   ,  23.1  ,  46.5    ,  51.7   ,  53.6  ,  43.1  ,  36.1   ,  55.2   ,  65.2    ,  65.0   ,  57.2  ,  46.2   ,  55.4   ,  41.6   ,  61.6  ,  50.0  ]) *f_CIII / 100.0 # ,  36.1  ])
    f_NII_b19         = np.array([ np.nan  , np.nan  ,  9.8    , np.nan  , np.nan , np.nan   , np.nan  ,  4.4   , np.nan , np.nan  , np.nan  , np.nan   , np.nan  , np.nan , np.nan  , np.nan  , np.nan  , np.nan , np.nan ]) *f_CIII / 100.0 # , np.nan ])

    ferr_CIV1_b19     = np.array([    5.8    ,  3.4    , np.nan ,  4.0   ,  14.6    ,  6.6     ,  13.7  , np.nan , np.nan  ,  4.3    ,  8.8   ,  4.1    ,  4.8     ,  5.0    ,  4.8    ,  9.4   , np.nan , np.nan  ,  2.7  ]) *f_CIII / 100.0 # , np.nan])
    ferr_CIV2_b19     = np.array([    5.9    ,  3.0    ,  5.1   ,  3.7   ,  11.7    ,  6.7     ,  9.5   , np.nan , np.nan  ,  4.8    ,  6.4   ,  3.9    ,  4.5     ,  4.9    ,  4.8    ,  10.3  ,  5.9   ,  7.3    ,  2.7  ]) *f_CIII / 100.0 # , np.nan])
    ferr_HeII_b19     = np.array([   np.nan  ,  4.2    ,  3.6   , np.nan , np.nan   ,  5.8     ,  12.2  ,  4.1   , np.nan  ,  4.9    ,  5.4   , np.nan  ,  5.4     , np.nan  , np.nan  , np.nan , np.nan , np.nan  ,  5.5  ]) *f_CIII / 100.0 # , np.nan])
    ferr_OIII1_b19    = np.array([    4.9    ,  3.5    ,  6.4   ,  3.7   ,  6.1     ,  6.0     , np.nan ,  4.1   ,  9.8    ,  4.9    ,  5.7   ,  3.1    , np.nan   , np.nan  , np.nan  , np.nan , np.nan ,  5.2    , np.nan]) *f_CIII / 100.0  # , np.nan])
    ferr_OIII2_b19    = np.array([    5.5    ,  5.6    ,  7.3   ,  3.9   ,  13.4    ,  7.2     ,  12.2  ,  4.5   ,  14.9   ,  5.2    ,  10.5  ,  3.2    ,  6.9     ,  5.8    ,  4.1    ,  6.5   ,  5.9   , np.nan  ,  5.5  ]) *f_CIII / 100.0 # ,  25.7 ])
    ferr_SiIII1_b19   = np.array([    8.3    ,  4.9    ,  7.6   ,  6.9   , np.nan   ,  13.4    , np.nan ,  6.2   , np.nan  ,  7.8    ,  11.0  ,  5.3    ,  11.4    ,  17.1   ,  9.7    ,  13.0  ,  7.9   ,  8.3    ,  7.0  ]) *f_CIII / 100.0 # , np.nan])
    ferr_SiIII2_b19   = np.array([    7.9    ,  4.7    ,  7.4   , np.nan , np.nan   ,  12.2    , np.nan , np.nan , np.nan  ,  7.8    ,  10.4  ,  5.2    ,  11.1    , np.nan  , np.nan  ,  12.8  , np.nan ,  8.4    ,  6.3  ]) *f_CIII / 100.0 # ,  42.0 ])
    ferr_CIII1_b19    = np.array([    10.5   ,  10.3   ,  9.2   ,  7.7   ,  26.8    ,  15.0    ,  13.8  ,  7.1   ,  20.8   ,  10.0   ,  11.9  ,  5.8    ,  12.2    ,  19.7   ,  11.2   ,  14.8  ,  9.9   ,  9.3    ,  7.7  ]) *f_CIII / 100.0 # ,  54.0 ])
    ferr_CIII2_b19    = np.array([    9.0    ,  9.0    ,  8.8   ,  9.0   ,  19.1    ,  14.3    ,  14.1  ,  7.5   ,  19.0   ,  8.3    ,  12.8  ,  7.0    ,  14.9    ,  21.7   ,  10.7   ,  15.9  ,  8.9   ,  10.8   ,  7.7  ]) *f_CIII / 100.0 # ,  44.9 ])
    ferr_NII_b19      = np.array([   np.nan  , np.nan  ,  22.3  , np.nan , np.nan   , np.nan   , np.nan ,  3.8   , np.nan  , np.nan  , np.nan , np.nan  , np.nan   , np.nan  , np.nan  , np.nan , np.nan , np.nan  , np.nan]) *f_CIII / 100.0  # , np.nan])

    ### Line fluxes from Berg+16 (scaled by Hbeta flux) ###
    #Ion 'J082555','J104457','J120122','J124159','J122622','J122436','J124827']
    f_Hbeta = np.array([   230.8,413.7,114.3,98.3,8131.0,138.4,78.1 ])

    f_CIV1_b16       = np.array([ 0.55 ,  2.07  ,  np.nan ,  np.nan ,  np.nan  ,  np.nan ,  0.82 ]) *f_Hbeta
    f_CIV2_b16       = np.array([ 0.42 ,  1.52  ,  np.nan ,  np.nan ,  np.nan  ,  np.nan ,  0.40 ]) *f_Hbeta
    f_HeII_b16       = np.array([ 0.38 ,  0.70  ,  np.nan ,  np.nan ,  np.nan  ,  np.nan ,  0.96 ]) *f_Hbeta
    f_OIII1_b16      = np.array([ 0.49 ,  0.85  ,  0.57   ,   0.20   , np.nan  ,  0.53   , np.nan]) *f_Hbeta
    f_OIII2_b16      = np.array([ 0.83 ,  1.49  ,  1.05   ,   1.72   ,  0.010  ,  1.28   ,  0.77 ]) *f_Hbeta
    f_SiIII1_b16     = np.array([ 1.16 ,  0.75  ,  np.nan ,   1.10   ,  np.nan ,  np.nan , np.nan]) *f_Hbeta
    f_SiIII2_b16     = np.array([ 1.25 ,  0.66  ,  np.nan ,   1.38   ,  np.nan ,  1.29   , np.nan]) *f_Hbeta
    f_CIII1_b16      = np.array([ 2.63 ,  2.83  ,  3.89   ,   2.94   ,  0.022  ,  3.38   ,  2.12 ]) *f_Hbeta
    f_CIII2_b16      = np.array([ 3.44 ,  1.12  ,  2.04   ,   1.80   ,  0.012  ,  2.96   ,  1.05 ]) *f_Hbeta

    ferr_CIV1_b16    = np.array([ 0.02 ,  0.07  ,  np.nan ,  np.nan  ,  np.nan , np.nan ,   0.24]) *f_Hbeta
    ferr_CIV2_b16    = np.array([ 0.01 ,  0.06  ,  np.nan ,  np.nan  ,  np.nan , np.nan ,   0.24]) *f_Hbeta
    ferr_HeII_b16    = np.array([ 0.01 ,  0.05  ,  np.nan ,  np.nan  ,  np.nan , np.nan ,   0.24]) *f_Hbeta
    ferr_OIII1_b16   = np.array([ 0.02 ,  0.05  ,   0.20  ,   0.20   ,  np.nan , 0.11   , np.nan]) *f_Hbeta
    ferr_OIII2_b16   = np.array([ 0.02 ,  0.06  ,   0.20  ,   0.28   ,  0.001  , 0.11   ,   0.22]) *f_Hbeta
    ferr_SiIII1_b16  = np.array([ 0.03 ,  0.06  ,  np.nan ,   0.25   ,  np.nan , np.nan , np.nan]) *f_Hbeta
    ferr_SiIII2_b16  = np.array([ 0.04 ,  0.05  ,  np.nan ,   0.25   ,  np.nan , 0.04   , np.nan]) *f_Hbeta
    ferr_CIII1_b16   = np.array([ 0.08 ,  0.08  ,   0.22  ,   0.32   ,  0.001  , 0.13   ,   0.18]) *f_Hbeta
    ferr_CIII2_b16   = np.array([ 0.10 ,  0.06  ,   0.20  ,   0.27   ,  0.001  , 0.12   ,   0.16]) *f_Hbeta


    datadic['f_CIV1']        = np.append(f_CIV1_b19     ,f_CIV1_b16      )
    datadic['f_CIV2']        = np.append(f_CIV2_b19     ,f_CIV2_b16      )
    datadic['f_HeII']        = np.append(f_HeII_b19     ,f_HeII_b16      )
    datadic['f_OIII1']       = np.append(f_OIII1_b19    ,f_OIII1_b16     )
    datadic['f_OIII2']       = np.append(f_OIII2_b19    ,f_OIII2_b16     )
    datadic['f_SiIII1']      = np.append(f_SiIII1_b19   ,f_SiIII1_b16    )
    datadic['f_SiIII2']      = np.append(f_SiIII2_b19   ,f_SiIII2_b16    )
    datadic['f_CIII1']       = np.append(f_CIII1_b19    ,f_CIII1_b16     )
    datadic['f_CIII2']       = np.append(f_CIII2_b19    ,f_CIII2_b16     )
    #datadic['f_NII']        =           f_NII_b19

    datadic['ferr_CIV1']     = np.append(ferr_CIV1_b19  ,ferr_CIV1_b16   )
    datadic['ferr_CIV2']     = np.append(ferr_CIV2_b19  ,ferr_CIV2_b16   )
    datadic['ferr_HeII']     = np.append(ferr_HeII_b19  ,ferr_HeII_b16   )
    datadic['ferr_OIII1']    = np.append(ferr_OIII1_b19 ,ferr_OIII1_b16  )
    datadic['ferr_OIII2']    = np.append(ferr_OIII2_b19 ,ferr_OIII2_b16  )
    datadic['ferr_SiIII1']   = np.append(ferr_SiIII1_b19,ferr_SiIII1_b16 )
    datadic['ferr_SiIII2']   = np.append(ferr_SiIII2_b19,ferr_SiIII2_b16 )
    datadic['ferr_CIII1']    = np.append(ferr_CIII1_b19 ,ferr_CIII1_b16  )
    datadic['ferr_CIII2']    = np.append(ferr_CIII2_b19 ,ferr_CIII2_b16  )
    #datadic['ferr_NII']     =           ferr_NII_b19

    linename = 'CIV'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    linename = 'OIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    linename = 'SiIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

    linename = 'CIII'
    datadic['f_'+linename], datadic['ferr_'+linename], \
    datadic['FR_'+linename+'1'+linename+'2'], datadic['FRerr_'+linename+'1'+linename+'2'], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])


    # ---------------------------------------------------------------------------------
    if verbose: print('   Converting fluxes to 1e-20 erg/s/cm2 using fluxscale = '+str(fluxscale))
    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    dataarray = lce.build_dataarray(catreference, datadic, S2Nlim=3.0,verbose=False)
    if verbose: print('   Returning catalog reference and data array')
    return catreference, dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
