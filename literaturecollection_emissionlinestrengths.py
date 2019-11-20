# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Collection of literature emission line measurements combined into single table
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import sys
import os
import astropy.io.fits as afits
from fits2ascii import ascii2fits
import MiGs
import numpy as np
import collections
import literaturecollection_emissionlinestrengths as lce
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def emissionlinelist(verbose=True):
    """

    """
    linelistdic = MiGs.linelistdic(listversion='full')
    if verbose: print(' - Loaded the line dictionary from MiGs containing:\n   '+linelistdic.keys())
    return linelistdic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_reference(idlist,verbose=True):
    """
    Function to return the reference(s) given a list of literature IDs

    --- INPUT ---
    idlist       List of "literature ids". They should be 11 digits long where the first 3 digits refers
                 to the reference the IDs belong to as provided by lce.referencedictionary()
    verbose      Toggle verbosity of function

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce

    idlist     = [1000000022,1004567890,2000000888,3000000888]
    references = lce.get_reference(idlist)


    """
    refdic     = lce.referencedictionary()
    returnrefs = []

    if type(idlist) is not list:
        idlist = [int(idlist)]

    if verbose: print(' - Looping over the '+str(len(idlist))+' IDs provide in list\n')
    for id in idlist:
        idstr  = str(id)

        if (len(idstr) != 11) or idstr.startswith('0'):
            if verbose: print(' - The id '+idstr+' is not on the right format. It should be 11 digits long and not start with 0')
        else:
            baseid = int(idstr[0:3].ljust(11,'0'))

            foundref = False
            for key in refdic.keys():
                if int(refdic[key][0]) == baseid:
                    returnrefs.append([idstr,key,refdic[key][0],refdic[key][1],refdic[key][2]])
                    foundref = True

            if not foundref:
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

    """
    refdic = collections.OrderedDict()
    #                  baseid   reference                     plotsymbol
    refdic['sen17'] = [1e10,    'Senchyna et al. 2017',       '+']
    refdic['nan19'] = [2e10,    'Nanaykkara et al. 2019',     's']
    refdic['eor']   = [3e10,    'EoR objects',                'v']

    return refdic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def lineinfo(verbose=True):
    """

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce
    linetypedic, linelistdic = lineinfo()

    """
    linetypedic = collections.OrderedDict()
    linetypedic['Lya']   = 'single'
    linetypedic['NV']    = 'doublet'
    linetypedic['CIV']   = 'doublet'
    linetypedic['HeII']  = 'doublet'
    linetypedic['OIII']  = 'doublet'
    linetypedic['SiIII'] = 'doublet'
    linetypedic['CIII']  = 'doublet'
    linetypedic['MgII']  = 'single'

    linelistdic = MiGs.linelistdic(listversion='full')

    return linetypedic, linelistdic
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
        print(' Something went wrong in lce.set_ratios() - setting trace to investigate')
        pdb.set_trace()

    return ratio, ratioerr

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def collect_data(outputfile,verbose=True):
    """
    Collecting literature values and storing them in a single file

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce

    outdir     = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/'
    outputfile = outdir+'literaturecollection_emissionlinestrengths.txt'
    dataarray  = lce.collect_data(outputfile)

    """

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Initializing the output file: \n   '+outputfile)
    fout = open(outputfile,'w')
    fout.write('# Flux, EW and line ratios collected from the literature:\n')
    fout.write('# \n')
    fout.write('# Upper and lower limits are given as values with uncertainty of +99 or -99, respectively. \n')
    fout.write('# \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Determine columns to fill in output')
    fluxratiodic = build_master_datadictionary()

    Ncols      = len(fluxratiodic.keys())
    # a simple dictionary containing the column locations in the output array (indexes)
    colents    = {}
    for oo, colname in enumerate(fluxratiodic.keys()):
        colents[colname] = oo
    if verbose: print('   The output file will contain '+str(Ncols)+' columns ')
    fout.write('# This file contains the following '+str(Ncols)+' columns:\n')
    fout.write('# '+' '.join(fluxratiodic.keys())+'  \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Collecting literature data arrays')

    appenddatadic = collections.OrderedDict()
    #------------------------------------------------
    dataref, appenddata     = lce.data_sen17(fluxscale=1e5)
    appenddatadic[dataref]  = appenddata

    dataref, appenddata     = lce.data_nan19(fluxscale=1.0)
    appenddatadic[dataref]  = appenddata

    dataref, appenddata     = lce.data_eor(fluxscale=1.0)
    appenddatadic[dataref]  = appenddata
    #------------------------------------------------
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Filling the columns of the output file with data ')
    for appendkey in appenddatadic.keys():
        appenddata = appenddatadic[appendkey]
        Nobjappend     = len(appenddata)
        fluxratioarray = np.zeros([Nobjappend,Ncols])*np.nan
        for objent, objappend in enumerate(appenddata):
            for acol in appenddata.dtype.names:
                if (acol in fluxratiodic.keys()) & (acol != 'reference'):
                    fluxratioarray[objent,colents[acol]] = appenddata[acol][objent]
                else:
                    if (acol != 'reference'):
                        if verbose: print(' WARNING: The column "'+str(acol)+'" has no match in output file')

        for ll in np.arange(Nobjappend):
            outstr = str(int(fluxratioarray[ll,0]))+'  '+appendkey+' '+' '.join([str("%10.4f" % ff) for ff in fluxratioarray[ll,2:]])
            fout.write(outstr+' \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Done filling array. Writing the data to \n   '+outputfile)
    fout.close()
    fmt = '20a,20a,'+','.join((Ncols-2)*['d'])
    dataarray = np.genfromtxt(outputfile,skip_header=5,dtype=fmt,comments='#',names=True)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Creating fits version of output: '+outputfile.replace('.txt','.fits'))
    fitsformat = ['D','20A'] + ['D']*(Ncols-2)
    fitsoutput = ascii2fits(outputfile,asciinames=True,skip_header=5,fitsformat=fitsformat,verbose=False)

    return dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_dataarray(catreference, linetypedic, datadic, verbose=True):
    """
    Function building data array from literature data

    --- INPUT ---
    ids             The ids of the objects
    redshifts       The redshifts of the objects
    catreference    Catalog reference (should exist in lce.referencedictionary())
    linetypedic     Dictionary with line-types as keys providing keyword on whether to treat them as doublets or
                    single lines when populating flux ratio columns.
    f_lines         Fluxes of objects for each of the lines. Expects list of arrays, i.e. [f_line1,f_line2,...,f_lineN]
                    where each f_line* is a 1D numpy array with data
    f_errlines      Errors on fluxe of objects for each of the line. Format similar to f_lines
    EW0_lines       EW0 of objects for each of the line. Format similar to f_lines
    EW0err_lines    Error on EW0 of objects for each of the line. Format similar to f_lines
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---
    dataarray = lce.build_dataarray(catreference, linetypedic, datadic) # see example in lce.data_sen17()

    """
    Ndataobj      = len(datadic['id'])

    masterdic     = lce.build_master_datadictionary()
    columns       = masterdic.keys()

    dtypebuild    = [(colname, 'd') for colname in columns]
    dtypebuild[1] = ('reference', '20a')
    dataarray     = np.array(np.zeros(Ndataobj)*np.nan,dtype=dtypebuild)

    if verbose: print('   Inserting data from literature into master dictionary ')
    for ll, mastercol in enumerate(masterdic.keys()):
        if mastercol in datadic.keys():
            dataarray[mastercol] = datadic[mastercol]

    if verbose: print('   Estimating further entries based on data from literature')

    for ii, id in enumerate(dataarray['id']):
        for ll, numerator_line in enumerate(linetypedic.keys()):
            fnum    = dataarray['f_'+numerator_line][ii]
            ferrnum = dataarray['ferr_'+numerator_line][ii]
            dataarray['vshift_'+numerator_line][ii] = 0.0
            dataarray['sigma_'+numerator_line][ii]  = np.nan

            if ferrnum == 99:
                dataarray['s2n_'+numerator_line][ii]     = np.nan
            else:
                dataarray['s2n_'+numerator_line][ii]     = fnum / ferrnum

            for kk, denominator_line in enumerate(linetypedic.keys()):
                fdenom    = dataarray['f_'+denominator_line][ii]
                ferrdenom = dataarray['ferr_'+denominator_line][ii]

                if numerator_line == denominator_line:
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

                    dataarray['FR_'+numerator_line+denominator_line][ii]     = FR
                    dataarray['FRerr_'+numerator_line+denominator_line][ii]  = FRerr
                    if np.abs(FRerr) == 99:
                        dataarray['FRs2n_'+numerator_line+denominator_line][ii]  = np.nan
                    else:
                        dataarray['FRs2n_'+numerator_line+denominator_line][ii]  = FR/FRerr

    return dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_master_datadictionary():
    """

    """
    datadictionary = collections.OrderedDict()
    datadictionary['id']        = np.array([])
    datadictionary['reference'] = np.array([])
    datadictionary['redshift']  = np.array([])

    linetypedic, linelistdic = lineinfo()

    for ll, numerator_line in enumerate(linetypedic.keys()):
        datadictionary['f_'+numerator_line]       = np.array([])
        datadictionary['ferr_'+numerator_line]    = np.array([])
        datadictionary['s2n_'+numerator_line]     = np.array([])
        datadictionary['EW0_'+numerator_line]      = np.array([])
        datadictionary['EW0err_'+numerator_line]   = np.array([])
        datadictionary['sigma_'+numerator_line]   = np.array([])
        datadictionary['vshift_'+numerator_line]  = np.array([])

        if  linetypedic[numerator_line] == 'doublet':
            datadictionary['f_'+numerator_line+'1']      = np.array([])
            datadictionary['ferr_'+numerator_line+'1']   = np.array([])
            datadictionary['f_'+numerator_line+'2']      = np.array([])
            datadictionary['ferr_'+numerator_line+'2']   = np.array([])
            datadictionary['FR_'+numerator_line+'1'+numerator_line+'2']     = np.array([])
            datadictionary['FRerr_'+numerator_line+'1'+numerator_line+'2']  = np.array([])
            datadictionary['FRs2n_'+numerator_line+'1'+numerator_line+'2']  = np.array([])

        for kk, denominator_line in enumerate(linetypedic.keys()):
            if numerator_line == denominator_line:
                continue
            else:
                datadictionary['FR_'+numerator_line+denominator_line]     = np.array([])
                datadictionary['FRerr_'+numerator_line+denominator_line]  = np.array([])
                datadictionary['FRs2n_'+numerator_line+denominator_line]  = np.array([])

                if (linetypedic[denominator_line] == 'single'):
                    datadictionary['FR_'+numerator_line+'1'+denominator_line]     = np.array([])
                    datadictionary['FRerr_'+numerator_line+'1'+denominator_line]  = np.array([])
                    datadictionary['FRs2n_'+numerator_line+'1'+denominator_line]  = np.array([])
                    datadictionary['FR_'+numerator_line+'2'+denominator_line]     = np.array([])
                    datadictionary['FRerr_'+numerator_line+'2'+denominator_line]  = np.array([])
                    datadictionary['FRs2n_'+numerator_line+'2'+denominator_line]  = np.array([])
                if (linetypedic[denominator_line] == 'doublet'):
                    datadictionary['FR_'+numerator_line+denominator_line+'1']     = np.array([])
                    datadictionary['FRerr_'+numerator_line+denominator_line+'1']  = np.array([])
                    datadictionary['FRs2n_'+numerator_line+denominator_line+'1']  = np.array([])
                    datadictionary['FR_'+numerator_line+denominator_line+'2']     = np.array([])
                    datadictionary['FRerr_'+numerator_line+denominator_line+'2']  = np.array([])
                    datadictionary['FRs2n_'+numerator_line+denominator_line+'2']  = np.array([])
                if (linetypedic[numerator_line] == 'doublet') & \
                        (linetypedic[denominator_line] == 'doublet'):
                    datadictionary['FR_'+numerator_line+'1'+denominator_line+'1']     = np.array([])
                    datadictionary['FRerr_'+numerator_line+'1'+denominator_line+'1']  = np.array([])
                    datadictionary['FRs2n_'+numerator_line+'1'+denominator_line+'1']  = np.array([])
                    datadictionary['FR_'+numerator_line+'1'+denominator_line+'2']     = np.array([])
                    datadictionary['FRerr_'+numerator_line+'1'+denominator_line+'2']  = np.array([])
                    datadictionary['FRs2n_'+numerator_line+'1'+denominator_line+'2']  = np.array([])
                    datadictionary['FR_'+numerator_line+'2'+denominator_line+'1']     = np.array([])
                    datadictionary['FRerr_'+numerator_line+'2'+denominator_line+'1']  = np.array([])
                    datadictionary['FRs2n_'+numerator_line+'2'+denominator_line+'1']  = np.array([])
                    datadictionary['FR_'+numerator_line+'2'+denominator_line+'2']     = np.array([])
                    datadictionary['FRerr_'+numerator_line+'2'+denominator_line+'2']  = np.array([])
                    datadictionary['FRs2n_'+numerator_line+'2'+denominator_line+'2']  = np.array([])
    return datadictionary
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_nan19(fluxscale=1.0,verbose=True):
    """
    Data collected from Nanaykkara et al. 19

    --- INPUT ---
    baseid    Number to add to ids to keep them distinct from other literature data sets
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'nan19'
    refdic              = lce.referencedictionary()
    if verbose: print(' - Assembling the data from '+refdic[catreference][1])
    fluxratiodic        = collections.OrderedDict()
    fluxratiodic['id']  = np.array([])
    baseid              = lce.referencedictionary()[catreference][0]

    linetypedic = collections.OrderedDict()
    linetypedic['CIII']  = 'single'
    linetypedic['OIII']  = 'single'
    linetypedic['HeII']  = 'single'

    datadic = {}
    datadic['id']        = np.array([1024, 1036, 1045, 1079, 1273, 3621, 87,   109, 144,    97,   39,  84, 161]) + baseid
    datadic['redshift']  = np.array([2.87, 2.69, 2.61, 2.68, 2.17, 3.07, 2.67, 2.2, 4.02, 2.11, 3.96, 3.1, 3.1])
    datadic['reference'] = [catreference]*len(datadic['id'])

    # # EW - from Anna
    # datadic['EW0_HeII']   = np.array([0.79487179, 0.81081081, 1., 2.5625, 2.51219512, 1.2972973, 1., 2.14,
    #                                    5.77419355, 4.94, 6.56097561, 12.26829268])
    # datadic['EW0err_HeII'] = np.array([0.15384615, 0.2972973, 0.25, 0.625, +99.0, 0.21621622, 0.125, 0.44,
    #                                    0.77419355, +99.0, 0.85365854, +99.0])
    # datadic['EW0_OIII']   = np.array([0.71794872, 1.2972973, 1.27777778, 2.25, 2.46341463, 1.02702703,
    #                                    1.34375, 7.44, 5.4516129, 0.26, 1.53658537, 1.68292683])
    # datadic['EW0err_OIII'] = np.array([0.15384615, 0.18918919, 0.22222222, 0.4375, +99.0, 0.13513514,
    #                                     0.15625, 0.86, 0.61290323, +99.0, 0.48780488, +99.0])
    # datadic['EW0_CIII']    = np.array([2.56410256, 5.64864865, 4.94444444, 11.15625, 4.7804878,
    #                                     2.72972973, 3.34375, np.nan, 15.83870968, np.nan, 7.46341463, 2.6097561])
    # datadic['EW0err_CIII'] = np.array([0.26399052, 0.34611482, 0.45218946, 0.73221091, +99.0,
    #                                     0.30577591, 0.22097087, +99.0,  1.58196127, +99.0,  +99.0, +99.0])
    # datadic['EW0_CIII']    = np.array([1.5384615385, 3.3513, 3.1944, 6.65625, 1.756, 1.9189, 1.78125,
    #                                     np.nan,9.32258, np.nan, 4.53659, 2.5366])
    # datadic['EW0err_CIII'] = np.array([0.23,0.2162,  0.33, 0.5625, +99.0, 0.2162, 0.15625, +99.0, 1.,
    #                                     +99.0, +99.0, +99.0 ])
    # line intensities - from Anna
    datadic['f_HeII']    = np.array([177., 142., 156., 290., 217., 213., 59., 54., 48., 306., 153., 161., 318.])
    datadic['ferr_HeII'] = np.array([52., 53., 58., 91., 79., 45., 12., 13., 12., 55., 33., 29., 39.])
    datadic['f_OIII']    = np.array([161., 230., 200., 81., 195., 105., 49., 72., 150., 284., 66., 54., 58])
    datadic['ferr_OIII'] = np.array([48., 54., 56., +99.0, 56., +99.0, 11., 13., 21., 44., +99.0, +99.0, 19.])
    datadic['f_CIII']    = np.array([514., 736., 599., 222., 673., 464., 111., 134., np.nan, 626., np.nan, 246., 205.])
    datadic['ferr_CIII'] = np.array([64.54,54.56, 72.92, +99.0, 67.18, +99.0,  15.55, 16.97, +99.0, 79.65, +99.0, +99.0, +99.0])

    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    dataarray = lce.build_dataarray(catreference, linetypedic, datadic)

    return catreference, dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_eor(fluxscale=1.0,verbose=True):
    """

    --- INPUT ---
    baseid    Number to add to ids to keep them distinct from other literature data sets
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'eor'
    refdic              = lce.referencedictionary()
    if verbose: print(' - Assembling the data from '+refdic[catreference][1])
    fluxratiodic        = collections.OrderedDict()
    fluxratiodic['id']  = np.array([])
    baseid              = lce.referencedictionary()[catreference][0]

    linetypedic = collections.OrderedDict()
    linetypedic['CIV']   = 'single'
    linetypedic['CIII']  = 'single'
    linetypedic['OIII']  = 'single'
    linetypedic['HeII']  = 'single'

    datadic = {}
    datadic['id']          = np.array([2248, 2, 3]) + baseid
    datadic['redshift']    = np.array([6.1, 7.0, 7.1])
    datadic['reference']   = [catreference]*len(datadic['id'])

    # line intensities
    datadic['f_HeII']    = np.array([150,  210,    126])
    datadic['ferr_HeII'] = np.array([+99,  +99,    29])
    datadic['f_OIII']    = np.array([440,  180,    np.nan])
    datadic['ferr_OIII'] = np.array([60,   70,     99])
    datadic['f_CIII']    = np.array([500,  np.nan, 175])
    datadic['ferr_CIII'] = np.array([+99,  +99,    +99])
    datadic['f_CIV']     = np.array([1140, 790,    270])
    datadic['ferr_CIV']  = np.array([90,   90,     +99])

    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    dataarray = lce.build_dataarray(catreference, linetypedic, datadic)

    return catreference, dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_sen17(fluxscale=1e5,verbose=True):
    """
    Data collected from Senchyna et al. 17

    --- INPUT ---
    baseid      Number to add to ids to keep them distinct from other literature data sets
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce
    catreference, dataarray = lce.data_sen17()

    """
    catreference        = 'sen17'
    refdic              = lce.referencedictionary()
    if verbose: print(' - Assembling the data from '+refdic[catreference][1])
    fluxratiodic        = collections.OrderedDict()
    fluxratiodic['id']  = np.array([])
    baseid              = refdic[catreference][0]

    linetypedic = collections.OrderedDict()
    linetypedic['CIII']  = 'single'
    linetypedic['OIII']  = 'single'
    linetypedic['HeII']  = 'single'

    datadic = {}
    datadic['id']          = np.array([2, 36, 80, 82, 110, 111, 179, 182, 191, 198]) + baseid
    datadic['redshift']    = np.array([0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003]) # nearby dwarfs
    datadic['reference']   = [catreference]*len(datadic['id'])

    #XXX CIV exists...

    # EW - from Anna
    datadic['EW0_HeII']    = np.array([1.7, 0.4, 0.7, 0.44, 0.4, 1.42, 0.6, 0.94, 0.5, 0.5])
    datadic['EW0err_HeII'] = np.array([0.06, +99.0, +99.0, 0.02, +99.0, 0.03, +99.0, 0.04, +99.0, +99.0])
    datadic['EW0_OIII']    = np.array([5.05, 0.84, 1.68, 1.89, 0.48, 1.73, 1.17, 1.73, np.nan, 0.9])
    datadic['EW0err_OIII'] = np.array([0.1, 0.03, 0.08, 0.02, 0.02, 0.05, 0.06, 0.04, np.nan, +99.0])
    datadic['EW0_CIII']    = np.array([14.86, 4.98, 4.0, 12.09, 2.8, 3.3, 8.71, 13.35, 11.33, 3.38])
    datadic['EW0err_CIII'] = np.array([1.07, 0.59, +99.0, 0.3, +99.0, +99.0, 0.42, 0.52, 0.34, 0.31])

    # line intensities - from Anna
    datadic['f_HeII']    = np.array([1.45, 1.1, 0.8, 1.22, 1.2,2.37, 0.8, 1.42, 1.0, 0.8])
    datadic['ferr_HeII'] = np.array([0.06, +99.0, +99.0, 0.05, +99.0, 0.05, +99.0, 0.05, +99.0, +99.0])
    datadic['f_OIII']    = np.array([3.35, 2.16, 1.69, 5.05, 1.28, 2.68, 1.36, 2.27, np.nan, 1.1])
    datadic['ferr_OIII'] = np.array([0.04, 0.08, 0.07, 0.05, 0.05, 0.07, 0.06, 0.05, np.nan, +99.0])
    datadic['f_CIII']    = np.array([11.06, 7.60, 4.1, 22.5, 5.2, 3.9, 6.80, 12.3, 15.9, 3.79])
    datadic['ferr_CIII'] = np.array([0.63, 0.88, +99.0, 0.40, +99.0, +99.0, 0.23, 0.33, 0.38, 0.30])

    for key in datadic.keys():
        if key.startswith('f'):
            datadic[key][np.abs(datadic[key]) != 99] = datadic[key][np.abs(datadic[key]) != 99]*fluxscale

    dataarray = lce.build_dataarray(catreference, linetypedic, datadic)

    return catreference, dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
