# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Collection of literature emission line measurements combined into single table
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import sys
import os
import astropy.io.fits as afits
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
def lineinfo(lines,verbose=True):
    """

    """
    linelistdic = MiGs.linelistdic(listversion='full')

    return linelistdic

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
def restUVline_literaturecollection(outputfile,verbose=True):
    """
    Collecting literature values and storing them in single file

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce

    outputfile = '/Users/kschmidt/Desktop/photplots/literature_UVlinemeasurements.txt'
    dataarray  = lce.restUVline_literaturecollection(outputfile)

    """

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Initializing the output file: \n   '+outputfile)
    fout = open(outputfile,'w')
    fout.write('# Flux, EW and line ratios collected from the literature:\n')
    fout.write('# \n')
    fout.write('# Upper and lower limits are given as values with uncertainty of +99 or -99, respectively. \n')
    fout.write('# \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Determine columns to fill in output')
    fluxratiodic = collections.OrderedDict()
    fluxratiodic['id']        = np.array([])
    fluxratiodic['reference'] = np.array([])
    fluxratiodic['redshift']  = np.array([])

    linetypedic = collections.OrderedDict()
    linetypedic['CIII']  = 'doublet'
    linetypedic['SiIII'] = 'doublet'
    linetypedic['CIV']   = 'doublet'
    linetypedic['MgII']  = 'single'
    linetypedic['OIII']  = 'doublet'
    linetypedic['NV']    = 'doublet'
    linetypedic['HeII']  = 'doublet'

    for ll, numerator_line in enumerate(linetypedic.keys()):
        fluxratiodic['f_'+numerator_line]       = np.array([])
        fluxratiodic['ferr_'+numerator_line]    = np.array([])
        fluxratiodic['s2n_'+numerator_line]     = np.array([])
        fluxratiodic['ew_'+numerator_line]      = np.array([])
        fluxratiodic['ewerr_'+numerator_line]   = np.array([])
        fluxratiodic['sigma_'+numerator_line]   = np.array([])
        fluxratiodic['vshift_'+numerator_line]  = np.array([])

        if  linetypedic[numerator_line] == 'doublet':
            fluxratiodic['f_'+numerator_line+'1']      = np.array([])
            fluxratiodic['ferr_'+numerator_line+'1']   = np.array([])
            fluxratiodic['f_'+numerator_line+'2']      = np.array([])
            fluxratiodic['ferr_'+numerator_line+'2']   = np.array([])
            fluxratiodic['FR_'+numerator_line+'1'+numerator_line+'2']     = np.array([])
            fluxratiodic['FRerr_'+numerator_line+'1'+numerator_line+'2']  = np.array([])
            fluxratiodic['FRs2n_'+numerator_line+'1'+numerator_line+'2']  = np.array([])

        for kk, denominator_line in enumerate(linetypedic.keys()):
            if numerator_line == denominator_line:
                continue
            else:
                fluxratiodic['FR_'+numerator_line+denominator_line]     = np.array([])
                fluxratiodic['FRerr_'+numerator_line+denominator_line]  = np.array([])
                fluxratiodic['FRs2n_'+numerator_line+denominator_line]  = np.array([])

                if (linetypedic[denominator_line] == 'single'):
                    fluxratiodic['FR_'+numerator_line+'1'+denominator_line]     = np.array([])
                    fluxratiodic['FRerr_'+numerator_line+'1'+denominator_line]  = np.array([])
                    fluxratiodic['FRs2n_'+numerator_line+'1'+denominator_line]  = np.array([])
                    fluxratiodic['FR_'+numerator_line+'2'+denominator_line]     = np.array([])
                    fluxratiodic['FRerr_'+numerator_line+'2'+denominator_line]  = np.array([])
                    fluxratiodic['FRs2n_'+numerator_line+'2'+denominator_line]  = np.array([])
                if (linetypedic[denominator_line] == 'doublet'):
                    fluxratiodic['FR_'+numerator_line+denominator_line+'1']     = np.array([])
                    fluxratiodic['FRerr_'+numerator_line+denominator_line+'1']  = np.array([])
                    fluxratiodic['FRs2n_'+numerator_line+denominator_line+'1']  = np.array([])
                    fluxratiodic['FR_'+numerator_line+denominator_line+'2']     = np.array([])
                    fluxratiodic['FRerr_'+numerator_line+denominator_line+'2']  = np.array([])
                    fluxratiodic['FRs2n_'+numerator_line+denominator_line+'2']  = np.array([])
                if (linetypedic[numerator_line] == 'doublet') & \
                        (linetypedic[denominator_line] == 'doublet'):
                    fluxratiodic['FR_'+numerator_line+'1'+denominator_line+'1']     = np.array([])
                    fluxratiodic['FRerr_'+numerator_line+'1'+denominator_line+'1']  = np.array([])
                    fluxratiodic['FRs2n_'+numerator_line+'1'+denominator_line+'1']  = np.array([])
                    fluxratiodic['FR_'+numerator_line+'1'+denominator_line+'2']     = np.array([])
                    fluxratiodic['FRerr_'+numerator_line+'1'+denominator_line+'2']  = np.array([])
                    fluxratiodic['FRs2n_'+numerator_line+'1'+denominator_line+'2']  = np.array([])
                    fluxratiodic['FR_'+numerator_line+'2'+denominator_line+'1']     = np.array([])
                    fluxratiodic['FRerr_'+numerator_line+'2'+denominator_line+'1']  = np.array([])
                    fluxratiodic['FRs2n_'+numerator_line+'2'+denominator_line+'1']  = np.array([])
                    fluxratiodic['FR_'+numerator_line+'2'+denominator_line+'2']     = np.array([])
                    fluxratiodic['FRerr_'+numerator_line+'2'+denominator_line+'2']  = np.array([])
                    fluxratiodic['FRs2n_'+numerator_line+'2'+denominator_line+'2']  = np.array([])

    Ncols      = len(fluxratiodic.keys())
    # a simple dictionary containing the column locations in the output array (indexes)
    colents    = {}
    for oo, colname in enumerate(fluxratiodic.keys()):
        colents[colname] = oo
    if verbose: print('   The output file will contain '+str(Ncols)+' columns ')
    fout.write('# This file contains the following '+str(Ncols)+' columns:\n')
    fout.write('# '+' '.join(fluxratiodic.keys())+'  \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Filling the columns with data ')

    appenddatadic = {}
    #------------------------------------------------
    dataref, appenddata     = lce.restUVline_literaturecollection_Senchyna17(baseid=1e10)
    appenddatadic[dataref]  = appenddata

    dataref, appenddata     = lce.restUVline_literaturecollection_Nanaykkara19(baseid=2e10)
    appenddatadic[dataref]  = appenddata

    dataref, appenddata     = lce.restUVline_literaturecollection_EoR(baseid=3e10)
    appenddatadic[dataref]  = appenddata
    #------------------------------------------------

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
                        if verbose: print(' WARNING: The Senchyna17 column "'+str(acol)+'" has no match in output file')

        for ll in np.arange(Nobjappend):
            outstr = str(int(fluxratioarray[ll,0]))+'  '+appendkey+' '+' '.join([str("%10.4f" % ff) for ff in fluxratioarray[ll,2:]])
            fout.write(outstr+' \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Done filleing array. Writing the flux ratio output to \n   '+outputfile)
    fout.close()
    fmt = '20a,20a,'+','.join((Ncols-2)*['d'])
    dataarray = np.genfromtxt(outputfile,skip_header=5,dtype=fmt,comments='#',names=True)
    return dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def restUVline_literaturecollection_Senchyna17(fluxscale=1e5,verbose=True):
    """

    --- INPUT ---
    baseid    Number to add to ids to keep them distinct from other literature data sets
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'sen17'
    fluxratiodic        = collections.OrderedDict()
    fluxratiodic['id']  = np.array([])
    baseid              = lce.referencedictionary[catreference][0]

    ids        = np.array([2, 36, 80, 82, 110, 111, 179, 182, 191, 198]) + baseid
    redshifts  = np.array([0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003]) # nearby dwarfs

    dataarray = np.array(np.zeros(len(ids))*np.nan,
                         dtype=[('id', 'd'), ('reference', '20a'), ('redshift', 'd'),
                                ('f_HeII', 'd'), ('ferr_HeII', 'd'), ('s2n_HeII', 'd'), ('vshift_HeII', 'd'), ('sigma_HeII', 'd'),
                                ('ew_HeII', 'd'), ('ewerr_HeII', 'd'),
                                ('f_OIII', 'd'), ('ferr_OIII', 'd'), ('s2n_OIII', 'd'), ('vshift_OIII', 'd'), ('sigma_OIII', 'd'),
                                ('ew_OIII', 'd'), ('ewerr_OIII', 'd'),
                                ('f_CIII', 'd'), ('ferr_CIII', 'd'), ('s2n_CIII', 'd'), ('vshift_CIII', 'd'), ('sigma_CIII', 'd'),
                                ('ew_CIII', 'd'), ('ewerr_CIII', 'd'),
                                ('FR_HeIICIII', 'd'), ('FRerr_HeIICIII', 'd'), ('FRs2n_HeIICIII', 'd'),
                                ('FR_HeIIOIII', 'd'), ('FRerr_HeIIOIII', 'd'), ('FRs2n_HeIIOIII', 'd'),
                                ('FR_CIIIHeII', 'd'), ('FRerr_CIIIHeII', 'd'), ('FRs2n_CIIIHeII', 'd'),
                                ('FR_CIIIOIII', 'd'), ('FRerr_CIIIOIII', 'd'), ('FRs2n_CIIIOIII', 'd'),
                                ('FR_OIIICIII', 'd'), ('FRerr_OIIICIII', 'd'), ('FRs2n_OIIICIII', 'd'),
                                ('FR_OIIIHeII', 'd'), ('FRerr_OIIIHeII', 'd'), ('FRs2n_OIIIHeII', 'd') ])

    if verbose: print(' - Filling output for '+catreference)
    if verbose: print('   Inserting columns from literature ')
    dataarray['id']        = ids
    dataarray['reference'] = [catreference]*len(ids)
    dataarray['redshift']  = redshifts

    #XXX CIV exists...


    # EW - from Anna
    dataarray['ew_HeII']    = np.array([1.7, 0.4, 0.7, 0.44, 0.4, 1.42, 0.6, 0.94, 0.5, 0.5])
    dataarray['ewerr_HeII'] = np.array([0.06, +99.0, +99.0, 0.02, +99.0, 0.03, +99.0, 0.04, +99.0, +99.0])
    dataarray['ew_OIII']    = np.array([5.05, 0.84, 1.68, 1.89, 0.48, 1.73, 1.17, 1.73, np.nan, 0.9])
    dataarray['ewerr_OIII'] = np.array([0.1, 0.03, 0.08, 0.02, 0.02, 0.05, 0.06, 0.04, np.nan, +99.0])
    dataarray['ew_CIII']    = np.array([14.86, 4.98, 4.0, 12.09, 2.8, 3.3, 8.71, 13.35, 11.33, 3.38])
    dataarray['ewerr_CIII'] = np.array([1.07, 0.59, +99.0, 0.3, +99.0, +99.0, 0.42, 0.52, 0.34, 0.31])

    # line intensities - from Anna
    dataarray['f_HeII']    = np.array([1.45, 1.1, 0.8, 1.22, 1.2,2.37, 0.8, 1.42, 1.0, 0.8])
    dataarray['ferr_HeII'] = np.array([0.06, +99.0, +99.0, 0.05, +99.0, 0.05, +99.0, 0.05, +99.0, +99.0])
    dataarray['f_OIII']    = np.array([3.35, 2.16, 1.69, 5.05, 1.28, 2.68, 1.36, 2.27, np.nan, 1.1])
    dataarray['ferr_OIII'] = np.array([0.04, 0.08, 0.07, 0.05, 0.05, 0.07, 0.06, 0.05, np.nan, +99.0])
    dataarray['f_CIII']    = np.array([11.06, 7.60, 4.1, 22.5, 5.2, 3.9, 6.80, 12.3, 15.9, 3.79])
    dataarray['ferr_CIII'] = np.array([0.63, 0.88, +99.0, 0.40, +99.0, +99.0, 0.23, 0.33, 0.38, 0.30])

    for colname in dataarray.dtype.names:
        if colname.startswith('f_') or colname.startswith('ferr_'):
            goodent = np.where(dataarray[colname] != 99)[0]
            if len(goodent) > 0:
                dataarray[colname][goodent] = dataarray[colname][goodent] * fluxscale

    if verbose: print('   Calculating further columns from literature')
    # fluxratios
    linetypedic = collections.OrderedDict()
    linetypedic['CIII']  = 'single'
    linetypedic['OIII']  = 'single'
    linetypedic['HeII']  = 'single'

    for ii, id in enumerate(ids):
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
                        FR, FRerr = uves.set_ratios('None','denominator_exists',
                                                    fnum,  'dummy', fdenom, ferrdenom)
                    elif (ferrnum != 99) & (ferrdenom == 99):
                        FR, FRerr = uves.set_ratios('numerator_exists','none',
                                                    fnum,  ferrnum, fdenom, 'dummy')
                    else:
                        FR, FRerr = uves.set_ratios('numerator_exists','denominator_exists',
                                                    fnum,  ferrnum, fdenom, ferrdenom)

                    dataarray['FR_'+numerator_line+denominator_line][ii]     = FR
                    dataarray['FRerr_'+numerator_line+denominator_line][ii]  = FRerr
                    if np.abs(FRerr) == 99:
                        dataarray['FRs2n_'+numerator_line+denominator_line][ii]  = np.nan
                    else:
                        dataarray['FRs2n_'+numerator_line+denominator_line][ii]  = FR/FRerr

    return catreference, dataarray

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def restUVline_literaturecollection_Nanaykkara19(fluxscale=1.0,verbose=True):
    """

    --- INPUT ---
    baseid    Number to add to ids to keep them distinct from other literature data sets
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'nan19'
    fluxratiodic        = collections.OrderedDict()
    fluxratiodic['id']  = np.array([])
    baseid              = lce.referencedictionary[catreference][0]

    ids        = np.array([1024, 1036, 1045, 1079, 1273, 3621, 87,   109, 144,    97,   39,  84, 161]) + baseid
    redshifts  = np.array([2.87, 2.69, 2.61, 2.68, 2.17, 3.07, 2.67, 2.2, 4.02, 2.11, 3.96, 3.1, 3.1])

    dataarray = np.array(np.zeros(len(ids))*np.nan,
                         dtype=[('id', 'd'), ('reference', '20a'), ('redshift', 'd'),
                                ('f_HeII', 'd'), ('ferr_HeII', 'd'), ('s2n_HeII', 'd'), ('vshift_HeII', 'd'), ('sigma_HeII', 'd'),
                                ('ew_HeII', 'd'), ('ewerr_HeII', 'd'),
                                ('f_OIII', 'd'), ('ferr_OIII', 'd'), ('s2n_OIII', 'd'), ('vshift_OIII', 'd'), ('sigma_OIII', 'd'),
                                ('ew_OIII', 'd'), ('ewerr_OIII', 'd'),
                                ('f_CIII', 'd'), ('ferr_CIII', 'd'), ('s2n_CIII', 'd'), ('vshift_CIII', 'd'), ('sigma_CIII', 'd'),
                                ('ew_CIII', 'd'), ('ewerr_CIII', 'd'),
                                ('FR_HeIICIII', 'd'), ('FRerr_HeIICIII', 'd'), ('FRs2n_HeIICIII', 'd'),
                                ('FR_HeIIOIII', 'd'), ('FRerr_HeIIOIII', 'd'), ('FRs2n_HeIIOIII', 'd'),
                                ('FR_CIIIHeII', 'd'), ('FRerr_CIIIHeII', 'd'), ('FRs2n_CIIIHeII', 'd'),
                                ('FR_CIIIOIII', 'd'), ('FRerr_CIIIOIII', 'd'), ('FRs2n_CIIIOIII', 'd'),
                                ('FR_OIIICIII', 'd'), ('FRerr_OIIICIII', 'd'), ('FRs2n_OIIICIII', 'd'),
                                ('FR_OIIIHeII', 'd'), ('FRerr_OIIIHeII', 'd'), ('FRs2n_OIIIHeII', 'd') ])

    if verbose: print(' - Filling output for '+catreference)
    if verbose: print('   Inserting columns from literature ')
    dataarray['id']        = ids
    dataarray['reference'] = [catreference]*len(ids)
    dataarray['redshift']  = redshifts

    # # EW - from Anna
    # dataarray['ew_HeII']   = np.array([0.79487179, 0.81081081, 1., 2.5625, 2.51219512, 1.2972973, 1., 2.14,
    #                                    5.77419355, 4.94, 6.56097561, 12.26829268])
    # dataarray['ewerr_HeII'] = np.array([0.15384615, 0.2972973, 0.25, 0.625, +99.0, 0.21621622, 0.125, 0.44,
    #                                    0.77419355, +99.0, 0.85365854, +99.0])
    # dataarray['ew_OIII']   = np.array([0.71794872, 1.2972973, 1.27777778, 2.25, 2.46341463, 1.02702703,
    #                                    1.34375, 7.44, 5.4516129, 0.26, 1.53658537, 1.68292683])
    # dataarray['ewerr_OIII'] = np.array([0.15384615, 0.18918919, 0.22222222, 0.4375, +99.0, 0.13513514,
    #                                     0.15625, 0.86, 0.61290323, +99.0, 0.48780488, +99.0])
    # dataarray['ew_CIII']    = np.array([2.56410256, 5.64864865, 4.94444444, 11.15625, 4.7804878,
    #                                     2.72972973, 3.34375, np.nan, 15.83870968, np.nan, 7.46341463, 2.6097561])
    # dataarray['ewerr_CIII'] = np.array([0.26399052, 0.34611482, 0.45218946, 0.73221091, +99.0,
    #                                     0.30577591, 0.22097087, +99.0,  1.58196127, +99.0,  +99.0, +99.0])
    # dataarray['ew_CIII']    = np.array([1.5384615385, 3.3513, 3.1944, 6.65625, 1.756, 1.9189, 1.78125,
    #                                     np.nan,9.32258, np.nan, 4.53659, 2.5366])
    # dataarray['ewerr_CIII'] = np.array([0.23,0.2162,  0.33, 0.5625, +99.0, 0.2162, 0.15625, +99.0, 1.,
    #                                     +99.0, +99.0, +99.0 ])
    # line intensities - from Anna
    dataarray['f_HeII']    = np.array([177., 142., 156., 290., 217., 213., 59., 54., 48., 306., 153., 161., 318.])
    dataarray['ferr_HeII'] = np.array([52., 53., 58., 91., 79., 45., 12., 13., 12., 55., 33., 29., 39.])
    dataarray['f_OIII']    = np.array([161., 230., 200., 81., 195., 105., 49., 72., 150., 284., 66., 54., 58])
    dataarray['ferr_OIII'] = np.array([48., 54., 56., +99.0, 56., +99.0, 11., 13., 21., 44., +99.0, +99.0, 19.])
    dataarray['f_CIII']    = np.array([514., 736., 599., 222., 673., 464., 111., 134., np.nan, 626., np.nan, 246., 205.])
    dataarray['ferr_CIII'] = np.array([64.54,54.56, 72.92, +99.0, 67.18, +99.0,  15.55, 16.97, +99.0, 79.65, +99.0, +99.0, +99.0])

    for colname in dataarray.dtype.names:
        if colname.startswith('f_') or colname.startswith('ferr_'):
            goodent = np.where(dataarray[colname] != 99)[0]
            if len(goodent) > 0:
                dataarray[colname][goodent] = dataarray[colname][goodent] * fluxscale

    if verbose: print('   Calculating further columns from literature')
    # fluxratios
    linetypedic = collections.OrderedDict()
    linetypedic['CIII']  = 'single'
    linetypedic['OIII']  = 'single'
    linetypedic['HeII']  = 'single'

    for ii, id in enumerate(ids):
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
                        FR, FRerr = uves.set_ratios('None','denominator_exists',
                                                    fnum,  'dummy', fdenom, ferrdenom)
                    elif (ferrnum != 99) & (ferrdenom == 99):
                        FR, FRerr = uves.set_ratios('numerator_exists','none',
                                                    fnum,  ferrnum, fdenom, 'dummy')
                    else:
                        FR, FRerr = uves.set_ratios('numerator_exists','denominator_exists',
                                                    fnum,  ferrnum, fdenom, ferrdenom)

                    dataarray['FR_'+numerator_line+denominator_line][ii]     = FR
                    dataarray['FRerr_'+numerator_line+denominator_line][ii]  = FRerr
                    if np.abs(FRerr) == 99:
                        dataarray['FRs2n_'+numerator_line+denominator_line][ii]  = np.nan
                    else:
                        dataarray['FRs2n_'+numerator_line+denominator_line][ii]  = FR/FRerr

    return catreference, dataarray

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def restUVline_literaturecollection_EoR(fluxscale=1.0,verbose=True):
    """

    --- INPUT ---
    baseid    Number to add to ids to keep them distinct from other literature data sets
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'EoR'
    fluxratiodic        = collections.OrderedDict()
    fluxratiodic['id']  = np.array([])
    baseid              = lce.referencedictionary[catreference][0]

    ids        = np.array([2248, 2, 3]) + baseid
    redshifts  = np.array([6.1, 7.0, 7.1])

    dataarray = np.array(np.zeros(len(ids))*np.nan,
                         dtype=[('id', 'd'), ('reference', '20a'), ('redshift', 'd'),
                                ('f_HeII', 'd'), ('ferr_HeII', 'd'), ('s2n_HeII', 'd'), ('vshift_HeII', 'd'), ('sigma_HeII', 'd'),
                                ('ew_HeII', 'd'), ('ewerr_HeII', 'd'),
                                ('f_OIII', 'd'), ('ferr_OIII', 'd'), ('s2n_OIII', 'd'), ('vshift_OIII', 'd'), ('sigma_OIII', 'd'),
                                ('ew_OIII', 'd'), ('ewerr_OIII', 'd'),
                                ('f_CIII', 'd'), ('ferr_CIII', 'd'), ('s2n_CIII', 'd'), ('vshift_CIII', 'd'), ('sigma_CIII', 'd'),
                                ('ew_CIII', 'd'), ('ewerr_CIII', 'd'),
                                ('f_CIV', 'd'), ('ferr_CIV', 'd'), ('s2n_CIV', 'd'), ('vshift_CIV', 'd'), ('sigma_CIV', 'd'),
                                ('ew_CIV', 'd'), ('ewerr_CIV', 'd'),

                                ('FR_HeIICIII', 'd'), ('FRerr_HeIICIII', 'd'), ('FRs2n_HeIICIII', 'd'),
                                ('FR_HeIIOIII', 'd'), ('FRerr_HeIIOIII', 'd'), ('FRs2n_HeIIOIII', 'd'),
                                ('FR_HeIICIV', 'd'), ('FRerr_HeIICIV', 'd'), ('FRs2n_HeIICIV', 'd'),
                                ('FR_CIIIHeII', 'd'), ('FRerr_CIIIHeII', 'd'), ('FRs2n_CIIIHeII', 'd'),
                                ('FR_CIIIOIII', 'd'), ('FRerr_CIIIOIII', 'd'), ('FRs2n_CIIIOIII', 'd'),
                                ('FR_CIIICIV', 'd'), ('FRerr_CIIICIV', 'd'), ('FRs2n_CIIICIV', 'd'),
                                ('FR_CIVHeII', 'd'), ('FRerr_CIVHeII', 'd'), ('FRs2n_CIVHeII', 'd'),
                                ('FR_CIVOIII', 'd'), ('FRerr_CIVOIII', 'd'), ('FRs2n_CIVOIII', 'd'),
                                ('FR_CIVCIII', 'd'), ('FRerr_CIVCIII', 'd'), ('FRs2n_CIVCIII', 'd'),
                                ('FR_OIIICIII', 'd'), ('FRerr_OIIICIII', 'd'), ('FRs2n_OIIICIII', 'd'),
                                ('FR_OIIIHeII', 'd'), ('FRerr_OIIIHeII', 'd'), ('FRs2n_OIIIHeII', 'd'),
                                ('FR_OIIICIV', 'd'), ('FRerr_OIIICIV', 'd'), ('FRs2n_OIIICIV', 'd') ])

    if verbose: print(' - Filling output for '+catreference)
    if verbose: print('   Inserting columns from literature ')
    dataarray['id']        = ids
    dataarray['reference'] = [catreference]*len(ids)
    dataarray['redshift']  = redshifts

    # EW
    # dataarray['ew_HeII']   = np.array([0.79487179, 0.81081081, 1., 2.5625, 2.51219512, 1.2972973, 1., 2.14,
    #                                    5.77419355, 4.94, 6.56097561, 12.26829268])
    # dataarray['ewerr_HeII'] = np.array([0.15384615, 0.2972973, 0.25, 0.625, +99.0, 0.21621622, 0.125, 0.44,
    #                                    0.77419355, +99.0, 0.85365854, +99.0])
    # dataarray['ew_OIII']   = np.array([0.71794872, 1.2972973, 1.27777778, 2.25, 2.46341463, 1.02702703,
    #                                    1.34375, 7.44, 5.4516129, 0.26, 1.53658537, 1.68292683])
    # dataarray['ewerr_OIII'] = np.array([0.15384615, 0.18918919, 0.22222222, 0.4375, +99.0, 0.13513514,
    #                                     0.15625, 0.86, 0.61290323, +99.0, 0.48780488, +99.0])
    # dataarray['ew_CIII']    = np.array([2.56410256, 5.64864865, 4.94444444, 11.15625, 4.7804878,
    #                                     2.72972973, 3.34375, np.nan, 15.83870968, np.nan, 7.46341463, 2.6097561])
    # dataarray['ewerr_CIII'] = np.array([0.26399052, 0.34611482, 0.45218946, 0.73221091, +99.0,
    #                                     0.30577591, 0.22097087, +99.0,  1.58196127, +99.0,  +99.0, +99.0])
    # dataarray['ew_CIII']    = np.array([1.5384615385, 3.3513, 3.1944, 6.65625, 1.756, 1.9189, 1.78125,
    #                                     np.nan,9.32258, np.nan, 4.53659, 2.5366])
    # dataarray['ewerr_CIII'] = np.array([0.23,0.2162,  0.33, 0.5625, +99.0, 0.2162, 0.15625, +99.0, 1.,
    #                                     +99.0, +99.0, +99.0 ])
    # line intensities
    dataarray['f_HeII']    = np.array([150,  210,    126])
    dataarray['ferr_HeII'] = np.array([+99,  +99,    29])
    dataarray['f_OIII']    = np.array([440,  180,    np.nan])
    dataarray['ferr_OIII'] = np.array([60,   70,     99])
    dataarray['f_CIII']    = np.array([500,  np.nan, 175])
    dataarray['ferr_CIII'] = np.array([+99,  +99,    +99])
    dataarray['f_CIV']     = np.array([1140, 790,    270])
    dataarray['ferr_CIV']  = np.array([90,   90,     +99])


    for colname in dataarray.dtype.names:
        if colname.startswith('f_') or colname.startswith('ferr_'):
            goodent = np.where(dataarray[colname] != 99)[0]
            if len(goodent) > 0:
                dataarray[colname][goodent] = dataarray[colname][goodent] * fluxscale

    if verbose: print('   Calculating further columns from literature')
    # fluxratios
    linetypedic = collections.OrderedDict()
    linetypedic['CIV']   = 'single'
    linetypedic['CIII']  = 'single'
    linetypedic['OIII']  = 'single'
    linetypedic['HeII']  = 'single'

    for ii, id in enumerate(ids):
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
                        FR, FRerr = uves.set_ratios('None','denominator_exists',
                                                    fnum,  'dummy', fdenom, ferrdenom)
                    elif (ferrnum != 99) & (ferrdenom == 99):
                        FR, FRerr = uves.set_ratios('numerator_exists','none',
                                                    fnum,  ferrnum, fdenom, 'dummy')
                    else:
                        FR, FRerr = uves.set_ratios('numerator_exists','denominator_exists',
                                                    fnum,  ferrnum, fdenom, ferrdenom)

                    dataarray['FR_'+numerator_line+denominator_line][ii]     = FR
                    dataarray['FRerr_'+numerator_line+denominator_line][ii]  = FRerr
                    if np.abs(FRerr) == 99:
                        dataarray['FRs2n_'+numerator_line+denominator_line][ii]  = np.nan
                    else:
                        dataarray['FRs2n_'+numerator_line+denominator_line][ii]  = FR/FRerr

    return catreference, dataarray


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =