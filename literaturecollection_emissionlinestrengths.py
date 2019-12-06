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
import matplotlib.pyplot as plt
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
    refdic['sen17'] = [1e10,    'Senchyna et al. 2017',       'v']
    refdic['nan19'] = [2e10,    'Nanaykkara et al. 2019',     '^']
    refdic['eor']   = [3e10,    'EoR objects',                '<']
    refdic['dummy'] = [99e10,   'dummy',                      '>']
    refdic['dummy'] = [99e10,   'dummy',                      '8']
    refdic['dummy'] = [99e10,   'dummy',                      's']
    refdic['dummy'] = [99e10,   'dummy',                      'p']
    refdic['dummy'] = [99e10,   'dummy',                      'P']
    refdic['dummy'] = [99e10,   'dummy',                      '*']
    refdic['dummy'] = [99e10,   'dummy',                      'h']
    refdic['dummy'] = [99e10,   'dummy',                      'H']
    refdic['dummy'] = [99e10,   'dummy',                      '+']
    refdic['dummy'] = [99e10,   'dummy',                      'x']
    refdic['dummy'] = [99e10,   'dummy',                      'D']
    refdic['dummy'] = [99e10,   'dummy',                      'd']
    refdic['dummy'] = [99e10,   'dummy',                      '1']
    refdic['dummy'] = [99e10,   'dummy',                      '2']
    refdic['dummy'] = [99e10,   'dummy',                      '3']
    refdic['dummy'] = [99e10,   'dummy',                      '4']
    refdic['dummy'] = [99e10,   'dummy',                      '$\\alpha$']
    refdic['dummy'] = [99e10,   'dummy',                      (5, 0, 180)] # pentagon rotated 180 degrees

    # --- MUSE-Wide def: ---
    # CDFS and COSMOS:  'o'
    # UDF:              'D'
    # UDF10:            'X'

    # --- potentially outer symbols def: ---
    #  D    GLASS
    #  s    stack/composite
    #  *    AGN
    #  o    MUSE data

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
def generate_literature_fitscatalog(verbose=True):
    """
    Collecting literature values and storing them in a single file

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import literaturecollection_emissionlinestrengths as lce
    lce.generate_literature_fitscatalog()

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
    # # a simple dictionary containing the column locations in the output array (indexes)
    # colents    = {}
    # for oo, colname in enumerate(fluxratiodic.keys()):
    #     colents[colname] = oo
    if verbose: print('   The output file will contain '+str(Ncols)+' columns ')
    fout.write('# This file contains the following '+str(Ncols)+' columns:\n')
    fout.write('# id      reference '+' '.join([str("%20s" % colname) for colname in fluxratiodic.keys()[2:]])+'  \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Collecting literature data arrays and appending them')
    refdic                 = lce.referencedictionary()

    dataref, dataarray     = lce.data_nan19(fluxscale=1.0,verbose=True)
    outputdataarray        = dataarray
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0]][1])

    dataref, dataarray     = lce.data_nan19(fluxscale=1.0,verbose=True)
    outputdataarray        = np.append(outputdataarray,dataarray)
    if verbose: print('   Added data from   '+refdic[dataarray['reference'][0]][1])

    #
    # appenddatadic = collections.OrderedDict()
    # #------------------------------------------------
    # dataref, appenddata     = lce.data_sen17(fluxscale=1e5,verbose=True)
    # appenddatadic[dataref]  = appenddata
    #
    # dataref, appenddata     = lce.data_nan19(fluxscale=1.0,verbose=True)
    # appenddatadic[dataref]  = appenddata
    #
    # dataref, appenddata     = lce.data_eor(fluxscale=1.0,verbose=True)
    # appenddatadic[dataref]  = appenddata
    # #------------------------------------------------
    # # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # if verbose: print('\n - Filling the columns of the output file with data. Adding data for: ')
    # for appendkey in appenddatadic.keys():
    #     appenddata = appenddatadic[appendkey]
    #     if verbose: print('      '+appendkey+'      ('+refdic[appendkey][1])+')'
    #     Nobjappend     = len(appenddata)
    #     fluxratioarray = np.zeros([Nobjappend,Ncols])*np.nan
    #     for objent, objappend in enumerate(appenddata):
    #         for acol in appenddata.dtype.names:
    #             if (acol in fluxratiodic.keys()) & (acol != 'reference'):
    #                 fluxratioarray[objent,colents[acol]] = appenddata[acol][objent]
    #             else:
    #                 if (acol != 'reference'):
    #                     if verbose: print(' WARNING: The column "'+str(acol)+'" has no match in output file')
    #
    #     for ll in np.arange(Nobjappend):
    #         outstr = str(int(fluxratioarray[ll,0]))+'  '+appendkey+' '+' '.join([str("%10.4f" % ff) for ff in fluxratioarray[ll,2:]])
    #         fout.write(outstr+' \n')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Writing output data array to ascii file \n   '+outputfile)
    for oo, id in enumerate(outputdataarray['id']):
        outstr = str(int(id))+' '+outputdataarray['reference'][oo]+' '+\
                 ' '.join([str("%20.4f" % ff) for ff in outputdataarray[oo].tolist()[2:]])
        fout.write(outstr+' \n')
    fout.close()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # fmt = '20a,20a,'+','.join((Ncols-2)*['d'])
    # dataarray = np.genfromtxt(outputfile,skip_header=5,dtype=fmt,comments='#',names=True)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Creating fits version of output: \n   '+outputfile.replace('.txt','.fits'))
    fitsformat = ['K','20A'] + ['D']*(Ncols-2)
    fitsoutput = ascii2fits(outputfile,asciinames=True,skip_header=5,fitsformat=fitsformat,verbose=False)
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
    dtypebuild[1] = ('reference', '20a')
    dataarray     = np.array(np.zeros(Ndataobj)*np.nan,dtype=dtypebuild)

    if verbose: print('   Inserting data from literature into master dictionary ')
    for ll, mastercol in enumerate(masterdic.keys()):
        if mastercol in datadic.keys():
            dataarray[mastercol] = datadic[mastercol]

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

                # if (numerator_line == 'CIII1') & (denominator_line == 'HeII'): pdb.set_trace()

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

                    # only attempt to fill empty columns to avoid overwriting values from input
                    if np.isnan(dataarray['FR_'+numerator_line+denominator_line][ii]) & \
                            np.isnan(dataarray['FRerr_'+numerator_line+denominator_line][ii]):
                        dataarray['FR_'+numerator_line+denominator_line][ii]     = FR
                        dataarray['FRerr_'+numerator_line+denominator_line][ii]  = FRerr
                        if np.abs(FRerr) == 99:
                            dataarray['FRs2n_'+numerator_line+denominator_line][ii]  = S2Nlim
                        else:
                            dataarray['FRs2n_'+numerator_line+denominator_line][ii]  = FR/FRerr


    return dataarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_master_datadictionary(verbose=True):
    """
    Defining the master dictionary to fill with data from the literature
    """
    datadictionary = collections.OrderedDict()
    datadictionary['id']        = np.array([])
    datadictionary['reference'] = np.array([])
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

        # if  linetypedic[numerator_line] == 'doublet':
        #     datadictionary['f_'+numerator_line+'1']      = np.array([])
        #     datadictionary['ferr_'+numerator_line+'1']   = np.array([])
        #     datadictionary['s2n_'+numerator_line+'1']   = np.array([])
        #     datadictionary['f_'+numerator_line+'2']      = np.array([])
        #     datadictionary['ferr_'+numerator_line+'2']   = np.array([])
        #     datadictionary['s2n_'+numerator_line+'2']      = np.array([])
        #     datadictionary['FR_'+numerator_line+'1'+numerator_line+'2']     = np.array([])
        #     datadictionary['FRerr_'+numerator_line+'1'+numerator_line+'2']  = np.array([])
        #     datadictionary['FRs2n_'+numerator_line+'1'+numerator_line+'2']  = np.array([])

        for kk, denominator_line in enumerate(componentlist):
            if numerator_line == denominator_line:
                continue
            elif numerator_line in denominator_line:
                continue
            elif denominator_line in numerator_line:
                continue
            else:
                datadictionary['FR_'+numerator_line+denominator_line]     = np.array([])
                datadictionary['FRerr_'+numerator_line+denominator_line]  = np.array([])
                datadictionary['FRs2n_'+numerator_line+denominator_line]  = np.array([])

                # if (linetypedic[denominator_line] == 'single'):
                #     datadictionary['FR_'+numerator_line+'1'+denominator_line]     = np.array([])
                #     datadictionary['FRerr_'+numerator_line+'1'+denominator_line]  = np.array([])
                #     datadictionary['FRs2n_'+numerator_line+'1'+denominator_line]  = np.array([])
                #     datadictionary['FR_'+numerator_line+'2'+denominator_line]     = np.array([])
                #     datadictionary['FRerr_'+numerator_line+'2'+denominator_line]  = np.array([])
                #     datadictionary['FRs2n_'+numerator_line+'2'+denominator_line]  = np.array([])
                # if (linetypedic[denominator_line] == 'doublet'):
                #     datadictionary['FR_'+numerator_line+denominator_line+'1']     = np.array([])
                #     datadictionary['FRerr_'+numerator_line+denominator_line+'1']  = np.array([])
                #     datadictionary['FRs2n_'+numerator_line+denominator_line+'1']  = np.array([])
                #     datadictionary['FR_'+numerator_line+denominator_line+'2']     = np.array([])
                #     datadictionary['FRerr_'+numerator_line+denominator_line+'2']  = np.array([])
                #     datadictionary['FRs2n_'+numerator_line+denominator_line+'2']  = np.array([])
                # if (linetypedic[numerator_line] == 'doublet') & \
                #         (linetypedic[denominator_line] == 'doublet'):
                #     datadictionary['FR_'+numerator_line+'1'+denominator_line+'1']     = np.array([])
                #     datadictionary['FRerr_'+numerator_line+'1'+denominator_line+'1']  = np.array([])
                #     datadictionary['FRs2n_'+numerator_line+'1'+denominator_line+'1']  = np.array([])
                #     datadictionary['FR_'+numerator_line+'1'+denominator_line+'2']     = np.array([])
                #     datadictionary['FRerr_'+numerator_line+'1'+denominator_line+'2']  = np.array([])
                #     datadictionary['FRs2n_'+numerator_line+'1'+denominator_line+'2']  = np.array([])
                #     datadictionary['FR_'+numerator_line+'2'+denominator_line+'1']     = np.array([])
                #     datadictionary['FRerr_'+numerator_line+'2'+denominator_line+'1']  = np.array([])
                #     datadictionary['FRs2n_'+numerator_line+'2'+denominator_line+'1']  = np.array([])
                #     datadictionary['FR_'+numerator_line+'2'+denominator_line+'2']     = np.array([])
                #     datadictionary['FRerr_'+numerator_line+'2'+denominator_line+'2']  = np.array([])
                #     datadictionary['FRs2n_'+numerator_line+'2'+denominator_line+'2']  = np.array([])
    if verbose: print(' - Assembled and build a master dictionary with '+str(len(datadictionary.keys()))+' keys')
    return datadictionary
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def load_literaturecatalog(verbose=True):
    """
    Loading the literature catalog.

    """
    outdir      = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/'
    fitscatalog = outdir+'literaturecollection_emissionlinestrengths.fits'
    if verbose: print('\n - lce.load_literaturecatalog() Loading data from:\n   '+fitscatalog)
    literaturedata = afits.open(fitscatalog)[1].data

    return fitscatalog, literaturedata

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def data_nan19(fluxscale=1.0,verbose=True):
    """
    Data collected from Nanaykkara et al. (2019)

    Non-existing data is provided as NaNs, 3-sigma upper/lower limits are given in flux columns with errors of +/-99

    --- INPUT ---
    fluxscale   Flux scale to bring fluxes and flux errors to 1e-20 erg/s/cm2
    verbose     Toggle verbosity

    """
    catreference        = 'nan19'
    refdic              = lce.referencedictionary()
    if verbose: print('\n - Assembling the data from '+refdic[catreference][1])
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
    if verbose: print('   Putting together measurements from '+str(len(datadic['id']))+' objects ')

    # line intensities - from Anna
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
    datadic['FR_'+linename], datadic['FRrerr_'+linename], \
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
    datadic['FR_'+linename], datadic['FRrerr_'+linename], \
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
    datadic['FR_'+linename], datadic['FRrerr_'+linename], \
    datadic['EW0_'+linename], datadic['EW0err_'+linename] = \
        lce.calc_doubletValuesFromSingleComponents(datadic['f_'+linename+'1'],datadic['ferr_'+linename+'1'],
                                                   datadic['f_'+linename+'2'],datadic['ferr_'+linename+'2'],
                                                   EW1=datadic['EW0_'+linename+'1'], EW1err=datadic['EW0err_'+linename+'1'],
                                                   EW2=datadic['EW0_'+linename+'2'], EW2err=datadic['EW0err_'+linename+'2'])

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

    # id   z_TN       HeII_ew   HeII_ew_nmad     CIIIa_ew CIIIa_ew_nmad    CIIIb_ew  CIIIb_ew_nmad     OIIIa_ew OIIIa_ew_nmad    OIIIb_ew  OIIIb_ew_nmad   SiIIIa_ew   SiIIIa_ew_nmad   SiIIIb_ew   SiIIIb_ew_nmad
    # 1024.0   2.9           3.1   0.6              6.0   ,    0.9   ,           4.0   ,   0.5   ,            1.3 , 0.6  ,                2.8 , 0.6  ,              3.6   ,  0.7   ,                  3.2   ,  1.1,
    # 1036.0   2.7           3.0   1.1              12.4  ,    0.8   ,           8.5   ,   1.0   ,            2.1 , 0.7  ,                4.8 , 0.7  ,              4.6   ,  0.6   ,                  2.5   ,  0.8,
    # 1045.0   2.6           3.6   0.9              11.5  ,    1.2   ,           6.3   ,   1.1   ,            3.2 , 0.8  ,                4.6 , 0.8  ,              4.1   ,  0.7   ,                  3.8   ,  0.8,
    # 1079.0   2.7           9.1   2.6              0.7   ,    0.6   ,           1.5   ,   0.7   ,            0.2 , 0.7  ,                0.2 , 0.4  ,              0.0   ,  0.4   ,                  0.4   ,  0.6,
    # 1273.0   2.2           8.2   2.0              21.3  ,    1.8   ,           14.4  ,   1.5   ,            3.6 , 1.1  ,                7.2 , 1.4  ,              7.3   ,  1.9   ,                  3.0   ,  1.4,
    # 3621.0   3.1           10.3  -99.0            7.2   ,    -99.0 ,           12.4  ,   -99.0 ,            10.2, -99.0,                10.1, -99.0,              14.5  ,  -99.0 ,                    7.0 ,  -99.0,
    # 87.0     2.7           4.8   0.8              7.1   ,    0.8   ,           3.0   ,   0.8   ,            2.1 , 0.5  ,                3.8 , 0.5  ,              2.9   ,  0.5   ,                  4.1   ,  0.7,
    # 109.0    2.2           3.2   0.4              5.7   ,    0.5   ,           5.0   ,   0.5   ,            1.5 , 0.5  ,                4.3 , 0.5  ,              4.5   ,  0.7   ,                  2.2   ,  0.7,
    # 144.0    4.0           10.7  2.2              np.nan,    np.nan,           np.nan,   np.nan,            9.7 , 1.8  ,                37.2, 4.3  ,              np.nan,  np.nan,                  np.nan,  np.nan,
    # 97.0     2.1           17.9  2.4              28.9  ,    3.1   ,           20.2  ,   3.8   ,            8.8 , 1.4  ,                16.9, 1.9  ,              7.7   ,  1.9   ,                  3.0   ,  1.8,
    # 39.0     4.0           24.7  -99.0            np.nan,    np.nan,           np.nan,   np.nan,            5.1 , -99.0,                1.3 , -99.0,              np.nan,  np.nan,                  np.nan,  np.nan,
    # 84.0     3.1           26.9  3.5              18.6  ,    -99.0 ,           12.0  ,   -99.0 ,            8.5 , 2.4  ,                6.3 , 2.0  ,              39.4  ,  -99.0 ,                  23.3  ,  -99.0,
    # 161.0    3.1           50.3  -99.0            10.4  ,    -99.0 ,           0.3   ,   -99.0 ,            1.7 , -99.0,                6.9 , -99.0,              6.6   ,  -99.0 ,                  12.0  ,  -99.0,
    if verbose: print('   Returning catalog reference and data array')
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
