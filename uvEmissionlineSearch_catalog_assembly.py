# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Collectin and assembling catalogs for release with UVES paper
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import sys
import os
import numpy.lib.recfunctions as rfn
import astropy.io.fits as afits
import numpy as np
import uvEmissionlineSearch_catalog_assembly as uca
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def assemble_catalogs(verbose=True):
    """
    Collecting velocity offset measurements in a catalog to be released with the UVES paper

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch_catalog_assembly as uca
    uca.assemble_catalogs()

    """
    if verbose: print(' - Assembling the catalogs to acompony Schmidt et al. (2021)')
    basename = 'Schmidt_etal_2021_AandA_XXXX_YYYY'

    uca.uves_catalog_assemble(basename,verbose=verbose)
    uca.lit_catalog_assemble(basename,verbose=verbose)
    uca.dv_catalog_assemble(basename,verbose=verbose)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def dv_catalog_assemble(basename,verbose=True):
    """
    Collecting velocity offset measurements in a catalog to be released with the UVES paper

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch_catalog_assembly as uca

    """
    if verbose: print('========================================================\n - Assembling the dv_Lya literature catalog')

    dvcat = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/vshift_Lya_and_MUV.txt'
    dvdat = np.genfromtxt(dvcat,skip_header=0,names=True,dtype='d,d,d,40a,d,d,d,d,d,d,d,d,d,d,d,40a')

    collist_trans = {'id':'id','ra':'ra', 'dec':'dec','name':'objname','reference':'Reference','redshift':'redshift','magabs_uv':'magabsUV','magabserr_uv':'magabsUVerr','magapp_uv':'magappUV','magapp_uv_err':'magappUVerr','vshift_lya':'vshift_Lya','vshifterr_lya':'vshifterr_Lya','fwhm_lya':'FWHM_lya','fwhmerr_lya':'FWHM_lyaerr','peaksep_lya':'peaksep_lya','peakseperr_lya':'peaksep_lya_err'}

    collist_trans_back = {}
    for key in collist_trans.keys():
        collist_trans_back[collist_trans[key]] = key

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    collist_out  = []
    colnames     = dvdat.dtype.names
    for cc, colname in enumerate(colnames):
        if colname in collist_trans_back.keys():
            collist_out.append(collist_trans_back[colname])

    dvdat['id'] = np.arange(1,len(dvdat['id'])+1,1)

    ent_reorder       = [0, 1, 2, 3, len(collist_out)-1]+list(np.arange(4,len(collist_out)-1,1))
    collist_reordered = [collist_out[i] for i in ent_reorder]


    outarray = uca.build_dataarray(dvdat, collist_reordered, colnametranslation=collist_trans, verbose=True)

    outputfile = basename+'_dvLya_literature.txt'
    write_catalogs_to_disk(outputfile, outarray, literaturecontent=True, verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def lit_catalog_assemble(basename,verbose=True):
    """
    Collecting UV emission literature compilation in a catalog to be released with the UVES paper

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch_catalog_assembly as uca

    """
    if verbose: print('========================================================\n - Assembling the UV emission literature catalog')
    litcat = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/literaturecollection_emissionlinestrengths.fits'
    litdat = afits.open(litcat)[1].data
    colnames = litdat.dtype.names
    collist  = ['id','ra','dec','name','reference','redshift',
                'f','ferr','s2n',
                'fr','frerr','frs2n',
                'ew0','ew0err']

    collist_out   = []

    for cc, colname in enumerate(colnames):
        if colname.split('_')[0].lower() in collist:
            collist_out.append(colname)

    ent_reorder       = [0, 3, 4, 1, 2]+list(np.arange(5,len(collist_out),1))
    collist_reordered = [collist_out[i] for i in ent_reorder]

    outarray = uca.build_dataarray(litdat, collist_reordered, verbose=verbose)

    outputfile = basename+'_UVemissionlines_literature.txt'
    uca.write_catalogs_to_disk(outputfile, outarray, literaturecontent=True, verbose=verbose)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def uves_catalog_assemble(basename,verbose=True):
    """
    Collecting UVES MUSE sources in a catalog to be released with the UVES paper

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch_catalog_assembly as uca

    """
    if verbose: print('========================================================\n - Assembling the UVES main source catalog')
    datestr       = '200213'
    parentdir     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/back2backAnalysis_'+datestr+'/'
    infocat       = parentdir+'../objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    infodat       = afits.open(infocat)[1].data
    infodat       = infodat[np.where((infodat['id']<4.9e8) | (infodat['id']>5.9e8))[0]] # ignoring UDF MW mock ids
    mastercat     = parentdir+'results_master_catalog_version'+datestr+'.fits'
    masterdat     = afits.open(mastercat)[1].data

    for ii, id in enumerate(masterdat['id']): # double checking that order of objects is the same
        if infodat['id'][ii] != id:
            sys.exit('There was a mismatch in ids for entry '+str(ii))

    collist_i_trans = {'id':'id','ra':'ra','dec':'dec','redshift':'redshift','lead_line':'leadline',
                       'confidence':'leadlineConf','id_duplication':'duplicationID',
                       'id_kerutt':'ID_jk100',
                       'f_lya':'F_3KRON_jk100','ferr_lya':'F_3KRON_ERR_jk100',
                       'ew0_lya':'EW_0_beta_own_median_jk100','ew0err_lya':'EW_0_beta_own_median_error_jk100',
                       'magabs_uv':'abs_mag_UV_cont_own_median_jk100','magabserr_uv':'abs_mag_UV_cont_own_median_error_jk100',
                       'fwhm_lya':'fwhm_kms_jk100','fwhmerr_lya':'fwhm_kms_err_jk100',
                       'peaksep_lya':'peak_sep_kms_jk100','peakseperr_lya':'peak_sep_kms_err_jk100',
                       'id_guo':'id_guo','sep_guo':'sep_guo',
                       'id_skelton':'id_skelton','sep_skelton':'sep_skelton',
                       'id_rafelski':'id_rafelski','sep_rafelski':'sep_rafelski',
                       'id_laigle':'id_Laigle','sep_laigle':'sep_Laigle'}
    collist_i_trans_back = {}
    for key in collist_i_trans.keys():
        collist_i_trans_back[collist_i_trans[key]] = key

    collist_m       = ['f','ferr','s2n','sigma','vshift',
                       'fr','frerr','frs2n',
                       'ew0','ew0err',
                       'contband','contmagAB','contmagABerr','photref']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    collist_i_out  = []
    colnames_i     = infodat.dtype.names
    for cc, colname in enumerate(colnames_i):
        if colname in collist_i_trans_back.keys():
            collist_i_out.append(collist_i_trans_back[colname])
    outarray_i = uca.build_dataarray(infodat, collist_i_out, colnametranslation=collist_i_trans, verbose=True)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    collist_m_out  = []
    colnames_m     = masterdat.dtype.names

    for cc, colname in enumerate(colnames_m):
        if colname.split('_')[0].lower() in collist_m:
            collist_m_out.append(colname)
    outarray_m = uca.build_dataarray(masterdat, collist_m_out, verbose=True)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    datarray_combined = rfn.merge_arrays([outarray_i,outarray_m], flatten = True, usemask = False)
    collist_out       = outarray_i.dtype.names+outarray_m.dtype.names
    ent_reorder       = [0, 1, 2, 3]+[12,13,14]+list(np.arange(26,len(collist_out),1))\
                        +[15,16,17,18,19,20,21,22,23,24,25]\
                        +[4,5,6,7,8,9,10,11]
    collist_reordered = [collist_out[i] for i in ent_reorder]
    outarray          = uca.build_dataarray(datarray_combined, collist_reordered, verbose=True)

    outputfile = basename+'_MUSE_source_catalog.txt'
    uca.write_catalogs_to_disk(outputfile, outarray, literaturecontent=False, verbose=verbose)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_dataarray(dataarray, colnames, colnametranslation=None, verbose=True):
    """
    Function building structured data array from input


    --- EXAMPLE OF USE ---
    dataarray = lce.build_dataarray(catreference, linetypedic, datadic) # see example in lce.data_sen17()

    """
    if verbose: print(' - Bulding and filling structure data array')
    Ncol          = len(colnames)
    Ndataobj      = len(dataarray)
    dtypebuild    = [(colname.lower(), 'd') for colname in colnames]
    for cc, colname in enumerate(colnames):
        if 'id' in colname.lower():
            dtypebuild[cc] = (colname.lower(), '>i8')
        if colname in ['name','reference','lead_line']:
            dtypebuild[cc] = (colname.lower(), '30a')

    outarray     = np.array(np.zeros(Ndataobj)*np.nan,dtype=dtypebuild)

    if verbose: print('   Inserting data from input dataarray')
    for cc, colname in enumerate(colnames):
        if colnametranslation is not None:
            outarray[colname.lower()] = dataarray[colnametranslation[colname]]
        else:
            outarray[colname.lower()] = dataarray[colname]

    return outarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def write_catalogs_to_disk(outputfile, outputdataarray, literaturecontent=True, verbose=True):
    """

    """
    outputdataarray = np.sort(outputdataarray, order='id') # sorting array by ID
    outputdir = '/Users/kschmidt/work/publications/MUSE_UVemissionlines/catalog_releases/'
    outtxt    = outputdir+outputfile
    outfits   = outputdir+outputfile.replace('.txt','.fits')
    #-------------------------------------------------------------------------------------------------------------
    if verbose: print(' - Initializing the output file: \n   '+outtxt)
    fout = open(outtxt,'w')
    if verbose: print('   Putting together header for ascii output')
    fout.write('# Value-added source catalog released with Schmidt et al. (2021) A&A XXXX:YYYY \n')
    fout.write('# Upper and lower limits are given as values with uncertainty of +99 or -99, respectively. \n')
    if literaturecontent:
        fout.write('# When using this catalog please cite Schmidt et al. (2021) A&A XXXX:YYYY '
                   'and the original papers the data were assmbled from in the reference column (see Schmidt et al. 2021) \n')
    else:
        fout.write('# When using this catalog please cite Schmidt et al. (2021) A&A XXXX:YYYY \n')
    fout.write('# \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    colnames     = outputdataarray.dtype.names
    Ncols        = len(colnames)
    if verbose: print('   The output files will contain '+str(Ncols)+' columns ')
    fout.write('# The catalog contains the following '+str(Ncols)+' columns:\n')
    fout.write('# '+(' '.join([str("%20s" % colname) for colname in list(colnames)])).replace('reference','                    reference')+'  \n')
    if verbose: print('   The total number of objects in the output is '+str(len(outputdataarray['id'])))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Writing output data array to ascii file \n   '+outputfile)
    dtypetranslator = {'>i8':'%20i',  'float64':'%20.4f', '|S30':'%30s'}

    for oo, id in enumerate(outputdataarray['id']):
        outstr = ' '
        for cc, colval in enumerate(outputdataarray[oo].tolist()):
            strfmt = dtypetranslator[str(outputdataarray.dtype[cc])]
            if outputdataarray.dtype.names[cc] in  ['ra','dec']:
                strfmt = '%20.10f'

            if strfmt == '%30s':
                outstr = outstr +' '+ str(strfmt % colval.decode('UTF-8'))
            else:
                outstr = outstr +' '+ str(strfmt % colval)

        fout.write(outstr+' \n')
    fout.close()

    #-------------------------------------------------------------------------------------------------------------
    if verbose: print(' - Creating fits version of output: \n   '+outfits)
    if literaturecontent:
        fitsformat = ['K','D','D','20A','30A','D'] + ['D']*(Ncols-6)
    else:
        fitsformat      = ['K','D','D','D','10A','K','K'] + ['D']*(Ncols-15) + ['K', 'D']*4
        fitsformat[-19] = 'K' # The Kerutt ID
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('   Reading ascii file ')
    data    = np.genfromtxt(outtxt,names=colnames,skip_header=6,comments='#',dtype=outputdataarray.dtype)
    keys    = data.dtype.names
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('   Initialize and fill dictionary with data')
    datadic = {}
    for kk in keys:
        datadic[kk] = []
        try:
            lenarr = len(np.asarray(data[kk]))
            datadic[kk] = np.asarray(data[kk])
        except: # if only one row of data is to be written
            datadic[kk] = np.asarray([data[kk]])

    if verbose: print('   Found '+str(len(keys))+' columns to insert into fits binary table')

    if len(fitsformat) != len(keys):
        fitsformat = np.asarray([fitsformat]*len(keys))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('   Writing to fits table')
    columndefs = []
    for kk, key in enumerate(keys):
        columndefs.append(afits.Column(name=key  , format=fitsformat[kk], unit=uca.get_unit(key), array=datadic[key]))

    cols     = afits.ColDefs(columndefs)
    tbhdu    = afits.BinTableHDU.from_columns(cols)  # creating table header
    hdu      = afits.PrimaryHDU()                    # creating primary (minimal) header
    thdulist = afits.HDUList([hdu, tbhdu])           # combine primary and table header to hdulist
    thdulist.writeto(outfits,overwrite=True)      # write fits file

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_unit(colname):
    """
    Returning units for fits table columns. Formats described in PDF linked from https://fits.gsfc.nasa.gov/fits_standard.html
    """
    unitdic = {'id':'','ra':'deg','dec':'deg','redshift':'','lead_line':'','name':'','reference':'',
               'confidence':'','id_duplication':'','id_kerutt':'',
               'f_lya':'10**-20 erg s-1 cm-2','ferr_lya':'10**-20 erg s-1 cm-2',
               'ew0_lya':'Angstrom','ew0err_lya':'Angstrom',
               'magabs_uv':'mag','magabserr_uv':'mag',
               'fwhm_lya':'km s-1','fwhmerr_lya':'km s-1',
               'peaksep_lya':'km s-1','peakseperr_lya':'km s-1',
               'id_guo':'id_guo','sep_guo':'sep_guo',
               'id_skelton':'','sep_skelton':'arcsec',
               'id_rafelski':'','sep_rafelski':'arcsec',
               'id_laigle':'','sep_laigle':'arcsec',
               'contmagAB':'mag','contmagABerr':'mag',
               'magapp_uv':'mag','magapp_uv_err':'mag','vshift_lya':'km s-1',
               'vshifterr_lya':'km s-1'}

    if colname in unitdic.keys():
        returnunit = unitdic[colname]
    elif colname.startswith('f_') or colname.startswith('ferr_'):
        returnunit = '10**-20 erg s-1 cm-2'
    elif colname.startswith('fr_') or colname.startswith('frerr_') or colname.startswith('s2n_') or colname.startswith('frs2n_') or\
            colname.startswith('contband_') or colname.startswith('photref_'):
        returnunit = ''
    elif colname.startswith('ew0_') or colname.startswith('ew0err_') or colname.startswith('sigma_') or colname.startswith('sigmaerr_'):
        returnunit = 'Angstrom'
    elif colname.startswith('vshift_') or colname.startswith('vshifterr_'):
        returnunit = 'km s-1'
    else:
        sys.exit("uca.get_unit() couldn't find a unit to return for '"+colname+"'")

    return returnunit
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
