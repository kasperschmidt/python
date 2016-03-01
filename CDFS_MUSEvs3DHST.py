# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import os
import commands
import sys
import glob
import collections
import pyfits
import kbsutilities as kbs
import numpy as np
import manipulate_3DHSTdata as m3d
import CDFS_MUSEvs3DHST as cm3
import fits2ascii as f2a
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def extract_objectsinfo(MUSEidlist,cmfile='./MUSEcdfs_cm2_3DHSTgoods.fits',
                        crossmatchtol=0.1,output_basename='./objectinfo_output',clobber=False,
                        printcolnames=False,verbose=True):
    """

    Get information for a set MUSE objects from the 3D-HST catalogs. Assumes a crossmatch file has
    already been generated (with, e.g., kbsutilities.crossmatch2cat() and kbsutilities.crossmatch())

    --- INPUT ---
    MUSEidlist      List of MUSE IDs to extract information for
    cmfile          Fits fiel with crossmatch to 3D-HST catalogs
    crossmatchtol   Tolerance for crossmatch. Deafult is to only extract information of r_match <0.1 arc sec
    printcolnames   Set to True to print names of columns extracted from 3D-HST catalogs
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---
    import CDFS_MUSEvs3DHST as cm3
    cm3.extract_objectsinfo([10101006])

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading information from crossmatch file'
    cmdat = pyfits.open(cmfile)[1].data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading 3D-HST catalogs to extract information from'
    path3DHST = '/Volumes/DATABCKUP3/GOODSdata/3D-HST/GOODSS_WFC3_V4.1.5/goodss_3dhst_v4.1.5_catalogs/'

    #cat_phot    = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
    cat_dupspec = 'goodss_3dhst.v4.1.5.duplicates_2d.dat'        # list of duplicate spectra
    cat_fast    = 'goodss_3dhst.v4.1.5.zbest.fout'               # fast output using best redshift
    #cat_rf      = 'goodss_3dhst.v4.1.5.zbest.rf'                 # rest frame colors using best redshift
    cat_sfr     = 'goodss_3dhst.v4.1.5.zbest.sfr'                # SFR estimates using best redshift
    cat_zbest   = 'goodss_3dhst.v4.1.5.zbest.fits'               # information about how zbest was obtained
    cat_zinfo   = 'goodss_3dhst.v4.1.5.zfit.linematched.fits'    # redshift infpormation
    cat_ELflux  = 'goodss_3dhst.v4.1.5.linefit.linematched.fits' # line flux measurements for best spectra

    #dat_phot    = pyfits.open(cat_phot)[1].data
    dat_dupspec = np.genfromtxt(path3DHST+cat_dupspec,comments='#',
                                names=['id','s1','s2','s3','s4','s5','s6','s7'],
                                dtype='f,40a,40a,40a,40a,40a,40a,40a')
    dat_fast    = np.genfromtxt(path3DHST+cat_fast,names=True,comments='#')
    #dat_rf      = np.genfromtxt(path3DHST+cat_rf,names=True,comments='#')
    dat_sfr     = np.genfromtxt(path3DHST+cat_sfr,names=True,comments='#')
    dat_zbest   = pyfits.open(path3DHST+cat_zbest)[1].data
    dat_zinfo   = pyfits.open(path3DHST+cat_zinfo)[1].data
    dat_ELflux  = pyfits.open(path3DHST+cat_ELflux)[1].data

    if printcolnames:
        print   '   Loaded the following catalogs on columns'
        #print   '-------------',cat_phot   ,'-------------\n',dat_phot.dtype.names
        print   '-------------',cat_dupspec,'-------------\n',dat_dupspec.dtype.names
        print   '-------------',cat_fast   ,'-------------\n',dat_fast.dtype.names
        #print   '-------------',cat_rf     ,'-------------\n',dat_rf.dtype.names
        print   '-------------',cat_sfr    ,'-------------\n',dat_sfr.dtype.names
        print   '-------------',cat_zbest  ,'-------------\n',dat_zbest.dtype.names
        print   '-------------',cat_zinfo  ,'-------------\n',dat_zinfo.dtype.names
        print   '-------------',cat_ELflux ,'-------------\n',dat_ELflux.dtype.names

    collist = cm3.objectsinfo_columnlist()

    Ncols   = len([cc.split('___')[1] for cc in collist])
    Ncols_u = len(np.unique(np.asarray([cc.split('___')[1] for cc in collist])))
    if Ncols != Ncols_u:
        sys.exit(' The columns to be extracted from the 3D-HST catalos have multiple entries; fix...')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Setting up output ascii and fits files to contain object info'
    asciifile = output_basename+'.txt'
    fitsfile  = output_basename+'.fits'

    if os.path.isfile(asciifile) and (clobber == False):
        sys.exit(asciifile+' exists and clobber=False so no file generated --> ABORTING')
    else:
        if verbose: print ' - Generating '+asciifile
        fout = open(asciifile,'w')
        fout.write('# Output from CDFS_MUSEvs3DHST.extract_objectsinfo() produced on '+kbs.DandTstr2()+' \n')
        columnsting = '  '.join([cc.split('___')[1] for cc in collist])
        fout.write('# id_MUSE ra_MUSE dec_MUSE id_3DHST ra_3DHST dec_3DHST r_match_arcsec '+columnsting+' specinfo \n')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Nobj = len(MUSEidlist)
    if verbose: print ' - Looping over the '+str(Nobj)+' objects and extracting information'
    for MUSEid in MUSEidlist:
        objent = np.where(cmdat['ID'] == MUSEid)[0]
        if len(objent) == 0:
            if verbose: print '   No match top MUSEid = '+str(MUSEid)
        else:
            if cmdat['r_match_arcsec'][objent] > crossmatchtol:
                if verbose: print '   r_match ('+str(cmdat['r_match_arcsec'][objent])+\
                                  ') above crossmatch tolerance for MUSEid = '+str(MUSEid)
            else:
                id_MUSE   = cmdat['id'][objent]
                if id_MUSE != MUSEid: print ' ERROR ids dont mathc --> stopping ';pdb.set_trace()
                ra_MUSE   = cmdat['ra'][objent]
                dec_MUSE  = cmdat['dec'][objent]
                id_3DHST  = cmdat['id_match'][objent]
                ra_3DHST  = cmdat['ra_match'][objent]
                dec_3DHST = cmdat['dec_match'][objent]
                r_match   = cmdat['r_match_arcsec'][objent]
                objstr  = str("%20i"   % id_MUSE)    +'  '+\
                          str("%16.8f" % ra_MUSE)    +'  '+\
                          str("%16.8f" % dec_MUSE)   +'  '+\
                          str("%20i"   % id_3DHST)   +'  '+\
                          str("%16.8f" % ra_3DHST)   +'  '+\
                          str("%16.8f" % dec_3DHST)  +'  '+\
                          str("%16.8f" % r_match)
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                for col in collist:
                    if 'fast___' in col:
                        fastent = np.where(dat_fast['id'] == id_3DHST)[0]
                        objstr  = objstr+' '+str("%e" % dat_fast[col.split('___')[1]][fastent][0])
                    if 'sfr___' in col:
                        sfrent = np.where(dat_sfr['id'] == id_3DHST)[0]
                        objstr  = objstr+' '+str("%e" % dat_sfr[col.split('___')[1]][sfrent][0])
                    if 'zbest___' in col:
                        zbestent = np.where(dat_zbest['phot_id'] == id_3DHST)[0]
                        objstr  = objstr+' '+str("%e" % dat_zbest[col.split('___')[1]][zbestent][0])
                    if 'zinfo___' in col:
                        zinfoent = np.where(dat_zinfo['phot_id'] == id_3DHST)[0]

                        if 'grism_id' in col:
                            objstr  = objstr+' '+str("%30s" % dat_zinfo[col.split('___')[1]][zinfoent][0])
                        else:
                            objstr  = objstr+' '+str("%e" % dat_zinfo[col.split('___')[1]][zinfoent][0])
                    if 'ELflux___' in col:
                        ELfluxent = np.where(dat_ELflux['number'] == id_3DHST)[0]

                        if 'grism_id' in col:
                            objstr  = objstr+' '+str("%30s" % dat_ELflux[col.split('___')[1]][ELfluxent][0])
                        else:
                            objstr  = objstr+' '+str("%e" % dat_ELflux[col.split('___')[1]][ELfluxent][0])
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                dupspecent = np.where(dat_dupspec['id'] == id_3DHST)[0]
                specinfo   = 'xxx'.join(dat_dupspec[dupspecent][0].tolist()[1:]).replace('00000','None')
                objstr     = objstr+'  '+specinfo
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                fout.write(objstr+'\n')
        fout.close()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #if os.path.isfile(fitsfile) and (clobber == False):
    #    if verbose: print ' - WARNING: '+fitsfile+' exists and clobber=False so no file generated'
    #else:
    #    if verbose: print ' - Generating (coverting ascii output file to) '+fitsfile
    #    fitspath   = kbs.pathAname(asciifile)[0]
    #    outputfile = f2a.ascii2fits(asciifile,asciinames=True,skip_header=1,outpath=fitspath,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def objectsinfo_columnlist():
    """
    List with column and catalog information for extractions

    See /Users/kschmidt/work/MUSE/candelsCDFS_3DHST/3DHSTcolumninfo.txt for column name list

    """
    columnlist = ['fast___ltau',
    'fast___metal',
    'fast___lage',
    'fast___Av',
    'fast___lmass',
    'fast___lsfr',
    'fast___lssfr',
    'sfr___sfr',
    'sfr___sfr_IR',
    'sfr___sfr_UV',
    'sfr___beta',
    'zbest___z_best_s',
    'zbest___use_phot',
    'zbest___use_zgrism',
    'zbest___z_best',
    'zbest___z_best_l95',
    'zbest___z_best_l68',
    'zbest___z_best_u68',
    'zbest___z_best_u95',
    'zinfo___phot_id',
    'zinfo___grism_id',
    'zinfo___jh_mag',
    'zinfo___z_spec',
#    'ELflux___grism_id',
    'ELflux___Lya_FLUX',
    'ELflux___Lya_FLUX_ERR',
    'ELflux___Lya_SCALE',
    'ELflux___Lya_EQW',
    'ELflux___Lya_EQW_ERR',
    'ELflux___CIV_FLUX',
    'ELflux___CIV_FLUX_ERR',
    'ELflux___CIV_SCALE',
    'ELflux___CIV_EQW',
    'ELflux___CIV_EQW_ERR',
    'ELflux___MgII_FLUX',
    'ELflux___MgII_FLUX_ERR',
    'ELflux___MgII_SCALE',
    'ELflux___MgII_EQW',
    'ELflux___MgII_EQW_ERR',
    'ELflux___OII_FLUX',
    'ELflux___OII_FLUX_ERR',
    'ELflux___OII_SCALE',
    'ELflux___OII_EQW',
    'ELflux___OII_EQW_ERR',
    'ELflux___OIII_FLUX',
    'ELflux___OIII_FLUX_ERR',
    'ELflux___OIII_SCALE',
    'ELflux___OIII_EQW',
    'ELflux___OIII_EQW_ERR',
    'ELflux___Ha_FLUX',
    'ELflux___Ha_FLUX_ERR',
    'ELflux___Ha_SCALE',
    'ELflux___Ha_EQW',
    'ELflux___Ha_EQW_ERR']
    return columnlist
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =