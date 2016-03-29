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
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
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
    if os.path.isfile(fitsfile) and (clobber == False):
        if verbose: print ' - WARNING: '+fitsfile+' exists and clobber=False so no file generated'
    else:
        if verbose: print ' - Generating (coverting ascii output file to) '+fitsfile
        fitspath   = kbs.pathAname(asciifile)[0]
        asciidat   = np.genfromtxt(asciifile,names=True,skip_header=1,comments='#')
        keys       = asciidat.dtype.names
        fitsformat = np.asarray(['D']*(len(keys)),dtype='5a')
        fitsformat[np.asarray(keys) == 'id_MUSE'] = 'J' #32 bit integer
        fitsformat[np.asarray(keys) == 'id_3DHST'] = 'J' #32 bit integer
        fitsformat[np.asarray(keys) == 'specinfo'] = '200A'
        fitsformat[np.asarray(keys) == 'grism_id'] = '25A'
        outputfile = f2a.ascii2fits(asciifile,asciinames=True,skip_header=1,outpath=fitspath,
                                    fitsformat=fitsformat,verbose=verbose)

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
def objinfo_CDFSobj(MUSEcat,zcut=[2.9-10.0],output_basename='./objectinfo_output',cmfile=False,crossmatchtol=0.1
                    ,clobber=False,MUSEcatext=1,verbose=True):
    """

    Make cut on MUSE redshift and extract 3D-HST info

    --- INPUT ---
    MUSEcat          MUSE catalog to extract objects from
    output_basename  Basename to use for output files
    cmfile           Fits file with crossmatch to 3D-HST catalogs. If None the crossmatch will be performed
    crossmatchtol    Tolerance for crossmatch. Deafult is to only extract information of r_match <0.1 arc sec
    clobber          Overwrite files?
    MUSEcatext       Extension containing data in MUSEcat
    verbose          Toggle verbosity

    --- EXAMPLE OF USE ---
    import CDFS_MUSEvs3DHST as cm3
    MUSEcat = '/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_candels-cdfs_v0.2.fits'
    cm3.objinfo_CDFSobj(MUSEcat,zcut=[2.9,10.0],output_basename='./MUSECDFS_z2p9-10p0_cmtol1p0',crossmatchtol=1.0,clobber=True)

    cm3.objinfo_CDFSobj(MUSEcat,zcut=[2.9,10.0],output_basename='./MUSECDFS_z2p9-10p0_cmtol0p5',crossmatchtol=0.5,clobber=True)


    cm3.objinfo_CDFSobj(MUSEcat,zcut=[0.0,10.0],output_basename='./MUSECDFS_z0p0-10p0_cmtol1p0',crossmatchtol=1.0,clobber=True)

    cm3.objinfo_CDFSobj(MUSEcat,zcut=[0.0,10.0],output_basename='./MUSECDFS_z0p0-10p0_cmtol0p5',crossmatchtol=0.5,clobber=True)


    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if not cmfile:
        if verbose: print ' - No crossmatch output provide, so crossmatching input catalog to 3D_HST goodss catalog'
        cmfile = output_basename+'_crossmatch_3DHST'
        radeccat = MUSEcat
        matchcat = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
        cmout = kbs.crossmatch2cat(radeccat,matchcat=matchcat,clobber=clobber,
                                   writetofile=cmfile.split('.fit')[0],verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading MUSE catalog and extracting ids based on selection'
    musedat   = pyfits.open(MUSEcat)[MUSEcatext].data
    goodmask  = (musedat['redshift'] >= zcut[0]) & (musedat['redshift'] <= zcut[1])
    gooddat   = musedat[goodmask]
    goodids   = musedat['ID'][goodmask]
    Ngoodobj  = len(goodids)
    if verbose: print '   Found '+str(Ngoodobj)+' objects with z = ',zcut
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if Ngoodobj > 0:
        subcatname = output_basename+'_MUSEsubcat.fits'
        if verbose: print ' - Saving MUSE sub-catalog to '+subcatname
        pyfits.writeto(subcatname, gooddat, header=pyfits.open(MUSEcat)[MUSEcatext].header, clobber=clobber)

        if verbose: print ' - Extracting 3D-HST data for the objects'
        cm3.extract_objectsinfo(goodids,printcolnames=False,clobber=clobber,verbose=verbose,
                                crossmatchtol=crossmatchtol,output_basename=output_basename+'_3DHSTinfo')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def extract_3DHSTdata(infofits_3DHST,outputdir=False,verbose=True):
    """
    Extract list of 3DHST grism spectra for objects in 3DHST info cat and (if requested)
    copy data to seperate directory for further inspection.

    --- INFO ---

    infofits_3DHST   Fits file with 3DHST info containing column specinfo to extract file names from
                     (generated with extract_objectsinfo)
    outputdir        If data is to be copied to directory, provide directory path here. If False,
                     the filenames and IDs will just be returned
    verbose          Toggle verbosity

    --- EXAMPLE OF USE ---

    idM_out, id3_out, spec_out = cm3.extract_3DHSTdata('MUSECDFS_z2p9-10p0_cmtol0p5_3DHSTinfo.fits')

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading data'
    infodat   = pyfits.open(infofits_3DHST)[1].data
    id_MUSE   = infodat['id_MUSE']
    id_3DHST  = infodat['id_3DHST']
    specinfo  = infodat['specinfo']
    Nobj      = len(id_MUSE)
    Nspec     = len((' '.join(specinfo.tolist()).replace('xxx',' ').replace('None',' ')).split())
    if verbose: print '   Found a total of '+str(Nspec)+' for the '+str(Nobj)+' objects in \n   '+infofits_3DHST
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Putting output together'
    idM_out   = []
    id3_out   = []
    spec_out  = []
    for ii, idM in enumerate(id_MUSE):
        id3   = id_3DHST[ii]
        specs = specinfo[ii].replace('xxx',' ').replace('None',' ').split()
        for ss in specs:
            idM_out.append(idM)
            id3_out.append(id3)
            spec_out.append(ss)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if outputdir:
        if verbose: print ' - Copying data for objects to '+outputdir
        print 'ERROR - not enabled... stopping'
        pdb.set_trace()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return idM_out, id3_out, spec_out
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_fits2Dfilenames(specinfo,dataparentdir,verbose=True):
    """

    glob for and return absolute paths for existing fits 2D files

    """
    filelist = []
    namelist = specinfo.split('xxx')[0]

    for name in namelist:
        if name != 'None':
            filename = glob.glob(dataparentdir+'/'+name[:-11]+'/2D/FITS/'+name+'.2D.fits')
            Nfiles   = len(filename)
            if len(filename) == 1:
                filelist.append(os.path.abspath(filename[0]))
            elif Nfiles == 0:
                if verbose: print '   No 2D.fits files found in data sub-directories '
            else:
                if verbose: print '   NB! More than one 2D.fits file found; returning path for first file '
                filelist.append(os.path.abspath(filename[0]))

    return filelist
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_objectinfo_all(catMUSE,cat3DHST,outputdir='./',zrange=[0.0,10.0],dzcut=0.25,verbose=True):
    """
    run all plotting for a set of MUSE and 3DHST info catalogs

    --- EXAMPLE OF USE ---
    catMUSE  = 'MUSECDFS_z0p0-10p0_cmtol1p0_MUSEsubcat.fits'
    cat3DHST = 'MUSECDFS_z0p0-10p0_cmtol1p0_3DHSTinfo.fits'
    cm3.plot_objectinfo_all(catMUSE,cat3DHST,outputdir='./',verbose=True)

    """

    cm3.plot_objectinfo_zVSz(catMUSE,cat3DHST,outputdir=outputdir,zrange=zrange,verbose=verbose)
    cm3.plot_objectinfo_zVSline(catMUSE,cat3DHST,outputdir=outputdir,plotline='OII',zrange=zrange,
                                dzcut=dzcut,verbose=verbose)
    cm3.plot_objectinfo_zVSline(catMUSE,cat3DHST,outputdir=outputdir,plotline='OIII',zrange=zrange,
                                dzcut=dzcut,verbose=verbose)
    cm3.plot_objectinfo_zVSline(catMUSE,cat3DHST,outputdir=outputdir,plotline='MgII',zrange=zrange,
                                dzcut=dzcut,verbose=verbose)
    cm3.plot_objectinfo_zVSline(catMUSE,cat3DHST,outputdir=outputdir,plotline='CIV',zrange=zrange,
                                dzcut=dzcut,verbose=verbose)
    cm3.plot_objectinfo_zVSline(catMUSE,cat3DHST,outputdir=outputdir,plotline='Ha',zrange=zrange,
                                dzcut=dzcut,verbose=verbose)
    cm3.plot_objectinfo_zVSline(catMUSE,cat3DHST,outputdir=outputdir,plotline='Lya',zrange=zrange,
                                dzcut=dzcut,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_objectinfo_zVSz(catMUSE,cat3DHST,outputdir='./',zrange=[0.0,10.0],verbose=True):
    """
    Plotting content of a MUSE catalog, with corresponding 3DHST info catalog

    --- INPUT ---

    catMUSE       MUSE catalog to plot
    cat3DHST      3DHST info catalog corresponding to MUSE catalog
    outputdir     Directory to save plots in
    verbose       Toggle verbosity

    --- EXAMPLE OF USE ---

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading data in MUSE and 3DHST catalogs'
    dat_MUSE  = pyfits.open(catMUSE)[1].data
    dat_3DHST = pyfits.open(cat3DHST)[1].data

    id_MUSEcat     = dat_MUSE['ID']
    id_3DHSTcat    = dat_3DHST['id_MUSE']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = outputdir+'/'+catMUSE.split('MUSEsubcat.fit')[0]+'zMUSEvsz3DHSTBEST.pdf'
    if verbose: print ' - Setting up and generating plot'
    fig = plt.figure(figsize=(6, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.95, bottom=0.10, top=0.95)
    Fsize    = 10
    lthick   = 1
    marksize = 3
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('/')[-1].replace('_','\_'),fontsize=Fsize)

    for oo, obj in enumerate(id_MUSEcat):
        objent = np.where(id_3DHSTcat == obj)[0]

        zMUSE          = dat_MUSE['redshift'][oo]
        zMUSE_conf     = dat_MUSE['confidence'][oo]
        zMUSE_conf     = dat_MUSE['line_id'][oo]

        z3DHST         = dat_3DHST['z_best'][objent]
        z3DHST_l95     = dat_3DHST['z_best_l95'][objent]
        z3DHST_l68     = dat_3DHST['z_best_l68'][objent]
        z3DHST_u95     = dat_3DHST['z_best_u95'][objent]
        z3DHST_u68     = dat_3DHST['z_best_u68'][objent]
        z3DHST_source  = dat_3DHST['z_best_s'][objent] # 1=ground-based spec; 2=grism; 3=photometry; 0=star

        dz = np.abs(zMUSE - z3DHST)

        if dz < 0.25:
            yerrvals = None
        else:
            yerrvals = None #[z3DHST_l68,z3DHST_u68]


        if z3DHST_source == 1:
            plt.errorbar(zMUSE,z3DHST,xerr=None,yerr=None,
                         fmt='o',lw=lthick,ecolor='white', markersize=marksize,markerfacecolor='white',
                         markeredgecolor = 'k')
        elif z3DHST_source == 2:
            plt.errorbar(zMUSE,z3DHST,xerr=None,yerr=yerrvals,
                         fmt='o',lw=lthick,ecolor='red', markersize=marksize,markerfacecolor='red',
                         markeredgecolor = 'k')
        elif z3DHST_source == 3:
            plt.errorbar(zMUSE,z3DHST,xerr=None,yerr=yerrvals,
                         fmt='o',lw=lthick,ecolor='blue', markersize=marksize,markerfacecolor='blue',
                         markeredgecolor = 'k')
        elif z3DHST_source == 0:
            plt.errorbar(zMUSE,z3DHST,xerr=None,yerr=None,
                         fmt='*',lw=lthick,ecolor='orange', markersize=marksize,markerfacecolor='orange',
                         markeredgecolor = 'k')

    plt.xlabel('zMUSE', fontsize=Fsize)
    plt.ylabel('z3DHST', fontsize=Fsize)

    plt.xlim(zrange)
    plt.ylim(zrange)


    #--------- LEGEND ---------
    plt.errorbar(-1,-1,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='white', markersize=marksize*2,
                 markerfacecolor='white',markeredgecolor = 'k',label='Ground-based spec')
    plt.errorbar(-1,-1,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='red', markersize=marksize*2,
                 markerfacecolor='red',markeredgecolor = 'k',label='G141 Grism')
    plt.errorbar(-1,-1,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='blue', markersize=marksize*2,
                 markerfacecolor='blue',markeredgecolor = 'k',label='EA$z$Y Photo')
    plt.errorbar(-1,-1,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='orange', markersize=marksize*2,
                 markerfacecolor='orange',markeredgecolor = 'k',label='Star')


    leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize},ncol=1,numpoints=1)
                     #bbox_to_anchor=(1.25, 1.03))  # add the legend
    leg.get_frame().set_alpha(0.7)
    #--------------------------

    if verbose: print '   Saving plot to',plotname
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_objectinfo_zVSline(catMUSE,cat3DHST,outputdir='./',plotline='OII',zrange=[0.0,10.0],
                            linerange=False,EWrange=False,dzcut=0.25,verbose=True):
    """
    Plotting content of a MUSE catalog, with corresponding 3DHST info catalog

    --- INPUT ---

    catMUSE       MUSE catalog to plot
    cat3DHST      3DHST info catalog corresponding to MUSE catalog
    outputdir     Directory to save plots in
    verbose       Toggle verbosity

    --- EXAMPLE OF USE ---

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading data in MUSE and 3DHST catalogs'
    dat_MUSE  = pyfits.open(catMUSE)[1].data
    dat_3DHST = pyfits.open(cat3DHST)[1].data

    id_MUSEcat     = dat_MUSE['ID']
    id_3DHSTcat    = dat_3DHST['id_MUSE']
    lineflux       = dat_3DHST[plotline+'_FLUX']
    lineEW         = dat_3DHST[plotline+'_EQW']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = outputdir+'/'+catMUSE.split('MUSEsubcat.fit')[0]+'zMUSEvs'+plotline+'_Flux.pdf'
    if verbose: print ' - Setting up and generating plot'
    fig = plt.figure(figsize=(6, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.95, bottom=0.10, top=0.95)
    Fsize    = 10
    lthick   = 1
    marksize = 6
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('/')[-1].replace('_','\_'),fontsize=Fsize)

    if not linerange:
        try:
            linerange = [np.min(lineflux[lineflux > 0])*0.9,np.max(lineflux[lineflux > 0])*1.1]
        except:
            linerange = [0,10]

    for oo, obj in enumerate(id_MUSEcat):
        objent = np.where(id_3DHSTcat == obj)[0]

        zMUSE          = dat_MUSE['redshift'][oo]
        zMUSE_conf     = dat_MUSE['confidence'][oo]
        zMUSE_conf     = dat_MUSE['line_id'][oo]

        z3DHST         = dat_3DHST['z_best'][objent]
        f_3DHST        = dat_3DHST[plotline+'_FLUX'][objent]
        ferr_3DHST     = dat_3DHST[plotline+'_FLUX_ERR'][objent]

        dz = np.abs(zMUSE - z3DHST)


        if dz < dzcut:
            col = 'green'
        else:
            col = 'red'

        if f_3DHST > 0:
            plt.errorbar(zMUSE,f_3DHST,xerr=None,yerr=ferr_3DHST,
                         fmt='o',lw=lthick,ecolor=col, markersize=marksize,markerfacecolor=col,
                         markeredgecolor = 'k')#,label='Ground-based spec')

    plt.xlabel('zMUSE', fontsize=Fsize)
    plt.ylabel('Flux '+plotline+' [$10^{-17}$ ergs/s/cm$^2$]', fontsize=Fsize)

    plt.xlim(zrange)
    plt.ylim(linerange)

    #--------- LEGEND ---------
    plt.errorbar(-1,-1,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='red', markersize=marksize,
                 markerfacecolor='red',markeredgecolor = 'k',label='$|z$MUSE-$z$3DHST$| >$ '+str("%.2f" % dzcut))
    plt.errorbar(-1,-1,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='green', markersize=marksize,
                 markerfacecolor='green',markeredgecolor = 'k',label='$|z$MUSE-$z$3DHST$| <$ '+str("%.2f" % dzcut))

    leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize},ncol=1,numpoints=1)
                     #bbox_to_anchor=(1.25, 1.03))  # add the legend
    leg.get_frame().set_alpha(0.7)
    #--------------------------

    if verbose: print '   Saving plot to',plotname
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = outputdir+'/'+catMUSE.split('MUSEsubcat.fit')[0]+'zMUSEvs'+plotline+'_EW.pdf'
    if verbose: print ' - Setting up and generating plot'
    fig = plt.figure(figsize=(6, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.95, bottom=0.10, top=0.95)
    Fsize    = 10
    lthick   = 1
    marksize = 6
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('/')[-1].replace('_','\_'),fontsize=Fsize)

    if not EWrange:
        try:
            lineEW[lineEW > 10**4] = -99999
            EWrange = [np.min(lineEW[lineEW > 0])*0.9,np.max(lineEW[lineEW > 0])*1.1]
        except:
            EWrange = [0,200]

    for oo, obj in enumerate(id_MUSEcat):
        objent = np.where(id_3DHSTcat == obj)[0]

        zMUSE          = dat_MUSE['redshift'][oo]
        zMUSE_conf     = dat_MUSE['confidence'][oo]
        zMUSE_conf     = dat_MUSE['line_id'][oo]

        z3DHST         = dat_3DHST['z_best'][objent]
        EW_3DHST       = dat_3DHST[plotline+'_EQW'][objent]
        EWerr_3DHST    = dat_3DHST[plotline+'_EQW_ERR'][objent]

        dz = np.abs(zMUSE - z3DHST)


        if dz < dzcut:
            col = 'green'
        else:
            col = 'red'

        if EW_3DHST > 0:
            plt.errorbar(zMUSE,EW_3DHST,xerr=None,yerr=EWerr_3DHST,
                         fmt='o',lw=lthick,ecolor=col, markersize=marksize,markerfacecolor=col,
                         markeredgecolor = 'k')

    plt.xlabel('zMUSE', fontsize=Fsize)
    plt.ylabel('EW '+plotline+' [\AA]', fontsize=Fsize)

    plt.xlim(zrange)
    plt.ylim(EWrange)

    #--------- LEGEND ---------
    plt.errorbar(-1,-1,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='red', markersize=marksize,
                 markerfacecolor='red',markeredgecolor = 'k',label='$|z$MUSE-$z$3DHST$| >$ '+str("%.2f" % dzcut))
    plt.errorbar(-1,-1,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='green', markersize=marksize,
                 markerfacecolor='green',markeredgecolor = 'k',label='$|z$MUSE-$z$3DHST$| <$ '+str("%.2f" % dzcut))

    leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize},ncol=1,numpoints=1)
                     #bbox_to_anchor=(1.25, 1.03))  # add the legend
    leg.get_frame().set_alpha(0.7)
    #--------------------------

    if verbose: print '   Saving plot to',plotname
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_3dhst_SpecAndInfo(catMUSE,cat3DHST,plotID='all',outputdir='./',plotscale=[1,99],figxsize=13,cmap='gray',
                           left=0.1, bottom=0.07, right=0.02, top=0.03, wspace=0.02, hspace=0.04,
                           dataparentdir='/Volumes/DATABCKUP3/GOODSdata/3D-HST/GOODSS_WFC3_V4.1.5/',
                           catphoto='/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog'
                                    '/goodss_3dhst.v4.1.cat.FITS',verbose=True):
    """
    Plotting the available 3D-HST spectra and the 3D-HST (and MUSE) information
    for all (or a few) objects in an info-cat generated with, e.g., extract_objectsinfo()
    or objinfo_CDFSobj()

    --- INPUT ---

    cat3DHST        Output fits catalog from extract_objectsinfo() extraction of 3D-HST info
    catMUSE         MUSE fits catalog the 3D-HST info extraction was based on
    plotID          List of MUSE IDs to plot. If 'all' then all objects found in catMUSE are plotted.
    dataparentdir   Parent dir containing the 3D-HST spectra, e.g., GOODSS_WFC3_V4.1.5, which
                    contains subdirectories of data from individual pointings in 1D, 2D, BIG,
                    etc. sub-directories.
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---

    catMUSE  = 'MUSECDFS_z0p0-10p0_cmtol1p0_MUSEsubcat.fits'
    cat3DHST = 'MUSECDFS_z0p0-10p0_cmtol1p0_3DHSTinfo.fits'
    cm3.plot_3dhst_SpecAndInfo(catMUSE,cat3DHST,plotID='all',verbose=True,cmap='coolwarm')

    cm3.plot_3dhst_SpecAndInfo(catMUSE,cat3DHST,plotID=[10611042],verbose=True)

    cm3.plot_3dhst_SpecAndInfo(catMUSE,cat3DHST,plotID=[12341163],verbose=True)

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading data in MUSE and 3DHST catalogs'
    dat_MUSE  = pyfits.open(catMUSE)[1].data
    dat_3DHST = pyfits.open(cat3DHST)[1].data
    dat_photo = pyfits.open(catphoto)[1].data

    id_MUSEcat     = dat_MUSE['ID']
    id_3DHSTcat    = dat_3DHST['id_MUSE']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plotID == 'all':
        plotIDs = np.asarray(id_MUSEcat)
    else:
        plotIDs = np.asarray(plotID)
    Nplots = len(plotIDs)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Generating plots for '+str(Nplots)+' objects to plot'
    for oo in xrange(Nplots):
        skipplot   = False            # keep track of whether to generate plot or not
        objid      = plotIDs[oo]
        if verbose: print ' - Grabbing data for ID = '+str(objid)
        ent_muse   = np.where(id_MUSEcat == objid)[0]
        if len(ent_muse) != 1:
            if verbose: print '   WARNING - no match in MUSE catalog; no plot generated'
            skipplot = True

        ent_3dhst  = np.where(id_3DHSTcat == objid)[0]

        if len(ent_3dhst) != 1:
            if verbose: print '   WARNING - no match in 3DHST catalog; no plot generated'
            skipplot = True

        if not skipplot:
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if verbose: print '   Getting 2D fits file names'
            objdat_MUSE   = dat_MUSE[ent_muse]
            objdat_3DHST  = dat_3DHST[ent_3dhst]

            try:
                fits2Dfilenames = cm3.get_fits2Dfilenames(objdat_3DHST['specinfo'],dataparentdir,verbose=verbose)
            except:
                pdb.set_trace()
            N2Dfitsfiles    = len(fits2Dfilenames)
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if N2Dfitsfiles == 0:
                if verbose: print '   WARNING - no 3DHST grism spectra stored in specinfo (all "None" and "xxx")'
                continue
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            plotname = outputdir+'/'+str(objid)+'_3dhstSpecAndInfo.pdf'
            if verbose: print '   Setting up and generating plot'

            #----------------- Plot setup -----------------
            Fsize  = 10
            lwidth = 1
            Tcolor = 'black'
            xscale = 1e4
            g141range = [11100,16500]
            xvalrange = [1.10,1.67]
            plt.rc('text', usetex=True)                         # enabling LaTex rendering of text
            plt.rc('font', family='serif',size=Fsize)           # setting text font
            plt.rc('xtick', labelsize=Fsize)
            plt.rc('ytick', labelsize=Fsize)

            xs     = figxsize
            sh     = np.shape(pyfits.open(fits2Dfilenames[0])['SCI'].data)
            aspect = (N2Dfitsfiles+2)*sh[0]*1./(3*sh[1])
            lrbt   = np.array([left,right,bottom,top])*5./xs
            ys     = (1-lrbt[1]-lrbt[0])/(1-lrbt[3]-lrbt[2])*xs*aspect
            lrbt[[2,3]] /= aspect

            fig = Figure(figsize=(xs,ys))
            fig.subplots_adjust(left=lrbt[0], bottom=lrbt[2], right=1-lrbt[1], top=1-lrbt[3],
                                wspace=wspace, hspace=hspace)
            #fig.subplots_adjust(wspace=0.2, hspace=0.2,left=0.04, right=0.98, bottom=0.13, top=0.89)
            #plt.title(plotname.split('/')[-1].replace('_','\_'),fontsize=Fsize)

            #----------------- Load Data -----------------
            plotcounter = 0
            for row, filetwod in enumerate(fits2Dfilenames):
                spec_flux          = pyfits.open(filetwod)['SCI'].data
                spec_contam        = pyfits.open(filetwod)['CONTAM'].data
                spec_var           = (pyfits.open(filetwod)['WHT'].data)**2
                spec_flux_nocontam = spec_flux-spec_contam
                wave               = pyfits.open(filetwod)['WAVE'].data

                # - - - -  setting plotting scale - - - -
                if verbose: print '   Setting scales for plotting to',
                # diff = spec_flux_nocontam
                # if len(plotscale) == 0:
                #     plotscale = [plotscale, plotscale]
                #
                # ok = (spec_var > 0) & (spec_flux != 0) & (spec_flux < 0.01)
                # if ok.sum() > 1:
                #     vmin = np.percentile(diff[ok].flatten(), plotscale[0])
                #     vmax = np.percentile(diff[ok].flatten(), plotscale[1])
                # else:
                #     vmin, vmax = -1,1

                scalearr   = spec_flux_nocontam[10:-10,10:-10]
                scalearr   = scalearr[scalearr != 0.0]
                if len(scalearr) > 0:
                    vmaximg    = np.abs(np.max(scalearr)*0.7)
                    vmin, vmax = vmaximg*1e-10, vmaximg
                else:
                    vmin, vmax = 0, 1

                if verbose: print vmin,vmax

                # - - - -  setting xrange and interpolate - - - -
                xx = np.arange(spec_flux.shape[1])
                if 'G141' in filetwod:
                    ltick = np.array([1.1,1.2,1.3,1.4,1.5,1.6,1.7])
                elif 'G102' in filetwod:
                    ltick = np.array([0.8,0.9,1.0,1.1])
                elif 'G800L' in filetwod:
                    ltick = np.array([0.6,0.7,0.8,0.9,1.0])
                else:
                    sys.exit('Grism name (for xticks) not found in filename '+filetwod+' --> using G141')

                scale = 1.e4
                xtick = np.interp(ltick*scale, wave, xx)
                # - - - -  display images - - - -
                imglist  = [spec_flux,spec_contam,spec_flux_nocontam]
                textlist = ['[SCI]','[CONTAM]','[SCI-CONTAM]']
                for colno in [0,1,2]:
                    plotcounter = plotcounter+1
                    ax = fig.add_subplot(N2Dfitsfiles,3,plotcounter) #row,column,plotnumber
                    ax.imshow(imglist[colno], interpolation='nearest', vmin=vmin, vmax=vmax, cmap=cmap)
                    ax.set_xticks(xtick)
                    ax.set_xticklabels(ltick)
                    if plotcounter > int(N2Dfitsfiles*3-N2Dfitsfiles):
                        ax.set_xlabel(r'$\lambda / \mu\mathrm{m}$')

                    if colno == 0:
                        ax.set_ylabel(r' r [arc sec]')
                    else:
                        ax.set_yticklabels([])

                    ax.text(0.95, 0.94, textlist[colno],
                            backgroundcolor='w',ha='right', va='top',
                            transform=ax.transAxes, color=Tcolor, size=Fsize)

                    if colno == 1:
                        ax.text(0.5, 0.94, filetwod.split('/')[-1].replace('_','\_'),#+textlist[colno],
                                backgroundcolor='w',ha='center', va='top',
                                transform=ax.transAxes, color=Tcolor, size=Fsize)

            #----------------- Position Legend -----------------
            # box = ax.get_position()
            # ax.set_position([box.x0, box.y0, box.width, box.height * 0.83])
            # ax.plot(0,0,'orange',label='G102 PA=XX',linewidth=lwidth*2)
            # ax.plot(0,0,'red',label='G141 PA=XX',linewidth=lwidth*2)
            # ax.plot(0,0,'blue',label='Photo-$z$',linewidth=lwidth*2)
            # ax.plot(0,0,'k',label='$z$fit model',linewidth=lwidth)
            # leg = ax.legend(fancybox=True, loc='upper left',numpoints=1,prop={'size':Fsize-3.},
            #                 ncol=3,bbox_to_anchor=(-0.08, 1.27))
            # leg.get_frame().set_alpha(0.7)

            #----------------- Title with info -----------------
            id3dhst   = dat_3DHST['ID_3DHST'][ent_3dhst][0]
            titletext = 'idMUSE='+str(objid)+', id3dhst='+str(id3dhst)
            titletext = titletext+', rmatch='+str("%.2f" % dat_3DHST['r_match_arcsec'][ent_3dhst])+"''"
            titletext = titletext+', log(M/M$\odot$)='+str("%.2f" % dat_3DHST['lmass'][ent_3dhst][0])
            titletext = titletext+', SFR/[M$\odot$/yr]='+str("%.2f" % (10**dat_3DHST['lsfr'][ent_3dhst][0]))
            titletext = titletext+', zbest='+str("%.4f" % dat_3DHST['z_best'][ent_3dhst][0])

            fOIII     = dat_3DHST['OIII_FLUX'][ent_3dhst][0]
            eOIII     = dat_3DHST['OIII_FLUX_ERR'][ent_3dhst][0]
            if fOIII != -99:
                titletext = titletext+', f[OIII]='+str("%.2f" % fOIII)+'$\pm$'+str("%.2f" % eOIII)
                titletext = titletext+', SNR[OIII]='+str("%.2f" % (fOIII/eOIII))
            else:
                titletext = titletext+', f[OIII]=None, SNR[OIII]=None'


            # - - - - get photo info - - - -
            photid    = dat_photo['id'][id3dhst-1]
            f160      = dat_photo['f_f160w'][id3dhst-1]
            if f160 != -99:
                mf160w    = -2.5*np.log10(f160)+25.0
                ef160w    = (2.5/np.log(10)) * dat_photo['e_f160w'][id3dhst-1]/f160
                titletext = titletext+', m160='+str("%.2f" % mf160w)+'$\pm$'+str("%.2f" % ef160w)
            else:
                titletext = titletext+', m160=None'


            #titletext = titletext+', photid='+str(photid)

            #ax.text(0.0,1.1, titletext,transform=ax.transAxes, horizontalalignment='center')
            fig.suptitle(titletext, fontsize=Fsize)
            #----------------- Write plot to file -----------------
            if verbose: print '   Saving plot to',plotname
            canvas = FigureCanvasAgg(fig)
            canvas.print_figure(plotname, dpi=300, transparent=False)

            # plt.savefig(plotname)
            # plt.clf()
            # plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =