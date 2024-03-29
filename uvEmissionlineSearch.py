# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Scripts, functions and routines to (enable) search for UV emission lines in MUSE data
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import os
import sys
import glob
import subprocess
import collections


import scipy
import scipy.stats as ss


import pyregion
import pyfits as pyfitsOLD
import datetime
import numpy as np
from uncertainties import unumpy
import shutil
import time

import astropy
import astropy.io.fits as afits
from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.cosmology as acosmo
import astropy.coordinates as acoord
import astropy.units as u
from astropy.convolution import Gaussian2DKernel, convolve_fft

import MiGs
import fits2ascii as f2a
import MUSEWideUtilities as mu
import MUSEWidePlots as mwp
import kbsutilities as kbs
import tdose_utilities as tu
import uvEmissionlineSearch as uves
import NEOGALmodels as nm
import photoionizationPDFs as pp
#import rxj2248_BooneBalestraSource as bbs
import felis_build_template as fbt
import felis
import literaturecollection_emissionlinestrengths as lce
import stacking
from itertools import combinations
import pickle
import re
import pyneb as pn

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def buildANDgenerate(clobber=True):
    """
    Convenience wrapper to build and generate all the files needed for the TDOSE run

    --- Needs to be updated to be used as of 171019 ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    uves.buildANDgenerate()

    """
    LAEinfofile = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits'
    uves.build_LAEfitstable(fitsname=LAEinfofile,clobber=clobber)

    sourcecatdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_sourcecats/'
    uves.gen_LAEsourceCats(sourcecatdir,LAEinfofile,modelcoord=True)

    SETUPinfofile = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_setupfiles/MUSEWide_infofile_arche_PSFupdate_LAEs.txt'
    uves.gen_TDOSEsetupfiles(SETUPinfofile,clobber=clobber)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def run_TDOSEextraction():
    """
    Command (to copy-paste into arche) to run TDOSE on setup files generated with uves.gen_TDOSEsetupfiles()

    --- EXAMPLE OF USE ---
    copy-past into Max terminal (for copying over files to arche) and on arche (for running TDOSE)

    """
    # ---------------------------- Copying over files from Mac ----------------------------
    # scp /Users/kschmidt/work/MUSE/uvEmissionlineSearch/ref_image_galfit_models/*.fits kasper@arche.aip.de:/store/data/musewide/TDOSE/ref_image_galfit_models/

    # scp /Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_setupfiles/*candels*.txt kasper@arche.aip.de:/store/data/musewide/TDOSE/tdose_setupfiles/

    # scp /Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_sourcecats/*.fits kasper@arche.aip.de:/store/data/musewide/TDOSE/tdose_sourcecats/

    # ------------------------ Running TDOSE on Arche - FewFileRun ------------------------
    # mkdir tdose_models tdose_cutouts tdose_spectra
    # ur_setup
    # ipython
    import tdose, glob
    import numpy as np
    Nsessions = 1

    setupfiles = [glob.glob('/store/data/musewide/TDOSE/tdose_setupfiles/MUSEWide_tdose_setup_LAEs_candels-*.txt')[0]] # COSMOS 06
    setupfiles = [glob.glob('/store/data/musewide/TDOSE/tdose_setupfiles/MUSEWide_tdose_setup_LAEs_candels-*.txt')[1]] # CDFS 01

    bundles, paralleldic = tdose.perform_extractions_in_parallel(setupfiles,Nsessions=Nsessions,clobber=True,performcutout=True,store1Dspectra=True,plot1Dspectra=True,generateFullFoVmodel=True,generateOverviewPlots=True,skipextractedobjects=False,logterminaloutput=True,verbosePE=True,verbosefull=True)

    # -------------------------- Running TDOSE on Arche - Full Run -------------------------
    # mkdir tdose_models, tdose_cutouts, tdose_spectra
    # ur_setup
    # nice ipython
    import tdose, glob
    import numpy as np
    Nsessions = 30

    setupfiles = glob.glob('/store/data/musewide/TDOSE/tdose_setupfiles/MUSEWide_tdose_setup_LAEs_candels-*.txt')

    bundles, paralleldic = tdose.perform_extractions_in_parallel(setupfiles,Nsessions=Nsessions,clobber=True,performcutout=True,store1Dspectra=True,plot1Dspectra=True,generateFullFoVmodel=True,generateOverviewPlots=True,skipextractedobjects=True,logterminaloutput=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_LAEfitstable(fitsname='./LAEinfoRENAME.fits',genDS9region=False,clobber=False,verbose=True):
    """
    Building a fits table containing information on the sources.
    Generated by combining multiple sources of information.

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    objinfofile = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats.fits'
    uves.build_LAEfitstable(fitsname=objinfofile,clobber=False)

    uves.append_JK100fieldinfocat(objinfofile=objinfofile, overwrite=False)
    #uves.append_JK100fieldinfocatBAD(objinfofile=objinfofile, matchthresh=0.20, overwrite=False)
    #uves.append_JKthesisCat2maininfofile(objinfofile=objinfofile, objrmatch=0.20, overwrite=False)

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading fits catalogs for LAEs:')
    catPSF            = '/Users/kschmidt/work/catalogs/MUSE_GTO/psf_all_Converted_cleaned.fits'
    catE24eltab       = '/Users/kschmidt/work/catalogs/MUSE_GTO/MW_1-24_emline_table_v3.2.fits'

    if verbose: print('   '+catPSF)
    datPSF      = afits.open(catPSF)[1].data
    if verbose: print('   Columns: '+str(datPSF.dtype.names)+'\n')

    if verbose: print('   '+catE24eltab)
    datE24eltab = afits.open(catE24eltab)[1].data
    if verbose: print('   Columns: '+str(datE24eltab.dtype.names)+'\n')

    catE24main        = '/Users/kschmidt/work/catalogs/MUSE_GTO/MW_1-24_main_table_v3.2.fits'
    catE36main        = '/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_e36_v1.0.fits'
    catE40main        = '/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_e40_v0.9.fits'
    catUDFSmain       = '/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_mosaic_shallow_v0.9.fits'
    catU10main        = '/Users/kschmidt/work/catalogs/MUSE_GTO/object_catalog_udf-10_v0.9.fits'
    catUDFmain        = '/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_udf-mosaic_v0.9.fits'
    # IDs in UDF mosaic catalog consist of 4 digits id followed by 4 digits running line number

    if verbose: print('   '+catE24main)
    datE24main  = afits.open(catE24main)[1].data
    if verbose: print('   Columns: '+str(datE24main.dtype.names)+'\n')

    if verbose: print('   '+catE36main)
    datE36main  = afits.open(catE36main)[1].data
    if verbose: print('   Columns: '+str(datE36main.dtype.names)+'\n')

    if verbose: print('   '+catE40main)
    datE40main  = afits.open(catE40main)[1].data
    if verbose: print('   Columns: '+str(datE40main.dtype.names)+'\n')

    if verbose: print('   '+catUDFSmain)
    datUDFSmain  = afits.open(catUDFSmain)[1].data
    if verbose: print('   Columns: '+str(datUDFSmain.dtype.names)+'\n')

    if verbose: print('   '+catUDFmain)
    datUDFmain  = afits.open(catUDFmain)[1].data
    if verbose: print('   Columns: '+str(datUDFmain.dtype.names)+'\n')

    if verbose: print('   '+catU10main)
    datU10main= afits.open(catU10main)[1].data
    if verbose: print('   Columns: '+str(datU10main.dtype.names)+'\n')


    catE24lineprops   = '/Users/kschmidt/work/catalogs/MUSE_GTO/MW_1-24_v3.1_LAEs_line_props.fits'
    catE36lineprops   = '/Users/kschmidt/work/catalogs/MUSE_GTO/e36_emline_master_v1.0_LAEs_line_props.fits'
    catE40lineprops   = 'None'
    catUDFSlineprops  = 'None'
    catUDFlineprops   = 'None'
    catU10lineprops   = 'None'

    if verbose: print('   '+catE24lineprops)
    datE24lp  = afits.open(catE24lineprops)[1].data
    if verbose: print('   Columns: '+str(datE24lp.dtype.names)+'\n')

    if verbose: print('   '+catE36lineprops)
    datE36lp    = afits.open(catE36lineprops)[1].data
    if verbose: print('   Columns: '+str(datE36lp.dtype.names)+'\n')

    catLyaEW          = '/Users/kschmidt/work/catalogs/MUSE_GTO/fluxes_EWs_line_props.fits'
    if verbose: print('   '+catLyaEW)
    datLyaEW  = afits.open(catLyaEW)[1].data
    if verbose: print('   Columns: '+str(datLyaEW.dtype.names)+'\n')

    catLyaJKthesis    = '/Users/kschmidt/work/catalogs/MUSE_GTO/kerutt_LAEparameters190926_EWs_0_clumps_ratio_line_props.fits'
    if verbose: print('   '+catLyaJKthesis)
    datLyaJKthesis  = afits.open(catLyaJKthesis)[1].data
    if verbose: print('   Columns: '+str(datLyaJKthesis.dtype.names)+'\n')

    catGuo     = '/Users/kschmidt/work/catalogs/guo/CANDELS.GOODSS.F160W.v1.fits'
    if verbose: print('   '+catGuo)
    datGuo     = afits.open(catGuo)[1].data

    catSkeltonGS  = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
    if verbose: print('   '+catSkeltonGS)
    datSkeltonGS  = afits.open(catSkeltonGS)[1].data

    catSkeltonCOS = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat.FITS'
    if verbose: print('   '+catSkeltonCOS)
    datSkeltonCOS = afits.open(catSkeltonCOS)[1].data

    catRafelski = '/Users/kschmidt/work/catalogs/rafelski/uvudf_rafelski_2015.fits'
    if verbose: print('   '+catRafelski)
    datRafelski = afits.open(catRafelski)[1].data

    #catLaigle = '/Users/kschmidt/work/catalogs/COSMOS2015_Laigle_v1.1.fits'
    catLaigle = '/Users/kschmidt/work/catalogs/laigle/COSMOS2015_Laigle_v1.1_candelsregion.fits' # generated with TOPCAT
    if verbose: print('   '+catLaigle )
    datLaigle  = afits.open(catLaigle)[1].data

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Counting LAEs and putting together ID list')
    e24_ids  = datE24main['UNIQUE_ID'].astype(str)
    e36_ids  = datE36main['ID'].astype(str)
    e40_ids  = datE40main['ID'].astype(str)
    udfs_ids = datUDFSmain['ID'].astype(str)
    udf_ids  = np.asarray(["6"+str("%08d" % udfid) for udfid in datUDFmain['ID']])
    u10_ids  = np.asarray(["7"+str("%08d" % u10id) for u10id in datU10main['ID']])
    objids   = []

    zcut     = 2.7 # LAEs
    zcut     = 1.5 # UVemitters; objects with potential CIII]1909 (or below) and no OII
    # - - - - - - - - - - - - - E24 - - - - - - - - - - - - - -
    for ii,id in enumerate(e24_ids):
        if datE24main['Z'][ii] > zcut:
            objids.append( id )
    # - - - - - - - - - - - - - E36 - - - - - - - - - - - - - -
    for ii,id in enumerate(e36_ids):
        if datE36main['REDSHIFT'][ii] > zcut:
            objids.append( id )
    # - - - - - - - - - - - - - E40 - - - - - - - - - - - - - -
    for ii,id in enumerate(e40_ids):
        if datE40main['REDSHIFT'][ii] > zcut:
            objids.append( id )
    # - - - - - - - - - - - - - UDF Shallow - - - - - - - - - - - - - -
    for ii,id in enumerate(udfs_ids):
        if datUDFSmain['REDSHIFT'][ii] > zcut:
            objids.append( id )
    # - - - - - - - - - - - - - UDF Mosaic - - - - - - - - - - - - - -
    for ii,id in enumerate(udf_ids):
        if datUDFmain['REDSHIFT'][ii] > zcut:
            objids.append( id )
    # - - - - - - - - - - - - - UDF-10 - - - - - - - - - - - - - -
    for ii,id in enumerate(u10_ids):
        if datU10main['REDSHIFT'][ii] > zcut:
            objids.append( id )
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    objids = np.sort(np.asarray(objids).astype(int))
    # objids[0] = 1: CDFS
    # objids[0] = 2: COSMOS
    # objids[0] = 3: Parallel hudf09-1
    # objids[0] = 4: Parallel hudf09-2
    # objids[0] = 5: UDF-mosaic-shallow
    # objids[0] = 6: UDF-mosaic
    # objids[0] = 7: UDF-10
    #objids = np.array([106003018,131016105,153024080,206004030,302038138,404010192,509084195,600100628,614564367,720060067])
    NLAEs  = len(objids)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Assembling info for the '+str(NLAEs)+' LAEs found')
    galfitmodeldir  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/imgblocks_josieGALFITmodels/'
    redshifts       = []
    ras             = []
    decs            = []
    pointing        = []
    x_image         = []
    y_image         = []
    name_model      = []
    N_model_comp    = []
    ras_model       = []
    decs_model      = []
    delta_coords    = []
    x_image_model   = []
    y_image_model   = []

    # v v v    line props    v v v
    z_vac_red           = []
    z_vac_error         = []
    z_vac_mean          = []
    num_peaks           = []
    fwhm_A              = []
    fwhm_A_std          = []
    fwhm_kms            = []
    fwhm_kms_std        = []
    peak_sep_A          = []
    peak_sep_A_std      = []
    peak_sep_kms        = []
    peak_sep_kms_std    = []
    sum_fit             = []
    sum_fit_std         = []
    sum_fit_blue        = []
    sum_fit_blue_std    = []
    sum_fit_red         = []
    sum_fit_red_std     = []
    sum_lsdcat          = []
    sum_lsdcat_std      = []

    red_peak_shift_V18_kms      = []  # Red peak shift estimate based on A. Verhamme et al. (2017) relations
    red_peak_shift_V18_kms_err  = []
    z_sys_V18                   = []  # Systemic redshift estimate based on A. Verhamme et al. (2017) relations
    z_sys_V18_err               = []

    # v v v    Lya EW props  v v v
    EW_0                 = []
    EW_0_err             = []
    beta                 = []
    beta_err             = []
    flux_acs_606w        = []
    flux_err_acs_606w    = []
    flux_acs_775w        = []
    flux_err_acs_775w    = []
    flux_acs_814w        = []
    flux_err_acs_814w    = []
    flux_wfc3_125w       = []
    flux_err_wfc3_125w   = []
    flux_wfc3_160w       = []
    flux_err_wfc3_160w   = []

    # v v v    Phot Cat Match  v v v
    idsGuo      = []
    sepsGuo     = []
    rasGuo      = []
    decsGuo     = []

    idsSkelton  = []
    sepsSkelton = []
    rasSkelton  = []
    decsSkelton = []

    idsRafelski  = []
    sepsRafelski = []
    rasRafelski  = []
    decsRafelski = []

    idsLaigle  = []
    sepsLaigle = []
    rasLaigle  = []
    decsLaigle = []

    leadline    = []
    leadlineSN  = []
    leadlineCONF= []

    # v v v  Get duplication info v v v
    ids_ignore, duplicate_dictionary = uves.get_object_duplication_list(matchtol=0.25)

    # dupspecial0 = [135010177,135010178,    601621599,601621600,    608733593,608733594]
    # dupspecial1 = [608933623,722290970,    722290971]
    # dupspecial2 = [610943897,610943898,722661024,722661025]

    duplication_main = []
    for ii,id in enumerate(objids):
        if verbose:
            infostr = '  >Getting info for '+str(id)+' ('+str("%.5d" % (ii+1))+' / '+str("%.5d" % NLAEs)+')        '
            # if verbose: print(infostr)
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()
        imgignorstr = '_wht_'

        if id in duplicate_dictionary.keys():
            duplication_main.append(duplicate_dictionary[id])
        elif id in [135010177,135010178]:
            duplication_main.append(0)
        elif id in [601621599,601621600]:
            duplication_main.append(0)
        elif id in [608733593,608733594]:
            duplication_main.append(0)
        elif id in [608933623]:
            duplication_main.append(722290970)
        elif id in [610943897]:
            duplication_main.append(722661024)
        elif id in [610943898]:
            duplication_main.append(722661025)
        else:
            duplication_main.append(0)

        # - - - - - - - - - - GET LSDCAT COORDINATES - - - - - - - - - -
        if str(id) in e24_ids:
            pointingname = mu.gen_pointingname(id)
            pointing.append(pointingname)

            objent = np.where(datE24main['UNIQUE_ID'] == str(id))[0]
            redshifts.append(datE24main['Z'][objent][0])
            ras.append(datE24main['RA'][objent][0])
            decs.append(datE24main['DEC'][objent][0])
            ximg, yimg = mu.get_pixelpos(datE24main['RA'][objent],datE24main['DEC'][objent],pointingname,pixorigin=0,
                                         imgdir='/Users/kschmidt/work/images_MAST/MUSEWidePointings/*814*',
                                         ignorestr=imgignorstr,imgext=0,verbose=False)
            x_image.append(ximg)
            y_image.append(yimg)
            leadline.append(datE24main['LEAD_LINE'][objent][0])
            leadlineSN.append(datE24main['SN'][objent][0])
            leadlineCONF.append(datE24main['CONFIDENCE'][objent][0])
        elif str(id) in e36_ids:
            pointingname = mu.gen_pointingname(id)
            pointing.append(pointingname)

            objent = np.where(datE36main['ID'] == id)[0]
            redshifts.append(datE36main['REDSHIFT'][objent][0])
            ras.append(datE36main['RA'][objent][0])
            decs.append(datE36main['DEC'][objent][0])
            ximg, yimg = mu.get_pixelpos(datE36main['RA'][objent],datE36main['DEC'][objent],pointingname,
                                         imgdir='/Users/kschmidt/work/images_MAST/MUSEWidePointings/*814*',
                                         ignorestr=imgignorstr,imgext=0,verbose=False)
            x_image.append(ximg)
            y_image.append(yimg)
            leadline.append(datE36main['LINE_ID'][objent][0])
            leadlineSN.append(datE36main['S2N'][objent][0])
            leadlineCONF.append(datE36main['CONFIDENCE'][objent][0])
        elif str(id) in e40_ids:
            pointingname = mu.gen_pointingname(id)
            pointing.append(pointingname)

            objent = np.where(datE40main['ID'] == id)[0]
            redshifts.append(datE40main['REDSHIFT'][objent][0])
            ras.append(datE40main['RA'][objent][0])
            decs.append(datE40main['DEC'][objent][0])
            if str(id).startswith('1') or str(id).startswith('2'):
                imgdir='/Users/kschmidt/work/images_MAST/MUSEWidePointings/*814*'
            else:
                imgdir='/Users/kschmidt/work/images_MAST/MUSEWidePointings/*775*'
            ximg, yimg = mu.get_pixelpos(datE40main['RA'][objent],datE40main['DEC'][objent],pointingname,
                                         ignorestr=imgignorstr,imgdir=imgdir,imgext=0,verbose=False)
            x_image.append(ximg)
            y_image.append(yimg)
            leadline.append(datE40main['LINE_ID'][objent][0])
            leadlineSN.append(datE40main['S2N'][objent][0])
            leadlineCONF.append(datE40main['CONFIDENCE'][objent][0])
        elif str(id) in udfs_ids:
            objent = np.where(udfs_ids == str(id))[0]

            pointingname = mu.gen_pointingname(id)
            pointing.append(pointingname)

            redshifts.append(datUDFSmain['REDSHIFT'][objent][0])
            ras.append(datUDFSmain['RA'][objent][0])
            decs.append(datUDFSmain['DEC'][objent][0])
            ximg, yimg = mu.get_pixelpos(datUDFSmain['RA'][objent],datUDFSmain['DEC'][objent],pointingname+'*v1.0*',
                                         imgdir='/Users/kschmidt/work/images_MAST/MUSEWidePointings/*775*',
                                         ignorestr=imgignorstr,imgext=0,verbose=False)
            x_image.append(ximg)
            y_image.append(yimg)
            leadline.append(datUDFSmain['LINE_ID'][objent][0])
            leadlineSN.append(datUDFSmain['S2N'][objent][0])
            leadlineCONF.append(datUDFSmain['CONFIDENCE'][objent][0])
        elif str(id) in udf_ids:
            objent = np.where(udf_ids == str(id))[0]

            pointingname = 'udf-mosaic'
            pointing.append(pointingname)

            redshifts.append(datUDFmain['REDSHIFT'][objent][0])
            ras.append(datUDFmain['RA'][objent][0])
            decs.append(datUDFmain['DEC'][objent][0])
            ximg, yimg = mu.get_pixelpos(datUDFmain['RA'][objent],datUDFmain['DEC'][objent],pointingname,
                                         imgdir='/Users/kschmidt/work/images_MAST/MUSEWidePointings/*775*',
                                         ignorestr=imgignorstr,imgext=0,verbose=False)
            x_image.append(ximg)
            y_image.append(yimg)
            leadline.append(datUDFmain['LINE_ID'][objent][0])
            leadlineSN.append(datUDFmain['S2N'][objent][0])
            leadlineCONF.append(datUDFmain['CONFIDENCE'][objent][0])
        elif str(id) in u10_ids:
            objent = np.where(u10_ids == str(id))[0]

            pointingname = 'udf-10'
            pointing.append(pointingname)

            redshifts.append(datU10main['REDSHIFT'][objent][0])
            ras.append(datU10main['RA'][objent][0])
            decs.append(datU10main['DEC'][objent][0])
            ximg, yimg = mu.get_pixelpos(datU10main['RA'][objent],datU10main['DEC'][objent],'udf-mosaic_v1.5',
                                         imgdir='/Users/kschmidt/work/images_MAST/MUSEWidePointings/*775*',
                                         ignorestr=imgignorstr,imgext=0,verbose=False)
            x_image.append(ximg)
            y_image.append(yimg)
            leadline.append(datU10main['LINE_ID'][objent][0])
            leadlineSN.append(datU10main['S2N'][objent][0])
            leadlineCONF.append(datU10main['CONFIDENCE'][objent][0])
        else:
            print('\nWeird... ID= '+str(id)+' not found in E24, E36, E40, UDF10, UDF-shallow or UDF-mosaic id-list... #1\n')
            pdb.set_trace()
        # - - - - - - - - - - GET MODEL COORDINATES - - - - - - - - - -
        modelfile = glob.glob(galfitmodeldir+'imgblock_'+str("%.9d" % id)+'.fits')

        if len(modelfile) == 0:
            if verbose: print('   No model found; ')
            name_model.append("NoModelFoundIn_"+galfitmodeldir)
            N_model_comp.append(0)
            ras_model.append(0)
            decs_model.append(0)
            delta_coords.append(0)
            x_image_model.append(0)
            y_image_model.append(0)
        elif len(modelfile) > 1:
            sys.exit('Found more than one model file for '+str("%.9d" % id)+'; Found the models '+modelfile)
        else:
            refimg_hdr  = afits.open(modelfile[0])[1].header
            model_hdr   = afits.open(modelfile[0])[2].header
            comps       = []
            for hdrkey in model_hdr.keys():
                if ('COMP_' in hdrkey) & (model_hdr[hdrkey] != 'sky'):
                    comps.append(hdrkey)

            imgwcs      = wcs.WCS(tu.strip_header(refimg_hdr.copy()))

            pix_based_on_model = False
            if pix_based_on_model:
                xstr        = model_hdr['1_XC'].split(' ')
                ystr        = model_hdr['1_YC'].split(' ')

                if len(xstr) > 1:
                    xpix    = int(float(xstr[0]))
                else:
                    if verbose: print('   Model xpix has no err; ')
                    xpix    = int(float(xstr[0][1:-1]))

                if len(ystr) > 1:
                    ypix    = int(float(ystr[0]))
                else:
                    if verbose: print('   Model ypix has no err; ')
                    ypix    = int(float(ystr[0][1:-1]))
            else:
                fit_region     = model_hdr['FITSECT']
                cutrange_low_x = int(float(fit_region.split(':')[0].split('[')[-1]))
                cutrange_low_y = int(float(fit_region.split(',')[-1].split(':')[0]))
                xsize          = model_hdr['NAXIS1']
                ysize          = model_hdr['NAXIS2']

                xpix           = cutrange_low_x + int(xsize/2.)
                ypix           = cutrange_low_y + int(ysize/2.)

            if 'cdfs' in pointingname:
                skycoord    = wcs.utils.pixel_to_skycoord(xpix,ypix,imgwcs, origin=1)
            elif 'cosmos' in pointingname:
                skycoord    = wcs.utils.pixel_to_skycoord(xpix,ypix,imgwcs, origin=0)

            ra_model    = skycoord.ra.value
            dec_model   = skycoord.dec.value

            delta_coord = np.sqrt( (np.cos(np.deg2rad(dec_model))*(ras[ii]-ra_model))**2.0 + (decs[ii]-dec_model)**2.0 )

            name_model.append(modelfile[0])
            N_model_comp.append(len(comps))
            ras_model.append(ra_model)
            decs_model.append(dec_model)
            delta_coords.append(delta_coord*3600.)
            x_image_model.append(xpix)
            y_image_model.append(ypix)

        # - - - - - - - - - - ADD INFO FROM LINE PROPS TABLES - - - - - - - - - -
        if len(modelfile) == 0:
            z_vac_red.append(0.0)
            z_vac_error.append(0.0)
            z_vac_mean.append(0.0)
            num_peaks.append(0.0)
            fwhm_A.append(0.0)
            fwhm_A_std.append(0.0)
            fwhm_kms.append(0.0)
            fwhm_kms_std.append(0.0)
            peak_sep_A.append(0.0)
            peak_sep_A_std.append(0.0)
            peak_sep_kms.append(0.0)
            peak_sep_kms_std.append(0.0)
            sum_fit.append(0.0)
            sum_fit_std.append(0.0)
            sum_fit_blue.append(0.0)
            sum_fit_blue_std.append(0.0)
            sum_fit_red.append(0.0)
            sum_fit_red_std.append(0.0)
            sum_lsdcat.append(0.0)
            sum_lsdcat_std.append(0.0)
            red_peak_shift_V18_kms.append(0.0)
            red_peak_shift_V18_kms_err.append(0.0)
            z_sys_V18.append(0.0)
            z_sys_V18_err.append(0.0)
        else:
            if str(id) in e24_ids:
                objent = np.where(datE24lp['UNIQUE_ID'] == str(id))[0]
                z_vac_red.append(datE24lp['z_vac_red'][objent][0])
                z_vac_error.append(datE24lp['z_vac_error'][objent][0])
                z_vac_mean.append(datE24lp['z_vac_mean'][objent][0])
                num_peaks.append(datE24lp['num_peaks'][objent][0])
                fwhm_A.append(datE24lp['fwhm_A'][objent][0])
                fwhm_A_std.append(datE24lp['fwhm_A_std'][objent][0])
                fwhm_kms.append(datE24lp['fwhm_kms'][objent][0])
                fwhm_kms_std.append(datE24lp['fwhm_kms_std'][objent][0])
                peak_sep_A.append(datE24lp['peak_sep_A'][objent][0])
                peak_sep_A_std.append(datE24lp['peak_sep_A_std'][objent][0])
                peak_sep_kms.append(datE24lp['peak_sep_kms'][objent][0])
                peak_sep_kms_std.append(datE24lp['peak_sep_kms_std'][objent][0])
                sum_fit.append(datE24lp['sum_fit'][objent][0])
                sum_fit_std.append(datE24lp['sum_fit_std'][objent][0])
                sum_fit_blue.append(datE24lp['sum_fit_blue'][objent][0])
                sum_fit_blue_std.append(datE24lp['sum_fit_blue_std'][objent][0])
                sum_fit_red.append(datE24lp['sum_fit_red'][objent][0])
                sum_fit_red_std.append(datE24lp['sum_fit_red_std'][objent][0])
                sum_lsdcat.append(datE24lp['sum_lsdcat'][objent][0])
                sum_lsdcat_std.append(datE24lp['sum_lsdcat_std'][objent][0])
            elif str(id) in e36_ids:
                objent = np.where(datE36lp['UNIQUE_ID'] == str(id))[0]
                z_vac_red.append(datE36lp['z_vac_red'][objent][0])
                z_vac_error.append(datE36lp['z_vac_error'][objent][0])
                z_vac_mean.append(datE36lp['z_vac_mean'][objent][0])
                num_peaks.append(datE36lp['num_peaks'][objent][0])
                fwhm_A.append(datE36lp['fwhm_A'][objent][0])
                fwhm_A_std.append(datE36lp['fwhm_A_std'][objent][0])
                fwhm_kms.append(datE36lp['fwhm_kms'][objent][0])
                fwhm_kms_std.append(datE36lp['fwhm_kms_std'][objent][0])
                peak_sep_A.append(datE36lp['peak_sep_A'][objent][0])
                peak_sep_A_std.append(datE36lp['peak_sep_A_std'][objent][0])
                peak_sep_kms.append(datE36lp['peak_sep_kms'][objent][0])
                peak_sep_kms_std.append(datE36lp['peak_sep_kms_std'][objent][0])
                sum_fit.append(datE36lp['sum_fit'][objent][0])
                sum_fit_std.append(datE36lp['sum_fit_std'][objent][0])
                sum_fit_blue.append(datE36lp['sum_fit_blue'][objent][0])
                sum_fit_blue_std.append(datE36lp['sum_fit_blue_std'][objent][0])
                sum_fit_red.append(datE36lp['sum_fit_red'][objent][0])
                sum_fit_red_std.append(datE36lp['sum_fit_red_std'][objent][0])
                sum_lsdcat.append(datE36lp['sum_lsdcat'][objent][0])
                sum_lsdcat_std.append(datE36lp['sum_lsdcat_std'][objent][0])
            elif (str(id) in e40_ids) or (str(id) in udfs_ids) or (str(id) in udf_ids) or (str(id) in u10_ids):
                z_vac_red.append(-99)
                z_vac_error.append(-99)
                z_vac_mean.append(-99)
                num_peaks.append(-99)
                fwhm_A.append(-99)
                fwhm_A_std.append(-99)
                fwhm_kms.append(-99)
                fwhm_kms_std.append(-99)
                peak_sep_A.append(-99)
                peak_sep_A_std.append(-99)
                peak_sep_kms.append(-99)
                peak_sep_kms_std.append(-99)
                sum_fit.append(-99)
                sum_fit_std.append(-99)
                sum_fit_blue.append(-99)
                sum_fit_blue_std.append(-99)
                sum_fit_red.append(-99)
                sum_fit_red_std.append(-99)
                sum_lsdcat.append(-99)
                sum_lsdcat_std.append(-99)
            else:
                print('Weird... ID not found in E24, E36, E40, UDF-10, UDF-shallow or UDF-mosaic id-list... #2')
                pdb.set_trace()

            if peak_sep_kms[ii] != 0.0:
                rp_shift_V18_kms     = 1.00 * peak_sep_kms[ii]/2.
                rp_shift_V18_kms_err = np.abs(rp_shift_V18_kms) * \
                                        np.sqrt( (peak_sep_kms_std[ii]/peak_sep_kms[ii])**2 + (0.04/1.00)**2)
            elif peak_sep_kms[ii] == -99:
                rp_shift_V18_kms     = -99
                rp_shift_V18_kms_err = -99
            else:
                rp_shift_V18_kms     = 0.86 * fwhm_kms[ii]
                rp_shift_V18_kms_err = np.abs(rp_shift_V18_kms) * \
                                        np.sqrt( (fwhm_kms_std[ii]/fwhm_kms[ii])**2 + (0.04/0.86)**2)

            # Estimate systemic redshift using Lya offest from Verhamme+17 and Eq. (5) Erb+14 relating this to z_sys
            if rp_shift_V18_kms == -99:
                z_sys           = -99
                z_sys_err       = -99
            else:
                c_val           = astropy.constants.c.value/1000.
                numerator       = ( z_vac_red[ii] - rp_shift_V18_kms/c_val)
                numerator_err   = np.sqrt( z_vac_error[ii]**2.0 + (rp_shift_V18_kms_err/rp_shift_V18_kms)**2.0 )
                denominator     = (rp_shift_V18_kms/c_val + 1.0)
                denominator_err = rp_shift_V18_kms_err/np.abs(rp_shift_V18_kms)
                z_sys           = numerator / denominator
                z_sys_err       = np.abs(z_sys) * \
                                  np.sqrt( (numerator_err/numerator)**2.0 + (denominator_err/denominator)**2.0 )

            red_peak_shift_V18_kms.append(rp_shift_V18_kms)
            red_peak_shift_V18_kms_err.append(rp_shift_V18_kms_err)
            z_sys_V18.append(z_sys)
            z_sys_V18_err.append(z_sys_err)

        # - - - - - - - - - - ADD INFO FROM EW LINE PROPS TABLE - - - - - - - - - -
        if len(modelfile) == 0:
            EW_0.append(0.0)
            EW_0_err.append(0.0)
            beta.append(0.0)
            beta_err.append(0.0)
            flux_acs_606w.append(0.0)
            flux_err_acs_606w.append(0.0)
            flux_acs_775w.append(0.0)
            flux_err_acs_775w.append(0.0)
            flux_acs_814w.append(0.0)
            flux_err_acs_814w.append(0.0)
            flux_wfc3_125w.append(0.0)
            flux_err_wfc3_125w.append(0.0)
            flux_wfc3_160w.append(0.0)
            flux_err_wfc3_160w.append(0.0)
        else:
            objent = np.where(datLyaEW['IDs'] == str(id))[0]
            EW_0.append(datLyaEW['EW_0'][objent][0])
            EW_0_err.append(datLyaEW['EW_0_err'][objent][0])
            beta.append(datLyaEW['beta'][objent][0])
            beta_err.append(datLyaEW['beta_err'][objent][0])
            flux_acs_606w.append(datLyaEW['flux_acs_606w'][objent][0])
            flux_err_acs_606w.append(datLyaEW['flux_err_acs_606w'][objent][0])
            flux_acs_775w.append(datLyaEW['flux_acs_775w'][objent][0])
            flux_err_acs_775w.append(datLyaEW['flux_err_acs_775w'][objent][0])
            flux_acs_814w.append(datLyaEW['flux_acs_814w'][objent][0])
            flux_err_acs_814w.append(datLyaEW['flux_err_acs_814w'][objent][0])
            flux_wfc3_125w.append(datLyaEW['flux_wfc3_125w'][objent][0])
            flux_err_wfc3_125w.append(datLyaEW['flux_err_wfc3_125w'][objent][0])
            flux_wfc3_160w.append(datLyaEW['flux_wfc3_160w'][objent][0])
            flux_err_wfc3_160w.append(datLyaEW['flux_err_wfc3_160w'][objent][0])

        # - - - - - - - - - - ADD MATCHES TO PHOTOMETRIC CATALOGS - - - - - - - - - -
        coordMUSE = SkyCoord(ra=ras[ii]*u.degree, dec=decs[ii]*u.degree)

        if str(id)[0] in ['1','5','6','7']: # <------- Match to Guo catalog
            coordGuo  = SkyCoord(ra=datGuo['RA']*u.degree, dec=datGuo['DEC']*u.degree)
            entGuo, match_2DGuo, match_3DGuo = coordMUSE.match_to_catalog_sky(coordGuo)
            idsGuo.append(datGuo['ID'][entGuo])
            sepsGuo.append(match_2DGuo.value[0] * 3600.0)
            rasGuo.append(datGuo['RA'][entGuo])
            decsGuo.append(datGuo['DEC'][entGuo])
        else:
            idsGuo.append(-99)
            sepsGuo.append(-99)
            rasGuo.append(-99)
            decsGuo.append(-99)

        if str(id)[0] in ['1','3','4','5','6','7']: # <------- Match to Skelton GOODS-SOUTH catalog
            coordSkelton  = SkyCoord(ra=datSkeltonGS['RA']*u.degree, dec=datSkeltonGS['DEC']*u.degree)
            entSkelton, match_2DSkelton, match_3DSkelton = coordMUSE.match_to_catalog_sky(coordSkelton)
            idsSkelton.append(datSkeltonGS['ID'][entSkelton])
            sepsSkelton.append(match_2DSkelton.value[0] * 3600.0)
            rasSkelton.append(datSkeltonGS['RA'][entSkelton])
            decsSkelton.append(datSkeltonGS['DEC'][entSkelton])
        elif str(id)[0] in ['2']: # <------- Match to Skelton COSMOS catalog
            coordSkelton  = SkyCoord(ra=datSkeltonCOS['RA']*u.degree, dec=datSkeltonCOS['DEC']*u.degree)
            entSkelton, match_2DSkelton, match_3DSkelton = coordMUSE.match_to_catalog_sky(coordSkelton)
            idsSkelton.append(datSkeltonCOS['ID'][entSkelton])
            sepsSkelton.append(match_2DSkelton.value[0] * 3600.0)
            rasSkelton.append(datSkeltonCOS['RA'][entSkelton])
            decsSkelton.append(datSkeltonCOS['DEC'][entSkelton])
        else:
            idsSkelton.append(-99)
            sepsSkelton.append(-99)
            rasSkelton.append(-99)
            decsSkelton.append(-99)

        if str(id)[0] in ['1','5','6','7']: # <------- Match to Rafelski catalog
            coordRafelski  = SkyCoord(ra=datRafelski['RA']*u.degree, dec=datRafelski['DEC']*u.degree)
            entRafelski, match_2DRafelski, match_3DRafelski = coordMUSE.match_to_catalog_sky(coordRafelski)
            idsRafelski.append(datRafelski['ID'][entRafelski])
            sepsRafelski.append(match_2DRafelski.value[0] * 3600.0)
            rasRafelski.append(datRafelski['RA'][entRafelski])
            decsRafelski.append(datRafelski['DEC'][entRafelski])
        else:
            idsRafelski.append(-99)
            sepsRafelski.append(-99)
            rasRafelski.append(-99)
            decsRafelski.append(-99)

        if str(id)[0] in ['2']: # <------- Match to Laigle cosmos catalog
            coordLaigle  = SkyCoord(ra=datLaigle['ALPHA_J2000']*u.degree, dec=datLaigle['DELTA_J2000']*u.degree)
            entLaigle, match_2DLaigle, match_3DLaigle = coordMUSE.match_to_catalog_sky(coordLaigle)
            idsLaigle.append(datLaigle['NUMBER'][entLaigle])
            sepsLaigle.append(match_2DLaigle.value[0] * 3600.0)
            rasLaigle.append(datLaigle['ALPHA_J2000'][entLaigle])
            decsLaigle.append(datLaigle['DELTA_J2000'][entLaigle])
        else:
            idsLaigle.append(-99)
            sepsLaigle.append(-99)
            rasLaigle.append(-99)
            decsLaigle.append(-99)

    if verbose: print('\n   done...')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Defining fits table and filling it with data')
    c1  = pyfitsOLD.Column(name='id', format='J', unit='', array=objids)
    c2  = pyfitsOLD.Column(name='pointing', format='A30', unit='', array=pointing)
    c3  = pyfitsOLD.Column(name='ra', format='D', unit='DEF', array=ras)
    c4  = pyfitsOLD.Column(name='dec', format='D', unit='DEG', array=decs)
    c5  = pyfitsOLD.Column(name='redshift', format='D', unit='', array=redshifts)
    c6  = pyfitsOLD.Column(name='x_image', format='D', unit='PIXEL', array=x_image)
    c7  = pyfitsOLD.Column(name='y_image', format='D', unit='PIXEL', array=y_image)

    c8  = pyfitsOLD.Column(name='modelname', format='A110', unit='', array=name_model)
    c9  = pyfitsOLD.Column(name='Nmodelcomponents', format='D', unit='DEF', array=N_model_comp)
    c10 = pyfitsOLD.Column(name='ra_model', format='D', unit='DEF', array=ras_model)
    c11 = pyfitsOLD.Column(name='dec_model', format='D', unit='DEG', array=decs_model)
    c12 = pyfitsOLD.Column(name='deltacoord', format='D', unit='DEG', array=delta_coords)
    c13 = pyfitsOLD.Column(name='x_image_model', format='D', unit='PIXEL', array=x_image_model)
    c14 = pyfitsOLD.Column(name='y_image_model', format='D', unit='PIXEL', array=y_image_model)

    c15 = pyfitsOLD.Column(name='z_vac_red', format='D', unit='', array=z_vac_red)
    c16 = pyfitsOLD.Column(name='z_vac_error', format='D', unit='', array=z_vac_error)
    c17 = pyfitsOLD.Column(name='z_vac_mean', format='D', unit='', array=z_vac_mean)
    c18 = pyfitsOLD.Column(name='num_peaks', format='D', unit='', array=num_peaks)
    c19 = pyfitsOLD.Column(name='fwhm_A', format='D', unit='A', array=fwhm_A)
    c20 = pyfitsOLD.Column(name='fwhm_A_std', format='D', unit='A', array=fwhm_A_std)
    c21 = pyfitsOLD.Column(name='fwhm_kms', format='D', unit='KM/S', array=fwhm_kms)
    c22 = pyfitsOLD.Column(name='fwhm_kms_std', format='D', unit='KM/S', array=fwhm_kms_std)
    c23 = pyfitsOLD.Column(name='peak_sep_A', format='D', unit='A', array=peak_sep_A)
    c24 = pyfitsOLD.Column(name='peak_sep_A_std', format='D', unit='A', array=peak_sep_A_std)
    c25 = pyfitsOLD.Column(name='peak_sep_kms', format='D', unit='KM/S', array=peak_sep_kms)
    c26 = pyfitsOLD.Column(name='peak_sep_kms_std', format='D', unit='KM/S', array=peak_sep_kms_std)
    c27 = pyfitsOLD.Column(name='sum_fit', format='D', unit='1e-20*ERG/S/CM**2', array=sum_fit)
    c28 = pyfitsOLD.Column(name='sum_fit_std', format='D', unit='1e-20*ERG/S/CM**2', array=sum_fit_std)
    c29 = pyfitsOLD.Column(name='sum_fit_blue', format='D', unit='1e-20*ERG/S/CM**2', array=sum_fit_blue)
    c30 = pyfitsOLD.Column(name='sum_fit_blue_std', format='D', unit='1e-20*ERG/S/CM**2', array=sum_fit_blue_std)
    c31 = pyfitsOLD.Column(name='sum_fit_red', format='D', unit='1e-20*ERG/S/CM**2', array=sum_fit_red)
    c32 = pyfitsOLD.Column(name='sum_fit_red_std', format='D', unit='1e-20*ERG/S/CM**2', array=sum_fit_red_std)
    c33 = pyfitsOLD.Column(name='sum_lsdcat', format='D', unit='1e-20*ERG/S/CM**2', array=sum_lsdcat)
    c34 = pyfitsOLD.Column(name='sum_lsdcat_std', format='D', unit='1e-20*ERG/S/CM**2', array=sum_lsdcat_std)

    c35 = pyfitsOLD.Column(name='red_peak_shift_V18_kms', format='D', unit='KM/S', array=red_peak_shift_V18_kms)
    c36 = pyfitsOLD.Column(name='red_peak_shift_V18_kms_err', format='D', unit='KM/S', array=red_peak_shift_V18_kms_err)
    c37 = pyfitsOLD.Column(name='z_sys_V18', format='D', unit='', array=z_sys_V18)
    c38 = pyfitsOLD.Column(name='z_sys_V18_err', format='D', unit='', array=z_sys_V18_err)

    c39 = pyfitsOLD.Column(name='EW_0', format='D', unit='A', array=EW_0)
    c40 = pyfitsOLD.Column(name='EW_0_err', format='D', unit='A', array=EW_0_err)
    c41 = pyfitsOLD.Column(name='beta', format='D', unit='', array=beta)
    c42 = pyfitsOLD.Column(name='beta_err', format='D', unit='', array=beta_err)
    c43 = pyfitsOLD.Column(name='flux_acs_606w', format='D', unit='ERG/S/CM**2', array=flux_acs_606w)
    c44 = pyfitsOLD.Column(name='flux_err_acs_606w', format='D', unit='ERG/S/CM**2', array=flux_err_acs_606w)
    c45 = pyfitsOLD.Column(name='flux_acs_775w', format='D', unit='ERG/S/CM**2', array=flux_acs_775w)
    c46 = pyfitsOLD.Column(name='flux_err_acs_775w', format='D', unit='ERG/S/CM**2', array=flux_err_acs_775w)
    c47 = pyfitsOLD.Column(name='flux_acs_814w', format='D', unit='ERG/S/CM**2', array=flux_acs_814w)
    c48 = pyfitsOLD.Column(name='flux_err_acs_814w', format='D', unit='ERG/S/CM**2', array=flux_err_acs_814w)
    c49 = pyfitsOLD.Column(name='flux_wfc3_125w', format='D', unit='ERG/S/CM**2', array=flux_wfc3_125w)
    c50 = pyfitsOLD.Column(name='flux_err_wfc3_125w', format='D', unit='ERG/S/CM**2', array=flux_err_wfc3_125w)
    c51 = pyfitsOLD.Column(name='flux_wfc3_160w', format='D', unit='ERG/S/CM**2', array=flux_wfc3_160w)
    c52 = pyfitsOLD.Column(name='flux_err_wfc3_160w', format='D', unit='ERG/S/CM**2', array=flux_err_wfc3_160w)

    c53 = pyfitsOLD.Column(name='id_guo'     , format='D', unit='', array=idsGuo)
    c54 = pyfitsOLD.Column(name='sep_guo'    , format='D', unit='ARCSEC', array=sepsGuo)
    c55 = pyfitsOLD.Column(name='ra_guo'     , format='D', unit='DEG', array=rasGuo)
    c56 = pyfitsOLD.Column(name='dec_guo'    , format='D', unit='DEG', array=decsGuo)

    c57 = pyfitsOLD.Column(name='id_skelton' , format='D', unit='', array=idsSkelton)
    c58 = pyfitsOLD.Column(name='sep_skelton', format='D', unit='ARCSEC', array=sepsSkelton)
    c59 = pyfitsOLD.Column(name='ra_skelton' , format='D', unit='DEG', array=rasSkelton)
    c60 = pyfitsOLD.Column(name='dec_skelton', format='D', unit='DEG', array=decsSkelton)

    c61 = pyfitsOLD.Column(name='id_rafelski' , format='D', unit='', array=idsRafelski)
    c62 = pyfitsOLD.Column(name='sep_rafelski', format='D', unit='ARCSEC', array=sepsRafelski)
    c63 = pyfitsOLD.Column(name='ra_rafelski' , format='D', unit='DEG', array=rasRafelski)
    c64 = pyfitsOLD.Column(name='dec_rafelski', format='D', unit='DEG', array=decsRafelski)

    c65 = pyfitsOLD.Column(name='id_Laigle' , format='D', unit='', array=idsLaigle)
    c66 = pyfitsOLD.Column(name='sep_Laigle', format='D', unit='ARCSEC', array=sepsLaigle)
    c67 = pyfitsOLD.Column(name='ra_Laigle' , format='D', unit='DEG', array=rasLaigle)
    c68 = pyfitsOLD.Column(name='dec_Laigle', format='D', unit='DEG', array=decsLaigle)

    c69 = pyfitsOLD.Column(name='leadline', format='10A', unit='', array=leadline)
    c70 = pyfitsOLD.Column(name='leadlineS2N', format='D', unit='', array=leadlineSN)
    c71 = pyfitsOLD.Column(name='leadlineConf', format='D', unit='', array=leadlineCONF)

    c72 = pyfitsOLD.Column(name='duplicationID', format='D', unit='', array=duplication_main)


    coldefs = pyfitsOLD.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,
                                 c9,c10,c11,c12,c13,c14,
                                 c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,
                                 c35,c36,c37,c38,
                                 c39,c40,c41,c42,c43,c44,c45,c46,c47,c48,c49,c50,c51,c52,
                                 c53,c54,c55,c56,c57,c58,c59,c60,
                                 c61,c62,c63,c64,c65,c66,c67,c68,
                                 c69,c70,c71,c72])
    th      = pyfitsOLD.new_table(coldefs) # creating default header

    # writing hdrkeys:'---KEY--',                             '----------------MAX LENGTH COMMENT-------------'
    #th.header.append(('MAG     ' , spec2D[0].header['MAG']   ,'MAG_AUTO from interlaced catalog'),end=True)

    head    = th.header
    tbHDU   = pyfitsOLD.new_table(coldefs, header=head)
    tbHDU.writeto(fitsname, clobber=clobber)
    if verbose: print('   Fits table stored in \n   '+fitsname)

    if genDS9region:
        if verbose: print(' - Generating DS9 region file')
        regionname = fitsname.replace('.fits','.reg')
        kbs.create_DS9region(regionname,ras,decs,color='magenta',circlesize=0.5,textlist=objids.astype(str),clobber=clobber)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def append_JK100fieldinfocat(objinfofile, overwrite=False, verbose=True):
    """
    Script appending the LAE parameter information from Joesie Kerutt's thesis to the objecet infofile generated with
    uves.build_LAEfitstable().

    Doing this through a direct ID match

    --- INPUT ---

    --- EXAMPLE OF USE ---
    See header of uves.build_LAEfitstable()

    """
    out_file = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    if os.path.isfile(out_file) and (overwrite == False):
        sys.exit('The output file '+out_file+' already exists and overwrite=False ')

    if verbose: print(' - Loading data ')
    # cat_jk     = '/Users/kschmidt/work/catalogs/MUSE_GTO/kerutt_LAEparameters200709_EWs_all_fields_v0p9.fits'

    cat_jk     = '/Users/kschmidt/work/catalogs/MUSE_GTO/kerutt_LAEparameters210125_EWs_all_fields_v1.0.fits'


    hdu_jk     = afits.open(cat_jk)
    dat_jk     = hdu_jk[1].data
    col_jk     = hdu_jk[1].columns

    hdu_info   = afits.open(objinfofile)
    dat_info   = hdu_info[1].data
    col_info   = hdu_info[1].columns

    for cc, colname in enumerate(hdu_jk[1].columns.names):
        hdu_jk[1].columns[cc].name = colname+'_jk100'

    col_out     = col_info + col_jk
    hdu_out     = afits.BinTableHDU.from_columns(col_out,fill=True) # filling with zeros or blanks to start from empty table

    if verbose: print(' - Fill the columns in the new table looping through info file objects')
    JKids = hdu_jk[1].data['ID_jk100']
    JKids_uvesformat = []
    for jkid in JKids:
        if (jkid < 2e7): # MOSAIC ID
            JKids_uvesformat.append(int(jkid + 6e8))
        elif (jkid < 1e8) * (jkid > 2e7): # UDF10 ID
            JKids_uvesformat.append(int(jkid + 7e8))
        else:
            JKids_uvesformat.append(jkid)

    for ii, id in enumerate(dat_info['id']):
        if verbose:
            infostr = '   Matching data for id='+str(id)+' ('+str("%.5d" % (ii+1))+' / '+str("%.5d" % len(dat_info['id']))+')     '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        for colname in col_info.names:
            hdu_out.data[colname][ii] = hdu_info[1].data[colname][ii]

        if not str(id).startswith('5'): # ignore UDF mosaic shallow
            jkent   = np.where(JKids_uvesformat == id)[0]
            objz    = hdu_info[1].data['redshift'][ii]

            if len(jkent) == 0:
                if objz > 2.9:
                    print('\n WARNING: No match to infofile LAE (z>2.9) id for '+str(id))
            elif len(jkent) == 1:
                for colname in col_jk.names:
                    hdu_out.data[colname][ii] = hdu_jk[1].data[colname][jkent][0]
            elif len(jkent) > 1:
                print('\n WARNING: More than 1 match to '+str(id)+', namely '+str(JKids_uvesformat[jkent])+'arcsec')

    hdu_out.writeto(out_file, overwrite=overwrite)
    if verbose: print(' - Output written to:\n   '+out_file)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def append_JK100fieldinfocatBAD(objinfofile, matchthresh=0.2, overwrite=False, verbose=True):
    """
    Script appending the LAE parameter information from Joesie Kerutt's thesis to the objecet infofile generated with
    uves.build_LAEfitstable()

    --- INPUT ---

    --- EXAMPLE OF USE ---
    See header of uves.build_LAEfitstable()

    """

    out_file = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    if os.path.isfile(out_file) and (overwrite == False):
        sys.exit('The output file '+out_file+' already exists and overwrite=False ')

    if verbose: print(' - Loading data ')
    cat_jk     = '/Users/kschmidt/work/catalogs/MUSE_GTO/kerutt_LAEparameters200709_EWs_all_fields_v0p9.fits'
    hdu_jk     = afits.open(cat_jk)
    dat_jk     = hdu_jk[1].data
    col_jk     = hdu_jk[1].columns
    Ncol_jk    = len(dat_jk.dtype.names)

    hdu_info   = afits.open(objinfofile)
    dat_info   = hdu_info[1].data
    col_info   = hdu_info[1].columns
    Ncol_info  = len(dat_info.dtype.names)

    for cc, colname in enumerate(hdu_jk[1].columns.names):
        if colname.lower() in [cn.lower() for cn in hdu_info[1].columns.names]:
            hdu_jk[1].columns[cc].name = colname.lower()+'_jk'

    xtracol     = afits.ColDefs([afits.Column(name='sep_infoVSjk', format='D',array=np.zeros(len(dat_info))*np.nan)])

    col_out     = col_info + col_jk + xtracol
    hdu_out     = afits.BinTableHDU.from_columns(col_out)

    if verbose: print(' - Fill the columns in the new table looping through info file objects')
    JKids_all = hdu_jk[1].data['id_jk']
    for ii, id in enumerate(dat_info['id']):
        objz                   = hdu_info[1].data['redshift'][ii]
        objra                  = hdu_info[1].data['ra'][ii]
        objdec                 = hdu_info[1].data['dec'][ii]
        matchcoord             = SkyCoord(ra=objra, dec=objdec, unit='deg')
        catalogcoord           = SkyCoord(ra=hdu_jk[1].data['RA_Lya'], dec=hdu_jk[1].data['DEC_Lya'], unit='deg')
        threshcheck            = catalogcoord.separation(matchcoord) < matchthresh*u.arcsec
        goodent                = np.where(threshcheck == True)[0]

        if verbose:
            infostr = '   Matching data for id='+str(id)+' ('+str("%.5d" % (ii+1))+' / '+str("%.5d" % len(dat_info['id']))+')     '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        #--------------------------------------------------------------------------
        if len(goodent) == 0:
            if objz > 2.9:
                print('\n WARNING: No match to infofile LAE (z>2.9) coordinates for '+str(id)+'=('+str(objra)+','+str(objdec)+') below '+
                      str(matchthresh)+'arcsec')
                bestent, bestsep, d3dq = astropy.coordinates.match_coordinates_sky(matchcoord, catalogcoord, nthneighbor=1)
                print('          Best match is ID_JK='+str(JKids_all[bestent])+' at seperation '+str(bestsep.arcsec[0])+' arcsec')
                if JKids_all[bestent] == id:
                    print('          But IDs match (ID_JK='+str(JKids_all[bestent])+' = ID_info='+str(id)+') so keeping best match')
                    goodent = np.atleast_1d(bestent)
                    rsep    = bestsep.arcsec[0]
                else:
                    continue
        else:
            rsep = catalogcoord.separation(matchcoord)[goodent].arcsec
        #--------------------------------------------------------------------------
        if (len(goodent) == 1):
            if objz > 2.9:
                # print('           Redshift = '+str(objz)+', i.e. beyond 2.9 so keeping match in catalog ')
                for colname in col_info.names:
                    hdu_out.data[colname][ii] = hdu_info[1].data[colname][ii]
                for colname in col_jk.names:
                    hdu_out.data[colname][ii] = hdu_jk[1].data[colname][goodent][0]
                hdu_out.data['sep_infoVSjk'][ii] = rsep
            else:
                print('\n WARNING Redshift = '+str(objz)+', i.e. below 2.9 so non-LAE and moving on')
        #--------------------------------------------------------------------------
        elif len(goodent) > 1:
            print('\n WARNING: '+str(len(goodent))+' matches to coordinates for '+str(id)+'=('+str(objra)+','+str(objdec)+'). \n'
                  '          Matches are with entries '+str(goodent)+', i.e. IDs='+str(hdu_jk[1].data['id_jk'][goodent])+' of the JK catalog.')

            if objz > 2.9:
                print('           Redshift = '+str(objz)+', i.e. beyond 2.9 so keep match in catalog... ')
                JKids = JKids_all[goodent]
                goodent_use = None
                for jj, jkid in enumerate(JKids):
                    if (jkid < 2e7): # MOSAIC ID
                        if goodent_use is None:
                            goodent_use = goodent[jj]
                        elif (hdu_jk[1].data['id_jk'][goodent_use] < 1e8): # UDF10 ID
                            pass
                        elif (hdu_jk[1].data['id_jk'][goodent_use] > 1e8): # MUSE-Wide
                            goodent_use = goodent[jj]
                        else:
                            print('\n WARNING: There was already a MOSAIC ID kept as goodent_use...')

                    if (jkid < 1e8) * (jkid > 2e7): # UDF10 ID
                        if goodent_use is None:
                            goodent_use = goodent[jj]
                        elif (hdu_jk[1].data['id_jk'][goodent_use] < 2e7): # mosaic
                            goodent_use = goodent[jj]
                        elif (hdu_jk[1].data['id_jk'][goodent_use] > 1e8): # MUSE-Wide
                            goodent_use = goodent[jj]
                        else:
                            print('\n WARNING: There was already a UDF10 ID kept as goodent_use...')

                if goodent_use is None: # No MOSAIC or UDF10 ids in satisfying object; simply using closest match
                    print('           WARNING2: No MOSAIC or UDF10 id in list of objects - keeping best MUSE-Wide match')
                    goodent_use = goodent[np.where(rsep == np.min(rsep))[0]][0]

                for colname in col_info.names:
                    hdu_out.data[colname][ii] = hdu_info[1].data[colname][ii]
                for colname in col_jk.names:
                    hdu_out.data[colname][ii] = hdu_jk[1].data[colname][goodent_use]
                hdu_out.data['sep_infoVSjk'][ii] = rsep[np.where(goodent == goodent_use)[0]]
                print('           Kept match for ID = '+str(id)+' and matched it with '+
                      str(hdu_jk[1].data['id_jk'][goodent_use])+' at rsep='+str(hdu_out.data['sep_infoVSjk'][ii])+'arcsec')
            else:
                print('           Redshift = '+str(objz)+', i.e. below 2.9 so non-LAE and moving on')
        #--------------------------------------------------------------------------
    hdu_out.writeto(out_file, overwrite=overwrite)
    if verbose: print(' - Output written to:\n   '+out_file)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def append_JKthesisCat2maininfofile(objinfofile='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo_UVemitters_3timesUDFcats.fits', objrmatch=0.2, overwrite=False, verbose=True):
    """
    Script appending the LAE parameter information from Joesie Kerutt's thesis to the objecet infofile generated with
    uves.build_LAEfitstable()

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    uves.append_JKthesisCat2maininfofile()

    """
    out_file = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    if os.path.isfile(out_file) and (overwrite == False):
        sys.exit('The output file '+out_file+' already exists and overwrite=False ')

    if verbose: print(' - Loading data ')
    cat_jk     = '/Users/kschmidt/work/catalogs/MUSE_GTO/kerutt_LAEparameters190926_EWs_0_clumps_ratio_line_props.fits'
    hdu_jk     = afits.open(cat_jk)
    dat_jk     = hdu_jk[1].data
    col_jk     = hdu_jk[1].columns
    Ncol_jk    = len(dat_jk.dtype.names)

    hdu_info   = afits.open(objinfofile)
    dat_info   = hdu_info[1].data
    col_info   = hdu_info[1].columns
    Ncol_info  = len(dat_info.dtype.names)

    for cc, colname in enumerate(hdu_jk[1].columns.names):
        if colname.lower() in [cn.lower() for cn in hdu_info[1].columns.names]:
            hdu_jk[1].columns[cc].name = colname.lower()+'_jk'

    xtracol     = afits.ColDefs([afits.Column(name='sep_infoVSjk', format='D',array=np.zeros(len(dat_info))*np.nan)])

    col_out     = col_info + col_jk + xtracol
    hdu_out     = afits.BinTableHDU.from_columns(col_out)

    if verbose: print(' - Fill the columns in the new table looping though info file objects')
    for ii, id in enumerate(dat_info['id']):
        objra      = hdu_info[1].data['ra'][ii]
        objdec     = hdu_info[1].data['dec'][ii]
        rsep       = np.sqrt( (np.cos(np.deg2rad(objdec))*(hdu_jk[1].data['ra_jk']-objra))**2.0 +
                              (hdu_jk[1].data['dec_jk']-objdec)**2.0 )
        # minsep_ent = np.where(rsep == np.min(rsep))[0]
        minsep_ent = np.where(rsep*3600 < objrmatch)[0]
        objselids  = hdu_jk[1].data['id_jk'][minsep_ent]

        if verbose:
            infostr = '   Matching data for id='+str(id)+' ('+str("%.5d" % (ii+1))+' / '+str("%.5d" % len(dat_info['id']))+')     '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()
        jk_ent = np.where(dat_jk['id_jk'].astype(int) == int(id))[0]
        if len(jk_ent) == 1:
            for colname in col_info.names:
                hdu_out.data[colname][ii] = hdu_info[1].data[colname][ii]
            for colname in col_jk.names:
                if colname ==  'z_vac_red': pdb.set_trace()
                hdu_out.data[colname][ii] = hdu_jk[1].data[colname][jk_ent[0]]

            if (len(minsep_ent) > 1):
                useent   = minsep_ent[np.where(objselids == hdu_jk[1].data['id_jk'][jk_ent[0]])[0]]
                hdu_out.data['sep_infoVSjk'][ii] = rsep[useent]*3600
            else:
                hdu_out.data['sep_infoVSjk'][ii] = rsep[minsep_ent]*3600
        else:
            if (len(minsep_ent) > 1):
                print('\n   > WARNING: There were multiple best matches below the threshold ('+str(objrmatch)+'arcsec) for id='+
                      str(id)+', namely id_jk='+str(hdu_jk[1].data['id_jk'][minsep_ent])+' at '+str(rsep[minsep_ent]*3600)+' arcsec')
                bestent   = minsep_ent[np.where(objselids == np.min(objselids))[0]]
                print('   > Selected to match with id_jk='+str(hdu_jk[1].data['id_jk'][bestent[0]])+' (lowest id indicating UDF)\n')
                minsep_ent = bestent

            if (len(minsep_ent) == 1):
                for colname in col_info.names:
                    hdu_out.data[colname][ii] = hdu_info[1].data[colname][ii]
                for colname in col_jk.names:
                    hdu_out.data[colname][ii] = hdu_jk[1].data[colname][minsep_ent[0]]
                hdu_out.data['sep_infoVSjk'][ii] = rsep[minsep_ent]*3600
            elif (len(minsep_ent) == 0):
                for colname in col_jk.names:
                    if hdu_out.data[colname].dtype == int:
                        hdu_out.data[colname][ii] =  - 99
                    elif 'string' in hdu_out.data[colname].dtype.name:
                        hdu_out.data[colname][ii] = 'None'
                    else:
                        hdu_out.data[colname][ii] = np.nan
            else:
                print('Something went wrong in selecting one of the IDs below the threshold ')
                pdb.set_trace()


    hdu_out.writeto(out_file, overwrite=overwrite)
    if verbose: print(' - Output written to:\n   '+out_file)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_LAEidLists(sourcecatalog,skipids=True,includecomponentinfo=True,verbose=True):
    """
    Generate TDOSE setupfiles for the LAE extractions

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    sourcecatalog = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits'
    idlists = uves.get_LAEidLists(sourcecatalog)

    """
    sourcetab = afits.open(sourcecatalog)[1].data
    pointings = np.unique(np.sort(sourcetab['pointing']))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    counterNB = 0
    counterNC = 0
    ids2skip = []
    if skipids:
        ids2skip.append(121033078)  # Object with CIV being main line (conf=1) with potential Lya; no model available
        ids2skip.append(211049280)  # Potential CIII emitter (conf=1) at 2.76, i.e., no Lya in MUSE

        if includecomponentinfo:
            # IDs with no component corresponding to LAE in model, according to:
            compinfo  = open('/Users/kschmidt/work/MUSE/uvEmissionlineSearch/171012_LAEs_component_info.txt','r')
            for ll, line in enumerate(compinfo.readlines()):
                if not line.startswith('#'):
                    cols = line.split()
                    # - - - - - - - - - - - - Check for no assigned components - - - - - - - - - - - -
                    assignedcomponent = False
                    for col in cols:
                        if (len(col) == 3) & (':' in col):
                            if col.split(':')[1] == '1':
                                assignedcomponent = True

                    if not assignedcomponent:
                        counterNC = counterNC+1
                        ids2skip.append(int(cols[1]))
                        if verbose: print('   '+str(ll)+', '+str(counterNC)+
                                          '   No assigned component for '+cols[1]+' as '+
                                          ' '.join(cols[2:8])+'[...]')
                    # - - - - - - Check for neighbors, i.e., dublicate IDs in source catalogs - - - - - -
                    if ' NB ' in line:
                        counterNB = counterNB+1
                        ids2skip.append(int(cols[1]))
                        if verbose: print('   '+str(ll)+', '+str(counterNB)+
                                          "   Close neighbor causing duplicate ID'ing for "+cols[1])
                    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ids2skip  = np.unique(np.sort(np.asarray(ids2skip)))
    Nobj_skip = len(ids2skip)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Nobj      = len(sourcetab['id'])
    if verbose: print(' - Will put id lists together for the '+str(Nobj-Nobj_skip)+
                      '; (Nobj, Nobj_NB, Nobj_Ncomp, Nobj_skip) = ('+str(Nobj)+', '+str(counterNB)+', '+str(counterNC)+', 2)')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    idlists = {}

    for pp, pointing in enumerate(pointings):
        objents     = np.where(sourcetab['pointing'] == pointing)[0]
        idlist      = []
        for laeid in sourcetab['id'][objents]:
            if not laeid in ids2skip:
                idlist.append(laeid)

        if verbose: print(pointing+'    '+str(idlist).replace(', ',','))
        idlists[pointing] = idlist

    return idlists
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_LAEsourceCats(outputdir,sourcecatalog,modelcoord=False,verbose=True):
    """
    Generating MUSE-Wide pointing source catalogs for TDOSE extraction of LAEs

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    uves.gen_LAEsourceCats('/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_sourcecats/','/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits',modelcoord=True)

    """
    sourcetab = afits.open(sourcecatalog)[1].data
    pointings = np.unique(np.sort(sourcetab['pointing']))

    for pp, pointing in enumerate(pointings):
        objents     = np.where(sourcetab['pointing'] == pointing)[0]
        pointingcat = outputdir+'tdose_sourcecat_LAEs_'+pointing+'.txt'
        fout = open(pointingcat,'w')
        fout.write('# TDOSE Source catalog generated with uvEmissionlineSearch.gen_LAEsourceCats() \n')
        fout.write('# \n')
        fout.write('# parent_id id ra dec x_image y_image fluxscale \n')

        # if '02' in pointing:        pdb.set_trace()
        ids2skip = []
        ids2skip.append(121033078)  # Object with CIV being main line (conf=1) with potential Lya; no model available
        ids2skip.append(211049280)  # Potential CIII emitter (conf=1) at 2.76, i.e., no Lya in MUSE

        for objent in objents:
            objstr = ' -99  '

            if sourcetab['ID'][objent] in ids2skip:
                continue
            else:
                objstr = objstr + str(sourcetab['ID'][objent]) + ' '
                if modelcoord == True:
                    objstr = objstr + str(sourcetab['ra_model'][objent]) + ' '
                    objstr = objstr + str(sourcetab['dec_model'][objent]) + ' '
                    objstr = objstr + str(sourcetab['x_image_model'][objent]) + ' '
                    objstr = objstr + str(sourcetab['y_image_model'][objent]) + ' '
                else:
                    objstr = objstr + str(sourcetab['RA'][objent]) + ' '
                    objstr = objstr + str(sourcetab['DEC'][objent]) + ' '
                    objstr = objstr + str(sourcetab['x_image_F814W'][objent]) + ' '
                    objstr = objstr + str(sourcetab['y_image_F814W'][objent]) + ' '
                objstr = objstr + ' 1.0000 ' + ' \n'

                fout.write(objstr)
        fout.close()

        pointingcat_fits = f2a.ascii2fits(pointingcat,asciinames=True,skip_header=2,fitsformat='D',verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_LAEsourceCats_fromGALFITmodelCubeSourceCats(outputdir,sourcecatalog,modelsourcecatdir,ignore99s=False,verbose=True):
    """
    Generating MUSE-Wide pointing source catalogs for TDOSE extraction of LAEs where all objects from
    the source catalogs generated when converting GALFIT models into cubes are combined

    --- INPUT ---
    outputdir            Output directory to contain pointing source catalogs
    sourcecatalog        Source catalog of LAEs to get pointing names from
    modelsourcecatdir    Directory containing the source catalogs generated when converting
                         LAE galfit models into cubes that TDOSE can understand for the spectral extractions.
    ignore99s            Ignore objects with (parent)IDs of -99? These are the central coordinates of the models.

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves

    outputdir         = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_sourcecats/'
    sourcecat         = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits'
    modelsourcecatdir = '/Volumes/DATABCKUP1/TDOSEextractions/MW_LAEs_JKgalfitmodels/'

    uves.gen_LAEsourceCats_fromGALFITmodelCubeSourceCats(outputdir,sourcecat,modelsourcecatdir,ignore99s=False)

    """
    sourcetab = afits.open(sourcecatalog)[1].data
    pointings = np.unique(np.sort(sourcetab['pointing']))

    for pp, pointing in enumerate(pointings):
        modelsourcecats = glob.glob(modelsourcecatdir+'/*'+pointing+'*_sourcecatalog.txt')

        pointingcat     = outputdir+'tdose_sourcecat_LAEs_'+pointing+'.txt'
        fout = open(pointingcat,'w')
        fout.write('# TDOSE Source catalog generated with uvEmissionlineSearch.gen_LAEsourceCats_fromGALFITmodelCubeSourceCats() on '+tu.get_now_string()+'  \n')
        fout.write('# \n')
        fout.write('# parent_id id ra dec x_image y_image fluxscale \n')

        for modcat in modelsourcecats:
            catinfo = open(modcat,'r')
            for line in catinfo.readlines():
                if line.startswith('#'):
                    pass
                else:
                    if ignore99s:
                        if line.split()[0] == '-99':
                            pass
                        else:
                            fout.write(line)
                    else:
                        fout.write(line)
        fout.close()

        pointingcat_fits = f2a.ascii2fits(pointingcat,asciinames=True,skip_header=2,fitsformat='D',verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_LAEsourceCats_FromPhotCat(outputdir,MUSEIDlist,LAEinfo,sourcecatradius=15.0,photcat='skelton',
                                  refimgdir='/Volumes/DATABCKUP1/MUSE-Wide/hst_cutouts/',returnSeparations=False,
                                  clobber=False,verbose=True):
    """
    Generate source catalogs based on photometric catalogs for individual MUSE Wide objects
    but adding the LSDCat coordinates as seperate soource


    --- INPUT ---
    outputdir            Output directory to contain object source catalogs
    MUSEIDlist           List of objects to generate source catalogs for
    LAEinfo              LAE info file
    sourcecatradius      Radius to search and return in the source catalogs
    photcat              The photmetric cataog to match to
    refimgdir            The reference image directory used to get WCS info
    returnSeparations    Return the separations in the fluxscale column of the source catalog

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    outputdir         = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_sourcecats_noModelComponent/'
    LAEinfo           = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits'
    MUSEIDlist        = [121033078,211015198]
    #MUSEIDlist        = [101005016,117027076,123016117,125049122,141003075,144008046,146069355,201073224,202013030,202044085,203007099,204053120,206014089,207022169,209006108,211015198,212029067,213022109,215016042]
    uves.gen_LAEsourceCats_FromPhotCat(outputdir,MUSEIDlist,LAEinfo,photcat='skelton',returnSeparations=False)

    """
    LAEinfo = afits.open(LAEinfo)[1].data

    for Mid in MUSEIDlist:
        idstr = str(Mid)
        if idstr.startswith('1') or idstr.startswith('6'):
            if photcat == 'skelton': # GOODS-S (incl. UDF)
                fitscat = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
            elif photcat == 'guo':
                fitscat = '/Users/kschmidt/work/catalogs/guo/CANDELS.GOODSS.F160W.v1.fits'
            else:
                print(' WARNING photcat = '+photcat+' has no setup for GOODS-S source; using Skelton catalog ')
                fitscat = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
        elif idstr.startswith('2'): # COSMOS
            if photcat == 'skelton':
                fitscat = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat.FITS'
            else:
                print(' WARNING photcat = '+photcat+' has no setup for COSMOS source; using Skelton catalog ')
                fitscat = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat.FITS'
        # elif idstr.startswith('6'): # UDF
        #     if photcat == 'skelton':
        #         print(' WARNING photcat = '+photcat+' will be ignored, as object in UDF where Rafelski will be used')
        #         fitscat = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat.FITS'
        #     else:
        #         print(' WARNING photcat = '+photcat+' has no setup for COSMOS source; using Skelton catalog ')
        #         fitscat = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat.FITS'
        else:
            sys.exit(' MUSE-Wide id = '+idstr+' is neither in GOODS-S (starts with "1") nor in COSMOS (starts with "2")')

        objent  = np.where(LAEinfo['id'] == Mid)[0]
        objra   = LAEinfo['ra'][objent][0]
        objdec  = LAEinfo['dec'][objent][0]
        if len(objent) == 0:
            print(' WARNING No match in LAE info file for MUSE-Wide_ID = '+idstr)

        refimg       = refimgdir+'acs_814w_'+LAEinfo['pointing'][objent][0]+'_cut_v1.0.fits'
        if idstr.startswith('6'):
            refimg   = refimgdir+'hlsp_hlf_hst_acs-30mas_goodss_f775w_udf-mosaic_v1.5_sci.fits'
        imgheader    = afits.open(refimg)[0].header
        sourcelist   = [ [Mid, Mid, objra, objdec, 1]]
        outname      = outputdir+'tdose_sourcecat_from_fitscat_id'+idstr+'.txt'

        if returnSeparations:
            fluxfactor = 'separation'
        else:
            fluxfactor = 1.0

        sourcecat    = tu.gen_sourcecat_from_FitsCat(fitscat,'id','ra','dec',[objra,objdec],sourcecatradius,imgheader,
                                                     fluxfactor=fluxfactor,outname=outname,newsources=sourcelist,
                                                     clobber=clobber,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_GALFITmodelcubes(GALFITmodels,outputdir,PSFmodels=None,PSFmodelext=2,sourcecat_compinfo=None,
                         refnamebase='model_acs_814w_PPPP_cut_v1.0_idIIII_cutout2p0x2p0arcsec.fits',
                         pointsourcefile=None,pointsourcescale=1.0,ignore_radius=0.5,clobber=False,verbose=True):
    """

    Function loading galfit models from modelinputdir (assumed to be names as imgblock_ID.fits), renaming them,
    converting them to cubes and generating the corresponding source catalogs needed by TDOSE. It also generates
    a template component info file which can be edited (after copying to a new file) and be provided back to the
    script for a second run updating the the cubes and source catalogs accordingly.

    --- INPUT ---


    --- EXMAMPLE OF USE ---
    import glob
    import uvEmissionlineSearch as uves

    GALFITmodels    = glob.glob('/Users/kschmidt/work/MUSE/uvEmissionlineSearch/imgblocks_josieGALFITmodels/imgblock_*.fits')
    outputdir       = '/Volumes/DATABCKUP2/TDOSEextractions/MW_LAEs_JKgalfitmodels/'
    PSFmodels       = ['/Users/kschmidt/work/MUSE/uvEmissionlineSearch/F814Wpsfmodel_imgblock_6475.fits']*len(GALFITmodels)
    pointsourcefile = None #'/Users/kschmidt/work/MUSE/uvEmissionlineSearch/pointsourceobjects.txt'
    uves.gen_GALFITmodelcubes(GALFITmodels,outputdir,PSFmodels=PSFmodels,sourcecat_compinfo=None,pointsourcefile=pointsourcefile)

    """
    Nmodels = len(GALFITmodels)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Renaming files of '+str(Nmodels)+' models profived to GALFITmodels keyword')
    models_renamed = []
    model_ids      = []
    for modelname in GALFITmodels:
        objid    = modelname.split('block_')[-1].split('.fit')[0]
        pointing = mu.gen_pointingname(objid)
        newname  = outputdir+'/'+refnamebase.replace('IIII',str(objid)).replace('PPPP',pointing)
        cpcmd    = ' cp '+modelname+' '+newname
        if os.path.isfile(newname) & (clobber == False):
            if verbose: print(' clobber = False and '+newname+' already exists, so moving on to next file.')
        else:
            cpout = commands.getoutput(cpcmd)
        models_renamed.append(newname)
        model_ids.append(objid)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if pointsourcefile is not None:
        if verbose: print(' - Assembling list of models for objects to use point source extractions for ')
        pointsources      = np.genfromtxt(pointsourcefile,dtype=None,comments='#')
        try:
            pointsourcescales = [pointsourcescale]*len(pointsources)
            ignore_radii      = [ignore_radius]*2 # same radius in x and y dimension
        except:
            pointsourcescales = pointsourcescale
            ignore_radii      = ignore_radius
    else:
        pointsources      = None
        pointsourcescales = 'dummy'
        ignore_radii      = 'dummy'

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if PSFmodels is None:
        PSFlist = None
    else:
        if verbose: print(' - Loading PSF models ')
        if type(PSFmodels) is list:
            PSFlist = []
            if type(PSFmodelext) is not list:
                PSFmodelext = [PSFmodelext] * len(PSFmodels)
            for mm, PSFmodel in enumerate(PSFmodels):
                PSFlist.append(afits.open(PSFmodel)[PSFmodelext[mm]].data)
        else:
            PSFlist = [afits.open(PSFmodels)[PSFmodelext].data]*len(GALFITmodels)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    gen_compinfofile = True
    if sourcecat_compinfo is None:
        compinfofile = None
    else:
        if os.path.isfile(sourcecat_compinfo):
            if verbose: print(' - Will use existing point source component file provided:\n   '+sourcecat_compinfo)
            compinfofile     = sourcecat_compinfo
            gen_compinfofile = False
            if verbose: print('   (no new file/template will be generated)')
        else:
            compinfofile = None

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Building cubes from renamed GALFIT models')
    # newlist = []
    # for mod in models_renamed:
    #     if 'id12400' in mod: newlist.append(mod)
    # models_renamed = newlist

    tu.galfit_convertmodel2cube(models_renamed,includewcs=True,savecubesumimg=True,convkernels=PSFlist,
                                sourcecat_compinfo=compinfofile,normalizecomponents=True,pointsources=pointsources,
                                ignore_radius=ignore_radii,pointsourcescales=pointsourcescales,includesky=False,
                                clobber=clobber,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if gen_compinfofile:
        if verbose: print(' - Generating component info file template for source catalog updates')
        skip = False
        if sourcecat_compinfo is None:
            compinfofile = './component_info_template_RENAME_.txt'
            if os.path.isfile(compinfofile) & (clobber == False):
                if verbose: print('   ... but '+compinfofile+' exists and clobber=False, so skipping.')
                skip = True
        else:
            if os.path.isfile(compinfofile) & (clobber == False):
                if verbose: print('   ... but '+compinfofile+' exists and clobber=False, so skipping.')
                skip = True
            else:
                compinfofile = sourcecat_compinfo

        if not skip:
            fout = open(compinfofile,'w')
            fout.write("""# TDOSE source catalog components keys for J. Kerutt's 2x2 arcsec GALFIT models of the MUSE-Wide LAEs from
# the first 60 MUSE-Wide pointings.
#
# --- TEMPLATE --- generated with uvEmissionlineSearch.gen_GALFITmodelcubes() on %s
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# modelfilename        Path and name of model file
# id                   MUSE-Wide object ID
# componentinfo        Information on the model components given as ComponentNumber:InfoKey
#                      where the info keys are:  1 = object, 2 = contaminant and 3 = sky
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The following default notes are used for commenting after ">>>Notes>>>:" (appended by notes from JK's inspection)
#
# NoteKey      ShortRef           NoteExplanation
#
# ND           Non Detection      A (visual) non detection in all the optical bands (606, 775W, 814W, 850LP)
# PS           Point Source       Add a point source at central location to represent the soource
# BD           Blue detection     Detection in filters blue-wards of the Lyalpha line (filters not including the Lyalpha wavelength)
# WD           Weak drop          There is a week drop in the filters blue wards of the Lyalpha line
# MO           Model offset       The model appears offset compared to LSDCat (Lyalpha) location
# NBid         Neighbor           A neioghboring LAE (with id "id") which can potentially cause confusion/overlap of LAE spectra exists
# DFfilter     Detection Filter   Object only detected in (a few) filters. Indicate those with multiple DFfilter comments
# OCtext       Other Comment      Anything else to comment on? Follow comment by text
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# modefilename  id  componentinfo
""" % tu.get_now_string())

            for mm, GFmodel in enumerate(GALFITmodels):
                modelheader = afits.open(models_renamed[mm])[2].header
                compstring  = ' '
                for key in modelheader.keys():
                    if 'COMP_' in key:
                        compNo = key.split('OMP_')[-1]
                        if modelheader[key] == 'sky':
                            compstring = compstring + compNo + ':3  '
                        else:
                            compstring = compstring + compNo + ':?  '

                outstring = models_renamed[mm]+'  '+model_ids[mm]+'  '+compstring.ljust(50)+\
                            '     # >>>Notes>>>:  ND  PS  BD  WD  MO  NBid  DFfilter  OCtext    >>>JK notes>>>: '
                jknotes   = open('/Users/kschmidt/work/MUSE/uvEmissionlineSearch/imgblocks_josieGALFITmodels_all_ids.txt','r')
                for line in jknotes.readlines():
                    if str(model_ids[mm]) in line:
                        outstring = outstring+'  '+line.replace('\n','').replace('	','   ')+'  '
                jknotes.close()
                fout.write(outstring+' \n')
            fout.close()
            if verbose: print(' - Wrote component info to: '+compinfofile)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def inspect_GALFITmodels(modeldir='/Volumes/DATABCKUP2/TDOSEextractions/MW_LAEs_JKgalfitmodels/',
                         LAEinfofile='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits',
                         imgdir='/Volumes/DATABCKUP2/MUSE-Wide/hst_cutouts/',modelstart=1,showPhotRegions=True,
                         objids=None,verbose=True):
    """
    Script to put open DS9 windows showing galfit models so they can be inspected

    --- INPUT ---
    modeldir     Directory containing models to display
    modelstart   Where to start the inspection in list of models/objects. Useful to skip ahead in long object lists.
                 E.g. when all models in a directory are to be inspected. First model has modelstart=1
    objids       List of objects ids to display. If None, all objects found in modeldir will be displayed

    --- EXAMPLE OF USE ---

    models = glob.glob('/Volumes/DATABCKUP2/TDOSEextractions/MW_LAEs_JKgalfitmodels/model*arcsec.fits')
    tu.galfit_model_ds9region(models,clobber=True)

    uves.inspect_GALFITmodels(modelstart=3)

    """
    import commands # only works on Python2.7; has to convert to subprocess.Popen() for Python3
    LAEinfo = afits.open(LAEinfofile)[1].data

    if objids is None:
        GALFITmodels = glob.glob(modeldir+'model*arcsec.fits')
    else:
        GALFITmodels = []
        for objid in objids:
            GALFITmodels = GALFITmodels + [mod for mod in glob.glob(modeldir+'model*'+str(objid)+'*arcsec.fits')]
    GALFITmodels      = np.asarray(GALFITmodels)
    if verbose: print(' - Found '+str(len(GALFITmodels))+' GALFIT models')

    loopmodels        = GALFITmodels[modelstart-1:]

    MWregion_cosmos   = '/Users/kschmidt/work/catalogs/MUSE_GTO/MUSE-Wide_objects_cosmos.reg' #candels_cosmos_pointings-all.reg'
    MWregion_cdfs     = '/Users/kschmidt/work/catalogs/MUSE_GTO/MUSE-Wide_objects_cdfs.reg' #candels_cdfs_pointings-all.reg'

    if verbose: print(' - Will look through '+str(len(loopmodels))+' of the models starting with model number '+str(modelstart))

    ds9cmd       = "ds9 -view layout horizontal -lock frame wcs -height 650 -width 800 -tile grid layout 4 4 "
    pds9         = subprocess.Popen(ds9cmd,shell=True,executable=os.environ["SHELL"])
    time.sleep(1.1)# sleep to make sure ds9 appear in PIDlist
    for ii in np.arange(1,13):
        out = commands.getoutput('xpaset -p ds9 frame new')
    out = commands.getoutput('xpaset -p ds9 tile yes ')

    Guo_goodss_reg     = '/Users/kschmidt/work/catalogs/guo/CANDELS.GOODSS.F160W.v1_ds9.reg'
    Skelton_goodss_reg = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat_ds9.reg'
    Skelton_cosmos_reg = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat_ds9.reg'

    for mm, GFmodel in enumerate(loopmodels):
        if verbose:
            infostr = '   Displaying files for model '+str("%.5d" % (mm+1))+' / '+str("%.5d" % len(loopmodels))+' in DS9.'
            print(infostr)

        modelid      = GFmodel.split('_id')[-1][:9]
        compregion   = GFmodel.replace('.fits','_ds9region.reg')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        out = commands.getoutput('xpaset -p ds9 frame 1 ')
        out = commands.getoutput('xpaset -p ds9 file '+GFmodel+'[1]')
        out = commands.getoutput('xpaset -p ds9 regions '+compregion)
        out = commands.getoutput('xpaset -p ds9 zoom to fit ')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        out = commands.getoutput('xpaset -p ds9 frame 2 ')
        out = commands.getoutput('xpaset -p ds9 file '+GFmodel+'[2]')
        out = commands.getoutput('xpaset -p ds9 regions '+compregion)
        out = commands.getoutput('xpaset -p ds9 zoom to fit ')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        modelcubesum = GFmodel.replace('.fits','_cubesum.fits')
        if os.path.isfile(modelcubesum):
            out = commands.getoutput('xpaset -p ds9 frame 4 ')
            out = commands.getoutput('xpaset -p ds9 file '+modelcubesum)
            out = commands.getoutput('xpaset -p ds9 regions '+compregion)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        HSTcutouts   = glob.glob(imgdir+'*'+GFmodel.split('/')[-1][15:-37]+'*fits')
        for cc, HSTcutout in enumerate(HSTcutouts):
            out = commands.getoutput('xpaset -p ds9 frame '+str(5+cc))
            out = commands.getoutput('xpaset -p ds9 file '+HSTcutout)
            if 'cosmos' in HSTcutout:
                out = commands.getoutput('xpaset -p ds9 regions '+MWregion_cosmos)
                if showPhotRegions & ('814' in HSTcutout):
                    out = commands.getoutput('xpaset -p ds9 regions '+Skelton_cosmos_reg)
            else:
                out = commands.getoutput('xpaset -p ds9 regions '+MWregion_cdfs)
                if showPhotRegions & ('814' in HSTcutout):
                    out = commands.getoutput('xpaset -p ds9 regions '+Skelton_goodss_reg)
                    out = commands.getoutput('xpaset -p ds9 regions '+Guo_goodss_reg)

            out = commands.getoutput('xpaset -p ds9 scale log 1 10')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        modelcube = GFmodel.replace('.fits','_cube.fits')
        if os.path.isfile(modelcube):
            out = commands.getoutput('xpaset -p ds9 frame 3 ')
            out = commands.getoutput('xpaset -p ds9 file '+modelcube)
            out = commands.getoutput('xpaset -p ds9 regions '+compregion)
            out = commands.getoutput('xpaset -p ds9 zoom to fit ')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        narrowbandimages = glob.glob(modeldir+'*'+modelid+'*narrowbandimage*.fits')
        if len(narrowbandimages) > 0:
            for nn, nbimg in enumerate(narrowbandimages):
                out = commands.getoutput('xpaset -p ds9 frame '+str(13+nn))
                out = commands.getoutput('xpaset -p ds9 file '+nbimg)
                if 'cosmos' in GFmodel:
                    out = commands.getoutput('xpaset -p ds9 regions '+MWregion_cosmos)
                else:
                    out = commands.getoutput('xpaset -p ds9 regions '+MWregion_cdfs)
            out = commands.getoutput('xpaset -p ds9 zoom to fit ')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Printing info of object
        if verbose:
            objent =  np.where(LAEinfo['ID'] == int(modelid))
            print('   ID        =  '+modelid)
            print('   [ra,dec]  = ['+str(LAEinfo['RA'][objent][0])+','+str(LAEinfo['DEC'][objent][0])+']')
            print('   zMUSE     = '+str(LAEinfo['redshift'][objent][0]))
            lamLya    = (LAEinfo['redshift'][objent][0]+1.0) * 1216.0
            print('   lamdaLya  = '+str("%.2f" % lamLya))
            bandsLya  = uves.wavelength_in_bands(lamLya)
            print('   bandsLya  = '+str(bandsLya))
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print('\n Move on to the next model (y/n)? ')
        input = raw_input()
        if (input.lower() == 'y') or (input.lower() == 'yes'):
            continue
        else:
            if verbose: print('\n - Okay; then shutting down ')
            return
    if verbose: print('\n - Done; no more objects in loop')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def wavelength_in_bands(wavelength):
    """
    Returning band names containing a given wavelength
    Can be used to return bands where a certain emission line is included

    Band widths are taken from
    http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=HST&gname2=ACS_WFC

    Curves are at
    http://www.stsci.edu/hst/wfc3/ins_performance/UVIS_sensitivity/UVIS_Longx.jpg
    http://www.stsci.edu/hst/wfc3/ins_performance/UVIS_sensitivity/UVIS_Wide1.jpg
    http://www.stsci.edu/hst/wfc3/ins_performance/UVIS_sensitivity/UVIS_Wide2.jpg
    http://www.stsci.edu/hst/wfc3/ins_performance/IR_sensitivity/IR4_Wide1_single.jpg
    Linked from
    http://www.stsci.edu/hst/wfc3/ins_performance/ground/components/filters

    """
    infodic           = {}
    infodic['F275W']  = [2286,3120]
    infodic['F336W']  = [3014,3707]
    infodic['F435W']  = [3599,4861]
    infodic['F606W']  = [4634,7180]
    infodic['F775W']  = [6804,8632]
    infodic['F814W']  = [6885,9648]
    infodic['F850LP'] = [8007,10865]
    infodic['F105W']  = [8947,12129]
    infodic['F125W']  = [10845,14139]
    infodic['F160W']  = [13854,16999]

    bands = []
    for key in infodic.keys():
        if (wavelength >= infodic[key][0]) & (wavelength <= infodic[key][1]):
            bands.append(key)

    return bands
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_TDOSEsetupfiles(infofile,namebase='MUSEWide_tdose_setup_LAEs',clobber=False,
                        outputdir='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_setupfiles/',verbose=True):
    """
    Generate TDOSE setupfiles for the LAE extractions

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    infofile = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_setupfiles/MUSEWide_infofile_arche_PSFupdate_LAEs.txt'
    uves.gen_TDOSEsetupfiles(infofile)

    """
    tu.duplicate_setup_template(outputdir,infofile,namebase=namebase,clobber=clobber,loopcols='all',infofmt="S250",infohdr=2)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def rename_models(outputdir,sourcecatalog,cutoutsize=[2.0,2.0],clobber=False,
                  modeldir='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/imgblocks_josieGALFITmodels/',verbose=True):
    """
    Renmae GALFIT models to comply with TDOSE naming convention (i.e. so TDOSE can find the models when
    looking for them using model_*refimage+cutoutstring*)

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    outputdir     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/ref_image_galfit_models/'
    sourcecatalog = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits'
    uves.rename_models(outputdir,sourcecatalog,cutoutsize=[2.0,2.0],clobber=False)
    """
    modelfiles = glob.glob(modeldir+'/imgblock*.fits')
    sourcetab  = afits.open(sourcecatalog)[1].data
    if verbose: print(' - Found '+str(len(modelfiles))+' in modeldir to rename ')

    for oldname in modelfiles:
        id     = oldname.split('/')[-1].split('_')[-1].split('.fit')[0]
        objent = np.where(sourcetab['id'] == int(id))[0]

        if len(objent) != 1:
            print(' - No match in sourcecatalog to object '+id)
        else:
            pointing = sourcetab['pointing'][objent][0]

            if cutoutsize is None:
                cutoutstr = ''
            else:
                cutoutstr = ('_id'+str("%.9d" % float(id))+'_cutout'+str(cutoutsize[0])+
                             'x'+str(cutoutsize[1])+'arcsec').replace('.','p')

            if 'cdfs' in pointing:
                newname = outputdir+'model_acs_814w_'+pointing+'_cut_v1.0'+cutoutstr+'.fits'
            elif 'cosmos' in pointing:
                newname = outputdir+'model_acs_814w_'+pointing+'_cut_v1.0'+cutoutstr+'.fits'

            if os.path.isfile(newname) & (clobber == False):
                print(' - Clobber = False and '+newname+' already exists so no new copy made. Moving on')
            else:
                if verbose: print(' - Copying '+oldname+' to '+newname)
                shutil.copy(oldname,newname)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_ModelReferencePixelCoordinates(modeldir,pixpos='center',printcoords=True,verbose=True):
    """
    Extract the reference coordinates of the GALFIT models from the fits headers

    PROBLEM! GALFIT apparantly doesn't propogate the coordinates of the cutouts. It quotes the reference
             pixel coordinate from the image the cutout was generated from. In the case of the MUSE cutouts
             this means that the coordiantes are the reference position for the fullf-FoV muse pointing
             cutouts and not the individual object cutouts GALFIT is modeling.

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    modeldir     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/ref_image_galfit_models/'
    coordarray   = uves.get_ModelReferencePixelCoordinates(modeldir,printcoords=True,verbose=True)

    """
    modelfiles = glob.glob(modeldir+'*.fits')
    Nfiles     = len(modelfiles)
    if verbose: print(' - Found '+str(Nfiles)+' models to extract coordinates from ')

    if verbose: print(' - Looping over models and extracting coordinates from: ')
    coordarray = np.zeros(Nfiles, dtype={'names':['modelfile','xpix','ypix','ra','dec'],
                                         'formats':['a250', 'f8', 'f8', 'f8', 'f8']})
    for mm, modelfile in enumerate(modelfiles[0:5]):
        model_refimghdr = afits.open(modelfile)[1].header
        imgwcs    = wcs.WCS(tu.strip_header(model_refimghdr.copy()))

        if pixpos == 'center':
            model_shape     = afits.open(modelfile)[1].data.shape
            xpix      = int(model_shape[1]/2.)
            ypix      = int(model_shape[0]/2.)
        else:
            xpix      = pixpos[1]
            ypix      = pixpos[0]

        print(imgwcs)
        skycoord  = wcs.utils.pixel_to_skycoord(xpix,ypix,imgwcs, origin=0)
        ra        = skycoord.ra.value
        dec       = skycoord.dec.value

        if printcoords & verbose:
            print('   '+modelfile.split('/')[-1]+':  (ra,dec) = ('+str(ra)+','+str(dec)+')')

        coordarray['modelfile'][mm] = modelfile
        coordarray['xpix'][mm]      = xpix
        coordarray['ypix'][mm]      = ypix
        coordarray['ra'][mm]        = ra
        coordarray['dec'][mm]       = dec

    return coordarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_narrowbandimages(LAEinfofile,datacubestring,outputdir,linewaves=[1216,1549,1909],fwhmkey='FWHM',
                         clobber=False,verbose=True):
    """
    Generate narrow band images around the location for a set of emission lines.

    If FWHM value is found in LAEinfo file, the relation from Verhamme+17 is used to predict systemic redshift

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    LAEinfofile    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits'
    datacubestring = '/Volumes/DATABCKUP2/MUSE-Wide/datacubes_dcbgc_effnoised/DATACUBE_PPPP_v1.0_dcbgc_effnoised.fits'
    outputdir      = '/Volumes/DATABCKUP2/TDOSEextractions/MW_LAEs_JKgalfitmodels/'

    uves.gen_narrowbandimages(LAEinfofile,datacubestring,outputdir,linewaves=[1216,1549,1909],fwhmkey='FWHM',verbose=True)

    """
    LAEinfo = afits.open(LAEinfofile)[1].data

    pointings = LAEinfo['pointing']

    for pointing in np.unique(np.sort(pointings)):
        pointing_objs = np.where(pointings == pointing)[0]

        datacube = glob.glob(datacubestring.replace('PPPP',pointing))

        if len(datacube) == 0:
            if verbose: print(' -----> WARNING No data cube found globbing for ')
            if verbose: print('        '+datacubestring.replace('PPPP',pointing))
        elif len(datacube) > 1:
            if verbose: print(' -----> WARNING More than 1 data cube found globbing for ')
            if verbose: print('        '+datacubestring.replace('PPPP',pointing))
            if verbose: print('        Using the first found in the list, i.e., ')
            datacube = datacube[0]
            if verbose: print('        Extracting from: '+datacube)
        else:
            datacube = datacube[0]
            if verbose: print('\n - Extracting from: '+datacube)

        ras       = LAEinfo['ra'][pointing_objs]
        decs      = LAEinfo['dec'][pointing_objs]
        names     = LAEinfo['id'][pointing_objs].astype(str)
        redshifts = LAEinfo['redshift'][pointing_objs]

        if fwhmkey in  LAEinfo.columns.names:
            fwhms = LAEinfo[fwhmkey][pointing_objs]
        else:
            fwhms = []

        wcenters = []
        dwaves   = []
        for redshift in redshifts:
            wcen = []
            dwav = []
            if len(fwhms) != 0:
                if verbose: print(' - Estimating systemic redshift using Verhamme+17 z_sys vs Lya_FWHM relation ')
                zsys = None
            else:
                zsys = redshift

            for lw in linewaves:
                if lw == 1216:
                    wcen.append(lw*(redshift+1.0))
                    dwav.append(10)
                else:
                    wcen.append(lw*(zsys+1.0))
                    dwav.append(10)

            wcenters.append(wcen)
            dwaves.append(dwav)

        mu.create_narrowband_subcube(datacube,ras,decs,5.0,5.0,wcenters,dwaves,outputdir,names=names,clobber=clobber)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def estimate_limits(spectra,sourcecatalog,lines=['lya','civ','ciii'],deltalam=10,plot=True,verbose=True,printresults=False):
    """
    Get limits at line locations from 1D spectra

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves, glob
    spectra       = glob.glob('/Volumes/DATABCKUP1/TDOSEextractions/tdose_spectra/tdose_spectrum_candels*.fits')
    sourcecatalog = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits'
    limit_output  = uves.estimate_limits(spectra,sourcecatalog,lines=['lya','civ','ciii'],deltalam=3)

    plotbasename = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/estimatelimits'
    uves.plot_limits(sourcecatalog,plotbasename,limit_output)

    """
    if verbose: print(' - Loading source catalog ')
    sourcecat = afits.open(sourcecatalog)[1].data

    Nspec = len(spectra)
    if verbose: print(' - Will estimate limits for lines '+str(lines)+' for '+str(Nspec)+' spectra found')
    outputdic = {}
    outputdic['deltalam'] = deltalam


    for ll, line in enumerate(lines):
        if line.lower() == 'lya':
            line_lams  = [1215.6737]
            use_sys    = False
            keys       = ['lya']
        elif line.lower() == 'civ':
            line_lams  = [1548.195,1550.770]
            use_sys    = True
            keys       = ['civ1548','civ1551']
        elif line.lower() == 'ciii':
            line_lams  = [1907.00,1909.00]
            use_sys    = True
            keys       = ['ciii1907','ciii1909']
        elif line.lower() == 'heii':
            line_lams  = [1640.420]
            use_sys    = True
            keys       = ['heii1640']
        elif line.lower() == 'nv':
            line_lams  = [1238.821,1242.804]
            use_sys    = True
            keys       = ['nv1239','nv1243']

        else:
            sys.exit(' Did not find any setups for the line designated '+line)

        for ll, line_lam in enumerate(line_lams):
            ids                 = []
            # v v v   From uves.lineinfofromspec()   v v v
            fluxval             = []
            fluxerr             = []
            SNval               = []
            fluxval_Dlam        = []
            fluxstd_Dlam        = []
            fluxerr_Dlam        = []
            SNval_Dlam          = []
            fluxval_Dlam_max    = []
            SNval_Dlam_max      = []
            fluxval_Dlam_sum    = []
            SNval_Dlam_sum      = []
            # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^
            EW                  = []
            EW_Dlam             = []
            EW_Dlam_sum         = []

            if verbose: print(' - Looping over spectra for line = '+keys[ll])
            for ss, spec in enumerate(spectra):
                id       = spec.split('/')[-1].split('.fit')[0][-9:]
                objent   = np.where(sourcecat['id'] == int(id))[0]
                specdat  = afits.open(spec)[1].data
                spec_lam = specdat['wave']
                spec_f   = specdat['flux']
                spec_err = specdat['fluxerror']
                spec_s2n = specdat['s2n']

                if use_sys:
                    try:
                        z_lam    = sourcecat['z_sys_AV17'][objent] # using systemic redshift
                    except:
                        z_lam    = sourcecat['z_sys_V18'][objent] # using systemic redshift
                else:
                    z_lam    = sourcecat['z_vac_red'][objent]  # using Lya redshiftfrom red peak

                line_wave  = (z_lam+1)*line_lam
                lineinfo   = uves.lineinfofromspec(line_wave,spec_lam,spec_f,spec_err,spec_s2n,
                                                   deltalam=deltalam,verbose=verbose)

                ids.append(id)

                fluxval.append(lineinfo[0])
                fluxerr.append(lineinfo[1])
                SNval.append(lineinfo[2])
                fluxval_Dlam.append(lineinfo[3])
                fluxstd_Dlam.append(lineinfo[4])
                fluxerr_Dlam.append(lineinfo[5])
                SNval_Dlam.append(lineinfo[6])
                fluxval_Dlam_max.append(lineinfo[7])
                SNval_Dlam_max.append(lineinfo[8])
                fluxval_Dlam_sum.append(lineinfo[9])
                SNval_Dlam_sum.append(lineinfo[10])

                # - - - - -  Estimate EW  - - - - -
                beta          = sourcecat['beta'][objent][0]
                try:
                    bandswithLya  = uves.wavelength_in_bands( 1216 * (sourcecat['z_sys_AV17'][objent][0] + 1))
                except:
                    bandswithLya  = uves.wavelength_in_bands( 1216 * (sourcecat['z_sys_V18'][objent][0] + 1))
                bandswithline = uves.wavelength_in_bands( line_wave )

                f_cont        = []
                bands         = ['F606W','F775W','F814W','F125W']
                wave_refs     = [6034.0,7730.0,8140.2,12516.2]
                for bb, band in enumerate(bands):
                    if (band not in bandswithLya) & (band not in bandswithline):
                        if band == 'F125W':
                            f_ref  = sourcecat[ 'flux_wfc3_'+band.lower()[1:] ][objent][0] * 1e20
                        else:
                            f_ref  = sourcecat[ 'flux_acs_'+band.lower()[1:] ][objent][0] * 1e20
                        f_wave = uves.estimate_continuumlevel_viaBeta(line_wave, wave_refs[bb], f_ref, beta, verbose=False)
                        f_cont.append(f_wave)

                if (len(f_cont) > 0):
                    EW.append( fluxval[ss] / np.mean(np.asarray(f_cont)) )
                else:
                    EW.append( np.NaN )

                if (len(f_cont) > 0):
                    EW_Dlam.append( fluxval_Dlam[ss] / deltalam / np.mean(np.asarray(f_cont)) )
                else:
                    EW_Dlam.append( np.NaN )

                if (len(f_cont) > 0):
                    EW_Dlam_sum.append( fluxval_Dlam_sum[ss] / deltalam / np.mean(np.asarray(f_cont)) )
                else:
                    EW_Dlam_sum.append( np.NaN )

                # - - - - - - - - - - - - - - - - - -

            outputdic[keys[ll]] = ids, \
                                  fluxval, fluxerr, SNval, \
                                  fluxval_Dlam, fluxstd_Dlam, fluxerr_Dlam, SNval_Dlam, \
                                  fluxval_Dlam_max, SNval_Dlam_max, \
                                  fluxval_Dlam_sum, SNval_Dlam_sum, \
                                  EW, EW_Dlam, EW_Dlam_sum

    if printresults:
        for ss, spec in enumerate(spectra):
            id = spec.split('/')[-1].split('.fit')[0][-9:]
            print(' - - - - - - Object '+id+' (Dlam = lambda +/-'+str(outputdic['deltalam'])+'A)- - - - - -  ')
            for key in outputdic.keys():
                if key is not 'deltalam':
                    print(' - '+key+' flux           = '+str('%.4f' % outputdic[key][1][ss])+' +/- '+\
                          str('%.4f' % outputdic[key][2][ss]))
                    print(' - '+key+' S/N            = '+str('%.4f' % outputdic[key][3][ss]))
                    print(' - '+key+' flux_Dlam      = '+str('%.4f' % outputdic[key][4][ss])+' +/- '+\
                          str('%.4f' % outputdic[key][5][ss]))
                    print(' - '+key+' S/N_Dlam       = '+str('%.4f' % outputdic[key][7][ss]))
                    print(' - '+key+' flux_Dlam_max  = '+str('%.4f' % outputdic[key][8][ss]))
                    print(' - '+key+' S/N_Dlam_max   = '+str('%.4f' % outputdic[key][9][ss]))
                    print(' - '+key+' flux_Dlam_sum  = '+str('%.4f' % outputdic[key][10][ss]))
                    print(' - '+key+' S/N_Dlam_sum   = '+str('%.4f' % outputdic[key][11][ss]))
                    print(' - '+key+' EW             = '+str('%.4f' % outputdic[key][12][ss]))
                    print(' - '+key+' EW_Dlam        = '+str('%.4f' % outputdic[key][13][ss]))
                    print(' - '+key+' EW_Dlam_sum    = '+str('%.4f' % outputdic[key][14][ss]))

    return outputdic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def estimate_continuumlevel_viaBeta(wave, wave_ref, f_ref, beta, verbose=True):
    """
    Estiamting the continuum level assuming a spectrum with a fixed beta extrapolating ot the provided wavelength

    --- INPUT ---
    wave          The wavelength at whihc the continuum is estimated
    wave_ref      Referece wavelength at whihc the flux is known
    f_ref         Flux at wave_ref
    beta          The conitnuum slope to use for estrapolation to wave. Will assum f ~ wave**beta

    --- EXAMPLE OF USE ---

    sourcecatdat = afits.open('/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits')[1].data
    objent       = np.where(sourcecatdat['id'] == 102013086)[0]
    wave         = 1551 * (sourcecatdat['redshift'][objent][0] + 1)
    beta         = sourcecatdat['beta'][objent][0]
    bandswithLya = uves.wavelength_in_bands( 1216 * (sourcecatdat['redshift'][objent][0] + 1))

    wave_ref     = 8140.2
    f_ref        = sourcecatdat['flux_acs_814w'][objent][0]
    f_wave = uves.estimate_continuumlevel_viaBeta(wave, wave_ref, f_ref, beta)

    wave_ref     = 7730.0
    f_ref        = sourcecatdat['flux_acs_775w'][objent][0]
    f_wave = uves.estimate_continuumlevel_viaBeta(wave, wave_ref, f_ref, beta)

    wave_ref     = 6034.0
    f_ref        = sourcecatdat['flux_acs_606w'][objent][0]
    f_wave = uves.estimate_continuumlevel_viaBeta(wave, wave_ref, f_ref, beta)


    """
    if verbose: print(' - Estimating flux at wavelength '+str(wave)+' using:')
    if verbose: print('   reference wave = '+str(wave_ref))
    if verbose: print('   reference flux = '+str(f_ref))
    if verbose: print('   beta           = '+str(beta))
    if verbose: print('   ')

    f_wave = f_ref * (wave / wave_ref)**beta
    if verbose: print('   f_'+str(int(np.round(wave)))+'         = '+str(f_wave)+'\n')

    return f_wave
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_limits(sourcecatalog, namebase, limits_dictionary, colorcode=True, colortype='redshift', zoom = False,
                ignoreSNlt3=False, showids=False, verbose=True):
    """
    Plotting the output from uves.estimate_limits()

    """
    sourcedat = afits.open(sourcecatalog)[1].data

    for key in limits_dictionary.keys():
    #for key in ['lya']:
        if key == 'deltalam':
            continue

        ids, \
        fluxval, fluxerr, SNval, \
        fluxval_Dlam, fluxstd_Dlam, fluxerr_Dlam, SNval_Dlam, \
        fluxval_Dlam_max, SNval_Dlam_max, \
        fluxval_Dlam_sum, SNval_Dlam_sum, \
        EW, EW_Dlam, EW_Dlam_sum = limits_dictionary[key]

        # - - - - build vectors from source catalog using IDs - - - -
        z_sys      = []
        z_lya      = []
        EW_lya     = []
        EW_lya_err = []
        beta       = []

        for ii,id in enumerate(ids):
            objent = np.where(sourcedat['id'] == int(id))[0]
            try:
                z_sys.append(sourcedat['z_sys_AV17'][objent][0])
            except:
                z_sys.append(sourcedat['z_sys_V18'][objent][0])
            z_lya.append(sourcedat['z_vac_red'][objent][0])
            EW_lya.append(sourcedat['EW_0'][objent][0])
            EW_lya_err.append(sourcedat['EW_0_err'][objent][0])
            beta.append(sourcedat['beta'][objent][0])

        # - - - - - - - - - - - - - - - - - - - - - - PLOTTING - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print(' - Setting up and generating plot')
        plotname = namebase+'_'+key+'_EWlyaVSEW'+key+'.pdf'
        if zoom: plotname = plotname.replace('.pdf','_zoom.pdf')
        if ignoreSNlt3: plotname = plotname.replace('.pdf','_noSNlt3.pdf')
        fig = plt.figure(figsize=(7, 5))
        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.2, right=0.97, bottom=0.10, top=0.9)
        Fsize    = 10
        lthick   = 2
        marksize = 4
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        #plt.title(inforstr[:-2],fontsize=Fsize)

        xvalues = np.asarray(EW_lya)
        #yvalues = np.asarray(EW)
        #yvalues = np.asarray(EW_Dlam)
        yvalues = np.asarray(EW_Dlam_sum)
        xerr    = [None]*len(xvalues)
        yerr    = [None]*len(xvalues)

        if colorcode:
            cmap    = plt.cm.get_cmap('rainbow')

            if colortype == 'redshift':
                cmin    = 2.8
                cmax    = 6.2
            else:
                sys.exit(' Color type '+colortype+' not enabled ')

            colnorm = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
            cmaparr = np.linspace(cmin, cmax, num=50)
            m       = plt.cm.ScalarMappable(cmap=cmap)
            m.set_array(cmaparr)
            cb      = plt.colorbar(m)

            if colortype == 'redshift':
                cb.set_label('redshift')

            for ii,id in enumerate(ids):

                if colortype == 'redshift':
                    objcol = cmap(colnorm(z_sys[ii]))

                if np.isfinite(xvalues[ii]) & np.isfinite(yvalues[ii]):
                    marker = 'o'
                    ms     = marksize
                    mec    = objcol
                    mfc    = objcol

                    if SNval_Dlam_sum[ii] < 3.0:
                        mfc    = 'None'
                        if ignoreSNlt3: continue

                    if (beta[ii] == -2.0) & (SNval_Dlam_sum[ii] > 3.0):
                        #marker = r'$\downarrow$'
                        marker = r'$\nearrow$'
                        ms     = ms*2

                        plt.errorbar(xvalues[ii],yvalues[ii],xerr=xerr[ii],yerr=yerr[ii],
                                     marker=marker,lw=0, markersize=ms,alpha=1.0,
                                     markerfacecolor=mfc,ecolor=objcol,
                                     markeredgecolor=mec,zorder=10)
                    else:
                        if SNval_Dlam_sum[ii] > 0.0:
                            plt.errorbar(xvalues[ii],yvalues[ii],xerr=xerr[ii],yerr=yerr[ii],
                                         marker=marker,lw=0, markersize=ms,alpha=1.0,
                                         markerfacecolor=mfc,ecolor=objcol,
                                         markeredgecolor=mec,zorder=10)
                        else:
                            pass # not plotting S/N < 0 points

        else:
            plt.errorbar(xvalues,yvalues,xerr=xerr,yerr=yerr,
                         marker='o',lw=0, markersize=marksize,alpha=0.5,
                         markerfacecolor='gray',ecolor='k',
                         markeredgecolor='k',zorder=10)

        #marking AGN:
        AGN, AGNcand = uves.get_AGN_ids()
        # AGN     = ['104014050','115003085','214002011']
        # AGNcand = ['123048186','123501191','121033078']
        for ii,id in enumerate(ids):
            if np.isfinite(xvalues[ii]) & np.isfinite(yvalues[ii]) & (id in AGN):
                plt.errorbar(xvalues[ii],yvalues[ii],xerr=None,yerr=None,
                                 marker='*',lw=0, markersize=marksize*2,alpha=1.0,
                                 markerfacecolor='None',ecolor=objcol,
                                 markeredgecolor='black',zorder=20)

            if np.isfinite(xvalues[ii]) & np.isfinite(yvalues[ii]) & (id in AGNcand):
                plt.errorbar(xvalues[ii],yvalues[ii],xerr=None,yerr=None,
                                 marker='D',lw=0, markersize=marksize,alpha=1.0,
                                 markerfacecolor='None',ecolor=objcol,
                                 markeredgecolor='black',zorder=20)

        if showids:
            for ii,id in enumerate(ids):
                if np.isfinite(xvalues[ii]) & np.isfinite(yvalues[ii]):
                    plt.text(xvalues[ii],yvalues[ii],id,color='black',fontsize=Fsize/2.)

        if key.lower() == 'lya':
            plt.plot([-5000,5000],[-5000,5000],'--',color='gray',lw=lthick,zorder=5)
            plt.xscale('log')
            plt.yscale('log')

        plt.xlabel(' EW( Ly$\\alpha$ ) ')
        plt.ylabel(' EW( '+key+' )')

        #--------- RANGES ---------
        goodent = np.where(np.isfinite(xvalues) & np.isfinite(yvalues))
        xmin = np.min(xvalues[goodent])
        xmax = np.max(xvalues[goodent])
        dx   = xmax-xmin

        ymin = np.min(yvalues[goodent])
        ymax = np.max(yvalues[goodent])
        dy   = ymax-ymin

        plt.xlim([xmin-dx*0.05,xmax+dx*0.05])
        plt.ylim([ymin-dy*0.05,ymax+dy*0.05])
        if zoom:
            plt.xlim([0,400])
            plt.ylim([0,30])

        # if logx:
        #     plt.xscale('log')
        # if logy:
        #     plt.yscale('log')

        #--------- LEGEND ---------
        plt.errorbar(-5000,-5000,xerr=None,yerr=1,marker='o',lw=0, markersize=marksize,alpha=1.0,
                     markerfacecolor='k',ecolor='k',markeredgecolor='black',zorder=1,label='MW LAE (S/N\_'+key+' $>$ 3)')
        plt.errorbar(-5000,-5000,xerr=None,yerr=1,marker='o',lw=0, markersize=marksize,alpha=1.0,
                     markerfacecolor='None',ecolor='k',markeredgecolor='black',zorder=1,label='MW LAE (S/N\_'+key+' $<$ 3)')
        plt.errorbar(-5000,-5000,xerr=None,yerr=1,marker=r'$\nearrow$',lw=0, markersize=marksize*2,alpha=1.0,
                     markerfacecolor='None',ecolor='k',markeredgecolor='black',zorder=1,
                     label='MW LAE (HST non-det. lower limit)')


        plt.errorbar(-5000,-5000,xerr=None,yerr=None,marker='*',lw=0, markersize=marksize*2,alpha=1.0,
                     markerfacecolor='None',ecolor='None',markeredgecolor='black',zorder=1,label='AGN')
        plt.errorbar(-5000,-5000,xerr=None,yerr=None,marker='D',lw=0, markersize=marksize,alpha=1.0,
                     markerfacecolor='None',ecolor='None',markeredgecolor='black',zorder=1,label='AGN candidate')

        leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize/1.7},ncol=3,numpoints=1,
                         bbox_to_anchor=(0.5, 1.1),)  # add the legend
        leg.get_frame().set_alpha(0.7)
        #--------------------------

        if verbose: print('   Saving plot to',plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')


        # - - - - - - - - - - - - - - - - - - - - - - PLOTTING - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print(' - Setting up and generating plot')
        plotname = namebase+'_'+key+'_fluxVSs2n.pdf'
        if zoom: plotname = plotname.replace('.pdf','_zoom.pdf')
        fig = plt.figure(figsize=(7, 5))
        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.2, right=0.97, bottom=0.10, top=0.9)
        Fsize    = 10
        lthick   = 2
        marksize = 4
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        #plt.title(inforstr[:-2],fontsize=Fsize)

        xvalues = np.asarray(SNval_Dlam_sum)
        yvalues = np.asarray(fluxval)
        xerr    = None
        yerr    = np.asarray(fluxerr)

        if colorcode:
            cmap    = plt.cm.get_cmap('rainbow')

            if colortype == 'redshift':
                cmin    = 2.8
                cmax    = 6.2
            else:
                sys.exit(' Color type '+colortype+' not enabled ')

            colnorm = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
            cmaparr = np.linspace(cmin, cmax, num=50)
            m       = plt.cm.ScalarMappable(cmap=cmap)
            m.set_array(cmaparr)
            cb      = plt.colorbar(m)

            if colortype == 'redshift':
                cb.set_label('redshift')

            for ii,id in enumerate(ids):

                if colortype == 'redshift':
                    objcol = cmap(colnorm(z_sys[ii]))

                if np.isfinite(xvalues[ii]) & np.isfinite(yvalues[ii]):
                    plt.errorbar(xvalues[ii],yvalues[ii],xerr=xerr,yerr=yerr[ii],
                                 marker='o',lw=0, markersize=marksize,alpha=1.0,
                                 markerfacecolor=objcol,ecolor=objcol,
                                 markeredgecolor='None',zorder=10)
        else:
            plt.errorbar(xvalues,yvalues,xerr=xerr,yerr=yerr,
                         marker='o',lw=0, markersize=marksize,alpha=0.5,
                         markerfacecolor='gray',ecolor='k',
                         markeredgecolor='k',zorder=10)

        #marking AGN:
        AGN, AGNcand = uves.get_AGN_ids()
        # AGN     = ['104014050','115003085','214002011']
        # AGNcand = ['123048186','123501191','121033078']
        for ii,id in enumerate(ids):
            if np.isfinite(xvalues[ii]) & np.isfinite(yvalues[ii]) & (id in AGN):
                plt.errorbar(xvalues[ii],yvalues[ii],xerr=None,yerr=None,
                                 marker='*',lw=0, markersize=marksize*2,alpha=1.0,
                                 markerfacecolor='None',ecolor=objcol,
                                 markeredgecolor='black',zorder=20)

            if np.isfinite(xvalues[ii]) & np.isfinite(yvalues[ii]) & (id in AGNcand):
                plt.errorbar(xvalues[ii],yvalues[ii],xerr=None,yerr=None,
                                 marker='D',lw=0, markersize=marksize,alpha=1.0,
                                 markerfacecolor='None',ecolor=objcol,
                                 markeredgecolor='black',zorder=20)


        if showids:
            for ii,id in enumerate(ids):
                if np.isfinite(xvalues[ii]) & np.isfinite(yvalues[ii]):
                    plt.text(xvalues[ii],yvalues[ii],id,color='black',fontsize=Fsize/2.)

        plt.plot([3,3],[-5000,5000],'--',color='gray',lw=lthick,zorder=5)

        plt.xlabel('S/N of MUSE-Wide LAE 1D spectra at location of '+key)
        plt.ylabel('Flux [1e-20cgs] of MUSE-Wide LAE 1D spectra at location of '+key)


        #--------- RANGES ---------
        xmin = np.min(xvalues[np.isfinite(xvalues)])
        xmax = np.max(xvalues[np.isfinite(xvalues)])
        dx   = xmax-xmin

        ymin = np.min(yvalues[np.isfinite(yvalues)])
        ymax = np.max(yvalues[np.isfinite(yvalues)])
        dy   = ymax-ymin

        plt.xlim([xmin-dx*0.05,xmax+dx*0.05])
        plt.ylim([ymin-dy*0.05,ymax+dy*0.05])

        if zoom:
            plt.xlim([-10,10])
            plt.ylim([-100,100])
            if key.lower() == 'lya':
                plt.xlim([0,10])
                plt.ylim([0,300])

        # if logx:
        #     plt.xscale('log')
        # if logy:
        #     plt.yscale('log')

        #--------- LEGEND ---------
        plt.errorbar(-5000,-5000,xerr=None,yerr=1,marker='o',lw=0, markersize=marksize,alpha=1.0,
                     markerfacecolor='k',ecolor='k',markeredgecolor='black',zorder=1,label='MUSE-Wide LAE')
        plt.errorbar(-5000,-5000,xerr=None,yerr=None,marker='*',lw=0, markersize=marksize*2,alpha=1.0,
                     markerfacecolor='None',ecolor='None',markeredgecolor='black',zorder=1,label='AGN')
        plt.errorbar(-5000,-5000,xerr=None,yerr=None,marker='D',lw=0, markersize=marksize,alpha=1.0,
                     markerfacecolor='None',ecolor='None',markeredgecolor='black',zorder=1,label='AGN candidate')

        leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize/1.0},ncol=5,numpoints=1,
                         bbox_to_anchor=(0.5, 1.1),)  # add the legend
        leg.get_frame().set_alpha(0.7)
        #--------------------------

        if verbose: print('   Saving plot to '+plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')

        # # - - - - - - - - - - - - - - - - - - - - - - PLOTTING - - - - - - - - - - - - - - - - - - - - - -
        # if verbose: print(' - Setting up and generating plot'
        # plotname = namebase+'_'+key+'_LyaEWVSflux.pdf'
        # fig = plt.figure(figsize=(7, 5))
        # fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.2, right=0.97, bottom=0.10, top=0.9)
        # Fsize    = 10
        # lthick   = 2
        # marksize = 4
        # plt.rc('text', usetex=True)
        # plt.rc('font', family='serif',size=Fsize)
        # plt.rc('xtick', labelsize=Fsize)
        # plt.rc('ytick', labelsize=Fsize)
        # plt.clf()
        # plt.ioff()
        # #plt.title(inforstr[:-2],fontsize=Fsize)
        #
        # xvalues = np.asarray(LyaEW)
        # yvalues = np.asarray(fluxval)
        # xerr    = None
        # yerr    = np.asarray(fluxerr)
        #
        # if colorcode:
        #     cmap    = plt.cm.get_cmap('rainbow')
        #
        #     if colortype == 'redshift':
        #         cmin    = 2.8
        #         cmax    = 6.2
        #     else:
        #         sys.exit(' Color type '+colortype+' not enabled ')
        #
        #     colnorm = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
        #     cmaparr = np.linspace(cmin, cmax, num=50)
        #     m       = plt.cm.ScalarMappable(cmap=cmap)
        #     m.set_array(cmaparr)
        #     cb      = plt.colorbar(m)
        #
        #     if colortype == 'redshift':
        #         cb.set_label('redshift')
        #
        #     for ii,id in enumerate(ids):
        #
        #         if colortype == 'redshift':
        #             objcol = cmap(colnorm(z_sys[ii]))
        #
        #         if np.isfinite(xvalues[ii]) & np.isfinite(yvalues[ii]):
        #             plt.errorbar(xvalues[ii],yvalues[ii],xerr=xerr,yerr=yerr[ii],
        #                          marker='o',lw=0, markersize=marksize,alpha=1.0,
        #                          markerfacecolor=objcol,ecolor=objcol,
        #                          markeredgecolor='None',zorder=10)
        # else:
        #     plt.errorbar(xvalues,yvalues,xerr=xerr,yerr=yerr,
        #                  marker='o',lw=0, markersize=marksize,alpha=0.5,
        #                  markerfacecolor='gray',ecolor='k',
        #                  markeredgecolor='k',zorder=10)
        #
        # #marking AGN:
        # AGN, AGNcand = uves.get_AGN_ids()
        # #AGN     = ['104014050','115003085','214002011']
        # #AGNcand = ['123048186','123501191','121033078']
        # for ii,id in enumerate(ids):
        #     if np.isfinite(xvalues[ii]) & np.isfinite(yvalues[ii]) & (id in AGN):
        #         plt.errorbar(xvalues[ii],yvalues[ii],xerr=None,yerr=None,
        #                          marker='*',lw=0, markersize=marksize*2,alpha=1.0,
        #                          markerfacecolor='None',ecolor=objcol,
        #                          markeredgecolor='black',zorder=20)
        #
        #     if np.isfinite(xvalues[ii]) & np.isfinite(yvalues[ii]) & (id in AGNcand):
        #         plt.errorbar(xvalues[ii],yvalues[ii],xerr=None,yerr=None,
        #                          marker='D',lw=0, markersize=marksize,alpha=1.0,
        #                          markerfacecolor='None',ecolor=objcol,
        #                          markeredgecolor='black',zorder=20)
        #
        #
        # if showids:
        #     for ii,id in enumerate(ids):
        #         if np.isfinite(xvalues[ii]) & np.isfinite(yvalues[ii]):
        #             plt.text(xvalues[ii],yvalues[ii],id,color='black',fontsize=Fsize/2.,zorder=30)
        #
        # plt.plot([3,3],[-5000,5000],'--',color='gray',lw=lthick,zorder=5)
        #
        # plt.xlabel('S/N of MUSE-Wide LAE 1D spectra at location of '+key)
        # plt.ylabel('Flux [1e-20cgs] of MUSE-Wide LAE 1D spectra at location of '+key)
        #
        #
        # #--------- RANGES ---------
        # xmin = np.min(xvalues[np.isfinite(xvalues)])
        # xmax = np.max(xvalues[np.isfinite(xvalues)])
        # dx   = xmax-xmin
        #
        # ymin = np.min(yvalues[np.isfinite(yvalues)])
        # ymax = np.max(yvalues[np.isfinite(yvalues)])
        # dy   = ymax-ymin
        #
        # plt.xlim([xmin-dx*0.05,xmax+dx*0.05])
        # plt.ylim([ymin-dy*0.05,ymax+dy*0.05])
        #
        # # if logx:
        # #     plt.xscale('log')
        # # if logy:
        # #     plt.yscale('log')
        #
        # #--------- LEGEND ---------
        # plt.errorbar(-5000,-5000,xerr=None,yerr=1,marker='o',lw=0, markersize=marksize,alpha=1.0,
        #              markerfacecolor='k',ecolor='k',markeredgecolor='black',zorder=1,label='MUSE-Wide LAE')
        # plt.errorbar(-5000,-5000,xerr=None,yerr=None,marker='*',lw=0, markersize=marksize*2,alpha=1.0,
        #              markerfacecolor='None',ecolor='None',markeredgecolor='black',zorder=1,label='AGN')
        # plt.errorbar(-5000,-5000,xerr=None,yerr=None,marker='D',lw=0, markersize=marksize,alpha=1.0,
        #              markerfacecolor='None',ecolor='None',markeredgecolor='black',zorder=1,label='AGN candidate')
        #
        # leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize/1.0},ncol=5,numpoints=1,
        #                  bbox_to_anchor=(0.5, 1.1),)  # add the legend
        # leg.get_frame().set_alpha(0.7)
        # #--------------------------
        #
        # if verbose: print('   Saving plot to',plotname
        # plt.savefig(plotname)
        # plt.clf()
        # plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def lineinfofromspec(wavelength,spec_lam,spec_flux,spec_fluxerr,spec_s2n,deltalam=10,verbose=True):
    """
    return info at given wavelength based on spectrum

    """
    if (wavelength > np.min(spec_lam)) & (wavelength < np.max(spec_lam)):

        wavediff     = np.abs(spec_lam-wavelength)
        waveent      = np.where(wavediff == np.min(wavediff))

        if len(waveent) > 1:
            if verbose: print(' - multiple matches for '+str(wavelength)+'; returning info for the first match: '+
                              str(spec_lam[waveent[0]]))
        ent      = waveent[0]
        ent_dlam = np.where( (spec_lam > (spec_lam[ent]-deltalam)) & (spec_lam < (spec_lam[ent]+deltalam)) )

        fluxval       = spec_flux[ent]
        fluxerr       = spec_fluxerr[ent]
        SNval         = spec_s2n[ent]

        fluxval_Dlam  = np.median(spec_flux[ent_dlam])
        fluxstd_Dlam  = np.std(spec_fluxerr[ent_dlam])
        fluxerr_Dlam  = np.median(spec_fluxerr[ent_dlam])
        SNval_Dlam    = np.median(spec_s2n[ent_dlam])

        fluxval_Dlam_max  = np.max(spec_flux[ent_dlam])
        SNval_Dlam_max    = np.max(spec_s2n[ent_dlam])

        deltawave         = np.median(np.diff(spec_lam[ent_dlam]))
        fluxval_Dlam_sum  = np.sum(spec_flux[ent_dlam]) * deltawave
        SNval_Dlam_sum    = fluxval_Dlam_sum/fluxerr_Dlam

    else:
        if verbose: print(' - '+str(wavelength)+' not within spectral range; returning NaNs')
        fluxval       = [np.NaN]
        fluxerr       = [np.NaN]
        SNval         = [np.NaN]

        fluxval_Dlam  = np.NaN
        fluxstd_Dlam  = np.NaN
        fluxerr_Dlam  = np.NaN
        SNval_Dlam    = np.NaN

        fluxval_Dlam_max  = np.NaN
        SNval_Dlam_max    = np.NaN

        fluxval_Dlam_sum  = np.NaN
        SNval_Dlam_sum    = np.NaN

    return fluxval[0], fluxerr[0], SNval[0], \
           fluxval_Dlam, fluxstd_Dlam, fluxerr_Dlam, SNval_Dlam, \
           fluxval_Dlam_max, SNval_Dlam_max, \
           fluxval_Dlam_sum, SNval_Dlam_sum

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview_forsample(MUSEidlist,outputdir='./',yrangefullflux=[-400,1200],yrangefullSN=[-3,30],
                                  specdir  = '/Volumes/DATABCKUP1/TDOSEextractions/TDOSEext_171019/tdose_spectra/',
                                  tol3DHSTmatch=0.5,showPSFspec=True,clobber=False,verbose=True):
    """
    Wrapper to run mwp.plot_1DspecOverview() for a sample of objects collecting the relevant spectra

    --- INPUT ---
    MUSEidlist      List of MUSE ids of objects to plot
    outputdir       Directory to save figure to
    yrangefullflux  Yrange of full-spectra overview in flux figure
    yrangefullSN    Yrange of full-spectra overview in S/N figure
    specdir         Directory containing TDOSE spectra to plot
    tol3DHSTmatch   Tolerance of match to 3D-HST catalog (and spectra)
    showPSFspec     Plot the spectra extracted via PSF weighting (TU's extraction)
    clobber         Overwrite existing files?
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves

    MUSEidlist  = [103006046,119031070,208006149]
    MUSEidlist  = afits.open('/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits')[1].data['id']

    outputdir   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/171108_1DspecOverview/'
    plottedspec = uves.plot_1DspecOverview_forsample(MUSEidlist,outputdir=outputdir)

    """
    Nobj     = len(MUSEidlist)
    if verbose: print(' - Plotting spectra of '+str(Nobj)+' objects in "MUSEidlist"')
    LAEinfo        = afits.open('/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits')[1].data
    spectraplottet = collections.OrderedDict()

    if type(specdir) == str:
        specdir = [specdir]

    for oo, objID in enumerate(MUSEidlist):
        idstr        = str(objID)
        idstr_short  = str(objID)[:3]+str(objID)[-5:]

        outputfigure = outputdir+'id'+idstr+'_1DspecOverview.pdf'

        objent = np.where(LAEinfo['id'] == objID)[0]
        if len(objent) == 0:
            sys.exit(' object ID '+str(objID)+' not found in LAEinfo file ')

        sep3DHST = LAEinfo['sep_skelton'][objent]

        if (sep3DHST < tol3DHSTmatch) & (sep3DHST != 0.0):
            id3DHST = LAEinfo['id_skelton'][objent]
        else:
            id3DHST  = -99

        # - - - - - - - TDOSE extractions - - - - - - -
        TDOSEspecs = []
        TDOSElabel = []
        for sdir in specdir:
            TDOSEspec    = glob.glob(sdir+'tdose_spectrum_candels*'+idstr+'*.fits')+\
                           glob.glob(sdir+'tdose_spectrum_candels*'+idstr_short+'*.fits')
            if len(TDOSEspec) == 0:
                if verbose: print('----- No TDOSE spectrum found for '+str(objID)+' in \n      '+
                                  sdir+'\n      -> moving on to next object ----- ')
                continue

            labelstrings = ['TDOSE '+specname.split('/')[-1].split('_')[3] for specname in TDOSEspec]
            TDOSElabel   = TDOSElabel + labelstrings
            TDOSEspecs   = TDOSEspecs + TDOSEspec

        TDOSEwave    = ['wave']*len(TDOSEspecs)
        TDOSEflux    = ['flux']*len(TDOSEspecs)
        TDOSEferr    = ['fluxerror']*len(TDOSEspecs)
        TDOSEsky     = [None]*len(TDOSEspecs)
        TDOSEskyW    = [None]*len(TDOSEspecs)
        TDOSEskyF    = [None]*len(TDOSEspecs)
        TDOSEsky[0]  = '/Users/kschmidt/work/MUSE/skyspectra/SKY_SPECTRUM_'+LAEinfo['pointing'][objent][0]+'_av.fits'
        TDOSEskyW[0] = 'lambda'
        TDOSEskyF[0] = 'data'

        # - - - - - - - PSF weighted extraction - - - - - - -
        if showPSFspec:
            PSFextdir = '/Users/kschmidt/work/MUSE/spectra1D/Arche170127/spectra/'
            PSFspec   = glob.glob(PSFextdir+'spectrum_'+idstr+'.fits')+\
                        glob.glob(PSFextdir+'spectrum_'+idstr_short+'.fits')
            if len(PSFspec) > 1:
                sys.exit(' More than one PSF extraction spectrum found for '+idstr+' in '+PSFextdir)
            elif len(PSFspec) == 1:
                PSFlabel = ['PSFext']
                PSFwave  = ['WAVE_AIR']
                PSFflux  = ['FLUX']
                PSFferr  = ['FLUXERR']
                PSFsky   = [None]
                PSFskyW  = [None]
                PSFskyF  = [None]
            else:
                PSFlabel = []
                PSFwave  = []
                PSFflux  = []
                PSFferr  = []
                PSFsky   = []
                PSFskyW  = []
                PSFskyF  = []
        else:
            PSFspec  = []
            PSFlabel = []
            PSFwave  = []
            PSFflux  = []
            PSFferr  = []
            PSFsky   = []
            PSFskyW  = []
            PSFskyF  = []

        # - - - - - - - 3D-HST spectra - - - - - - -
        if id3DHST != -99:
            if idstr.startswith('1'):
                specdic = uves.get_3DHSTspecname([id3DHST[0]],field='goodss',spec1D=True,verbose=verbose)
            elif idstr.startswith('2'):
                specdic = uves.get_3DHSTspecname([id3DHST[0]],field='cosmos',spec1D=True,verbose=verbose)
            else:
                sys.exit(' Invalid MUSE id - does not start with "1" or "2"')

            grismspec  = specdic[str(int(id3DHST[0]))]

            if len(grismspec) > 0:
                grismlabel   = [gs.split('/')[-1].split('.')[0].replace('_','\_') for gs in grismspec]
                grismwave    = ['wave']*len(grismspec)
                grismflux    = ['flux']*len(grismspec)
                grismferr    = ['error']*len(grismspec)
                grismsky     = [None]*len(grismspec)
                grismskyW    = [None]*len(grismspec)
                grismskyF    = [None]*len(grismspec)
                grismsky[0]  = '/Users/kschmidt/work/MUSE/skytable.fits'
                grismskyW[0] = 'lam'
                grismskyF[0] = 'flux'

            else:
                grismlabel = []
                grismwave  = []
                grismflux  = []
                grismferr  = []
                grismsky   = []
                grismskyW  = []
                grismskyF  = []
        else:
            grismspec  = []
            grismlabel = []
            grismwave  = []
            grismflux  = []
            grismferr  = []
            grismsky   = []
            grismskyW  = []
            grismskyF  = []

        # - - - - - - - Assemble input for plot - - - - - - -
        spectra      = TDOSEspecs + PSFspec  + grismspec
        labels       = TDOSElabel + PSFlabel + grismlabel
        wavecols     = TDOSEwave  + PSFwave  + grismwave
        fluxcols     = TDOSEflux  + PSFflux  + grismflux
        fluxerrcols  = TDOSEferr  + PSFferr  + grismferr
        skyspectra   = TDOSEsky   + PSFsky   + grismsky
        wavecols_sky = TDOSEskyW  + PSFskyW  + grismskyW
        fluxcols_sky = TDOSEskyF  + PSFskyF  + grismskyF

        zLya         = LAEinfo['redshift'][objent][0]
        try:
            zsys         = LAEinfo['z_sys_AV17'][objent][0]
            voffset      = LAEinfo['red_peak_shift_AV17_kms'][objent][0]
        except:
            zsys         = LAEinfo['z_sys_V18'][objent][0]
            voffset      = LAEinfo['red_peak_shift_V18_kms'][objent][0]

        spectraplottet[idstr] = spectra
        if verbose:
            idno    = oo+1
            infostr = '----- Plot '+str(objID)+' at z = '+str("%.4f" % zLya)+' indicating voff = '+\
                      str("%.4f" % voffset)+'  ('+str(idno)+'/'+str(Nobj)+') -----                           '

            if (clobber == False) & (len(glob.glob(outputfigure.replace('.pdf','*.pdf'))) != 0):
                infostr = infostr.replace(') -----              ',
                                          ') -----> skip (clobber=False)')

                print(infostr)
            else:
                print(infostr)

                for plotSN in [True,False]:
                    if plotSN:
                        yrangefull = yrangefullSN
                    else:
                        yrangefull = yrangefullflux

                    mwp.plot_1DspecOverview(spectra, labels, wavecols, fluxcols, fluxerrcols, zLya, voffset=voffset,
                                            skyspectra=skyspectra, wavecols_sky=wavecols_sky, fluxcols_sky=fluxcols_sky,
                                            outputfigure=outputfigure, yrangefull=yrangefull, plotSN=plotSN,
                                            verbose=verbose)

    if verbose: print(' - Done...')
    return spectraplottet
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_3DHSTspecname(ids,spec1D=False,field='goodss',verbose=True):
    """

    Return absolute paths for existing 3D-HST 2D fits files. If

    --- INPUT ---
    ids           list of 3D-HST ids to return list of spectr for
    spec1D        if true, paths of the 1D spectra will be returned

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    spectradic = uves.get_3DHSTspecname([22684,21926,9033,7609])

    """
    path3DHST     = '/Volumes/DATABCKUP1/3DHST/'+field.lower()+'_WFC3_V4.1.5/'
    cat_spec      = path3DHST+field.lower()+'_3dhst_v4.1.5_catalogs/'+\
                    field.lower()+'_3dhst.v4.1.5.duplicates_2d.dat'    # list of (duplicate) spectra for each ID

    if field.lower() == 'goodss':
        dat_spec      = np.genfromtxt(cat_spec,comments='#',
                                      names=['id','s1','s2','s3','s4','s5','s6','s7'],
                                      dtype='d,40a,40a,40a,40a,40a,40a,40a')
    elif field.lower() == 'cosmos':
        dat_spec      = np.genfromtxt(cat_spec,comments='#',
                                      names=['id','s1','s2','s3','s4'],
                                      dtype='d,40a,40a,40a,40a')

    else:
        sys.exit(' No 3D-HST directory setup available for the field "'+field+'"')

    infodic = collections.OrderedDict()

    for id in ids:
        objent   = np.where(dat_spec['id'] == float(id))[0][0]
        namelist = [val for val in dat_spec[objent]][1:]
        filelist = []

        for name in namelist:
            if name != '00000':
                if spec1D:
                    namestr  = '/1D/FITS/'+name
                    filename = glob.glob(path3DHST+'/'+name[:-11]+namestr+'.1D.fits')
                    if os.path.isfile(filename[0]):
                        filelist.append(os.path.abspath(filename[0]))
                else:
                    # look for regular spectrum
                    namestr  = '/2D/FITS/'+name
                    filename = glob.glob(path3DHST+'/'+name[:-11]+namestr+'.2D.fits')
                    if os.path.isfile(filename[0]):
                        filelist.append(os.path.abspath(filename[0]))

                    # look for BIG spectrum
                    namestr  = '/BIG/2D/'+name.replace('G141','G141-big')
                    filename = glob.glob(path3DHST+'/'+name[:-11]+namestr+'.2D.fits')
                    if os.path.isfile(filename[0]):
                        filelist.append(os.path.abspath(filename[0]))

        infodic[str(int(id))] = filelist

    return infodic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def copy_singleobjsourcecats(outputdir='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_sourcecats_singleobjects/',
                             verbose=True):
    """

    --- EXAMPLE OF USE ---
    uves.copy_singleobjsourcecats()

    """
    modeldir      = '/Volumes/DATABCKUP1/TDOSEextractions/MW_LAEs_JKgalfitmodels/'

    sourcecats = ['model_acs_814w_candels-cdfs-01_cut_v1.0_id101005016_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-01_cut_v1.0_id101023043_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-01_cut_v1.0_id101024044_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-02_cut_v1.0_id102015088_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-03_cut_v1.0_id103006046_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-03_cut_v1.0_id103050126_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-04_cut_v1.0_id104024069_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-06_cut_v1.0_id106014046_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-06_cut_v1.0_id106035088_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-11_cut_v1.0_id111013028_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-22_cut_v1.0_id122002034_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-22_cut_v1.0_id122002035_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-23_cut_v1.0_id123016117_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-25_cut_v1.0_id125049122_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-28_cut_v1.0_id128038236_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-31_cut_v1.0_id131016105_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-31_cut_v1.0_id131016106_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-35_cut_v1.0_id135010177_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-35_cut_v1.0_id135010178_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-41_cut_v1.0_id141003075_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-41_cut_v1.0_id141036146_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-43_cut_v1.0_id143033113_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-44_cut_v1.0_id144008046_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-45_cut_v1.0_id145022065_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-45_cut_v1.0_id145034089_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-45_cut_v1.0_id145065132_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-46_cut_v1.0_id146053338_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cdfs-46_cut_v1.0_id146069355_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cosmos-01_cut_v1.0_id201073224_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cosmos-02_cut_v1.0_id202013030_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cosmos-02_cut_v1.0_id202044085_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cosmos-03_cut_v1.0_id203007099_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cosmos-04_cut_v1.0_id204053120_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cosmos-06_cut_v1.0_id206014089_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cosmos-07_cut_v1.0_id207022169_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cosmos-09_cut_v1.0_id209006108_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cosmos-11_cut_v1.0_id211015198_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cosmos-12_cut_v1.0_id212029067_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cosmos-13_cut_v1.0_id213022109_cutout2p0x2p0arcsec_sourcecatalog.*',
                  'model_acs_814w_candels-cosmos-15_cut_v1.0_id215016042_cutout2p0x2p0arcsec_sourcecatalog.*']

    for scat in sourcecats:
        cpout = commands.getoutput('cp '+modeldir+scat+' '+outputdir)
        if cpout != '':
            print(str(cpout))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_overviewdocument(outdir,outfile,clobber=False,
                         specoverview='/Volumes/DATABCKUP1/TDOSEextractions/171201_TDOSEextraction/overviewplots/idIIII*.pdf',
                         FoVoverview='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVoverviews/FoVoverview_IIII.png'):
    """
    Generate LaTeX document summarizing objects (using figures and text) to ease inspections

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    outdir  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/overviewdocument/'
    outfile = 'overviewdocument.tex'
    uves.gen_overviewdocument(outdir,outfile)

    """

    objdata   = afits.open('/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits')[1].data
    ids       = objdata['id']
    ras       = objdata['ra']
    decs      = objdata['dec']
    redshifts = objdata['redshift']
    pointings = objdata['pointing']

    # - - - - - Generate main document - - - - -
    if os.path.isfile(outdir+outfile) & (clobber==False):
        sys.exit('Document '+outdir+outfile+' exists and clobber = False')

    fmain = open(outdir+outfile,'w')
    fmain.write("""
\documentclass[a4paper,10pt]{article}
\\usepackage[latin1]{inputenc}
\\usepackage{float}
\\usepackage[pdftex]{graphicx}
\\usepackage[left=1cm,right=1cm,top=1cm,bottom=2cm]{geometry}
%===========================================================================
\\begin{document}
%===========================================================================
""")

    for ii, id in enumerate(ids[:]):
        pagename = 'page'+str(id)+'.tex'

        FoVfig   = FoVoverview.replace('IIII',str(id))
        specfigs = glob.glob(specoverview.replace('IIII',str(id)))
        fluxfig  = specfigs[0]
        snfig    = specfigs[1]

        fpage = open(outdir+pagename,'w')
        fpage.write("""
\section*{%s (%s)}
\small
\\begin{verbatim}
ID(MUSE-Wide) = %s     z  = %s     (ra,dec) = (%s,%s)
ID(Guo)       = %s     dr = %s     (ra,dec) = (%s,%s)
ID(Skelton)   = %s     dr = %s     (ra,dec) = (%s,%s)
\end{verbatim}
\\normalsize

\\begin{figure}[h]
\\begin{center}
\includegraphics[width=0.45\\textwidth]{%s}
\caption{Candels F814W 5$\\times$5 arcsec field of view. Circles have r=0.5arcsec. Red circles mark MUSE-Wide LAEs. White circles mark MUSE-Wide non-Ly$\\alpha$ EL sources (low-$z$).}
\label{fig:FoV%s}
\end{center}
\end{figure}

\\begin{figure*}
\\begin{center}
\includegraphics[width=0.9\\textwidth]{%s}
\includegraphics[width=0.9\\textwidth]{%s}
\caption{Spectral overview (incl. crossmatch to 3D-HST) with zoom-ins on rest-frame UV line regions. V-offset (gray region) estimated based on Verhamme+18 relations with peak seperation or FWHM.}
\label{fig:spec%s}
\end{center}
\end{figure*}

        """ % (id,pointings[ii],
               str("%11s" % int(id)),str("%.4f" % redshifts[ii]),ras[ii],decs[ii],
               str("%11s" % int(objdata['id_guo'][ii])),str("%.4f" % objdata['sep_guo'][ii]),objdata['ra_guo'][ii],objdata['dec_guo'][ii],
               str("%11s" % int(objdata['id_skelton'][ii])),str("%.4f" % objdata['sep_skelton'][ii]),objdata['ra_skelton'][ii],objdata['dec_skelton'][ii],
               FoVfig,id,fluxfig,snfig,id))
        fpage.close()

        fmain.write('\input{./'+pagename+'}\\newpage\n')

    fmain.write("""
%===========================================================================
\end{document}""")
    fmain.close()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_NEOGALmodels(limits_dictionary,plotnumbers=[99],cutSFmodels=False,smallsyms=False,verbose=True):
    """
    Commands to generate plots of the NEOGAL photo-ionization models with MUSE-Wide included

    --- EXAMPLE OF USE ---
    import rxj2248_BooneBalestraSource as bbs
    bbs.plot_feltregutkinmodels_cmds(plotnumbers=[1,2],cutSFmodels=False,smallsyms=True)

    """

    modeldata  = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/nebular_emission/')
    modeldata2 = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/AGN_NLR_nebular_feltre16/')

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

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (1 in plotnumbers) or (99 in plotnumbers):
        namebase  = './photomodels_CIVdCIIIvsCIVdHeII'
        line1     = 'CIV1550'
        line2     = 'CIII1908'
        line3     = 'CIV1550'
        line4     = 'HeII1640'
        xrange    = [2e-2,1e3]
        yrange    = [2e-4,1.5e2]
        boxranges = [0.7,1e99,0.5,1e99]

        bbs.plot_feltregutkinmodels(modeldata,line1=line1,line2=line2,line3=line3,line4=line4,verbose=verbose,
                                    varyparam='Zgas',logx=True,logy=True,logp1=True,
                                    fixxrange=xrange,fixyrange=yrange,cutSFmodels=cutSFmodels,
                                    plotname=namebase+'_Zgas.pdf', xid=xid, nh=nh, COratio=COratio, Mcutoff=Mcutoff,
                                    modeldata2=modeldata2,colormap='winter',boxranges=boxranges,legpos='lower right',
                                    smallsyms=smallsyms)

        bbs.plot_feltregutkinmodels(modeldata,line1=line1,line2=line2,line3=line3,line4=line4,verbose=verbose,
                                    varyparam='logUs',logx=True,logy=True,logp1=False,
                                    fixxrange=xrange,fixyrange=yrange,cutSFmodels=cutSFmodels,
                                    plotname=namebase+'_logUs.pdf', xid=xid, nh=nh, COratio=COratio, Mcutoff=Mcutoff,
                                    modeldata2=modeldata2,colormap='spring',boxranges=boxranges,legpos='lower right',
                                    smallsyms=smallsyms)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (2 in plotnumbers) or (99 in plotnumbers):
        line1='CIII1908'
        line2='HeII1640'
        line3='CIV1550'
        line4='HeII1640'
        xrange    = [1e-3,1e3]
        yrange    = [2e-4,1.5e2]
        boxranges = [1e-99,1e99,0.5,1e99]

        bbs.plot_feltregutkinmodels(modeldata,line1=line1,line2=line2,line3=line3,line4=line4,verbose=verbose,varyparam='Zgas',logx=True,logy=True,logp1=True,fixxrange=xrange,fixyrange=yrange,cutSFmodels=cutSFmodels,plotname='./photomodels_CIIIHeIIvsCIVdHeII_Zgas.pdf', xid=xid, nh=nh, COratio=COratio, Mcutoff=Mcutoff,modeldata2=modeldata2,colormap='winter',boxranges=boxranges,legpos='lower right',smallsyms=smallsyms)

        bbs.plot_feltregutkinmodels(modeldata,line1=line1,line2=line2,line3=line3,line4=line4,verbose=verbose,varyparam='logUs',logx=True,logy=True,logp1=False,fixxrange=xrange,fixyrange=yrange,cutSFmodels=cutSFmodels,plotname='./photomodels_CIIIHeIIvsCIVdHeII_logUs.pdf', xid=xid, nh=nh, COratio=COratio, Mcutoff=Mcutoff,modeldata2=modeldata2,colormap='spring',boxranges=boxranges,legpos='lower right',smallsyms=smallsyms)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (3 in plotnumbers) or (99 in plotnumbers):
        pass

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (4 in plotnumbers) or (99 in plotnumbers):
        line1='NV1240'
        line2='CIV1550'
        line3='CIV1550'
        line4='HeII1640'
        xrange    = [5e-5,1e1]
        yrange    = [2e-4,3e2]
        boxranges = [2.05,1e99,0.5,1e99]

        bbs.plot_feltregutkinmodels(modeldata,line1=line1,line2=line2,line3=line3,line4=line4,verbose=verbose,varyparam='Zgas',logx=True,logy=True,logp1=True,fixxrange=xrange,fixyrange=yrange,cutSFmodels=cutSFmodels,
plotname='./photomodels_NVCIVdvsCIVdHeII_Zgas.pdf', xid=xid, nh=nh, COratio=COratio, Mcutoff=Mcutoff,modeldata2=modeldata2,colormap='winter',boxranges=boxranges,legpos='lower right',smallsyms=smallsyms)

        bbs.plot_feltregutkinmodels(modeldata,line1=line1,line2=line2,line3=line3,line4=line4,verbose=verbose,varyparam='logUs',logx=True,logy=True,logp1=False,fixxrange=xrange,fixyrange=yrange,cutSFmodels=cutSFmodels,plotname='./photomodels_NVCIVdvsCIVdHeII_logUs.pdf', xid=xid, nh=nh, COratio=COratio, Mcutoff=Mcutoff,modeldata2=modeldata2,colormap='spring',boxranges=boxranges,legpos='lower right',smallsyms=smallsyms)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (5 in plotnumbers) or (99 in plotnumbers):
        line1='NV1240'
        line2='HeII1640'
        line3='CIV1550'
        line4='HeII1640'
        xrange    = [1e-6,1e3]
        yrange    = [1e-4,3e2]
        boxranges = [1.02,1e99,0.5,1e99]

        bbs.plot_feltregutkinmodels(modeldata,line1=line1,line2=line2,line3=line3,line4=line4,verbose=verbose,varyparam='Zgas',logx=True,logy=True,logp1=True,fixxrange=xrange,fixyrange=yrange,cutSFmodels=cutSFmodels,plotname='./photomodels_NVHeIIvsCIVdHeII_Zgas.pdf', xid=xid, nh=nh, COratio=COratio, Mcutoff=Mcutoff,modeldata2=modeldata2,colormap='winter',boxranges=boxranges,legpos='lower right',smallsyms=smallsyms)

        bbs.plot_feltregutkinmodels(modeldata,line1=line1,line2=line2,line3=line3,line4=line4,verbose=verbose,varyparam='logUs',logx=True,logy=True,logp1=False,fixxrange=xrange,fixyrange=yrange,cutSFmodels=cutSFmodels,plotname='./photomodels_NVHeIIvsCIVHeII_logUs.pdf', xid=xid, nh=nh, COratio=COratio, Mcutoff=Mcutoff,modeldata2=modeldata2,colormap='spring',boxranges=boxranges,legpos='lower right',smallsyms=smallsyms)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_singlelinetemplate(outfits='./felis_testing/uves_felis_template_singleline.fits',verbose=True):
    """
    Wrapper to generate spectral template with a single line

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    uves.gen_singlelinetemplate()
    """
    doubletlam  = [1907.0,1909.0]
    rangeDlam   = [np.min(doubletlam)-5.0,np.max(doubletlam)+15.0,0.1]
    tcdic = {}
    tcdic['LINE1'] = ['GAUSS', 1907.0, 0.5/4.5, 0.0, 10.0/4.5, 'CIII]1907A mimicking line']
    #tcdic['CONT']  = ['CONT', 1.0, 0.0, 1908.0,       'Continuum with flux 1.0 at 1908 + slope 0.0']
    valstring = '_Singleline_sig_0p5_flux_10p0'
    tempname = outfits.replace('.fits',valstring+'.fits')
    fbt.build_template(rangeDlam,tcdic,tempfile=tempname,overwrite=True)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_felistemplates(outfits='./uves_felis_template.fits',overwrite=False,addLSF=False,verbose=True):
    """
    Wrapper to generate spectral templates for cross-correlation with spectra.

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves

    outdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_testing/'
    uves.gen_felistemplates(outfits=outdir+'uves_felis_template.fits')

    """
    # - - - - - - - - - - - - - - - - - - LSF setup  - - - - - - - - - - - - - - - - - - -
    MUSELSFfwhm  = 1.0 #[A]
    MUSELSFsigma = MUSELSFfwhm/2.354
    LSFparam     = ['LSF', MUSELSFsigma,  'MUSE GAUSS LSF']

    # - - - - - - - - - - - - - - - - - - Width and raange setup - - - - - - - - - - - - - - - - - -
    dlam          = 0.05                                      # wavelength grid spacing in angstrom/pix
    tempwidth     = 10.0                                      # the wavelength range of template (half the width)
    sigmas        = np.array([ 0.2,  0.4,  0.6,  0.8,  1. ,  1.2]) # line width of emission lines in agnstrom
    sigmas_pix    = sigmas / dlam

    # - - - - - - - - - - - - - - - - - - CIII doublet - - - - - - - - - - - - - - - - - -
    doubletlam  = [1906.68,1908.73]
    rangeDlam   = [np.min(doubletlam)-tempwidth,np.max(doubletlam)+tempwidth,dlam]
    fluxCIII1   = 1.0
    fluxratios  = [0.1,0.4,0.7,1.0,1.3,1.6,1.9] # Osterbrock predicts CIII1/CIII2 < 1.6

    Ntemps      = len(fluxratios)*len(sigmas)
    if verbose: print(' - generating '+str(Ntemps)+' templates for the CIII doublet (varying sigma and flux ratio)')

    for fr in fluxratios:
        fluxCIII2 = fluxCIII1 / fr
        for sig in sigmas:
            tcdic = {}
            tcdic['CIII1']                 = ['GAUSS', doubletlam[0], sig, 0.0, fluxCIII1, 'CIII]1907A']
            tcdic['CIII2']                 = ['GAUSS', doubletlam[1], sig, 0.0, fluxCIII2,  'CIII]1909A']

            valstring = '_CIIIdoublet'+\
                        '_sig_'+str(sig).replace('.','p')+\
                        '_fluxratio_'+str(fr).replace('.','p')

            if addLSF:
                tcdic['LSF']        = LSFparam
                valstring = valstring+'_LSF_'+str(LSFparam[1]).replace('.','p')

            tempname = outfits.replace('.fits',valstring+'.fits')
            fbt.build_template(rangeDlam,tcdic,tempfile=tempname,overwrite=overwrite)

    # - - - - - - - - - - - - - - - - - - CIV doublet - - - - - - - - - - - - - - - - - -
    doubletlam  = [1548.195,1550.770]
    rangeDlam   = [np.min(doubletlam)-tempwidth,np.max(doubletlam)+tempwidth,dlam]
    fluxCIV1    = 1.0
    fluxratios  = [0.5,1.0,1.5,2.0,2.5,3.0] # Feibelman 1983 predicts CIV1/CIV2=2 from theory (Mainali+17 used CIV1/CIV2=1)

    Ntemps      = len(fluxratios)*len(sigmas)
    if verbose: print(' - generating '+str(Ntemps)+' templates for the CIV doublet (varying sigma and flux ratio)')

    for fr in fluxratios:
        fluxCIV2 = fluxCIV1 / fr
        for sig in sigmas:
            tcdic = {}
            tcdic['CIV1']                  = ['GAUSS', doubletlam[0], sig, 0.0, fluxCIV1, 'CIV1548A']
            tcdic['CIV2']                  = ['GAUSS', doubletlam[1], sig, 0.0, fluxCIV2,  'CIV1551A']

            valstring = '_CIVdoublet'+\
                        '_sig_'+str(sig).replace('.','p')+\
                        '_fluxratio_'+str(fr).replace('.','p')

            if addLSF:
                tcdic['LSF']        = LSFparam
                valstring = valstring+'_LSF_'+str(LSFparam[1]).replace('.','p')

            tempname = outfits.replace('.fits',valstring+'.fits')
            fbt.build_template(rangeDlam,tcdic,tempfile=tempname,overwrite=overwrite)

    # - - - - - - - - - - - - - - - - - - NV doublet - - - - - - - - - - - - - - - - - -
    doubletlam  = [1238.821,1242.804]
    rangeDlam   = [np.min(doubletlam)-tempwidth,np.max(doubletlam)+tempwidth,dlam]
    fluxNV1     = 1.0
    fluxratios  = [0.5,1.0,1.5,2.0,2.5,3.0] # Torres-Peimbert, S. & Pena, M. 1984; emissivity ratio of NV1/NV2~2 like CIV

    Ntemps      = len(fluxratios)*len(sigmas)
    if verbose: print(' - generating '+str(Ntemps)+' templates for the NV doublet (varying sigma and flux ratio)')

    for fr in fluxratios:
        fluxNV2 = fluxNV1 / fr
        for sig in sigmas:
            tcdic = {}
            tcdic['NV1']                   = ['GAUSS', doubletlam[0], sig, 0.0, fluxNV1, 'NV1239A']
            tcdic['NV2']                   = ['GAUSS', doubletlam[1], sig, 0.0, fluxNV2, 'NV1243A']

            valstring = '_NVdoublet'+\
                        '_sig_'+str(sig).replace('.','p')+\
                        '_fluxratio_'+str(fr).replace('.','p')

            if addLSF:
                tcdic['LSF']        = LSFparam
                valstring = valstring+'_LSF_'+str(LSFparam[1]).replace('.','p')

            tempname = outfits.replace('.fits',valstring+'.fits')
            fbt.build_template(rangeDlam,tcdic,tempfile=tempname,overwrite=overwrite)

    # - - - - - - - - - - - - - - - - - - OIII doublet - - - - - - - - - - - - - - - - - -
    doubletlam  = [1660.809,1666.150]
    rangeDlam   = [np.min(doubletlam)-tempwidth,np.max(doubletlam)+tempwidth,dlam]
    fluxOIII1   = 1.0
    fluxratios  = [0.1,0.3,0.5,0.7,0.9,1.1] # Morton1991tab2 OIII1/OIII2~0.5 (can be lower and higher but 1666 strongest, see eg. Mainali+17, Vanzella+16, )

    Ntemps      = len(fluxratios)*len(sigmas)
    if verbose: print(' - generating '+str(Ntemps)+' templates for the NV doublet (varying sigma and flux ratio)')

    for fr in fluxratios:
        fluxOIII2 = fluxOIII1 / fr
        for sig in sigmas:
            tcdic = {}
            tcdic['OIII1']  = ['GAUSS', doubletlam[0], sig, 0.0, fluxOIII1, 'OIII1661A']
            tcdic['OIII2']  = ['GAUSS', doubletlam[1], sig, 0.0, fluxOIII2, 'OIII1666A']

            valstring = '_OIIIdoublet'+\
                        '_sig_'+str(sig).replace('.','p')+\
                        '_fluxratio_'+str(fr).replace('.','p')

            if addLSF:
                tcdic['LSF']        = LSFparam
                valstring = valstring+'_LSF_'+str(LSFparam[1]).replace('.','p')

            tempname = outfits.replace('.fits',valstring+'.fits')
            fbt.build_template(rangeDlam,tcdic,tempfile=tempname,overwrite=overwrite)

    # - - - - - - - - - - - - - - - - - - HeII gauss - - - - - - - - - - - - - - - - - -
    linelam     = 1640.420
    rangeDlam   = [linelam-tempwidth,linelam+tempwidth,dlam]

    Ntemps      = len(sigmas)
    if verbose: print(' - generating '+str(Ntemps)+' templates for the HeII line (varying sigma)')

    for sig in sigmas:
        tcdic = {}
        tcdic['HeII']                  = ['GAUSS', linelam, sig, 0.0, 1.0, 'HeII1640A']

        valstring = '_HeII'+\
                    '_sig_'+str(sig).replace('.','p')

        if addLSF:
            tcdic['LSF']        = LSFparam
            valstring = valstring+'_LSF_'+str(LSFparam[1]).replace('.','p')

        tempname = outfits.replace('.fits',valstring+'.fits')
        fbt.build_template(rangeDlam,tcdic,tempfile=tempname,overwrite=overwrite)

    # - - - - - - - - - - - - - - - - - - SiIII doublet - - - - - - - - - - - - - - - - - -
    doubletlam  = [1882.71,1892.03] # See detection from Berg+19
    rangeDlam   = [np.min(doubletlam)-tempwidth,np.max(doubletlam)+tempwidth,dlam]
    fluxSiIII1  = 1.0
    fluxratios  = [0.1,0.4,0.7,1.0,1.3,1.6,1.9] # Osterbrock predicts SIII1/SIII2 < 1.7 similar to CIII doublet

    Ntemps      = len(fluxratios)*len(sigmas)
    if verbose: print(' - generating '+str(Ntemps)+' templates for the NV doublet (varying sigma and flux ratio)')

    for fr in fluxratios:
        fluxSiIII2 = fluxSiIII1 / fr
        for sig in sigmas:
            tcdic = {}
            tcdic['Si3_1']  = ['GAUSS', doubletlam[0], sig, 0.0, fluxSiIII1, 'SiIII1883A']
            tcdic['Si3_2']  = ['GAUSS', doubletlam[1], sig, 0.0, fluxSiIII2, 'SiIII1892A']

            valstring = '_SiIIIdoublet'+\
                        '_sig_'+str(sig).replace('.','p')+\
                        '_fluxratio_'+str(fr).replace('.','p')

            if addLSF:
                tcdic['LSF']        = LSFparam
                valstring = valstring+'_LSF_'+str(LSFparam[1]).replace('.','p')

            tempname = outfits.replace('.fits',valstring+'.fits')
            fbt.build_template(rangeDlam,tcdic,tempfile=tempname,overwrite=overwrite)


    # - - - - - - - - - - - - - - - - - - CII gauss - - - - - - - - - - - - - - - - - -
    linelam     = 1335.6627
    rangeDlam   = [linelam-tempwidth,linelam+tempwidth,dlam]

    Ntemps      = len(sigmas)
    if verbose: print(' - generating '+str(Ntemps)+' templates for the CII line (varying sigma)')

    for sig in sigmas:
        tcdic = {}
        tcdic['CII']                  = ['GAUSS', linelam, sig, 0.0, 1.0, 'CII1336A']

        valstring = '_CII'+\
                    '_sig_'+str(sig).replace('.','p')

        if addLSF:
            tcdic['LSF']        = LSFparam
            valstring = valstring+'_LSF_'+str(LSFparam[1]).replace('.','p')

        tempname = outfits.replace('.fits',valstring+'.fits')
        fbt.build_template(rangeDlam,tcdic,tempfile=tempname,overwrite=overwrite)

    # - - - - - - - - - - - - - - - - - - CII] gauss - - - - - - - - - - - - - - - - - -
    linelam     = 2326.113
    rangeDlam   = [linelam-tempwidth,linelam+tempwidth,dlam]

    Ntemps      = len(sigmas)
    if verbose: print(' - generating '+str(Ntemps)+' templates for the CII] line (varying sigma)')

    for sig in sigmas:
        tcdic = {}
        tcdic['CIIb']                  = ['GAUSS', linelam, sig, 0.0, 1.0, 'CII]2326A']

        valstring = '_CIIb'+\
                    '_sig_'+str(sig).replace('.','p')

        if addLSF:
            tcdic['LSF']        = LSFparam
            valstring = valstring+'_LSF_'+str(LSFparam[1]).replace('.','p')

        tempname = outfits.replace('.fits',valstring+'.fits')
        fbt.build_template(rangeDlam,tcdic,tempfile=tempname,overwrite=overwrite)

    # - - - - - - - - - - - - - - - - - - MgII doublet - - - - - - - - - - - - - - - - - -
    doubletlam  = [2795.528,2802.705]
    rangeDlam   = [np.min(doubletlam)-tempwidth,np.max(doubletlam)+tempwidth,dlam]
    fluxMgII1   = 1.0
    fluxratios  = [0.1,0.4,0.7,1.0,1.3,1.6,1.9] # ???

    Ntemps      = len(fluxratios)*len(sigmas)
    if verbose: print(' - generating '+str(Ntemps)+' templates for the NV doublet (varying sigma and flux ratio)')

    for fr in fluxratios:
        fluxMgII2 = fluxMgII1 / fr
        for sig in sigmas:
            tcdic = {}
            tcdic['MgII1']  = ['GAUSS', doubletlam[0], sig, 0.0, fluxMgII1, 'MgII2796A']
            tcdic['MgII2']  = ['GAUSS', doubletlam[1], sig, 0.0, fluxMgII2, 'MgII2803A']

            valstring = '_MgIIdoublet'+\
                        '_sig_'+str(sig).replace('.','p')+\
                        '_fluxratio_'+str(fr).replace('.','p')

            if addLSF:
                tcdic['LSF']        = LSFparam
                valstring = valstring+'_LSF_'+str(LSFparam[1]).replace('.','p')

            tempname = outfits.replace('.fits',valstring+'.fits')
            fbt.build_template(rangeDlam,tcdic,tempfile=tempname,overwrite=overwrite)


    # - - - - - - - - - - - - - - - - - - Width and raange setup - - - - - - - - - - - - - - - - - -
    #FWHM observed for MW LAEs <~10 which corresponds to 10/2.35 / (1+3) <~ 1.06 Angstrom rest-frame
    dlam          = 0.05                                       # wavelength grid spacing in angstrom/pix
    tempwidth     = 20.0                                       # the wavelength range of template (half the width)
    sigmas        = np.array([0.32,0.64,1.28,2.56,5.12,10.24]) # line width of emission lines in agnstrom (rest-frame)
    sigmas_pix    = sigmas/dlam
    # - - - - - - - - - - - - - - - - - - Lya skew gauss - - - - - - - - - - - - - - - - - -
    linelam     = 1215.6737
    rangeDlam   = [linelam-tempwidth,linelam+tempwidth,dlam]
    Lyaskew     = [0.0,3.0,6.0,9.0,12.0]

    Ntemps      = len(sigmas)*len(Lyaskew)
    if verbose: print(' - generating '+str(Ntemps)+' templates for the Lya line (varying sigma and skew)')

    for sig in sigmas:
        for skew in Lyaskew:
            tcdic = {}
            tcdic['Lya']                  = ['GAUSS', linelam, sig, skew, 1.0, 'Lya1216A']

            valstring = '_Lya'+\
                        '_sig_'+str(sig).replace('.','p')+\
                        '_skew_'+str(skew).replace('.','p')

            if addLSF:
                tcdic['LSF']        = LSFparam
                valstring = valstring+'_LSF_'+str(LSFparam[1]).replace('.','p')

            tempname = outfits.replace('.fits',valstring+'.fits')
            fbt.build_template(rangeDlam,tcdic,tempfile=tempname,overwrite=overwrite)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_felismockspec(outfits='./uves_felis_mock_MUSEspectrum.fits',redshift=3.5,
                      zoomxplot=None,verbose=True):
    """
    Wrapper to generate a set of mock MUSE spectra to test FELIS on

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves

    outdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_testing/'
    uves.gen_felismockspec(outfits=outdir+'uves_felis_mock_MUSEspectrum.fits',redshift=3.5,zoomxplot=np.array([1880,1940])*(1+3.5))

    """
    if verbose: print(' - Building mockspectra at redshift '+str(redshift)+' using FELIS tools')
    tcdic = {}
    #tcdic['CONT']                   = ['CONT',  1.0, 0.0, 0.0,          'Flat continuum at 1.0']

    voffsetCIV = -200 #km/s
    lam_obsCIV1, lam_offsetCIV1, dlamCIV1 = kbs.velocityoffset2dwave(redshift,voffsetCIV,1548.195,verbose=verbose)
    tcdic['CIV1']  = ['GAUSS', lam_offsetCIV1, 0.5, 0.0, 10.0, 'CIV1548A']

    lam_obsCIV2, lam_offsetCIV2, dlamCIV2 = kbs.velocityoffset2dwave(redshift,voffsetCIV,1550.770,verbose=verbose)
    tcdic['CIV2']  = ['GAUSS', lam_offsetCIV2, 0.5, 0.0, 5.0,  'CIV1551A']

    tcdic['CIII1'] = ['GAUSS', 1907.0 * (1+redshift), 0.5, 0.0, 10.0, 'CIII]1907A']

    tcdic['CIII2'] = ['GAUSS', 1909.0 * (1+redshift), 0.5, 0.0, 5.0,  'CIII]1909A']

    voffsetLYA = -500 #km/s
    lam_obsLYA, lam_offsetLYA, dlamLYA = kbs.velocityoffset2dwave(redshift,voffsetLYA,1216.0,verbose=verbose)
    tcdic['LYA']   = ['GAUSS', lam_offsetLYA, 20.0, 10.0, 100.0, 'Lya1216']

    noisesigmas = [0.05,0.5,1.0,3.0]
    for sigma in noisesigmas:
        mockspec = outfits.replace('.fits','_noisesigma'+str(sigma).replace('.','p')+'.fits')
        noise    = ['GAUSS', 0.0, sigma]
        fbt.build_template([4800,9300,0.25],tcdic,tempfile=mockspec,noise=noise,overwrite=True,zoomxplot=zoomxplot)
        if verbose: print(' - Wrote output spectrum to '+mockspec)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_felismockspec_fromsetupfile(specsetup,basename='./uves_felis_mock_spectrum_fromsetup.fits',
                                    plotspectra=False,overwrite=False,noisesetup='errspec',
                                    verbose=True,verbose_buildtemp=False):
    """
    Wrapper to generate a set of mock spectra to test FELIS on.
    The setup of each of the mock spectra are defined in the spec setup file

    Note; by setting noisesetup=None templates can be generated for the FELIS fitting using a similar setup file.

    --- INPUT ---
    specsetup           Setup file containing the components and spectra to generate
    basename            Main name (and dir) of mock fits spectra to generate
    plotspectra         To plot the generated spectra set this to true
    overwrite           Overwrite existing files?
    verbose             Toggle verbosity
    verbose_buildtemp   Toggle verbosity of the FELIS function building the template

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves

    basename   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra/uves_mock_spectrum_fromsetup.fits'
    specsetup  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_setup.txt'
    uves.gen_felismockspec_fromsetupfile(specsetup,basename=basename,verbose_buildtemp=True,overwrite=True,plotspectra=True)

    """
    if verbose: print(' - Building mockspectra using FELIS tools')
    templatedat = np.genfromtxt(specsetup,names=True,skip_header=7,comments='#',
                                dtype='80a,80a,80a,80a,80a,80a,80a,80a,80a,d,d,80a,80a')
    # ----- load error spectrum to use for noise simulations -----
    if noisesetup.lower() == 'errspec':
        skyspec     = None # should be noise spectrum '/Users/kschmidt/work/MUSE/skyspectra/SKY_SPECTRUM_candels-cdfs-36_av.fits'
        errspec     = '/Users/kschmidt/work/MUSE/spectra_noise/median_eff_noise_spectrum_70fields190819.fits'
        noisewave   = afits.open(errspec)[1].data['wave']
        noiseflux   = afits.open(errspec)[1].data['flux'] * 5.5 # 5.5 corresponds to r=0.6'' (30 pixel) aperture
                                                                # spectrum as total noise on such spectrum would
                                                                # be pix_noise x 30/sqrt(30)
        noise       = ['SPECTRUM', noisewave, noiseflux]
    elif noisesetup.lower() == 'fixedvalue':
        noisesigma  = 0.2
        noisemean   = 0.5
        noise       = ['GAUSS', noisemean, noisesigma]
    elif noisesetup.lower() == 'none':
        noise       = None
        noisesigma  = None

    for tt, tempname in enumerate(templatedat['namekey']):
        tempdic = {}
        linewaves  = np.asarray(templatedat['linewaves'][tt][1:-1].split(',')).astype(float)
        linesigmas = np.asarray(templatedat['sigmas'][tt][1:-1].split(',')).astype(float)
        lineskews  = np.asarray(templatedat['skew'][tt][1:-1].split(',')).astype(float)
        linefluxes = np.asarray(templatedat['scaling'][tt][1:-1].split(',')).astype(float)
        redshifts  = np.asarray(templatedat['redshift'][tt][1:-1].split(',')).astype(float)
        fratios    = np.asarray(templatedat['fratios'][tt][1:-1].split(',')).astype(float)

        for linesigma in linesigmas:
            for lineskew in lineskews:
                for lineflux in linefluxes:
                    for zval in redshifts:
                        for fratio in fratios:
                            tempdic['line1'] = [templatedat['type'][tt], linewaves[0] * (1+zval),
                                                linesigma, lineskew,
                                                lineflux, templatedat['headerinfo'][tt]+'_1']
                            if len(linewaves) == 2:
                                tempdic['line2'] = [templatedat['type'][tt], linewaves[1] * (1+zval),
                                                    linesigma, lineskew,
                                                    lineflux / fratio, templatedat['headerinfo'][tt]+'_2']
                                Ftotspec         = lineflux + lineflux / fratio
                            else:
                                Ftotspec         = lineflux

                            wavecen  = np.mean(linewaves)* (1+zval)
                            wavemin  = wavecen-templatedat['waverange'][tt]/2.* (1+zval)
                            wavemax  = wavecen+templatedat['waverange'][tt]/2.* (1+zval)

                            outstr   = '_'+tempname+\
                                       '_noisestd'+str(noisesigma).replace('.','p')+\
                                       '_sigma'+str("%.2f" % linesigma).replace('.','p')+\
                                       '_skew'+str("%.2f" % lineskew).replace('.','p')+\
                                       '_Ftot'+str("%.2f" % Ftotspec).replace('.','p')+\
                                       '_Fratio'+str("%.2f" % fratio).replace('.','p')+\
                                       '_z'+str("%.2f" % zval).replace('.','p')+\
                                       '.fits'
                            if noisesetup.lower() == 'errspec':
                                outstr = outstr.replace('_noisestd'+str(noisesigma).replace('.','p'),'_noisespec')
                            mockspec = basename.replace('.fits',outstr)
                            try:
                                fbt.build_template([wavemin,wavemax,templatedat['dwave'][tt]],
                                                   tempdic,tempfile=mockspec,
                                                   noise=noise,overwrite=overwrite,
                                                   plottemplate=plotspectra,zoomxplot=[wavemin,wavemax],
                                                   verbose=verbose_buildtemp,
                                                   waveunits=templatedat['waveunits'][tt],
                                                   fluxunits=templatedat['fluxunits'][tt])
                            except:
                                print('\n\n ERROR: fbt.build_template failed... stopping to enable invesitgations')
                                pdb.set_trace()

                            if verbose: print(' - Genreated spectrum:   '+mockspec)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_mockspeck_setup_parametertable(setupfile,skip_header=7,noisestr='_noisespec',basename='',verbose=True):
    """
    Function to write (and load) table with parameter sets of the individual spectra generated from setupfile with
    uves.gen_felismockspec_fromsetupfile(). Essentially just expanding the parameter lists.

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves

    basename   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra/uves_mock_spectrum_fromsetup.fits'
    specsetup  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_setup.txt'
    paramtable = uves.build_mockspeck_setup_parametertable(specsetup,basename=basename)

    """
    templatedat = np.genfromtxt(setupfile,names=True,skip_header=skip_header,
                                comments='#',dtype='80a,80a,80a,80a,80a,80a,80a,80a,80a,d,d,80a,80a')
    outputfile  = setupfile.replace('.txt','_parametertable.txt')

    fout = open(outputfile,'w')
    fout.write('# Parameter table of spectra generated with uves.gen_felismockspec_fromsetupfile() based on the setup\n')
    fout.write('# '+setupfile+'\n')
    fout.write('# \n')
    fout.write('#  linesigma     lineskew     lineflux     Ftotspec     redshift       Fratio       dwave              waveunits    fluxunits            specname \n')
    for tt, tempname in enumerate(templatedat['namekey']):
        linewaves  = np.asarray(templatedat['linewaves'][tt][1:-1].split(',')).astype(float)
        linesigmas = np.asarray(templatedat['sigmas'][tt][1:-1].split(',')).astype(float)
        lineskews  = np.asarray(templatedat['skew'][tt][1:-1].split(',')).astype(float)
        linefluxes = np.asarray(templatedat['scaling'][tt][1:-1].split(',')).astype(float)
        redshifts  = np.asarray(templatedat['redshift'][tt][1:-1].split(',')).astype(float)
        fratios    = np.asarray(templatedat['fratios'][tt][1:-1].split(',')).astype(float)
        dwave      = templatedat['dwave'][tt]
        fluxunits  = templatedat['fluxunits'][tt]
        waveunits  = templatedat['waveunits'][tt]

        for linesigma in linesigmas:
            for lineskew in lineskews:
                for lineflux in linefluxes:
                    for zval in redshifts:
                        for fratio in fratios:
                            if len(linewaves) == 2:
                                Ftotspec         = lineflux + lineflux/fratio
                            else:
                                Ftotspec         = lineflux

                            specext    = '_'+tempname+\
                                         noisestr+\
                                         '_sigma'+str("%.2f" % linesigma).replace('.','p')+\
                                         '_skew'+str("%.2f" % lineskew).replace('.','p')+\
                                         '_Ftot'+str("%.2f" % Ftotspec).replace('.','p')+\
                                         '_Fratio'+str("%.2f" % fratio).replace('.','p')+\
                                         '_z'+str("%.2f" % zval).replace('.','p')+\
                                         '.fits'
                            specname  = basename.replace('.fits',specext)

                            paramlist = str("%12.4f" % linesigma)+' '+\
                                        str("%12.4f" % lineskew)+' '+\
                                        str("%12.4f" % lineflux)+' '+\
                                        str("%12.4f" % Ftotspec)+' '+\
                                        str("%12.4f" % zval)+' '+\
                                        str("%12.4f" % fratio)+' '+\
                                        str("%12.4f" % dwave)+' '+\
                                        str("%20s" % waveunits)+' '+\
                                        str("%20s" % fluxunits)
                            fout.write(paramlist+'     '+specname+'\n')
    fout.close()

    paramtable = np.genfromtxt(outputfile,names=True,comments='#',skip_header=3,dtype='d,d,d,d,d,d,d,20a,20a,200a')
    return paramtable

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def pre190911_gen_mocspecFELISresults_summary(summaryfile,picklefiles,overwrite=False,verbose=True):
    """
    Generate a summary of the template characteristics FELIS determined to match the mock spectra
    the best, i.e. with the highest S/N.

    - - - - - - - - - - - - - - - - - - - - - -  NB - - - - - - - - - - - - - - - - - - - - - - -
    - - -  After correcting the FELIS normalization this function was depreciated on 190911 - - -
    - - -  Instead the summarizing can now be done with the new version of the function     - - -

    --- INPUT ---
    summaryfile        Path and name to summary file to generate.
    picklefiles        List of FELIS pickle files to summarize.
    overwrite          Overwrite the summary file if it already exists?
    verbose            Toggle verbosity

    --- EXAMPLE OF RUN ---
    import uvEmissionlineSearch as uves

    FELISoutputdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_CCresults/'
    summaryfile    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_CCresults_summary_all.txt'
    picklefiles    = glob.glob(FELISoutputdir+'*.pkl')
    summarydat     = uves.gen_mocspecFELISresults_summary(summaryfile,picklefiles)

    """
    if verbose: print(' - Generating a summary of the best-fit template in:\n   '+summaryfile)
    if os.path.isfile(summaryfile) & (not overwrite):
        sys.exit(' Summary file '+summaryfile+' already exists and overwrite=False ')

    fout = open(summaryfile,'w')
    fout.write('# Summary of '+str(len(picklefiles))+' FELIS pickle files provided \n')
    fout.write('# File contains the characteristics of the templates with max S/N from the FELIS template fits \n')
    fout.write('# The summary was generated with uves.compate_mockspec_to_FELISresults() on '+kbs.DandTstr2()+' \n')
    fout.write('# \n')
    fout.write('# Columns are:\n')
    fout.write('# z_spec                     Instrinspic redshift of matched mock spectrum \n')
    fout.write('# z_temp_S2Nmax              Estimated redshift from template match\n')
    fout.write('# sigma_spec_pix             Input sigma in pixels for mock spectrum (mean for multiple lines)\n')
    fout.write('# sigma_spec_ang_rf          Input sigma in rest-frame angstroms for mock spectrum (mean for multiple lines)\n')
    fout.write('# sigma_temp_pix             Input sigma in pixels for maxS/N template (mean for multiple lines)\n')
    fout.write('# sigma_temp_ang_rf          Input sigma in rest-frame angstroms for maxS/N template (mean for multiple lines)\n')
    fout.write('# Fratio_spec                Flux ratio of doublets in mock spectrum (line_lowwave/line_highwave) \n')
    fout.write('# Fratio_temp                Flux ratio of doublets in template      (line_lowwave/line_highwave) \n')
    fout.write('# Ftot_spec_intr             Intrinsic total flux of mock spectrum \n')
    fout.write('# Ftot_spec_trapz            Integreated total flux of mock spectrum after noise addition \n')
    fout.write('# Ftot_spec_trapz_err        Uncertainty on Ftot_spec_trapz (mock spectrum flux errors propogated) \n')
    fout.write('# Ftot_temp_trapz            Integreated total flux of normalized template scaled by fluxscale_S2Nmax\n')
    fout.write('# Ftot_temp_trapz_err        Uncertainty on Ftot_temp_trapz (mock spectrum flux errors propogated) \n')
    fout.write('# Ftot_temp_trapz_fsclaeerr  Uncertainty on Ftot_temp_trapz (using fluxscaleerr_S2Nmax for all pixels) \n')
    fout.write('# Ftot_temp_sum              Summed total flux of normalized template scaled by fluxscale_S2Nmax\n')
    fout.write('# Ftot_temp_sum_err          Uncertainty on Ftot_sum_trapz (mock spectrum flux errors propogated) \n')
    fout.write('# Ftot_temp_sum_fsclaeerr    Uncertainty on Ftot_sum_trapz (using fluxscaleerr_S2Nmax for all pixels) \n')
    fout.write('# vshift_spec                Known intrinsic velocity shift of mock spectrum \n')
    fout.write('# vshift_CCmatch             Estimated velocity shift from template match [ c*(z_spec-z_temp_S2Nmax)/(1+z_temp_S2Nmax) ]\n')
    fout.write('# fluxscale_S2Nmax           Flux scale applied to normalized template to obtain maxS/N match \n')
    fout.write('# fluxscaleerr_S2Nmax        Uncertainty on fluxscale_S2Nmax [sqrt(fluxscale_variance)]\n')
    fout.write('# S2Nmax                     The S/N value of the (scaled) template match to the mock spectrum \n')
    fout.write('# Ngoodent                   The number of good pixels used in the cross correlation \n')
    fout.write('# chi2                       Chi^2 value between the mock spectrum and the template match \n')
    fout.write('# lineS2N                    Estimated S/N of spectral feature within [lineS2Nwavemin,lineS2Nwavemin] \n')
    fout.write('# lineS2Nwavemin             Lower integration limit for S/N estimate \n')
    fout.write('# lineS2Nwavemax             Upper integration limit for S/N estimate \n')
    fout.write('# lineS2N_rf                 Estimated S/N (rest-frame) of spectral feature within [lineS2Nwavemin_rf,lineS2Nwavemin_rf] \n')
    fout.write('# lineS2Nwavemin_rf          Lower integration limit for rest-frame S/N estimate \n')
    fout.write('# lineS2Nwavemax_rf          Upper integration limit for rest-frame S/N estimate \n')
    fout.write('# Ftot                       Total flux of spectral feature (sum(f)*dwave) used to estimate line S/N \n')
    fout.write('# Ftot_sigma                 Squaroot of the variance/sqrt(Npix) of Ftot \n')
    fout.write('# spectrum                   The mock spectrum the templates were matched to \n')
    fout.write('# template                   The maxS/N template \n')
    fout.write('# \n')
    fout.write('# z_spec z_temp_S2Nmax    sigma_spec_pix sigma_spec_ang_rf sigma_temp_pix sigma_temp_ang_rf   '
               ' Fratio_spec Fratio_temp     '
               'Ftot_spec_intr Ftot_spec_trapz Ftot_spec_trapz_err '
               'Ftot_temp_trapz Ftot_temp_trapz_err Ftot_temp_trapz_fsclaeerr     '
               'Ftot_temp_sum Ftot_temp_sum_err Ftot_temp_sum_fsclaeerr '
               'vshift_spec vshift_CCmatch     fluxscale_S2Nmax fluxscaleerr_S2Nmax    '
               'S2Nmax Ngoodent chi2    '
               'lineS2N lineS2Nwavemin lineS2Nwavemax lineS2N_rf lineS2Nwavemin_rf lineS2Nwavemax_rf     '
               'Ftot Ftot_sigma     spectrum template \n')

    for pp, picklefile in enumerate(picklefiles):
        if verbose:
            infostr = ' - Summarizing picklefile  '+str("%.5d" % (pp+1))+' / '+str("%.5d" % len(picklefiles))+'     '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()
        loaddic       = felis.load_picklefile(picklefile)

        Nsigma_integration = 3.0

        for specname in loaddic.keys():
            keydic = loaddic[specname]

            # load info about max S/N template
            template, vshift_intr, vshift_match, fluxscale_S2Nmax, fluxscaleerr_S2Nmax, \
            S2Nmax, Ngoodent, chi2, z_spec, zS2Nmax =  \
                felis.getresult4maxS2N(loaddic,specname)

            # load matched spec and move to restframe
            s_wave   , s_flux   , s_df   , s_s2n    = felis.load_spectrum(specname,verbose=False)
            s_wave_rf, s_flux_rf, s_df_rf, s_s2n_rf = s_wave / (1+z_spec), s_flux * (1+z_spec), s_df * (1+z_spec), s_s2n

            # interpolate max S/N template to spec grid
            min_template_level = 1e-4
            t_wave_init, t_flux_init, t_df_init, t_s2n_init = felis.load_spectrum(template,verbose=False)
            func       = scipy.interpolate.interp1d(t_wave_init,t_flux_init,kind='linear',fill_value="extrapolate")
            t_flux     = func(s_wave_rf)
            t_flux[t_flux < min_template_level] = 0.0

            # Normalize template flux to 1
            if len(t_flux[t_flux != 0]) == 0:
                if verbose: print(' WARNING All interpolated template pixels are 0.0')
            else:
                temp_sum = np.sum(t_flux)
                t_flux   = t_flux / temp_sum

            # extract info on template and mockspec from fits headers
            spec_hdr        = afits.open(specname)[1].header
            spec_sigma_ang  = np.array([])
            spec_flux       = np.array([])
            spec_line_wave  = np.array([])
            for hdrkey in spec_hdr.keys():
                if ('noise' not in hdrkey.lower()) & ('err' not in hdrkey.lower()):
                    if '_1' in hdrkey: spec_line_wave = np.append(spec_line_wave,spec_hdr[hdrkey])
                    if '_2' in hdrkey: spec_sigma_ang = np.append(spec_sigma_ang,spec_hdr[hdrkey])
                    if '_4' in hdrkey: spec_flux      = np.append(spec_flux,spec_hdr[hdrkey])

            spec_sigma_ang_rf  = np.mean(spec_sigma_ang) / (1.0+z_spec)
            spec_sigma_pix     = spec_sigma_ang_rf       / np.median(np.diff(s_wave_rf))
            Ftot_spec_intr     = np.sum(spec_flux)
            if len(spec_flux) == 2:
                Fratio_spec = spec_flux[np.where(spec_line_wave == np.min(spec_line_wave))] / \
                              spec_flux[np.where(spec_line_wave == np.max(spec_line_wave))]
            else:
                Fratio_spec = 0.0


            lineS2Nwavemin = np.min(spec_line_wave)-Nsigma_integration*np.median(spec_sigma_ang)
            lineS2Nwavemax = np.max(spec_line_wave)+Nsigma_integration*np.median(spec_sigma_ang)
            waverange      = [lineS2Nwavemin,lineS2Nwavemax]
            goodent        = np.where((s_wave >= waverange[0]) & (s_wave <= waverange[1]))

            lineS2Nwavemin_rf = np.min(spec_line_wave/(1.0+z_spec))-Nsigma_integration*spec_sigma_ang_rf
            lineS2Nwavemax_rf = np.max(spec_line_wave/(1.0+z_spec))+Nsigma_integration*spec_sigma_ang_rf
            waverange_rf      = [lineS2Nwavemin_rf,lineS2Nwavemax_rf]
            goodent_rf        = np.where((s_wave_rf >= waverange_rf[0]) & (s_wave_rf <= waverange_rf[1]))

            # Ftot_trapz = np.trapz(fluxscale_S2Nmax * t_flux,s_wave_rf)
            datarr               = unumpy.uarray(fluxscale_S2Nmax * t_flux[goodent_rf], s_df_rf[goodent_rf])
            Ftot_trapz           = np.trapz(datarr,s_wave_rf[goodent_rf])
            Ftot_sum             = np.sum(datarr) * np.median(np.diff(s_wave_rf[goodent_rf]))

            datarr_fscaleerr     = unumpy.uarray(fluxscale_S2Nmax * t_flux[goodent_rf], fluxscaleerr_S2Nmax + t_flux[goodent_rf]*0.0)
            Ftot_trapz_fscaleerr = np.trapz(datarr_fscaleerr,s_wave_rf[goodent_rf])
            Ftot_sum_fscaleerr   = np.sum(datarr) * np.median(np.diff(s_wave_rf[goodent_rf]))

            datarr_spec        = unumpy.uarray(s_flux_rf[goodent_rf], s_df_rf[goodent_rf])
            Ftot_trapz_spec    = np.trapz(datarr_spec,s_wave_rf[goodent_rf])

            temp_hdr        = afits.open(template)[1].header
            temp_sigma_ang  = np.asarray([])
            temp_flux       = np.asarray([])
            temp_line_wave  = np.asarray([])
            for hdrkey in temp_hdr.keys():
                if ('noise' not in hdrkey.lower()) & ('err' not in hdrkey.lower()):
                    if '_1' in hdrkey: temp_line_wave = np.append(temp_line_wave,temp_hdr[hdrkey])
                    if '_2' in hdrkey: temp_sigma_ang = np.append(temp_sigma_ang,temp_hdr[hdrkey])
                    if '_4' in hdrkey: temp_flux      = np.append(temp_flux,temp_hdr[hdrkey])

            temp_sigma_ang_rf = np.mean(temp_sigma_ang) / (1.0+0.0)
            temp_sigma_pix    = np.mean(temp_sigma_ang) / np.median(np.diff(t_wave_init))
            if len(temp_flux) == 2:
                Fratio_temp = temp_flux[np.where(temp_line_wave == np.min(temp_line_wave))] / \
                              temp_flux[np.where(temp_line_wave == np.max(temp_line_wave))]
            else:
                Fratio_temp = 0.0

            #------------ Estimating S/N ------------
            Ftot, vartot, Npixgood, lineS2N = uves.calc_1Dspec_S2N(s_wave,s_flux,s_df**2.0,
                                                                   waverange,verbose=False)
            Ftot_sigma = np.sqrt(vartot)

            Ftot_rf, vartot_rf, Npixgood_rf, lineS2N_rf = uves.calc_1Dspec_S2N(s_wave_rf,s_flux_rf,s_df_rf**2.0,
                                                                               waverange_rf,verbose=False)

            if (Ftot-Ftot_rf) > 1.0:
                print(' - Ftot-Ftot_rf is less than 1e-20cgs; stopping to enable investigation')
                pdb.set_trace()

            #------------ Writing to output file ------------
            outstr = str("%7.4f" % z_spec)+'  '+\
                     str("%7.4f" % zS2Nmax)+'      '+\
                     str("%7.4f" % np.mean(spec_sigma_pix))+'  '+\
                     str("%7.4f" % spec_sigma_ang_rf)+'  '+\
                     str("%7.4f" % np.mean(temp_sigma_pix))+'  '+\
                     str("%7.4f" % temp_sigma_ang_rf)+'      '+\
                     str("%7.2f" % Fratio_spec)+'  '+\
                     str("%7.2f" % Fratio_temp)+'      '+\
                     str("%12.4f" % Ftot_spec_intr)+'  '+\
                     str("%12.4f" % Ftot_trapz_spec.nominal_value)+'  '+\
                     str("%12.4f" % Ftot_trapz_spec.std_dev)+'  '+\
                     str("%12.4f" % Ftot_trapz.nominal_value)+'  '+\
                     str("%12.4f" % Ftot_trapz.std_dev)+'  '+\
                     str("%12.4f" % Ftot_trapz_fscaleerr.std_dev)+'           '+\
                     str("%12.4f" % Ftot_sum.nominal_value)+'  '+\
                     str("%12.4f" % Ftot_sum.std_dev)+'  '+\
                     str("%12.4f" % Ftot_sum_fscaleerr.std_dev)+'           '+\
                     str("%12.4f" % vshift_intr)+'  '+\
                     str("%12.4f" % vshift_match)+'      '+\
                     str("%12.4f" % fluxscale_S2Nmax)+'  '+\
                     str("%12.4f" % fluxscaleerr_S2Nmax)+'  '+\
                     str("%12.4f" % S2Nmax)+'  '+\
                     str("%12.4f" % Ngoodent)+'  '+\
                     str("%12.4f" % chi2)+'  '+\
                     str("%12.4f" % lineS2N)+'  '+\
                     str("%12.4f" % lineS2Nwavemin)+'  '+\
                     str("%12.4f" % lineS2Nwavemax)+'  '+\
                     str("%12.4f" % lineS2N_rf)+'  '+\
                     str("%12.4f" % lineS2Nwavemin_rf)+'  '+\
                     str("%12.4f" % lineS2Nwavemax_rf)+'  '+\
                     str("%12.4f" % Ftot)+'  '+\
                     str("%12.4f" % Ftot_sigma)+'  '+\
                     specname+'  '+\
                     template+'  '
            fout.write(outstr+'\n')
    if verbose: print('\n   ...done')
    fout.close()

    fmt = 'd,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,300a,300a'
    summarydat = np.genfromtxt(summaryfile,skip_header=39,dtype=fmt,comments='#',names=True)
    return summarydat
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def pre190911_plot_mocspecFELISresults_summary(summaryfile,plotbasename,colortype='lineS2N_rf',histaxes=False,Nbins=50,
                                     overwrite=False,verbose=True):
    """
    plotting and evaluating the output from uves.gen_mocspecFELISresults_summary()

    - - - - - - - - - - - - - - - - - - - - - -  NB - - - - - - - - - - - - - - - - - - - - - - -
    - - -  After correcting the FELIS normalization this function was depreciated on 190911 - - -
    - - -  Instead the plotting can now be done with the new version of the function        - - -


    --- INPUT ---
    summaryfile        Path and name to summary file to evaluate
    plotbasename       The based name for the plots to generate (incl. output directory)
    overwrite          Overwrite the plots if they already exist?
    verbose            Toggle verbosity

    --- EXAMPLE OF RUN ---
    import uvEmissionlineSearch as uves

    summaryfile    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_CCresults_summary.txt'
    plotbasename   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_CCresults_summary_plots/190815test_'
    uves.plot_mocspecFELISresults_summary(summaryfile,plotbasename)

    """
    if verbose: print(' - Loading and plotting the content of \n   '+summaryfile)
    fmt = 'd,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,300a,300a'
    summarydat = np.genfromtxt(summaryfile,skip_header=39,dtype=fmt,comments='#',names=True)
    specnumber = np.arange(len(summarydat))+1.0

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'S2NvsS2N'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['lineS2N_rf'])
    yvalues  = np.asarray(summarydat['S2Nmax'])
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = 'Rest-frame emission line S/N'
    ylabel   = 'Template match S/N'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,colortype='Ftot',
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[1,1000],yrange=[1,1000],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'LineSigma'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['sigma_spec_ang_rf'])
    yvalues  = np.asarray(summarydat['sigma_temp_ang_rf'])
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = 'Mock spectrum $\sigma_\\textrm{Gauss, restframe}$ [\AA]'
    ylabel   = 'Best-fit template $\sigma_\\textrm{Gauss, restframe}$ [\AA]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='onetoone',
                                                   colortype='s2n',cdatvec = summarydat['lineS2N_rf'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'Ftot_intrinsic'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['Ftot_spec_intr'])
    yvalues  = np.asarray(summarydat['fluxscale_S2Nmax'])
    xerr     = [None]*len(xvalues)
    yerr     = summarydat['fluxscaleerr_S2Nmax']
    xlabel   = 'Intrinsic (pre-noise) flux in mock spectrum [1e-20 erg/s/cm$^2$]'
    ylabel   = '$\\alpha$ FELIS template flux estimate [1e-20 erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colortype=colortype,colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'Ftot_observed'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['Ftot_spec_trapz'])
    yvalues  = np.asarray(summarydat['fluxscale_S2Nmax'])
    xerr     = summarydat['Ftot_spec_trapz_err']
    yerr     = summarydat['fluxscaleerr_S2Nmax']
    xlabel   = 'Observed (post-noise) flux in mock spectrum [1e-20 erg/s/cm$^2$]'
    ylabel   = '$\\alpha$ FELIS template flux estimate [1e-20 erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colortype=colortype,colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'Ftot_intrinsic_sum'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['Ftot_spec_intr'])
    yvalues  = np.asarray(summarydat['Ftot_temp_sum'])
    xerr     = [None]*len(xvalues)
    yerr     = summarydat['Ftot_temp_sum_err']
    xlabel   = 'Total flux in mock spectrum (intrinsic)\\\\(no noise) [1e-20 erg/s/cm$^2$]'
    ylabel   = 'Total flux in best-fit template (sum)\\\\(mock spec noise) [1e-20 erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colortype='s2n',cdatvec = summarydat['lineS2N_rf'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'Ftot_observed_sum'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['Ftot'])
    yvalues  = np.asarray(summarydat['Ftot_temp_sum'])
    xerr     = summarydat['Ftot_sigma']
    yerr     = summarydat['Ftot_temp_sum_fsclaeerr']
    xlabel   = 'Total flux in mock spectrum (sum)\\\\("observed"; with noise) [1e-20 erg/s/cm$^2$]'
    ylabel   = 'Total flux in best-fit template (sum)\\\\(integrated; flux scale "noise") [1e-20 erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colortype='s2n',cdatvec = summarydat['lineS2N_rf'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'Redshift'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['z_spec'])
    yvalues  = np.asarray(summarydat['z_temp_S2Nmax'])
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = 'Intrinsic redshift mock spectrum'
    ylabel   = 'Redshift of best-fit template'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='onetoone',
                                                   colortype='vshift',colorcode=True,cdatvec = summarydat['vshift_CCmatch'],
                                                   overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'dFtot_vs_S2N'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['S2Nmax'])
    yvalues  = np.asarray(summarydat['Ftot_spec_intr']) - np.asarray(summarydat['Ftot_temp_sum'])
    xerr     = [None]*len(xvalues)
    yerr     = summarydat['Ftot_temp_sum_err']
    xlabel   = 'S/N of template cross match to mock spectrum'
    ylabel   = '$\Delta$Total flux; intrinsic mock spec - best-fit template [1e-20 erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='horizontal',
                                                   # yrange=[-3,3],colortype='redshift',
                                                   yrange=None,colortype='redshift',cdatvec = summarydat['z_spec'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'dsigma_vs_S2N'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['S2Nmax'])
    yvalues  = np.asarray(summarydat['sigma_spec_ang_rf']) - np.asarray(summarydat['sigma_temp_ang_rf'])
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(yvalues)
    xlabel   = 'S/N of template cross match to mock spectrum'
    ylabel   = '$\Delta\sigma_\\textrm{Gauss, restframe}$; intrinsic mock spec - best-fit template [\AA]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='horizontal',
                                                   # yrange=[-0.05,0.05],colortype='redshift',
                                                   yrange=None,colortype='redshift',cdatvec = summarydat['z_spec'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'dsigma_vs_sigmaspec'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['sigma_spec_ang_rf'])
    yvalues  = np.asarray(summarydat['sigma_spec_ang_rf']) - np.asarray(summarydat['sigma_temp_ang_rf'])
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(yvalues)
    xlabel   = 'Mock spec intrinsic $\sigma_\\textrm{Gauss, restframe}$ [\AA]'
    ylabel   = '$\Delta\sigma_\\textrm{Gauss, restframe}$; intrinsic mock spec - best-fit template [\AA]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='horizontal',
                                                   # yrange=[-0.05,0.05],colortype=colortype,
                                                   yrange=None,colortype='s2n',cdatvec = summarydat['lineS2N_rf'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'dFtot_vs_specno'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = specnumber
    yvalues  = np.asarray(summarydat['Ftot_spec_intr']) - np.asarray(summarydat['Ftot_temp_sum'])
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = 'Spectrum number - according to summary file \n'+summaryfile.split('/')[-1].replace('_','\_')
    ylabel   = '$\Delta$Ftot; mock spec - temp match [1e-20erg/s/cm$^2$]'


    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='horizontal',colortype='Sigma',cdatvec = summarydat['sigma_spec_ang_rf'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'dFtot_vs_Ftot'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['Ftot_spec_intr'])
    yvalues  = np.asarray(summarydat['Ftot_spec_intr']) - np.asarray(summarydat['Ftot_temp_sum'])
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = 'Ftot mock spec [1e-20erg/s/cm$^2$]'
    ylabel   = '$\Delta$Ftot; mock spec - temp match [1e-20erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='horizontal',colortype='Sigma',cdatvec = summarydat['sigma_spec_ang_rf'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'Ratio_Ftot_vs_Ftot'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['Ftot_spec_intr'])
    yvalues  = (np.asarray(summarydat['Ftot_temp_sum'])/np.asarray(summarydat['Ftot_spec_intr'])) -1.0
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = 'Ftot mock spec [1e-20erg/s/cm$^2$]'
    ylabel   = '(Ftot temp match / Ftot mock spec) - 1 '

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='horizontal',colortype='Sigma',cdatvec = summarydat['sigma_spec_ang_rf'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # nameext  = 'Ratio_Ftot_vs_lineS2N'
    # plotname = plotbasename+nameext+'.pdf'
    # xvalues  = np.asarray(summarydat['lineS2N_rf'])
    # yvalues  = (np.asarray(summarydat['Ftot_temp_trapz'])/np.asarray(summarydat['Ftot_spec_intr'])) -1.0
    # xerr     = [None]*len(xvalues)
    # yerr     = [None]*len(xvalues)
    # xlabel   = 'Line S/N'
    # ylabel   = '(Ftot temp match / Ftot mock spec) - 1 '
    #
    # uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
    #                                                histaxes=histaxes,Nbins=Nbins,
    #                                                linetype='horizontal',colortype='Ftot_spec_intr',
    #                                                colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'Ratio_Ftot_vs_lineS2N_sum_temp'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['lineS2N_rf'])
    yvalues  = (np.asarray(summarydat['Ftot_temp_sum'])/np.asarray(summarydat['Ftot_spec_intr'])) -1.0
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = 'Line S/N'
    ylabel   = '(Ftot sum tempalte / Ftot intrinsic) - 1 '

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,yrange=[-0.7,0.7],
                                                   linetype='horizontal',colortype='Ftot_spec_intr',
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'Ratio_Ftot_vs_lineS2N_sum_spec'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['lineS2N_rf'])
    yvalues  = (np.asarray(summarydat['Ftot'])/np.asarray(summarydat['Ftot_spec_intr'])) -1.0
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = 'Line S/N'
    ylabel   = '(Ftot sum mock spec w. noise / Ftot mock spec intrinsic) - 1 '

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,yrange=[-0.7,0.7],
                                                   linetype='horizontal',colortype='Ftot_spec_intr',
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'dsigma_vs_specno'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = specnumber
    yvalues  = np.asarray(summarydat['sigma_spec_ang_rf']) - np.asarray(summarydat['sigma_temp_ang_rf'])
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = 'Spectrum number - according to summary file \n'+summaryfile.split('/')[-1].replace('_','\_')
    ylabel   = '$\Delta\sigma_\\textrm{rest}$; mock spec - temp match [\AA]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='horizontal',
                                                   colortype='s2n',cdatvec = summarydat['lineS2N_rf'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'dsigma_vs_sigma'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['sigma_spec_ang_rf'])
    yvalues  = np.asarray(summarydat['sigma_spec_ang_rf']) - np.asarray(summarydat['sigma_temp_ang_rf'])
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = '$\sigma_\\textrm{rest}$ mock spec'
    ylabel   = '$\Delta\sigma_\\textrm{rest}$; mock spec - temp match [\AA]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='horizontal',
                                                   colortype='s2n',cdatvec = summarydat['lineS2N_rf'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'Ratio_sigma_vs_sigma'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = np.asarray(summarydat['sigma_spec_ang_rf'])
    yvalues  = (np.asarray(summarydat['sigma_temp_ang_rf'])/np.asarray(summarydat['sigma_spec_ang_rf'])) -1.0
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = '$\sigma_\\textrm{rest}$ mock spec'
    ylabel   = '($\sigma_\\textrm{rest}$ temp match / $\sigma_\\textrm{rest}$ mock spec) - 1'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='horizontal',
                                                   colortype='s2n',cdatvec = summarydat['lineS2N_rf'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'dz_vs_specno'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = specnumber
    yvalues  = np.asarray(summarydat['z_spec'])-np.asarray(summarydat['z_temp_S2Nmax'])
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = 'Spectrum number - according to summary file \n'+summaryfile.split('/')[-1].replace('_','\_')
    ylabel   = '$\Delta z$; mock spec - temp match'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='horizontal',colortype='redshift',cdatvec = summarydat['z_spec'],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    #-------------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------------
    goodFratio = np.where(summarydat['Fratio_spec'] > 0)
    if len(goodFratio[0]) > 0:
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'FluxRatio'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = np.asarray(summarydat['Fratio_spec'][goodFratio])
        yvalues  = np.asarray(summarydat['Fratio_temp'][goodFratio])
        xerr     = [None]*len(xvalues)
        yerr     = [None]*len(xvalues)
        xlabel   = 'Flux ratio mock spectrum doublet lines'
        ylabel   = 'Flux ratio best-fit template doublet lines'

        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       linetype='onetoone',
                                                       colortype='s2n',cdatvec = summarydat['lineS2N_rf'][goodFratio],
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'dFluxRatio_vs_specno'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = specnumber[goodFratio]
        yvalues  = np.asarray(summarydat['Fratio_spec'][goodFratio])-np.asarray(summarydat['Fratio_temp'][goodFratio])
        xerr     = [None]*len(xvalues)
        yerr     = [None]*len(xvalues)
        xlabel   = 'Spectrum number - according to summary file \n'+summaryfile.split('/')[-1].replace('_','\_')
        ylabel   = '$\Delta$Flux ratio; mock spec- temp match'
        cdatvec  = summarydat['sigma_spec_ang_rf'][goodFratio]

        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       linetype='horizontal',colortype='Sigma',cdatvec=cdatvec,
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'dFluxRatio_vs_Fluxratio'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = np.asarray(summarydat['Fratio_spec'][goodFratio])
        yvalues  = np.asarray(summarydat['Fratio_spec'][goodFratio])-np.asarray(summarydat['Fratio_temp'][goodFratio])
        xerr     = [None]*len(xvalues)
        yerr     = [None]*len(xvalues)
        xlabel   = 'Flux ratio mock spec'
        ylabel   = '$\Delta$Flux ratio; mock spec - temp match'
        cdatvec  = summarydat['sigma_spec_ang_rf'][goodFratio]

        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       linetype='horizontal',colortype='Sigma',cdatvec=cdatvec,
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'Ratio_FluxRatio_vs_Fluxratio'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = np.asarray(summarydat['Fratio_spec'][goodFratio])
        yvalues  = (np.asarray(summarydat['Fratio_temp'][goodFratio])/np.asarray(summarydat['Fratio_spec'][goodFratio])) - 1.0
        xerr     = [None]*len(xvalues)
        yerr     = [None]*len(xvalues)
        xlabel   = 'Flux ratio mock spec '
        ylabel   = '(Flux ratio temp match / Flux ratio mock spec) - 1'
        cdatvec  = summarydat['sigma_spec_ang_rf'][goodFratio]

        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       linetype='horizontal',colortype='Sigma',cdatvec=cdatvec,
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'Ratio_FluxRatio_vs_specno'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = specnumber[goodFratio]
        yvalues  = (np.asarray(summarydat['Fratio_temp'][goodFratio])/np.asarray(summarydat['Fratio_spec'][goodFratio])) - 1.0
        xerr     = [None]*len(xvalues)
        yerr     = [None]*len(xvalues)
        xlabel   = 'Spectrum number - according to summary file \n'+summaryfile.split('/')[-1].replace('_','\_')
        ylabel   = '(Flux ratio temp match / Flux ratio mock spec) - 1'
        cdatvec  = summarydat['sigma_spec_ang_rf'][goodFratio]

        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       linetype='horizontal',colortype='Sigma',cdatvec=cdatvec,
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'Ratio_FluxRatio_vs_specno'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = specnumber[goodFratio]
        yvalues  = (np.asarray(summarydat['Fratio_temp'][goodFratio])/np.asarray(summarydat['Fratio_spec'][goodFratio])) - 1.0
        xerr     = [None]*len(xvalues)
        yerr     = [None]*len(xvalues)
        xlabel   = 'Spectrum number - according to summary file \n'+summaryfile.split('/')[-1].replace('_','\_')
        ylabel   = '(Flux ratio temp match / Flux ratio mock spec) - 1'
        cdatvec  = summarydat['sigma_spec_ang_rf'][goodFratio]

        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       linetype='horizontal',colortype='Sigma',cdatvec=cdatvec,
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'Ratio_FluxRatio_vs_lineS2N'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = summarydat['lineS2N_rf'][goodFratio]
        yvalues  = (summarydat['Fratio_temp'][goodFratio]/summarydat['Fratio_spec'][goodFratio]) - 1.0
        xerr     = [None]*len(xvalues)
        yerr     = [None]*len(xvalues)
        xlabel   = 'Line S/N'
        ylabel   = '(Flux ratio temp match / Flux ratio mock spec) - 1'
        cdatvec  = summarydat['sigma_spec_ang_rf'][goodFratio]
        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       linetype='horizontal',colortype='Sigma',
                                                       cdatvec=cdatvec,
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'Ratio_FluxRatio_vs_Ratio_sigma'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = (summarydat['sigma_temp_ang_rf'][goodFratio]/summarydat['sigma_spec_ang_rf'][goodFratio]) -1.0
        yvalues  = (summarydat['Fratio_temp'][goodFratio]/summarydat['Fratio_spec'][goodFratio]) - 1.0
        xerr     = [None]*len(xvalues)
        yerr     = [None]*len(xvalues)
        xlabel   = '($\sigma_\\textrm{rest}$ temp match / $\sigma_\\textrm{rest}$ mock spec) - 1'
        ylabel   = '(Flux ratio temp match / Flux ratio mock spec) - 1'
        cdatvec  = summarydat['lineS2N_rf'][goodFratio]
        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       linetype='plus',
                                                       colortype='s2n',cdatvec=cdatvec,
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)


        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'Ratio_Ftot_vs_Ratio_FluxRatio'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = (summarydat['Fratio_temp'][goodFratio]/summarydat['Fratio_spec'][goodFratio]) - 1.0
        yvalues  = (summarydat['Ftot_temp_sum'][goodFratio]/summarydat['Ftot_spec_intr'][goodFratio]) - 1.0
        xerr     = [None]*len(xvalues)
        yerr     = [None]*len(xvalues)
        xlabel   = '(Flux ratio temp match / Flux ratio mock spec) - 1'
        ylabel   = '(Ftot temp match / Ftot mock spec) - 1 '
        cdatvec  = summarydat['sigma_spec_ang_rf'][goodFratio]
        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       linetype='plus',colortype='Sigma',
                                                       cdatvec=cdatvec,
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'Ratio_Ftot_vs_Ratio_sigma'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = (summarydat['sigma_temp_ang_rf']/summarydat['sigma_spec_ang_rf']) - 1.0
    yvalues  = (summarydat['Ftot_temp_sum']/summarydat['Ftot_spec_intr']) - 1.0
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = '($\sigma_\\textrm{rest}$ temp match / $\sigma_\\textrm{rest}$ mock spec) - 1'
    ylabel   = '(Ftot temp match / Ftot mock spec) - 1 '
    cdatvec  = summarydat['lineS2N_rf']
    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   linetype='plus',colortype='s2n',
                                                   cdatvec=cdatvec,
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def match_mockspectra_to_templates(outputdir,CCwavewindow=25.0,plot_allCCresults=False,noisefree=False,
                                   spec2match='all',verbose=True):
    """
    Wrapper around match_mockspectrum_to_templates() to match templates to all mockspectra.


    --- INPUT ---
    outputdir          Directory to store picklefiles with cross-correlation results to
    CCwavewindow       Window around line to perform cross-correlations over (rest frame)
    plot_allCCresults  To plot all CC results set this to True. Plots will be stored in outputdir
    noisefree          To matcht the noise-free version of the mock spectra set to True
    spec2match         To only match a a subsample of spectra provide their parameter tabel line index here
    verbose            Toggle verbosity

    --- EXAMPLE OF RUN ---
    import uvEmissionlineSearch as uves

    outputdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_CCresults/'
    uves.match_mockspectra_to_templates(outputdir)

    """
    if verbose: print(' - Decide which templates to fit to what mock spectra (through a dictionary)')
    mockVStemp_lines = {}
    mockVStemp_lines['Lya']             = ['Lya']
    mockVStemp_lines['CIIIdoublet']     = ['CIII']
    mockVStemp_lines['CIVdoublet']      = ['CIV']
    mockVStemp_lines['NVdoublet']       = ['NV']
    mockVStemp_lines['HeII']            = ['HeII']
    mockVStemp_lines['OIII1663doublet'] = ['OIII']
    mockVStemp_lines['MgIIdoublet']     = ['MgII']

    mockVStemp_lines['testsinglet'] = ['HeII']
    mockVStemp_lines['testdoublet'] = ['CIII']

    waverest     = {'Lya':1215.67, 'CIII':1908.0, 'CIV':1550.0, 'NV':1241.0, 'HeII':1640.0, 'OIII':1663.0, 'MgII':2800.}

    if verbose: print(' - Defining files to perform matches between (hardcoded)')
    parentdir    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'
    # templatedir  = parentdir+'felis_templates_190816/'
    templatedir  = parentdir+'felis_templates_fromsetup/'
    paramfile    = parentdir+'mockspectra_setup_parametertable.txt'
    paramtable   = np.genfromtxt(paramfile,names=True,comments='#',skip_header=3,dtype='d,d,d,d,d,d,d,20a,20a,200a')

    if verbose: print(' - Crosscorrelating templates to spectra using FELIS')
    for ss, mockspec in enumerate(paramtable['specname']):
        if spec2match != 'all':
            if ss not in spec2match:
                continue
            else:
                print('------ spec2match: following spectrum in list so matching:')
                print('                   '+mockspec)

        if noisefree:
            mockspec = mockspec.replace('noisespec','noisestdNone')
        mockline  = mockspec.split('fromsetup_')[-1].split('_')[0]
        templines = mockVStemp_lines[mockline]

        for templine in templines:
            picklefile = outputdir+mockspec.split('/')[-1].replace('.fits',
                                                                   '_CCresults_template'+templine+
                                                                   '_matchto_spectrum'+mockline+'.pkl')
            templates  = glob.glob(templatedir+'uves_felis_template_*'+templine+'*.fits')

            if not plot_allCCresults:
                plotdir  = None

            windowcen  = CCwavewindow * (1+paramtable['redshift'][ss])
            ccdic      = felis.match_templates2specs(templates,[mockspec],[paramtable['redshift'][ss]],
                                                     picklefile,wavewindow=[windowcen],
                                                     plotdir=outputdir,wavecen_restframe=[waverest[templine]],
                                                     vshift=None,min_template_level=1e-4,
                                                     plot_allCCresults=plot_allCCresults,
                                                     subtract_spec_median=False)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_mocspecFELISresults_summary(summaryfile,picklefiles,overwrite=False,verbose=True):
    """
    Generate a summary of the template characteristics FELIS determined to match the mock spectra
    the best, i.e. with the highest S/N.

    --- INPUT ---
    summaryfile        Path and name to summary file to generate.
    picklefiles        List of FELIS pickle files to summarize.
    overwrite          Overwrite the summary file if it already exists?
    verbose            Toggle verbosity

    --- EXAMPLE OF RUN ---
    import uvEmissionlineSearch as uves

    FELISoutputdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_CCresults/'
    summaryfile    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_CCresults_summary_all.txt'
    picklefiles    = glob.glob(FELISoutputdir+'*.pkl')
    summarydat     = uves.gen_mocspecFELISresults_summary(summaryfile,picklefiles)

    """
    if verbose: print(' - Generating a summary of the best-fit template in:\n   '+summaryfile)
    if os.path.isfile(summaryfile) & (not overwrite):
        sys.exit(' Summary file '+summaryfile+' already exists and overwrite=False ')

    fout = open(summaryfile,'w')
    fout.write('# Summary of '+str(len(picklefiles))+' FELIS pickle files provided \n')
    fout.write('# File contains the characteristics of the templates with max S/N from the FELIS template fits \n')
    fout.write('# The summary was generated with uves.compate_mockspec_to_FELISresults() on '+kbs.DandTstr2()+' \n')
    fout.write('# \n')
    fout.write('# Columns are:\n')
    fout.write('# z_spec                     Instrinspic redshift of matched mock spectrum \n')
    fout.write('# z_temp_S2Nmax              Estimated redshift from template match\n')
    fout.write('# sigma_spec_ang_obs         Observed line width in angstroms for mock spectrum \n')
    fout.write('# sigma_spec_ang_rf          Rest-frame line width in angstroms for mock spectrum \n')
    fout.write('# sigma_temp_ang_rf          Rest-frame line width in angstroms for maxS/N template \n')
    fout.write('# Fratio_spec                Flux ratio of doublets in mock spectrum (line_lowwave/line_highwave) \n')
    fout.write('# Fratio_temp                Flux ratio of doublets in template      (line_lowwave/line_highwave) \n')
    fout.write('# Ftot_spec_intr             Intrinsic (input) total flux of mock spectrum \n')
    fout.write('# Ftot_spec_trapz            Integreated (trapz) total flux of noisy mock spectrum \n')
    fout.write('# Ftot_spec_trapz_err        Uncertainty on Ftot_spec_trapz (mock spectrum flux errors propogated) \n')
    fout.write('# Ftot_spec_sum              Integreated (sum*dwave) total flux of noisy mock spectrum \n')
    fout.write('# Ftot_spec_sum_err          Integreated (sum*dwave) total flux of noisy mock spectrum \n')
    fout.write('# Ftot_FELIS_S2Nmax          Estimated total flux from FELIS template match for maxS/N match \n')
    fout.write('# Ftot_FELIS_S2Nmax_err      Uncertainty on FELISflux_S2Nmax [sqrt(Ftot_FELIS_S2Nmax_variance)]\n')
    fout.write('# FELIS_S2Nmax               The S/N value of the (scaled) template match to the mock spectrum \n')
    fout.write('# Ngoodent                   The number of good pixels used in the cross correlation \n')
    fout.write('# chi2                       Chi^2 value between the mock spectrum and the template match \n')
    fout.write('# vshift_spec                Known intrinsic velocity shift of mock spectrum \n')
    fout.write('# vshift_CCmatch             Estimated velocity shift from template match '
               ' [ c*(z_spec-z_temp_S2Nmax)/(1+z_temp_S2Nmax) ]\n')
    fout.write('# lineS2N                    Estimated S/N of spectral feature within [lineS2Nwavemin,lineS2Nwavemin] \n')
    fout.write('# lineS2Nwavemin             Lower integration limit for S/N estimate \n')
    fout.write('# lineS2Nwavemax             Upper integration limit for S/N estimate \n')
    fout.write('# lineS2N_rf                 Estimated S/N (rest-frame) of spectral feature within [lineS2Nwavemin_rf,lineS2Nwavemin_rf] \n')
    fout.write('# lineS2Nwavemin_rf          Lower integration limit for rest-frame S/N estimate \n')
    fout.write('# lineS2Nwavemax_rf          Upper integration limit for rest-frame S/N estimate \n')
    fout.write('# Ftot_lineS2N               Total flux of spectral feature (sum(f)*dwave) used to estimate line S/N \n')
    fout.write('# Ftot_lineS2N_sigma         Square root of the variance/sqrt(Npix) of Ftot \n')
    fout.write('# spectrum                   The mock spectrum the templates were matched to \n')
    fout.write('# template                   The maxS/N template \n')
    fout.write('# \n')
    fout.write('# z_spec z_temp_S2Nmax sigma_spec_ang_obs sigma_spec_ang_rf sigma_temp_ang_rf Fratio_spec Fratio_temp Ftot_spec_intr Ftot_spec_trapz Ftot_spec_trapz_err Ftot_spec_sum Ftot_spec_sum_err Ftot_FELIS_S2Nmax Ftot_FELIS_S2Nmax_err FELIS_S2Nmax Ngoodent chi2 vshift_spec vshift_CCmatch lineS2N lineS2Nwavemin lineS2Nwavemax lineS2N_rf lineS2Nwavemin_rf lineS2Nwavemax_rf Ftot_lineS2N Ftot_lineS2N_sigma spectrum template \n')

    for pp, picklefile in enumerate(picklefiles):
        if verbose:
            infostr = ' - Summarizing picklefile  '+str("%.5d" % (pp+1))+' / '+str("%.5d" % len(picklefiles))+'     '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()
        pkldic  = felis.load_picklefile(picklefile)

        Nsigma_integration = 3.0

        for specname in pkldic.keys():
            tempdic = pkldic[specname]

            #------ load info about max S/N template ------
            template, vshift_intr, vshift_match, Ftot_FELIS_S2Nmax, Ftot_FELIS_S2Nmax_err, \
            FELIS_S2Nmax, Ngoodent, chi2, z_spec, zS2Nmax =  \
                felis.getresult4maxS2N(pkldic,specname)

            #------ load matched spec and move to restframe ------
            s_wave   , s_flux   , s_df   , s_s2n    = felis.load_spectrum(specname,verbose=False)
            s_wave_rf, s_flux_rf, s_df_rf, s_s2n_rf = s_wave / (1+z_spec), s_flux * (1+z_spec), s_df * (1+z_spec), s_s2n

            #------ extract info on mock spectrum from fits headers ------
            spec_hdr        = afits.open(specname)[1].header
            spec_sigma_ang  = np.array([])
            spec_flux       = np.array([])
            spec_line_wave  = np.array([])
            for hdrkey in spec_hdr.keys():
                if ('noise' not in hdrkey.lower()) & ('err' not in hdrkey.lower()):
                    if '_1' in hdrkey: spec_line_wave = np.append(spec_line_wave,spec_hdr[hdrkey])
                    if '_2' in hdrkey: spec_sigma_ang = np.append(spec_sigma_ang,spec_hdr[hdrkey])
                    if '_4' in hdrkey: spec_flux      = np.append(spec_flux,spec_hdr[hdrkey])

            spec_sigma_ang_obs = np.mean(spec_sigma_ang)
            spec_sigma_ang_rf  = spec_sigma_ang_obs / (1.0+z_spec)
            Ftot_spec_intr     = np.sum(spec_flux)
            if len(spec_flux) == 2:
                Fratio_spec = spec_flux[np.where(spec_line_wave == np.min(spec_line_wave))] / \
                              spec_flux[np.where(spec_line_wave == np.max(spec_line_wave))]
            else:
                Fratio_spec = 0.0

            #------ extract info on template from fits headers ------
            temp_hdr        = afits.open(template)[1].header
            temp_sigma_ang  = np.array([])
            temp_flux       = np.array([])
            temp_line_wave  = np.array([])
            for hdrkey in temp_hdr.keys():
                if ('noise' not in hdrkey.lower()) & ('err' not in hdrkey.lower()):
                    if '_1' in hdrkey: temp_line_wave = np.append(temp_line_wave,temp_hdr[hdrkey])
                    if '_2' in hdrkey: temp_sigma_ang = np.append(temp_sigma_ang,temp_hdr[hdrkey])
                    if '_4' in hdrkey: temp_flux      = np.append(temp_flux,temp_hdr[hdrkey])

            temp_sigma_ang_rf  = np.mean(temp_sigma_ang)
            if len(temp_flux) == 2:
                Fratio_temp = temp_flux[np.where(temp_line_wave == np.min(temp_line_wave))] / \
                              temp_flux[np.where(temp_line_wave == np.max(temp_line_wave))]
            else:
                Fratio_temp = 0.0

            #------ integrate mock spectrum to obtain observed total fluxes (rest-frame) ------
            lineS2Nwavemin = np.min(spec_line_wave)-Nsigma_integration*spec_sigma_ang_obs
            lineS2Nwavemax = np.max(spec_line_wave)+Nsigma_integration*spec_sigma_ang_obs
            waverange      = [lineS2Nwavemin,lineS2Nwavemax]
            goodent        = np.where((s_wave >= waverange[0]) & (s_wave <= waverange[1]))

            lineS2Nwavemin_rf = np.min(spec_line_wave/(1.0+z_spec))-Nsigma_integration*spec_sigma_ang_rf
            lineS2Nwavemax_rf = np.max(spec_line_wave/(1.0+z_spec))+Nsigma_integration*spec_sigma_ang_rf
            waverange_rf      = [lineS2Nwavemin_rf,lineS2Nwavemax_rf]
            goodent_rf        = np.where((s_wave_rf >= waverange_rf[0]) & (s_wave_rf <= waverange_rf[1]))

            datarr_spec       = unumpy.uarray(s_flux_rf[goodent_rf], s_df_rf[goodent_rf])
            Ftot_trapz_spec   = np.trapz(datarr_spec,s_wave_rf[goodent_rf])
            Ftot_sum_spec     = np.sum(datarr_spec) * np.median(np.diff(s_wave_rf[goodent_rf]))

            #------ estimate signal to noise of emission feature ------
            Ftot_lineS2N, Ftot_lineS2N_var, Npixgood, lineS2N = uves.calc_1Dspec_S2N(s_wave,s_flux,s_df**2.0,
                                                                   waverange,verbose=False)
            Ftot_lineS2N_sigma = np.sqrt(Ftot_lineS2N_var)

            Ftot_rf, vartot_rf, Npixgood_rf, lineS2N_rf = uves.calc_1Dspec_S2N(s_wave_rf,s_flux_rf,s_df_rf**2.0,
                                                                               waverange_rf,verbose=False)

            if np.abs(Ftot_lineS2N-Ftot_rf) > 1.0:
                print(' - Ftot-Ftot_rf is larger than 1e-20cgs; stopping to enable investigation')
                pdb.set_trace()

            if np.abs(lineS2N-lineS2N_rf) > 0.1:
                print(' - linsS2N-lineS2N_ref is larger than 10%; stopping to enable investigation')
                pdb.set_trace()


            #------------ Writing to output file ------------
            outstr = str("%7.4f" % z_spec)+'  '+\
                     str("%7.4f" % zS2Nmax)+'      '+\
                     str("%7.4f" % spec_sigma_ang_obs)+'  '+\
                     str("%7.4f" % spec_sigma_ang_rf)+'  '+\
                     str("%7.4f" % temp_sigma_ang_rf)+'      '+\
                     str("%7.2f" % Fratio_spec)+'  '+\
                     str("%7.2f" % Fratio_temp)+'      '+\
                     str("%12.4f" % Ftot_spec_intr)+'  '+\
                     str("%12.4f" % Ftot_trapz_spec.nominal_value)+'  '+\
                     str("%12.4f" % Ftot_trapz_spec.std_dev)+'  '+\
                     str("%12.4f" % Ftot_sum_spec.nominal_value)+'  '+\
                     str("%12.4f" % Ftot_sum_spec.std_dev)+'  '+\
                     str("%12.4f" % Ftot_FELIS_S2Nmax)+'  '+\
                     str("%12.4f" % Ftot_FELIS_S2Nmax_err)+'  '+\
                     str("%12.4f" % FELIS_S2Nmax)+'  '+\
                     str("%12.4f" % Ngoodent)+'  '+\
                     str("%12.4f" % chi2)+'  '+\
                     str("%12.4f" % vshift_intr)+'  '+\
                     str("%12.4f" % vshift_match)+'      '+\
                     str("%12.4f" % lineS2N)+'  '+\
                     str("%12.4f" % lineS2Nwavemin)+'  '+\
                     str("%12.4f" % lineS2Nwavemax)+'  '+\
                     str("%12.4f" % lineS2N_rf)+'  '+\
                     str("%12.4f" % lineS2Nwavemin_rf)+'  '+\
                     str("%12.4f" % lineS2Nwavemax_rf)+'  '+\
                     str("%12.4f" % Ftot_lineS2N)+'  '+\
                     str("%12.4f" % Ftot_lineS2N_sigma)+'  '+\
                     specname+'  '+\
                     template+'  '
            fout.write(outstr+'\n')
    if verbose: print('\n   ...done')
    fout.close()

    fmt = 'd,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,300a,300a'
    summarydat = np.genfromtxt(summaryfile,skip_header=34,dtype=fmt,comments='#',names=True)
    return summarydat

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_mocspecFELISresults_summary(summaryfile,plotbasename,colortype='lineS2N_rf',histaxes=False,Nbins=50,
                                     overwrite=False,verbose=True):
    """
    plotting and evaluating the output from uves.gen_mocspecFELISresults_summary()

    --- INPUT ---
    summaryfile        Path and name to summary file to evaluate
    plotbasename       The based name for the plots to generate (incl. output directory)
    overwrite          Overwrite the plots if they already exist?
    verbose            Toggle verbosity

    --- EXAMPLE OF RUN ---
    import uvEmissionlineSearch as uves

    summaryfile    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_CCresults_summary.txt'
    plotbasename   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_CCresults_summary_plots/190815test_'
    uves.plot_mocspecFELISresults_summary(summaryfile,plotbasename)

    """
    if verbose: print(' - Loading and plotting the content of \n   '+summaryfile)
    fmt = 'd,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,300a,300a'
    summarydat = np.genfromtxt(summaryfile,skip_header=34,dtype=fmt,comments='#',names=True)
    specnumber = np.arange(len(summarydat))+1.0


    line         = summaryfile.split('cgs_')[-1].split('19')[0].lower()

    sigmaerrval  = {'ciii':0.1, 'civ':0.1, 'siiii':0.1, 'nv':0.1, 'mgii':0.1, 'oiii':0.1, 'heii':0.1, 'lya':0.3, 'all':0.1}
    fratioerrval = {'ciii':0.1, 'civ':0.2, 'siiii':0.1, 'nv':0.2, 'mgii':0.2, 'oiii':0.1, 'all':0.2}


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_z_lineVSfelis'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['z_spec']
    yvalues  = summarydat['z_temp_S2Nmax']
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = '$z$(mock spectrum)'
    ylabel   = '$z$(FELIS)'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[.1,600],yrange=[.1,600],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_S2N_lineVSfelis'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['lineS2N_rf']
    yvalues  = summarydat['FELIS_S2Nmax']
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = 'Rest-frame emission line S/N'
    ylabel   = 'FELIS template match S/N'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='Ftot_spec_intr',cdatvec=summarydat['Ftot_spec_intr'],
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[.1,600],yrange=[.1,600],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_Ftot_intrVStrapz'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['Ftot_spec_intr']
    xerr     = [None]*len(xvalues)
    yvalues  = summarydat['Ftot_spec_trapz']
    yerr     = summarydat['Ftot_spec_trapz_err']
    xlabel   = 'Intrinsic mock spectrum line flux [1e-20erg/s/cm$^2$]'
    ylabel   = 'Trapz int mock spectrum line flux [1e-20erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_Ftot_intrVStrapz_sigma'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['Ftot_spec_intr']
    xerr     = [None]*len(xvalues)
    yvalues  = summarydat['Ftot_spec_trapz']
    yerr     = summarydat['Ftot_spec_trapz_err']
    xlabel   = 'Intrinsic mock spectrum line flux [1e-20erg/s/cm$^2$]'
    ylabel   = 'Trapz int mock spectrum line flux [1e-20erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='sigma',cdatvec=summarydat['sigma_spec_ang_rf'],
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_Ftot_intrVSsum'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['Ftot_spec_intr']
    xerr     = [None]*len(xvalues)
    yvalues  = summarydat['Ftot_spec_sum']
    yerr     = summarydat['Ftot_spec_sum_err']
    xlabel   = 'Intrinsic mock spectrum line flux [1e-20erg/s/cm$^2$]'
    ylabel   = 'Sum*dwave mock spectrum line flux [1e-20erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_Ftot_intrVSline'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['Ftot_spec_intr']
    xerr     = [None]*len(xvalues)
    yvalues  = summarydat['Ftot_lineS2N']
    yerr     = summarydat['Ftot_lineS2N_sigma']
    xlabel   = 'Intrinsic mock spectrum line flux [1e-20erg/s/cm$^2$]'
    ylabel   = 'Trapz int +/- 3 sigma mock spectrum line flux [1e-20erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_Ftot_intrVSfelis'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['Ftot_spec_intr']
    xerr     = [None]*len(xvalues)
    yvalues  = summarydat['Ftot_FELIS_S2Nmax']
    yerr     = summarydat['Ftot_FELIS_S2Nmax_err']
    xlabel   = 'Intrinsic mock spectrum line flux [1e-20erg/s/cm$^2$]'
    ylabel   = 'FELIS line flux estimate [1e-20erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_Ftot_trapzVSfelis'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['Ftot_spec_trapz']
    xerr     = summarydat['Ftot_spec_trapz_err']
    yvalues  = summarydat['Ftot_FELIS_S2Nmax']
    yerr     = summarydat['Ftot_FELIS_S2Nmax_err']

    xlabel   = 'Trapz int mock spectrum line flux [1e-20erg/s/cm$^2$]'
    ylabel   = 'FELIS line flux estimate [1e-20erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_Ftot_trapzVSfelis_sigma'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['Ftot_spec_trapz']
    xerr     = summarydat['Ftot_spec_trapz_err']
    yvalues  = summarydat['Ftot_FELIS_S2Nmax']
    yerr     = summarydat['Ftot_FELIS_S2Nmax_err']

    xlabel   = 'Trapz int mock spectrum line flux [1e-20erg/s/cm$^2$]'
    ylabel   = 'FELIS line flux estimate [1e-20erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='sigma',cdatvec=summarydat['sigma_spec_ang_rf'],
                                                   linetype='onetoone',
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_sigma_intrVSfelis'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['sigma_spec_ang_rf']
    xerr     = [None]*len(xvalues)
    yvalues  = summarydat['sigma_temp_ang_rf']
    yerr     = [None]*len(yvalues)

    xlabel   = '$\sigma$(mock spectrum) [\AA]'
    ylabel   = '$\sigma$(FELIS) [\AA]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='onetoone',
                                                   xlog=False,ylog=False,xrange=[0.0,2.7],yrange=[0.0,2.7],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'horizontal_FtotVSFtot'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['Ftot_spec_intr']
    yvalues  = summarydat['Ftot_spec_trapz']/summarydat['Ftot_FELIS_S2Nmax']-1.0
    xerr     = [None]*len(xvalues)
    yerr     = np.sqrt( (summarydat['Ftot_spec_trapz_err']   / summarydat['Ftot_spec_trapz'])**2 +
                        (summarydat['Ftot_FELIS_S2Nmax_err'] / summarydat['Ftot_FELIS_S2Nmax'])**2 ) * np.abs(yvalues)
    xlabel   = 'F(Intrinsic) [1e-20erg/s/cm$^2$]'
    ylabel   = 'F(mock spectrum)/F(FELIS) - 1'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='horizontal',
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'horizontal_FtotVSFtot_sigma'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['Ftot_spec_intr']
    yvalues  = summarydat['Ftot_spec_trapz']/summarydat['Ftot_FELIS_S2Nmax']-1.0
    xerr     = [None]*len(xvalues)
    yerr     = np.sqrt( (summarydat['Ftot_spec_trapz_err']   / summarydat['Ftot_spec_trapz'])**2 +
                        (summarydat['Ftot_FELIS_S2Nmax_err'] / summarydat['Ftot_FELIS_S2Nmax'])**2 ) * np.abs(yvalues)
    xlabel   = 'F(Intrinsic) [1e-20erg/s/cm$^2$]'
    ylabel   = 'F(mock spectrum)/F(FELIS) - 1'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='sigma',cdatvec=summarydat['sigma_spec_ang_rf'],
                                                   linetype='horizontal',
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'horizontal_sigmaVSsigma'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['sigma_spec_ang_rf']
    yvalues  = summarydat['sigma_spec_ang_rf']/summarydat['sigma_temp_ang_rf']-1.0
    xerr     = [None]*len(xvalues)
    yerr     = np.sqrt( (summarydat['sigma_temp_ang_rf']*0.0+sigmaerrval[line] / summarydat['sigma_temp_ang_rf'])**2 +
                        (summarydat['sigma_spec_ang_rf']*0.0+0.0 / summarydat['sigma_spec_ang_rf'])**2 ) * np.abs(yvalues)

    xlabel   = '$\sigma$(mock spectrum) [\AA]'
    ylabel   = '$\sigma$(mock spectrum)/$\sigma$(FELIS) - 1'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='horizontal',
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'horizontal_zVSz'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['z_spec']
    yvalues  = summarydat['z_spec']/summarydat['z_temp_S2Nmax']-1.0
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(yvalues)

    xlabel   = '$z$(mock spectrum)'
    ylabel   = '$z$(mock spectrum)/$z$(FELIS) - 1'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='horizontal',
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'horizontal_zVSvshift'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = summarydat['z_spec']
    yvalues  = summarydat['vshift_CCmatch']
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(yvalues)

    xlabel   = '$z$(mock spectrum)'
    ylabel   = '$\Delta v$/[km/s] = $c$[$z$(mock spectrum) - $z$(FELIS)] / [1+ $z$(FELIS)]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='horizontal',
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    #-------------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------------
    goodFratio = np.where(summarydat['Fratio_spec'] > 0)
    if len(goodFratio[0]) > 0:
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'onetoone_Fratio_intrVSfelis'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = summarydat['Fratio_spec'][goodFratio]
        xerr     = [None]*len(xvalues)
        yvalues  = summarydat['Fratio_temp'][goodFratio]
        yerr     = [fratioerrval[line]]*len(yvalues)

        xlabel   = 'Doublet flux ratio mock spec'
        ylabel   = 'Doublet flux ratio FELIS match'

        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'][goodFratio],
                                                       linetype='onetoone',
                                                       xlog=False,ylog=False,xrange=[0.0,3.3],yrange=[0.0,3.3],
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'horizontal_fluxratioVSfluxratio'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = summarydat['Fratio_spec'][goodFratio]
        yvalues  = summarydat['Fratio_spec'][goodFratio]/summarydat['Fratio_temp'][goodFratio]-1.0
        xerr     = [None]*len(xvalues)
        yerr     = np.sqrt( (summarydat['Fratio_temp'][goodFratio]*0.0+fratioerrval[line] / summarydat['Fratio_temp'][goodFratio])**2 +
                            (summarydat['Fratio_spec'][goodFratio]*0.0+0.0 / summarydat['Fratio_spec'][goodFratio])**2 ) * np.abs(yvalues)

        xlabel   = 'Doublet flux ratio (FR) of mock spectrum'
        ylabel   = 'FR(mock spectrum)/FR(FELIS) - 1'

        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'][goodFratio],
                                                       linetype='horizontal',
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'horizontal_sigmaVSfluxratio'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = summarydat['sigma_spec_ang_rf'][goodFratio]
        yvalues  = summarydat['Fratio_spec'][goodFratio]/summarydat['Fratio_temp'][goodFratio]-1.0
        xerr     = [None]*len(xvalues)
        yerr     = np.sqrt( (summarydat['Fratio_temp'][goodFratio]*0.0+fratioerrval[line] / summarydat['Fratio_temp'][goodFratio])**2 +
                            (summarydat['Fratio_spec'][goodFratio]*0.0+0.0 / summarydat['Fratio_spec'][goodFratio])**2 ) * np.abs(yvalues)

        xlabel   = '$\sigma$(mock spectrum) [\AA]'
        ylabel   = 'FR(mock spectrum)/FR(FELIS) - 1'

        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'][goodFratio],
                                                       linetype='horizontal',
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'horizontal_sigmafelisVSfluxratio'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = summarydat['sigma_temp_ang_rf'][goodFratio]
        yvalues  = summarydat['Fratio_spec'][goodFratio]/summarydat['Fratio_temp'][goodFratio]-1.0
        xerr     = [None]*len(xvalues)
        yerr     = np.sqrt( (summarydat['Fratio_temp'][goodFratio]*0.0+fratioerrval[line] / summarydat['Fratio_temp'][goodFratio])**2 +
                            (summarydat['Fratio_spec'][goodFratio]*0.0+0.0 / summarydat['Fratio_spec'][goodFratio])**2 ) * np.abs(yvalues)

        xlabel   = '$\sigma$(FELIS) [\AA]'
        ylabel   = 'FR(mock spectrum)/FR(FELIS) - 1'

        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'][goodFratio],
                                                       linetype='horizontal',
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                              yrange=None,xrange=None,linetype='onetoone',ylog=False,xlog=False,
                                              colortype=None,colorcode=True,cdatvec=None,point_text=None,ids=None,
                                              overwrite=False,verbose=True,title=None,MUSEsymbolblackedge=True,
                                              photoionizationplotparam=None, showgraylimits=True,
                                              histaxes=False,Nbins=50):
    """

    """
    if summarydat == 'dummydat':
        finitevalues_ent = np.where(np.isfinite(xvalues) & np.isfinite(yvalues))[0]
        xvalues = np.asarray(xvalues)[finitevalues_ent]
        yvalues = np.asarray(yvalues)[finitevalues_ent]
        if xerr is not None:
            xerr    = np.asarray(xerr)[finitevalues_ent]
        if yerr is not None:
            yerr    = np.asarray(yerr)[finitevalues_ent]
        if point_text is not None:
            point_text = np.asarray(point_text)[finitevalues_ent]
        if ids is not None:
            ids = np.asarray(ids)[finitevalues_ent]
        if cdatvec is not None:
            cdatvec = np.asarray(cdatvec)[finitevalues_ent]

    if verbose: print(' - Setting up and generating plot')
    if os.path.isfile(plotname) & (not overwrite):
        if verbose: print('\n - WARNING: the plot '+plotname+' exists and overwrite=False so moving on \n')
    else:
        if histaxes:
            fig = plt.figure(1, figsize=(6, 6))
            # partially based on https://matplotlib.org/examples/pylab_examples/scatter_hist.html
            # definitions for the axes
            left, width = 0.20, 0.60
            bottom, height = 0.15, 0.60

            bottom_h = bottom + height + 0.01
            left_h   = left + width + 0.01

            fig.subplots_adjust(wspace=0.1, hspace=0.1,left=left, right=left+width, bottom=bottom, top=bottom+height)
            rect_histx = [left, bottom_h, width, 0.2]
            rect_histy = [left_h, bottom, 0.2, height]
        else:
            fig = plt.figure(2, figsize=(6, 5))
            fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.16, right=0.97, bottom=0.15, top=0.95)
            if photoionizationplotparam:
                fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=1.2, bottom=0.15, top=0.93)

        Fsize    = 16
        lthick   = 2
        marksize = 6
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()

        if colorcode:
            cmap    = plt.cm.get_cmap('viridis_r') # 'autumn_r'

            if cdatvec is None:
                cdatvec = summarydat[colortype]

            if colortype.lower() == 'redshift':
                # clabel  = '$z$'
                # cmin    = 0.01 # 0.0, 1.4
                # cmax    = 7.5 # 10.2, 6.2
                # cextend = 'neither'
                clabel  = '$z$'
                cmin    = 1.5
                cmax    = 6.6
                cextend = 'both'
            elif colortype.lower() == 's2nfelis':
                clabel  = 'S/N(FELIS)'
                cmin    = 3.0
                cmax    = 10.0
                cextend = 'both'
            elif colortype.lower() == 's2n':
                clabel  = 'S/N'
                cmin    = 3.0
                cmax    = 10.0
                cextend = 'both'
            elif colortype.lower() == 's2n_ciii':
                clabel  = 'CIII S/N(FELIS)'
                cmin    = 3.0
                cmax    = 10.0
                cextend = 'both'
            elif colortype.lower() == 'vshift':
                clabel  = 'Velocity shift (spec vs. template match) [km/s]'
                cmin    = 0.0
                cmax    = 1000.0
                cextend = 'both'
            elif colortype == 'f(CIII) [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]':
                clabel  = colortype
                cmin    = np.min(cdatvec[np.isfinite(cdatvec)])
                cmax    = np.max(cdatvec[np.isfinite(cdatvec)])
                cextend = 'neither'
            elif colortype.lower() == 'zmanual':
                clabel  = '$z$'
                try:
                    cmin    = np.min(cdatvec[np.isfinite(cdatvec)])
                    cmax    = np.max(cdatvec[np.isfinite(cdatvec)])
                    cextend = 'neither'
                except:
                    cmin    = 1.5
                    cmax    = 6.6
                    cextend = 'both'
            elif colortype.lower() == 'ew_0':
                clabel  = 'EW$_0$(Ly$\\alpha$) [\AA]'
                cmin    = np.min(cdatvec[np.isfinite(cdatvec)])
                cmax    = np.max(cdatvec[np.isfinite(cdatvec)])
                cextend = 'neither'
            elif colortype.lower() == 'sigma':
                clabel  = '$\sigma$ [\AA]'
                cmin    = np.min(cdatvec[np.isfinite(cdatvec)])
                cmax    = np.max(cdatvec[np.isfinite(cdatvec)])
                cextend = 'neither'
            elif colortype.lower() == 'ew0_ciii':
                clabel  = 'EW$_0$(CIII)'
                cmin    = np.min(cdatvec[np.isfinite(cdatvec)])
                cmax    = np.max(cdatvec[np.isfinite(cdatvec)])
                cextend = 'neither'
            elif colortype in summarydat.dtype.names:
                # if colortype == 'Fratio_spec': pdb.set_trace()
                clabel  = colortype.replace('_','\_')
                cmin    = np.min(cdatvec[np.isfinite(cdatvec)])
                cmax    = np.max(cdatvec[np.isfinite(cdatvec)])
                cextend = 'neither'
            else:
                sys.exit(' Color type '+colortype+' not enabled ')

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
                    colanchor   = (-1.78,0.0)
                    colshrink   = colshrink/colbarscale
                    colaspect   = colaspect/colbarscale

                cb      = plt.colorbar(m,extend=cextend,orientation='vertical',
                                       pad=0.01,aspect=colaspect,shrink=colshrink,anchor=colanchor,use_gridspec=False)
                cb.set_label(clabel)

            colvec   = []
            for ii,xval in enumerate(xvalues):
                colvec.append(cmap(colnorm(cdatvec[ii])))
            facecol  = colvec
            alphaval = 1.0
        else:
            colvec   = ['k']*len(xvalues)
            facecol  = ['gray']*len(xvalues)
            alphaval = 0.5

        #--------- RANGES ---------
        if not xrange:
            xmin   = np.min(xvalues[np.isfinite(xvalues)])
            xmax   = np.max(xvalues[np.isfinite(xvalues)])
            dx     = xmax-xmin
            xrange = [xmin-dx*0.05,xmax+dx*0.05]
        plt.xlim(xrange)
        xminsys, xmaxsys = plt.xlim() # use to get automatically expanded axes if xmin = xmax

        if not yrange:
            ymin   = np.min(yvalues[np.isfinite(yvalues)])
            ymax   = np.max(yvalues[np.isfinite(yvalues)])
            dy     = ymax-ymin
            yrange = [ymin-dy*0.05,ymax+dy*0.05]
        plt.ylim(yrange)
        yminsys, ymaxsys = plt.ylim() # use to get automatically expanded axes if xmin = xmax

        #--------- X and Y limits ---------
        if yerr is not None:
            y_uplimarr = (np.asarray(yerr) == +99)  # .astype(int) # KBS removed on 200217
            y_lolimarr = (np.asarray(yerr) == -99)  # .astype(int) # KBS removed on 200217
        else:
            y_uplimarr = [False]*len(xvalues)
            y_lolimarr = y_uplimarr

        if xerr is not None:
            x_uplimarr = (np.asarray(xerr) == +99)  # .astype(int) # KBS removed on 200217
            x_lolimarr = (np.asarray(xerr) == -99)  # .astype(int) # KBS removed on 200217
        else:
            x_uplimarr = [False]*len(xvalues)
            x_lolimarr = x_uplimarr

        for ii,xval in enumerate(xvalues): # loop necessary for coloring and upper/lower limits markers
            # change color of limits
            ecol         = colvec[ii]
            mecol        = colvec[ii]
            fcol         = facecol[ii]

            # checking for upper/lower limits
            if ids is not None:
                mfc         = True
                if (ids[ii] < 6e8): # CDFS and COSMOS
                    markersym   = 'o'
                    markerzorder = 25
                    if MUSEsymbolblackedge:
                        mecol        = 'black'
                elif (ids[ii] < 7e8) & (ids[ii] > 6e8): # UDF
                    markersym   = 'D'
                    markerzorder = 25
                    if MUSEsymbolblackedge:
                        mecol        = 'black'
                elif (ids[ii] < 9e8) & (ids[ii] > 7e8): # UDF10
                    markersym   = 'X'
                    markerzorder = 25
                    if MUSEsymbolblackedge:
                        mecol        = 'black'
                elif (ids[ii] > 1e9): # Literature objects
                    if ids[ii] == 990000000000:
                        markersym   = '.'
                        markerzorder = 20
                    else:
                        markersym   = lce.get_reference_fromID(ids[ii],verbose=False)[4]
                        mfc         = False
                        markerzorder = 22
                else:
                    print(' WARNING - stopped as could not assing a marker symbol to the id '+str(ids[ii]))
                    pdb.set_trace()
            else:
                markersym   = 'o'
                mfc         = True
            ms          = marksize
            limsizefrac = 0.05

            if yerr is not None:
                if y_uplimarr[ii].all():
                    if ylog:
                        dlog     = np.abs(np.diff(np.log10(plt.ylim()))) * limsizefrac
                        yerr[ii] = np.abs(yvalues[ii] - 10.**(np.log10(yvalues[ii])-dlog))
                    else:
                        yerr[ii] = np.abs(np.diff(plt.ylim())) * limsizefrac
                if y_lolimarr[ii].all():
                    if ylog:
                        dlog     = np.abs(np.diff(np.log10(plt.ylim()))) * limsizefrac
                        yerr[ii] = np.abs(yvalues[ii] - 10.**(np.log10(yvalues[ii])+dlog))
                    else:
                        yerr[ii] = np.abs(np.diff(plt.ylim())) * limsizefrac

            if xerr is not None:
                if x_uplimarr[ii].all():
                    if xlog:
                        dlog     = np.abs(np.diff(np.log10(plt.xlim()))) * limsizefrac
                        xerr[ii] = np.abs(xvalues[ii] - 10.**(np.log10(xvalues[ii])-dlog))
                    else:
                        xerr[ii] = np.abs(np.diff(plt.xlim())) * limsizefrac
                if x_lolimarr[ii].all():
                    if xlog:
                        dlog     = np.abs(np.diff(np.log10(plt.xlim()))) * limsizefrac
                        xerr[ii] = np.abs(xvalues[ii] - 10.**(np.log10(xvalues[ii])+dlog))
                    else:
                        xerr[ii] = np.abs(np.diff(plt.xlim())) * limsizefrac

            if (xerr is not None) & showgraylimits:
                if x_uplimarr[ii].all() or x_lolimarr[ii].all():
                    ecol         = 'darkgray'
                    mecol        = 'darkgray'
                    fcol         = 'darkgray'
                    markerzorder = 18

            if (xerr is not None):
                if len(np.atleast_1d(xerr[ii])) == 2:
                    xerrshow = [xerr[ii]]
                else:
                    xerrshow = xerr[ii]
            else:
                xerrshow = xerr

            if (yerr is not None) & showgraylimits:
                if y_uplimarr[ii].all() or y_lolimarr[ii].all():
                    ecol         = 'darkgray'
                    mecol        = 'darkgray'
                    fcol         = 'darkgray'
                    markerzorder = 18

            if (yerr is not None):
                if len(np.atleast_1d(yerr[ii])) == 2:
                    yerrshow = [yerr[ii]]
                else:
                    yerrshow = yerr[ii]
            else:
                yerrshow = yerr


            if mfc:
                markerfacecolor = fcol
            else:
                markerfacecolor = 'None'

            plt.errorbar(xvalues[ii],yvalues[ii],xerr=xerrshow,yerr=yerrshow,capthick=0.5,
                         uplims=y_uplimarr[ii],lolims=y_lolimarr[ii],xuplims=x_uplimarr[ii],xlolims=x_lolimarr[ii],
                         marker=markersym,lw=lthick/2., markersize=ms,alpha=alphaval,
                         markerfacecolor=markerfacecolor,ecolor=ecol,
                         markeredgecolor=mecol,zorder=markerzorder)

        #------------------ Drawing lines to guide the eye ------------------
        Zsun = 8.69
        xminsys, xmaxsys = plt.xlim()
        yminsys, ymaxsys = plt.ylim()
        tot_maxval = np.max([xmaxsys, ymaxsys])
        tot_minval = np.min([xminsys, yminsys])
        if linetype == 'horizontal':
            plt.plot([xminsys,xmaxsys],[0,0],'--',color='black',lw=lthick,zorder=10)
        elif linetype == 'horizontalWlya':
            plt.plot([xminsys,xmaxsys],[0,0],'--',color='black',lw=lthick,zorder=10)
            plt.plot([2.9,2.9],[yminsys,ymaxsys],'--',color='gray',lw=lthick,zorder=10)
        elif linetype == 'threesigma_y':
            plt.plot([xminsys,xmaxsys],[3.0,3.0],'--',color='black',lw=lthick,zorder=10)
        elif linetype == 'threesigma_x':
            plt.plot([3.0,3.0],[yminsys,ymaxsys],'--',color='black',lw=lthick,zorder=10)
        elif linetype == 'onetoone':
            plt.plot([-1e10,1e10],[-1e10,1e10],'--',color='black',lw=lthick,zorder=10)
        elif linetype == 'plus':
            plt.plot([-1e10,1e10],[0,0],'--',color='black',lw=lthick,zorder=10)
            plt.plot([0,0],[-1e10,1e10],'--',color='black',lw=lthick,zorder=10)
        elif linetype == 'onetooneWZsun':
            plt.plot([-1,1e5],[-1,1e5],'--',color='black',lw=lthick,zorder=10)
            plt.fill_between([5.0,Zsun,Zsun,10],[10.0,10.0,10.0,10.0],[Zsun,Zsun,5.0,5.0],color='lightgray')
            plt.text(7.0,8.75,'super-solar',color='gray',zorder=10)
        elif linetype == 'onethreeten':
            plt.plot([tot_minval,tot_maxval],[tot_minval,tot_maxval],'-',color='black',lw=lthick,zorder=10)
            plt.plot([tot_minval,tot_maxval],[tot_minval/10.,tot_maxval/10.],':',color='black',lw=lthick,zorder=10)
            plt.plot([tot_minval,tot_maxval],[tot_minval*10,tot_maxval*10],':',color='black',lw=lthick,zorder=10)
            plt.plot([tot_minval,tot_maxval],[tot_minval/3.,tot_maxval/3.],'--',color='black',lw=lthick,zorder=10)
            plt.plot([tot_minval,tot_maxval],[tot_minval*3,tot_maxval*3],'--',color='black',lw=lthick,zorder=10)
        elif linetype == 'yZsun':
            plt.fill_between([-1.0,10],[10.0,10.0],[Zsun,Zsun],color='lightgray')
            plt.text(1.0,8.75,'super-solar',color='gray',zorder=10)
        elif linetype == 'AV18_peaksep_wHorizontal':
            plt.plot([xminsys,xmaxsys],[0,0],'--',color='black',lw=lthick,zorder=10)
            Npoints = 100.0
            x         = np.arange(xminsys,xmaxsys,(xmaxsys-xminsys)/Npoints)
            yarr      = np.zeros([4,int(Npoints)])
            yarr[0,:] = (1.05-0.11) * x - 12 - 37
            yarr[1,:] = (1.05-0.11) * x - 12 + 37
            yarr[2,:] = (1.05+0.11) * x - 12 - 37
            yarr[3,:] = (1.05+0.11) * x - 12 + 37
            ylow      = np.min(yarr,axis=0)
            yhigh     = np.max(yarr,axis=0)
            plt.fill_between(x,ylow,yhigh,color='lightgray')
        elif linetype == 'AV18_fwhm_wHorizontal':
            plt.plot([xminsys,xmaxsys],[0,0],'--',color='black',lw=lthick,zorder=10)
            Npoints = 100.0
            x         = np.arange(xminsys,xmaxsys,(xmaxsys-xminsys)/Npoints)
            yarr      = np.zeros([4,int(Npoints)])
            yarr[0,:] = (0.9-0.14) * x - 36 - 60
            yarr[1,:] = (0.9-0.14) * x - 36 + 60
            yarr[2,:] = (0.9+0.14) * x - 36 - 60
            yarr[3,:] = (0.9+0.14) * x - 36 + 60
            ylow      = np.min(yarr,axis=0)
            yhigh     = np.max(yarr,axis=0)
            plt.fill_between(x,ylow,yhigh,color='lightgray')
        elif linetype == 'horizontal_and_nakajima18EWvsDv':
            plt.plot([-1e10,1e10],[0,0],'--',color='black',lw=lthick,zorder=10)

            # Nakajima+2018 Eq.1
            x = np.array([70,xmaxsys])
            y = np.array([150,150])
            plt.plot(x,y,'-.',color='gray',lw=lthick,zorder=10)

            x = np.array([20,70])
            y = 360-3.0*x
            plt.plot(x,y,'-.',color='gray',lw=lthick,zorder=10)

            x = np.array([xminsys,20])
            y = 600-15.0*x
            plt.plot(x,y,'-.',color='gray',lw=lthick,zorder=10)

            # Adelberger+2003 Eq.2
            x = np.array([0,xmaxsys])
            y = 670.0 - 8.9 * x
            plt.plot(x,y,':',color='gray',lw=lthick,zorder=10)
        elif linetype == 'CM18_wHorizontal':
            plt.plot([xminsys,xmaxsys],[0,0],'--',color='black',lw=lthick,zorder=10)

            plotTabulatedcurves = True
            if plotTabulatedcurves:
                CM18tabulated = np.genfromtxt('/Users/kschmidt/work/catalogs/Mason18_DV_Muv_z.txt',
                                              skip_header=1,comments='#',dtype=None,names=True)

                CM18_MUV      = CM18tabulated['Muv']
                zvec          = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]
                for zz, zvalue in enumerate(zvec):
                    lcolor  = cmap(colnorm(zvalue))
                    CM18_Dv = CM18tabulated['DvATz'+str(int(zvalue))]
                    plt.plot(CM18_MUV,CM18_Dv,'--',color=lcolor,lw=lthick,zorder=10)

            plotEq3curves = False
            if plotEq3curves:
                x = np.arange(xminsys,xmaxsys,0.01)
                for zvals in [[2.0,':','tab:orange'],[7.0,'-.','tab:blue']]:
                    zfct, lstyle, lcolor = zvals

                    gamma_trans = -20 - 0.26*zfct
                    ent_low     = np.where(x < gamma_trans)[0] # gamma = -0.7
                    ent_high    = np.where(x > gamma_trans)[0] # gamma = -0.3

                    gammas = [-0.7,-0.3]
                    for ee, ent in enumerate([ent_low,ent_high]):
                        gamma = gammas[ee]
                        Dv    = 10.0 ** (0.32*gamma*(x[ent] + 20.0 + 0.26*zfct) + 2.34)
                        plt.plot(x[ent],Dv,lstyle,color=lcolor,lw=lthick,zorder=10)

        elif (str(linetype).lower == 'none') or (linetype is None) :
            pass
        else:
            sys.exit(' Unknown value of linetype = "'+linetype+'"')

        #------------------ PHOTOIONIZATION GRIDS ------------------
        titleaddition = ''
        if photoionizationplotparam is not None:
            # titleaddition = uves.add_photoionization_models_to_lineratioplot(photoionizationplotparam)
            titleaddition = lce.add_photoionization_models_to_plot(photoionizationplotparam)

        if (title is not None) & (histaxes == False):
            plt.title(title+titleaddition,fontsize=Fsize-4)

        for ii,xval in enumerate(xvalues): # loop necessary for coloring
            if point_text is not None:
                if xlog:
                    xtextval = xvalues[ii]*1.03
                else:
                    xminsys, xmaxsys = plt.xlim()
                    xtextval = xvalues[ii] + (xmaxsys-xminsys)*0.02

                plt.text(xtextval,yvalues[ii],
                         point_text[ii],color='white',fontsize=Fsize*0.2,zorder=30,
                         bbox=dict(boxstyle="round",edgecolor='k',facecolor=colvec[ii],linewidth=lthick*0.2))

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if ylog:
            plt.yscale('log')
        if xlog:
            plt.xscale('log')
        #------------------ HISTOGRAM AXES ------------------
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

            if ids is not None:
                axHistx.hist(xvalues[ids.astype(int) < 1e9][np.isfinite(xvalues[ids.astype(int) < 1e9])],linestyle='-',
                             bins=bindefs,histtype='step',color='k')
                if len(np.where(ids.astype(int) >= 1e9)[0]) > 0: # checking for literature IDs
                    axHistx.hist(xvalues[np.isfinite(xvalues)],linestyle=':',
                                 bins=bindefs,histtype='step',color='k')
            else:
                axHistx.hist(xvalues[np.isfinite(xvalues)], bins=bindefs,histtype='step',color='k')
            axHistx.set_xticks([])
            axHistx.set_xlim([xminsys,xmaxsys])

            binwidth_y = np.diff([yminsys,ymaxsys])/Nbins
            bindefs    = np.arange(yminsys, ymaxsys+binwidth_y, binwidth_y)
            if ylog:
                bindefs = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))
                axHisty.set_yscale('log')

            if ids is not None:
                axHisty.hist(yvalues[ids.astype(int) < 1e9][np.isfinite(yvalues[ids.astype(int) < 1e9])],linestyle='-',
                             bins=bindefs,histtype='step',color='k', orientation='horizontal')
                if len(np.where(ids.astype(int) >= 1e9)[0]) > 0: # checking for literature IDs
                    axHisty.hist(yvalues[np.isfinite(yvalues)],linestyle=':',
                                 bins=bindefs,histtype='step',color='k', orientation='horizontal')
            else:
                axHisty.hist(yvalues[np.isfinite(yvalues)], bins=bindefs,histtype='step',color='k', orientation='horizontal')
            axHisty.set_yticks([])
            axHisty.set_ylim([yminsys,ymaxsys])

            if colorcode:
                cb      = plt.colorbar(m,extend=cextend,orientation='vertical',
                                       pad=0.01,aspect=10,shrink=0.35,anchor=(-15.0,1.58),use_gridspec=False)
                cb.set_label(clabel)
        else:
            pass

        #--------- LEGEND ---------
        # plt.errorbar(-5000,-5000,xerr=None,yerr=1,marker='o',lw=0, markersize=marksize,alpha=1.0,
        #              markerfacecolor='k',ecolor='k',markeredgecolor='black',zorder=1,label='MUSE-Wide LAE')
        # plt.errorbar(-5000,-5000,xerr=None,yerr=None,marker='*',lw=0, markersize=marksize*2,alpha=1.0,
        #              markerfacecolor='None',ecolor='None',markeredgecolor='black',zorder=1,label='AGN')
        # plt.errorbar(-5000,-5000,xerr=None,yerr=None,marker='D',lw=0, markersize=marksize,alpha=1.0,
        #              markerfacecolor='None',ecolor='None',markeredgecolor='black',zorder=1,label='AGN candidate')
        #
        # leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize/1.0},ncol=5,numpoints=1,
        #                  bbox_to_anchor=(0.5, 1.1),)  # add the legend
        # leg.get_frame().set_alpha(0.7)
        #--------------------------

        if verbose: print('   Saving plot to '+plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def add_photoionization_models_to_lineratioplot(piplotparam,verbose=True):
    """
    Wrapper to add NEOGAL photoionization model grids to flux ratio plots generated with
    uves.plot_mocspecFELISresults_summary_plotcmds() - see e.g. uves.plot_lineratios_fromsummaryfiles_wrapper()

    plotting based on rxj2248_BooneBalestraSource.plot_feltregutkinmodels()
    which is based on NEOGALmodels.plot_lineratios()

    NB! as of 200622 use literaturecollection_emissionlinestrengths.add_photoionization_models_to_plot() instead.


    --- INPUT ---
    piplotparam            The photoionization plot parameters.

    --- EXAMPLE OF RUN ---

    """
    modeldataSF  = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/nebular_emission/')
    modeldataAGN = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/AGN_NLR_nebular_feltre16/')

    linestrings, doubletratios, varyparam, cutSFmodels, markersize, SFmarker, AGNmarker = piplotparam
    # varyparam options: 'Zgas','logUs','xid','nh','COratio','Mcutoff'
    logcolors = ['Zgas']

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
    NgoodentSF  = len(modeldataSF)

    if NgoodentSF > 1:
        if verbose: print(' - Getting data for '+str(NgoodentSF)+' data points satisfying (SFR)model selection ')
        varydatSF  = modeldataSF[varyparam]
        if varyparam in logcolors:
            varydatSF  = np.log10(modeldataSF[varyparam])

        linedataSF = []
        for ll, linestr in enumerate(linestrings):
            if linestr == 'CIV1550':
                linedataSF.append(modeldataSF['CIV1548']+modeldataSF['CIV1551'])
            elif linestr == 'CIII1907':
                linedataSF.append(modeldataSF['CIII1908']/(1.+1./doubletratios[ll]))
            elif linestr == 'CIII1910':
                ciii1907 = modeldataSF['CIII1908']/(1.+1./doubletratios[ll])
                linedataSF.append(ciii1907/doubletratios[ll])
            elif linestr == 'OIII1663':
                linedataSF.append(modeldataSF['OIII1661']+modeldataSF['OIII1666'])
            elif linestr == 'SiIII1892':
                linedataSF.append(modeldataSF['SiIII1888']-modeldataSF['SiIII1883'])
            else:
                linedataSF.append(modeldataSF[linestr])

        # --------- Line Ratios ---------
        ratioSF_x   = linedataSF[0]/linedataSF[1]
        ratioSF_y   = linedataSF[2]/linedataSF[3]

    else:
        print(' WARNING: uves.add_photoionization_models_to_lineratioplot >>>'
              ' Less than 2 (SFR)model grid points to plot; no data added to plot')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    NgoodentAGN  = len(modeldataAGN)

    if NgoodentAGN > 1:
        if verbose: print(' - Getting data for '+str(NgoodentAGN)+' data points satisfying (AGN)model selection ')
        varydatAGN  = modeldataAGN[varyparam]
        if varyparam in logcolors:
            varydatAGN = np.log10(modeldataAGN[varyparam])

        linedataAGN = []
        for ll, linestr in enumerate(linestrings):
            if linestr == 'CIV1550':
                linedataAGN.append(modeldataAGN['CIV1548']+modeldataAGN['CIV1551'])
            elif linestr == 'OIII1663':
                linedataAGN.append(modeldataAGN['OIII1661']+modeldataAGN['OIII1666'])
            elif linestr == 'CIII1908':
                linedataAGN.append(modeldataAGN['CIII1907']+modeldataAGN['CIII1910'])
            elif linestr == 'SiIII1892':
                linedataAGN.append(modeldataAGN['SiIII1888']-modeldataAGN['SiIII1883'])
            else:
                linedataAGN.append(modeldataAGN[linestr])

        # --------- Line Ratios ---------
        ratioAGN_x   = linedataAGN[0]/linedataAGN[1]
        ratioAGN_y   = linedataAGN[2]/linedataAGN[3]
    else:
        print(' WARNING: uves.add_photoionization_models_to_lineratioplot >>>'
              ' Less than 2 (AGN)model grid points to plot; no data added to plot')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    edgecol  = 'None'

    cmap    = plt.cm.get_cmap('Reds') #'plasma', 'Greys' 'copper_r'
    cmin    = np.min(np.append(varydatSF,varydatAGN))
    cmax    = np.max(np.append(varydatSF,varydatAGN))
    colnorm = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, num=50) #cmax-cmin)
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
    colanchor   = (0.0,1.0)
    colshrink   = colshrink/colbarscale
    colaspect   = colaspect/colbarscale
    cextend     = 'neither'

    cb1      = plt.colorbar(mm,extend=cextend,orientation='vertical',ticks=colortickvals,
                            pad=0.01,aspect=colaspect,shrink=colshrink,anchor=colanchor,use_gridspec=False)
    cb1.ax.set_yticklabels(colorlabels,rotation=0)
    cb1.set_label(cbarlegend)     #Zgas is unitless as it is a mass ratio (see Gutkin et al. 2016 Eq. 10)

    for vdSF in np.unique(varydatSF):
        SFcol    = cmap(colnorm(vdSF))
        SFcolent = np.where(varydatSF == vdSF)

        plt.scatter(ratioSF_x[SFcolent],ratioSF_y[SFcolent],s=markersize,
                    marker=SFmarker,lw=0.2, facecolor='None',edgecolor=SFcol, zorder=5)

    for vdAGN in np.unique(varydatAGN):
        AGNcol    = cmapAGN(colnorm(vdAGN)) # 'gray'
        AGNcolent = np.where(varydatAGN == vdAGN)

        plt.scatter(ratioAGN_x[AGNcolent],ratioAGN_y[AGNcolent],s=markersize,
                    marker=AGNmarker,lw=0.2, facecolor='None',edgecolor=AGNcol, zorder=5)

    titleaddition = infostrSFcut
    return titleaddition
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def match_tdosespectra_to_templates(spectra,outputdir,outkeystr='FELISmatch_RENAME',
                                    subtract_spec_median=True,
                                    lines2find=['CIII','CIV','HeII','NV','OIII','SiIII','MgII'],
                                    CCwavewindow=10.0,plot_allCCresults=False,
                                    overwrite=False,verbose_FELIS=False,verbose=True):
    """
    Wrapper around felis.match_templates2specs() to match FELIS templates to a sample of TDOSE spectra.

    --- INPUT ---
    spectra                List of spectra to search for emission features via the FELIS template matches.
    outputdir              Directory to contiain the output (plots and *pkl files).
    outkeystr              String used in naming output.
    subtract_spec_median   Subtract the median during the template match to approximate the continuum in the
                           matched region of the spectrum.
    lines2find             List of lines to look for; determines what templates to match to each spectrum.
    CCwavewindow           Rest-frame width (central wavelengths +/- CCwavewindow) to match templates in.
    plot_allCCresults      Plot all template matches; timeconsuming so default is False. Can always be done afterwards.
    overwrite              Overwrite the output pickle file if it already exists?
    verbose_FELIS          Set to true to get vebosity of FELIS matching
    verbose                Set to true to get vebpsity of progress of the matching by the wrapper.

    --- EXAMPLE OF RUN ---
    import uvEmissionlineSearch as uves

    specdir   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/spectra_aperture/'
    spectra   = glob.glob(specdir+'tdose_spectrum*aperture_07209*.fits')
    outputdir = ' /Users/kschmidt/work/MUSE/uvEmissionlineSearch/FELIStemplatematch2uvesobjects/'
    uves.match_tdosespectra_to_templates(spectra,outputdir,outkeystr='FELISmatch2udf10_07209starUVESobj190913')

    """
    templatedir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_templates_fromsetup/'
    uvesobjinfo = afits.open('/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo_UVemitters_3timesUDFcats.fits')[1].data
    MUSEwaves   = np.array([4790,9310])
    waverest    = {'lya':1215.67, 'ciii':1907.7, 'civ':1549.5, 'nv':1240.8,
                   'heii':1640.4, 'oiii':1663.5, 'mgii':2799.1, 'siiii':1887.4}
    matchno     = 0
    if verbose: print('\n - Will match the '+str(len(spectra))+' spectra to templates of the lines '+str(lines2find))
    for ss, spectrum in enumerate(spectra):
        obj_id       = int(spectrum.split('_')[-1].split('.fit')[0])
        obj_ent      = np.where(uvesobjinfo['id'] == obj_id)[0]
        if len(obj_ent) != 1:
            sys.exit('ERROR - There are '+str(len(obj_ent))+' matches in the UVES info file to object '+str(obj_id))

        obj_z         = uvesobjinfo['redshift'][obj_ent]
        MUSEwaves_rf  = MUSEwaves / (1.0 + obj_z)

        for templine in lines2find:
            Nmatch      = len(spectra)*len(lines2find)
            matchno     = matchno+1
            if verbose:
                infostr = '   FELIS match '+str("%.5d" % (matchno))+' / '+str("%.5d" % Nmatch)+'     '
                sys.stdout.write("%s\r" % infostr)
                sys.stdout.flush()

            linewave_rf  = waverest[templine.lower()]
            if (MUSEwaves_rf[0]+CCwavewindow < linewave_rf) & (linewave_rf < MUSEwaves_rf[1]-CCwavewindow):
                templates  = glob.glob(templatedir+'uves_felis_template_fromsetup_'+templine+'*.fits')
                if len(templates) == 0:
                    print('\n WARNING: No templates found for "lines2find"='+str(templine))
                else:
                    picklefile = outputdir+spectrum.split('/')[-1].replace('.fits','_CCresults_template'+
                                                                           templine+'_'+outkeystr+'.pkl')

                    wavewindow = CCwavewindow * (1+obj_z) # turn CCwavewindow into observed frame
                    ccdic      = felis.match_templates2specs(templates,[spectrum],[obj_z],
                                                             picklefile,wavewindow=[wavewindow],
                                                             plotdir=outputdir,
                                                             wavecen_restframe=[linewave_rf],
                                                             vshift=None,min_template_level=1e-4,
                                                             plot_allCCresults=plot_allCCresults,
                                                             subtract_spec_median=subtract_spec_median,
                                                             overwrite=overwrite,
                                                             verbose=verbose_FELIS)
            else:
                pass
    if verbose: print('\n   done...\n')


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_tdosespecFELISresults_summary(summaryfile,picklefiles,overwrite=False,verbose=True):
    """
    Generate a summary of the template characteristics FELIS determined to match the tdose spectra
    the best, i.e. with the highest S/N.

    --- INPUT ---
    summaryfile        Path and name of summary file to generate.
    picklefiles        List of FELIS pickle files to summarize.
    overwrite          Overwrite the summary file if it already exists?
    verbose            Toggle verbosity

    --- EXAMPLE OF RUN ---
    import uvEmissionlineSearch as uves

    FELISoutputdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_CCresults/'
    summaryfile    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra_CCresults_summary_all.txt'
    picklefiles    = glob.glob(FELISoutputdir+'*.pkl')
    summarydat     = uves.gen_mocspecFELISresults_summary(summaryfile,picklefiles)

    """
    if verbose: print(' - Generating a summary of the best-fit template in:\n   '+summaryfile)
    if os.path.isfile(summaryfile) & (not overwrite):
        sys.exit(' Summary file '+summaryfile+' already exists and overwrite=False ')

    fout = open(summaryfile,'w')
    fout.write('# Summary of '+str(len(picklefiles))+' FELIS pickle files provided \n')
    fout.write('# File contains the characteristics of the templates with max S/N from the FELIS template fits \n')
    fout.write('# The summary was generated with uves.gen_tdosespecFELISresults_summary() on '+kbs.DandTstr2()+' \n')
    fout.write('# \n')
    fout.write('# Columns are:\n')
    fout.write('# id                         ID of object matched \n')
    fout.write('# z_spec                     Instrinspic redshift of matched mock spectrum \n')
    fout.write('# z_temp_S2Nmax              Estimated redshift from template match\n')
    fout.write('# sigma_temp_ang_rf          Rest-frame line width in angstroms for maxS/N template \n')
    fout.write('# Fratio_temp                Flux ratio of doublets in template      (line_lowwave/line_highwave) \n')
    fout.write('# Ftot_FELIS_S2Nmax          Estimated total flux from FELIS template match for maxS/N match \n')
    fout.write('# Ftot_FELIS_S2Nmax_err      Uncertainty on FELISflux_S2Nmax [sqrt(Ftot_FELIS_S2Nmax_variance)]\n')
    fout.write('# FELIS_S2Nmax               The S/N value of the (scaled) template match to the mock spectrum \n')
    fout.write('# Ngoodent                   The number of good pixels used in the cross correlation \n')
    fout.write('# chi2                       Chi^2 value between the mock spectrum and the template match \n')
    fout.write('# vshift_spec                Known intrinsic velocity shift of mock spectrum \n')
    fout.write('# vshift_CCmatch             Estimated velocity shift from template match '
               ' [ c*(z_spec-z_temp_S2Nmax)/(1+z_temp_S2Nmax) ]\n')
    fout.write('# lineS2N_rf                 Estimated S/N (rest-frame) of spectral feature within [lineS2Nwavemin_rf,lineS2Nwavemin_rf] \n')
    fout.write('# lineS2Nwavemin_rf          Lower integration limit for rest-frame S/N estimate \n')
    fout.write('# lineS2Nwavemax_rf          Upper integration limit for rest-frame S/N estimate \n')
    fout.write('# Ftot_lineS2N_rf            Total flux of spectral feature (sum(f)*dwave) used to estimate line S/N \n')
    fout.write('# Ftot_lineS2N_sigma_rf      Square root of the variance/sqrt(Npix) of Ftot \n')
    fout.write('# spectrum                   The mock spectrum the templates were matched to \n')
    fout.write('# template                   The maxS/N template \n')
    fout.write('# \n')
    fout.write('# id z_spec z_temp_S2Nmax sigma_temp_ang_rf Fratio_temp Ftot_FELIS_S2Nmax Ftot_FELIS_S2Nmax_err FELIS_S2Nmax Ngoodent chi2 vshift_spec vshift_CCmatch lineS2N_rf lineS2Nwavemin_rf lineS2Nwavemax_rf Ftot_lineS2N_rf Ftot_lineS2N_sigma_rf spectrum template \n')

    for pp, picklefile in enumerate(picklefiles):
        if verbose:
            infostr = ' - Summarizing picklefile  '+str("%.5d" % (pp+1))+' / '+str("%.5d" % len(picklefiles))+'     '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()
        pkldic  = felis.load_picklefile(picklefile)

        Nsigma_integration = 3.0

        for specname in pkldic.keys():
            tempdic = pkldic[specname]

            #------ load info about max S/N template ------
            template, vshift_intr, vshift_match, Ftot_FELIS_S2Nmax, Ftot_FELIS_S2Nmax_err, \
            FELIS_S2Nmax, Ngoodent, chi2, z_spec, zS2Nmax =  \
                felis.getresult4maxS2N(pkldic,specname)

            #------ load matched spec and move to restframe ------
            s_wave   , s_flux   , s_df   , s_s2n    = felis.load_spectrum(specname,verbose=False)
            s_wave_rf, s_flux_rf, s_df_rf, s_s2n_rf = s_wave / (1+z_spec), s_flux * (1+z_spec), s_df * (1+z_spec), s_s2n

            #------ extract info on mock spectrum from fits headers ------
            spec_hdr        = afits.open(specname)[1].header
            spec_sigma_ang  = np.array([])
            spec_flux       = np.array([])
            spec_line_wave  = np.array([])
            for hdrkey in spec_hdr.keys():
                if ('noise' not in hdrkey.lower()) & ('err' not in hdrkey.lower()):
                    if '_1' in hdrkey: spec_line_wave = np.append(spec_line_wave,spec_hdr[hdrkey])
                    if '_2' in hdrkey: spec_sigma_ang = np.append(spec_sigma_ang,spec_hdr[hdrkey])
                    if '_4' in hdrkey: spec_flux      = np.append(spec_flux,spec_hdr[hdrkey])

            spec_sigma_ang_obs = np.mean(spec_sigma_ang)
            spec_sigma_ang_rf  = spec_sigma_ang_obs / (1.0+z_spec)
            Ftot_spec_intr     = np.sum(spec_flux)
            if len(spec_flux) == 2:
                Fratio_spec = spec_flux[np.where(spec_line_wave == np.min(spec_line_wave))] / \
                              spec_flux[np.where(spec_line_wave == np.max(spec_line_wave))]
            else:
                Fratio_spec = 0.0

            #------ extract info on template from fits headers ------
            temp_hdr        = afits.open(template)[1].header
            temp_sigma_ang  = np.array([])
            temp_flux       = np.array([])
            temp_line_wave  = np.array([])
            for hdrkey in temp_hdr.keys():
                if ('noise' not in hdrkey.lower()) & ('err' not in hdrkey.lower()):
                    if '_1' in hdrkey: temp_line_wave = np.append(temp_line_wave,temp_hdr[hdrkey])
                    if '_2' in hdrkey: temp_sigma_ang = np.append(temp_sigma_ang,temp_hdr[hdrkey])
                    if '_4' in hdrkey: temp_flux      = np.append(temp_flux,temp_hdr[hdrkey])

            temp_sigma_ang_rf  = np.mean(temp_sigma_ang)
            if len(temp_flux) == 2:
                Fratio_temp = temp_flux[np.where(temp_line_wave == np.min(temp_line_wave))] / \
                              temp_flux[np.where(temp_line_wave == np.max(temp_line_wave))]
            else:
                Fratio_temp = 0.0

            #------ estimate signal to noise by integrating observed spectrum in restfram ------
            lineS2Nwavemin_rf = np.min(temp_line_wave)-Nsigma_integration*temp_sigma_ang_rf
            lineS2Nwavemax_rf = np.max(temp_line_wave)+Nsigma_integration*temp_sigma_ang_rf
            waverange_rf      = [lineS2Nwavemin_rf,lineS2Nwavemax_rf]
            Ftot_lineS2N_rf, Ftot_lineS2N_var_rf, Npixgood_rf, lineS2N_rf = \
                uves.calc_1Dspec_S2N(s_wave_rf,s_flux_rf,s_df_rf**2.0,waverange_rf,verbose=False)
            Ftot_lineS2N_sigma_rf = np.sqrt(Ftot_lineS2N_var_rf)

            #------------ Writing to output file ------------
            outstr = specname.split('_')[-1].split('.fit')[0]+'  '+\
                     str("%7.8f" % z_spec)+'  '+\
                     str("%7.8f" % zS2Nmax)+'      '+\
                     str("%7.4f" % temp_sigma_ang_rf)+'      '+\
                     str("%7.2f" % Fratio_temp)+'      '+\
                     str("%12.4f" % Ftot_FELIS_S2Nmax)+'  '+\
                     str("%12.4f" % Ftot_FELIS_S2Nmax_err)+'  '+\
                     str("%12.4f" % FELIS_S2Nmax)+'  '+\
                     str("%12.4f" % Ngoodent)+'  '+\
                     str("%12.4f" % chi2)+'  '+\
                     str("%12.4f" % vshift_intr)+'  '+\
                     str("%12.4f" % vshift_match)+'      '+\
                     str("%12.4f" % lineS2N_rf)+'  '+\
                     str("%12.4f" % lineS2Nwavemin_rf)+'  '+\
                     str("%12.4f" % lineS2Nwavemax_rf)+'  '+\
                     str("%12.4f" % Ftot_lineS2N_rf)+'  '+\
                     str("%12.4f" % Ftot_lineS2N_sigma_rf)+'  '+\
                     specname+'  '+\
                     template+'  '
            fout.write(outstr+'\n')
    if verbose: print('\n   ...done')
    fout.close()

    fmt = '12a,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,300a,300a'
    summarydat = np.genfromtxt(summaryfile,skip_header=25,dtype=fmt,comments='#',names=True)
    return summarydat


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_tdosespecFELISresults_summary(summaryfile,plotbasename,colortype='lineS2N_rf',obj2show='all',
                                       histaxes=True,Nbins=50,S2Ncut=[0.0,1000.0],point_text=None,
                                       overwrite=False,verbose=True):
    """
    plotting and evaluating the output from uves.gen_tdosespecFELISresults_summary()

    --- INPUT ---
    summaryfile        Path and name to summary file to evaluate
    plotbasename       The based name for the plots to generate (incl. output directory)
    colortype          The type of color bar to show
    obj2show           To only show a sumbsample of objects from the summary file provide the ids here.
    overwrite          Overwrite the plots if they already exist?
    verbose            Toggle verbosity

    --- EXAMPLE OF RUN ---
    import uvEmissionlineSearch as uves

    outdir         = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FELIStemplatematch2uvesobjects/CCresults_summary/'
    summaryfile    = outdir+'CCresults_summary_templateCIII_FELISmatch2udf10_07209starUVESobj190913.txt'
    plotbasename   = outdir+'plottest_'
    uves.plot_mocspecFELISresults_summary(summaryfile,plotbasename)

    """
    if verbose: print(' - Loading and plotting the content of \n   '+summaryfile)
    fmt = '12a,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,300a,300a'
    summarydat  = np.genfromtxt(summaryfile,skip_header=25,dtype=fmt,comments='#',names=True)
    Nspecin     = len(summarydat['spectrum'])

    if obj2show is not 'all':
        showent = np.array([])
        for objid in obj2show:
            objent = np.where( summarydat['id'].astype(int) == objid)[0]
            if len(objent) > 0:
                showent = np.append(showent,objent)
    else:
        showent = np.arange(len(summarydat))

    if verbose: print(' - Plotting FELIS matches in summary file\n   '+summaryfile+'\n   where the following holds:')
    if verbose: print('    S/N(FELIS)          = ['+str(S2Ncut[0])+','+str(S2Ncut[1])+']    (both ends included)')
    selectionAll  = np.where( (summarydat['FELIS_S2Nmax'] >= S2Ncut[0]) & (summarydat['FELIS_S2Nmax'] <= S2Ncut[1]) )[0]

    selection     = []
    for selent in selectionAll:
        if selent in showent:
            selection.append(selent)

    Nselspec    = len(selection)
    if Nselspec > 0:
        selecteddat = summarydat[selection]
        if verbose: print(' - '+str(Nselspec)+'/'+str(Nspecin)+' matched spectra in summary satisfies the cuts\n')
    else:
        if verbose: print(' WARNING No FELIS matches found in summary file satisfying cuts; returning...')
        return

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_z_lineVSfelis'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = selecteddat['z_spec']
    yvalues  = selecteddat['z_temp_S2Nmax']
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = '$z$(spectrum)'
    ylabel   = '$z$(FELIS)'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,selecteddat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=selecteddat['FELIS_S2Nmax'],
                                                   linetype='onetoone',point_text=point_text,
                                                   xlog=False,ylog=False,xrange=[1,6.5],yrange=[1,6.5],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_S2N_lineVSfelis'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = selecteddat['lineS2N_rf']
    yvalues  = selecteddat['FELIS_S2Nmax']
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(xvalues)
    xlabel   = 'Integrated spectrum flux S/N (+/-3$\sigma$)'
    ylabel   = 'S/N(FELIS)'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,selecteddat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='sigma',cdatvec=selecteddat['sigma_temp_ang_rf'],
                                                   linetype='onetoone',point_text=point_text,
                                                   xlog=True,ylog=True,xrange=[1.,200],yrange=[1.,200],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_Ftot_lineVSfelis'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = selecteddat['Ftot_lineS2N_rf']
    xerr     = selecteddat['Ftot_lineS2N_sigma_rf']
    yvalues  = selecteddat['Ftot_FELIS_S2Nmax']
    yerr     = selecteddat['Ftot_FELIS_S2Nmax_err']

    xlabel   = 'Integrated spectrum flux (+/-3$\sigma$) [1e-20erg/s/cm$^2$]'
    ylabel   = 'F(FELIS) [1e-20erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,selecteddat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=selecteddat['FELIS_S2Nmax'],
                                                   linetype='onetoone',point_text=point_text,
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_Ftot_lineVSfelis_sigma'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = selecteddat['Ftot_lineS2N_rf']
    xerr     = selecteddat['Ftot_lineS2N_sigma_rf']
    yvalues  = selecteddat['Ftot_FELIS_S2Nmax']
    yerr     = selecteddat['Ftot_FELIS_S2Nmax_err']

    xlabel   = 'Integrated spectrum flux (+/-3$\sigma$) [1e-20erg/s/cm$^2$]'
    ylabel   = 'F(FELIS) [1e-20erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,selecteddat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='sigma',cdatvec=selecteddat['sigma_temp_ang_rf'],
                                                   linetype='onetoone',point_text=point_text,
                                                   xlog=True,ylog=True,xrange=[10,2e4],yrange=[10,2e4],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)




    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_sigmaVSfelisflux'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = selecteddat['sigma_temp_ang_rf']
    xerr     = [None]*len(xvalues)
    yvalues  = selecteddat['Ftot_FELIS_S2Nmax']
    yerr     = selecteddat['Ftot_FELIS_S2Nmax_err']

    xlabel   = '$\sigma$(FELIS) [\AA]'
    ylabel   = 'F(FELIS) [1e-20erg/s/cm$^2$]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,selecteddat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=selecteddat['FELIS_S2Nmax'],
                                                   linetype='onetoone',point_text=point_text,
                                                   xlog=False,ylog=True,xrange=[0.0,2.7],yrange=[10.0,2e4],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'horizontal_zVSvshift'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = selecteddat['z_spec']
    yvalues  = selecteddat['vshift_CCmatch']
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(yvalues)

    xlabel   = '$z$(spectrum)'
    ylabel   = '$\Delta v$/[km/s] = $c$[$z$(spectrum)-$z$(FELIS)] / [1+$z$(FELIS)]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='horizontal',point_text=point_text,
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext  = 'onetoone_vshift_literatureVSfelis'
    plotname = plotbasename+nameext+'.pdf'
    xvalues  = selecteddat['vshift_spec']
    yvalues  = selecteddat['vshift_CCmatch']
    xerr     = [None]*len(xvalues)
    yerr     = [None]*len(yvalues)

    xlabel   = '$\Delta v$(Literature)/[km/s]'
    ylabel   = '$\Delta v$/[km/s] = $c$[$z$(spectrum)-$z$(FELIS)] / [1+$z$(FELIS)]'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,summarydat,
                                                   histaxes=histaxes,Nbins=Nbins,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'],
                                                   linetype='horizontal',point_text=point_text,
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    #-------------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------------
    goodFratio = np.where(selecteddat['Fratio_temp'] > 0)
    if len(goodFratio[0]) > 0:
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nameext  = 'onetoone_FratioVSs2n'
        plotname = plotbasename+nameext+'.pdf'
        xvalues  = selecteddat['Fratio_temp'][goodFratio]
        xerr     = [None]*len(xvalues)
        yvalues  = selecteddat['Ftot_FELIS_S2Nmax'][goodFratio]
        yerr     = selecteddat['Ftot_FELIS_S2Nmax_err'][goodFratio]

        xlabel   = 'FR(FELIS)'
        ylabel   = 'F(FELIS) [1e-20erg/s/cm$^2$]'

        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,selecteddat,
                                                       histaxes=histaxes,Nbins=Nbins,
                                                       colortype='s2nfelis',cdatvec=selecteddat['FELIS_S2Nmax'][goodFratio],
                                                       linetype='onetoone',point_text=point_text,
                                                       xlog=False,ylog=True,xrange=[0.0,3.3],yrange=[10.0,2e4],
                                                       colorcode=True,overwrite=overwrite,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def calc_1Dspec_S2N(wavelengths,fluxes,variances,waverange,verbose=True):
    """
    Estimating the signal to noise ratio of a defined region in a 1D spectrum.
    Signal and noise is obtained by trapexoidal integration of the flux propogating errors.

    --- INPUT ---
    wavelengths      Wavelength vector of 1D spec
    fluxes           Flux values of pixels in 1D spec
    variance         Variances for fluxes in 1D spec
    waverange        Wavelenght range to estimate S/N over
    verbose          Toggle verbosity

    --- EXAMPLE OF USE ---
    import astropy.io.fits as afits
    import uvEmissionlineSearch as uves
    import glob

    specs = glob.glob('/Users/kschmidt/work/MUSE/uvEmissionlineSearch/mockspectra/uves_mock_spectrum_fromsetup_CIIIdoublet_noisespec_sigma0p50_skew0p00_Ftot85p71_Fratio1p40_z*fits')

    for spec in specs:
        redshift    = float(spec.split('io1p40_z')[-1].split('.fit')[0].replace('p','.'))
        wavelengths = afits.open(spec)[1].data['wave']
        fluxes      = afits.open(spec)[1].data['flux']
        variances   = afits.open(spec)[1].data['fluxerror']**2.0
        dwave_rest  = 5.0
        waverange   = [(1908.0-dwave_rest/2.)*(1.0+redshift),(1908.0+dwave_rest/2.)*(1.0+redshift)]

        print(' --- Calc for '+spec.split('/')[-1]+';\n --- waverange='+str(waverange))
        Ftot, vartot, Npix, S2N = uves.calc_1Dspec_S2N(wavelengths,fluxes,variances,waverange)

    """
    if verbose: print(' - Estimating S/N of 1D spectral range '+str(waverange))
    goodent = np.where((wavelengths >= waverange[0]) & (wavelengths <= waverange[1]) & np.isfinite(fluxes) & (variances != 0))

    if len(goodent[0]) == 0.0:
        if verbose: print(' - No good (finite) pixels in wavelength range '+str(waverange))
        Ftot, vartot, Npix, S2N = 0.0, 0.0, 0.0, 0.0
    else:
        Npix   = len(goodent[0])
        if Npix == 1:
            dwave  = np.median(np.diff(wavelengths))
            Ftot   = fluxes[goodent] * dwave
            vartot = variances[goodent] * dwave**2
            S2N    = Ftot/np.sqrt(vartot)
        else:
            datarr = unumpy.uarray(fluxes[goodent], np.sqrt(variances[goodent]))
            Ftot   = np.trapz(datarr,wavelengths[goodent])
            S2N    = Ftot.nominal_value/Ftot.std_dev
            Ftot, vartot = Ftot.nominal_value, Ftot.std_dev**2

    if verbose: print(' - Returning values  Ftot(trapz), vartot, Npix, S/N = '+
                      str(Ftot)+', '+str(vartot)+', '+str(Npix)+', '+str(S2N)+'')
    return Ftot, vartot, Npix, S2N

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def match_MUSEWideLAEs(templatedir,zrange=[1.516,3.874],datestr='dateofrun',line='CIII',
                       wave_restframe=1908.0,generateplots=False,specificobj=None,
                       lamwidth_restframe='dvoffset',runonallspecs=False,subtract_spec_median=True,verbose=True):
    """
    Wrapper around felis.match_templates2specs()

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves

    #list of IDs of test objects
    specificobj = [214002011,123048186,115003085]

    tempdir  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_testing/'

    picklefileCIII = uves.match_MUSEWideLAEs(tempdir,zrange=[1.516,3.874],line='CIII',wave_restframe=1908.0,generateplots=False)
    picklefileCIV  = uves.match_MUSEWideLAEs(tempdir,zrange=[2.100,4.996],line='CIV',wave_restframe=1549.0,generateplots=False)

    ccdicCIII      = uves.load_picklefile(picklefileCIII)
    ccdicCIV       = uves.load_picklefile(picklefileCIV)


    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Grabbing TDOSE spectra and loading object info ')
    #specdir    = '/Volumes/DATABCKUP1/TDOSEextractions/171201_TDOSEextraction/Modelimg/tdose_spectra/'
    specdir    = '/Volumes/DATABCKUP1/TDOSEextractions/180824_TDOSEextraction_LAEs60fields/modelimg/tdose_spectra/'
    if runonallspecs:
        specs_all  = glob.glob(specdir+'tdose_spectrum_candels-*.fits')
    else:
        specs_all  = [specdir+'tdose_spectrum_candels-cdfs-04_modelimg_0104014050-0104014050.fits',
                      specdir+'tdose_spectrum_candels-cdfs-15_modelimg_0115003085-0115003085.fits',
                      specdir+'tdose_spectrum_candels-cdfs-06_modelimg_0106004019-0106004019.fits',
                      specdir+'tdose_spectrum_candels-cdfs-25_modelimg_0125042115-0125042115.fits']

    uvlinesdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'
    z_all      = afits.open(uvlinesdir+'LAEinfo.fits')[1].data['redshift']
    id_all     = afits.open(uvlinesdir+'LAEinfo.fits')[1].data['id']
    try:
        vshift_all = afits.open(uvlinesdir+'LAEinfo.fits')[1].data['red_peak_shift_AV17_kms']
    except:
        vshift_all = afits.open(uvlinesdir+'LAEinfo.fits')[1].data['red_peak_shift_V18_kms']

    specs     = []
    objzs     = []
    vshift    = []

    for spec in specs_all:
        specid  = int(spec.split('_')[-1].split('-')[0])
        goodent = np.where(id_all == specid)
        if len(goodent) != 1:
            sys.exit(' ERROR: found '+str(len(goodent))+' matches to id = '+str(specid))

        zobj = z_all[goodent][0]
        if (zobj > zrange[0]) & (zobj < zrange[1]):

            if specificobj is not None:
                if specid not in specificobj: # skipping objects not in "specificobj" list
                    continue

            specs.append(spec)
            objzs.append(zobj)
            vshift.append(vshift_all[goodent][0])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Grabbing templates ')
    #temps    = glob.glob(templatedir+'uves_felis_template_'+line+'doublet_sig_*_flux'+line+'1_1p0_flux*.fits')
    temps    = glob.glob(templatedir+'uves_felis_template_'+line+'*_sig_*.fits')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if lamwidth_restframe == 'dvoffset':
        if verbose: print(' - Getting wavelength range width (obs frame) corresponding to systemic velocity offset from AV15')
        lamwidth  = []
        wave_rest = [wave_restframe]*len(specs)
        for ss, spec in enumerate(specs):
            lam_obs, lam_offset, dlam = kbs.velocityoffset2dwave(objzs[ss],vshift[ss],wave_rest[ss])
            if dlam > 2.0:
                lamwidth.append(dlam*5.0)
            else:
                lamwidth.append(60.0)
    else:
        if verbose: print(' - Convert fixed rest-frame wavelength width to obs frame widths')
        lamwidth           = []
        wave_rest = [wave_restframe]*len(specs)
        for ss, spec in enumerate(specs):
            waveobs_low  = (wave_rest[ss]-lamwidth_restframe) * (1.0 + objzs[ss-1])
            waveobs_high = (wave_rest[ss]+lamwidth_restframe) * (1.0 + objzs[ss-1])
            lamwidth.append( (waveobs_high - waveobs_low) / 2. )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Crosscorrelating templates to spectra using FELIS')
    picklefile = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/MUSEWideLAEs_CCresults'+datestr+'_'+line+'_RENAME_.pkl'
    if generateplots:
        plotdir            = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/MUSEwideLAE_FELISplots/'
        plot_allCCresults  = True
    else:
        plotdir            = None
        plot_allCCresults  = False
    ccdic      = felis.match_templates2specs(temps,specs,objzs,picklefile,wavewindow=lamwidth,plotdir=plotdir,
                                             wavecen_restframe=wave_rest,vshift=vshift,min_template_level=1e-4,
                                             plot_allCCresults=plot_allCCresults,subtract_spec_median=subtract_spec_median)
    return picklefile
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_FELISmatchOutput_OLD(picklefile,line='CIII',verbose=True,
                          plotdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/MUSEwideLAE_FELISplots/'):
    """
    Producing analytic plots based on a pickled output from match_templates2specs()

    --- EXAMPLE OF USE ---
    uves.plot_FELISmatchOutput('MUSEWideLAEs_CCresults180325_CIII_9templaterun.pkl',line='CIII')
    uves.plot_FELISmatchOutput('MUSEWideLAEs_CCresults180325_CIV_9templaterun.pkl',line='CIV')

    """
    CCdic   = felis.load_picklefile(picklefile)

    Nspecs  = len(CCdic.keys())
    if verbose: print(' - Loaded pickle file and found cross-correlation results for '+str(Nspecs)+' spectra')

    cr_arrmax        = []
    cn_arrmax        = []
    besttemp_fratios = []
    besttemp_sigmas  = []
    besttemp_zmax    = []
    v_offsetCC       = []
    v_offsetAH       = []

    NaNcount     = 0
    for ss, spec in enumerate(CCdic.keys()):
        cr_arr = CCdic[spec]['ccresultsarray']
        cr_arr = cr_arr[np.isfinite(cr_arr)]

        if len(cr_arr) == 0:
            if verbose: print(' WARNING No finite values in cross-correlation array for \n   '+spec)
            NaNcount = NaNcount + 1.0
            continue
        else:
            cr_arr = CCdic[spec]['ccresultsarray'] # need to re-define to get 2D
            cn_arr = CCdic[spec]['ccnormarray']    # need to re-define to get 2D

            try:
                cr_maxent = np.where(cr_arr == np.max(cr_arr))
                cn_maxent = np.where(cn_arr == np.max(cn_arr))

                if cn_arr[cn_maxent][0] < 0.6:
                    continue
                else:
                    cr_arrmax.append(cr_arr[cr_maxent][0])
                    cn_arrmax.append(cn_arr[cn_maxent][0])
            except:
                print('---->problems... cr array has no proper max (even though NaN arrays are igored)')
                pdb.set_trace()

            cr_besttemp  = cr_maxent[0][0]
            cn_besttemp  = cn_maxent[0][0]

            if cr_besttemp != cn_besttemp:
                print('---->problems... best temp do not match ')
                print('     '+spec)
                print('     Bestent: cr '+str(cr_maxent)+'  and cn '+str(cn_maxent)+'\n')
                #pdb.set_trace()

            besttemp = CCdic[spec]['templatevec'][cn_besttemp]

            besttempt_hdr = afits.open(besttemp)[1].header

            besttemp_fratios.append(besttempt_hdr['F'+line+'1_4']/besttempt_hdr['F'+line+'2_4'])
            besttemp_sigmas.append(besttempt_hdr['F'+line+'1_2'])

            zSysCC = CCdic[spec]['zCCmaxvec'][cn_besttemp]
            zLya    = CCdic[spec]['zLya']
            cc      = 299792.458 # km/s
            v_off   = ( (zLya - zSysCC) / (zSysCC + 1.0) ) * cc

            besttemp_zmax.append(zSysCC)
            v_offsetCC.append(v_off)
            v_offsetAH.append(CCdic[spec]['vshift'])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = plotdir+picklefile.split('/')[-1].replace('.pkl','_CCnormmax_hist.pdf')
    if verbose: print(' - Setting up and generating plot')
    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.3,left=0.1, right=0.97, bottom=0.10, top=0.95)
    Fsize    = 12
    lthick   = 2
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    hist = plt.hist(cn_arrmax,color="r",bins=np.arange(-1,1,0.05),histtype="step",lw=1,label=r'')

    plt.xlabel(' Maximum value of normalized cross-correlation for '+str(Nspecs-NaNcount)+' spectra ('+str(NaNcount)+
               ' spectra with all-NaN cross correlation solutions)')
    # plt.ylabel(' S/N ')

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = plotdir+picklefile.split('/')[-1].replace('.pkl','_Fratio_hist.pdf')
    if verbose: print(' - Setting up and generating plot')
    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.3,left=0.1, right=0.97, bottom=0.10, top=0.95)
    Fsize    = 12
    lthick   = 2
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    hist = plt.hist(besttemp_fratios,color="r",bins=30,histtype="step",lw=1,label=r'')

    if line == 'CIII':
        plt.xlabel(' Best matched template line flux ratio (CIII1907/CIII1909) ')
    elif line == 'CIV':
        plt.xlabel(' Best matched template line flux ratio (CIV1548/CIV1551) ')
    # plt.ylabel(' S/N ')

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = plotdir+picklefile.split('/')[-1].replace('.pkl','_linesigma_hist.pdf')
    if verbose: print(' - Setting up and generating plot')
    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.3,left=0.1, right=0.97, bottom=0.10, top=0.95)
    Fsize    = 12
    lthick   = 2
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    hist = plt.hist(besttemp_sigmas,color="r",bins=30,histtype="step",lw=1,label=r'')

    plt.xlabel(' Best matched template line width (Gauss sigma) of '+line+' components [A]')

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = plotdir+picklefile.split('/')[-1].replace('.pkl','_voffset_hist.pdf')
    if verbose: print(' - Setting up and generating plot')
    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.3,left=0.1, right=0.97, bottom=0.10, top=0.95)
    Fsize    = 12
    lthick   = 2
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    hist = plt.hist(v_offsetCC,color="r",bins=np.arange(-1000,1000,10.0),histtype="step",lw=1,label=r'Cross-Corr. prediction')
    hist = plt.hist(v_offsetAH,color="k",bins=np.arange(-1000,1000,10.0),histtype="step",lw=1,label=r'Verhamme prediction')

    plt.xlabel(' Velocity shift [km/s] ')

    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)
    leg.get_frame().set_alpha(0.7)

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_FELISmatchOutput(picklefile,line='CIII',verbose=True,S2Ncut=3,  # only consider CC detections with S/N>S2Ncut
                          plotdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/MUSEwideLAE_FELISplots/',
                          zspecISzLya=False):
    """
    Producing analytic plots based on a pickled output from match_templates2specs()

    --- EXAMPLE OF USE ---
    uves.plot_FELISmatchOutput('MUSEWideLAEs_CCresults180325_CIII_9templaterun.pkl',line='CIII')
    uves.plot_FELISmatchOutput('MUSEWideLAEs_CCresults180325_CIV_9templaterun.pkl',line='CIV')

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Ensure compatibility with FELIS output dictionaries from before 180912
    zkey = 'zspec'
    if zspecISzLya: zkey = 'zLya'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    CCdic   = felis.load_picklefile(picklefile)

    Nspecs  = len(CCdic.keys())
    if verbose: print(' - Loaded pickle file and found cross-correlation results for '+str(Nspecs)+' spectra')

    spec_vec         = []
    max_S2N_vec      = []
    flux_vec         = []
    variance_vec     = []
    besttemp_fratios = []
    besttemp_sigmas  = []
    besttemp_zmax    = []
    v_offsetCC       = []
    v_offsetAH       = []

    for ss, spec in enumerate(CCdic.keys()):
        maxS2N   = np.max(CCdic[spec]['S2NCCmaxvec'])

        if maxS2N > S2Ncut:
            best_ent = np.where(CCdic[spec]['S2NCCmaxvec'] == maxS2N)[0]
            if len(best_ent) > 1:
                print('------> WARNING: More than one pixel with maximum S/N = '+str(maxS2N)+' -> choosing "first" template')
                best_ent = best_ent[0]
            besttemp = CCdic[spec]['templatevec'][best_ent][0]

            cc_best_S2N      = CCdic[spec]['ccresultsarr_S2N'][best_ent][0]
            cc_best_flux     = CCdic[spec]['ccresultsarray_flux'][best_ent][0]
            cc_best_variance = CCdic[spec]['ccresultsarray_variance'][best_ent][0]

            SNmax_ent        = np.where(cc_best_S2N == maxS2N)[0]

            if len(SNmax_ent) == 0:
                print('------> WARNING: No match to S/N in S/N vector. skipping data for key:')
                print('       '+spec)
            else:
                if len(np.atleast_1d(cc_best_flux)) == 1:
                    print('WARNING: CC flux vector only contains 1 value for '+spec)
                    continue
                else:
                    best_flux        = cc_best_flux[SNmax_ent]
                    best_variance    = cc_best_variance[SNmax_ent]

                besttempt_hdr = afits.open(besttemp)[1].header
                if 'fluxratio' in besttemp:
                    besttemp_fratios.append(besttempt_hdr['F'+line+'1_4']/besttempt_hdr['F'+line+'2_4'])
                    besttemp_sigmas.append(besttempt_hdr['F'+line+'1_2'])
                else:
                    besttemp_sigmas.append(besttempt_hdr['F'+line+'_2'])

                zSysCC  = CCdic[spec]['zCCmaxvec'][best_ent]
                zLya    = CCdic[spec][zkey]
                cc      = 299792.458 # km/s
                v_off   = ( (zLya - zSysCC) / (zSysCC + 1.0) ) * cc

                spec_vec.append(spec)
                max_S2N_vec.append(maxS2N)
                flux_vec.append(best_flux[0])
                variance_vec.append(best_variance[0])
                besttemp_zmax.append(zSysCC[0])
                v_offsetCC.append(v_off[0])
                v_offsetAH.append(CCdic[spec]['vshift'])

    if len(max_S2N_vec) == 0:
        sys.exit(' No detections passing the S/N cut of '+str(S2Ncut))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if 'fluxratio' in besttemp:
        plotname = plotdir+picklefile.split('/')[-1].replace('.pkl','_Fratio_hist.pdf')
        if verbose: print(' - Setting up and generating plot')
        fig = plt.figure(figsize=(5, 5))
        fig.subplots_adjust(wspace=0.1, hspace=0.3,left=0.1, right=0.97, bottom=0.10, top=0.95)
        Fsize    = 12
        lthick   = 2
        marksize = 4
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()

        hist = plt.hist(besttemp_fratios,color="r",bins=30,histtype="step",lw=1,label=r'')

        if line == 'CIII':
            lineratiostring = '(CIII1907/CIII1909)'
        elif line == 'CIV':
            lineratiostring = '(CIV1548/CIV1551)'
        elif line == 'OIII':
            lineratiostring = '(OIII1661/OIII1666)'
        elif line == 'NV':
            lineratiostring = '(NV1239/NV1243)'

        plt.xlabel(' Template line flux ratio '+lineratiostring)

        if verbose: print('   Saving plot to '+plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = plotdir+picklefile.split('/')[-1].replace('.pkl','_linesigma_hist.pdf')
    if verbose: print(' - Setting up and generating plot')
    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.3,left=0.1, right=0.97, bottom=0.10, top=0.95)
    Fsize    = 12
    lthick   = 2
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    hist = plt.hist(besttemp_sigmas,color="r",bins=30,histtype="step",lw=1,label=r'')

    plt.xlabel(' Template line width (Gauss sigma) of '+line+' components [A]')

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = plotdir+picklefile.split('/')[-1].replace('.pkl','_voffset_hist.pdf')
    if verbose: print(' - Setting up and generating plot')
    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.3,left=0.1, right=0.97, bottom=0.10, top=0.95)
    Fsize    = 12
    lthick   = 2
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    hist = plt.hist(v_offsetCC,color="r",bins=np.arange(-1000,1000,10.0),histtype="step",lw=1,label=r'Cross-Corr. prediction')
    hist = plt.hist(v_offsetAH,color="k",bins=np.arange(-1000,1000,10.0),histtype="step",lw=1,label=r'Verhamme prediction')

    plt.xlabel(' Velocity shift wrt. Ly$\\alpha$ [km/s] ')

    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)
    leg.get_frame().set_alpha(0.7)

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = plotdir+picklefile.split('/')[-1].replace('.pkl','_S2Nmax_hist.pdf')
    if verbose: print(' - Setting up and generating plot')
    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.3,left=0.1, right=0.97, bottom=0.10, top=0.95)
    Fsize    = 12
    lthick   = 2
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    hist = plt.hist(max_S2N_vec,color="r",bins=30,histtype="step",lw=1,label=r'')

    plt.xlabel(' S/N ($>$'+str(S2Ncut)+' only)')

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = plotdir+picklefile.split('/')[-1].replace('.pkl','_FluxmaxVSS2Nmax.pdf')
    if verbose: print(' - Setting up and generating plot')
    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.3,left=0.15, right=0.97, bottom=0.10, top=0.95)
    Fsize    = 12
    lthick   = 2
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    yerrbar = np.sqrt(variance_vec)
    plt.errorbar(max_S2N_vec,flux_vec, yerr=yerrbar,fmt='o', markerfacecolor='red', markeredgecolor='black',
                 ecolor='r', capthick=2)

    plt.xlabel(' S/N ($>$'+str(S2Ncut)+' only)')

    if line == 'CIII':
        tempindicator = 'CIII1907 + CIII1909'
    elif line == 'CIV':
        tempindicator = 'CIV1548 + CIV1551'
    elif line == 'OIII':
        tempindicator = 'OIII1661 + OIII1666'
    elif line == 'NV':
        tempindicator = 'NV1239 + NV1243'
    elif line == 'HEII':
        tempindicator = 'HeII1640'

    plt.ylabel(' $F_\\textrm{tot}$/[1e-20cgs] scaling for ('+tempindicator+')')

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def check_neighbors(ids=[214063213],
                    sourcecatdir='/Volumes/DATABCKUP1/TDOSEextractions/180822_TDOSEextraction_LAEs60fields/tdose_sourcecats/',
                    modeldir='/Volumes/DATABCKUP1/TDOSEextractions/MW_LAEs_JKgalfitmodels/'):
    """

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    uves.check_neighbors()

    """

    for id in ids:

        if str(id).startswith('2'):
            field = 'cosmos'
        else:
            field = 'cdfs'

        fieldno = str(id)[1:3]

        sc_generated = sourcecatdir+'tdose_sourcecat_LAEs_candels-'+field+'-'+str(fieldno)+'_id'+\
                       str(id)+'_cutout2p0x2p0arcsec.fits'
        sc_JKmodel   = modeldir+'model_acs_814w_candels-'+field+'-'+str(fieldno)+'_cut_v1.0_id'+\
                       str(id)+'_cutout2p0x2p0arcsec_sourcecatalog.fits'

        sc_gen = afits.open(sc_generated)[1].data
        sc_mod = afits.open(sc_JKmodel)[1].data

        print(' - - - - - - - - - - '+str(id)+' - - - - - - - - - - - ')
        print(" - SOURCE CATALOG GENERATED FROM CROSSMATCH TO MASTER SOURCE CAT")
        print('   '+str(sc_gen['PARENT_ID']))
        print('   '+str(sc_gen['ID']))
        print(" - SOURCE CATALOG CONTENT BASED ON JOSIE'S MODEL")
        print('   '+str(sc_mod['PARENT_ID']))
        print('   '+str(sc_mod['ID']))
        print(' ')
        try:
            print('   diff(PARENT_ID): '+str(sc_gen['PARENT_ID']-sc_mod['PARENT_ID']))
        except:
            print('   WARNING: mismatch in PARENT_ID source cats')
        try:
            print('   diff(ID): '+str(sc_gen['ID']-sc_mod['ID']))
        except:
            print('   WARNING: mismatch in ID source cats')
        print(' ')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def calc_lineratios_fromsummaryfiles(summaryfiles,lineindicators,outputfile, Nsigmalimits=3, vetfelis_output=None, verbose=True):
    """
    Function to calculate the flux and line ratios for a set of summary files
    containing the results from FELIS template matches to TDOSE spectra.

    Based (partially) on uves.calculatelineratios() below.

    --- INPUT ---
    summaryfiles     A list of summaryfiles to generate flux ratio output for
    lineindicators   List of strings indicating the content of the summary files (i.e. the lines that were
                     searched for with the FELIS template match). These are used to name output columns.
    outputfile       The ascii file to write results to.
    verbose          Toggle verbosity

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves

    outfile        = outputdir+'fluxratios190919_'+outkeystr+'.txt'
    summaryfiles   = glob.glob(outputdir+'CCresults_summary/CCresults_summary_template*_'+outkeystr+'.txt')
    lineindicators = [sf.split('template')[-1].split('_')[0] for sf in summaryfiles]

    fluxratiodat   = uves.calc_lineratios_fromsummaryfiles(summaryfiles,lineindicators,outfile)

    """
    if len(summaryfiles) != len(lineindicators):
        sys.exit('\nThe '+str(len(summaryfiles))+' provided do not match the '+str(len(lineindicators))+' provided \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading content of the '+str(lineindicators)+' summary files provided')
    dic_summarydat = {}
    fmt = '12a,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,300a,300a'
    for ss, sfile in enumerate(summaryfiles):
        summarydat  = np.genfromtxt(sfile,skip_header=25,dtype=fmt,comments='#',names=True)
        dic_summarydat[lineindicators[ss]] = summarydat

    spectra   = []
    for ll in lineindicators:
        spectra = spectra+list(dic_summarydat[ll]['spectrum'])
    spectra = np.unique(np.array(spectra))

    pointings = []
    ids       = []
    for spec in spectra:
        pointings.append(spec.split('/tdose_spectrum_')[-1].split('-full')[0])
        ids.append(int(spec.split('_')[-1].split('.fit')[0]))
    pointings = np.asarray(pointings)
    ids       = np.asarray(ids)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Initializing the output file: \n   '+outputfile)
    fout = open(outputfile,'w')
    fout.write('# Flux and line ratios estimated based on FELIS template match results summarized in:\n')
    fout.write('# '+str(summaryfiles)+'\n')
    if vetfelis_output is None:
        fout.write('# \n')
        fout.write('# No limits are stored in file - s2n values are provided instead to determine limits in post-processing. \n')
        fout.write('# See for instance uves.plot_lineratios_fromsummaryfiles_wrapper() \n')
    else:
        fout.write('# Using FELIS vetting info from '+vetfelis_output+'\n')
        fout.write('# Limits are provided according to FELIS vetting, i.e., "untrustworthy" lines and non-detections '
                   'are quoted as '+str(Nsigmalimits)+'sigma upper limits. '
                   'Note that all other error bars are also '+str(Nsigmalimits)+'sigma errors. \n')
        fout.write('# The s2n values above this threshold are also provided for post-processing. See for instance uves.plot_lineratios_fromsummaryfiles_wrapper() \n')
    fout.write('# \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Determine columns to fill in output')
    fluxratiodic = collections.OrderedDict()
    fluxratiodic['id'] = np.array([])
    fluxratiodic['pointing'] = pointings
    for ll, numerator_line in enumerate(lineindicators):
        fluxratiodic['f_'+numerator_line]       = np.array([])
        fluxratiodic['ferr_'+numerator_line]    = np.array([])
        fluxratiodic['s2n_'+numerator_line]     = np.array([])
        fluxratiodic['sigma_'+numerator_line]   = np.array([])
        fluxratiodic['vshift_'+numerator_line]  = np.array([])

        if len(dic_summarydat[numerator_line]['id']) == 0:
            continue
        elif dic_summarydat[numerator_line]['Fratio_temp'][0] != 0:
            fluxratiodic['f_'+numerator_line+'1']      = np.array([])
            fluxratiodic['ferr_'+numerator_line+'1']   = np.array([])
            fluxratiodic['f_'+numerator_line+'2']      = np.array([])
            fluxratiodic['ferr_'+numerator_line+'2']   = np.array([])
            fluxratiodic['FR_'+numerator_line+'1'+numerator_line+'2']     = np.array([])
            fluxratiodic['FRerr_'+numerator_line+'1'+numerator_line+'2']  = np.array([])
            fluxratiodic['FRs2n_'+numerator_line+'1'+numerator_line+'2']  = np.array([])

        for kk, denominator_line in enumerate(lineindicators):
            if numerator_line == denominator_line:
                continue
            elif len(dic_summarydat[denominator_line]['id']) == 0:
                continue
            else:
                fluxratiodic['FR_'+numerator_line+denominator_line]     = np.array([])
                fluxratiodic['FRerr_'+numerator_line+denominator_line]  = np.array([])
                fluxratiodic['FRs2n_'+numerator_line+denominator_line]  = np.array([])

                if (dic_summarydat[numerator_line]['Fratio_temp'][0] != 0):
                    fluxratiodic['FR_'+numerator_line+'1'+denominator_line]     = np.array([])
                    fluxratiodic['FRerr_'+numerator_line+'1'+denominator_line]  = np.array([])
                    fluxratiodic['FRs2n_'+numerator_line+'1'+denominator_line]  = np.array([])
                    fluxratiodic['FR_'+numerator_line+'2'+denominator_line]     = np.array([])
                    fluxratiodic['FRerr_'+numerator_line+'2'+denominator_line]  = np.array([])
                    fluxratiodic['FRs2n_'+numerator_line+'2'+denominator_line]  = np.array([])
                if (dic_summarydat[denominator_line]['Fratio_temp'][0] != 0):
                    fluxratiodic['FR_'+numerator_line+denominator_line+'1']     = np.array([])
                    fluxratiodic['FRerr_'+numerator_line+denominator_line+'1']  = np.array([])
                    fluxratiodic['FRs2n_'+numerator_line+denominator_line+'1']  = np.array([])
                    fluxratiodic['FR_'+numerator_line+denominator_line+'2']     = np.array([])
                    fluxratiodic['FRerr_'+numerator_line+denominator_line+'2']  = np.array([])
                    fluxratiodic['FRs2n_'+numerator_line+denominator_line+'2']  = np.array([])
                if (dic_summarydat[numerator_line]['Fratio_temp'][0] != 0) & \
                        (dic_summarydat[denominator_line]['Fratio_temp'][0] != 0):
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
    if vetfelis_output is not None:
        if verbose: print(' - FELIS vetting results provided so loading these - only vetted objects stored in output.')
        dat_vetfelis  = np.genfromtxt(vetfelis_output,  dtype=None,comments='#',names=True,skip_header=12)
        Nobj_fvet     = len(np.unique(dat_vetfelis['id']))
        if verbose: print('   Found '+str(Nobj_fvet)+' objects in the FELIS vetting output to be stored to line ratios output ')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Filling the columns with data ')
    fluxratioarray = np.zeros([len(ids),Ncols])*np.nan
    for ii, id in enumerate(ids):
        if verbose:
            infostr = '   Filling output for object '+str(id)+'   (id '+str(ii+1)+'/'+str(len(ids))+')'
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()
        fluxratioarray[ii,0] = float(id)

        if vetfelis_output is not None:
            if 'UVESselectio' in pointings[ii]:
                ent_vet_felis = np.where((dat_vetfelis['id'] == id))[0]
            else:
                ent_vet_felis = np.where((dat_vetfelis['id'] == id) & (dat_vetfelis['pointing'] == pointings[ii]))[0]

            if len(ent_vet_felis) > 1:
                sys.exit('\n WARNING: object '+str(id)+' in '+pointings[ii]+
                         ' exisits multiple times in FELIS vetting output; sort this out befor calculating lineratios...')
            elif len(ent_vet_felis) == 0:
                if verbose: print('\n           Object '+str(id)+' in '+pointings[ii]+' not found in FELIS vetting output')
                # creating dummy obj_vetfelis entry filled with 9999s
                obj_vetfelis = dat_vetfelis[0].copy()#*0.0 (101009024, 'cdfs-01', 99, 99, 99, 99, 99, 99, -99)
                obj_vetfelis['id']       = id
                obj_vetfelis['pointing'] = pointings[ii]
                for col in obj_vetfelis.dtype.names[2:]:
                    obj_vetfelis[col] = 9999 # no inspection exists; obj_vetfelis made by hand
            else:
                obj_vetfelis = dat_vetfelis[ent_vet_felis]

        for ll, numerator_line in enumerate(lineindicators):
            numerator_dat = dic_summarydat[numerator_line]
            ent_num       = np.where((numerator_dat['id'].astype(int) == id) & (numerator_dat['spectrum'] == spectra[ii]))[0]
            if len(ent_num) == 0:
                # if verbose: print('       object '+str(id)+' in '+pointings[ii]+' not found in '+
                #                   numerator_line+' summary (numerator)')
                continue

            f_num    = numerator_dat['Ftot_FELIS_S2Nmax'][ent_num]
            ferr_num = numerator_dat['Ftot_FELIS_S2Nmax_err'][ent_num]

            fluxratioarray[ii,colents['f_'+numerator_line]]      = f_num
            fluxratioarray[ii,colents['ferr_'+numerator_line]]   = ferr_num
            fluxratioarray[ii,colents['s2n_'+numerator_line]]    = numerator_dat['FELIS_S2Nmax'][ent_num]
            fluxratioarray[ii,colents['sigma_'+numerator_line]]  = numerator_dat['sigma_temp_ang_rf'][ent_num]
            fluxratioarray[ii,colents['vshift_'+numerator_line]] = numerator_dat['vshift_CCmatch'][ent_num]

            trust_numerator = 'yes'

            if vetfelis_output is not None:
                ferr_num = ferr_num * Nsigmalimits

                if (obj_vetfelis['trust'+numerator_line] == 0) or \
                        (obj_vetfelis['trust'+numerator_line] == 9) or \
                        (obj_vetfelis['trust'+numerator_line] == 99) or \
                        (obj_vetfelis['trust'+numerator_line] == 9999) or \
                        (fluxratioarray[ii,colents['s2n_'+numerator_line]] < Nsigmalimits):

                    if (obj_vetfelis['trust'+numerator_line] == 9999) & \
                            (fluxratioarray[ii,colents['s2n_'+numerator_line]] > Nsigmalimits):
                        if verbose: print('\033[91m        ... and '+numerator_line+' has S/N > '+
                                          str(Nsigmalimits)+' (vshift = '+
                                          str(fluxratioarray[ii,colents['vshift_'+numerator_line]])+
                                          '); make sure to include it in the vetting results as it will be set to '+
                                              str(Nsigmalimits)+'sigma limits here ...\033[0m')

                    fluxratioarray[ii,colents['f_'+numerator_line]]      = ferr_num
                    fluxratioarray[ii,colents['ferr_'+numerator_line]]   = 99
                    if ferr_num != 0.0:
                        fluxratioarray[ii,colents['s2n_'+numerator_line]] = Nsigmalimits
                    fluxratioarray[ii,colents['sigma_'+numerator_line]]  = 99
                    fluxratioarray[ii,colents['vshift_'+numerator_line]] = 99

                    trust_numerator = 'None'

            # if (id == 102014087) & (pointings[ii] == 'cdfs-08') & (numerator_line == 'CIV'):
            #     pdb.set_trace()

            if (dic_summarydat[numerator_line]['Fratio_temp'][0] != 0) & (trust_numerator != 'None'):
                # error on q=|B|x where |B| is known exact is just dq=|B|dx
                f1_num      = numerator_dat['Ftot_FELIS_S2Nmax'][ent_num]     / (1 + 1/numerator_dat['Fratio_temp'][ent_num])
                f1err_num   = numerator_dat['Ftot_FELIS_S2Nmax_err'][ent_num] / (1 + 1/numerator_dat['Fratio_temp'][ent_num])
                f2_num      = f1_num    / numerator_dat['Fratio_temp'][ent_num]
                f2err_num   = f1err_num / numerator_dat['Fratio_temp'][ent_num]
                FR12, FR12err = lce.set_ratios('gooddoublet','gooddoublet',f1_num,f1err_num,f2_num,f2err_num)

                fluxratioarray[ii,colents['f_'+numerator_line+'1']]                         = f1_num
                fluxratioarray[ii,colents['ferr_'+numerator_line+'1']]                      = f1err_num
                fluxratioarray[ii,colents['f_'+numerator_line+'2']]                         = f2_num
                fluxratioarray[ii,colents['ferr_'+numerator_line+'2']]                      = f2err_num
                fluxratioarray[ii,colents['FR_'+numerator_line+'1'+numerator_line+'2']]     = FR12
                fluxratioarray[ii,colents['FRerr_'+numerator_line+'1'+numerator_line+'2']]  = FR12err
                fluxratioarray[ii,colents['FRs2n_'+numerator_line+'1'+numerator_line+'2']]  = FR12/FR12err

            for kk, denominator_line in enumerate(lineindicators):
                denominator_dat = dic_summarydat[denominator_line]

                ent_denom       = np.where((denominator_dat['id'].astype(int) == id) & (denominator_dat['spectrum'] == spectra[ii]))[0]
                if len(ent_denom) > 1:
                    sys.exit(' id pointing combination '+str(id)+' '+pointings[ii]+' appears '+str(len(ent_denom))+
                             ' times in '+numerator_line+' summary -> that should not happen!')

                if len(ent_denom) == 0:
                    # if verbose: print('       object '+str(id)+' in '+pointings[ii]+' not found in '+
                    #                   denominator_line+' summary (numerator)')
                    continue

                if numerator_line == denominator_line:
                    continue
                else:
                    trust_denominator = 'yes'
                    f_denom    = denominator_dat['Ftot_FELIS_S2Nmax'][ent_denom]
                    ferr_denom = denominator_dat['Ftot_FELIS_S2Nmax_err'][ent_denom]

                    if vetfelis_output is not None:
                        ferr_denom = ferr_denom * Nsigmalimits

                        if (obj_vetfelis['trust'+denominator_line] == 0) or \
                                (obj_vetfelis['trust'+denominator_line] == 9) or \
                                (obj_vetfelis['trust'+denominator_line] == 99) or \
                                (f_denom/ferr_denom < Nsigmalimits):
                            f_denom           = ferr_denom
                            ferr_denom        = 99
                            trust_denominator = 'None'

                    if (trust_numerator == 'None') & (trust_denominator == 'None'):
                        continue

                    FR, FRerr = lce.set_ratios(trust_numerator,trust_denominator,f_num,ferr_num,f_denom,ferr_denom)
                    fluxratioarray[ii,colents['FR_'+numerator_line+denominator_line]]     = FR
                    fluxratioarray[ii,colents['FRerr_'+numerator_line+denominator_line]]  = FRerr
                    fluxratioarray[ii,colents['FRs2n_'+numerator_line+denominator_line]]  = FR/FRerr

                    if (dic_summarydat[numerator_line]['Fratio_temp'][0] != 0) & (trust_numerator != 'None'):
                        FR, FRerr = lce.set_ratios(trust_numerator,trust_denominator,f1_num,f1err_num,f_denom,ferr_denom)
                        fluxratioarray[ii,colents['FR_'+numerator_line+'1'+denominator_line]]     = FR
                        fluxratioarray[ii,colents['FRerr_'+numerator_line+'1'+denominator_line]]  = FRerr
                        fluxratioarray[ii,colents['FRs2n_'+numerator_line+'1'+denominator_line]]  = FR/FRerr

                        FR, FRerr = lce.set_ratios(trust_numerator,trust_denominator,f2_num,f2err_num,f_denom,ferr_denom)
                        fluxratioarray[ii,colents['FR_'+numerator_line+'2'+denominator_line]]     = FR
                        fluxratioarray[ii,colents['FRerr_'+numerator_line+'2'+denominator_line]]  = FRerr
                        fluxratioarray[ii,colents['FRs2n_'+numerator_line+'2'+denominator_line]]  = FR/FRerr

                    if (dic_summarydat[denominator_line]['Fratio_temp'][0] != 0)  & (trust_denominator != 'None'):
                        f1_denom    = denominator_dat['Ftot_FELIS_S2Nmax'][ent_denom]     / \
                                      (1 + 1/denominator_dat['Fratio_temp'][ent_denom])
                        f1err_denom = denominator_dat['Ftot_FELIS_S2Nmax_err'][ent_denom] / \
                                      (1 + 1/denominator_dat['Fratio_temp'][ent_denom])
                        f2_denom    = f1_denom    / denominator_dat['Fratio_temp'][ent_denom]
                        f2err_denom = f1err_denom / denominator_dat['Fratio_temp'][ent_denom]

                        FR, FRerr = lce.set_ratios(trust_numerator,'gooddoublet',f_num,ferr_num,f1_denom,f1err_denom)
                        fluxratioarray[ii,colents['FR_'+numerator_line+denominator_line+'1']]     = FR
                        fluxratioarray[ii,colents['FRerr_'+numerator_line+denominator_line+'1']]  = FRerr
                        fluxratioarray[ii,colents['FRs2n_'+numerator_line+denominator_line+'1']]  = FR/FRerr

                        FR, FRerr = lce.set_ratios(trust_numerator,'gooddoublet',f_num,ferr_num,f2_denom,f2err_denom)
                        fluxratioarray[ii,colents['FR_'+numerator_line+denominator_line+'2']]     = FR
                        fluxratioarray[ii,colents['FRerr_'+numerator_line+denominator_line+'2']]  = FRerr
                        fluxratioarray[ii,colents['FRs2n_'+numerator_line+denominator_line+'2']]  = FR/FRerr

                    if (dic_summarydat[numerator_line]['Fratio_temp'][0] != 0) & \
                        (dic_summarydat[denominator_line]['Fratio_temp'][0] != 0) & \
                            (trust_numerator != 'None') & (trust_denominator != 'None'):
                        FR, FRerr = lce.set_ratios('gooddoublet','gooddoublet',f1_num,f1err_num,f1_denom,f1err_denom)
                        fluxratioarray[ii,colents['FR_'+numerator_line+'1'+denominator_line+'1']]     = FR
                        fluxratioarray[ii,colents['FRerr_'+numerator_line+'1'+denominator_line+'1']]  = FRerr
                        fluxratioarray[ii,colents['FRs2n_'+numerator_line+'1'+denominator_line+'1']]  = FR/FRerr

                        FR, FRerr = lce.set_ratios('gooddoublet','gooddoublet',f1_num,f1err_num,f2_denom,f2err_denom)
                        fluxratioarray[ii,colents['FR_'+numerator_line+'1'+denominator_line+'2']]     = FR
                        fluxratioarray[ii,colents['FRerr_'+numerator_line+'1'+denominator_line+'2']]  = FRerr
                        fluxratioarray[ii,colents['FRs2n_'+numerator_line+'1'+denominator_line+'2']]  = FR/FRerr

                        FR, FRerr = lce.set_ratios('gooddoublet','gooddoublet',f2_num,f2err_num,f1_denom,f1err_denom)
                        fluxratioarray[ii,colents['FR_'+numerator_line+'2'+denominator_line+'1']]     = FR
                        fluxratioarray[ii,colents['FRerr_'+numerator_line+'2'+denominator_line+'1']]  = FRerr
                        fluxratioarray[ii,colents['FRs2n_'+numerator_line+'2'+denominator_line+'1']]  = FR/FRerr

                        FR, FRerr = lce.set_ratios('gooddoublet','gooddoublet',f2_num,f2err_num,f2_denom,f2err_denom)
                        fluxratioarray[ii,colents['FR_'+numerator_line+'2'+denominator_line+'2']]     = FR
                        fluxratioarray[ii,colents['FRerr_'+numerator_line+'2'+denominator_line+'2']]  = FRerr
                        fluxratioarray[ii,colents['FRs2n_'+numerator_line+'2'+denominator_line+'2']]  = FR/FRerr

    for ll in np.arange(len(ids)):
        outstr = str(int(fluxratioarray[ll,0]))+'  '+pointings[ll]+' '+' '.join([str("%10.4f" % ff) for ff in fluxratioarray[ll,2:]])
        fout.write(outstr+' \n')
    fout.close()
    if verbose: print('\n - Wrote the flux ratio output to \n   '+outputfile)
    fmt = 'i,12a,'+','.join((Ncols-2)*['d'])
    fluxratiodat = np.genfromtxt(outputfile,skip_header=7,dtype=fmt,comments='#',names=True)
    return fluxratiodat

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_lineratios_fromsummaryfiles(lineratiofile, plotbasename, infofile, colorvar_obj='s2n_CIII', point_text=None,
                                     Nsigma=3.0, colorvar_pi='logUs', vshiftmax=1e5, obj2show='all', showlimits=True,
                                     addliteraturevalues = False, litsymboldot=False,
                                     overwrite=False, verbose=True):
    """
    Function to plot the output containing flux ratios generated with uves.calc_lineratios_fromsummary()

    --- INPUT ---


    --- EXAMPLE OF USE ---
    uves.plot_lineratios_fromsummaryfiles(lineratiofile, plotbasename, overwrite=False, verbose=True)

    """
    if verbose: print(' - Loading flux ratio data to plot ')
    fluxratiodatALL = np.genfromtxt(lineratiofile,skip_header=7,dtype='d',comments='#',names=True)

    if verbose: print(' - Loading infofile data for expanded selection and plotting ')
    infofiledat      = afits.open(infofile)[1].data
    infofiledat      = infofiledat[np.where((infofiledat['id']<4.9e8) | (infofiledat['id']>5.9e8))[0]] # ignoring UDF MW mock ids

    if obj2show is 'all':
        showent = np.arange(len(fluxratiodatALL))
    elif obj2show is 'all_nodup':
        idlist2show = infofiledat[(infofiledat['duplicationID'] == 0.0) & (infofiledat['redshift'] <= 4.955)]['id']
        showent     = np.array([])
        for objid in idlist2show:
            objent = np.where( fluxratiodatALL['id'].astype(int) == objid)[0]
            if len(objent) > 0:
                showent = np.append(showent,objent)
            # else:
            #     print(str(objid)+'  '+str(infofiledat[(infofiledat['id'] == objid)]['redshift']))
    elif obj2show.lower() == 'none':
        showent = np.array([0])
    else:
        sample = 'udf10'
        if obj2show is 'goodspec_only':
            ids_badTDOSEspec, ids_goodTDOSEspec = uves.summarize_tdosevetting(returnsample=sample,verbose=verbose)
            idlist2show = ids_goodTDOSEspec
        elif obj2show is 'badspec_only':
            ids_badTDOSEspec, ids_goodTDOSEspec = uves.summarize_tdosevetting(returnsample=sample,verbose=verbose)
            idlist2show = ids_badTDOSEspec
        else:
            idlist2show = obj2show

        showent = np.array([])
        for objid in idlist2show:
            objent = np.where( fluxratiodatALL['id'].astype(int) == objid)[0]
            if len(objent) > 0:
                showent = np.append(showent,objent)

    Nselspec    = len(showent)
    if Nselspec > 0:
        fluxratiodat = fluxratiodatALL[showent.astype(int)]
        if verbose: print(' - '+str(Nselspec)+'/'+str(len(fluxratiodatALL))+' spectra in flux ratio summary satisfies the cuts\n')
    else:
        if verbose: print(' WARNING No flux ratio matches found in summary file satisfying cuts; returning...')
        return

    if colorvar_obj in fluxratiodat.dtype.names:
        cdatvec   = fluxratiodat[colorvar_obj]
    elif colorvar_obj in infofiledat.columns.names:
        cdatvec = np.zeros(len(fluxratiodat['id']))*np.nan
        for ii, id in enumerate(fluxratiodat['id']):
            infoent     = np.where(infofiledat['id'] == int(id))
            cdatvec[ii] = infofiledat[colorvar_obj][infoent]

    if colorvar_obj.lower() == 's2n_ciii':
        cdattype = 's2n_ciii'
    elif colorvar_obj.lower() == 'redshift':
        cdattype = colorvar_obj.lower()
    else:
        cdattype = None

    #------------------------------- Append literature measurements -------------------------------
    if addliteraturevalues:
        litcat  = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/' \
                  'literaturecollection_emissionlinestrengths.fits'
        litdat  = afits.open(litcat)[1].data
        # litdat  = np.genfromtxt(litcat,names=True,skip_header=5,comments='#',dtype=None)
        Nlitobj = len(litdat)
        litarr  = np.array(np.zeros(Nlitobj)*np.nan,fluxratiodat.dtype)
        for lo, litobject in enumerate(litdat):
            for litcol in litdat.dtype.names:
                if litcol in fluxratiodat.dtype.names:
                    litarr[litcol][lo] = litdat[litcol][lo]

        fluxratiodat = np.hstack((fluxratiodat,litarr))

        if colorvar_obj.lower() == 'redshift':
            cdatvec  = np.hstack((cdatvec,litdat['redshift']))
    # pdb.set_trace() #  litdat['FR_CIVOIII'][-3:], litdat['FR_CIVHeII'][-3:],  litdat['FRerr_CIVHeII'][-3:],  litdat['FRerr_CIVOIII'][-3:]
    # fluxratiodat['FR_CIVOIII'][-3:], fluxratiodat['FR_CIVHeII'][-3:],  fluxratiodat['FRerr_CIVHeII'][-3:],  fluxratiodat['FRerr_CIVOIII'][-3:]

    #----------------------------------------------------------------------------------------------

    if point_text is not None:
        point_text = fluxratiodat['id'].astype(str)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # NOTE: See uves.get_infodat_plotcols() for column definitions and plot labels

    infocols        = ['lyaew_JKmed',
                       'lyaew_b2',
                       'lyafwhm_a',
                       'lyafwhm_kms',
                       # 'beta_beta2',
                       'absmagUV_-2',
                       'magUV_-2',
                       'absmagUV_median',
                       'magUV_median',
                       'lyaflux',
                       'zsys',
                       'redshift',
                       'peaksep_kms',
                       'peaksep_A',
                       'loglLlya',
                       'R_e']

    info_ranges     = [[-100,440],
                       [-100,440],
                       [-5,30],
                       [-100,990],
                       # [-4,3],
                       [-25,-10],
                       [20,35],
                       [-25,-10],
                       [20,35],
                       [1e-19,1e-15],
                       [0.0,8.0],
                       [0.0,8.0],
                       [-100.,990.],
                       [-5.,25.],
                       [39.8,43.8],
                       [0.,9.]]

    addvaryingbetavalues = True
    if addvaryingbetavalues:
        infocols        = infocols+['lyaew_many',
                           'beta_many',
                           'absmagUV_many',
                           'magUV_many']

        info_ranges     = info_ranges+[[-40,440],
                           [-3,-1],#[-4,3],
                           [-25,-10],
                           [20,35]]

    s2n_range       = [2.0,55.]
    if addliteraturevalues:
        # fluxes_range    = [10,9e7]
        fluxes_range    = [10,8e3]
    else:
        fluxes_range    = [10,8e3]
    ratios_range    = [1e-4,1e3]
    FR_range        = [0.0,3.75]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    linesetlist_lya = []
    for ll, infocol in enumerate(infocols):
        linesetlist_lya.append([infocol, None, 'CIV',   None ,info_ranges[ll], fluxes_range, None])
        linesetlist_lya.append([infocol, None, 'CIII',  None ,info_ranges[ll], fluxes_range, None])
        linesetlist_lya.append([infocol, None, 'OIII',  None ,info_ranges[ll], fluxes_range, None])
        linesetlist_lya.append([infocol, None, 'HeII',  None ,info_ranges[ll], fluxes_range, None])
        # linesetlist_lya.append([infocol, None, 'MgII',  None ,info_ranges[ll], fluxes_range, None])
        linesetlist_lya.append([infocol, None, 'SiIII', None ,info_ranges[ll], fluxes_range, None])
        linesetlist_lya.append([infocol, None, 'NV',    None ,info_ranges[ll], fluxes_range, None])

        linesetlist_lya.append([infocol, None, 'CIV',   'CIII',info_ranges[ll], ratios_range, None])
        linesetlist_lya.append([infocol, None, 'HeII',  'CIII',info_ranges[ll], ratios_range, None])
        linesetlist_lya.append([infocol, None, 'SiIII', 'CIII',info_ranges[ll], ratios_range, None])
        linesetlist_lya.append([infocol, None, 'OIII',  'CIII',info_ranges[ll], ratios_range, None])
        linesetlist_lya.append([infocol, None, 'NV',    'CIII',info_ranges[ll], ratios_range, None])

    Nhistbins = 30
    histaxes  = True
    for lineset in linesetlist_lya:
        uves.plot_lineratios_fromsummaryfiles_vsInfofile(plotbasename,fluxratiodat,lineset,histaxes,Nhistbins,cdatvec,'zmanual',
                                                         Nsigma=Nsigma,point_text=point_text,vshiftmax=vshiftmax,performlinearfit=True,
                                                         ylog=True,xlog=False,addliteraturevalues=addliteraturevalues,
                                                         overwrite=overwrite,verbose=verbose,showlimits=showlimits,litsymboldot=litsymboldot)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    linesetlist_lyaFR = []

    for ll, infocol in enumerate(infocols):
        linesetlist_lyaFR.append([infocol, None, 'FR_NV1NV2',       None ,info_ranges[ll], FR_range, None])
        linesetlist_lyaFR.append([infocol, None, 'FR_CIV1CIV2',     None ,info_ranges[ll], FR_range, None])
        linesetlist_lyaFR.append([infocol, None, 'FR_OIII1OIII2',   None ,info_ranges[ll], FR_range, None])
        linesetlist_lyaFR.append([infocol, None, 'FR_SIIII1SIIII2', None ,info_ranges[ll], FR_range, None])
        linesetlist_lyaFR.append([infocol, None, 'FR_CIII1CIII2',   None ,info_ranges[ll], FR_range, None])

    Nhistbins = 30
    histaxes  = False
    for lineset in linesetlist_lyaFR:
        uves.plot_lineratios_fromsummaryfiles_vsInfofile(plotbasename,fluxratiodat,lineset,histaxes,Nhistbins,cdatvec,'zmanual',
                                                         Nsigma=Nsigma,point_text=point_text,vshiftmax=vshiftmax,performlinearfit=True,
                                                         ylog=False,xlog=False,addliteraturevalues=addliteraturevalues,
                                                         overwrite=overwrite,verbose=verbose,showlimits=showlimits,litsymboldot=litsymboldot)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    linesetlist_fluxesFR = []

    linesetlist_fluxesFR.append(['CIV1'  ,'CIV2'    ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxesFR.append(['OIII1' ,'OIII2'   ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxesFR.append(['MgII1' ,'MgII2'   ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxesFR.append(['NV1'   ,'NV2'     ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxesFR.append(['SiIII1','SiIII2'  ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxesFR.append(['CIII1' ,'CIII2'   ,None,None,fluxes_range, fluxes_range,   None])

    Nhistbins = 30
    histaxes  = True
    for lineset in linesetlist_fluxesFR:
        plot_lineratios_fromsummaryfiles_wrapper(plotbasename,fluxratiodat,lineset,histaxes,Nhistbins,cdatvec,cdattype,
                                                 Nsigma=Nsigma,point_text=point_text,vshiftmax=vshiftmax,
                                                 ylog=True, xlog=True,
                                                 literaturevaluesadded=addliteraturevalues,performlinearfit=True,
                                                 showlimits=showlimits,overwrite=overwrite,verbose=verbose,litsymboldot=litsymboldot)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    linesetlist_fluxes = []
    linesetlist_fluxes.append(['NV'   ,'s2n_NV'   ,None,None,fluxes_range, s2n_range,   None])
    linesetlist_fluxes.append(['CIV'  ,'s2n_CIV'  ,None,None,fluxes_range, s2n_range,   None])
    linesetlist_fluxes.append(['HeII' ,'s2n_HeII' ,None,None,fluxes_range, s2n_range,   None])
    linesetlist_fluxes.append(['OIII' ,'s2n_OIII' ,None,None,fluxes_range, s2n_range,   None])
    linesetlist_fluxes.append(['SiIII','s2n_SiIII',None,None,fluxes_range, s2n_range,   None])
    linesetlist_fluxes.append(['CIII' ,'s2n_CIII' ,None,None,fluxes_range, s2n_range,   None])
    linesetlist_fluxes.append(['MgII' ,'s2n_MgII' ,None,None,fluxes_range, s2n_range,   None])

    linesetlist_fluxes.append(['CIII','CIV'   ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['CIII','OIII'  ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['CIII','HeII'  ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['CIII','MgII'  ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['CIII','NV'    ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['CIII','SiIII' ,None,None,fluxes_range, fluxes_range,   None])

    linesetlist_fluxes.append(['CIV','OIII'   ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['CIV','HeII'   ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['CIV','MgII'   ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['CIV','NV'     ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['CIV','SiIII'  ,None,None,fluxes_range, fluxes_range,   None])
    #
    linesetlist_fluxes.append(['OIII','HeII'  ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['OIII','MgII'  ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['OIII','NV'    ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['OIII','SiIII' ,None,None,fluxes_range, fluxes_range,   None])
    #
    linesetlist_fluxes.append(['HeII','MgII'  ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['HeII','NV'    ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['HeII','SiIII' ,None,None,fluxes_range, fluxes_range,   None])
    #
    # No good values linesetlist_fluxes.append(['MgII','NV'    ,None,None,fluxes_range, fluxes_range,   None])
    linesetlist_fluxes.append(['MgII','SiIII' ,None,None,fluxes_range, fluxes_range,   None])
    #
    linesetlist_fluxes.append(['NV','SiIII'   ,None,None,fluxes_range, fluxes_range,   None])

    Nhistbins = 30
    histaxes  = True
    for lineset in linesetlist_fluxes:
        plot_lineratios_fromsummaryfiles_wrapper(plotbasename,fluxratiodat,lineset,histaxes,Nhistbins,cdatvec,cdattype,
                                                 Nsigma=Nsigma,point_text=point_text,vshiftmax=vshiftmax,
                                                 literaturevaluesadded=addliteraturevalues,performlinearfit=True,
                                                 showlimits=showlimits,overwrite=overwrite,verbose=verbose,litsymboldot=litsymboldot)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    linesetlist  = []
    linesetlist.append(['CIV','CIII','CIV','HeII',ratios_range,ratios_range  , 'Schmidt+17 fig. 7 top,   Feltre+16 fig A2a'])
    linesetlist.append(['CIII','HeII','CIV','HeII',ratios_range,ratios_range , 'Schmidt+17 fig. 7 center                  '])
    linesetlist.append(['CIV','OIII','CIV','HeII',ratios_range,ratios_range  , 'Schmidt+17 fig. 7 bottom                  '])
    linesetlist.append(['SiIII','CIII','OIII','CIII',ratios_range,ratios_range,'Byler+20 fig. 3                           '])
    linesetlist.append(['HeII','CIII','OIII','CIII',ratios_range,ratios_range, 'Byler+20 fig. 4                           '])
    linesetlist.append(['OIII','CIII','CIV','OIII',ratios_range,ratios_range,  'Byler+20 fig. 5                           '])
    linesetlist.append(['CIII','HeII','NV','HeII',ratios_range,ratios_range  ,'Plat+19 fig. 6d                           '])
    linesetlist.append(['CIII','OIII','CIV','CIII',ratios_range,ratios_range ,'Plat+19 fig. 6f                           '])
    linesetlist.append(['CIV','HeII','CIV','CIII',ratios_range,ratios_range  ,'Feltre+16 fig 5                           '])
    linesetlist.append(['CIV','HeII','CIII','HeII',ratios_range,ratios_range ,'Feltre+16 fig 6                           '])
    linesetlist.append(['NV','HeII','CIII','HeII',ratios_range,ratios_range  ,'Feltre+16 fig 8, fig A1b                  '])
    linesetlist.append(['CIV','CIII','CIII','HeII',ratios_range,ratios_range ,'Feltre+16 fig A1a                         '])
    linesetlist.append(['NV','CIV','CIII','HeII',ratios_range,ratios_range   ,'Feltre+16 fig A1c                         '])
    linesetlist.append(['NV','HeII','CIV','HeII',ratios_range,ratios_range   ,'Feltre+16 fig A2b                         '])
    linesetlist.append(['NV','CIV','CIV','HeII',ratios_range,ratios_range    ,'Feltre+16 fig A2c                         '])
    linesetlist.append(['OIII','HeII','CIV','HeII',ratios_range,ratios_range ,'Feltre+16 fig A2e                         '])
    linesetlist.append(['OIII','HeII','CIII','HeII',ratios_range,ratios_range ,'Hirschmann+19 fig 6                         '])
    linesetlist.append(['SiIII','HeII','CIII','HeII',ratios_range,ratios_range ,'Hirschmann+19 fig 6                         '])
    #
    linesetlist.append(['SiIII','HeII','CIV','HeII',ratios_range,ratios_range,'Feltre+16 fig A2i                         '])
    linesetlist.append(['CIII','OIII','HeII','CIII',ratios_range,ratios_range  , None])
    linesetlist.append(['CIII','CIV','OIII','HeII',ratios_range,ratios_range  , None])
    linesetlist.append(['CIII','CIV','OIII','SiIII',ratios_range,ratios_range  , None])
    linesetlist.append(['CIV','SiIII','OIII','HeII',ratios_range,ratios_range, None])

    linesetlist.append(['MgII','SiIII','OIII','HeII',ratios_range,ratios_range, None])

    Nhistbins = 30
    histaxes  = False
    for lineset in linesetlist:
        if 'MgII' in lineset:                # No MgII columns in the NEOGAL photoionisation models
            photoionizationplotparam = None
        # elif 'NV' in lineset:                # No columns with NV in Maseda line ratio summaryfile
        #     continue
        else:
            #varyparam, cutSFmodels, markersize, SFmarker, AGNmarker, linestrings, doubletratios = piplotparam # for uves. version of add-function
            # photoionizationplotparam = colorvar_pi, False, 1.5, 's', 'D', \
            #                            [uves.linenameUVES2NEOGAL(lineset[0]),uves.linenameUVES2NEOGAL(lineset[1]),
            #                             uves.linenameUVES2NEOGAL(lineset[2]),uves.linenameUVES2NEOGAL(lineset[3])], \
            #                            [None,None,None,None]

            # x2plot, y2plot, varyparam, cutSFmodels, markersize, SFmarker, AGNmarker = piplotparam  # for lce. version of add-function
            photoionizationplotparam = uves.linenameUVES2NEOGAL(lineset[0])+'/'+uves.linenameUVES2NEOGAL(lineset[1]),\
                                       uves.linenameUVES2NEOGAL(lineset[2])+'/'+uves.linenameUVES2NEOGAL(lineset[3]), \
                                       colorvar_pi, False, 1.5, 's', 'D'


        plot_lineratios_fromsummaryfiles_wrapper(plotbasename,fluxratiodat,lineset,histaxes,Nhistbins,cdatvec,cdattype,
                                                 Nsigma=Nsigma,point_text=point_text,vshiftmax=vshiftmax,
                                                 literaturevaluesadded=addliteraturevalues,performlinearfit=True,
                                                 overwrite=overwrite,verbose=verbose,showlimits=showlimits,
                                                 photoionizationplotparam=photoionizationplotparam,litsymboldot=litsymboldot)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_lineratios_fromsummaryfiles_wrapper(plotbasename,fluxratiodat,lineset,histaxes,Nhistbins,cdatvec,cdattype,
                                             photoionizationplotparam=None,point_text=None,showlimits=True,
                                             literaturevaluesadded=False,performlinearfit=False, ylog=True, xlog=True,
                                             overwrite=False,Nsigma=3.0,vshiftmax=1e4,verbose=True,litsymboldot=False):
    """
    Wrapper to define input data and excecute plot command

    """
    line1,line2,line3,line4,xrange,yrange,title = lineset

    linetype='onetoone'
    if 's2n' in line2:
        linetype='threesigma_y'

    if (line3 is None) & (line4 is None):
        if 's2n' in line2:
            plotname = plotbasename+'_linefluxes_'+line1+'vs'+line2+\
                       '_Nsigma'+str(Nsigma).replace('.','p')+\
                       '_vshiftLT'+str(vshiftmax).replace('.','p')+'.pdf'
            if ('f_'+line1 in fluxratiodat.dtype.names) & (line2 in fluxratiodat.dtype.names) & showlimits:
                goodent  = np.where(np.isfinite(fluxratiodat['f_'+line1]) & np.isfinite(fluxratiodat['f_'+line2]) &
                                    (np.abs(fluxratiodat['vshift_'+line1]) < vshiftmax) &
                                    (fluxratiodat['id'].astype(float) < 1e9))[0]
            elif ('f_'+line1 in fluxratiodat.dtype.names) & (line2 in fluxratiodat.dtype.names) & (not showlimits):
                goodent  = np.where((np.abs(fluxratiodat['ferr_'+line1]) != 99) &
                                    np.isfinite(fluxratiodat['f_'+line1]) & np.isfinite(fluxratiodat[line2]) &
                                    (np.abs(fluxratiodat['vshift_'+line1]) < vshiftmax) &
                                    (fluxratiodat['id'].astype(float) < 1e9))[0]
                plotname = plotname.replace('.pdf','_nolimits.pdf')
            else:
                goodent  = []
        elif ('1' in line1) & ('2' in line2):
            plotname = plotbasename+'_linefluxes_'+line1+'vs'+line2+\
                       '_Nsigma'+str(Nsigma).replace('.','p')+\
                       '_vshiftLT'+str(vshiftmax).replace('.','p')+'.pdf'

            if ('f_'+line1 in fluxratiodat.dtype.names) & ('f_'+line2 in fluxratiodat.dtype.names) & showlimits:
                goodent  = np.where(np.isfinite(fluxratiodat['f_'+line1]) & np.isfinite(fluxratiodat['f_'+line2]) &
                                    (np.abs(fluxratiodat['vshift_'+line1[:-1]]) < vshiftmax) &
                                    (np.abs(fluxratiodat['vshift_'+line2[:-1]]) < vshiftmax) &
                                    (fluxratiodat['id'].astype(float) < 1e9))[0]
            elif ('f_'+line1 in fluxratiodat.dtype.names) & ('f_'+line2 in fluxratiodat.dtype.names) & (not showlimits):
                goodent  = np.where((np.abs(fluxratiodat['ferr_'+line1]) != 99) & (np.abs(fluxratiodat['ferr_'+line2]) != 99) &
                                    np.isfinite(fluxratiodat['f_'+line1]) & np.isfinite(fluxratiodat['f_'+line2]) &
                                    (np.abs(fluxratiodat['vshift_'+line1[:-1]]) < vshiftmax) &
                                    (np.abs(fluxratiodat['vshift_'+line2[:-1]]) < vshiftmax) &
                                    (fluxratiodat['id'].astype(float) < 1e9))[0]
                plotname = plotname.replace('.pdf','_nolimits.pdf')
            else:
                goodent  = []

        else:
            plotname = plotbasename+'_linefluxes_'+line1+'vs'+line2+\
                       '_Nsigma'+str(Nsigma).replace('.','p')+\
                       '_vshiftLT'+str(vshiftmax).replace('.','p')+'.pdf'

            if ('f_'+line1 in fluxratiodat.dtype.names) & ('f_'+line2 in fluxratiodat.dtype.names) & showlimits:
                goodent  = np.where(np.isfinite(fluxratiodat['f_'+line1]) & np.isfinite(fluxratiodat['f_'+line2]) &
                                    (np.abs(fluxratiodat['vshift_'+line1]) < vshiftmax) &
                                    (np.abs(fluxratiodat['vshift_'+line2]) < vshiftmax) &
                                    (fluxratiodat['id'].astype(float) < 1e9))[0]
            elif ('f_'+line1 in fluxratiodat.dtype.names) & ('f_'+line2 in fluxratiodat.dtype.names) & (not showlimits):
                goodent  = np.where((np.abs(fluxratiodat['ferr_'+line1]) != 99) & (np.abs(fluxratiodat['ferr_'+line2]) != 99) &
                                    np.isfinite(fluxratiodat['f_'+line1]) & np.isfinite(fluxratiodat['f_'+line2]) &
                                    (np.abs(fluxratiodat['vshift_'+line1]) < vshiftmax) &
                                    (np.abs(fluxratiodat['vshift_'+line2]) < vshiftmax) &
                                    (fluxratiodat['id'].astype(float) < 1e9))[0]
                plotname = plotname.replace('.pdf','_nolimits.pdf')
            else:
                goodent  = []

        if 's2n' in line1:
            xlabel   = 'S/N('+line1.split('_')[1]+')'
        else:
            xlabel   = line1+' [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]'

        if 's2n' in line2:
            ylabel   = 'S/N('+line2.split('_')[1]+')'
        else:
            ylabel   = line2+' [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]'

        if len(goodent) == 0:
            pass
            # if verbose: print('\n - OLD WARNING No good values found for the plot: \n           '+plotname.split('/')[-1]+'\n')
            # goodent  = np.asarray([0,1])
            # xvalues  = [1e10]*2
            # xerr     = [1.0]*2
            # yvalues  = [1e10]*2
            # yerr     = [1.0]*2
            # cdatvec  = np.asarray([0.0]*2)
        else:
            if 's2n' in line1:
                xvalues  = fluxratiodat[line1][goodent]
                xerr     = fluxratiodat[line1][goodent]*np.nan
            elif ('1' in line1) & ('2' in line2):
                xvalues  = fluxratiodat['f_'+line1][goodent]
                xerr     = fluxratiodat['ferr_'+line1][goodent]

                xlimits_ent  = np.where(fluxratiodat['s2n_'+line1[:-1]][goodent] < Nsigma)[0]
                if len(xlimits_ent) > 0:
                    xvalues[xlimits_ent] = xerr[xlimits_ent] * Nsigma
                    xerr[xlimits_ent]    = +99 # upper limit
            else:
                xvalues  = fluxratiodat['f_'+line1][goodent]
                xerr     = fluxratiodat['ferr_'+line1][goodent]

                xlimits_ent  = np.where(fluxratiodat['s2n_'+line1][goodent] < Nsigma)[0]
                if len(xlimits_ent) > 0:
                    xvalues[xlimits_ent] = xerr[xlimits_ent] * Nsigma
                    xerr[xlimits_ent]    = +99 # upper limit


            if 's2n' in line2:
                yvalues  = fluxratiodat[line2][goodent]
                yerr     = fluxratiodat[line2][goodent]*np.nan
            elif ('1' in line1) & ('2' in line2):
                yvalues  = fluxratiodat['f_'+line2][goodent]
                yerr     = fluxratiodat['ferr_'+line2][goodent]

                ylimits_ent  = np.where(fluxratiodat['s2n_'+line2[:-1]][goodent] < Nsigma)[0]
                if len(ylimits_ent) > 0:
                    yvalues[ylimits_ent] = yerr[ylimits_ent] * Nsigma
                    yerr[ylimits_ent]    = +99 # upper limit
            else:
                yvalues  = fluxratiodat['f_'+line2][goodent]
                yerr     = fluxratiodat['ferr_'+line2][goodent]

                ylimits_ent  = np.where(fluxratiodat['s2n_'+line2][goodent] < Nsigma)[0]
                if len(ylimits_ent) > 0:
                    yvalues[ylimits_ent] = yerr[ylimits_ent] * Nsigma
                    yerr[ylimits_ent]    = +99 # upper limit
    else:
        plotname = plotbasename+'_lineratios_'+line1+line2+'vs'+line3+line4+\
                   '_Nsigma'+str(Nsigma).replace('.','p')+\
                   '_vshiftLT'+str(vshiftmax).replace('.','p')+'.pdf'

        if ('FR_'+line1+line2 in fluxratiodat.dtype.names) & ('FR_'+line3+line4 in fluxratiodat.dtype.names) & showlimits:
            goodent  = np.where(np.isfinite(fluxratiodat['FR_'+line1+line2]) &
                                np.isfinite(fluxratiodat['FR_'+line3+line4]) &
                                (np.abs(fluxratiodat['vshift_'+line1]) < vshiftmax) &
                                (np.abs(fluxratiodat['vshift_'+line2]) < vshiftmax) &
                                (np.abs(fluxratiodat['vshift_'+line3]) < vshiftmax) &
                                (np.abs(fluxratiodat['vshift_'+line4]) < vshiftmax) &
                                (fluxratiodat['id'].astype(float) < 1e9))[0]
        elif ('FR_'+line1+line2 in fluxratiodat.dtype.names) & ('FR_'+line3+line4 in fluxratiodat.dtype.names) & (not showlimits):
            goodent  = np.where(np.isfinite(fluxratiodat['FR_'+line1+line2]) &
                                np.isfinite(fluxratiodat['FR_'+line3+line4]) &
                                (np.abs(fluxratiodat['FRerr_'+line1+line2]) != 99) &
                                (np.abs(fluxratiodat['FRerr_'+line3+line4]) != 99) &
                                (np.abs(fluxratiodat['vshift_'+line1]) < vshiftmax) &
                                (np.abs(fluxratiodat['vshift_'+line2]) < vshiftmax) &
                                (np.abs(fluxratiodat['vshift_'+line3]) < vshiftmax) &
                                (np.abs(fluxratiodat['vshift_'+line4]) < vshiftmax) &
                                (fluxratiodat['id'].astype(float) < 1e9))[0]
            plotname = plotname.replace('.pdf','_nolimits.pdf')
        else:
            goodent  = []

        xlabel   = line1+'/'+line2
        ylabel   = line3+'/'+line4

        if len(goodent) == 0:
            pass
            # if verbose: print('\n - OLD WARNING No good values found for the plot \n           '+plotname.split('/')[-1]+'\n')
            # goodent  = np.asarray([0,1])
            # xvalues  = [1e10]*2
            # xerr     = [1.0]*2
            # yvalues  = [1e10]*2
            # yerr     = [1.0]*2
            # cdatvec  = np.asarray([0.0]*2)
        else:
            xvalues  = fluxratiodat['FR_'+line1+line2][goodent]
            xerr     = fluxratiodat['FRerr_'+line1+line2][goodent]
            yvalues  = fluxratiodat['FR_'+line3+line4][goodent]
            yerr     = fluxratiodat['FRerr_'+line3+line4][goodent]

            xlimits_ent_num  = np.where((fluxratiodat['s2n_'+line1][goodent] < Nsigma) &
                                        (fluxratiodat['s2n_'+line1][goodent] > 0.0))[0]
            xlimits_ent_den  = np.where((fluxratiodat['s2n_'+line2][goodent] < Nsigma) &
                                        (fluxratiodat['s2n_'+line2][goodent] > 0.0))[0]

            for xent, xval in enumerate(xvalues):
                if (xent in xlimits_ent_num) & (xent in xlimits_ent_den):
                    xvalues[xent]     = np.nan
                    xerr[xent] = 0      # no constraint
                elif (xent in xlimits_ent_num) & (xent not in xlimits_ent_den):
                    xvalues[xent] = xerr[xent] * Nsigma
                    xerr[xent]    = +99 # upper limit
                elif (xent not in xlimits_ent_num) & (xent in xlimits_ent_den):
                    xvalues[xent] = xerr[xent] * Nsigma
                    xerr[xent]    = -99 # lower limit

            ylimits_ent_num  = np.where((fluxratiodat['s2n_'+line3][goodent] < Nsigma) &
                                        (fluxratiodat['s2n_'+line3][goodent] > 0.0))[0]
            ylimits_ent_den  = np.where((fluxratiodat['s2n_'+line4][goodent] < Nsigma) &
                                        (fluxratiodat['s2n_'+line4][goodent] > 0.0))[0]
            for yent, yval in enumerate(yvalues):
                if (yent in ylimits_ent_num) & (yent in ylimits_ent_den):
                    yvalues[yent]   = np.nan
                    yerr[yent]      = 0      # no constraint
                elif (yent in ylimits_ent_num) & (yent not in ylimits_ent_den):
                    yvalues[yent] = yerr[yent] * Nsigma
                    yerr[yent]    = +99 # upper limit
                elif (yent not in ylimits_ent_num) & (yent in ylimits_ent_den):
                    yvalues[yent] = yerr[yent] * Nsigma
                    yerr[yent]    = -99 # lower limit

            if (xerr == 0).all():
                if verbose: print('\n - WARNING all values for '+line1+' and '+line2+' were below '+str(Nsigma)+
                                  'sigma for the plot \n           '+plotname.split('/')[-1]+'\n')
                #return
            if (yerr == 0).all():
                if verbose: print('\n - WARNING all values for '+line3+' and '+line4+' were below '+str(Nsigma)+
                                  'sigma for the plot \n           '+plotname.split('/')[-1]+'\n')
                #return

    if point_text is not None:
        point_textALL = point_text[goodent]
    else:
        point_textALL = point_text

    IDsALL     = fluxratiodat['id'][goodent]
    cdatvecALL = cdatvec[goodent]

    if len(cdatvecALL) == 0: # no objects found in FELIS output
        xvalues = np.array([])
        xerr    = np.array([])
        yvalues = np.array([])
        yerr    = np.array([])

    # - - - - - - - - - Literature - - - - - - - - -
    if (line3 is None) & (line4 is None):
        if 's2n' in line2:
            if ('f_'+line1 in fluxratiodat.dtype.names) & (line2 in fluxratiodat.dtype.names) & showlimits:
                litent  = np.where(np.isfinite(fluxratiodat['f_'+line1]) &
                                   np.isfinite(fluxratiodat[line2]) &
                                   (fluxratiodat['id'].astype(float) > 1e9))[0]

            elif ('f_'+line1 in fluxratiodat.dtype.names) & (line2 in fluxratiodat.dtype.names) & (not showlimits):
                litent  = np.where(np.isfinite(fluxratiodat['f_'+line1]) &
                                   np.isfinite(fluxratiodat[line2]) &
                                   np.abs((fluxratiodat['ferr_'+line1]) != 99) &
                                   (fluxratiodat['id'].astype(float) > 1e9))[0]
            else:
                litent = []

            if len(litent) > 0:
                xvalues    = np.append(xvalues,fluxratiodat['f_'+line1][litent])
                xerr       = np.append(xerr,fluxratiodat['ferr_'+line1][litent])
                yvalues    = np.append(yvalues,fluxratiodat[line2][litent])
                yerr       = np.append(yerr,fluxratiodat[line2][litent]*np.nan)
                cdatvecALL = np.append(cdatvecALL,cdatvec[litent])
                if litsymboldot:
                    IDsALL     = np.append(IDsALL,fluxratiodat['id'][litent]*0.0+990000000000)
                else:
                    IDsALL     = np.append(IDsALL,fluxratiodat['id'][litent])
                if point_text is not None:
                    point_textALL = np.append(point_textALL,point_text[litent])
        else:
            if ('f_'+line1 in fluxratiodat.dtype.names) & ('f_'+line2 in fluxratiodat.dtype.names) & showlimits:
                litent  = np.where(np.isfinite(fluxratiodat['f_'+line1]) &
                                   np.isfinite(fluxratiodat['f_'+line2]) &
                                   (fluxratiodat['id'].astype(float) > 1e9))[0]

            elif ('f_'+line1 in fluxratiodat.dtype.names) & ('f_'+line2 in fluxratiodat.dtype.names) & (not showlimits):
                litent  = np.where(np.isfinite(fluxratiodat['f_'+line1]) &
                                   np.isfinite(fluxratiodat['f_'+line2]) &
                                   np.abs((fluxratiodat['ferr_'+line1]) != 99) &
                                   np.abs((fluxratiodat['ferr_'+line2]) != 99) &
                                   (fluxratiodat['id'].astype(float) > 1e9))[0]
            else:
                litent = []

            if len(litent) > 0:
                xvalues    = np.append(xvalues,fluxratiodat['f_'+line1][litent])
                xerr       = np.append(xerr,fluxratiodat['ferr_'+line1][litent])
                yvalues    = np.append(yvalues,fluxratiodat['f_'+line2][litent])
                yerr       = np.append(yerr,fluxratiodat['ferr_'+line2][litent])
                cdatvecALL = np.append(cdatvecALL,cdatvec[litent])
                if litsymboldot:
                    IDsALL     = np.append(IDsALL,fluxratiodat['id'][litent]*0.0+990000000000)
                else:
                    IDsALL     = np.append(IDsALL,fluxratiodat['id'][litent])
                if point_text is not None:
                    point_textALL = np.append(point_textALL,point_text[litent])
    else:
        if ('FR_'+line1+line2 in fluxratiodat.dtype.names) & ('FR_'+line3+line4 in fluxratiodat.dtype.names) & showlimits:
            litent  = np.where(np.isfinite(fluxratiodat['FR_'+line1+line2]) &
                               np.isfinite(fluxratiodat['FR_'+line3+line4]) &
                               (fluxratiodat['id'].astype(float) > 1e9))[0]
        elif ('FR_'+line1+line2 in fluxratiodat.dtype.names) & ('FR_'+line3+line4 in fluxratiodat.dtype.names) & (not showlimits):
            litent  = np.where(np.isfinite(fluxratiodat['FR_'+line1+line2]) &
                               np.isfinite(fluxratiodat['FR_'+line3+line4]) &
                               np.abs((fluxratiodat['FRerr_'+line1+line2]) != 99) &
                               np.abs((fluxratiodat['FRerr_'+line3+line4]) != 99) &
                               (fluxratiodat['id'].astype(float) > 1e9))[0]
        else:
            litent = []

        if len(litent) > 0:
            xvalues    = np.append(xvalues,fluxratiodat['FR_'+line1+line2][litent])
            xerr       = np.append(xerr,fluxratiodat['FRerr_'+line1+line2][litent])
            yvalues    = np.append(yvalues,fluxratiodat['FR_'+line3+line4][litent])
            yerr       = np.append(yerr,fluxratiodat['FRerr_'+line3+line4][litent])
            cdatvecALL = np.append(cdatvecALL,cdatvec[litent])
            if litsymboldot:
                IDsALL     = np.append(IDsALL,fluxratiodat['id'][litent]*0.0+990000000000)
            else:
                IDsALL     = np.append(IDsALL,fluxratiodat['id'][litent])
            if point_text is not None:
                point_textALL = np.append(point_textALL,point_text[litent])

    # - - - - - - - - - - - - - - - - - - - - - - - -

    if len(cdatvecALL) == 0:
        if verbose: print('\n - WARNING No good values found for the plot: \n           '+plotname.split('/')[-1]+'\n')
        xvalues    = [1e10]*2
        xerr       = [1.0]*2
        yvalues    = [1e10]*2
        yerr       = [1.0]*2
        cdatvecALL = np.asarray([0.0]*2)
        IDsALL     = np.asarray([0.0]*2)

    if literaturevaluesadded:
        plotname = plotname.replace('.pdf','_wLit.pdf')

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,
                                                   'dummydat',linetype=linetype,title=title,ids=IDsALL,
                                                   ylog=ylog,xlog=xlog,yrange=yrange,xrange=xrange,
                                                   colortype=cdattype,colorcode=True,cdatvec=cdatvecALL,
                                                   point_text=point_textALL,photoionizationplotparam=photoionizationplotparam,
                                                   histaxes=histaxes,Nbins=Nhistbins,
                                                   overwrite=overwrite,verbose=verbose)

    if performlinearfit & (len(cdatvecALL) > 0) & (len(xvalues[np.isfinite(xvalues)]) > 2) & (len(yvalues[np.isfinite(yvalues)]) > 2):
        plotname   = plotname.replace('.pdf','_ODRfit2data.pdf')

        if ylog:
            yval_fit = np.log10(yvalues)
            yerr_fit = 0.434 * yerr/yvalues  # https://faculty.washington.edu/stuve/log_error.pdf
        else:
            yval_fit = yvalues
            yerr_fit = yerr

        if xlog:
            xval_fit = np.log10(xvalues)
            xerr_fit = 0.434 * xerr/xvalues  # https://faculty.washington.edu/stuve/log_error.pdf
        else:
            xval_fit = xvalues
            xerr_fit = xerr

        if line1[:-1] == line2[:-1]:
            if ylog & xlog:
                ratios = 10**xval_fit/10**yval_fit
            else:
                ratios = xval_fit/yval_fit
            meanval   = np.mean(ratios[np.isfinite(ratios)])
            medianval = np.median(ratios[np.isfinite(ratios)])
            stdval    = np.std(ratios[np.isfinite(ratios)])
            if verbose: print('   ----> Ratio '+line1+'/'+line2+' (mean,median,std) = '+
                              str("%.2f" % meanval)+','+str("%.2f" % medianval)+','+str("%.2f" % stdval))
        fitres,rp,pp,rs,ps = kbs.fit_function_to_data_with_errors_on_both_axes(xval_fit,yval_fit,xerr_fit,yerr_fit,
                                                                               fitfunction='linear',plotresults=plotname,
                                                                               returnCorrelationCoeffs=True)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_lineratios_fromsummaryfiles_vsInfofile(plotbasename,fluxratiodat,lineset,histaxes,Nhistbins,cdatvec,cdattype,
                                                addliteraturevalues=True,performlinearfit=False, litsymboldot=False,
                                                point_text=None,showlimits=True,ylog=False,xlog=False,
                                                overwrite=False,Nsigma=3.0,vshiftmax=1e4,verbose=True):
    """
    Wrapper to define input data and excecute plot command

    """
    # infofile         = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile         = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    infofiledat      = afits.open(infofile)[1].data
    infofiledat      = infofiledat[np.where((infofiledat['id']<4.9e8) | (infofiledat['id']>5.9e8))[0]] # ignoring UDF MW mock ids

    line1,line2,line3,line4,xrange,yrange,title = lineset
    if line2 is not None:
        sys.exit(' Expects line2 to be None for uves.plot_lineratios_fromsummaryfiles_vsInfofile() linelist input; it was "'+line2+'"')

    infocols = uves.get_infodat_plotcols()

    xlabel    = infocols[line1][2]

    if (line4 is None):
        plotname = plotbasename+'_linefluxes_'+line1+'vs'+line3+\
                   '_Nsigma'+str(Nsigma).replace('.','p')+\
                   '_vshiftLT'+str(vshiftmax).replace('.','p')+'.pdf'

        if 'vshift_'+line3 in fluxratiodat.dtype.names:
            vshiftcol3 = 'vshift_'+line3
        elif 'vshift_'+line3.split('R_')[1].split('1')[0] in fluxratiodat.dtype.names:
            vshiftcol3 = 'vshift_'+line3.split('R_')[1].split('1')[0]
        else:
            vshiftcol3 = 'vshift_'+line3[:-1]

        if ('f_'+line3 in fluxratiodat.dtype.names) & showlimits:
            goodent  = np.where(np.isfinite(fluxratiodat['f_'+line3]) &
                                (np.abs(fluxratiodat[vshiftcol3]) < vshiftmax) &
                                (fluxratiodat['id'].astype(float) < 1e9))[0]
        elif ('f_'+line3 in fluxratiodat.dtype.names) & (not showlimits):
            goodent  = np.where((np.abs(fluxratiodat['ferr_'+line3]) != 99) &
                                np.isfinite(fluxratiodat['f_'+line3]) &
                                (np.abs(fluxratiodat[vshiftcol3]) < vshiftmax) &
                                (fluxratiodat['id'].astype(float) < 1e9))[0]
            plotname = plotname.replace('.pdf','_nolimits.pdf')
        elif (line3 in fluxratiodat.dtype.names) & showlimits:
            goodent  = np.where(np.isfinite(fluxratiodat[line3]) &
                                (np.abs(fluxratiodat[vshiftcol3]) < vshiftmax) &
                                (fluxratiodat['id'].astype(float) < 1e9))[0]
        elif (line3 in fluxratiodat.dtype.names) & (not showlimits):
            goodent  = np.where((np.abs(fluxratiodat[line3.replace('_','err_')]) != 99) &
                                np.isfinite(fluxratiodat[line3]) &
                                (np.abs(fluxratiodat[vshiftcol3]) < vshiftmax) &
                                (fluxratiodat['id'].astype(float) < 1e9))[0]
            plotname = plotname.replace('.pdf','_nolimits.pdf')
        else:
            goodent  = []

        if 'FR_' in line3:
            linename = line3.split('R_')[1].split('1')[0]
            ylabel   = linename+' doublet flux ratio'
        else:
            ylabel   = line3+' [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]'

        if len(goodent) == 0:
            pass
        else:
            if 'FR_' in line3:
                yvalues  = fluxratiodat[line3][goodent]
                yerr     = fluxratiodat[line3.replace('_','err_')][goodent]

            else:
                yvalues  = fluxratiodat['f_'+line3][goodent]
                yerr     = fluxratiodat['ferr_'+line3][goodent]

                try:
                    ylimits_ent  = np.where(fluxratiodat['s2n_'+line3][goodent] < Nsigma)[0]
                except:
                    ylimits_ent  = np.where(fluxratiodat['f_'+line3][goodent]/fluxratiodat['ferr_'+line3][goodent] < Nsigma)[0]

                if len(ylimits_ent) > 0:
                    yvalues[ylimits_ent] = yerr[ylimits_ent] * Nsigma
                    yerr[ylimits_ent]    = +99 # upper limit
    else:
        plotname = plotbasename+'_lineratios_'+line1+'vs'+line3+line4+\
                   '_Nsigma'+str(Nsigma).replace('.','p')+\
                   '_vshiftLT'+str(vshiftmax).replace('.','p')+'.pdf'

        if 'vshift_'+line3 in fluxratiodat.dtype.names:
            vshiftcol3 = 'vshift_'+line3
        else:
            vshiftcol3 = 'vshift_'+line3[:-1]
        if 'vshift_'+line4 in fluxratiodat.dtype.names:
            vshiftcol4 = 'vshift_'+line4
        else:
            vshiftcol4 = 'vshift_'+line4[:-1]


        if ('FR_'+line3+line4 in fluxratiodat.dtype.names) & showlimits:
            goodent  = np.where(np.isfinite(fluxratiodat['FR_'+line3+line4]) &
                                (np.abs(fluxratiodat[vshiftcol3]) < vshiftmax) &
                                (np.abs(fluxratiodat[vshiftcol4]) < vshiftmax) &
                                (fluxratiodat['id'].astype(float) < 1e9))[0]
        elif ('FR_'+line3+line4 in fluxratiodat.dtype.names) & (not showlimits):
            goodent  = np.where(np.isfinite(fluxratiodat['FR_'+line3+line4]) &
                                (np.abs(fluxratiodat['FRerr_'+line3+line4]) != 99) &
                                (np.abs(fluxratiodat[vshiftcol3]) < vshiftmax) &
                                (np.abs(fluxratiodat[vshiftcol4]) < vshiftmax) &
                                (fluxratiodat['id'].astype(float) < 1e9))[0]
            plotname = plotname.replace('.pdf','_nolimits.pdf')
        else:
            goodent  = []

        ylabel   = line3+'/'+line4

        if len(goodent) == 0:
            pass
        else:
            yvalues  = fluxratiodat['FR_'+line3+line4][goodent]
            yerr     = fluxratiodat['FRerr_'+line3+line4][goodent]

            try:
                ylimits_ent_num  = np.where((fluxratiodat['s2n_'+line3][goodent] < Nsigma) &
                                            (fluxratiodat['s2n_'+line3][goodent] > 0.0))[0]
                ylimits_ent_den  = np.where((fluxratiodat['s2n_'+line4][goodent] < Nsigma) &
                                            (fluxratiodat['s2n_'+line4][goodent] > 0.0))[0]
            except:
                ylimits_ent_num  = np.where((fluxratiodat['f_'+line3][goodent]/fluxratiodat['ferr_'+line3][goodent] < Nsigma) &
                                            (fluxratiodat['f_'+line3][goodent]/fluxratiodat['ferr_'+line3][goodent] > 0.0))[0]
                ylimits_ent_den  = np.where((fluxratiodat['f_'+line4][goodent]/fluxratiodat['ferr_'+line4][goodent] < Nsigma) &
                                            (fluxratiodat['f_'+line4][goodent]/fluxratiodat['ferr_'+line4][goodent] > 0.0))[0]

            for yent, yval in enumerate(yvalues):
                if (yent in ylimits_ent_num) & (yent in ylimits_ent_den):
                    yvalues[yent]   = np.nan
                    yerr[yent]      = 0      # no constraint
                elif (yent in ylimits_ent_num) & (yent not in ylimits_ent_den):
                    yvalues[yent] = yerr[yent] * Nsigma
                    yerr[yent]    = +99 # upper limit
                elif (yent not in ylimits_ent_num) & (yent in ylimits_ent_den):
                    yvalues[yent] = yerr[yent] * Nsigma
                    yerr[yent]    = -99 # lower limit

            if (yerr == 0).all():
                if verbose: print('\n - WARNING all values for '+line3+' and '+line4+' were below '+str(Nsigma)+
                                  'sigma for the plot \n           '+plotname.split('/')[-1]+'\n')
                #return

    lyaents   = []
    for objid in fluxratiodat['id'][goodent]:
        lyaent = np.where(infofiledat['id'] == int(objid))[0]
        if len(lyaent) == 0:
            pdb.set_trace()
            sys.exit(' Weird... no match found in uves.plot_lineratios_fromsummaryfiles_vsInfofile()!')
        else:
            lyaents.append(lyaent[0])

    if len(goodent) > 0:
        xvalues  = infofiledat[infocols[line1][0]][np.asarray(lyaents)]
        if infocols[line1][1] is None:
            xerr     = np.asarray([np.nan]*len(xvalues))
        else:
            xerr     = infofiledat[infocols[line1][1]][np.asarray(lyaents)]

        if len(xvalues[xvalues == 0.0]) > 0:
            xvalues[xvalues == 0.0] = np.nan
            xerr[xvalues == 0.0]    = np.nan

    if point_text is not None:
        point_textALL = point_text[goodent]
    else:
        point_textALL = point_text

    IDsALL     = fluxratiodat['id'][goodent]
    cdatvecALL = cdatvec[goodent]

    if len(cdatvecALL) == 0: # no objects found in FELIS output
        xvalues = np.array([])
        xerr    = np.array([])
        yvalues = np.array([])
        yerr    = np.array([])

    # - - - - - - - - - Literature - - - - - - - - -
    if addliteraturevalues & ('FR_' not in line3):
        litfile = '/Users/kschmidt/work/catalogs/' \
                  'literaturecollection_emissionlinestrengths/literaturecollection_emissionlinestrengths.fits'
        fluxratiodat_lit = afits.open(litfile)[1].data

        coltranslationdic= {'f(CIII) [1e-20 erg/s/cm$^2$]':'f_CIII', 'zmanual':'redshift', 'ew_0':'EW0_Lya',
                            'redshift':'redshift'}
        cdatvec_lit = fluxratiodat_lit[coltranslationdic[cdattype]]

        if ('lyaew' in line1):
            if (line4 is None):
                if ('EW0_Lya' in fluxratiodat_lit.dtype.names) & showlimits:
                    litent  = np.where(np.isfinite(fluxratiodat_lit['EW0_Lya']) &
                                       np.isfinite(fluxratiodat_lit['f_'+line3]) &
                                       (fluxratiodat_lit['id'].astype(float) > 1e9))[0]

                elif ('EW0_Lya' in fluxratiodat_lit.dtype.names) & (not showlimits):
                    litent  = np.where(np.isfinite(fluxratiodat_lit['EW0_Lya']) &
                                       np.isfinite(fluxratiodat_lit['f_'+line3]) &
                                       np.abs((fluxratiodat_lit['EW0err_Lya']) != 99) &
                                       np.abs((fluxratiodat_lit['ferr_'+line3]) != 99) &
                                       (fluxratiodat_lit['id'].astype(float) > 1e9))[0]
                else:
                    litent = []

                if len(litent) > 0:
                    xvalues    = np.append(xvalues,fluxratiodat_lit['EW0_Lya'][litent])
                    xerr       = np.append(xerr,fluxratiodat_lit['EW0err_Lya'][litent])
                    yvalues    = np.append(yvalues,fluxratiodat_lit['f_'+line3][litent])
                    yerr       = np.append(yerr,fluxratiodat_lit['ferr_'+line3][litent])
                    cdatvecALL = np.append(cdatvecALL,cdatvec_lit[litent])
                    if litsymboldot:
                        IDsALL     = np.append(IDsALL,fluxratiodat_lit['id'][litent]*0.0+990000000000)
                    else:
                        IDsALL     = np.append(IDsALL,fluxratiodat_lit['id'][litent])
                    if point_text is not None:
                        point_textALL = np.append(point_textALL,point_text[litent])
            else:
                if ('EW0_Lya' in fluxratiodat_lit.dtype.names)  & showlimits:
                    litent  = np.where(np.isfinite(fluxratiodat_lit['EW0_Lya']) &
                                       np.isfinite(fluxratiodat_lit['FR_'+line3+line4]) &
                                       (fluxratiodat_lit['id'].astype(float) > 1e9))[0]
                elif ('EW0_Lya' in fluxratiodat_lit.dtype.names) & ('FR_'+line3+line4 in fluxratiodat_lit.dtype.names) & (not showlimits):
                    litent  = np.where(np.isfinite(fluxratiodat_lit['EW0_Lya']) &
                                       np.isfinite(fluxratiodat_lit['FR_'+line3+line4]) &
                                       np.abs((fluxratiodat_lit['EW0err_Lya']) != 99) &
                                       np.abs((fluxratiodat_lit['FRerr_'+line3+line4]) != 99) &
                                       (fluxratiodat_lit['id'].astype(float) > 1e9))[0]
                else:
                    litent = []

                if len(litent) > 0:
                    xvalues    = np.append(xvalues,fluxratiodat_lit['EW0_Lya'][litent])
                    xerr       = np.append(xerr,fluxratiodat_lit['EW0err_Lya'][litent])
                    yvalues    = np.append(yvalues,fluxratiodat_lit['FR_'+line3+line4][litent])
                    yerr       = np.append(yerr,fluxratiodat_lit['FRerr_'+line3+line4][litent])
                    cdatvecALL = np.append(cdatvecALL,cdatvec_lit[litent])
                    if litsymboldot:
                        IDsALL     = np.append(IDsALL,fluxratiodat_lit['id'][litent]*0.0+990000000000)
                    else:
                        IDsALL     = np.append(IDsALL,fluxratiodat_lit['id'][litent])
                    if point_text is not None:
                        point_textALL = np.append(point_textALL,point_text[litent])
        if ('redshift' in line1):
            if (line4 is None):
                if ('redshift' in fluxratiodat_lit.dtype.names) & showlimits:
                    litent  = np.where(np.isfinite(fluxratiodat_lit['redshift']) &
                                       np.isfinite(fluxratiodat_lit['f_'+line3]) &
                                       (fluxratiodat_lit['id'].astype(float) > 1e9))[0]

                elif ('redshift' in fluxratiodat_lit.dtype.names) & (not showlimits):
                    litent  = np.where(np.isfinite(fluxratiodat_lit['redshift']) &
                                       np.isfinite(fluxratiodat_lit['f_'+line3]) &
                                       np.abs((fluxratiodat_lit['EW0err_Lya']) != 99) &
                                       np.abs((fluxratiodat_lit['ferr_'+line3]) != 99) &
                                       (fluxratiodat_lit['id'].astype(float) > 1e9))[0]
                else:
                    litent = []

                if len(litent) > 0:
                    xvalues    = np.append(xvalues,fluxratiodat_lit['redshift'][litent])
                    xerr       = np.append(xerr,fluxratiodat_lit['redshift'][litent]*0.0)
                    yvalues    = np.append(yvalues,fluxratiodat_lit['f_'+line3][litent])
                    yerr       = np.append(yerr,fluxratiodat_lit['ferr_'+line3][litent])
                    cdatvecALL = np.append(cdatvecALL,cdatvec_lit[litent])
                    if litsymboldot:
                        IDsALL     = np.append(IDsALL,fluxratiodat_lit['id'][litent]*0.0+990000000000)
                    else:
                        IDsALL     = np.append(IDsALL,fluxratiodat_lit['id'][litent])
                    if point_text is not None:
                        point_textALL = np.append(point_textALL,point_text[litent])
            else:
                if ('redshift' in fluxratiodat_lit.dtype.names)  & showlimits:
                    litent  = np.where(np.isfinite(fluxratiodat_lit['redshift']) &
                                       np.isfinite(fluxratiodat_lit['FR_'+line3+line4]) &
                                       (fluxratiodat_lit['id'].astype(float) > 1e9))[0]
                elif ('EW0_Lya' in fluxratiodat_lit.dtype.names) & ('FR_'+line3+line4 in fluxratiodat_lit.dtype.names) & (not showlimits):
                    litent  = np.where(np.isfinite(fluxratiodat_lit['redshift']) &
                                       np.isfinite(fluxratiodat_lit['FR_'+line3+line4]) &
                                       np.abs((fluxratiodat_lit['FRerr_'+line3+line4]) != 99) &
                                       (fluxratiodat_lit['id'].astype(float) > 1e9))[0]
                else:
                    litent = []

                if len(litent) > 0:
                    xvalues    = np.append(xvalues,fluxratiodat_lit['redshift'][litent])
                    xerr       = np.append(xerr,fluxratiodat_lit['redshift'][litent]*0.0)
                    yvalues    = np.append(yvalues,fluxratiodat_lit['FR_'+line3+line4][litent])
                    yerr       = np.append(yerr,fluxratiodat_lit['FRerr_'+line3+line4][litent])
                    cdatvecALL = np.append(cdatvecALL,cdatvec_lit[litent])
                    if litsymboldot:
                        IDsALL     = np.append(IDsALL,fluxratiodat_lit['id'][litent]*0.0+990000000000)
                    else:
                        IDsALL     = np.append(IDsALL,fluxratiodat_lit['id'][litent])
                    if point_text is not None:
                        point_textALL = np.append(point_textALL,point_text[litent])

    # - - - - - - - - - - - - - - - - - - - - - - - -

    if len(cdatvecALL) == 0:
        if verbose: print('\n - WARNING No good values found for the plot: \n           '+plotname.split('/')[-1]+'\n')
        xvalues    = np.asarray([1e10]*2)
        xerr       = np.asarray([1.0]*2)
        yvalues    = np.asarray([1e10]*2)
        yerr       = np.asarray([1.0]*2)
        cdatvecALL = np.asarray([0.0]*2)
        IDsALL     = np.asarray([0.0]*2)

    if 'FR_' in line3:
        linetype=None
    else:
        linetype='onetoone'

    if addliteraturevalues:
        plotname = plotname.replace('.pdf','_wLit.pdf')
    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,
                                                   'dummydat',linetype=linetype,title=title,ids=IDsALL,
                                                   ylog=ylog,xlog=xlog,yrange=yrange,xrange=xrange,
                                                   colortype=cdattype,colorcode=True,cdatvec=cdatvecALL,
                                                   point_text=point_textALL,photoionizationplotparam=None,
                                                   histaxes=histaxes,Nbins=Nhistbins,
                                                   overwrite=overwrite,verbose=verbose)

    if performlinearfit & (len(cdatvecALL) > 0) & (len(xvalues[np.isfinite(xvalues)]) > 2) & (len(yvalues[np.isfinite(yvalues)]) > 2):
        plotname   = plotname.replace('.pdf','_ODRfit2data.pdf')

        if ylog:
            yval_fit = np.log10(yvalues)
            yerr_fit = 0.434 * yerr/yvalues  # https://faculty.washington.edu/stuve/log_error.pdf
        else:
            yval_fit = yvalues
            yerr_fit = yerr

        if xlog:
            xval_fit = np.log10(xvalues)
            xerr_fit = np.log10(xerr)
        else:
            xval_fit = xvalues
            xerr_fit = 0.434 * xerr/xvalues  # https://faculty.washington.edu/stuve/log_error.pdf

        fitres,rp,pp,rs,ps = kbs.fit_function_to_data_with_errors_on_both_axes(xval_fit,yval_fit,xerr_fit,yerr_fit,
                                                                               fitfunction='linear',plotresults=plotname,
                                                                               returnCorrelationCoeffs=True)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_infodat_plotcols():
    """
    Function returning columns from info file to plot
    """
    colinfo  = {}
    # from JK's thesis catalog: /Users/kschmidt/work/catalogs/MUSE_GTO/kerutt_LAEparameters190926_EWs_0_clumps_ratio_line_props.fits
    # colinfo['lyaew_b2']      = ['EW_0_beta_beta2','EW_0_beta_beta2_error',             'EW(Ly$\\alpha$) [\AA] ($\\beta = -2$)']
    # colinfo['lyaew_many']    = ['EW_0_beta_linear_many','EW_0_beta_linear_many_error', 'EW(Ly$\\alpha$) [\AA] ($\\beta$ linear multiband fit)']
    # colinfo['lyafwhm_a']     = ['fwhm_a_jk','fwhm_a_std_jk',                           'FWHM(Ly$\\alpha$) [\AA]']
    # colinfo['lyafwhm_kms']   = ['fwhm_kms_jk','fwhm_kms_std_jk',                       'FWHM(Ly$\\alpha$) [km/s]']
    # colinfo['beta_many']     = ['beta_linear_many','beta_linear_many_error',           '$\\beta$ (linear multiband fit)']
    # colinfo['beta_beta2']    = ['beta_beta2','beta_beta2_error',                       '$\\beta$ (fixed)']
    # colinfo['absmagUV_fix']  = ['abs_mag_UV_cont_beta2',None,                          'M(UV) ($\\beta = -2$)']
    # colinfo['absmagUV_many'] = ['abs_mag_UV_cont_linear_many',None,                    'M(UV) ($\\beta$ linear multiband fit)']
    # colinfo['magUV_fix']     = ['mag_UV_cont_beta2',None,                              'm(UV) ($\\beta = -2$)']
    # colinfo['magUV_many']    = ['mag_UV_cont_linear_many',None,                        'm(UV) ($\\beta$ linear multiband fit)']
    # colinfo['redshift']      = ['redshift',None,                                       '$z$']

    # from JK's 100 field catalog: /Users/kschmidt/work/catalogs/MUSE_GTO/kerutt_LAEparameters200709_EWs_all_fields_v0p9.fits
    colinfo['lyaew_JKmed']     = ['EW_0_beta_own_median','EW_0_beta_own_median_error',               'EW$_0$(Ly$\\alpha$) [\AA]'] # ($\\beta = -1.97$)
    colinfo['lyaew_b2']        = ['EW_0_beta_beta2','EW_0_beta_beta2_error',                         'EW$_0$(Ly$\\alpha$) [\AA] ($\\beta = -2$)']
    colinfo['lyaew_many']      = ['EW_0_beta_linear_many','EW_0_beta_linear_many_error',             'EW$_0$(Ly$\\alpha$) [\AA] ($\\beta$ linear multiband fit)']
    colinfo['lyafwhm_a']       = ['fwhm_a','fwhm_A_err',                                             'FWHM(Ly$\\alpha$) [\AA]']
    colinfo['lyafwhm_kms']     = ['fwhm_kms','fwhm_kms_err',                                         'FWHM(Ly$\\alpha$) [km/s]']
    colinfo['beta_many']       = ['beta_linear_many','beta_linear_many_error',                       '$\\beta$'] # (linear multiband fit)
    colinfo['absmagUV_-2']     = ['abs_mag_UV_cont_beta2','abs_mag_UV_cont_beta2_error',             'M(UV,1500) ($\\beta = -2$)']
    colinfo['absmagUV_many']   = ['abs_mag_UV_cont_linear_many','abs_mag_UV_cont_linear_many_error', 'M(UV,1500) ($\\beta$ linear multiband fit)']
    colinfo['absmagUV_median'] = ['abs_mag_UV_cont_own_median','abs_mag_UV_cont_own_median_error',   'M(UV,1500)'] #($\\beta = -1.97$)
    colinfo['magUV_-2']        = ['mag_UV_cont_beta2','mag_UV_cont_beta2_error',                     'm(UV,1500) ($\\beta = -2$)']
    colinfo['magUV_many']      = ['mag_UV_cont_linear_many','mag_UV_cont_linear_many_error',         'm(UV,1500) ($\\beta$ linear multiband fit)']
    colinfo['magUV_median']    = ['mag_UV_cont_own_median','mag_UV_cont_own_median_error',           'm(UV,1500)'] # ($\\beta = -1.97$)
    colinfo['zsys']            = ['z_vac','z_vac_err',                                               '$z_\\textrm{sys}$ (Verhamme et al. (2018) approx.)']
    colinfo['lyaflux']         = ['line_flux','line_flux_error',                'F(Ly$\\alpha$) [erg/s/cm$^2$]']
    colinfo['peaksep_kms']     = ['peak_sep_kms','peak_sep_kms_err',            'Peak separation [km/s]']
    colinfo['peaksep_A']       = ['peak_sep_A','peak_sep_A_err',                'Peak separation [\AA]']
    colinfo['loglLlya']        = ['logLLya','logLLya_err',                      'log$_{10}$(L(Ly$\\alpha$)/[erg/s])']
    colinfo['R_e']             = ['R_e','R_e_error',                            'R$_e$ [kpc]']
    colinfo['zsysFWHM']        = ['z_vac_fwhm','z_vac_fwhm_err',                '$z_\\textrm{sys}$ from FWHM(Ly$\\alpha$) ']
    colinfo['zsysPS']          = ['z_vac_peak_sep','z_vac_peak_sep_err',        '$z_\\textrm{sys}$ from Ly$\\alpha$ Peak separation ']
    colinfo['zfit']            = ['z_vac_fit','z_vac_fit_err',                  '$z_\\textrm{Ly$\\alpha$ fit}$']

    for key in colinfo.keys():
        for ent in [0,1]:
            if colinfo[key][ent] is not None:
                colinfo[key][ent] = colinfo[key][ent]+'_jk100'

    colinfo['redshift']        = ['redshift',None,                                       '$z_\\textrm{lead line}$']

    return colinfo

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def estimate_EW0(lineratiofile,infofile,outputfile='default', vetfelis_included=False, fixbeta=False,
                 defaultbeta=-1.97, overwrite=False, s2nlimit=3.0,fcontverbose=False, verbose=True):
    """
    Estimate the rest-frame EWs based on the the Kerutt+20 counterpart assignments and GALFIT magnitudes.
    For objects without matches to the Kerutt+20 catalog, 3D-HST or Rafelski photometry is used instead.

    --- INPUT ---

    lineratiofile        File containing the line ratio and flux estimates to use for the EW0 estimates.
    infofile             The general info catalog for the UVES sample including the 100 fields results from
                         Kerutt+20.
    outputfile           File name of output. If 'default' the output wil be stored to a file named like
                         and at the location of the lineratiofile appended "_EW0estimates"
    vetfelis_included    Set to True if the FELIS vetting is accounted for in lineartiofile in which case
                         the fluxes in the input catalog are Nsigma limits already.
    fixbeta              To fix the beta slope instead of using the values from the infofile provide fixed
                         beta value here.
    defaultbeta          If fixbeta=False and no beta value was found in the infofile, then use this default
                         beta value to obtain flux at emission line location in EW estimate.
    overwrite            Overwrite output if it already exists?
    s2nlimit             The S/N of the limits to provide in the output catalog.
    fcontverbose         Toggle verbosity on the continuume estimates
    verbose              Toggle verbosity


    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves

    parentdir     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'
    lineratiofile = parentdir+'back2backAnalysis_200213/fluxratios_FELISmatch2all_200213_postFELISvetting.txt'
    infofile      = parentdir+'objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    EW0data       = uves.estimate_EW0(lineratiofile,infofile,fixbeta=-1.97,vetfelis_included=True,overwrite=False)


    """
    if verbose: print(' - Estimating the continuum flux level for each source in \n   '+lineratiofile)
    EW0dic = collections.OrderedDict()
    EW0dic['id'] = np.array([])

    if verbose: print(' - Loading flux ratio data and infofile with Kerutt+20 flux estimates ')
    infofiledat     = afits.open(infofile)[1].data
    frdatBadFMT     = np.genfromtxt(lineratiofile,skip_header=7,dtype=None,comments='#',names=True)
    fmt             = 'd,12a,'+','.join((len(frdatBadFMT.dtype.names)-2)*['d'])
    fluxratiodatALL = np.genfromtxt(lineratiofile,skip_header=7,dtype=fmt,comments='#',names=True)
    ids             = fluxratiodatALL['id']
    pointings       = fluxratiodatALL['pointing']
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Initializing output ')
    if outputfile == 'default':
        outputfile = lineratiofile.replace('.txt','_EW0estimates.txt')
        if outputfile == lineratiofile:
            sys.exit(' The default output name is the same as the lineratiofile (missing .txt?) which is \n '+lineratiofile)
    if os.path.isfile(outputfile) & (overwrite == False):
        sys.exit(' Overwrite=False and the expected output '+outputfile+' already exists so exiting')
    else:
        fout = open(outputfile,'w')
        fout.write('# Estimated rest-frame EWs for the objects and fluxes in \n#  '+lineratiofile+'\n')
        fout.write('# Output generated with uves.estimate_EW0() on  '+kbs.DandTstr2()+'\n')
        if fixbeta:
            strout = '# Fixing beta = '+str(fixbeta)+'\n'
            if verbose: print(strout.replace('# ',' - '))
            fout.write(strout)
        else:
            strout = '# beta from Kerutt et al (2020) if available, otherwise beta = '+str(defaultbeta)+' is used.\n'
            if verbose: print(strout.replace('# ',' - '))
            fout.write(strout)
        fout.write('# EW0 upper(lower) limits are indicated by EW0err = +99(-99) and are '+str(s2nlimit)+'sigma limits. 0s indicate that both line flux and continuum flux were <'+str(s2nlimit)+'sigma \n')

        bandskey = {'F275W':275,'F336W':336,'F435W':435,
                    'F606W':606,'F775W':775,'F814W':814,'F850LP':850,
                    'F105W':105,'F125W':125,'F160W':160}
        fout.write('# The reference continuum band (before extrapolation to line location using beta) is based on the band indicated in the cont_band columns and use the keys: \n# '+str(bandskey)+'\n')
        fout.write('# Only the bands F435W, F606W, F775W, F814W, F125W, and F160W have fluxes estimated in the catalog by Kerutt et al. (2020). \n')
        fout.write('# The photref columns provides a reference to the flux catalog used: Kerutt (1), 3D-HST (2) or Rafelski (3).\n#\n')
    prefkey = {'Kerutt':1, '3D-HST':2, 'Rafelski':3}
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading photometric catalogs for continuum estimates when not match to Kerutt+20 catalog')
    catSkeltonGS  = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
    datSkeltonGS  = afits.open(catSkeltonGS)[1].data
    catSkeltonCOS = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat.FITS'
    datSkeltonCOS = afits.open(catSkeltonCOS)[1].data
    catRafelski   = '/Users/kschmidt/work/catalogs/rafelski/uvudf_rafelski_2015.fits'
    datRafelski   = afits.open(catRafelski)[1].data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fluxcols = []
    for colname in fluxratiodatALL.dtype.names:
        if colname.startswith('f_'):
            fluxcols.append(colname)
    Nfluxcol   = len(fluxcols)
    outcolumns = '# id  pointing  beta  '+'    '.join([('EW0_'+fc[2:]+' '+'EW0err_'+fc[2:]+' '+'contband_'+fc[2:]+' '+'contmagAB_'+fc[2:]+' '+'contmagABerr_'+fc[2:]+' '+'photref_'+fc[2:]) for fc in fluxcols])#+'  photref  contmag '
    fout.write(outcolumns+'\n')

    sizeColSet = 6 # size of set of columns: EW0, EW0err, contband, magAB, magABerr, photref
    EW0array = np.zeros([len(ids),Nfluxcol*sizeColSet])*np.nan
    magarray = np.zeros([len(ids),Nfluxcol*sizeColSet])*np.nan

    if verbose: print(' - Estimating continuum level based on band fluxes and from that EW0 for all lines for object: ')
    betavals = np.zeros(len(ids))
    for ii, id in enumerate(ids):
        if verbose:
            infostr = '   '+str(int(id))+'   ( '+str("%.5d" % (ii+1))+' / '+str("%.5d" % len(ids))+' ) '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()
        obj_infoent = np.where(infofiledat['id'].astype(int) == int(id))[0]
        redshift    = infofiledat['redshift'][obj_infoent][0]
        for ff, fluxcol in enumerate(fluxcols):
            line_name    = fluxcol[2:]
            line_wave    = uves.linewavesUVES(line_name) * (1.0 + redshift)
            line_flux    = fluxratiodatALL[fluxcol][ii]
            line_fluxerr = fluxratiodatALL[fluxcol.replace('f_','ferr_')][ii]

            if vetfelis_included: # FELIS vetting accounted for meaning fluxes are Nsigma limits already
                if line_fluxerr == 99:
                    line_fluxerr = line_flux / s2nlimit
                else:
                    line_fluxerr = line_fluxerr / s2nlimit

            if np.isfinite(line_flux) & (line_flux > 0):
                if fixbeta:
                    beta    = fixbeta
                else:
                    if infofiledat['beta_linear_many_jk100'][obj_infoent] != 0.0:
                        beta  = infofiledat['beta_linear_many_jk100'][obj_infoent]
                    else:
                        beta  = defaultbeta

                f_conts, photrefs, magABs = uves.estimate_fcont(infofiledat,obj_infoent,line_wave,fixbeta=beta,verbose=fcontverbose,
                                                                datSkeltonGS=datSkeltonGS, datSkeltonCOS=datSkeltonCOS, datRafelski=datRafelski)

                betavals[ii] = beta
                # selecting continuum band to extrapolate from as:
                # The band with effective wavelength closest to the line that does not contain Lya or
                # the line itself and is not from WFC3_UVIS
                bands        = f_conts.keys()
                goodphotrefs = []
                goodbands    = []
                for bb, band in enumerate(bands):
                    if (band not in ['F275W', 'F336W']) & (f_conts[band] !=  (-99, -99)) & (f_conts[band] !=  (np.nan, np.nan)):
                        goodbands.append(band)
                        goodphotrefs.append(photrefs[bb])
                goodbands    = np.asarray(goodbands)
                goodphotrefs = np.asarray(goodphotrefs)

                if len(goodbands) > 0:
                    if 'Kerutt' in goodphotrefs: # Choose Kerutt et al. (2020) photometry if available over the others
                        selectent    = np.where(goodphotrefs == 'Kerutt')[0]
                        goodphotrefs = goodphotrefs[selectent]
                        goodbands    = goodbands[selectent]
                    elif 'Rafelski' in goodphotrefs: # Chose Rafelski over 3D-HST and LyaOrLineInBand
                        selectent    = np.where(goodphotrefs == 'Rafelski')[0]
                        goodphotrefs = goodphotrefs[selectent]
                        goodbands    = goodbands[selectent]
                    elif '3D-HST' in goodphotrefs: # Chose 3D-HST over LyaOrLineInBand
                        selectent    = np.where(goodphotrefs == '3D-HST')[0]
                        goodphotrefs = goodphotrefs[selectent]
                        goodbands    = goodbands[selectent]

                    wavediff   = np.asarray([np.abs(uves.band_waveeff(gb) - line_wave) for gb in goodbands])
                    mindiffent = np.where(wavediff == np.min(wavediff))[0][0]
                    cont_band  = goodbands[mindiffent]
                    photoref   = goodphotrefs[mindiffent]

                    f_cont, f_cont_err   = f_conts[cont_band]
                    f_cont_restframe     = f_cont     * (1.0 + redshift)
                    f_cont_err_restframe = f_cont_err * (1.0 + redshift)

                    magAB, magABerr      = magABs[cont_band]
                    if magABerr == 0.0: # denote lower limit (fainter than) magnitude by -99
                        magABerr = -99

                    if (line_flux/line_fluxerr > s2nlimit) & (f_cont/f_cont_err > s2nlimit):
                        EW0array[ii,ff*sizeColSet]   = line_flux / f_cont_restframe
                        EW0array[ii,ff*sizeColSet+1] = np.sqrt((line_fluxerr/line_flux)**2 +
                                                      (f_cont_err_restframe/f_cont_restframe)**2) * EW0array[ii,ff*sizeColSet]
                        EW0array[ii,ff*sizeColSet+3] = magAB
                        EW0array[ii,ff*sizeColSet+4] = magABerr
                    elif (line_flux/line_fluxerr <= s2nlimit) & (f_cont/f_cont_err > s2nlimit):
                        EW0array[ii,ff*sizeColSet]   = s2nlimit * line_fluxerr / f_cont_restframe
                        EW0array[ii,ff*sizeColSet+1] = +99 # indicate that EW0 value is "s2nlimit sigma upper limit"
                        EW0array[ii,ff*sizeColSet+3] = magAB
                        EW0array[ii,ff*sizeColSet+4] = magABerr
                    elif (line_flux/line_fluxerr > s2nlimit) & (f_cont/f_cont_err <= s2nlimit):
                        EW0array[ii,ff*sizeColSet]   = line_flux / (s2nlimit * f_cont_err_restframe)
                        EW0array[ii,ff*sizeColSet+1] = -99 # indicate that EW0 value is "s2nlimit sigma lower limit"
                        EW0array[ii,ff*sizeColSet+3] = magAB
                        EW0array[ii,ff*sizeColSet+4] = -99   # quote magAB as lower limit (fainter than) if S/N < s2nlimit
                    else:
                        EW0array[ii,ff*sizeColSet]   = 0.0
                        EW0array[ii,ff*sizeColSet+1] = 0.0
                        EW0array[ii,ff*sizeColSet+3] = magAB
                        EW0array[ii,ff*sizeColSet+4] = -99   # quote magAB as upper limit (fainter than) if S/N < s2nlimit

                    EW0array[ii,ff*sizeColSet+2] = bandskey[cont_band.upper()]
                    EW0array[ii,ff*sizeColSet+5] = prefkey[photoref]

        outstr = str(int(id))+' '+pointings[ii]+' '+str("%10.2f" % betavals[ii])+' '+\
                 ' '.join([str("%10.4f" % ff) for ff in EW0array[ii,:]])
        fout.write(outstr.replace('.0000','     ')+' \n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fout.close()
    if verbose: print('\n - Wrote the EW0 estimates to \n   '+outputfile)
    fmt = ','.join((2+Nfluxcol*3)*['d'])
    EW0dat = np.genfromtxt(outputfile,skip_header=10,dtype=fmt,comments='#',names=True)
    return EW0dat

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def estimate_fcont(infofiledat,obj_infoent,wavelength,fixbeta=False,photmatchlim=0.5,
                   datSkeltonGS=None, datSkeltonCOS=None, datRafelski=None, verbose=True):
    """
    Estimate the rest-frame continuum flux at a given wavelength for an object in the infofile using the
    observed broad band magnitudes and the beta estimate.

    """
    objid         = infofiledat['id'][obj_infoent][0]
    redshift      = infofiledat['redshift'][obj_infoent][0]
    if fixbeta:
        beta      = fixbeta
    else:
        beta      = infofiledat['beta'][obj_infoent][0]

    if datSkeltonGS is None:
        catSkeltonGS  = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
        datSkeltonGS  = afits.open(catSkeltonGS)[1].data

    if datSkeltonCOS is None:
        catSkeltonCOS = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat.FITS'
        datSkeltonCOS = afits.open(catSkeltonCOS)[1].data

    if datRafelski is None:
        catRafelski   = '/Users/kschmidt/work/catalogs/rafelski/uvudf_rafelski_2015.fits'
        datRafelski   = afits.open(catRafelski)[1].data

    bands         = ['F275W','F336W','F435W','F606W','F775W','F814W','F880LP','F105W','F125W','F160W']

    try:
        bandswithLya  = uves.wavelength_in_bands( 1216 * (redshift + 1))
    except:
        bandswithLya  = uves.wavelength_in_bands( 1216 * (redshift + 1))

    bandswithline = uves.wavelength_in_bands( wavelength )

    f_cont        = collections.OrderedDict()
    photrefs      = ['']*len(bands)
    magABs        = collections.OrderedDict()
    for bb, band in enumerate(bands):
        if verbose: print(' - Getting flux estimate for band '+band)
        if (band not in bandswithLya) & (band not in bandswithline):
            # - - - - - first check infofile for flux estimate - - - - -
            if band in ['F105W','F125W','F160W']:
                hstinst = 'wfc3'
            else:
                hstinst = 'acs'

            try:
                if infofiledat['use_'+hstinst+'_'+band.lower()[1:]+'_jk100' ][obj_infoent][0]:
                    f_ref         = infofiledat[ 'flux_'+hstinst+'_'+band.lower()[1:]+'_jk100' ][obj_infoent][0] * 1e20
                    f_ref_err     = infofiledat[ 'flux_error_'+hstinst+'_'+band.lower()[1:]+'_jk100' ][obj_infoent][0] * 1e20
                    photrefs[bb]  = 'Kerutt'
                    magAB         = infofiledat[ 'mag_'+hstinst+'_'+band.lower()[1:]+'_jk100' ][obj_infoent][0]
                    magABerr      = infofiledat[ 'mag_error_'+hstinst+'_'+band.lower()[1:]+'_jk100' ][obj_infoent][0]
                else:
                    f_ref     = 0.0
                    f_ref_err = 0.0
            except:
                f_ref     = 0.0
                f_ref_err = 0.0

            # - - - - - if not there try to get flux from photometric catalogs - - - - -
            if f_ref == 0.0:
                if verbose: print('   was not in the infofile; trying photmetric catalogs')
                if str(objid)[0] in ['6', '7']:
                    phot_id    = infofiledat['id_rafelski'][obj_infoent][0]
                    phot_match = infofiledat['sep_rafelski'][obj_infoent][0]
                    phot_ent   = np.where(datRafelski['id'] == phot_id)[0]

                    if phot_match <= photmatchlim:
                        try:
                            f_cont_jy     =  datRafelski['FLUX_'+band.upper()][phot_ent] * 1e-6
                            f_cont_jy_err =  datRafelski['FLUXERR_'+band.upper()][phot_ent] * 1e-6

                            # ------ from http://www.stsci.edu/~strolger/docs/UNITS.txt: ------
                            # [Y Jy]            = 3.33564095E+04 * [X1 erg/cm^2/s/A] * [X2 A]^2
                            # [Y erg/cm^2/s/A]  = 2.99792458E-05 * [X1 Jy] / [X2 A]^2
                            f_ref     = f_cont_jy     * 2.99792458E-05 / uves.band_waveeff(band)**2.0 / 1e-20
                            f_ref_err = f_cont_jy_err * 2.99792458E-05 / uves.band_waveeff(band)**2.0 / 1e-20
                            photrefs[bb]  = 'Rafelski'
                            magAB         =   datRafelski['MAG_'+band.upper()][phot_ent]
                            magABerr      =   datRafelski['MAGERR_'+band.upper()][phot_ent]
                        except:
                            f_ref     = 0.0
                            f_ref_err = 0.0
                else:
                    if str(objid)[0] in ['1','3','4']:
                        datSkel = datSkeltonGS
                    elif str(objid)[0] in ['2']:
                        datSkel = datSkeltonCOS
                    else:
                        sys.exit(' - No Skelton data for ID = '+str(objid))

                    phot_id    = infofiledat['id_skelton'][obj_infoent][0]
                    phot_match = infofiledat['sep_skelton'][obj_infoent][0]
                    phot_ent   = np.where(datSkel['id'] == phot_id)[0][0]

                    if phot_match <= photmatchlim:
                        try:
                            # magAB    = -2.5*log10(f_nu/Jy) + 8.90    <- https://en.wikipedia.org/wiki/AB_magnitude
                            skelflux   = datSkel['f_'+band.lower()][phot_ent]
                            skelerr    = datSkel['e_'+band.lower()][phot_ent]
                            magAB      = 25.0-2.5*np.log10(skelflux)
                            magABerr   = np.abs(-2.5* np.log10(np.e) / skelflux * skelerr)

                            f_cont_jy      = 10**( (magAB-8.90)/-2.5 )
                            f_cont_jy_err  = np.abs(10**((magAB-8.90)/-2.5) * np.log(10)/-2.5 * magABerr) # deriv. from Schaums 15.28
                            # ------ from http://www.stsci.edu/~strolger/docs/UNITS.txt: ------
                            # [Y Jy]            = 3.33564095E+04 * [X1 erg/cm^2/s/A] * [X2 A]^2
                            # [Y erg/cm^2/s/A]  = 2.99792458E-05 * [X1 Jy] / [X2 A]^2
                            f_ref     = f_cont_jy     * 2.99792458E-05 / uves.band_waveeff(band)**2.0 / 1e-20
                            f_ref_err = f_cont_jy_err * 2.99792458E-05 / uves.band_waveeff(band)**2.0 / 1e-20
                            photrefs[bb]  = '3D-HST'
                        except:
                            f_ref     = 0.0
                            f_ref_err = 0.0
            else:
                phot_match = 0.0
                if verbose: print('   found flux estimate in infofile')
            # - - - - - if still not available set to NaNs - - - - -
            if f_ref == 0.0:
                if phot_match > photmatchlim:
                    if verbose: print('   Photometric match further than '+str(photmatchlim)+' arcsec from object coordinates.')
                else:
                    if verbose: print('   Hmm; no value in photmetric catalogs either. Returning NaNs')
                f_cont[band] = np.nan, np.nan
                photrefs[bb]  = 'None'
                magABs[band]  = np.nan, np.nan
            else:
                if verbose: print('   Using flux value to estimate the continuum strength at lambda='+
                                  str(wavelength)+' using beta='+str(beta))
                f_wave       = uves.estimate_continuumlevel_viaBeta(wavelength, uves.band_waveeff(band), f_ref, beta, verbose=False)
                f_wave_err   = f_ref_err * np.abs((wavelength / uves.band_waveeff(band))**beta)
                f_cont[band] = f_wave, f_wave_err
                magABs[band]  = magAB, magABerr
        else:
            if verbose: print('   Band contains either Lya or the wavelength itself, so returning -99s')
            f_cont[band] = -99, -99
            photrefs[bb] = 'LyaOrLineInBand'
            magABs[band] = np.nan, np.nan

    if '' in photrefs:
        print(' ---> There should not be empty "photrefs" entries; stopping to enable further investigation ')
        pdb.set_trace()
    return f_cont, photrefs, magABs
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_SED(infofiledat,objent,photmatchlim=0.5,
            datSkeltonGS=None, datSkeltonCOS=None, datRafelski=None, # Enable providing data from outside function
            verbose=True):
    """
    Return flux and ABmag SEDs from photometric catalogs of given object
    Similar to uves.estimate_fcont() but returns all photometric values

    """
    objid         = infofiledat['id'][objent]

    if datSkeltonGS is None:
        catSkeltonGS  = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
        datSkeltonGS  = afits.open(catSkeltonGS)[1].data

    if datSkeltonCOS is None:
        catSkeltonCOS = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat.FITS'
        datSkeltonCOS = afits.open(catSkeltonCOS)[1].data

    if datRafelski is None:
        catRafelski   = '/Users/kschmidt/work/catalogs/rafelski/uvudf_rafelski_2015.fits'
        datRafelski   = afits.open(catRafelski)[1].data

    bands         = ['F275W','F336W','F435W','F606W','F775W','F814W','F880LP','F105W','F125W','F160W']

    f_out         = collections.OrderedDict()
    photrefs      = ['']*len(bands)
    magABs        = collections.OrderedDict()
    for bb, band in enumerate(bands):
        # if verbose: print(' - Getting flux estimate for band '+band)
        if band in ['F105W','F125W','F160W']:
            hstinst = 'wfc3'
        else:
            hstinst = 'acs'

        try:
            if infofiledat['use_'+hstinst+'_'+band.lower()[1:]+'_jk100' ][objent][0]:
                f_ref         = infofiledat[ 'flux_'+hstinst+'_'+band.lower()[1:]+'_jk100' ][objent][0] * 1e20
                f_ref_err     = infofiledat[ 'flux_error_'+hstinst+'_'+band.lower()[1:]+'_jk100' ][objent][0] * 1e20
                photrefs[bb]  = 'Kerutt'
                magAB         = infofiledat[ 'mag_'+hstinst+'_'+band.lower()[1:]+'_jk100' ][objent][0]
                magABerr      = infofiledat[ 'mag_error_'+hstinst+'_'+band.lower()[1:]+'_jk100' ][objent][0]
            else:
                f_ref     = 0.0
                f_ref_err = 0.0
        except:
            f_ref     = 0.0
            f_ref_err = 0.0

        # - - - - - if not there try to get flux from photometric catalogs - - - - -
        if f_ref == 0.0:
            # if verbose: print('   was not in the infofile; trying photometric catalogs')
            if str(objid)[0] in ['6', '7']:
                phot_id    = infofiledat['id_rafelski'][objent]
                phot_match = infofiledat['sep_rafelski'][objent]
                phot_ent   = np.where(datRafelski['id'] == phot_id)[0]

                if phot_match <= photmatchlim:
                    try:
                        f_cont_jy     =  datRafelski['FLUX_'+band.upper()][phot_ent] * 1e-6
                        f_cont_jy_err =  datRafelski['FLUXERR_'+band.upper()][phot_ent] * 1e-6

                        # ------ from http://www.stsci.edu/~strolger/docs/UNITS.txt: ------
                        # [Y Jy]            = 3.33564095E+04 * [X1 erg/cm^2/s/A] * [X2 A]^2
                        # [Y erg/cm^2/s/A]  = 2.99792458E-05 * [X1 Jy] / [X2 A]^2
                        f_ref     = f_cont_jy     * 2.99792458E-05 / uves.band_waveeff(band)**2.0 / 1e-20
                        f_ref_err = f_cont_jy_err * 2.99792458E-05 / uves.band_waveeff(band)**2.0 / 1e-20
                        photrefs[bb]  = 'Rafelski'
                        magAB         =   datRafelski['MAG_'+band.upper()][phot_ent]
                        magABerr      =   datRafelski['MAGERR_'+band.upper()][phot_ent]
                    except:
                        f_ref     = 0.0
                        f_ref_err = 0.0
            else:
                if str(objid)[0] in ['1','3','4']:
                    datSkel = datSkeltonGS
                elif str(objid)[0] in ['2']:
                    datSkel = datSkeltonCOS
                else:
                    sys.exit(' - No Skelton data for ID = '+str(objid))

                phot_id    = infofiledat['id_skelton'][objent]
                phot_match = infofiledat['sep_skelton'][objent]
                phot_ent   = np.where(datSkel['id'] == phot_id)[0]

                if phot_match <= photmatchlim:
                    try:
                        # magAB    = -2.5*log10(f_nu/Jy) + 8.90    <- https://en.wikipedia.org/wiki/AB_magnitude
                        skelflux   = datSkel['f_'+band.lower()][phot_ent]
                        skelerr    = datSkel['e_'+band.lower()][phot_ent]
                        magAB      = 25.0-2.5*np.log10(skelflux)
                        magABerr   = np.abs(-2.5* np.log10(np.e) / skelflux * skelerr)

                        f_cont_jy      = 10**( (magAB-8.90)/-2.5 )
                        f_cont_jy_err  = np.abs(10**((magAB-8.90)/-2.5) * np.log(10)/-2.5 * magABerr) # deriv. from Schaums 15.28
                        # ------ from http://www.stsci.edu/~strolger/docs/UNITS.txt: ------
                        # [Y Jy]            = 3.33564095E+04 * [X1 erg/cm^2/s/A] * [X2 A]^2
                        # [Y erg/cm^2/s/A]  = 2.99792458E-05 * [X1 Jy] / [X2 A]^2
                        f_ref     = f_cont_jy     * 2.99792458E-05 / uves.band_waveeff(band)**2.0 / 1e-20
                        f_ref_err = f_cont_jy_err * 2.99792458E-05 / uves.band_waveeff(band)**2.0 / 1e-20
                        photrefs[bb]  = '3D-HST'
                    except:
                        f_ref     = 0.0
                        f_ref_err = 0.0
        else:
            phot_match = 0.0
            if verbose: print('   found flux estimate in infofile')
        # - - - - - if still not available set to NaNs - - - - -
        if f_ref == 0.0:
            # if phot_match > photmatchlim:
            #     if verbose: print('   Photometric match further than '+str(photmatchlim)+' arcsec from object coordinates.')
            # else:
            #     if verbose: print('   Hmm; no value in photometric catalogs either. Returning NaNs')
            f_out[band]    = np.nan, np.nan
            photrefs[bb]   = 'None'
            magABs[band]   = np.nan, np.nan
        else:
            f_out[band]    = f_ref, f_ref_err
            magABs[band]   = magAB, magABerr

    if '' in photrefs:
        print(' ---> There should not be empty "photrefs" entries; stopping to enable further investigation ')
        pdb.set_trace()
    return f_out, photrefs, magABs
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def band_waveeff(bandname):
    """
    Returning the effective wavelength of an HST band

    """
    #  WFC3-UVIS filters: F225W, F275W, and F336W
    #  ACS-WFC  optical filters: F435W, F606W, F775W, and F850LP
    #  WFC3-IR filters: F105W, F125W, F140W, and F160W
    #
    # effective wavelengths of ACS-WFC, WFC3-UVIS1 and WFC3-IR from
    # http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=HST&gname2=ACS_WFC

    bands_wavecen = collections.OrderedDict()
    bands_wavecen['F275W' ] = 2720.9
    bands_wavecen['F336W' ] = 3359.4
    bands_wavecen['F435W' ] = 4341.9
    bands_wavecen['F606W' ] = 5810.8
    bands_wavecen['F775W' ] = 7652.5
    bands_wavecen['F814W' ] = 7972.9
    bands_wavecen['F880LP'] = 9008.7
    bands_wavecen['F105W' ] = 10431.7
    bands_wavecen['F125W' ] = 12364.6
    bands_wavecen['F160W' ] = 15279.1

    if bandname == 'all':
        return bands_wavecen
    else:
        return bands_wavecen[bandname.upper()]

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_EW0estimates(lineratiofile, plotbasename, infofile, EW0file, colorvar_obj='EW_0', point_text=None, showlimits=True,
                      addliteraturevalues=True, ErrNsigma=1.0, vshiftmax=1e5, obj2show='all', xlog=True, ylog=True,
                      lyaEWtype='lyaew_b2', overwrite=False, verbose=True, litsymboldot=False):
    """
    Function to plot the output containing EW0 estimates generated with uves.estimate_EW0()

    --- INPUT ---


    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves

    B2Bdir          = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/back2backAnalysis_200213/'
    FRfile_postvet  = B2Bdir+'fluxratios_FELISmatch2all_200213_postFELISvetting.txt'
    plotbase_EWest  = B2Bdir+'plots_FRaEW/fluxratios_FELISmatch2all_200213_postFELISvetting_EW0estimates'
    infofile        = B2Bdir+'../objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    EWest_postvet   = B2Bdir+'fluxratios_FELISmatch2all_200213_postFELISvetting_EW0estimates.txt'

    addliteraturevalues = True
    showlimits          = False
    lyaEWtype           = 'lyaew_JKmed' # from
    ErrNsigma           = 3.0 # the size of the errorbars (limits are kept at 3sigma limits)
    obj2show            = 'all_nodup' # [102014087,122022112,157001017,602922055,606813186,721590802] #

    uves.plot_EW0estimates(FRfile_postvet, plotbase_EWest, infofile, EWest_postvet, lyaEWtype=lyaEWtype, ErrNsigma=ErrNsigma, vshiftmax=1000.0, overwrite=True, colorvar_obj='redshift', xlog=True,ylog=True, addliteraturevalues=addliteraturevalues, showlimits=showlimits, obj2show=obj2show)


    """
    if verbose: print(' - Loading flux ratio and EW data to plot ')
    fluxratiodatALL = np.genfromtxt(lineratiofile,skip_header=7,dtype='d',comments='#',names=True)
    EW0datALL       = np.genfromtxt(EW0file,skip_header=10,dtype='d',comments='#',names=True)

    if verbose: print(' - Loading infofile data for expanded selection and plotting ')
    infofiledat      = afits.open(infofile)[1].data
    infofiledat      = infofiledat[np.where((infofiledat['id']<4.9e8) | (infofiledat['id']>5.9e8))[0]] # ignoring UDF MW mock ids

    if verbose: print(' - Selecting objects to show ')
    if str(obj2show).lower() == 'all':
        showent = np.arange(len(fluxratiodatALL))
    elif str(obj2show).lower() == 'all_nodup':
        idlist2show = infofiledat[(infofiledat['duplicationID'] == 0.0) & (infofiledat['redshift'] <= 4.955)]['id']
        showent     = np.array([])
        for objid in idlist2show:
            objent = np.where( fluxratiodatALL['id'].astype(int) == objid)[0]
            if len(objent) > 0:
                showent = np.append(showent,objent)
    elif str(obj2show).lower() == 'none':
        showent = np.array([0])
    elif str(obj2show).lower() == 'udf10_goodspec_only':
        sample = 'udf10'
        ids_badTDOSEspec, ids_goodTDOSEspec = uves.summarize_tdosevetting(returnsample=sample,verbose=verbose)
        idlist2show = ids_goodTDOSEspec
    elif str(obj2show).lower() == 'udf10_badspec_only':
        sample = 'udf10'
        ids_badTDOSEspec, ids_goodTDOSEspec = uves.summarize_tdosevetting(returnsample=sample,verbose=verbose)
        idlist2show = ids_badTDOSEspec
    else:
        idlist2show = obj2show

        showent = np.array([])
        for objid in idlist2show:
            objent = np.where( fluxratiodatALL['id'].astype(int) == objid)[0]
            if len(objent) > 0:
                showent = np.append(showent,objent)
            else:
                if verbose: print(' - The ID '+str(objid)+' was not found in the lineratiofile ')

    Nselspec    = len(showent)
    if Nselspec > 0:
        fluxratiodat = fluxratiodatALL[showent.astype(int)]
        EW0dat       = EW0datALL[showent.astype(int)]
        if verbose: print(' - '+str(Nselspec)+'/'+str(len(fluxratiodatALL))+
                          ' spectra in flux ratio summary (and EW0 file) satisfies the cuts')
    else:
        if verbose: print(' WARNING No flux ratio matches found in summary file satisfying cuts; returning...')
        #return

    if verbose: print(' - Pulling out Lya measurements from Kerutt et al. (2020) columns in infofile ')
    LyaEW         = np.zeros(len(fluxratiodat['id']))*np.nan
    LyaEWerr      = np.zeros(len(fluxratiodat['id']))*np.nan

    LyaFWHM       = np.zeros(len(fluxratiodat['id']))*np.nan
    LyaFWHMerr    = np.zeros(len(fluxratiodat['id']))*np.nan

    LyaPeaksep    = np.zeros(len(fluxratiodat['id']))*np.nan
    LyaPeakseperr = np.zeros(len(fluxratiodat['id']))*np.nan

    MUV           = np.zeros(len(fluxratiodat['id']))*np.nan
    MUVerr        = np.zeros(len(fluxratiodat['id']))*np.nan
    beta          = np.zeros(len(fluxratiodat['id']))*np.nan
    betaerr       = np.zeros(len(fluxratiodat['id']))*np.nan
    Llya          = np.zeros(len(fluxratiodat['id']))*np.nan
    Llyaerr       = np.zeros(len(fluxratiodat['id']))*np.nan

    colinfo = uves.get_infodat_plotcols()
    for ii, id in enumerate(fluxratiodat['id']):
        infoent            = np.where(infofiledat['id'] == int(id))
        LyaEW[ii]          = infofiledat[colinfo[lyaEWtype][0]][infoent]
        LyaEWerr[ii]       = infofiledat[colinfo[lyaEWtype][1]][infoent]
        LyaFWHM[ii]        = infofiledat[colinfo['lyafwhm_kms'][0]][infoent]
        LyaFWHMerr[ii]     = infofiledat[colinfo['lyafwhm_kms'][1]][infoent]
        LyaPeaksep[ii]     = infofiledat[colinfo['peaksep_kms'][0]][infoent]
        LyaPeakseperr[ii]  = infofiledat[colinfo['peaksep_kms'][1]][infoent]
        MUV[ii]            = infofiledat[colinfo['absmagUV_median'][0]][infoent]
        MUVerr[ii]         = infofiledat[colinfo['absmagUV_median'][1]][infoent]
        beta[ii]           = infofiledat[colinfo['beta_many'][0]][infoent]
        betaerr[ii]        = infofiledat[colinfo['beta_many'][1]][infoent]
        Llya[ii]           = infofiledat[colinfo['loglLlya'][0]][infoent]
        Llyaerr[ii]        = infofiledat[colinfo['loglLlya'][1]][infoent]


    if colorvar_obj in fluxratiodat.dtype.names:
        cdatvec   = fluxratiodat[colorvar_obj]
    elif colorvar_obj in infofiledat.columns.names:
        cdatvec = np.zeros(len(fluxratiodat['id']))*np.nan
        for ii, id in enumerate(fluxratiodat['id']):
            infoent     = np.where(infofiledat['id'] == int(id))
            cdatvec[ii] = infofiledat[colorvar_obj][infoent]
    else:
        cdatvec = None

    if colorvar_obj.lower() == 's2n_ciii':
        cdattype = 's2n_ciii'
    elif colorvar_obj.lower() == 'redshift':
        cdattype = colorvar_obj.lower()
    elif colorvar_obj.lower() == 'ew_0':
        cdattype = colorvar_obj.lower()

    if point_text is not None:
        point_text = fluxratiodat['id'].astype(str)

    #------------------------------- Append literature measurements -------------------------------
    if addliteraturevalues:
        plotbasename = plotbasename+'_wLit'

        if verbose: print(' - Loading and appending data from catalog of literature observations ')
        litcat  = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/' \
                  'literaturecollection_emissionlinestrengths.fits'
        litdat  = afits.open(litcat)[1].data
        # litdat  = np.genfromtxt(litcat,names=True,skip_header=5,comments='#',dtype=None)
        Nlitobj = len(litdat)

        litarr  = np.array(np.zeros(Nlitobj)*np.nan,fluxratiodat.dtype)
        for lo, litobject in enumerate(litdat):
            for litcol in litdat.dtype.names:
                if litcol in fluxratiodat.dtype.names:
                    litarr[litcol][lo] = litdat[litcol][lo]
        fluxratiodat = np.hstack((fluxratiodat,litarr))

        litarr  = np.array(np.zeros(Nlitobj)*np.nan,EW0dat.dtype)
        for lo, litobject in enumerate(litdat):
            for litcol in litdat.dtype.names:
                if litcol in EW0dat.dtype.names:
                    litarr[litcol][lo] = litdat[litcol][lo]

        if litsymboldot:
            litarr['id'] = litarr['id']*0.0+990000000000

        EW0dat = np.hstack((EW0dat,litarr))

        LyaEW         = np.hstack((LyaEW,        litdat['EW0_Lya']))
        LyaEWerr      = np.hstack((LyaEWerr,     litdat['EW0err_Lya']))

        LyaFWHM       = np.hstack((LyaFWHM,      litdat['sigma_Lya']*2.355))
        LyaFWHMerr    = np.hstack((LyaFWHMerr,   litdat['sigmaerr_Lya']*2.355))

        LyaPeaksep    = np.hstack((LyaPeaksep,   litdat['sigma_Lya']*0.0))
        LyaPeakseperr = np.hstack((LyaPeakseperr,litdat['sigma_Lya']*0.0))

        if colorvar_obj.lower() == 's2n_ciii':
            cdatvec = np.hstack((cdatvec,litdat['s2n_CIII']))
        elif colorvar_obj.lower() == 'redshift':
            cdatvec = np.hstack((cdatvec,litdat['redshift']))
        elif colorvar_obj.lower() == 'ew_0':
            cdatvec = np.hstack((cdatvec,litdat['EW0_Lya']))
        else:
            sys.exit(' Selected color vector not enabled for literature values ')
        if verbose: print('   Appended data for '+str(Nlitobj)+' literature objects \n')

    #----------------------------------------------------------------------------------------------s

    if verbose: print(' - Defining plotting ranges and sets of parameters to plot ')
    if ylog:
        EW0_range_y = [0.1,900]
    else:
        EW0_range_y = [-5,42]

    if xlog:
        EW0_range_x = [0.1,900]
    else:
        EW0_range_x = [-5,42]

    if xlog:
        LyaEW_range = [0.1,900]
    else:
        LyaEW_range = [-100,320]

    if xlog:
        LyaPS_range = [10,10000]
    else:
        LyaPS_range = [0,1000]

    if xlog:
        LyaFWHM_range = [1,3000]
    else:
        LyaFWHM_range = [0,1000]

    if xlog:
        beta_range = [0.1,3]
    else:
        beta_range = [-3,-1] #[-4,3]

    if xlog:
        MUV_range = [10,25]
    else:
        MUV_range = [-25,-10]

    if xlog:
        Llya_range = [10,25]
    else:
        Llya_range = [39.8,43.8]

    linesetlist_EWs = []

    for yline in ['MgII', 'CIII','CIV', 'OIII', 'HeII', 'SiIII', 'NV']:
        linesetlist_EWs.append([['LyaEW',     'EW$_0$(Ly$\\alpha$) [\AA]',LyaEW,LyaEWerr],
                                yline   ,LyaEW_range, EW0_range_y,   None])
        linesetlist_EWs.append([['LyaFWHM',   'FWHM(Ly$\\alpha$) [km/s]',LyaFWHM,LyaFWHMerr],
                                yline   ,LyaFWHM_range, EW0_range_y,   None])
        linesetlist_EWs.append([['LyaPeaksep','Ly$\\alpha$ Peak Seperation [km/s]',LyaPeaksep,LyaPeakseperr],
                                yline   ,LyaPS_range, EW0_range_y,   None])

        if not addliteraturevalues:
            cdattype = 'zmanual'
            histaxes = False
            linesetlist_EWs.append([['logLlya','log$_{10}$(L(Ly$\\alpha$)/[erg/s])',Llya,Llyaerr],
                                    yline   ,Llya_range, EW0_range_y,   None])
            linesetlist_EWs.append([['MUVlya','M(UV,1500)',MUV,MUVerr],
                                    yline   ,MUV_range, EW0_range_y,   None])
            linesetlist_EWs.append([['betalya','$\\beta$',beta,betaerr],
                                    yline   ,beta_range, EW0_range_y,   None])
        else:
            histaxes  = True

    linesetlist_EWs.append(['CIII','CIV'   ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['CIII','OIII'  ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['CIII','HeII'  ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['CIII','MgII'  ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['CIII','NV'    ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['CIII','SiIII' ,EW0_range_x, EW0_range_y,   None])

    linesetlist_EWs.append(['CIV','OIII'   ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['CIV','HeII'   ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['CIV','MgII'   ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['CIV','NV'     ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['CIV','SiIII'  ,EW0_range_x, EW0_range_y,   None])

    linesetlist_EWs.append(['OIII','HeII'  ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['OIII','MgII'  ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['OIII','NV'    ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['OIII','SiIII' ,EW0_range_x, EW0_range_y,   None])

    linesetlist_EWs.append(['HeII','MgII'  ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['HeII','NV'    ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['HeII','SiIII' ,EW0_range_x, EW0_range_y,   None])

    linesetlist_EWs.append(['MgII','NV'    ,EW0_range_x, EW0_range_y,   None])
    linesetlist_EWs.append(['MgII','SiIII' ,EW0_range_x, EW0_range_y,   None])

    linesetlist_EWs.append(['NV','SiIII'   ,EW0_range_x, EW0_range_y,   None])

    Nhistbins = 30
    for lineset in linesetlist_EWs:
        plot_EW0estimates_wrapper(plotbasename,EW0dat,fluxratiodat,lineset,histaxes,Nhistbins,cdatvec,cdattype,
                                  ErrNsigma=ErrNsigma,point_text=point_text,vshiftmax=vshiftmax,performlinearfit=True,
                                  xlog=xlog,ylog=ylog,showlimits=showlimits,overwrite=overwrite,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_EW0estimates_wrapper(plotbasename,EWdat,fluxratiodat,EWset,histaxes,Nhistbins,cdatvec,cdattype,
                              photoionizationplotparam=None,point_text=None,showlimits=True,performlinearfit=False,
                              xlog=True,ylog=True,overwrite=False,ErrNsigma=1.0,vshiftmax=1e4,verbose=True):
    """
    Wrapper to define input data and excecute plot command

    """
    lineEW1,lineEW2,xrange,yrange,title = EWset

    if type(lineEW1) != str:
        str1  = lineEW1[0]
    else:
        str1  = lineEW1

    if type(lineEW2) != str:
        str2 = lineEW2[0]
    else:
        str2 = lineEW2

    plotname = plotbasename+'_EW0_'+str1+'vs'+str2+\
               '_ErrNsigma'+str(ErrNsigma).replace('.','p')+\
               '_vshiftLT'+str(vshiftmax).replace('.','p')+'.pdf'
    if not showlimits:
        plotname = plotname.replace('.pdf','_nolimits.pdf')

    if 'lya' in str1.lower():
        if ('f_'+str2 in fluxratiodat.dtype.names) & showlimits:
            goodent  = np.where(np.isfinite(EWdat['EW0_'+str2]) &  np.isfinite(lineEW1[2]) &
                                np.isfinite(lineEW1[2]) & (lineEW1[2] != 0) )[0]
        elif ('f_'+str2 in fluxratiodat.dtype.names) & (not showlimits):
            goodent  = np.where((np.abs(EWdat['EW0err_'+str2]) != 99.) & (lineEW1[2] != 0) & np.isfinite(lineEW1[2]) &
                                (EWdat['EW0err_'+str2] != 0) & (np.abs(lineEW1[3]) != 99.) &
                                np.isfinite(EWdat['EW0_'+str2]))[0]
        else:
            goodent  = []
    elif 'lya' in str2.lower():
        if ('f_'+str1 in fluxratiodat.dtype.names) & showlimits:
            goodent  = np.where(np.isfinite(EWdat['EW0_'+str1]) & np.isfinite(lineEW2[2]) &
                                np.isfinite(lineEW2[2]) & (lineEW2[2] != 0))[0]
        elif ('f_'+str1 in fluxratiodat.dtype.names) & (not showlimits):
            goodent  = np.where((np.abs(EWdat['EW0err_'+str1]) != 99.) & (lineEW2[2] != 0) & np.isfinite(lineEW2[2]) &
                                (EWdat['EW0err_'+str1] != 0) & (np.abs(lineEW2[3]) != 99.) &
                                np.isfinite(EWdat['EW0_'+str1]) )[0]
        else:
            goodent  = []
    else:
        if ('f_'+str1 in fluxratiodat.dtype.names) & ('f_'+str2 in fluxratiodat.dtype.names) & showlimits:
            goodent  = np.where(np.isfinite(EWdat['EW0_'+str1]) & np.isfinite(EWdat['EW0_'+str2]))[0]
        elif ('f_'+str1 in fluxratiodat.dtype.names) & ('f_'+str2 in fluxratiodat.dtype.names) & (not showlimits):
            goodent  = np.where((np.abs(EWdat['EW0err_'+str1]) != 99.) & (np.abs(EWdat['EW0err_'+str2]) != 99.) &
                                (EWdat['EW0err_'+str1] != 0) & (EWdat['EW0err_'+str2] != 0) &
                                np.isfinite(EWdat['EW0_'+str1]) & np.isfinite(EWdat['EW0_'+str2]))[0]
        else:
            goodent  = []

    if len(goodent) == 0:
        if verbose: print('\n - WARNING No good values found for the plot: \n           '+plotname.split('/')[-1]+'\n')
        goodent  = np.asarray([0,1])
        xvalues  = np.asarray([1e10]*2)
        xerr     = np.asarray([1.0]*2)
        xlabel   = 'EW$_0$('+str1+') [\AA]'
        yvalues  = np.asarray([1e10]*2)
        yerr     = np.asarray([1.0]*2)
        ylabel   = 'EW$_0$('+str2+') [\AA]'
        cdatvec  = np.asarray([0.0]*2)
        IDsALL   = np.asarray([0.0]*2)
    else:
        if 'lya' in str1.lower():
            xlabel   = lineEW1[1]
            xvalues  = lineEW1[2][goodent]
            xerr     = lineEW1[3][goodent]
        else:
            xlabel   = 'EW$_0$('+str1+') [\AA]'
            xvalues  = EWdat['EW0_'+str1][goodent]
            xerr     = EWdat['EW0err_'+str1][goodent]

        if 'lya' in str2.lower():
            ylabel   = lineEW2[1]
            yvalues  = lineEW2[2][goodent]
            yerr     = lineEW2[3][goodent]
        else:
            ylabel   = 'EW$_0$('+str2+') [\AA]'
            yvalues  = EWdat['EW0_'+str2][goodent]
            yerr     = EWdat['EW0err_'+str2][goodent]

        IDsALL = EWdat['id'][goodent]

    xerr[np.abs(xerr) != 99] = xerr[np.abs(xerr) != 99]*ErrNsigma
    yerr[np.abs(yerr) != 99] = yerr[np.abs(yerr) != 99]*ErrNsigma

    if point_text is not None:
        point_text = point_text[goodent]

    if ('lya' in str1.lower()) or ('lya' in str2.lower()):
        if 'EW' in xlabel:
            lines2show = 'onethreeten'

            if verbose:
                print(' - Ratio ('+ylabel+'/'+xlabel+')    [xlog='+str(xlog)+'; ylog='+str(ylog)+']')
                goodval = np.where((np.abs(xerr) != 99) & (np.abs(yerr) != 99) &
                                   np.isfinite(xvalues) & np.isfinite(yvalues) &
                                   (xvalues != 0) & (yvalues != 0))[0]
                yxratio = yvalues[goodval]/xvalues[goodval]
                print('   all       (mean,median,std) = '+str((np.mean(yxratio),np.median(yxratio),np.std(yxratio))))

                goodval = np.where((np.abs(xerr) != 99) & (np.abs(yerr) != 99) &
                                   np.isfinite(xvalues) & np.isfinite(yvalues) &
                                   (xvalues != 0) & (yvalues != 0) & (yvalues<xvalues))[0]
                yxratio = yvalues[goodval]/xvalues[goodval]
                print('   y<x       (mean,median,std) = '+str((np.mean(yxratio),np.median(yxratio),np.std(yxratio))))

                goodval = np.where((np.abs(xerr) != 99) & (np.abs(yerr) != 99) &
                                   np.isfinite(xvalues) & np.isfinite(yvalues) &
                                   (xvalues != 0) & (yvalues != 0) & (xvalues<yvalues))[0]
                yxratio = yvalues[goodval]/xvalues[goodval]
                print('   x<y       (mean,median,std) = '+str((np.mean(yxratio),np.median(yxratio),np.std(yxratio))))

        else:
            lines2show = None
    else:
        lines2show = 'onetoone'

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,
                                                   'dummydat',linetype=lines2show,title=title,ids=IDsALL,
                                                   ylog=ylog,xlog=xlog,yrange=yrange,xrange=xrange,
                                                   colortype=cdattype,colorcode=True,cdatvec=cdatvec[goodent],
                                                   point_text=point_text,photoionizationplotparam=photoionizationplotparam,
                                                   histaxes=histaxes,Nbins=Nhistbins,
                                                   overwrite=overwrite,verbose=verbose)
    if performlinearfit & (len(goodent) > 0) & (len(xvalues[np.isfinite(xvalues)]) > 2) & (len(yvalues[np.isfinite(yvalues)]) > 2):
        plotname   = plotname.replace('.pdf','_ODRfit2data.pdf')

        if ylog:
            yval_fit = np.log10(yvalues)
            yerr_fit = 0.434 * yerr/yvalues  # https://faculty.washington.edu/stuve/log_error.pdf
        else:
            yval_fit = yvalues
            yerr_fit = yerr

        if xlog:
            xval_fit = np.log10(xvalues)
            xerr_fit = 0.434 * xerr/xvalues  # https://faculty.washington.edu/stuve/log_error.pdf
        else:
            xval_fit = xvalues
            xerr_fit = xerr

        fitres,rp,pp,rs,ps = kbs.fit_function_to_data_with_errors_on_both_axes(xval_fit,yval_fit,xerr_fit,yerr_fit,
                                                                               initguess = [1.,1.],verbose=verbose,
                                                                               fitfunction='linear',plotresults=plotname,
                                                                               show3sigmainterval=True,
                                                                               returnCorrelationCoeffs=True)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def linenameUVES2NEOGAL(uvesname):
    """
    Translateror between the UVES names of emission lines and the NEOGAL emission line names

    """
    translatedic = {}
    # - - - - - - - UVES -> NEOGAL - - - - - - -
    translatedic['NV']        = 'NV1240'
    translatedic['CIV1']      = 'CIV1549'
    translatedic['CIV']       = 'CIV1550'# Not a NEOGAL column but calculated in uves.add_photoionization_models_to_lineratioplot()
    translatedic['CIV2']      = 'CIV1551'
    translatedic['CIII']      = 'CIII1908'
    translatedic['CIII1']     = 'CIII1907'
    translatedic['CIII2']     = 'CIII1910'
    translatedic['HeII']      = 'HeII1640'
    translatedic['OIII1']     = 'OIII1661'
    translatedic['OIII']      = 'OIII1663'
    translatedic['OIII2']     = 'OIII1666'
    translatedic['SiIII1']    = 'SiIII1883'
    translatedic['SiIII']     = 'SiIII1888'
    translatedic['SiIII2']    = 'SiIII1892'# Not a NEOGAL column but calculated in uves.add_photoionization_models_to_lineratioplot()

    return translatedic[uvesname]
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def linewavesUVES(uvesname):
    """
    Return line wavelength for a UVES line

    """
    linewavedic = {}
    # - - - - - - - UVES -> NEOGAL - - - - - - -
    linewavedic['nv1']       = 1238.821
    linewavedic['nv']        = 1240
    linewavedic['nv2']       = 1242.804
    linewavedic['civ1']      = 1548.195
    linewavedic['civ']       = 1550.00
    linewavedic['civ2']      = 1550.770
    linewavedic['ciii1']     = 1906.68
    linewavedic['ciii']      = 1908.00
    linewavedic['ciii2']     = 1908.73
    linewavedic['heii']      = 1640.420
    linewavedic['oiii1']     = 1660.809
    linewavedic['oiii']      = 1663.00
    linewavedic['oiii2']     = 1666.150
    linewavedic['siiii1']    = 1882.71
    linewavedic['siiii']     = 1888.00
    linewavedic['siiii2']    = 1892.03
    linewavedic['mgii1']     = 2795.528
    linewavedic['mgii']      = 2799.00
    linewavedic['mgii2']     = 2802.705

    return linewavedic[uvesname.lower()]
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def calculatelineratios(outputfile='./fluxratioresults.txt', S2Nmaxrange=[5.0,100.0], zspecrange=[0.0,10.0],
                        voffsetrange=[-1500.0,1500.0], onesigmalimit=100.0, plotdir=None, verbose=True):
    """
    Function to calculate the flux and line ratios for a selection of objects with FELIS template matches
    satisfying a set of criteria based on the FELIS output pickle.

    --- INPUT ---
    outputfile       The ascii file to write results to
    S2Nmaxrange      Range of S/N template matches to unclude
    zspecrange       The range of redshifts for objects to include
    voffsetrange     Estimated velocity offset of line to include
    onesigmalimit    1 sigma flux limit to use for ratio limits (in units of spectra 10.0)
    plotdir          Provide path to output directory, to plot selected tempalte matches.
    verbose          Toggle verbosity

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    uves.calculatelineratios(outputfile='./fluxratioresults.txt', S2Nmaxrange=[5.0,100.0], zspecrange=[0.0,10.0], voffsetrange=[-1500.0,1500.0], onesigmalimit=100.0)

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - The estimated line flux ratios will be based on the following FELIS outputs ')
    picklepath  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'
    pickleCIII  = picklepath+'MUSEWideLAEs_CCresultsdateofrun_CIII_all575specX180914templates.pkl'
    pickleCIV   = picklepath+'MUSEWideLAEs_CCresultsdateofrun_CIV_all575specX180914templates.pkl'
    pickleHeII  = picklepath+'MUSEWideLAEs_CCresultsdateofrun_HEII_all575specX180914templates.pkl'
    pickleOIII  = picklepath+'MUSEWideLAEs_CCresultsdateofrun_OIII_all575specX180914templates.pkl'
    pickleNV    = picklepath+'MUSEWideLAEs_CCresultsdateofrun_NV_all575specX180914templates.pkl'
    zspecISzLya = False

    if verbose: print('   '+pickleCIII)
    if verbose: print('   '+pickleCIV)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Performing selection on FELIS output ')

    goodkeysCIII    = felis.selection_from_picklefile(pickleCIII, S2Nmaxrange=S2Nmaxrange, zspecrange=zspecrange,
                                                    voffsetrange=voffsetrange, zspecISzLya=zspecISzLya)
    goodkeysCIV     = felis.selection_from_picklefile(pickleCIV, S2Nmaxrange=S2Nmaxrange, zspecrange=zspecrange,
                                                    voffsetrange=voffsetrange, zspecISzLya=zspecISzLya)
    goodkeysHeII    = felis.selection_from_picklefile(pickleHeII, S2Nmaxrange=S2Nmaxrange, zspecrange=zspecrange,
                                                    voffsetrange=voffsetrange, zspecISzLya=zspecISzLya)
    goodkeysOIII    = felis.selection_from_picklefile(pickleOIII, S2Nmaxrange=S2Nmaxrange, zspecrange=zspecrange,
                                                    voffsetrange=voffsetrange, zspecISzLya=zspecISzLya)
    goodkeysNV      = felis.selection_from_picklefile(pickleNV, S2Nmaxrange=S2Nmaxrange, zspecrange=zspecrange,
                                                    voffsetrange=voffsetrange, zspecISzLya=zspecISzLya)

    goodkeys_unique = np.unique(goodkeysCIV+goodkeysCIII+goodkeysHeII+goodkeysOIII+goodkeysNV)

    loaddicCIII     = felis.load_picklefile(pickleCIII)
    loaddicCIV      = felis.load_picklefile(pickleCIV)
    loaddicHeII     = felis.load_picklefile(pickleHeII)
    loaddicOIII     = felis.load_picklefile(pickleOIII)
    loaddicNV       = felis.load_picklefile(pickleNV)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if os.path.isfile(outputfile):
        sys.exit('Ouptut file '+outputfile+' already exists')

    if verbose: print(' - Initializing the output file '+outputfile)
    fout = open(outputfile,'w')
    fout.write('# Flux and line ratios estimated based on FELIS template match results.\n')
    fout.write('# \n')
    fout.write('# Each template match was narrowed down to only contain matches with \n')
    fout.write('#   - z_spec         :  ['+str(zspecrange[0])+','+str(zspecrange[1])+']\n')
    fout.write('#   - max(S/N)       :  ['+str(S2Nmaxrange[0])+','+str(S2Nmaxrange[1])+']\n')
    fout.write('#   - voffset[km/s]  :  ['+str(voffsetrange[0])+','+str(voffsetrange[1])+']\n')
    fout.write('# \n')
    fout.write('# The selection was performed on the following pickle files outputted by FELIS \n')
    fout.write('#   - CIII: '+pickleCIII+'\n')
    fout.write('#     ('+str(len(goodkeysCIII))+' satisfying the cuts)\n')
    fout.write('#   - CIV: '+pickleCIV+'\n')
    fout.write('#     ('+str(len(goodkeysCIV))+' satisfying the cuts)\n')
    fout.write('#   - HeII: '+pickleHeII+'\n')
    fout.write('#     ('+str(len(goodkeysHeII))+' satisfying the cuts)\n')
    fout.write('#   - OIII: '+pickleOIII+'\n')
    fout.write('#     ('+str(len(goodkeysOIII))+' satisfying the cuts)\n')
    fout.write('#   - NV: '+pickleNV+'\n')
    fout.write('#     ('+str(len(goodkeysNV))+' satisfying the cuts)\n')
    fout.write('# \n')
    fout.write('# In total '+str(len(goodkeys_unique))+' template matches satisfied the cuts and are included in this file \n')
    fout.write('# \n')
    fout.write('# Upper and lower limits are given as negative values with uncertainty of +99 or -99, respectively. \n')
    fout.write('# For the limits a 1sigma limit of '+str(onesigmalimit)+' was used \n')
    fout.write('# \n')
    fout.write('# This file contains the following columns:\n')

    fluxratiodic = collections.OrderedDict()
    fluxratiodic['vshift_AV18']          = 999
    fluxratiodic['vshift_ciii1908']      = 999
    fluxratiodic['vshift_civ1550']       = 999
    fluxratiodic['vshift_heii1640']      = 999
    fluxratiodic['vshift_oiii1663']      = 999
    fluxratiodic['vshift_nv1241']        = 999

    fluxratiodic['f_ciii1908']           = 999
    fluxratiodic['ferr_ciii1908']        = 999
    fluxratiodic['S2N_ciii1908']         = 999
    fluxratiodic['sigma_ciii1908']       = 999

    fluxratiodic['f_civ1550']            = 999
    fluxratiodic['ferr_civ1550']         = 999
    fluxratiodic['S2N_civ1550']          = 999
    fluxratiodic['sigma_civ1550']        = 999

    fluxratiodic['f_heii1640']           = 999
    fluxratiodic['ferr_heii1640']        = 999
    fluxratiodic['S2N_heii1640']         = 999
    fluxratiodic['sigma_heii1640']       = 999

    fluxratiodic['f_oiii1663']           = 999
    fluxratiodic['ferr_oiii1663']        = 999
    fluxratiodic['S2N_oiii1663']         = 999
    fluxratiodic['sigma_oiii1663']       = 999

    fluxratiodic['f_nv1241']             = 999
    fluxratiodic['ferr_nv1241']          = 999
    fluxratiodic['S2N_nv1241']           = 999
    fluxratiodic['sigma_nv1241']         = 999

    fluxratiodic['ciii1907ciii1909']     = 999
    fluxratiodic['ciii1907ciii1909err']  = 999
    fluxratiodic['civ1549civ1551']       = 999
    fluxratiodic['civ1549civ1551err']    = 999
    fluxratiodic['oiii1661oiii1666']     = 999
    fluxratiodic['oiii1661oiii1666err']  = 999
    fluxratiodic['nv1239nv1243']         = 999
    fluxratiodic['nv1239nv1243err']      = 999

    fluxratiodic['civ1550ciii1908']      = 999
    fluxratiodic['civ1550ciii1908err']   = 999
    fluxratiodic['ciii1908heii1640']     = 999
    fluxratiodic['ciii1908heii1640err']  = 999
    fluxratiodic['civ1550heii1640']      = 999
    fluxratiodic['civ1550heii1640err']   = 999
    fluxratiodic['oiii1663heii1640']     = 999
    fluxratiodic['oiii1663heii1640err']  = 999
    fluxratiodic['nv1241heii1640']       = 999
    fluxratiodic['nv1241heii1640err']    = 999
    fluxratiodic['nv1241civ1550']        = 999
    fluxratiodic['nv1241civ1550err']     = 999

    hdrstr = '# id_musewide '
    for col in fluxratiodic.keys():
        hdrstr = hdrstr+col+' '
    hdrstr = hdrstr+'\n'
    fout.write(hdrstr)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Calculating line ratios based on FELIS template matches')
    for goodkey in goodkeys_unique:
        template_ciii1908, vshift_V18_ciii1908, vshift_ciii1908, f_ciii1908, \
        ferr_ciii1908, S2Nmax_ciii1908, Ngoodent_ciii1908, chi2_ciii1908, zspec, zS2Nmax = \
            felis.getresult4maxS2N(loaddicCIII,goodkey)
        template_civ1550, vshift_V18_civ1550, vshift_civ1550, f_civ1550, \
        ferr_civ1550, S2Nmax_civ1550, Ngoodent_civ1550, chi2_civ1550, zspec, zS2Nmax = \
            felis.getresult4maxS2N(loaddicCIV,goodkey)
        template_heii1640, vshift_V18_heii1640, vshift_heii1640, f_heii1640, \
        ferr_heii1640, S2Nmax_heii1640, Ngoodent_heii1640, chi2_heii1640, zspec, zS2Nmax = \
            felis.getresult4maxS2N(loaddicHeII,goodkey)
        template_oiii1663, vshift_V18_oiii1663, vshift_oiii1663, f_oiii1663, \
        ferr_oiii1663, S2Nmax_oiii1663, Ngoodent_oiii1663, chi2_oiii1663, zspec, zS2Nmax = \
            felis.getresult4maxS2N(loaddicOIII,goodkey)
        template_nv1241, vshift_V18_nv1241, vshift_nv1241, f_nv1241, \
        ferr_nv1241, S2Nmax_nv1241, Ngoodent_nv1241, chi2_nv1241, zspec, zS2Nmax = \
            felis.getresult4maxS2N(loaddicNV,goodkey)

        fluxratiodic['vshift_AV18']          = np.max(np.array([vshift_V18_ciii1908, vshift_V18_civ1550,
                                                                vshift_V18_heii1640, vshift_V18_oiii1663,
                                                                vshift_V18_nv1241])) # max to handle -99 values
        fluxratiodic['vshift_ciii1908']      = vshift_ciii1908
        fluxratiodic['vshift_civ1550']       = vshift_civ1550
        fluxratiodic['vshift_heii1640']      = vshift_heii1640
        fluxratiodic['vshift_oiii1663']      = vshift_oiii1663
        fluxratiodic['vshift_nv1241']        = vshift_nv1241


        if template_ciii1908 is 'None':
            f_ciii1908    = -1.0 * onesigmalimit
            ferr_ciii1908 = +99
        if template_civ1550 is 'None':
            f_civ1550     = -1.0 * onesigmalimit
            ferr_civ1550  = +99
        if template_heii1640 is 'None':
            f_heii1640    = -1.0 * onesigmalimit
            ferr_heii1640 = +99
        if template_oiii1663 is 'None':
            f_oiii1663    = -1.0 * onesigmalimit
            ferr_oiii1663 = +99
        if template_nv1241 is 'None':
            f_nv1241      = -1.0 * onesigmalimit
            ferr_nv1241   = +99

        fluxratiodic['f_ciii1908']           = f_ciii1908
        fluxratiodic['ferr_ciii1908']        = ferr_ciii1908
        # fluxratiodic['S2N_ciii1908_calc']    = f_ciii1908/ferr_ciii1908
        fluxratiodic['S2N_ciii1908']         = S2Nmax_ciii1908

        fluxratiodic['f_civ1550']            = f_civ1550
        fluxratiodic['ferr_civ1550']         = ferr_civ1550
        # fluxratiodic['S2N_civ1550_calc']     = f_civ1550/ferr_civ1550
        fluxratiodic['S2N_civ1550']          = S2Nmax_civ1550

        fluxratiodic['f_heii1640']           = f_heii1640
        fluxratiodic['ferr_heii1640']        = ferr_heii1640
        # fluxratiodic['S2N_heii1640_calc']    = f_heii1640/ferr_heii1640
        fluxratiodic['S2N_heii1640']         = S2Nmax_heii1640

        fluxratiodic['f_oiii1663']           = f_oiii1663
        fluxratiodic['ferr_oiii1663']        = ferr_oiii1663
        # fluxratiodic['S2N_oiii1663_calc']    = f_oiii1663/ferr_oiii1663
        fluxratiodic['S2N_oiii1663']         = S2Nmax_oiii1663

        fluxratiodic['f_nv1241']           = f_nv1241
        fluxratiodic['ferr_nv1241']        = ferr_nv1241
        # fluxratiodic['S2N_nv1241_calc']    = f_nv1241/ferr_nv1241
        fluxratiodic['S2N_nv1241']         = S2Nmax_nv1241

        if template_ciii1908 is 'None':
            fluxratiodic['ciii1907ciii1909']    = np.nan
            fluxratiodic['ciii1907ciii1909err'] = np.nan
            fluxratiodic['sigma_ciii1908']      = np.nan
        else:
            fluxratiodic['ciii1907ciii1909']    = float(template_ciii1908.split('fluxratio_')[-1].split('.')[0].replace('p','.'))
            fluxratiodic['ciii1907ciii1909err'] = ferr_ciii1908
            temp_sigma_ciii1908                 = float(template_ciii1908.split('sig_')[-1].split('_')[0].replace('p','.'))
            fluxratiodic['sigma_ciii1908']      = 299792.458 * 2.354 * temp_sigma_ciii1908 / 1908.0

        if template_civ1550 is 'None':
            fluxratiodic['civ1549civ1551']     = np.nan
            fluxratiodic['civ1549civ1551err']  = np.nan
            fluxratiodic['sigma_civ1550']      = np.nan
        else:
            fluxratiodic['civ1549civ1551']     = float(template_civ1550.split('fluxratio_')[-1].split('.')[0].replace('p','.'))
            fluxratiodic['civ1549civ1551err']  = ferr_civ1550
            temp_sigma_civ1550                 = float(template_civ1550.split('sig_')[-1].split('_')[0].replace('p','.'))
            fluxratiodic['sigma_civ1550']      = 299792.458 * 2.354 * temp_sigma_civ1550 / 1550.0

        if template_oiii1663 is 'None':
            fluxratiodic['oiii1661oiii1666']    = np.nan
            fluxratiodic['oiii1661oiii1666err'] = np.nan
            fluxratiodic['sigma_oiii1663']      = np.nan
        else:
            fluxratiodic['oiii1661oiii1666']    = float(template_oiii1663.split('fluxratio_')[-1].split('.')[0].replace('p','.'))
            fluxratiodic['oiii1661oiii1666err'] = ferr_oiii1663
            temp_sigma_oiii1663                 = float(template_oiii1663.split('sig_')[-1].split('_')[0].replace('p','.'))
            fluxratiodic['sigma_oiii1663']      = 299792.458 * 2.354 * temp_sigma_oiii1663 / 1550.0

        if template_nv1241 is 'None':
            fluxratiodic['nv1239nv1243']     = np.nan
            fluxratiodic['nv1239nv1243err']  = np.nan
            fluxratiodic['sigma_nv1241']     = np.nan
        else:
            fluxratiodic['nv1239nv1243']     = float(template_nv1241.split('fluxratio_')[-1].split('.')[0].replace('p','.'))
            fluxratiodic['nv1239nv1243err']  = ferr_nv1241
            temp_sigma_nv1241                = float(template_nv1241.split('sig_')[-1].split('_')[0].replace('p','.'))
            fluxratiodic['sigma_nv1241']     = 299792.458 * 2.354 * temp_sigma_nv1241 / 1550.0

        if template_heii1640 is 'None':
            fluxratiodic['sigma_heii1640']       = np.nan
        else:
            temp_sigma_heii1640                  = float(template_heii1640.split('sig_')[-1].split('.fi')[0].replace('p','.'))
            fluxratiodic['sigma_heii1640']       = 299792.458 * 2.354 * temp_sigma_heii1640 / 1640.0

        fluxratiodic['civ1550ciii1908'], fluxratiodic['civ1550ciii1908err'] = \
            lce.set_ratios(template_civ1550,template_ciii1908,f_civ1550,ferr_civ1550,f_ciii1908,ferr_ciii1908)

        fluxratiodic['ciii1908heii1640'], fluxratiodic['ciii1908heii1640err'] = \
            lce.set_ratios(template_ciii1908,template_heii1640,f_ciii1908,ferr_ciii1908,f_heii1640,ferr_heii1640)

        fluxratiodic['civ1550heii1640'], fluxratiodic['civ1550heii1640err'] = \
            lce.set_ratios(template_civ1550,template_heii1640,f_civ1550,ferr_civ1550,f_heii1640,ferr_heii1640)

        fluxratiodic['oiii1663heii1640'], fluxratiodic['oiii1663heii1640err'] = \
            lce.set_ratios(template_oiii1663,template_heii1640,f_oiii1663,ferr_oiii1663,f_heii1640,ferr_heii1640)

        fluxratiodic['nv1241heii1640'], fluxratiodic['nv1241heii1640err'] = \
            lce.set_ratios(template_nv1241,template_heii1640,f_nv1241,ferr_nv1241,f_heii1640,ferr_heii1640)

        fluxratiodic['nv1241civ1550'], fluxratiodic['nv1241civ1550err'] = \
            lce.set_ratios(template_nv1241,template_civ1550,f_nv1241,ferr_nv1241,f_civ1550,ferr_civ1550)

        outstr = ' '
        outstr = outstr+goodkey.split('-')[-1].split('.')[0][1:]+' '
        for col in fluxratiodic.keys():
            outstr = outstr+str("%12.2f" % fluxratiodic[col])+' '
        outstr = outstr+'\n'
        fout.write(outstr)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Wrote results to '+outputfile)
    fout.close()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plotdir is not None:
        if verbose: print(' - plotobjects=True so plotting selections to the directory:\n    '+plotdir)
        goodkeysall = [goodkeysCIII,goodkeysCIV,goodkeysHeII,goodkeysOIII,goodkeysNV]
        picklefiles = [pickleCIII,pickleCIV,pickleHeII,pickleOIII,pickleNV]
        lines = ['CIII1908','CIV1549','HeII1640','OIII1663','NV1241']

        for gg, goodkeys in enumerate(goodkeysall):
            plotnames = []
            plotnameinput = []
            for key in goodkeys:
                fname = key.split('/')[-1]
                plotnameinput.append(key.replace('.fits','_maxS2Ntemplatematch_'+lines[gg]+'.pdf'))
                plotnames.append(plotdir+fname.replace('.fits','_maxS2Ntemplatematch_'+lines[gg]+'.pdf'))

            felis.plot_picklefilecontent(goodkeys, picklefiles[gg], plotdir=plotdir, plotnames=plotnameinput,
                                         showspecerr=False, zspecISzLya=zspecISzLya)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def TDOSE_sourcecat_from_infofile(infofile,outputdir,minRaper=0.5,minCutwidth=4.5,goodmatchsep=0.25,overwrite=False,verbose=True):
    """
    Generate a TDOSE source catatalog (and catalog with intitial guesses for Gaussian modeling) from
    an infofile generated with uves.build_LAEfitstable().
    Based on the MUSEWideUtilities.TDOSE_sourcecat_from_*() scripts

    --- INPUT ---
    infofile          infofile from uves.build_LAEfitstable() to based source catalogs on
    outputdir         location of output
    goodmatchsep      Any match to the photometric catalog below goodmatchsep is replaced by MW object
    overwrite         Overwrite output files if they exists?
    verbose           Toggle verbosity

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo_UVemitters_3timesUDFcats.fits'
    outputdir     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_sourcecatalogs_100fields/'
    uves.TDOSE_sourcecat_from_infofile(infofile,outputdir,minRaper=0.25,minCutwidth=4.5,goodmatchsep=0.25,overwrite=True,verbose=True)

    """
    if verbose: print(' - Generating TDOSE source catalog from the infofile:\n'+infofile+' '
                      '   Restricting source to FoV of individual reference images\n'
                      '   Output will be saved to '+outputdir)

    infodat              = afits.open(infofile)[1].data
    ids_all              = infodat['id']
    ras_all              = infodat['ra']
    decs_all             = infodat['dec']
    aimage_arr_inarcsec  = np.array([])

    skelcosmos     = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat.FITS'
    skelcosmosdat  = afits.open(skelcosmos)[1].data

    skelcdfs       = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
    skelcdfsdat    = afits.open(skelcdfs)[1].data

    rafelskicat    = '/Users/kschmidt/work/catalogs/rafelski/uvudf_rafelski_2015.fits'
    rafelskidat    = afits.open(rafelskicat)[1].data

    imgpath        = '/Users/kschmidt/work/images_MAST/MUSEWidePointings/'
    refimages      = glob.glob(imgpath+'*_814w_*cdfs*.fits')
    refimages      = refimages + glob.glob(imgpath+'*_814w_*cosmos*.fits')
    refimages      = refimages + glob.glob(imgpath+'*_160w_*hudf09*.fits')
    refimages      = refimages + glob.glob(imgpath+'*_775w_*udf-0*rot.fits')
    refimages      = refimages + ['/Users/kschmidt/work/images_MAST/MUSEWidePointings/acs_775w_udf-10_cut.fits']

    for refimage in refimages:
        if '_wht_' in refimage:
            continue
        else:
            if 'udf-' in refimage:
                outputnamebase = outputdir+refimage.split('/')[-1].\
                    replace('acs_775w','tdose_sourcecat_MWuves_acs_775w').replace('.fits','')
            elif 'hudf09' in refimage:
                outputnamebase = outputdir+refimage.split('/')[-1].\
                    replace('wfc3_160w','tdose_sourcecat_MWuves_wfc3_160w').replace('.fits','')
            else:
                outputnamebase = outputdir+refimage.split('/')[-1].\
                    replace('acs_814w','tdose_sourcecat_MWuves_acs_814w').replace('.fits','')

        refimgdata       = afits.open(refimage)[0].data
        refimghdr        = afits.open(refimage)[0].header

        if '-cosmos-' in refimage:
            skeltondat       = skelcosmosdat
        else:
            skeltondat       = skelcdfsdat

        if '_udf-' in refimage:
            match_sep          = infodat['sep_rafelski']
            match_id           = infodat['id_rafelski']
            phot_id_all        = rafelskidat['ID']
            phot_ra_all        = rafelskidat['RA']
            phot_dec_all       = rafelskidat['DEC']
            phot_flux_all      = rafelskidat['FLUX_ISO_F775W']
            phot_theta_all     = rafelskidat['THETA']
            phot_a_image_all   = np.sqrt(rafelskidat['AREAF']/np.pi/(1-rafelskidat['ELLIPTICITY']))
            phot_b_image_all   = (1.0-rafelskidat['ELLIPTICITY'])*phot_a_image_all
            arcsecPerPix_phot  = 0.03
            fivesigma_maglimit = 29.5  # f775w 5sigma limit
            onesigma_flux     = 10**((8.90-fivesigma_maglimit)/2.5) / 5.0 * 1e6 # flux in muJy
        else:
            match_sep          = infodat['sep_skelton']
            match_id           = infodat['id_skelton']
            phot_id_all        = skeltondat['id']
            phot_ra_all        = skeltondat['ra']
            phot_dec_all       = skeltondat['dec']
            phot_flux_all      = skeltondat['f_f160w']
            phot_theta_all     = skeltondat['theta_j2000']
            phot_a_image_all   = skeltondat['a_image']
            phot_b_image_all   = skeltondat['b_image']
            arcsecPerPix_phot  = 0.06
            if '-cosmos-' in refimage:
                fivesigma_maglimit = 25.8  # f160w 5sigma limit COSMOS
            else:
                fivesigma_maglimit = 26.4  # f160w 5sigma limit GOODS-S
            onesigma_flux     = 10**((25-fivesigma_maglimit)/2.5) / 5.

        fluxf_all   = []
        theta_all   = []
        a_image_all = []
        b_image_all = []

        matchedids = []
        for mm, mid in enumerate(match_id):
            MWid = ids_all[mm]
            firstdigit = int(str(MWid)[0])
            if '-cosmos-' in refimage:
                if firstdigit in [2]:
                    includeobj = True
                    if match_sep[mm] < goodmatchsep: matchedids.append(mid)
                else:
                    includeobj = False
            elif '_udf-0' in refimage:
                if firstdigit in [6]: # Include the mosaic full-deth ids; not the MWmock ids (firstdigit=5).
                    includeobj = True
                    if match_sep[mm] < goodmatchsep: matchedids.append(mid)
                else:
                    includeobj = False
            elif '_udf-10' in refimage:
                if firstdigit in [7]:
                    includeobj = True
                    if match_sep[mm] < goodmatchsep: matchedids.append(mid)
                else:
                    includeobj = False
            else:
                if firstdigit in [1,3,4]:
                    includeobj = True
                    if match_sep[mm] < goodmatchsep: matchedids.append(mid)
                else:
                    includeobj = False

            if includeobj:
                if  (match_sep[mm] < goodmatchsep) & (match_sep[mm] > 0.0):
                    phot_ent = np.where(phot_id_all== mid)[0]
                    fluxf_all.append(phot_flux_all[phot_ent][0])
                    theta_all.append(phot_theta_all[phot_ent][0])
                    a_image_all.append(phot_a_image_all[phot_ent][0])
                    b_image_all.append(phot_b_image_all[phot_ent][0])
                else:
                    fluxf_all.append(onesigma_flux)
                    theta_all.append(0.0)
                    a_image_all.append(2.0)
                    b_image_all.append(2.0)

                if len(aimage_arr_inarcsec) == 0:
                    aimage_arr_inarcsec = np.array([str(MWid),str(np.asarray(a_image_all)[-1]*arcsecPerPix_phot)])
                else:
                    aimage_arr_inarcsec = np.vstack([aimage_arr_inarcsec,
                                                     np.array([str(MWid),str(np.asarray(a_image_all)[-1]*arcsecPerPix_phot)])])
            else:
                fluxf_all.append(-99)
                theta_all.append(-99)
                a_image_all.append(-99)
                b_image_all.append(-99)

        # if '_udf-10' in refimage: pdb.set_trace()
        fluxf_all    = np.asarray(fluxf_all)
        theta_all    = np.asarray(theta_all)
        a_image_all  = np.asarray(a_image_all)
        b_image_all  = np.asarray(b_image_all)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print(' - Finding objects within reference image field-of-view')
        if verbose: print('   (Estimating pixel positions using wcs info from header)')
        striphdr   = tu.strip_header(refimghdr,verbose=verbose,delkeys=['COMMENT','HISTORY','','A_ORDER','B_ORDER'])
        wcs_in     = wcs.WCS(striphdr)
        arcsecPerPix_refimg = wcs.utils.proj_plane_pixel_scales(wcs_in)*60*60

        skycoord   = SkyCoord(ras_all, decs_all, frame='fk5', unit='deg')
        pixcoord   = wcs.utils.skycoord_to_pixel(skycoord,wcs_in,origin=1)
        xpos       = pixcoord[0]
        ypos       = pixcoord[1]
        goodent    = np.where((xpos < refimghdr['NAXIS1']) & (xpos > 0) &
                              (ypos < refimghdr['NAXIS2']) & (ypos > 0) &
                              (a_image_all > 0.0))[0]

        if verbose: print('   (Make sure no 0s exist in a 6x6 pixel region around position, '
                          'i.e., ignoring 0-edges of ref images.)')
        Ngoodinit  = len(goodent)
        intxpos    = np.round(xpos).astype(int)
        intypos    = np.round(ypos).astype(int)
        for ent in goodent:
            if (refimgdata[np.max([intypos[ent]-5,0]):np.min([intypos[ent]+5,refimghdr['NAXIS2']-1]),
                np.max([intxpos[ent]-5,0]):np.min([intxpos[ent]+5,refimghdr['NAXIS1']-1])] == 0).any():
                goodent[np.where(goodent == ent)[0]] = -99
        Nedge     = len(np.where(goodent == -99)[0])
        goodent   = goodent[np.where(goodent != -99)[0]]
        if verbose: print('   (Ended up removing '+str(Nedge)+'/'+str(Ngoodinit)+
                          ' objects that fall in the edge region but are within the image FoV)')

        ids             = ids_all[goodent]
        ras             = ras_all[goodent]
        decs            = decs_all[goodent]
        x_image         = xpos[goodent]
        y_image         = ypos[goodent]
        fluxscale       = fluxf_all[goodent]
        a_image         = a_image_all[goodent] * arcsecPerPix_phot / np.mean(arcsecPerPix_refimg)
        b_image         = b_image_all[goodent] * arcsecPerPix_phot / np.mean(arcsecPerPix_refimg)
        theta           = theta_all[goodent]

        MWidsinfield    = ids
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print(' - Add photometric objects around MUSE-Wide objects of interest')

        skycoord   = SkyCoord(phot_ra_all, phot_dec_all, frame='fk5', unit='deg')
        pixcoord   = wcs.utils.skycoord_to_pixel(skycoord,wcs_in,origin=1)
        xpos       = pixcoord[0]
        ypos       = pixcoord[1]
        goodent    = np.where((xpos < refimghdr['NAXIS1']) & (xpos > 0) &
                              (ypos < refimghdr['NAXIS2']) & (ypos > 0) &
                              (phot_a_image_all > 0.0))[0]

        if verbose: print('   (Make sure no 0s exist in a 6x6 pixel region around position, '
                          'i.e., ignoring 0-edges of ref images.)')
        Ngoodinit  = len(goodent)
        intxpos    = np.round(xpos).astype(int)
        intypos    = np.round(ypos).astype(int)
        for ent in goodent:
            if (refimgdata[np.max([intypos[ent]-5,0]):np.min([intypos[ent]+5,refimghdr['NAXIS2']-1]),
                np.max([intxpos[ent]-5,0]):np.min([intxpos[ent]+5,refimghdr['NAXIS1']-1])] == 0).any() or \
                    (phot_id_all[ent] in matchedids):
                goodent[np.where(goodent == ent)[0]] = -99
        Nedge     = len(np.where(goodent == -99)[0])
        goodent   = goodent[np.where(goodent != -99)[0]]
        if verbose: print('   (Ended up removing '+str(Nedge)+'/'+str(Ngoodinit)+
                          ' objects that fall in the edge region but are within the image FoV)')

        ids        = np.append(ids,phot_id_all[goodent])
        ras        = np.append(ras,phot_ra_all[goodent])
        decs       = np.append(decs,phot_dec_all[goodent])
        x_image    = np.append(x_image,xpos[goodent])
        y_image    = np.append(y_image,ypos[goodent])
        fluxscale  = np.append(fluxscale,phot_flux_all[goodent])
        a_image    = np.append(a_image,phot_a_image_all[goodent]* arcsecPerPix_phot / np.mean(arcsecPerPix_refimg))
        b_image    = np.append(b_image,phot_b_image_all[goodent]* arcsecPerPix_phot / np.mean(arcsecPerPix_refimg))
        theta      = np.append(theta,phot_theta_all[goodent])
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        outtxt          = outputnamebase+'.txt'
        if (overwrite == False) & os.path.isfile(outtxt):
            sys.exit('Output ('+outtxt+') already exists and clobber=False')
        else:
            if verbose: print(' - Will save source catalog to '+outtxt+' (overwriting any existing files)')
            fout = open(outtxt,'w')
            fout.write('# TDOSE Source catalog generated with '
                       'uvEmissionlineSearch.TDOSE_sourcecat_from_infofile() on '+kbs.DandTstr2()+'\n')
            fout.write('# see objects with ds9 '+refimage+' -regions '+outputnamebase+'.reg \n')
            fout.write('# parent_id id ra dec x_image y_image fluxscale a_image b_image theta \n')
            for ii, id in enumerate(ids):
                fout.write(str(ids[ii])+' '+str(ids[ii])+' '+str(ras[ii])+' '+str(decs[ii])+' '+
                           str(x_image[ii])+' '+str(y_image[ii])+' '+str(fluxscale[ii])+' '+
                           str(a_image[ii])+' '+str(b_image[ii])+' '+str(theta[ii])+'  \n')

            fout.close()

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            regionfile = outtxt.replace('.txt','.reg')
            if verbose: print(' - Storing DS9 region file to '+regionfile)
            idsstr     = [str(id) for id in ids]
            tu.create_simpleDS9region(regionfile,ras,decs,color='cyan',
                                      circlesize=2.0 * a_image * np.mean(arcsecPerPix_refimg),
                                      textlist=idsstr,clobber=overwrite)
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            outnamefits = outtxt.replace('.txt','.fits')
            if verbose: print(' - Saving fits version of source catalog to '+outnamefits)
            fitsfmt       = ['D']*10
            sourcecatfits = tu.ascii2fits(outtxt,asciinames=True,skip_header=2,fitsformat=fitsfmt,verbose=verbose)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            extidcat      = outtxt.replace('.txt','_objects2extract.txt')
            if verbose: print(' - Will save ids of objects to extract in '+extidcat+' (overwriting any existing files)')
            fout = open(extidcat,'w')
            fout.write('# id \n')
            for mm, mwid in enumerate(MWidsinfield):
                fout.write(str(MWidsinfield[mm])+'  \n')
            fout.close()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Creating cutoutsizes.txt and aperturesizes.txt for the TDOSE extractions along the way ')
    uid, uindex = np.unique(aimage_arr_inarcsec[:,0].astype(int), return_index=True)
    aimage_arr_inarcsec = aimage_arr_inarcsec[uindex,:]

    # ===== CTUOUT SIZES =====
    cut_outtxt  = outputdir+'cutoutsizes_MWuves_4xA_IMAGE.txt'
    fout        = open(cut_outtxt,'w')
    fout.write('# Cutout sizes of 4*A_IMAGE on each side (corresponding to twice the aperture diameter '
               'used for the aperture extractions) estimated with '
               'uvEmissionlineSearch.TDOSE_sourcecat_from_infofile() on '+kbs.DandTstr2()+'\n')
    fout.write('# \n')
    fout.write('# id xsize ysize \n')
    for aa, arrid in enumerate(aimage_arr_inarcsec[:,0]):
        cutoutwidth = 4.0 * float(aimage_arr_inarcsec[aa,1])
        if (cutoutwidth >= 0.0) & (cutoutwidth < minCutwidth):
            cutoutwidth = minCutwidth
        if cutoutwidth > 10.0:
            cutoutwidth = 10.0000

        fout.write(str(arrid)+' '+
                   str(cutoutwidth)+' '+
                   str(cutoutwidth)+'  \n')
    fout.close()

    # ===== APERTURE SIZES =====
    aper_outtxt = outputdir+'apertureradii_MWuves_2xA_IMAGE.txt'
    fout        = open(aper_outtxt,'w')
    fout.write('# Aperturesizes of 2*A_IMAGE estimated with '
               'uvEmissionlineSearch.TDOSE_sourcecat_from_infofile() on '+kbs.DandTstr2()+'\n')
    fout.write('# \n')
    fout.write('# id aperturesize \n')
    for aa, arrid in enumerate(aimage_arr_inarcsec[:,0]):
        Raper_arcsec = 2.0 * float(aimage_arr_inarcsec[aa,1])
        if (Raper_arcsec >= 0.0) & (Raper_arcsec < minRaper):
            Raper_arcsec = minRaper
        fout.write(str(arrid)+' '+str(Raper_arcsec)+'  \n')
    fout.close()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_infofile_nondetections(infofile,goodmatchsep=0.25,outdir=None,magcuts=None,withheader=False,
                               magtocuton=['f814wcand','f160w','f814w','f775w'],verbose=True):
    """
    Generate object id lists of objects from the infofile which are (likely) not detected in the imaging
    and which should therfore be extracted as point sources with the "nondetections" keyword in the TDOSE setup files.

    --- INPUT ---
    infofile          infofile from uves.build_LAEfitstable() to base source catalogs on
    goodmatchsep      Any match to the photometric catalog below goodmatchsep is considered detectable
    outdir            To save results to output ascii files instead of printing IDs to screen provide an output directory
    magcuts           Applying magnitude cuts to the selection. Provide a list of 4 cuts to be applied
                      to [cdfs,cdfs-parallel,cosmos,udf] photometry.
    withheader        TDOSE expects no header in id list, but setting withheader=True will add on with info to each file
                      either way (the summary file always contains a header with the info)
    magtocuton        Corresponding to the magcuts, provide the names of the filters to cut on
    verbose           Toggle verbosity

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo_UVemitters_3timesUDFcats.fits'
    outdir        = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/nondetection_lists/'
    uves.get_infofile_nondetections(infofile,goodmatchsep=0.25,outdir=outdir,magcuts=[27.2,26.4,26.5,29.5],verbose=True) # mag cuts are Skelton and Rafelski 5sigma limiting depths

    """
    if outdir is not None:
        out_summary = outdir+'uves_nondetections_summary.txt'
        fout = open(out_summary,'w')
        fout.write('# List of MUSE ids to treat as non-detections for TDOSE extractions generated with uves.get_infofile_nondetections() on '+kbs.DandTstr2()+'\n')
        fout.write('# Selection based on: \n')
        fout.write('# goodmatchsep                                = '+str(goodmatchsep)+' \n')
        fout.write('# magtocuton[cdfs,cdfs-parallel,cosmos,udf]   = '+str(magtocuton)+' \n')
        fout.write('# magcuts[cdfs,cdfs-parallel,cosmos,udf]      = '+str(magcuts)+' \n# \n')

    infodat        = afits.open(infofile)[1].data

    skelcosmos     = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat.FITS'
    skelcosmosdat  = afits.open(skelcosmos)[1].data

    skelcdfs       = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
    skelcdfsdat    = afits.open(skelcdfs)[1].data

    rafelskicat    = '/Users/kschmidt/work/catalogs/rafelski/uvudf_rafelski_2015.fits'
    rafelskidat    = afits.open(rafelskicat)[1].data

    pointings      = np.unique(infodat['pointing'])

    if magcuts is None:
        if verbose: print(' - WARNING: No magnitude cuts provided')
        magcuts = [35,35,35,35]
    print(' - Including objects with matches brighter than: \n   [cdfs_'+str(magtocuton[0])+',cdfs-parallels_'+str(magtocuton[1])+',cosmos_'+str(magtocuton[2])+',udf_'+str(magtocuton[3])+'] = '+str(magcuts))

    if verbose & (outdir is None):
        print(' - - - - - - - - - - - - - - - OBJECTS FOR NONDETECTION KEYWORD - - - - - - - - - - - - - - - ')
    for pp, pointing in enumerate(pointings):
        pointingent     = np.where(infodat['pointing'] == pointing)
        ids_point       = infodat['id'][pointingent]
        match_sep_Skel  = infodat['sep_skelton'][pointingent]
        match_id_Skel   = infodat['id_skelton'][pointingent]
        match_sep_Raf   = infodat['sep_rafelski'][pointingent]
        match_id_Raf    = infodat['id_rafelski'][pointingent]

        pointing_nondetections   = []
        pointing_nondetections_s = []  # seperate list for UDF skelton objects
        for ii, museid in enumerate(ids_point):

            # magcuts = [cdfs,cdfs-parallel,cosmos,udf,udf10]
            if 'cdfs' in pointing:
                match_id  = match_id_Skel[ii]
                match_sep = match_sep_Skel[ii]
                phot_ent  = np.where(skelcdfsdat['id'] == match_id)[0]
                phot_mag  = 25.0 - 2.5**np.log10(skelcdfsdat['f_'+magtocuton[0]][phot_ent])
                mag_cut   = magcuts[0]

            elif 'hudf09' in pointing:
                match_id  = match_id_Skel[ii]
                match_sep = match_sep_Skel[ii]
                phot_ent  = np.where(skelcdfsdat['id'] == match_id)[0]
                phot_mag  = 25.0 - 2.5**np.log10(skelcdfsdat['f_'+magtocuton[1]][phot_ent])
                mag_cut   = magcuts[1]

            elif 'cosmos' in pointing:
                match_id  = match_id_Skel[ii]
                match_sep = match_sep_Skel[ii]
                phot_ent  = np.where(skelcosmosdat['id'] == match_id)[0]
                phot_mag  = 25.0 - 2.5*np.log10(skelcosmosdat['f_'+magtocuton[2]][phot_ent])
                mag_cut   = magcuts[2]

            elif ('mosaic' in pointing) or ('udf-' in pointing):
                match_id  = match_id_Raf[ii]
                match_sep = match_sep_Raf[ii]
                phot_ent  = np.where(rafelskidat['ID'] == match_id)[0]
                phot_mag  = rafelskidat['MAG_'+magtocuton[1].upper()][phot_ent]
                mag_cut   = magcuts[3]

                match_id_s  = match_id_Skel[ii]
                match_sep_s = match_sep_Skel[ii]
                phot_ent_s  = np.where(skelcdfsdat['id'] == match_id_s)[0]
                phot_mag_s  = 25.0 - 2.5**np.log10(skelcdfsdat['f_'+magtocuton[0]][phot_ent_s])
                mag_cut_s   = magcuts[0]

            else:
                sys.exit(' The pointing '+pointing+' was not recognised')


            if match_sep > goodmatchsep:
                pointing_nondetections.append(museid)
            else:
                if phot_mag > mag_cut:
                    pointing_nondetections.append(museid)

            if ('mosaic' in pointing) or ('udf-' in pointing):
                if match_sep_s > goodmatchsep:
                    pointing_nondetections_s.append(museid)
                else:
                    if phot_mag_s > mag_cut_s:
                        pointing_nondetections.append(museid)

        if verbose & (outdir is None):
            print('   '+str("%30s" % pointing)+' (Nobj='+str("%5i" % len(pointing_nondetections))+
                  ')  '+str(pointing_nondetections))

            if ('mosaic' in pointing) or ('udf-' in pointing):
                if verbose: print('   '+str("%30s" % (pointing+'_SKELTON'))+' (Nobj='+str("%5i" % len(pointing_nondetections_s))+
                                  ')  '+str(pointing_nondetections_s))
        else:
            fout.write('   '+str("%30s" % pointing)+' (Nobj='+str("%5i" % len(pointing_nondetections))+
                       ')  '+str(pointing_nondetections)+'\n')

            if ('mosaic' in pointing) or ('udf-' in pointing):
                fout.write('   '+str("%30s" % (pointing+'_SKELTON'))+' (Nobj='+str("%5i" % len(pointing_nondetections_s))+
                           ')  '+str(pointing_nondetections_s)+'\n')

        out_pointing = outdir+'uves_nondetections_'+pointing+'.txt'
        fout_p = open(out_pointing,'w')
        if withheader:
            fout_p.write('# List of MUSE ids to treat as non-detections for TDOSE extractions generated with uves.get_infofile_nondetections() on '+kbs.DandTstr2()+'\n')
            fout_p.write('# Selection based on: \n')
            fout_p.write('# goodmatchsep           = '+str(goodmatchsep)+' \n')
        if 'cdfs' in pointing:
            if withheader:
                fout_p.write('# magtocuton[cdfs]   = '+str(magtocuton[0])+' \n')
                fout_p.write('# magcuts[cdfs]      = '+str(magcuts[0])+' \n# \n')
            for nonid in pointing_nondetections: fout_p.write(str(nonid)+'\n')
        elif 'hudf09' in pointing:
            if withheader:
                fout_p.write('# magtocuton[cdfs-parallel]   = '+str(magtocuton[1])+' \n')
                fout_p.write('# magcuts[cdfs-parallel]      = '+str(magcuts[1])+' \n# \n')
            for nonid in pointing_nondetections: fout_p.write(str(nonid)+'\n')
        elif 'cosmos' in pointing:
            if withheader:
                fout_p.write('# magtocuton[cosmos]   = '+str(magtocuton[2])+' \n')
                fout_p.write('# magcuts[cosmos]      = '+str(magcuts[2])+' \n# \n')
            for nonid in pointing_nondetections: fout_p.write(str(nonid)+'\n')
        elif ('mosaic' in pointing) or ('udf-' in pointing):
            if withheader:
                fout_p.write('# magtocuton[udf]   = '+str(magtocuton[3])+' \n')
                fout_p.write('# magcuts[udf]      = '+str(magcuts[3])+' \n# \n')
            for nonid in pointing_nondetections: fout_p.write(str(nonid)+'\n')

            out_pointing_skelton = outdir+'uves_nondetections_'+pointing+'_skeltonbased.txt'
            fout_pS = open(out_pointing_skelton,'w')
            if withheader:
                fout_pS.write('# List of MUSE ids to treat as non-detections for TDOSE extractions generated with uves.get_infofile_nondetections() on '+kbs.DandTstr2()+'\n')
                fout_pS.write('# Selection based on: \n')
                fout_pS.write('# goodmatchsep           = '+str(goodmatchsep)+' \n')
                fout_pS.write('# magtocuton[udf]   = '+str(magtocuton[3])+' \n')
                fout_pS.write('# magcuts[udf]      = '+str(magcuts[3])+' \n# \n')
            for nonid in pointing_nondetections_s: fout_pS.write(str(nonid)+'\n')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_noise_spectrum(outfile='/Users/kschmidt/work/MUSE/spectra_noise/median_eff_noise_spectrum_RENAME.fits',
                         overwrite=False,verbose=True):
    """
    Build a median noise spectrum for the MUSE-Wide fields

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    uves.build_noise_spectrum(outfile='/Users/kschmidt/work/MUSE/median_eff_noise_spectrum.fits')

    """
    if os.path.isfile(outfile) & (not overwrite):
        sys.exit(' - The output file '+outfile+' exists and overwrite=False')

    #### UDF10 ####
    datacubes = ['/Users/kschmidt/work/MUSE/QtClassify/UDF10/udf-10_mfs-and-effvar-cube.fits',
                 '/Users/kschmidt/work/MUSE/QtClassify/UDF10/udf-10_mfs-and-effvar-cube.fits']
    noiseext  = 'EFFVAR'
    waveunits = 'Angstrom'
    fluxunits = '1e-20 erg/s/cm2/A'

    #### UDF ####
    # datacubes = glob.glob('/Users/kschmidt/work/MUSE/QtClassify/*/udf-*_mfs-and-effvar-cube.fits')
    # noiseext  = 'EFFVAR'
    # waveunits = 'Angstrom'
    # fluxunits = '1e-20 erg/s/cm2/A'

    #### CDFS/COSMOS MUSE-Wide ####
    # datacubes = glob.glob('/Volumes/DATABCKUP1/MUSE-Wide/DATACUBES/DATACUBE_candels-*_v1.0_dcbgc_effnoised.fits')
    # noiseext  = 'EFF_STAT'
    # waveunits = 'Angstrom'
    # fluxunits = '1e-20 erg/s/cm2/A'

    if verbose: print(' - Generating median vec for: ')
    for dd, cube in enumerate(datacubes):
        if verbose: print('   '+cube+'  (spec '+str(dd+1)+'/'+str(len(datacubes))+')')
        effstatarr = afits.open(cube)[noiseext].data
        mediannoisevec = np.sqrt(np.nanmedian(np.nanmedian(effstatarr[:,150:250,150:250],axis=1),axis=1))
        if dd == 0:
            noisearr               = mediannoisevec
            cubehdr_0              = afits.open(cube)[noiseext].header
            wavevec_0              = np.arange(cubehdr_0['NAXIS3'])*cubehdr_0['CD3_3']+cubehdr_0['CRVAL3']
        else:
            cubehdr                = afits.open(cube)[noiseext].header
            wavevec                = np.arange(cubehdr['NAXIS3'])*cubehdr['CD3_3']+cubehdr['CRVAL3']
            func                   = scipy.interpolate.interp1d(wavevec,mediannoisevec,kind='linear',fill_value="extrapolate")
            mediannoisevec_interp  = func(wavevec_0)

            noisearr  = np.vstack([noisearr,mediannoisevec_interp])

    noisevec = np.nanmedian(noisearr,axis=0)

    felis.save_spectrum(outfile,wavevec,noisevec,noisevec*0.0,
                        headerinfo=None,waveunits=waveunits,fluxunits=fluxunits,
                        overwrite=overwrite,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def strip_cube_from_spectra(specdir,outputdir,overwrite=False,verbose=True):
    """
    Remoce source cube extension from spectra.

    --- INPUT ---
    specdir            Directory containing spectra to strip soource cube extension from
    outputdir          Directory to store the stripped spectra to.
    overwrite          Overwrite output if it already exists.
    verbose            Toggle verbosity

    --- EXAMPLE OF RUN ---
    import uvEmissionlineSearch as uves

    outputdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/spectra_aperture/'
    specdir   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/spectra_test/original/'

    uves.strip_cube_from_spectra(specdir,outputdir,overwrite=False,verbose=True)

    """
    spectra   = glob.glob(specdir+'tdose_spectrum_*.fits')
    for spectrum in spectra:
        tu.strip_extension_from_fitsfile(spectrum,outputdir,removeextension='SOURCECUBE',overwrite=overwrite,verbose=verbose)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def checkfluxscales(specWnoise=True):
    """
    function to check conversion of output flux scales to integrated line fluxes


    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    uves.checkfluxscales()

    """
    Nsigma  = 3.0

    testdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FELIStesting/felistest190910/'
    if specWnoise:
        noisestr = 'noisespec'
        pklfile  = testdir+ 'uves_mock_spectrum_fromsetup_CIIIdoublet_noisespec_sigma2p00_skew0p00_' \
                            'Ftot1028p57_Fratio1p40_z2p70_CCresults_templateCIII_matchto_spectrumCIIIdoublet.pkl'
    else:
        noisestr = 'noisestdNone'
        pklfile  = testdir+ 'uves_mock_spectrum_fromsetup_CIIIdoublet_noisestdNone_sigma2p00_skew0p00_' \
                            'Ftot1028p57_Fratio1p40_z2p70_CCresults_templateCIII_matchto_spectrumCIIIdoublet.pkl'

    loaddic  = felis.load_picklefile(pklfile)
    spectrum = loaddic.keys()[0]

    template, vshift_intr, vshift_match, alpha, alphaerr, S2Nmax, Ngoodent, chi2, zspec, zS2Nmax = \
        felis.getresult4maxS2N(loaddic,spectrum,zspecISzLya=False)

    # - - - - - - - - - - - Load spec and temp info - - - - - - - - - - - - -
    specdat   = afits.open(spectrum)[1].data
    spechdr   = afits.open(spectrum)[1].header
    sigma_obs = spechdr['FLINE1_2']
    fluxline1 = spechdr['FLINE1_4']
    fluxline2 = spechdr['FLINE2_4']
    fluxratio = fluxline1/fluxline2

    threesigma_obs = sigma_obs * Nsigma
    wavemin = spechdr['FLINE1_1']-threesigma_obs
    wavemax = spechdr['FLINE2_1']+threesigma_obs

    tempdat = afits.open(template)[1].data
    temphdr = afits.open(template)[1].header
    temp_sigma     = temphdr['FLINE1_2']
    temp_fluxline1 = temphdr['FLINE1_4']
    temp_fluxline2 = temphdr['FLINE2_4']
    temp_fluxratio = temp_fluxline1/temp_fluxline2

    threesigma_rf = temp_sigma * Nsigma
    wavemin_rf = spechdr['FLINE1_1']/(1.0+zspec)-threesigma_rf
    wavemax_rf = spechdr['FLINE2_1']/(1.0+zspec)+threesigma_rf
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    specent   = np.where((specdat['wave'] > wavemin) & (specdat['wave'] < wavemax))[0]
    tempent   = np.where((tempdat['wave'] > wavemin_rf) & (tempdat['wave'] < wavemax_rf))[0]

    print(' --------- Observed frame: ----------')
    Ftot_spec_full = np.sum(specdat['flux'])*np.median(np.diff(specdat['wave']))
    Ftot_temp_full = np.sum(tempdat['flux']/(1+zspec))*np.median(np.diff(tempdat['wave']*(1+zspec)))
    print('Ftot_spec_full      = '+str(Ftot_spec_full))
    print('Ftot_temp_full      = '+str(Ftot_temp_full))

    Ftot_spec_cut = np.sum(specdat['flux'][specent])*np.median(np.diff(specdat['wave'][specent]))
    Ftot_temp_cut = np.sum(tempdat['flux'][tempent]/(1+zspec))*np.median(np.diff(tempdat['wave'][tempent]*(1+zspec)))
    print('Ftot_spec_cut       = '+str(Ftot_spec_cut))
    print('Ftot_temp_cut       = '+str(Ftot_temp_cut))

    Ftot_spec_cut_trapz = np.trapz(specdat['flux'][specent],specdat['wave'][specent])
    Ftot_temp_cut_trapz = np.trapz(tempdat['flux'][tempent]/(1+zspec),tempdat['wave'][tempent]*(1+zspec))
    print('Ftot_spec_cut_trapz = '+str(Ftot_spec_cut_trapz))
    print('Ftot_temp_cut_trapz = '+str(Ftot_temp_cut_trapz))

    #
    # print(' --------- Rest frame: ----------')
    # Ftot_spec_full = np.sum(specdat['flux']*(1+zspec))*np.median(np.diff(specdat['wave']/(1+zspec)))
    # Ftot_temp_full = np.sum(tempdat['flux'])*np.median(np.diff(tempdat['wave']))
    # Ftot_spec_cut = np.sum(specdat['flux'][specent]*(1+zspec))*np.median(np.diff(specdat['wave'][specent]/(1+zspec)))
    # Ftot_temp_cut = np.sum(tempdat['flux'][tempent])*np.median(np.diff(tempdat['wave'][tempent]))
    # print('Ftot_spec_full = '+str(Ftot_spec_full))
    # print('Ftot_temp_full = '+str(Ftot_temp_full))
    # print('Ftot_spec_cut = '+str(Ftot_spec_cut))
    # print('Ftot_temp_cut = '+str(Ftot_temp_cut))

    print(' --------- Values from headers: ----------')
    template_scaled               = tempdat['flux'] / (temp_fluxline1+temp_fluxline2) * alpha
    print('alpha +/- alphaerr     = '+str(alpha)+' +/- '+str(alphaerr))
    print('spec total flux        = '+str(fluxline1+fluxline2))
    print('temp total flux        = '+str(temp_fluxline1+temp_fluxline2))
    print('spec flux ratio        = '+str(fluxratio))
    print('temp flux ratio        = '+str(temp_fluxratio))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' --------- Plotting overview: ----------')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6, 4))
    fig.subplots_adjust(wspace=0.1, hspace=0.3,left=0.15, right=0.99, bottom=0.15, top=0.99)
    Fsize    = 14
    lthick   = 2
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    plt.plot(specdat['wave'],specdat['flux'],color='gray')
    plt.plot(tempdat['wave']*(1+zspec),template_scaled/(1+zspec),color='pink')

    plt.plot([wavemin,wavemin],[np.min(specdat['flux'][specent]),np.max(specdat['flux'][specent])],color='k',ls=':')
    plt.plot([wavemax,wavemax],[np.min(specdat['flux'][specent]),np.max(specdat['flux'][specent])],color='k',ls=':')
    plt.plot(specdat['wave'][specent],specdat['flux'][specent],color='k')

    plt.plot([wavemin_rf*(1+zspec),wavemin_rf*(1+zspec)],
             [np.min(template_scaled[tempent]/(1+zspec)),np.max(template_scaled[tempent]/(1+zspec))],color='red',ls=':')
    plt.plot([wavemax_rf*(1+zspec),wavemax_rf*(1+zspec)],
             [np.min(template_scaled[tempent]/(1+zspec)),np.max(template_scaled[tempent]/(1+zspec))],color='red',ls=':')
    plt.plot(tempdat['wave'][tempent]*(1+zspec),tempdat['flux'][tempent]/(1+zspec),color='green')

    plt.plot(tempdat['wave'][tempent]*(1+zspec),template_scaled[tempent]/(1+zspec),color='red')
    plt.ylim([-50,140])

    plt.xlabel('$\\lambda$ [\AA]')
    plt.ylabel('flux [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]')

    plotname = testdir+'fluxscalecheck_'+noisestr+'.pdf'
    plt.savefig(plotname)
    # plt.savefig('/Users/kschmidt/Desktop/fluxscalecheck_'+noisestr+'.pdf')
    print('Saved plot to '+plotname)
    plt.clf()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_UDF10masedaobjcomparison(gaussspec=False,overwrite=False,verbose=True):
    """
    plotting information on FELIS match to Maseda UDF-10 CIII emitters

    --- INPUT ---
    overwrite          Overwrite the plots if they already exist?
    verbose            Toggle verbosity

    --- EXAMPLE OF RUN ---
    import uvEmissionlineSearch as uves
    uves.plot_UDF10masedaobjcomparison(gaussspec=True,overwrite=True)
    uves.plot_UDF10masedaobjcomparison(gaussspec=False,overwrite=True)

    """
    Nsigmaplot   = 3.0 # the size of the error bars to show
    outdir       = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FELIStemplatematch2uvesobjects/CCresults_summary/'

    summaryfile  = outdir+'CCresults_summary_templateCIII_FELISmatch2udf10masedaobj190913.txt'
    if gaussspec:
        summaryfile = summaryfile.replace('.txt','_gauss.txt')
    fmt          = 'd,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,300a,300a'
    summarydat   = np.genfromtxt(summaryfile,skip_header=24,dtype=fmt,comments='#',names=True)
    sortindex_S  = np.argsort(summarydat['id'])
    Nspecin      = len(summarydat['spectrum'])

    masedainfo   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/masedaUDF10emitters.txt'
    fmt          = '12a,d,d,d,d,d,d,d,d,d,d,d'
    masedadat    = np.genfromtxt(masedainfo,skip_header=1,dtype=fmt,comments='#',names=True)
    sortindex_M  = np.argsort(masedadat['id_uves'])

    plotbasename = outdir+'UDF10_CIIIemitters_Maseda17comparison'
    if verbose: print(' - Plotting FELIS matches in summary file\n   '+summaryfile+'\n   where the following holds:')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext    = 'Comparison_onetoone_Fciii'
    plotname   = plotbasename+nameext+'.pdf'
    if gaussspec:
        plotname = plotname.replace('.pdf','_gauss.pdf')
    xvalues    = masedadat['f_ciii'][sortindex_M]
    xerr       = masedadat['df_ciii'][sortindex_M]*Nsigmaplot
    yvalues    = summarydat['Ftot_FELIS_S2Nmax'][sortindex_S]
    yerr       = summarydat['Ftot_FELIS_S2Nmax_err'][sortindex_S]*Nsigmaplot
    xlabel     = 'Flux(CIII,Maseda)'
    ylabel     = 'Flux(CIII,FELIS)'
    point_text = (summarydat['id'][sortindex_S].astype(int)).astype(str)


    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,'dummydat',
                                                   histaxes=True,Nbins=30,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'][sortindex_S],
                                                   linetype='onetoone',point_text=point_text,
                                                   xlog=True,ylog=True,xrange=[10,1200],yrange=[10,1200],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext    = 'Comparison_horizontal_FciiiVSz'
    plotname   = plotbasename+nameext+'.pdf'
    if gaussspec:
        plotname = plotname.replace('.pdf','_gauss.pdf')
    xvalues    = summarydat['z_temp_S2Nmax'][sortindex_S]
    xerr       = [None]*len(xvalues)
    yvalues    = (summarydat['Ftot_FELIS_S2Nmax'][sortindex_S]/masedadat['f_ciii'][sortindex_M]) - 1
    yerr       = np.sqrt((summarydat['Ftot_FELIS_S2Nmax_err'][sortindex_S] /
                          summarydat['Ftot_FELIS_S2Nmax'][sortindex_S])**2.0+
                         (masedadat['df_ciii'][sortindex_M] /
                          masedadat['f_ciii'][sortindex_M])**2.0) * Nsigmaplot
    xlabel     = '$z$(FELIS)'
    ylabel     = '[Flux(CIII,FELIS)-Flux(CIII,Maseda)] - 1'
    point_text = (summarydat['id'][sortindex_S].astype(int)).astype(str)


    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,'dummydat',
                                                   histaxes=True,Nbins=30,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'][sortindex_S],
                                                   linetype='horizontal',point_text=point_text,
                                                   xlog=False,ylog=False,xrange=[1.0,3.0],yrange=[-0.5,0.5],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nameext    = 'Comparison_horizontal_FciiiVSsigma'
    plotname   = plotbasename+nameext+'.pdf'
    if gaussspec:
        plotname = plotname.replace('.pdf','_gauss.pdf')
    xvalues    = summarydat['sigma_temp_ang_rf'][sortindex_S]
    xerr       = [None]*len(xvalues)
    yvalues    = (summarydat['Ftot_FELIS_S2Nmax'][sortindex_S]/masedadat['f_ciii'][sortindex_M]) - 1
    yerr       = np.sqrt((summarydat['Ftot_FELIS_S2Nmax_err'][sortindex_S] /
                          summarydat['Ftot_FELIS_S2Nmax'][sortindex_S])**2.0+
                         (masedadat['df_ciii'][sortindex_M]/
                          masedadat['f_ciii'][sortindex_M])**2.0) * Nsigmaplot
    xlabel     = '$\sigma$(FELIS)'
    ylabel     = '[Flux(CIII,FELIS)-Flux(CIII,Maseda)] - 1'
    point_text = (summarydat['id'][sortindex_S].astype(int)).astype(str)


    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,'dummydat',
                                                   histaxes=True,Nbins=30,
                                                   colortype='s2nfelis',cdatvec=summarydat['FELIS_S2Nmax'][sortindex_S],
                                                   linetype='horizontal',point_text=point_text,
                                                   xlog=False,ylog=False,xrange=[0.2,0.8],yrange=[-0.5,0.5],
                                                   colorcode=True,overwrite=overwrite,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_FELISmatches(objectids,pickledir,summaryfiles,outputdir,S2Nmin=3.0,vshiftmax=1e4,
                      plot_overview=False,skipplotting=False,
                      showFELISresults=False,FELISvetting=None,verbose=True,verboseplots=False):
    """
    Wrapper to search for objects with line detections and then plot
    the corresponding FELIS template overview.

    --- INPUT ---
    objectids          The ids of the objects to check and plot. Either provide a list or a "glob string" to append pickledir.
    pickledir          The picklefiles to plot information from. Generated with uves.match_tdosespectra_to_templates()
    summaryfiles       The summaryfile generate with uves.gen_tdosespecFELISresults_summary()
    outputdir          Directory to store plots in
    S2Nmin             Minimum FELIS S2N of line match to plot results for
    vshiftmax          Maximum velocity shift of line match to plot results for
    plot_overview      Plot the overview of the spectrum with zoom-in on lines using mwp.
    verbose            Toggle verbosity
    verboseplots       Toggle verbosity of plotting functions

    --- EXAMPLE OF RUN ---
    import uvEmissionlineSearch as uves

    parentdir   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FELIStemplatematch2uvesobjects/'
    objectids   = '7201*'
    summaryfile = glob.glob(parentdir+'CCresults_summary/CCresults_summary_template*_FELISmatch2udf10_190913.txt')
    outputdir   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FELIStemplatematch2uvesobjects/FELISmatches_plots/'
    uves.plot_FELISmatches(objectids,parentdir,summaryfile,outputdir,S2Nmin=3.0,plot_overview=True)

    """

    if type(objectids) is str:
        pfiles    = glob.glob(pickledir+'*'+objectids+'*.pkl')
        idlist    = [int(pfile.split('_CCresults_')[0].split('_')[-1]) for pfile in pfiles]
        objidlist = np.unique(np.asarray(idlist))
    else:
        objidlist = np.unique(np.asarray(objectids))

    if verbose: print('# - Loading the data of the '+str(len(summaryfiles))+' summary files into memory')
    summarfiledat_dic = {}
    for summaryfile in summaryfiles:
        fmt                             = '12a,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,300a,300a'
        summarfiledat_dic[summaryfile]  = np.genfromtxt(summaryfile,skip_header=25,dtype=fmt,comments='#',names=True)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('# - Looping over the '+str(len(objidlist))+
                      ' (unique) object IDs to look for FELIS matches with S/N(FELIS) > '+str(S2Nmin)+
                      ' and vshift < '+str(vshiftmax)+' km/s ')
    matchcount = 0
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for oo, objid in enumerate(objidlist):
        if verbose & (not skipplotting): print('---- Checking object  '+str(objid)+' ----')
        NlinesMatchCuts  = 0 # keeping track of number of lines above threshold
        zobj             = np.array([])
        matchline        = np.array([])
        pfilenumber      = np.array([])
        pointings        = np.array([])
        spectra          = np.array([])
        templates        = np.array([])
        matchS2N         = np.array([])
        vshift           = np.array([])
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for summaryfile in summaryfiles:
            summarydat       = summarfiledat_dic[summaryfile]
            objsummary_ent   = np.where(summarydat['id'].astype(int) == objid)[0]
            zobj             = np.append(zobj,summarydat['z_spec'][objsummary_ent])
            linesummarized   = summaryfile.split('_template')[-1].split('_')[0]
            outkeystr        = summaryfile.split('_template'+linesummarized)[-1].split('.txt')[0]

            if len(objsummary_ent) != 0:
                if len(objsummary_ent) > 1:
                    if verbose & (not skipplotting):
                        print('     NOTE: Found '+str(len(objsummary_ent))+
                              ' spectra (template matches) in '+linesummarized+' summaryfile ')

                picklefiles = glob.glob(pickledir+'*'+str(objid)+'*template'+linesummarized+'*'+outkeystr+'*.pkl')
                for pp, pfile in enumerate(picklefiles):
                    spec             = summarydat['spectrum'][objsummary_ent[pp]]
                    pfilenumber      = np.append(pfilenumber,pp+1)
                    pointings        = np.append(pointings,spec.split('/tdose_spectrum_')[-1].split('-full')[0])
                    spectra          = np.append(spectra,spec)
                    templates        = np.append(templates,summarydat['template'][objsummary_ent[pp]])
                    matchline        = np.append(matchline,linesummarized)
                    matchS2N         = np.append(matchS2N,summarydat['FELIS_S2Nmax'][objsummary_ent[pp]])
                    vshift           = np.append(vshift,summarydat['vshift_CCmatch'][objsummary_ent[pp]])

                    if (summarydat['FELIS_S2Nmax'][objsummary_ent[pp]] > S2Nmin) & \
                            (np.abs(summarydat['vshift_CCmatch'][objsummary_ent[pp]]) < vshiftmax):
                        NlinesMatchCuts  = NlinesMatchCuts  + 1
                        for pp, pickle in enumerate(picklefiles):
                            plotname  = pickle.replace('.pkl','_CCwith_maxS2N.pdf')
                            if not skipplotting:
                                felis.plot_picklefilecontent([summarydat['spectrum'][objsummary_ent[pp]]],
                                                             pickle,showspecerr=False,plotnames=[plotname],
                                                             plotdir=outputdir,verbose=verboseplots)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        zobj = np.unique(zobj)[0]
        if verbose & (not skipplotting):
            print('     Found '+str(NlinesMatchCuts )+' lines with S2N(FELIS) > '+str(S2Nmin)+
                  ' and vshift < '+str(vshiftmax)+' km/s at z='+str("%.6f" % zobj))
        if NlinesMatchCuts  > 0:
            matchcount = matchcount + 1
            outstr = ''
            for ii, mline in enumerate(matchline):
                if verbose & (not skipplotting):
                    print('     '+mline+' @ S/N = '+str(matchS2N[ii])+'  (vshift='+str(vshift[ii])+'km/s)')
                if skipplotting:
                    outstr = outstr+' S2N('+mline+'_'+str(pointings[ii])+\
                             ')='+str(matchS2N[ii])+' w. Dv='+str(vshift[ii])+'km/s | '
            if skipplotting:
                print(str(objid)+' '+str("%12.6f" % zobj)+'   # '+outstr)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if plot_overview & (NlinesMatchCuts  > 0) & (not skipplotting):
            uspec      = np.unique(spectra)
            gaussspec  = []
            aperspec   = []
            for spec in uspec:
                gfile = spec.replace('aperture','gauss')
                if os.path.isfile(gfile):
                    gaussspec.append(gfile)

                afile = spec.replace('gauss','aperture')
                if os.path.isfile(afile):
                    aperspec.append(afile)
            if gfile != afile:
                specoverview = gaussspec + aperspec
                labels       = ['TDOSE gauss']*len(gaussspec) +['TDOSE aper']*len(aperspec)
            else:
                specoverview = gaussspec
                labels       = ['TDOSE spectrum']*len(gaussspec)

            wavecols     = ['wave']*len(specoverview)
            fluxcols     = ['flux']*len(specoverview)
            fluxerrcols  = ['fluxerror']*len(specoverview)
            redshift     = zobj
            voffset      = 0.0

            if str(objid)[0] == '1':
                skyspec     = '/Users/kschmidt/work/MUSE/spectra_noise/median_eff_noise_spectrum_70fields190819.fits'
                wavecol_sky = 'wave'
                fluxcol_sky = 'flux'
                # skyspec     = '/Users/kschmidt/work/MUSE/spectra_sky/SKY_SPECTRUM_candels-cdfs-'+str(objid)[1:3]+'_av.fits'
                # wavecol_sky = 'lambda'
                # fluxcol_sky = 'data'
            elif str(objid)[0] == '2':
                skyspec     = '/Users/kschmidt/work/MUSE/spectra_noise/median_eff_noise_spectrum_70fields190819.fits'
                wavecol_sky = 'wave'
                fluxcol_sky = 'flux'
                # skyspec     = '/Users/kschmidt/work/MUSE/spectra_sky/SKY_SPECTRUM_candels-cosmos-'+str(objid)[1:3]+'_av.fits'
                # wavecol_sky = 'lambda'
                # fluxcol_sky = 'data'
            elif str(objid)[0] == '6':
                skyspec     = '/Users/kschmidt/work/MUSE/spectra_noise/median_eff_noise_spectrum_UDF190819.fits'
                wavecol_sky = 'wave'
                fluxcol_sky = 'flux'
            elif (str(objid)[0] == '3') or (str(objid)[0] == '4'):
                skyspec     = '/Users/kschmidt/work/MUSE/spectra_noise/median_eff_noise_spectrum_70fields190819.fits'
                wavecol_sky = 'wave'
                fluxcol_sky = 'flux'
                # skyspec     = '/Users/kschmidt/work/MUSE/spectra_sky/SKY_SPECTRUM_candels-cdfs-20_av.fits'
                # wavecol_sky = 'lambda'
                # fluxcol_sky = 'data'
            elif str(objid)[0] == '7':
                skyspec     = '/Users/kschmidt/work/MUSE/spectra_noise/median_eff_noise_spectrum_UDF10_200925.fits'
                wavecol_sky = 'wave'
                fluxcol_sky = 'flux'
            else:
                skyspec     = 'None'
                wavecol_sky = 'None'
                fluxcol_sky = 'None'

            skyspectra   = [None]*len(specoverview)
            wavecols_sky = [wavecol_sky]*len(specoverview)
            fluxcols_sky = [fluxcol_sky]*len(specoverview)

            if os.path.isfile(skyspec):
                skyspectra[0] = skyspec

            yrangefull   = None
            xrangefull   = [4600,9400]
            #linenames    = ['Lyb+OVI','Lya+NV','CIV','HeII','OIII1663','SiIII','CIII','MgII','OII']
            linenames    = ['Lyb+OVI','Lya','NV','CIV','HeII','OIII1663','SiIII','CIII','MgII']
            outputfigure = outputdir+'overview_1DspecWzooms_'+str(objid)+'.pdf'

            mwp.plot_1DspecOverview(specoverview, labels, wavecols, fluxcols, fluxerrcols, redshift,
                                    voffset=voffset, skyspectra=skyspectra, wavecols_sky=wavecols_sky,
                                    fluxcols_sky=fluxcols_sky, outputfigure=outputfigure,linenames=linenames,
                                    yrangefull=yrangefull, xrangefull=xrangefull, plotSN=False,
                                    showFELISresults=showFELISresults,FELISvetting=FELISvetting,verbose=verboseplots)

            mwp.plot_1DspecOverview(specoverview, labels, wavecols, fluxcols, fluxerrcols, redshift,
                                    voffset=voffset, skyspectra=skyspectra, wavecols_sky=wavecols_sky,
                                    fluxcols_sky=fluxcols_sky, outputfigure=outputfigure,linenames=linenames,
                                    yrangefull=yrangefull, xrangefull=xrangefull, plotSN=True,verbose=verboseplots)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose:
        print('# - Done; so that was '+str(matchcount)+' objects with at least one line matching the cuts: ')
        print('#   S2N(FELIS) > '+str(S2Nmin))
        print('#   vshift     < '+str(vshiftmax)+' km/s ')


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_FELISs2nMap(felisresults,outputdir='./',S2Nmarkerthreshold=3.0,verbose=True):
    """
    Function plotting a map of the S/N obtained from fitting a given set of templates

    --- INPUT ---
    felisresults    = The FELIS output pickle file to grab information for results from

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    plotdir      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FELIStemplatematch2uvesobjects/all_gauss190926/s2nmaps/'
    felisresults = plotdir+'../tdose_spectrum_udf-10-full-v1p0-MWuves_gauss_0720680577_CCresults_templateCIII_FELISmatch2uves190926_gauss.pkl'
    uves.plot_FELISs2nMap(felisresults,outputdir=plotdir)


    """
    resultsdic = felis.load_picklefile(felisresults)

    if verbose: print(' --------- Loading content of results dictionary key: ---------')
    for dickey in resultsdic.keys():
        if verbose: print(' - Loading data from \n   '+dickey)
        zvals      = resultsdic[dickey]['zCCmaxvec']
        S2Nvals    = resultsdic[dickey]['S2NCCmaxvec']
        templates  = resultsdic[dickey]['templatevec']
        sigmavals  = np.asarray([str(tt).split('sigma')[-1].split('_skew')[0].replace('p','.') for tt in templates]).astype(float)
        fratiovals = np.asarray([str(tt).split('Fratio')[-1].split('_z')[0].replace('p','.') for tt in templates]).astype(float)

        tempbase = templates[0].split('/')[-1].split('_noise')[0]
        plotname = outputdir+dickey.split('/')[-1].replace('.fits','_fitto_'+tempbase+'_templates_FELISresultsS2Nmap.pdf')

        if verbose: print('   Setting up and generating plot:\n   '+plotname)
        fig = plt.figure(figsize=(5, 4))
        fig.subplots_adjust(wspace=0.1, hspace=0.5,left=0.13, right=0.99, bottom=0.13, top=0.98)
        Fsize    = 11
        lthick   = 2
        marksize = 10
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()

        cmin        = np.min(S2Nvals)
        cmax        = np.max(S2Nvals)

        cmap        = plt.cm.viridis_r
        colnorm     = plt.Normalize(vmin=cmin,vmax=cmax)
        cmaparr     = np.linspace(cmin, cmax, num=50)
        tempcolors  = [cmap(colnorm(S2Nval)) for S2Nval in S2Nvals]
        m           = plt.cm.ScalarMappable(cmap=cmap)
        m.set_array(cmaparr)
        cb          = plt.colorbar(m)
        cb.set_label('Max S/N of template cross correlation')

        for cc, tempcol in enumerate(tempcolors):
            if S2Nvals[cc] == np.max(S2Nvals):
                tempmarker = 'X'
            elif S2Nvals[cc] < S2Nmarkerthreshold:
                tempmarker = 'o'
            else:
                tempmarker = 's'

            plt.plot(sigmavals[cc],fratiovals[cc],tempmarker,color=tempcol,markersize=marksize,zorder=10)

        plt.xlim([np.min(sigmavals)-0.1,np.max(sigmavals)+0.1])
        plt.ylim([np.min(fratiovals)-0.1,np.max(fratiovals)+0.1])
        plt.xlabel(' Template $\sigma$ [\AA] ')
        plt.ylabel(' Template doublet flux ratio [F(blue)/F(red)]')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print('   succesfully saved plot')
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def count_SpecOnArche(countfields):
    """
    Small script that can be copy-pasted into ipython on Arche to count the spectra and keep track of the reductions.
    Before this script I did copy-paste from muse notes Latex file.

    """
    import glob, sys
    setupglobstr       = 'tdose_setupfiles/tdose_setupfile_*_gauss.txt'
    setupglobstr_reext = 'tdose_setupfiles/tdose_setupfile_*_gauss_reext.txt'
    setupfiledir = {}
    setupfiledir['rafelski']   = list(np.sort(glob.glob('Rafelski-UDF-*/'+setupglobstr)))
    setupfiledir['rafudfmock'] = list(np.sort(glob.glob('Rafelski-UDF-MWmock-v1p0/'+setupglobstr)))
    setupfiledir['rafudffull'] = list(np.sort(glob.glob('Rafelski-UDF-full-v1p0/'+setupglobstr)))

    setupfiledir['skelton']         = list(np.sort(glob.glob('Skelton-*/'+setupglobstr)))
    setupfiledir['skelcdfsfull']    = list(np.sort(glob.glob('Skelton-CDFS-full-v1p0/'+setupglobstr)))
    setupfiledir['skelcdfsparfull'] = list(np.sort(glob.glob('Skelton-CDFSparallel-full-v1p0/'+setupglobstr)))
    setupfiledir['skeludffull']     = list(np.sort(glob.glob('Skelton-UDF-full-v1p0/'+setupglobstr)))
    setupfiledir['skeludf10full']   = list(np.sort(glob.glob('Skelton-UDF10-full-v1p0/'+setupglobstr)))
    setupfiledir['skeludfmock']     = list(np.sort(glob.glob('Skelton-UDF-MWmock-v1p0/'+setupglobstr)))
    setupfiledir['skelcosfull']     = list(np.sort(glob.glob('Skelton-COSMOS-full-v1p0/'+setupglobstr)))

    setupfiledir['uves']        = list(np.sort(glob.glob('MWuves100full/MWuves-*/'+setupglobstr)))
    setupfiledir['uvescos']     = list(np.sort(glob.glob('MWuves100full/MWuves-COSMOS*/'+setupglobstr)))
    setupfiledir['uvescdfs']    = list(np.sort(glob.glob('MWuves100full/MWuves-CDFS-*/'+setupglobstr)))
    setupfiledir['uvescdfspar'] = list(np.sort(glob.glob('MWuves100full/MWuves-CDFSpar*/'+setupglobstr)))
    setupfiledir['uvesudf']     = list(np.sort(glob.glob('MWuves100full/MWuves-UDF-*/'+setupglobstr)))
    setupfiledir['uvesudf10']   = list(np.sort(glob.glob('MWuves100full/MWuves-UDF10*/'+setupglobstr)))
    setupfiledir['uvesreext']   = list(np.sort(glob.glob('MWuves100full/*/'+setupglobstr_reext)))

    setupfiledir['laigle']      = list(np.sort(glob.glob('Laigle-*/'+setupglobstr)))

    setupfiledir['guo']            = list(np.sort(glob.glob('Guo-*/'+setupglobstr)))
    setupfiledir['guocdfsfull']    = list(np.sort(glob.glob('Guo-CDFS-full-v1p0/'+setupglobstr)))
    setupfiledir['guoudffull']     = list(np.sort(glob.glob('Guo-UDF-full-v1p0/'+setupglobstr)))
    setupfiledir['guoudf10full']   = list(np.sort(glob.glob('Guo-UDF10-full-v1p0/'+setupglobstr)))
    setupfiledir['guoudfmock']     = list(np.sort(glob.glob('Guo-UDF-MWmock-v1p0/'+setupglobstr)))

    setupfiledir['whitaker']    = list(np.sort(glob.glob('Whitaker-*/'+setupglobstr)))
    setupfiledir['whitcdfs']    = list(np.sort(glob.glob('Whitaker-CDFS-full-v1p0/'+setupglobstr)))
    setupfiledir['whitcdfspar'] = list(np.sort(glob.glob('Whitaker-CDFSparallel-full-v1p0/'+setupglobstr)))
    setupfiledir['whitudf']     = list(np.sort(glob.glob('Whitaker-UDF-full-v1p0/'+setupglobstr)))
    setupfiledir['whitudf10']   = list(np.sort(glob.glob('Whitaker-UDF10-full-v1p0/'+setupglobstr)))

    if countfields == 'all':
        setupfiles = setupfiledir['rafelski']+setupfiledir['skelton']+setupfiledir['uves']
    else:
        try:
            setupfiles = setupfiledir[countfields]
        except:
            sys.exit('"'+countfields+'" is not a valid setup file selection. Choices are: '+str(setupfiledir.keys()))

    print('                                 ExtractionDir               '
          ' Fieldname           Ngaussspec        Naperturespec     missingGaussspec')
    completefields   = []
    incompletefields = []
    emptyfields      = []
    for sf in setupfiles :

        extdir     = sf.split('/tdose')[0]
        fieldname  = sf.split('_gauss')[0].split('v1p0-')[-1].split('v1p0_')[-1].split('andels-')[-1]

        specsearchstr = sf.split('/tdose_setupfiles/')[0]+'/tdose_spectra/tdose_spec*'+fieldname+'*gauss*.fits'
        if '_reext' in sf:
            specsearchstr = sf.split('/tdose_setupfiles/')[0]+'/tdose_spectra_reext/tdose_spec*'+fieldname+'*gauss*.fits'

        Ng = len(glob.glob('/store/data/musewide/TDOSE/'+specsearchstr))
        Na = len(glob.glob('/store/data/musewide/TDOSE/'+specsearchstr.replace('gauss*.fits','aperture*.fits')))
        if '_reext' in sf:
            import sys
            sys.path.append('/store/data/musewide/TDOSE/TDOSE')
            import tdose_utilities as tu
            try:
                sfdat = tu.load_setup(sf,verbose=False)
                Na    = len(sfdat['sources_to_extract'])
            except:
                Na    = np.nan

        Ngmissing = Na-Ng
        outstr    = str("%50s" % extdir)+' '+str("%20s" % fieldname)+' '+\
                    str("%20s" % Ng)+' '+str("%20s" % Na)+' '+str("%20s" % Ngmissing)

        if Ngmissing > 0:
            print('\033[91m'+outstr+'\033[0m') # coloring output red if spectra are missing
            incompletefields.append(fieldname)
        else:
            print(outstr)
            if Ng > 0:
                completefields.append(fieldname)
        if Na == 0:
            emptyfields.append(fieldname)

    print('\n - Hence, the '+str(len(completefields))+'/'+str(len(setupfiles))+' completed fields are: \n'+str(completefields))

    print('\n The fields with no aperture spectra at all are: '+str(emptyfields))
    return completefields, incompletefields, emptyfields
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def vet_felisdetection(idlist,plotdirs,outputfile,lineratiosummary,S2Nmincheck=3.0,
                       emlines=['NV', 'CIV', 'HeII', 'OIII', 'SiIII', 'CIII', 'MgII'],
                       infofile='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo_UVemitters_3timesUDFcats.fits',
                       FELISsummary=False,EWestimates=False,overwrite=False,performinspections=True,verbose=True):
    """
    Script to automize the vetting of the supposed FELIS emission line detections.

    --- INPUT ---

    --- EXAMPLE OF RUN ---
    import uvEmissionlineSearch as uves

    parentdir    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FELIStemplatematch2uvesobjects/all_aperture190926/'

    idlist = [101005016,720190222,101011026,77]
    lineratiosummary = parentdir+'fluxratios/fluxratios_FELISmatch2uves190926_aperture.txt'
    EWestimates      = parentdir+'fluxratios/fluxratios_FELISmatch2uves190926_aperture_EW0estimates_191015run.txt'
    plotdirs         = [parentdir+'../FELIS_bestmatches_plots_aperture/']

    outputfile = parentdir+'vet_felisdetection_outputRENAME.txt'
    uves.vet_felisdetection(idlist,plotdirs,outputfile,lineratiosummary,S2Nmincheck=3.0,EWestimates=EWestimates,overwrite=True)

    """
    pversion = sys.version_info[0]
    Nobj   = len(idlist)
    if verbose: print(' - '+str(Nobj)+' IDs provided for vetting')
    if verbose: print(' - Preparing output file '+outputfile)
    if os.path.isfile(outputfile) & (not overwrite):
        sys.exit('The output ('+outputfile+') alreaduy exists and overwrite=False so exiting.')
    else:
        fout = open(outputfile,'w')
        fout.write('# Result from vetting the '+str(Nobj)+' IDs provided to uves.vet_felisdetection() on '+kbs.DandTstr2()+'\n')
        fout.write('# \n')
        fout.write('# Each of the "trustLINE" columns provide \n'
                   '#        1  = yes   I would trust the line detection. \n'
                   '#        0  = no    I would not the line detection. \n'
                   '#        9  = maybe I would trust the line detection. \n'
                   '#      -99  = no spectral coverage of the given line. \n'
                   '#       99  = estimated S/N < '+str(S2Nmincheck)+'  (limit from "S2Nmincheck")  \n'
                   '#      NaN  = ID missing in line flux ratio summary.\n')
        fout.write('# The columns are followed by notes on the object \n')
        fout.write('# \n')
        fout.write('# Columns are:\n')
        fout.write('#  id               pointing'+''.join([str("%13s" % ('trust'+el)) for el in emlines])+'\n')

    if verbose: print(' - Loading main info file: '+infofile)
    dat_maininfo = afits.open(infofile)[1].data

    if verbose: print(' - Loading ancillary information to display for each object ')

    frdatBadFMT     = np.genfromtxt(lineratiosummary,skip_header=7,dtype='d',comments='#',names=True)
    fmt             = 'd,12a,'+','.join((len(frdatBadFMT.dtype.names)-2)*['d'])
    dat_lineratio   = np.genfromtxt(lineratiosummary,skip_header=7,dtype=fmt,comments='#',names=True)

    if EWestimates:
        ew0BadFMT = np.genfromtxt(EWestimates,skip_header=10,dtype='d',comments='#',names=True)
        fmt       = 'd,12a,'+','.join((len(ew0BadFMT.dtype.names)-2)*['d'])
        dat_ew0   = np.genfromtxt(EWestimates,skip_header=10,dtype=fmt,comments='#',names=True)
    else:
        dat_ew0 = None

    answerkeys = {'y':1, 'n':0, 'm':9, 'nocov':-99, 'lows2n':99, 'idmissing':np.nan}
    if verbose: print(' - Loop over objects while opening figures and printing info ')
    for ii, objid in enumerate(idlist):
        objent_info      = uves.return_objent(objid,dat_maininfo,idcol='id',verbose=False)
        objent_lineratio = uves.return_objent(objid,dat_lineratio,idcol='id',verbose=False)
        if objent_lineratio is None:
            if verbose: print('--------- OBJECT  '+str("%.10d" % objid)+
                              ' ('+str("%.5d" % (ii+1))+'/'+str("%.5d" % len(idlist))+') --------- ')
            if verbose: print(' WARNING - No objects found in line ratio summary')
            if verbose: print('           Continuing to next object ')
            continue
        objpointings     = dat_lineratio['pointing'][objent_lineratio]

        mainplots = ' '
        for plotdir in plotdirs:
            plolist   = ' '.join(glob.glob(plotdir+'overview_1DspecWzooms*'+str("%.9d" % objid)+'*.pdf'))
            mainplots = mainplots + plolist
        opencommand = 'open -n -F '
        if performinspections:
            if mainplots != '':
                pipe_mainplots = subprocess.Popen(opencommand+mainplots,shell=True,executable=os.environ["SHELL"])
                time.sleep(1.1) # sleep to make sure process appears in PIDlist
                pid_mainplots  = MiGs.getPID('Preview.app',verbose=False) # get PID of png process

        for pp, objpoint in enumerate(objpointings):
            objent_lr        = objent_lineratio[pp]
            objent_ew0       = uves.return_objent([objid,objpoint],dat_ew0,idcol=['id','pointing'],verbose=False)
            if verbose: print('\n--------- OBJECT  '+str("%.10d" % objid)+' ('+str("%.5d" % (ii+1))+'/'+str("%.5d" % len(idlist))+') POINTING '+str(objpoint)+' --------- ')
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if verbose:
                print(' z(info)         = '+str("%7.4f" % dat_maininfo['redshift'][objent_info]))
                print(' FWHM(Lya)       = '+str("%8.2f" % dat_maininfo['fwhm_kms'][objent_info])+' +\-'+
                      str("%8.2f" % dat_maininfo['fwhm_kms_std'][objent_info])+'  km/s ')
                print(' Dv_redpeak(V18) = '+str("%8.2f" % dat_maininfo['red_peak_shift_V18_kms'][objent_info])+' +\-'+
                      str("%8.2f" % dat_maininfo['red_peak_shift_V18_kms_err'][objent_info])+'  km/s ')
                print(' EW0(Lya)        = '+str("%8.2f" % dat_maininfo['EW_0'][objent_info])+' +\-'+
                      str("%8.2f" % dat_maininfo['EW_0_err'][objent_info])+'  A ')
                print(' beta            = '+str("%8.2f" % dat_maininfo['beta'][objent_info])+' +\-'+
                      str("%8.2f" % dat_maininfo['beta_err'][objent_info]))
                print(' f606wJK         = '+str("%10.4f" % (1e20*dat_maininfo['flux_acs_606w'][objent_info]))+' +\-'+
                      str("%10.4f" % (1e20*dat_maininfo['flux_err_acs_606w'][objent_info]))+'  1e-20cgs')
                print(' f814wJK         = '+str("%10.4f" % (1e20*dat_maininfo['flux_acs_814w'][objent_info]))+' +\-'+
                      str("%10.4f" % (1e20*dat_maininfo['flux_err_acs_814w'][objent_info]))+'  1e-20cgs')
                print(' id_Laigle       = '+str("%8.2f" % dat_maininfo['id_Laigle'][objent_info])+
                      '         sep = '+str("%8.2f" % dat_maininfo['sep_Laigle'][objent_info]))
                print(' id_Skelton      = '+str("%8.2f" % dat_maininfo['id_skelton'][objent_info])+
                      '         sep = '+str("%8.2f" % dat_maininfo['sep_skelton'][objent_info]))
                print(' id_Guo          = '+str("%8.2f" % dat_maininfo['id_guo'][objent_info])+
                      '         sep = '+str("%8.2f" % dat_maininfo['sep_guo'][objent_info]))
                print(' id_Rafelski     = '+str("%8.2f" % dat_maininfo['id_rafelski'][objent_info])+
                      '         sep = '+str("%8.2f" % dat_maininfo['sep_rafelski'][objent_info])+'\n')
                for el in emlines:
                    print(' EW0('+str("%5s" % el)+') = '+str("%8.2f" % dat_ew0['EW0_'+el][objent_ew0])+' +\-'+
                          str("%8.2f" % dat_ew0['EW0err_'+el][objent_ew0])+'   for beta = '+
                          str("%8.4f" % dat_ew0['beta'][objent_ew0]))
                print('\n')
                for el in emlines:
                    printstr = ' f('+str("%5s" % el)+')   = '+str("%8.2f" % dat_lineratio['f_'+el][objent_lr])+\
                               ' @ S/N ='+str("%8.2f" % dat_lineratio['s2n_'+el][objent_lr])+\
                               '   with Dv = '+str("%8.2f" % dat_lineratio['vshift_'+el][objent_lr])+\
                               ' km/s,  sigma = '+str("%8.2f" % dat_lineratio['sigma_'+el][objent_lr]+' .')

                    if el.lower()  != 'heii':
                        FRstring = ' and double flux ratio = '+str("%8.2f" % dat_lineratio['FR_'+el+'1'+el+'2'][objent_lr]+' .')
                        printstr = printstr.replace(' .',FRstring)

                    if (dat_lineratio['s2n_'+el][objent_lr] > 3.0) & (np.abs(dat_lineratio['vshift_'+el][objent_lr]) < 1000.0):
                        print('\033[94m'+printstr+'\033[0m') # color potential detections blue
                    else:
                        print('\033[91m'+printstr+'\033[0m') # color low-S/N or high-vshift detections red
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if verbose: print('\n (Info from summaries; -99 = no data file provided; None = ID missing) ')
            outstr = str("%12s" % objid)+' '+str("%15s" % objpoint)+' '

            for emline in emlines:
                if verbose: print(' - Checking template matches to '+emline)
                answer = ''
                if objent_lr is None:
                    answer = 'idmissing'
                    outstr = outstr+str("%12s" % answerkeys[answer.lower()])+' '
                else:
                    lineS2N    = dat_lineratio['s2n_'+emline][objent_lr]
                    ELquestion = '   -> Do you trust FELIS LLLL (y)es/(n)o/(m)aybe/(e)exit? '.replace('LLLL',str("%5s" % emline))

                    if lineS2N < S2Nmincheck:
                        answer = 'lows2n'
                    elif ~np.isfinite(lineS2N):
                        answer = 'nocov'
                    elif np.abs(dat_lineratio['vshift_'+emline][objent_lr]) > 1000.0:
                        answer = 'n'
                    else:
                        lp_all = ''
                        for plotdir in plotdirs:
                            lp_all = lp_all+' '+' '.join(glob.glob(plotdir+'*'+str(objpoint)+'*'+
                                                                   str("%.9d" %  objid)+'*'+str(emline)+'*.pdf'))
                        lineplots = ''
                        for lp in lp_all.split():
                            if (lp.split('/')[-1]).startswith("overview_1DspecWzooms") == False:
                                lineplots = lineplots+' '+lp
                        if performinspections:
                            if lineplots != '':
                                pipe_lineplots = subprocess.Popen(opencommand+lineplots,shell=True,executable=os.environ["SHELL"])
                                time.sleep(1.1) # sleep to make sure process appears in PIDlist
                                pid_lineplots  = MiGs.getPID('Preview.app',verbose=False) # get PID of png process
                            else:
                                if verbose: print('   WARNING did not find a plot of the FELIS match to the line')

                            while answer.lower() not in ['y','n','m']:
                                if pversion == 2:
                                    answer = raw_input(ELquestion) # raw_input for python 2.X
                                elif pversion == 3:
                                    answer = input(ELquestion)     # input for python 3.X
                                else:
                                    sys.exit(' Unknown version of python: version = '+str(pversion))
                                if answer == 'e':
                                    fout.close()
                                    sys.exit('   Exiting as answer provided was "e". Vetting summarized in\n   '+outputfile)
                            if lineplots != '':
                                killsignal = 1
                                try:
                                    os.kill(pipe_lineplots.pid+1,killsignal)
                                except:
                                    print(' WARNING: Was unable to close lineplots ')
                        else:
                            answer = 'n' # no lines trusted when not performing inspections
                    outstr = outstr+str("%12i" % answerkeys[answer.lower()])+' '

            if performinspections:
                if pversion == 2:
                    notes  = raw_input('   -> Anything to add for this object?\n'
                                       '      (AGN? Marginal detections? Bad spec?)    ') # raw_input for python 2.X
                elif pversion == 3:
                    notes  = input('   -> Anything to add for this object?\n'
                                   '      (AGN? Marginal detections? Bad spec?)    ') # raw_input for python 3.X
                else:
                    sys.exit(' Unknown version of python: version = '+str(pversion))
            else:
                notes = ' Automatically generated entry from setting performinspections=True'

            if notes == 'e':
                fout.close()
                sys.exit('   Exiting as input for notes was "e". Vetting summarized in\n   '+outputfile)
            else:
                outstr = outstr+'  #Notes: '+notes+' \n'

            fout.write(outstr)
            fout.close()
            fout = open(outputfile,'a')

        if performinspections:
            if mainplots != '':
                killsignal = 1
                try:
                    os.kill(pid_mainplots,killsignal)
                except:
                    print(' WARNING: Was unable to close mainplots ')
    if verbose: print('\n--------- Done --------- ')
    if verbose: print(' Wrote output to:\n '+outputfile)
    fout.close()
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def return_objent(id,dataarray,idcol='id',verbose=True):
    """
    Little script to match and return object index in data array
    """
    if dataarray is not None:
        if type(id) is list: # list provided so ID and Pointing to be usef for selection
            objents = np.where((dataarray[idcol[0]].astype(int) == int(id[0])) &
                               (dataarray[idcol[1]] == str(id[1])))[0]
        else:
            objents = np.where(dataarray[idcol].astype(int) == int(id))[0]
        if len(objents) > 1:
            if verbose: print('   WARNING: Multiple matches ('+str(len(objents))+') to '+str(id)+' found by uves.return_objdat()')
        elif len(objents) == 0:
            if verbose: print('   WARNING: No matches to '+str(id)+' found in uves.return_objdat()')
            objents = None
    else:
        objents = -99

    return objents

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def summarize_tdosevetting(returnsample='udf10',verbose=True):
    """
    Summarizing the content of the vetting of the TDOSE spectra

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    idbad, idgood = uves.summarize_tdosevetting(returnsample='udf10')
    """

    vetdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/vet_tdose_extractions_outputs/'
    vetcat = vetdir+'MWuves-full-v1p0_vet_tdose_extractions_output_manuallycombined.txt'

    vetdat = np.genfromtxt(vetcat,names=True,comments='#',dtype=None,skip_header=28)

    selcuts = collections.OrderedDict()
    selcuts['all']        = [0.0,   1e12]
    selcuts['cdfs']       = [0.9e8, 1.9e8 ]
    selcuts['parallels']  = [2.9e8, 4.9e8  ]
    selcuts['cosmos']     = [1.9e8, 2.9e8  ]
    selcuts['udfmosaic']  = [5.9e8, 6.9e8  ]
    selcuts['udf10']      = [6.9e8, 7.9e8  ]

    if returnsample not in selcuts.keys():
        sys.exit(' - Invalid sample to return ("'+returnsample+'") provided to uves.summarize_tdosevetting(). '
                                                               '\n   Choose between: '+str(selcuts.keys()))

    # - - - -  count - - - -
    print(' ---- vetting results in '+vetcat.split('/')[-1]+' ----')
    for key in selcuts.keys():
        idmin, idmax = selcuts[key]
        Nobj  = len(np.where(( vetdat['id'].astype(float) > idmin) & ( vetdat['id'].astype(float) < idmax ) )[0] )
        zeros = np.where((vetdat['vetresult'] == 0) &
                         ( vetdat['id'].astype(float) > idmin) & ( vetdat['id'].astype(float) < idmax ) )[0]
        ones = np.where((vetdat['vetresult'] == 1) &
                        ( vetdat['id'].astype(float) > idmin) & ( vetdat['id'].astype(float) < idmax ) )[0]
        twos  = np.where((vetdat['vetresult'] == 2) &
                         ( vetdat['id'].astype(float) > idmin) & ( vetdat['id'].astype(float) < idmax ) )[0]
        threes = np.where((vetdat['vetresult'] == 3) &
                          ( vetdat['id'].astype(float) > idmin) & ( vetdat['id'].astype(float) < idmax ) )[0]
        fours = np.where((vetdat['vetresult'] == 4) &
                         ( vetdat['id'].astype(float) > idmin) & ( vetdat['id'].astype(float) < idmax ) )[0]

        ids_goodspec = vetdat['id'][np.append(ones,twos)]
        ids_badspec = vetdat['id'][np.append(threes,fours)]

        if Nobj>0:
            fracgood = float(len(ids_goodspec))/float(Nobj)
        else:
            fracgood = 0

        if verbose:
            print(' - The fraction of good spectra (vetresult = 1 or 2) among the '+
                  str("%5s" % Nobj)+' spectra from '+str("%10s" % key)+' is:        '+str(fracgood))
            print('   This corresponds to '+str(len(ids_badspec))+'/'+str(len(ids_badspec)+len(ids_goodspec))+' bad spectra')
            print('   The number of spectra marked for revisiting was: '+str(len(zeros)))

        if key == returnsample:
            ids_badspec_return, ids_goodspec_return = ids_badspec, ids_goodspec
            print('   >>>>>>>>>>>> Note: The '+key+' IDs were returned <<<<<<<<<<<<\n')
        else:
            print('   ')

    return ids_badspec_return, ids_goodspec_return
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def summarize_felisvetting(vetoutput,verbose=True):
    """
    Function summarizing the content of an output file from uves.vet_felisdetection()

    It also looks for close neighbors which can be used to clean out dublicate entries of objects in overlap regions.

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    vetoutput = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FELIStemplatematch2uvesobjects/all_aperture190926/vet_felisdetection_output_for_object_match_selection_all191029_aperture_s2nGT3_vshiftLT1000_191030vettingUDFonly.txt'
    uves.summarize_felisvetting(vetoutput)

    """
    fmt    = ['d','12a','d','d','d','d','d','d','d']
    vetdat = np.genfromtxt(vetoutput,skip_header=12,dtype=fmt,names=True,comments='#')

    lines = [nn.replace('trust','') for nn in vetdat.dtype.names[2:]]
    if verbose: print(' - Summarizing content of '+str(len(vetdat))+' entries in \n   '+vetoutput)
    if verbose: print('   (entries for '+str(len(np.unique(vetdat['id'])))+' unique ids)\n')
    for cc, col in enumerate(vetdat.dtype.names[2:]):
        N_cov   = float(len(np.where(vetdat[col] >= 0)[0]))
        N_yes   = float(len(np.where(vetdat[col] == 1)[0]))
        N_no    = float(len(np.where(vetdat[col] == 0)[0]))
        N_maybe = float(len(np.where(vetdat[col] == 9)[0]))
        N_lim   = float(len(np.where(vetdat[col] == 99)[0]))
        if verbose:
            print('  '+str("%5s" % lines[cc])+' (numbers)  : '+str("%8i" % N_yes)+' '+str("%8i" % N_no)+' '+
                  str("%8i" % N_maybe)+' '+str("%8i" % N_lim)+'   of '+str("%5i" % N_cov)+' entries')
            print('        (fraction) : '+str("%8.4f" % (N_yes/N_cov))+' '+
                  str("%8.4f" % (N_no/N_cov))+' '+str("%8.4f" % (N_maybe/N_cov))+' '+str("%8.4f" % (N_lim/N_cov)))

    N_systemic           = 0.0
    N_systemic_wCIV      = 0.0
    N_systemic_wCIVaMgII = 0.0
    id_sys               = []
    id_sys_wCIV          = []
    id_sys_wCIVaMgII     = []

    for ii, objid in enumerate(vetdat['id']):
        if  (vetdat['trustNV'][ii]    == 1) or (vetdat['trustHeII'][ii]  == 1) or (vetdat['trustOIII'][ii]  == 1) or (vetdat['trustSiIII'][ii] == 1) or (vetdat['trustCIII'][ii]  == 1):
            N_systemic = N_systemic + 1.0
            id_sys.append(int(objid))

        if  (vetdat['trustNV'][ii]    == 1) or (vetdat['trustCIV'][ii]   == 1) or  (vetdat['trustHeII'][ii]  == 1) or (vetdat['trustOIII'][ii]  == 1) or (vetdat['trustSiIII'][ii] == 1) or (vetdat['trustCIII'][ii]  == 1):
            N_systemic_wCIV = N_systemic_wCIV + 1.0
            id_sys_wCIV.append(int(objid))

        if  (vetdat['trustNV'][ii]    == 1) or (vetdat['trustCIV'][ii]   == 1) or (vetdat['trustHeII'][ii]  == 1) or (vetdat['trustOIII'][ii]  == 1) or (vetdat['trustSiIII'][ii] == 1) or (vetdat['trustCIII'][ii]  == 1) or (vetdat['trustMgII'][ii]  == 1):
            N_systemic_wCIVaMgII = N_systemic_wCIVaMgII + 1.0
            id_sys_wCIVaMgII.append(int(objid))

    if verbose:
        print('\n - Number of systemic redshifst = '+str(N_systemic)+'      (excl. CIV and MgII detections)')
        print(' - Number of systemic redshifst = '+str(N_systemic_wCIV)+'      (incl. CIV but excl. MgII detections)')
        print(' - Number of systemic redshifst = '+str(N_systemic_wCIVaMgII)+'      (incl. CIV and MgII detections)')

    coordsepsel = 0.5
    if verbose: print(' - Checking for dublicates/close neighbors within '+str(coordsepsel)+' arcsec:')
    infofile ='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo_UVemitters_3timesUDFcats.fits'
    datinfo  = afits.open(infofile)[1].data
    ralist   = []
    declist  = []
    zlist    = []

    for objid in np.unique(vetdat['id']):
        ralist.append( datinfo['ra'][np.where(datinfo['id'].astype(int) == int(objid))[0][0]] )
        declist.append( datinfo['dec'][np.where(datinfo['id'].astype(int) == int(objid))[0][0]] )
        zlist.append( datinfo['redshift'][np.where(datinfo['id'].astype(int) == int(objid))[0][0]] )

    for oo, objid in enumerate(np.unique(vetdat['id'])):
        coordiff_deg = np.sqrt( (np.cos(np.deg2rad(np.asarray(declist)))*(np.asarray(ralist)-np.asarray(ralist)[oo]))**2.0 +
                                (np.asarray(declist)-np.asarray(declist)[oo])**2.0 )
        coordiff = coordiff_deg * 3600.
        matchent = np.where((0 < coordiff) & (coordiff < coordsepsel))[0]
        if len(matchent) > 0:
            matchIDs = vetdat['id'][matchent]
            coorsep  = coordiff[matchent]
            if verbose: print('   matches to '+str(int(objid))+': '+str(matchIDs.astype(int))+'      @ sep = '+str("%.4f" % coorsep)+' arcsec')

    regionfile = vetoutput.replace('.txt','.reg')
    if regionfile != vetoutput:
        if verbose: print('\n - Storing region file to\n'+regionfile)
        colors = ['red']*len(ralist)
        for oo, objid in enumerate(vetdat['id'].astype(int)):
            if objid in id_sys_wCIVaMgII:
                colors[oo] = 'magenta'
            if objid in id_sys_wCIV:
                colors[oo] = 'orange'
            if objid in id_sys:
                colors[oo] = 'cyan'

        kbs.create_DS9region(regionfile,ralist,declist,color=colors,textlist=vetdat['id'].astype(str),clobber=True)

        ent_wLya = np.where(np.asarray(zlist) > 2.9)[0]
        regionfileLya = regionfile.replace('.reg','_zLT2p9.reg')
        if verbose: print('\n - Storing LAEs to region file \n'+regionfileLya)
        kbs.create_DS9region(regionfileLya,np.asarray(ralist)[ent_wLya],np.asarray(declist)[ent_wLya],
                             color=np.asarray(colors)[ent_wLya],textlist=vetdat['id'].astype(str)[ent_wLya]
                             ,circlesize=0.3,clobber=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_mastercat_v2(outputfits, file_info, file_fluxratio, file_EWestimates,
                       printwarning=True, overwrite=False, verbose=True):
    """
    Function assembling a master output catalog from the flux ratio measurements and EW estimates.
    Based on uves.build_mastercat() but now assumes that the TDOSE vetting has already been performed
    prior to FELIS vetting.

    Catalog also accounts for duplicate objects in the sence that that the information is included

    Ideally, this should be the catalog plotted, analyzed and released with any publication.

    This can also form the base of generating LaTeX tables

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    outputfits = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/results_master_catalog_testversion2002XX.fits'
    # file_info        = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    file_info        = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    file_fluxratio   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/back2backAnalysis_200213/fluxratios_FELISmatch2all_200213_preFELISvetting.txt'
    file_EWestimates = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/back2backAnalysis_200213/fluxratios_FELISmatch2all_200213_preFELISvetting_EW0estimates.txt'

    uves.build_mastercat_v2(outputfits,file_info,file_fluxratio,file_EWestimates)

    """
    if os.path.isfile(outputfits) and (overwrite == False):
        sys.exit('The output file '+outputfits+' already exists and overwrite=False ')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' --- Building the master catalog for the UVES study --- ')

    if verbose: print(' - Load infofile data ')
    dat_info         = afits.open(file_info)[1].data
    dat_info         = dat_info[np.where((dat_info['id']<4.9e8) | (dat_info['id']>5.9e8))[0]] # ignoring UDF MW mock ids

    if verbose: print(' - Load fluxratio data ')
    dat_fluxratio    = np.genfromtxt(file_fluxratio,    dtype=None,comments='#',names=True,skip_header=7)

    if verbose: print(' - Load EW estiamtes ')
    dat_EWestimates  = np.genfromtxt(file_EWestimates,  dtype=None,comments='#',names=True,skip_header=10)

    pointingIDdic            = {}
    pointingIDdic[106004019] = 'stackof01and06'
    pointingIDdic[109014056] = 'stackof04and09'
    pointingIDdic[121004014] = 'stackof20and21'
    pointingIDdic[122022112] = 'stackof16and22'
    pointingIDdic[124038073] = 'stackof24and18'
    pointingIDdic[130012014] = 'stackof11and30'
    pointingIDdic[131021114] = 'stackof31and41'
    pointingIDdic[133021057] = 'stackof33and28'
    pointingIDdic[139032271] = 'stackof39and40'
    pointingIDdic[143029102] = 'stackof43and45'
    pointingIDdic[154037129] = 'stackof54and50'
    pointingIDdic[158007026] = 'stackof57and58'
    pointingIDdic[159006015] = 'stackof02and59'
    pointingIDdic[202010025] = 'stackof02and03'
    pointingIDdic[215027065] = 'stackof08and15'

    file_specsel     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/' \
                       'MWuves-full-v1p0_spectra_paperselection200213.txt'
    dat_specsel      = np.genfromtxt(file_specsel,    dtype=None,comments='#',names=True,skip_header=4)
    pointing_sel     = [ss.split('/')[-1].split('-full-')[0].split('spectrum_')[-1] for ss in dat_specsel['spectrum']]
    for ss, selid in enumerate(dat_specsel['id']):
        if selid not in pointingIDdic.keys():
            pointingIDdic[selid] = pointing_sel[ss]
        else:
            print(' ID already in dictionary with pointings... that should not happen; stopping to investigate...')
            pdb.set_trace()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Setting up output structure. \n   Definiing columns: ')
    outcolnames_init = dat_fluxratio.dtype.names + dat_EWestimates.dtype.names[2:]

    outcolnames = []
    for outcolname in outcolnames_init: # looping over columns to ignore unnescessary entries (too keep below 1000 columns)
        if outcolname.startswith('FR') & (outcolname.split('_')[-1].endswith('1')):
            continue
        elif outcolname.startswith('FR') & (outcolname.split('_')[-1].endswith('2')):
            FRname = outcolname.split('_')[-1]
            if ('1' in FRname):
                if FRname.split('1')[0] == FRname.split('1')[-1][:-1]:
                    outcolnames.append(outcolname)
                else:
                    continue
        elif outcolname.startswith('FR') & ('1' in outcolname.split('_')[-1]):
            continue
        elif outcolname.startswith('FR') & ('2' in outcolname.split('_')[-1]):
            continue
        else:
            outcolnames.append(outcolname)

    print(outcolnames)

    Nrows       = len(dat_info)
    Ncols       = len(outcolnames)
    for cc, colname in enumerate(outcolnames):
        if verbose:
            infostr = '   '+str("%20s" % colname)+'   ( '+str("%.5d" % (cc+1))+' / '+str("%.5d" % Ncols)+' ) '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()
        try:
            npstr   = dat_fluxratio.dtype[colname].str
        except:
            npstr   = dat_EWestimates.dtype[colname].str

        if '|S' in npstr:
            fmt = npstr.split('S')[1]+'A'
            if colname == 'pointing': fmt = '30A'
            dat = np.asarray(['None'*Nrows])
        elif npstr == '<i8':
            fmt     = 'K'
            dat     = (np.zeros(Nrows)-99).astype(int)
        elif npstr == '<f8':
            fmt     = 'D'
            dat     = np.zeros(Nrows)*np.nan
        else:
            sys.exit(' No setup for the numpy array format "'+npstr+'"')

        col_def  = afits.ColDefs([afits.Column(name=colname, format=fmt, array=dat)])
        try:
            col_out = col_out + col_def
        except:
            col_out = col_def

    col_def  = afits.ColDefs([afits.Column(name='duplicationID', format='K', array=(np.zeros(Nrows)-99).astype(int))])
    col_out  = col_out[:2] + col_def + col_out[2:]

    col_def  = afits.ColDefs([afits.Column(name='redshift', format='D', array=dat_info['redshift'])])
    col_out  = col_out[:2] + col_def + col_out[2:]


    Ncols    = len(col_out)
    if verbose: print('\n - After adding additional columns the total number of columns to output are: '+str(Ncols))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('   Creating fits file HDU with column definitions ')
    hdu_out     = afits.BinTableHDU.from_columns(col_out)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Filling output with data')
    for ii, id in enumerate(dat_info['id']):
        if verbose:
            infostr = '   Checking data for id='+str(id)+' ('+str("%.5d" % (ii+1))+' / '+str("%.5d" % len(dat_info['id']))+')     '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        hdu_out.data['duplicationID'][ii] = dat_info['duplicationID'][ii]

        # ------------- Pulling out data from catalogs and writing it to master cat -------------
        objent_flux      = np.where(dat_fluxratio['id']   == id)[0]
        objent_ew        = np.where(dat_EWestimates['id'] == id)[0]

        if len(objent_flux) == 1:
            for colname in dat_fluxratio.dtype.names:
                if colname in outcolnames:
                    if colname == 'pointing':
                        hdu_out.data[colname][ii] = pointingIDdic[id]
                    else:
                        hdu_out.data[colname][ii] = dat_fluxratio[colname][objent_flux][0]

        if len(objent_ew) == 1:
            for colname in dat_EWestimates.dtype.names[2:]:
                if colname in outcolnames:
                    hdu_out.data[colname][ii] = dat_EWestimates[colname][objent_ew][0]

        if hdu_out.data['id'][ii] == -99: # checking that data for object was inserted
            hdu_out.data['id'][ii]       = id
            hdu_out.data['pointing'][ii] = 'NoDataToInsert'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Storing final catalog to '+outputfits)
    hdu_primary = afits.PrimaryHDU()
    hdulist     = afits.HDUList([hdu_primary, hdu_out])
    hdulist.writeto(outputfits, overwrite=overwrite)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_mastercat(outputfits, printwarning=True, overwrite=False, verbose=True):
    """
    --------------------------------------------------------------------------------------------------------
    ---> See uves.build_mastercat_v2() making this function obsolete for ouputs generated after 200218 <---
    --------------------------------------------------------------------------------------------------------

    Function assembling a master output catalog from the flux ratio measurements and EW estimates taking
    the results of the FELIS and TDOSE vetting into account.

    Ideally, this should be the catalog plotted, analyzed and released with any publication.

    This can also form the base of generating LaTeX tables

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    outputfits = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/results_master_catalog_testversion2002XX.fits'
    uves.build_mastercat(outputfits)

    """
    if os.path.isfile(outputfits) and (overwrite == False):
        sys.exit('The output file '+outputfits+' already exists and overwrite=False ')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' --- Building the master catalog for the UVES study --- ')
    dir_main         = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'

    # file_info        = dir_main+'objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    file_info        = dir_main+'objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    file_fluxratio   = dir_main+'FELIStemplatematch2uvesobjects/all_gauss190926/fluxratios/' \
                                'fluxratios_FELISmatch2uves190926_gauss_wpointings.txt'
    file_EWestimates = file_fluxratio.replace('.txt','_EW0estimates_191028run.txt')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Load infofile data ')
    dat_info         = afits.open(file_info)[1].data

    if verbose: print(' - Load fluxratio data ')
    dat_fluxratio    = np.genfromtxt(file_fluxratio,    dtype=None,comments='#',names=True,skip_header=7)

    if verbose: print(' - Load EW estiamtes ')
    dat_EWestimates  = np.genfromtxt(file_EWestimates,  dtype=None,comments='#',names=True,skip_header=8)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Setting up output structure. \n   Definiing columns: ')
    outcolnames = dat_fluxratio.dtype.names + dat_EWestimates.dtype.names[2:]
    Nrows       = len(dat_info)
    Ncols       = len(outcolnames)
    Nskip       = 0
    for cc, colname in enumerate(outcolnames):
        if verbose:
            infostr = '   '+str("%20s" % colname)+'   ( '+str("%.5d" % (cc+1))+' / '+str("%.5d" % Ncols)+' ) '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()
        try:
            npstr   = dat_fluxratio.dtype[colname].str
        except:
            npstr   = dat_EWestimates.dtype[colname].str

        if '|S' in npstr:
            fmt = npstr.split('S')[1]+'A'
            if colname == 'pointing': fmt = '30A'
            dat = np.asarray(['None'*Nrows])
        elif npstr == '<i8':
            fmt     = 'K'
            dat     = (np.zeros(Nrows)-99).astype(int)
        elif npstr == '<f8':
            fmt     = 'D'
            dat     = np.zeros(Nrows)*np.nan
        else:
            sys.exit(' No setup for the numpy array format "'+npstr+'"')

        if 'NV' in colname: # igrnoring columns incl. NV as this is not going to be sensible due the Lya either way;
                            # helps stay below Ncol = 999 which is a fundamental limitation of the fits format.
            Nskip = Nskip + 1
            continue
        else:
            col_def  = afits.ColDefs([afits.Column(name=colname, format=fmt, array=dat)])
            try:
                col_out = col_out + col_def
            except:
                col_out = col_def

    Ncols = Ncols - Nskip
    if verbose: print('\n   (Skipped the '+str(Nskip)+' rows containing NV so Ncols='+str(Ncols))

    if verbose: print('   Creating fits file HDU with column definitions ')
    hdu_out     = afits.BinTableHDU.from_columns(col_out)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Filling output with data adjusting for vetting results ')
    ids_badspec, ids_goodspec = uves.summarize_tdosevetting(returnsample='all')
    pointing_selector_dic     = uves.pointing_selector()

    for ii, id in enumerate(dat_info['id'][:]):
        if (len(str(id)) == 9) & (str(id).startswith('5')): # skipping UDF MWmock depth objects
            hdu_out.data['id'][ii]       = id
            hdu_out.data['pointing'][ii] = 'SkippingUDFMWmock'

            continue
        if verbose:
            infostr = '   Checking data for id='+str(id)+' ('+str("%.5d" % (ii+1))+' / '+str("%.5d" % len(dat_info['id']))+')     '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        if id in ids_badspec:
            hdu_out.data['id'][ii]       = id
            hdu_out.data['pointing'][ii] = 'TDOSEvetBadSpec'
            continue

        objent_flux      = np.where(dat_fluxratio['id']   == id)[0]
        objent_ew        = np.where(dat_EWestimates['id'] == id)[0]
        # ------------- check what pointing to use if more than 1 appearance of ID -------------
        if (len(objent_flux) > 1) or (len(objent_ew) > 1):
            if verbose & printwarning: print('\n   WARNING Multiple pointings for object id='+str(id))
            if str(id) in pointing_selector_dic.keys():
                objent_flux      = np.where((dat_fluxratio['id'] == id) &
                                            (dat_fluxratio['pointing'] == pointing_selector_dic[str(id)]) )[0]
                objent_ew        = np.where((dat_EWestimates['id'] == id) &
                                            (dat_EWestimates['pointing'] == pointing_selector_dic[str(id)]) )[0]
                objent_vet_felis = np.where((dat_vetfelis['id'] == id) &
                                            (dat_vetfelis['pointing'] == pointing_selector_dic[str(id)]) )[0]
            else:
                if verbose & printwarning:
                    print('   ID not found in dictionary from uves.pointing_selector().')
                    if len(objent_flux) > 1:
                        print('   Check the entries in flux ratio file to make a decision (id appears '+
                              str(len(objent_flux))+' times there)')
                    if len(objent_ew) > 1:
                        print('   Check the entries in flux ratio file to make a decision (id appears '+
                              str(len(objent_ew))+' times there)')
                    print('   Skipping object until add to uves.pointing_selector() dictionary ')

                hdu_out.data['id'][ii]       = id
                hdu_out.data['pointing'][ii] = 'MultiPointingNoDecision'
                continue

        # ------------- Pulling out data from catalogs and writing it to master cat -------------
        if len(objent_flux) == 1:
            for colname in dat_fluxratio.dtype.names:
                if 'NV' in colname:
                    continue
                else:
                    hdu_out.data[colname][ii] = dat_fluxratio[colname][objent_flux][0]

        if len(objent_ew) == 1:
            for colname in dat_EWestimates.dtype.names[2:]:
                if 'NV' in colname:
                    continue
                else:
                    hdu_out.data[colname][ii] = dat_EWestimates[colname][objent_ew][0]

        if hdu_out.data['id'][ii] == -99: # checking that data for object was inserted
            hdu_out.data['id'][ii]       = id
            hdu_out.data['pointing'][ii] = 'NoDataToInsert'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Storing final catalog to '+outputfits)
    hdu_primary = afits.PrimaryHDU()
    hdulist     = afits.HDUList([hdu_primary, hdu_out])
    hdulist.writeto(outputfits, overwrite=overwrite)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def pointing_selector():
    """
    Function returning a dictionary with the pointing name to use for analysis for all the objects (IDs) appearing
    in multiple pointings.

    Note: As of 191129 this does not account for overlaps between MUSE-Wide, UDF and UDF-10 as the IDs are
    different in these cases


    """
    pointing_selector_dic = {}
    #                    ['_____id_____'] = pointing                # bad pointings
    pointing_selector_dic['_____id_____'] = 'udf-01'                # bad pointings

    return pointing_selector_dic

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def prepare_reextractionPostVetting(verbose=True,printVetComment=False):
    """
    Function to prepare setupfiles (by editing exiting files) and collecting objects to perform re-extractions for
    based on the TDOSE vetting summary generated with tdose_utilities.vet_tdose_extractions() collected in
    parentdir+'vet_tdose_extractions_outputs/MWuves-full-v1p0_vet_tdose_extractions_output_manuallycombined.txt'

    --- Example of use ---
    import uvEmissionlineSearch as uves
    uves.prepare_reextractionPostVetting()

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading vetting results and grabbing list of original TDOSE setup file to modify.')
    parentdir     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/'
    # vetresults    = parentdir+'vet_tdose_extractions_outputs/MWuves-full-v1p0_vet_tdose_extractions_output_manuallycombined.txt'
    # vetresults    = parentdir+'vet_tdose_extractions_outputs/MWuves-full-v1p0_vet_tdose_re-extractions_output_manuallycombined.txt'
    vetresults    = parentdir+'vet_tdose_extractions_outputs/MWuves-full-v1p0_vet_tdose_re-extractions_2ndround_manuallycombined.txt'
    if verbose: print('   Vetting results from: '+vetresults)
    outdir        = parentdir+'tdose_reextractionPostVetting/'
    orig_setups   = glob.glob(parentdir+'tdose_setupfiles/*tdose_setupfile_MWuves*_gauss.txt')
    vetdat        = np.genfromtxt(vetresults,names=True,comments='#',dtype=None,skip_header=28)
    Nvets         = len(vetdat['id'])
    ent_reext     = np.where(vetdat['vetresult'] > 2)[0]
    id_reext      = vetdat['id'][ent_reext]
    #--------------------------------------------------------------------
    # speclistfile = parentdir+'MWuves-full-v1p0_spectra_paperselection.txt'
    # if os.path.isfile(speclistfile):
    #     speclistdat = np.genfromtxt(speclistfile,names=True,comments='#',dtype=None,skip_header=4)
    #     if verbose: print(' - Loading output from uves.print_specs to get list of objects with existing good specs')
    #     IDswGoodSpec = speclistdat['id']
    #
    #     for ii, id_re in enumerate(id_reext):
    #         if id_re in IDswGoodSpec:
    #             id_reext[ii]  = -99
    #             ent_reext[ii] = -99
    #
    #     id_reext  = id_reext[id_reext > -1]
    #     ent_reext = ent_reext[ent_reext > -1]
    # else:
    #     if verbose: print(' - Did not find uves.print_specs() list of objects with existing good specs; including all ids')
    #--------------------------------------------------------------------
    Nreext        = len(ent_reext)
    point_reext   = [spec.split('tdose_spectrum_')[-1].split('-full')[0] for spec in vetdat['spectrum'][ent_reext]]
    point_reext_u = np.unique(point_reext)
    if verbose: print(' - Found '+str(Nvets)+' objects that had been vetted; '+str(Nreext)+
                      ' of those should be re-extracted (vet results > 2)')
    if verbose: print('   These are spread out over '+str(len(point_reext_u))+' individual MUSE pointings.')

    if printVetComment:
        print('   The comments for these objects from the vetting file are:')
        print('   - - - - - - - - - - - - - - - - - - - - - - - - ')
        fvet      = open(vetresults,'r')
        for line in fvet.readlines():
            if line.split()[0] in id_reext.astype(str):
                notesstr = line.split('#Notes:')[-1].split('\n')[0]
                if ('point source' in notesstr.lower()):
                    print('\033[91m'+'  '+line.split()[0]+':  '+notesstr+'\033[0m')
                elif ('fov' in notesstr.lower()):
                    print('\033[94m'+'  '+line.split()[0]+':  '+notesstr+'\033[0m')
                else:
                    print('  '+line.split()[0]+':  '+notesstr)

        print('   - - - - - - - - - - - - - - - - - - - - - - - - ')
        fvet.close()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Initializing setup files for re-extractions based on original TDOSE setup files.')
    reext_setups  = []
    for orig_setup in orig_setups:
        point_setup = orig_setup.split('v1p0-')[-1].split('_gauss')[0].replace('candels-','')
        if point_setup in point_reext_u:
            reext_setup = outdir+'tdose_setupfiles/'+orig_setup.split('/')[-1].replace('.txt','_reext.txt')
            shutil.copyfile(orig_setup, reext_setup)
            reext_setups.append(reext_setup)
        else:
            continue
    if len(point_reext_u) != len(reext_setups):
        print(' ERROR: The number of unique pointings from the vetting results does not '
              'match the number of re-extraction setups initialized')
        pdb.set_trace()
    if verbose: print('   '+str(len(reext_setups))+' re-extraction setups ready for modifications.')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Looping over setupfiles and editing content')
    idlistcheck  = []
    for rr, reext_setup in enumerate(reext_setups):
        pointing = reext_setup.split('v1p0-')[-1].split('_gauss')[0].replace('candels-','')
        fin      = open(reext_setup,'r')
        tmpfile  = reext_setup.replace('.txt','_tmp.txt')
        ftmp     = open(tmpfile,'w')

        for line in fin.readlines():
            if line.startswith('models_directory'):
                line   = line.replace('tdose_models/','tdose_models_reext/')
            if line.startswith('spec1D_directory'):
                line   = line.replace('tdose_spectra/','tdose_spectra_reext/')

            if line.startswith('max_centroid_shift'):
                line   = line.replace('  10  ','  2  ')

            if line.startswith('sources_to_extract'):
                sourceextfile = line.split()[1]
                macloc     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/' \
                             'tdose_extraction_MWuves_100fields_maxdepth190808/tdose_sourcecatalogs/'
                origsource = np.genfromtxt(macloc+sourceextfile.split('/')[-1],names=True,dtype=None)
                idlist = uves.get_idlist_for_reext_pointing(id_reext,pointing,idmaster=origsource['id'])
                line   = line.replace(sourceextfile,str(idlist).replace(' ',''))
                idlistcheck = idlistcheck + list(idlist)

            if line.startswith('spec1D_name'):
                line   = line.replace('MWuves','MWuves_reext')

            ftmp.write(line)
        ftmp.close()
        fin.close()
        shutil.move(tmpfile,reext_setup)   # replacing original file with temporary file with edits

    if verbose: print(' - Chekcing that all IDs to re-extract were found in edited setups')
    idsinsetups = np.unique(np.asarray(idlistcheck))
    for id in id_reext:
        if id not in idsinsetups:
            print('   setups do not appear to contain the id '+str(id))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_idlist_for_reext_pointing(idsinput,pointing,idmaster=None):
    """
    Function returning the IDs in a given pointing provided a full list of IDs.

    """
    idsout = []
    if 'cdfs-' in pointing:
        pointno = pointing.split('cdfs-')[-1][:2]
        for idcheck in idsinput:
            if str(int(idcheck)).startswith('1'+pointno):
                idsout.append(idcheck)
    elif 'cosmos-' in pointing:
        pointno = pointing.split('cosmos-')[-1][:2]
        for idcheck in idsinput:
            if str(int(idcheck)).startswith('2'+pointno):
                idsout.append(idcheck)
    elif 'hudf09-1-' in pointing:
        pointno = pointing.split('hudf09-1-')[-1][:2]
        for idcheck in idsinput:
            if str(int(idcheck)).startswith('3'+pointno):
                idsout.append(idcheck)
    elif 'hudf09-2-' in pointing:
        pointno = pointing.split('hudf09-2-')[-1][:2]
        for idcheck in idsinput:
            if str(int(idcheck)).startswith('4'+pointno):
                idsout.append(idcheck)
    elif 'udf-0' in pointing:
        for idcheck in idsinput:
            if str(int(idcheck)).startswith('6') & (int(idcheck) in idmaster.astype(int)):
                idsout.append(idcheck)
    elif 'udf-1' in pointing:
        for idcheck in idsinput:
            if str(int(idcheck)).startswith('7'):
                idsout.append(idcheck)
    else:
        sys.exit(' Did not find a matching string in pointing '+pointing)

    return idsout
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def print_specs(outputfile, overwrite=True, verbose=True):
    """
    Small function to print ID and spec list from TDOSE vetting results

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    uves.print_specs('./outfile.txt')

    """
    parentdir     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/'
    vr_main       = parentdir+'vet_tdose_extractions_outputs/MWuves-full-v1p0_vet_tdose_extractions_output_manuallycombined.txt'
    vr_reext1     = parentdir+'vet_tdose_extractions_outputs/MWuves-full-v1p0_vet_tdose_re-extractions_output_manuallycombined.txt'
    vr_reext2     = parentdir+'vet_tdose_extractions_outputs/MWuves-full-v1p0_vet_tdose_re-extractions_2ndround_manuallycombined.txt'
    vr_reext3     = parentdir+'vet_tdose_extractions_outputs/MWuves-full-v1p0_vet_tdose_re-extractions_3rdround_manuallycombined.txt'
    vetresults    = [vr_main,vr_reext1,vr_reext2,vr_reext3]

    if (overwrite == False) & os.path.isfile(outputfile):
        sys.exit('The file '+outputfile+' already exists but overwrite=False ')

    fout = open(outputfile,'w')
    fout.write('# List of spectra to search for UV emission with FELIS. \n')
    fout.write('# List generated with uves.print_specs() on '+kbs.DandTstr2()+' based on the TDOSE vetting results in: \n')
    fout.write('# '+' and '.join(vetresults)+'\n')
    fout.write('# \n')
    fout.write('# id           spectrum \n')

    if verbose: print(' - Looping over vetting results and writing output to \n   '+outputfile)
    Nmultistack = 0
    stackids    = []
    for vetresult in vetresults:
        fout.write('### IDs and Spec from vetting results collected in: '+vetresult.split('/')[-1]+' ###\n')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - -
        vetinfo       = open(vetresult,'r')
        stackent      = []
        lineindex     = 0
        for line in vetinfo.readlines():
            if line.startswith('#'):
                continue
            else:
                if 'n6-' in line:
                    stackent.append(lineindex)
                    stackids.append(int(line.split()[0]))
            lineindex = lineindex+1
        Nmultistack = Nmultistack + len(stackent)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - -
        vetdat        = np.genfromtxt(vetresult,names=True,comments='#',dtype=None,skip_header=28)
        for ii, id in enumerate(vetdat['id']):
            if vetdat['vetresult'][ii] > 2:
                fout.write("# Needs re-extraction: "+str("%s" % id)+'    '+
                           vetdat['spectrum'][ii].replace('..','/store/data/musewide/TDOSE/')+' \n')
            elif ii in stackent:
                fout.write("# Replaced by multi-field stack: "+str("%s" % id)+'    '+
                           vetdat['spectrum'][ii].replace('..','/store/data/musewide/TDOSE/')+' \n')
            else:
                fout.write(str("%s" % id)+'    '+
                           vetdat['spectrum'][ii].replace('..','/store/data/musewide/TDOSE/')+' \n')
    fout.close()

    if verbose: print(' - Loading output and summarizing: ')
    outdat = np.genfromtxt(outputfile,names=True,comments='#',dtype=None,skip_header=4)
    Nobj   = len(outdat['id'])
    Nobj_u = len(np.unique(np.sort(outdat['id'])))
    if Nobj != Nobj_u:
        if verbose: print('    The following objects appear N times in output ')
        if verbose: print(' objid       N')
        for objid in np.unique(np.sort(outdat['id'])):
            objent = np.where(outdat['id'] == objid)[0]
            if len(objent) > 1:
                if verbose: print(str(objid)+'   '+str(len(objent)))
        if verbose: print(' -> But there should be no duplications... go fix. ')
    else:
        if verbose:
            print('   Found '+str(Nobj_u)+' unique objects in output')
            print('   To be compared with '+str(2197-Nmultistack)+' expected (2197 unique IDs minus '+
                  str(Nmultistack)+' multifield-stacks)')
            print('   (not accounted for duplicates UDF-CDFS-UDF10 as of 200206)')

    if verbose:
        print('\n - The following IDs were not found among the vetting results and hence, do not exist in output list')
        # infofile     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
        infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
        infodat       = afits.open(infofile)[1].data
        for infoid in infodat['id']:
            if (infoid not in np.unique(np.sort(outdat['id']))) & (~str(infoid).startswith('5')):
                if infoid in stackids:
                    print('   '+str(infoid)+'   (part of multi-field staced saple)')
                else:
                    print('   '+str(infoid))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def collectAndRenamArcheSpec(speclist='/store/data/musewide/TDOSE/MWuves100full/MWuves-full-v1p0_speclist12XXXX.txt',
                             outputdir='/store/data/musewide/TDOSE/MWuves100full/MWuvesSpecs20XXXX/',
                             namestring='20XXXXselection',verbose=True):
    """
    As the name says, this function collects and renames the spectra to search for UV emission lines.
    In other words it assembles the final sample including selecting pointings, including stacks,
    assembling re-extractions, etc.

    NB: This has to be run with copy-past on arche

    --- Example of use ---
    import uvEmissionlineSearch as uves
    uves.collectAndRenamArcheSpec(namestring='200121selection',outputdir='/store/data/musewide/TDOSE/MWuves100full/test200121/')

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    import sys, numpy as np, os, tdose_utilities as tu
    if os.path.isdir(outputdir) != True:
        sys.exit(' Did not find output directory '+outputdir)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    specdat  = np.genfromtxt(speclist,dtype=None,names=True,comments='#',skip_header=4)
    Nobjects = len(np.unique(specdat['id']))
    if Nobjects != 2197:
        answer   = 'no answer'
        print('WARNING There should be spectra associated with a total of 2197 objects (IDs). Only found '+str(Nobjects))
        question = '        Are you sure you want to proceed (y/n)? '
        pversion = sys.version_info[0]

        while str(answer).lower() not in ['y','n']:
            if pversion == 2:
                answer = raw_input(question) # raw_input for python 2.X
            elif pversion == 3:
                answer =     input(question) # input for python 3.X
            else:
                sys.exit(' Unknown version of python: version = '+str(pversion))
        if answer.lower() == 'n':
            sys.exit(' okay - then exited as requested')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for ss, spec in enumerate(specdat['spectrum']):
        newspec = outputdir+'tdose_spectrum_'+namestring+'_'+str("%.10d" % int(specdat['id'][ss]))+'.fits'
        if verbose:
            infostr = ' - Removing SOURCECUBE and renaming spectrum for id='+str("%.10d" % int(specdat['id'][ss]))+\
                      '  ('+str("%.5d" % (ss+1))+' / '+str("%.5d" % len(specdat['id']))+')     '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()
        try:
            tu.strip_extension_from_fitsfile(spec,outputdir,removeextension='SOURCECUBE',overwrite=False,verbose=False)
            os.rename(outputdir+spec.split('/')[-1],newspec)
        except:
            print('\nWARNING Attempt to copy '+spec+' to '+newspec+' failed.')
    if verbose: print('\n   ... done')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose:
        print(' - To tar up the directory and its content excecute: ')
        print('   bash> cd '+outputdir+'..')
        if outputdir.endswith('/'):
            tarname = outputdir.split('/')[-2]
        else:
            tarname = outputdir.split('/')[-1]
        print('   bash> tar -zcvf ./'+tarname+'.tar.gz '+outputdir.split('/')[-2]+'/tdose_spectrum_'+namestring+'*.fits')
        print('   kbs>  scp kasper@arche.aip.de:'+outputdir+'../'+tarname+'.tar.gz  /Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def stack_IndividualObjectsWithMultiSpec(plotstackoverview=True,verbose=True):
    """
    Function collecting spectra of individual objects with multiple good spectra and stacking them.
    Basing selection on TDOSE vetting comments containing n6-xx-yy notes.

    --- Example of use ---
    import uvEmissionlineSearch as uves
    uves.stack_IndividualObjectsWithMultiSpec(plotstackoverview=False)

    """
    parentdir     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/'
    outdir        = parentdir+'stacks1D_individual_objects/'
    specdir       = parentdir+'spectra_all/'
    vetresults    = parentdir+'vet_tdose_extractions_outputs/MWuves-full-v1p0_vet_tdose_extractions_output_manuallycombined.txt'
    vetdat        = np.genfromtxt(vetresults,names=True,comments='#',dtype=None,skip_header=28)
    Nvets         = len(vetdat['id'])

    if verbose: print(' - Looking through vetting results ('+vetresults.split('/')[-1]+') to find objects to stack')
    vetinfo       = open(vetresults,'r')
    stackids      = collections.OrderedDict()
    ll            = 0
    for line in vetinfo.readlines():
        if line.startswith('#'):
            continue
        else:
            if 'n6-' in line:
                objid       = vetdat['id'][ll]
                stackfields = line.split('n6-')[-1].split()[0].split('-')
                objspec     = specdir+'*/tdose_spectra_cubestripped/'+vetdat['spectrum'][ll].split('/')[-1]
                objspecid   = objspec.split('-full')[0][-2:]
                if verbose: print('   Stack '+str(objid)+' spectra in fields: '+str(stackfields))
                stackids[objid] = [glob.glob(objspec.replace(objspecid+'-full',sf+'-full'))[0] for sf in stackfields],stackfields
            ll = ll+1

    if verbose: print(' - Stacking spectra in observed frame and saving output to:\n   '+outdir)
    # infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    infodat       = afits.open(infofile)[1].data

    for oo, objid in enumerate(stackids.keys()):
        spectra, fields = stackids[objid]
        outfile         = outdir+'tdose_spectrum_stackInclFields'+'-'.join(fields)+'_'+str("%.10d" % int(objid))+'.fits'
        wavelengths     = []
        fluxes          = []
        variances       = []

        for spectrum in spectra:
            data        = afits.open(spectrum)[1].data
            wavelengths.append(data['wave'])
            fluxes.append(data['flux'])
            variances.append(data['fluxerror']**2.0)

        if verbose: print('   Generating '+outfile.split('/')[-1]+'   (stack '+str(oo+1)+'/'+str(len(stackids.keys()))+')')
        stacktype ='mean'   # 'median'
        errtype   ='varsum'
        wave_out, flux_out, variance_out, Nspecstack = \
            stacking.stack_1D(wavelengths, fluxes, variances, z_systemic=np.zeros(len(spectra)),
                              stacktype=stacktype, errtype=errtype, wavemin=4750, wavemax=9350, Nsigmaclip=None,
                              deltawave=1.25, outfile=outfile, verbose=False)

        if plotstackoverview:
            plotspecs    = spectra+[outfile]
            labels       = ['Spec from field '+str(ff) for ff in fields]+['Stack of individual spec']
            wavecols     = ['wave']*len(plotspecs)
            fluxcols     = ['flux']*len(plotspecs)
            fluxerrcols  = ['fluxerror']*len(plotspecs)

            infoent      = np.where(infodat['id'] == objid)[0]
            zLya         = infodat['redshift'][infoent]
            voffset      = 250.0
            skyspectra   = [None]*len(spectra)+['/Users/kschmidt/work/MUSE/spectra_sky/SKY_SPECTRUM_candels-cdfs-06_av.fits']
            wavecols_sky = [None]*len(spectra)+['lambda']
            fluxcols_sky = [None]*len(spectra)+['data']
            yrangefull   = [-300,1000]
            xrangefull   = [4700,9400]

            for plotSN in [True,False]:
                mwp.plot_1DspecOverview(plotspecs, labels, wavecols, fluxcols, fluxerrcols, zLya[0], voffset=voffset,
                                        skyspectra=skyspectra, wavecols_sky=wavecols_sky, fluxcols_sky=fluxcols_sky,
                                        outputfigure=outfile.replace('.fits','_overview.pdf'),
                                        yrangefull=yrangefull, xrangefull=xrangefull,
                                        plotSN=plotSN,verbose=False)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def stack_composites_generate_setup(outputfile,equalsizebins=False,overwrite=False,verbose=True):
    """
    Function automatically putting together a file containing the setups for generating composites
    needed by

    --- Example of use ---
    import uvEmissionlineSearch as uves
    outputfile = './stacks1D_sample_selection.txt'
    outarray   = uves.stack_composites_generate_setup(outputfile,overwrite=True)

    """
    if os.path.isfile(outputfile) & (overwrite == False):
        sys.exit(' Overwrite was set to "False" and found existing copy of the file \n '+outputfile)

    if verbose: print(' - Loading infofile to enable cutting composites ')
    # infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    infodat       = afits.open(infofile)[1].data
    infodat       = infodat[np.where((infodat['id']<4.9e8) | (infodat['id']>5.9e8))[0]] # ignoring UDF MW mock ids

    columns       = ['zmin','zmax',
                     'Llyamin','Llyamax',
                     'FWHMlyamin','FWHMlyamax',
                     'm814wmin','m814wmax',
                     'EW0lyamin','EW0lyamax',
                     'betamin','betamax']

    outarray    = np.array([],dtype=[('id', 'i4'), ('label', 'U50'), ('nspec', 'i4'), ('ztype', 'U5'), ('bintype', 'U12')]+
                                    [(cn,'>f8') for cn in columns])
    emptyval    = 9999
    for equalsizebins in [False,True]:
        if not equalsizebins:
            eqbinstring = 'EqNumberBins'
            binid       = 100000
        else:
            eqbinstring = 'EqSizeBins'
            binid       = 200000

        for ztype in ['zcat','zv18']:
            if verbose: print('----------------- Generate setups for '+eqbinstring+' and '+ztype+' -----------------')
            coltrans  = uves.stack_composite_col_translator(ztype)
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if verbose: print(' - Setup for composite cintaining all spectra')
            if ztype == 'zcat':
                specid   = 10001+binid
            elif ztype == 'zv18':
                specid   = 20001+binid
            else:
                sys.exit(' Invalid entry for ztype (='+ztype+')')

            outrow   = np.asarray((specid,'all',99,ztype,eqbinstring)+
                                  (-emptyval,emptyval)*(len(columns)/2),dtype=outarray.dtype)
            outarray = np.append(outarray,outrow)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if verbose: print(' - Generating z-binning setups')
            for Nbins in [4,10]:
                for colbase in ['z']:
                    binranges = uves.get_vector_intervals(infodat[coltrans['zmin']][infodat[coltrans['zmin']]>2.9],
                                                          Nbins,verbose=False,
                                                          equalsizebins=equalsizebins)
                    for bb, br in enumerate(binranges):
                        specid   = specid+1
                        label    = 'zGT2p9_bin'+str(bb+1)+'of'+str(Nbins)
                        outrow   = np.asarray((specid,label,99,ztype,eqbinstring)+
                                              (-emptyval,emptyval)*(len(columns)/2),dtype=outarray.dtype)
                        outrow[colbase+'min'] = br[0]
                        outrow[colbase+'max'] = br[1]
                        outarray = np.append(outarray,outrow)

            Nbins     = 3
            for colbase in ['z']:
                if ztype != 'zv18':
                    binranges = uves.get_vector_intervals(infodat[coltrans['zmin']][infodat[coltrans['zmin']]<2.9],
                                                          Nbins,verbose=False,equalsizebins=equalsizebins)
                    for bb, br in enumerate(binranges):
                        specid   = specid+1
                        label    = 'zLT2p9_bin'+str(bb+1)+'of'+str(Nbins)
                        outrow   = np.asarray((specid,label,99,ztype,eqbinstring)+
                                              (-emptyval,emptyval)*(len(columns)/2),dtype=outarray.dtype)
                        outrow[colbase+'min'] = br[0]
                        outrow[colbase+'max'] = br[1]
                        outarray = np.append(outarray,outrow)
                else:
                    for bb, br in enumerate(np.arange(Nbins)):
                        specid   = specid+1
                        label    = 'zLT2p9_bin'+str(bb+1)+'of'+str(Nbins)
                        outrow   = np.asarray((specid,label,99,ztype,eqbinstring)+
                                              (np.nan,np.nan)*(len(columns)/2),dtype=outarray.dtype)
                        outarray = np.append(outarray,outrow)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if verbose: print(' - Generating m814w binning (seperating non-detections out) ')
            for Nbins in [4,10]:
                for colbase in ['m814w']:
                    datvec    = infodat[coltrans[colbase+'min']][infodat[coltrans[colbase+'min']]<29.39]
                    binranges = uves.get_vector_intervals(datvec,Nbins,verbose=False,equalsizebins=equalsizebins)

                    for bb, br in enumerate(binranges):
                        specid   = specid+1
                        label    = colbase+'_bin'+str(bb+1)+'of'+str(Nbins)
                        outrow   = np.asarray((specid,label,99,ztype,eqbinstring)+
                                              (-emptyval,emptyval)*(len(columns)/2),dtype=outarray.dtype)
                        outrow[colbase+'min'] = br[0]
                        outrow[colbase+'max'] = br[1]
                        outarray = np.append(outarray,outrow)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if verbose: print(' - Selection of objects with m814w non-edetections (limiting mag) ')
            datvec    = infodat[coltrans['m814wmin']][infodat[coltrans['m814wmin']]>29.39]
            binranges = uves.get_vector_intervals(datvec,1,verbose=True,equalsizebins=equalsizebins)
            specid    = specid+1
            label     = colbase+'nondet_bin1of1'
            outrow    = np.asarray((specid,label,99,ztype,eqbinstring)+
                                   (-emptyval,emptyval)*(len(columns)/2),dtype=outarray.dtype)
            outrow[colbase+'min'] = binranges[0][0]-0.0001 # manually expand range as all values are identical
            outrow[colbase+'max'] = binranges[0][1]+0.0001 # manually expand range as all values are identical
            outarray = np.append(outarray,outrow)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if verbose: print(' - Generating 4-bin selections (cut on just one parameter) ')
            coltrans  = uves.stack_composite_col_translator(ztype)
            Nbins     = 4
            for colbase in ['Llya','FWHMlya','EW0lya','beta']:
                binranges = uves.get_vector_intervals(infodat[coltrans[colbase+'min']],Nbins,verbose=False,
                                                      equalsizebins=equalsizebins)
                for bb, br in enumerate(binranges):
                    specid   = specid+1
                    label    = colbase+'_bin'+str(bb+1)+'of'+str(Nbins)
                    outrow   = np.asarray((specid,label,99,ztype,eqbinstring)+
                                          (-emptyval,emptyval)*(len(columns)/2),dtype=outarray.dtype)
                    outrow[colbase+'min'] = br[0]
                    outrow[colbase+'max'] = br[1]
                    outarray = np.append(outarray,outrow)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if verbose: print(' - Generating 3-bin selections (cut on two parameters) ')
            coltrans  = uves.stack_composite_col_translator(ztype)
            Nbins     = 3

            for colbase1,colbase2 in combinations( ['z','m814w','Llya','FWHMlya','EW0lya','beta'],2):
                if colbase1 == 'z':
                    datvec1 = infodat[coltrans['zmin']][infodat[coltrans['zmin']]>2.9]
                elif colbase1 == 'm814w':
                    datvec1 = infodat[coltrans[colbase1+'min']][infodat[coltrans[colbase1+'min']]<29.39]
                else:
                    datvec1 = infodat[coltrans[colbase1+'min']]
                binranges1 = uves.get_vector_intervals(datvec1,Nbins,verbose=False,equalsizebins=equalsizebins)

                if colbase2 == 'z':
                    datvec2 = infodat[coltrans['zmin']][infodat[coltrans['zmin']]>2.9]
                elif colbase2 == 'm814w':
                    datvec2 = infodat[coltrans[colbase2+'min']][infodat[coltrans[colbase2+'min']]<29.39]
                else:
                    datvec2 = infodat[coltrans[colbase2+'min']]
                binranges2 = uves.get_vector_intervals(datvec2,Nbins,verbose=False,equalsizebins=equalsizebins)

                for bb1, br1 in enumerate(binranges1):
                    for bb2, br2 in enumerate(binranges2):
                        specid   = specid+1
                        label    = colbase1+'_bin'+str(bb1+1)+'of'+str(Nbins)+'_'+colbase2+'_bin'+str(bb2+1)+'of'+str(Nbins)
                        outrow   = np.asarray((specid,label,99,ztype,eqbinstring)+
                                              (-emptyval,emptyval)*(len(columns)/2),dtype=outarray.dtype)
                        outrow[colbase1+'min'] = br1[0]
                        outrow[colbase1+'max'] = br1[1]
                        outrow[colbase2+'min'] = br2[0]
                        outrow[colbase2+'max'] = br2[1]
                        outarray = np.append(outarray,outrow)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('----------------- Counting objects for each setup and generating output -----------------')
    Nrows        = outarray.shape[0]
    # ppp = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/'
    # ccc = ppp + 'stacks1D_sample_selection.txt'
    # sss = np.genfromtxt(ccc,names=True,skip_header=2,comments='#',dtype=None)
    # print(' comparison: outarray VS fileload')
    if os.path.isfile(outputfile):
        dat_outfile     = np.genfromtxt(outputfile,names=True,skip_header=2,comments='#',dtype=None)

    for rr in np.arange(Nrows):
        if os.path.isfile(outputfile):
            ent_sample_file = uves.stack_composites_objselection(dat_outfile[rr],infodat)
            outarray[rr][2] = len(ent_sample_file)
        else:
            ent_sample = uves.stack_composites_objselection(outarray[rr],infodat)
            outarray[rr][2] = len(ent_sample)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Writing output to \n   '+outputfile)
    fout = open(outputfile, 'w')
    fout.write('# Setup file for generating coposite spectra with uves.stack_composites() \n')
    fout.write('# Created with uves.stack_composites_generate_setup() on '+kbs.DandTstr2()+' \n')
    fout.write('#  id                                             label                '
               'nspec              ztype              bintype'+
               ' '.join([str("%14s" % cc) for cc in columns])+'\n')
    for rr in np.arange(Nrows):
        outstr = str('%.5d' % outarray[rr][0])+\
                 str('%50s' % outarray[rr][1])+\
                 str('%20s' % outarray[rr][2])+\
                 str('%20s' % outarray[rr][3])+\
                 str('%20s' % outarray[rr][4])+\
                 ' '.join([str('%14.4f' % val) for val in outarray[rr].tolist()[5:]])
        if '     all     ' in outstr:
            fout.write('#'+''.join(['-------------------']*16)+'\n')
        fout.write(outstr+'\n')

    fout.close()
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def stack_composites(compositesetup,plotstackoverview=True,inloopload=False,verbose=True):
    """
    Function collecting spectra and stacking them based on various collection cuts defined in "stackdefs"

    --- Example of use ---
    import uvEmissionlineSearch as uves

    parentdir       = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/'
    compositesetups = parentdir+'stacks1D_sample_selection_manual.txt'
    uves.stack_composites(compositesetups,plotstackoverview=False)

    """
    parentdir     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/'
    outdir        = parentdir+compositesetup.split('/')[-1].split('.tx')[0]+'/'
    specdir       = parentdir+'MWuves-full-v1p0_spectra_paperselection200213/'
    # infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    infodat       = afits.open(infofile)[1].data
    infodat       = infodat[np.where((infodat['id']<4.9e8) | (infodat['id']>5.9e8))[0]] # ignoring UDF_MWmock

    stackinfo     = np.genfromtxt(compositesetup,names=True,skip_header=2,comments='#',dtype=None)
    Nstacks       = len(stackinfo['id'])

    if not inloopload:
        specdatadic   = {}
        if verbose: print(' - Assembling master dictionary with data from spectra (instead of loading in-loop)')
        for ii, id in enumerate(infodat['id']):
            searchstr = specdir+'*_0'+str(id)+'.fits'
            objspecs  = glob.glob(searchstr)
            if len(objspecs) != 1:
                print('WARNING Found '+str(len(objspecs))+' spectra (expected 1) globbing for:\n        '+searchstr)
                specdatadic[str(id)] = np.nan, np.nan, np.nan
            else:
                data                     = afits.open(objspecs[0])[1].data
                specdatadic[objspecs[0]] = data['wave'], data['flux'], data['fluxerror']**2.0

    if verbose: print(' - Stacking spectra in observed frame and saving output to:\n   '+outdir)
    for ii, stackid in enumerate(stackinfo['id']):
        outfile         = outdir+'composite_spectrum_'+str("%.10d" % stackid)+'_'+stackinfo['label'][ii]+'.fits'
        if verbose: print('   === Generating '+outfile.split('/')[-1]+'   (stack '+str(ii+1)+'/'+str(Nstacks)+') === ')
        wavelengths     = []
        fluxes          = []
        variances       = []

        if verbose: print('       > selecting indexes')
        ent_sample      = uves.stack_composites_objselection(stackinfo[ii],infodat)
        spectra         = []
        for ee, ent_s in enumerate(ent_sample):
            searchstr = specdir+'*_0'+str(infodat['id'][ent_s])+'.fits'
            objspecs  = glob.glob(searchstr)
            if len(objspecs) != 1:
                # print('WARNING Found '+str(len(objspecs))+' spectra (expected 1) globbing for:\n        '+searchstr)
                ent_sample[ee] = -99
            else:
                spectra.append(objspecs[0])
        spectra   = np.asarray(spectra)
        Nspec     = len(spectra)
        ent_sample = ent_sample[ent_sample != -99]

        if stackinfo[ii]['ztype'].lower() == 'zcat':
            objredshift  = infodat['redshift'][ent_sample.astype(int)]
        elif stackinfo[ii]['ztype'].lower() == 'zv18':
            objredshift  = infodat['z_sys_V18'][ent_sample.astype(int)]
        else:
            sys.exit(' Invalid entry for ztype (='+stackinfo[ii]['ztype']+') in '+compositesetup)

        # only passing on objects with determined redshifts to stacking (other-wise they cant be moved to rest-frame)
        spectra     = spectra[objredshift > 0]
        objredshift = objredshift[objredshift > 0]
        noz_ent     = np.where(objredshift <= 0)[0]
        if len(noz_ent) > 0:
            print('   Objects removed because they have no good redshift: \n   '+str(infodat['id'][noz_ent]))
        Nspec       = len(spectra)

        if verbose: print('       > building input structures for stacking the selected '+str(Nspec)+' objects ')
        for spectrum in spectra:
            if inloopload:
                data        = afits.open(spectrum)[1].data
                wavelengths.append(data['wave'])
                fluxes.append(data['flux'])
                variances.append(data['fluxerror']**2.0)
            else:
                datwave, datflux, datvariance = specdatadic[spectrum]
                wavelengths.append(datwave)
                fluxes.append(datflux)
                variances.append(datvariance)

        if verbose: print('       > performing stacking of the selected '+str(Nspec)+' objects: ')

        stacktype ='mean'   # 'median'
        errtypes  =['varsum']#,'stdonmean','fspread','std','medianvar']
        for errtype in errtypes:
            if verbose: print('          > estimating errors with the method "'+errtype+'" ')
            outfile_errappend = outfile.replace('.fits','_ERR'+errtype+'.fits')
            wave_out, flux_out, variance_out, Nspecstack = \
                stacking.stack_1D(wavelengths, fluxes, variances, z_systemic=objredshift, Nsigmaclip=1,
                                  stacktype=stacktype, errtype=errtype, wavemin=600, wavemax=3800,
                                  deltawave=0.1, outfile=outfile_errappend, verbose=False)

            if plotstackoverview:
                if verbose: print('            plotting the result ')
                plotspecs    = [outfile_errappend]
                labels       = [stackinfo['label'][ii].replace('_','\_')+' (stack of '+str(Nspec)+' spectra)']
                wavecols     = ['wave']
                fluxcols     = ['flux']
                fluxerrcols  = ['fluxerror']

                plotz        = 0.0
                voffset      = 0.0
                skyspectra   = [None]*len(plotspecs)
                wavecols_sky = [None]*len(plotspecs)
                fluxcols_sky = [None]*len(plotspecs)
                yrangefull   = [-30,100]
                xrangefull   = [600,3800]

                mwp.plot_1DspecOverview(plotspecs, labels, wavecols, fluxcols, fluxerrcols, plotz, voffset=voffset,
                                        skyspectra=skyspectra, wavecols_sky=wavecols_sky, fluxcols_sky=fluxcols_sky,
                                        outputfigure=plotspecs[0].replace('.fits','_overview.pdf'),
                                        yrangefull=yrangefull, xrangefull=xrangefull,
                                        show_error=False,verbose=False,plotSN=False)
                mwp.plot_1DspecOverview(plotspecs, labels, wavecols, fluxcols, fluxerrcols, plotz, voffset=voffset,
                                        skyspectra=skyspectra, wavecols_sky=wavecols_sky, fluxcols_sky=fluxcols_sky,
                                        outputfigure=plotspecs[0].replace('.fits','_overview.pdf'),
                                        yrangefull=yrangefull, xrangefull=xrangefull,
                                        show_error=False,verbose=False,plotSN=True)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def stack_composites_objselection(selectinfo,selectdata,verbose=True):
    """
    Function returning the entries of a given composite stack selection based on data infofiles

    --- Example of use ---
    see uves.stack_composites above

    """
    coltranslationdic = uves.stack_composite_col_translator(selectinfo['ztype'].lower())

    indexlist = np.arange(len(selectdata)) # all entries returned by default
    colsets   = np.unique(np.asarray([colname[:-3] for colname in selectinfo.dtype.names[5:]]))

    for colset in colsets:
        mincol = colset+'min'
        maxcol = colset+'max'
        if (selectinfo[mincol] != -9999) & (selectinfo[maxcol] != 9999):
            indexsel   = np.where((selectdata[coltranslationdic[mincol]] >= selectinfo[mincol]) &
                                  (selectdata[coltranslationdic[maxcol]] <= selectinfo[maxcol]) &
                                  (selectdata[coltranslationdic[mincol]] != 0.0) &
                                  (selectdata[coltranslationdic['zmin']] != 0.0))[0]
        else:
            indexsel   = np.where((selectdata[coltranslationdic['zmin']] != 0.0))[0]

        indexlist  = np.asarray([ii for ii in indexlist if ii in indexsel])
    return indexlist

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def stack_composite_col_translator(ztype,verbose=True):
    """

    coltranslationdic = uves.stack_composite_col_translator('zcat')

    """
    if ztype.lower() == 'zcat':
        zcol = 'redshift'
    elif ztype.lower() == 'zv18':
        zcol = 'z_sys_V18'
    else:
        print('\nWARNING: uves.stack_composite_col_translator() got unknown redshift type: '+ztype+'\n')
        zcol = None

    coltranslationdic = {}
    coltranslationdic['zmin']       = zcol
    coltranslationdic['zmax']       = zcol
    coltranslationdic['Llyamin']    = 'F_3KRON'
    coltranslationdic['Llyamax']    = 'F_3KRON'
    coltranslationdic['FWHMlyamin'] = 'fwhm_a_jk'
    coltranslationdic['FWHMlyamax'] = 'fwhm_a_jk'
    coltranslationdic['m814wmin']   = 'mag_acs_814w'
    coltranslationdic['m814wmax']   = 'mag_acs_814w'
    coltranslationdic['EW0lyamin']  = 'EW_0_beta_linear_many'
    coltranslationdic['EW0lyamax']  = 'EW_0_beta_linear_many'
    coltranslationdic['betamin']    = 'beta_linear_many'
    coltranslationdic['betamax']    = 'beta_linear_many'

    return coltranslationdic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def stack_composite_plotNxNspecs(param1,param2,param1range,param2range,spectra,outname,
                                 ztype='zcat',bintype='EqNumberBins',verbose=True):
    """
    function plotting the composite spectra from the binning of the objects displaying the binning
    in a seperate panel of the overview plot.

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves, glob
    specdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/stacks1D_sample_selection/'

    # --- two parameters binned ---
    param1  = 'z'
    param2  = 'm814w'
    globstr = specdir+'*'+param1+'_bin*'+param2+'_bin*.fits'
    spectra = np.sort(glob.glob(globstr))

    outname = specdir+'testoverviewWcols.pdf'
    uves.stack_composite_plotNxNspecs(param1,param2,[2.9,6.7],[23.79,28.82],spectra,outname)

    # --- one parameter binned ---
    param1  = 'z'
    globstr = specdir+'*'+param1+'LT2p9_bin*of3.fits'
    spectra = np.sort(glob.glob(globstr))
    outname = specdir+'composite_overview_'+param1+'LT2p9.pdf'
    uves.stack_composite_plotNxNspecs(param1,None,[1.5,2.9],None,spectra,outname)

    """
    parentdir      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/'
    compositesetup = parentdir+'stacks1D_sample_selection.txt'
    setupinfo      = np.genfromtxt(compositesetup,names=True,skip_header=2,comments='#',dtype=None)

    Nspec   = len(spectra)
    setuplabels    = [(param1+ss.split('.fit')[0].split('_'+param1)[-1].split('_ERR')[0]) for ss in spectra]

    # plot overview figure with 9 spectra in
    if verbose: print(' - Plotting the '+str(Nspec)+' spectra provided \n   ')
    plotspecs        = spectra

    if len(spectra) == 9:
        red1, red2, red3      = (205/255.,0/255.,0/255.), (255/255.,76/255.,76/255.), (255/255.,153/255.,153/255.)
        green1, green2, gree3 = (0/255.,205/255.,0/255.), (76/255.,255/255.,76/255.), (153/255.,255/255.,153/255.)
        blue1, blue2, blue3   = (0/255.,0/255.,205/255.), (76/255.,76/255.,255/255.), (153/255.,153/255.,255/255.)
        speccolors            = [blue1, blue2, blue3, green1, green2, gree3, red1, red2, red3]
        # speccolors       = ['orange','blue','red','green','brown','cyan','pink','purple','yellow']
    elif len(spectra) == 3:
        speccolors            = ['blue','green','red']
    elif len(spectra) == 4:
        speccolors            = ['blue','green','orange','red']
    elif len(spectra) == 10:
        cmap         = plt.cm.viridis
        cmap_norm    = plt.Normalize(vmin=1.0, vmax=Nspec)
        speccolors   = [cmap(cmap_norm(cc))[0:3] for cc in np.arange(Nspec)+1]

    labels           = [' ']*len(spectra) #[ll.replace('_','\_') for ll in setuplabels]
    wavecols         = ['wave']*len(spectra)
    fluxcols         = ['flux']*len(spectra)
    fluxerrcols      = ['fluxerror']*len(spectra)
    col_matrix       = True
    col_matrix_title = '' #'The Color Matrix'
    # col_matrix_text  =  [ss.split('.fit')[0].split(param1+'_')[-1].replace('bin','').replace('of','/').replace('_'+param2+'_','-')
    #                      for ss in spectra]
    setupinfo_ents   = np.unique(np.asarray( [np.where((setupinfo['label'] == ll) &
                                                       (setupinfo['ztype'] == ztype) &
                                                       (setupinfo['bintype'] == bintype))[0] for ll in setuplabels] ))
    col_matrix_text  = [str(setupinfo[np.where((setupinfo['label'] == ll) &
                                               (setupinfo['ztype'] == ztype) &
                                               (setupinfo['bintype'] == bintype))[0]]['nspec'][0]) for ll in setuplabels]
    col_matrix_labels=[param1,param2]
    if param2 is None:
        col_matrix_labels=[param1,'']
    col_matrix_ranges=[param1range,param2range]

    transdic      = uves.stack_composite_col_translator(ztype)
    # infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    infodat       = afits.open(infofile)[1].data
    infodat       = infodat[np.where((infodat['id']<4.9e8) | (infodat['id']>5.9e8))[0]] # ignoring UDF_MWmock

    if param2 is not None:
        datent        = np.where((infodat[transdic[param1+'min']] >= param1range[0]) &
                                 (infodat[transdic[param1+'max']] <= param1range[1]) &
                                 (infodat[transdic[param1+'max']] != 0) &
                                 (infodat[transdic[param2+'min']] >= param2range[0]) &
                                 (infodat[transdic[param2+'max']] <= param2range[1]) &
                                 (infodat[transdic[param2+'max']] != 0))[0]

        col_matrix_p1dat = infodat[transdic[param1+'min']][datent]
        col_matrix_p2dat = infodat[transdic[param2+'min']][datent]

        minvals1     = np.unique(setupinfo[setupinfo_ents][param1+'min'])
        maxvals1     = np.unique(setupinfo[setupinfo_ents][param1+'max'])
        param1ranges = [[minvals1[ii],maxvals1[ii]] for ii in [0,1,2]]

        minvals2     = np.unique(setupinfo[setupinfo_ents][param2+'min'])
        maxvals2     = np.unique(setupinfo[setupinfo_ents][param2+'max'])
        param2ranges = [[minvals2[ii],maxvals2[ii]] for ii in [0,1,2]]
    else:
        datent        = np.where((infodat[transdic[param1+'min']] >= param1range[0]) &
                                 (infodat[transdic[param1+'max']] <= param1range[1]) &
                                 (infodat[transdic[param1+'max']] != 0.0))[0]

        col_matrix_p1dat = infodat[transdic[param1+'min']][datent]
        col_matrix_p2dat = None

        minvals1     = np.unique(setupinfo[setupinfo_ents][param1+'min'])
        maxvals1     = np.unique(setupinfo[setupinfo_ents][param1+'max'])
        param1ranges = [[minvals1[ii],maxvals1[ii]] for ii in np.arange(len(spectra))]
        param2ranges = [[0,np.max(setupinfo['nspec'][setupinfo_ents])]]*len(spectra)

    col_matrix_binranges = [param1ranges,param2ranges]

    plotz        = 0.0
    voffset      = 0.0
    skyspectra   = [None]*len(plotspecs)
    wavecols_sky = [None]*len(plotspecs)
    fluxcols_sky = [None]*len(plotspecs)
    xrangefull   = [600,3800]

    for plotSNval in [True,False]:
        if plotSNval is False:
            yrangefull   = [-30,150]
        else:
            yrangefull   = 'dummy'

        mwp.plot_1DspecOverview(plotspecs, labels, wavecols, fluxcols, fluxerrcols, plotz, voffset=voffset,
                                skyspectra=skyspectra, wavecols_sky=wavecols_sky, fluxcols_sky=fluxcols_sky,
                                outputfigure=outname, speccols=speccolors, show_error=False,
                                yrangefull=yrangefull, xrangefull=xrangefull,
                                col_matrix=col_matrix, col_matrix_title=col_matrix_title, col_matrix_text=col_matrix_text,
                                col_matrix_labels=col_matrix_labels, col_matrix_ranges=col_matrix_ranges,
                                col_matrix_binranges=col_matrix_binranges,
                                col_matrix_p1dat=col_matrix_p1dat,col_matrix_p2dat=col_matrix_p2dat,
                                plotSN=plotSNval,verbose=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_compositespec_wrapper(specdir='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/stacks1D_sample_selection/'):
    """

    wrapper to stack_composite_plotNxNspecs() setting up and plotting the overveiws of the composite spectra

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    uves.plot_compositespec_wrapper(ztype='zcat')

    """
    bintypes = ['EqNumberBins','EqSizeBins']
    ztypes   = ['zcat','zv18']

    for bintype in bintypes:
        if bintype == 'EqNumberBins':
            binid = 10
        elif bintype == 'EqSizeBins':
            binid = 20
        else:
            sys.exit('Uknown bintype = "'+bintype+'"')
        for ztype in ztypes:
            if ztype == 'zcat':
                zid = 1
            elif ztype == 'zv18':
                zid = 2
            else:
                sys.exit('Uknown ztype = "'+ztype+'"')

            baseid = '_0000'+str(binid+zid)

            # --- two parameters binned ---
            rangedic = {}
            rangedic['z']        = [2.8,6.7]
            rangedic['m814w']    = [23.76,28.82]
            rangedic['Llya']     = [-12420,24500]
            rangedic['FWHMlya']  = [0,17]
            rangedic['EW0lya']   = [0,2000]
            rangedic['beta']     = [-15.5,8.2]

            for param1,param2 in combinations( ['z','m814w','Llya','FWHMlya','EW0lya','beta'],2):
                globstr = specdir+'*'+baseid+'*'+param1+'_bin*'+param2+'_bin*.fits'
                spectra = np.sort(glob.glob(globstr))
                outname = specdir+'composite_overview_'+ztype+'_'+bintype+'_'+param1+'_VS_'+param2+'.pdf'
                if len(spectra) > 0:
                    uves.stack_composite_plotNxNspecs(param1,param2,rangedic[param1],rangedic[param2],spectra,outname,
                                                      ztype=ztype,bintype=bintype)

            # --- one parameter binned ---
            if ztype == 'zcat':
                param1  = 'z'
                globstr = specdir+'*'+baseid+'*'+param1+'LT2p9_bin*of3*.fits'
                spectra = np.sort(glob.glob(globstr))
                outname = specdir+'composite_overview_'+ztype+'_'+bintype+'_'+param1+'LT2p9.pdf'
                if len(spectra) > 0:
                    uves.stack_composite_plotNxNspecs(param1,None,[1.45,2.95],None,spectra,outname,
                                                      ztype=ztype,bintype=bintype)

            for Nbin in [4,10]:
                param1  = 'z'
                globstr = specdir+'*'+baseid+'*'+param1+'GT2p9_bin*of'+str(Nbin)+'*.fits'
                spectra = np.sort(glob.glob(globstr))
                outname = specdir+'composite_overview_'+ztype+'_'+bintype+'_'+param1+'GT2p9_'+str(Nbin)+'bins.pdf'
                if len(spectra) > 0:
                    uves.stack_composite_plotNxNspecs(param1,None,[2.8,6.7],None,spectra,outname,
                                                      ztype=ztype,bintype=bintype)

                param1  = 'm814w'
                globstr = specdir+'*'+baseid+'*'+param1+'_bin*of'+str(Nbin)+'*.fits'
                spectra = np.sort(glob.glob(globstr))
                outname = specdir+'composite_overview_'+ztype+'_'+bintype+'_'+param1+'Detections_'+str(Nbin)+'bins.pdf'
                if len(spectra) > 0:
                    uves.stack_composite_plotNxNspecs(param1,None,rangedic[param1],None,spectra,outname,
                                                      ztype=ztype,bintype=bintype)

            for plotparam in ['Llya','FWHMlya','EW0lya','beta']:
                globstr = specdir+'*'+baseid+'*'+plotparam+'_bin*of4*.fits'
                spectra = np.sort(glob.glob(globstr))
                outname = specdir+'composite_overview_'+ztype+'_'+bintype+'_'+plotparam+'_4bins.pdf'
                if len(spectra) > 0:
                    uves.stack_composite_plotNxNspecs(plotparam,None,rangedic[plotparam],None,spectra,outname,
                                                      ztype=ztype,bintype=bintype)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_vector_intervals(vector,Nsamples,equalsizebins=False,verbose=True):
    """
    Function to split vector in intervals for generating sub-samples.

    --- INPUT ---
    vector        The data vector to bin up
    Nsamples      Number of samples to generate, i.e., the number of bins to return
    equalsizebins To return bins of equal size instead of the default with bins
                  contianing equal number of objects set this keyword to True.
    verbose       Toggle verbosity

    --- Example of use ---
    import uvEmissionlineSearch as uves, astropy.io.fits as afits, numpy as np

    # infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    infodat       = afits.open(infofile)[1].data
    infodat       = infodat[np.where((infodat['id']<4.9e8) | (infodat['id']>5.9e8))[0]]

    binranges     = uves.get_vector_intervals(infodat['redshift'],4)

    """
    if verbose: print(' - Input vector has length '+str(len(vector))+
                      ' but will only consider finite and non-0 values in binning, hence...')
    vector      = np.asarray(vector)[np.isfinite(vector) & (vector != 0)]
    vector_s    = np.sort(vector)

    binranges = []
    if not equalsizebins:
        Nobj_perbin = int(np.floor(len(vector)/Nsamples))
        if verbose:
            print(' - Divding vector of length '+str(len(vector))+' with min and max values ['+str(np.min(vector))+','+str(np.max(vector))+'] into subsamples:')

        if verbose: print(' - The provided vector can be split into the following '+str(Nsamples)+' samples:')
        for bb in np.arange(Nsamples):
            if bb < Nsamples-1:
                binlen = len(vector_s[Nobj_perbin*bb:Nobj_perbin*(bb+1)])
                binmin = np.min(vector_s[Nobj_perbin*bb:Nobj_perbin*(bb+1)])
                binmax = np.max(vector_s[Nobj_perbin*bb:Nobj_perbin*(bb+1)])
            else:
                binlen = len(vector_s[Nobj_perbin*bb:])
                binmin = np.min(vector_s[Nobj_perbin*bb:])
                binmax = np.max(vector_s[Nobj_perbin*bb:])

            if verbose: print('   subsample '+str(bb+1)+'   ['+str("%12.4f" % binmin)+' '+str("%12.4f" % binmax)+
                              ']  of length '+str(binlen))
            binranges.append([binmin,binmax])
    else:
        vectorrange = np.max(vector)-np.min(vector)
        dbin = vectorrange/Nsamples

        for bb in np.arange(Nsamples):
            if bb == 0:
                binmin = np.min(vector) - np.abs(np.min(vector)*0.01)
            else:
                binmin = np.min(vector)+dbin*bb
            if bb == Nsamples-1:
                binmax = np.max(vector) +  np.abs(np.max(vector)*0.01)
            else:
                binmax = np.min(vector) +dbin*(bb+1)
            binlen = len(vector[(vector>=binmin) & (vector<=binmax)])

            if verbose: print('   subsample '+str(bb+1)+'   ['+str("%12.4f" % binmin)+' '+str("%12.4f" % binmax)+
                              ']  containing '+str(binlen)+' objects')
            binranges.append([binmin,binmax])

    return binranges
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_object_duplication_list(matchtol=0.1,verbose=True):
    """
    Function returning a dictionary of objects and their duplicates so they can be ignored/accounted for
    in the analysis.

    --- Example of use ---
    import uvEmissionlineSearch as uves
    ids_ignore, duplicate_dictionary = uves.get_object_duplication_list(matchtol=0.2)

    """
    # infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    infodat       = afits.open(infofile)[1].data
    infodat       = infodat[np.where((infodat['id']<4.9e8) | (infodat['id']>5.9e8))[0]] # ignoring UDF MW mock ids
    Nobj          = len(infodat['id'])
    matchtol_deg  = matchtol/3600.
    duplicate_dic = {}

    dupGT2_ids    = []
    dupGT2_ras    = []
    dupGT2_decs   = []

    if verbose: print(' - Looking for duplicates among the '+str(Nobj)+' objects in infofile ')
    for ii, objid in enumerate(infodat['id']):
        if verbose:
            infostr = '   checking duplicates for id='+str(objid)+\
                      ' ('+str("%.5d" % (ii+1))+' / '+str("%.5d" % len(infodat['id']))+')     '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        objra     = infodat['ra'][ii]
        objdec    = infodat['dec'][ii]
        rmatch    = np.sqrt( (np.cos(np.deg2rad(objdec))*(infodat['ra']-objra))**2.0 + (infodat['dec']-objdec)**2.0 )
        ent_dup   = np.where(rmatch < matchtol_deg)[0]

        if len(ent_dup) > 1:
            ids_dup    = infodat['id'][ent_dup]
            rmatch_dup = rmatch[ent_dup]
            if len(ent_dup) > 2:
                print('WARNING: Found more than 2 matched at location of '+str(objid)+' (within '+str(matchtol)+
                      ' arcsec) with distances and ids \n         '+str(rmatch_dup)+
                      '\n         '+str(ids_dup))
                dupGT2_ids.append(objid)
                dupGT2_ras.append(objra)
                dupGT2_decs.append(objdec)
            elif str(ids_dup[0])[0] == str(ids_dup[1])[0]:
                print('WARNING: The 2 matches to the location of '+str(objid)+' (within '+str(matchtol)+
                      ' arcsec) come from same survey: with distances and ids \n         '+str(rmatch_dup)+
                      '\n         '+str(ids_dup))
                dupGT2_ids.append(objid)
                dupGT2_ras.append(objra)
                dupGT2_decs.append(objdec)
            else:
                idkey = np.min(ids_dup)
                if idkey not in duplicate_dic.keys():
                    duplicate_dic[idkey] = ids_dup[ids_dup != idkey][0]


    if len(dupGT2_ids) > 0: # If there are objects with more than 2 duplications, make region file of those
        regionname = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objects_with_more_than_2_matches_within_'+\
                     str(matchtol).replace('.','p')+'_of_coordinate.reg'
        if verbose: print('\n - Generated DS9 region file with '+str(len(dupGT2_ids))+
                          ' case of >2 matches to coordinate. Saved to\n   '+regionname)
        kbs.create_DS9region(regionname,dupGT2_ras,dupGT2_decs,color='red',circlesize=matchtol,
                             textlist=np.asarray(dupGT2_ids).astype(str),clobber=True,point=None)

    ignoreids = np.sort(np.asarray([duplicate_dic[idkey] for idkey in duplicate_dic.keys()]))

    return ignoreids, duplicate_dic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def object_region_files(basename='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/MWuves-full-v1p0_DS9regions_200213selection',verbose=True):
    """
    Function to generate DS9 region files for the samples in the UVES study

    --- Example of use ---
    import uvEmissionlineSearch as uves
    uves.object_region_files()

    """
    matchtol      = 0.2
    # infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile      = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    infodat       = afits.open(infofile)[1].data
    infodat       = infodat[np.where((infodat['id']<4.9e8) | (infodat['id']>5.9e8))[0]] # ignoring UDF MW mock ids

    if verbose: print(' - Generating DS9 region file for excluded sources ')
    regionname = basename+'_sampleexclusion.reg'
    id_remove  = [158002004,208014258,213022109,600341002,601931670] # the five objects removed from sample analysis
    sampleent  = np.asarray([np.where(infodat['id'] == ir)[0][0] for ir in id_remove])
    ras        = infodat['ra'][sampleent]
    decs       = infodat['dec'][sampleent]
    kbs.create_DS9region(regionname,ras,decs,color='red',circlesize=20,textlist=None,clobber=True,point='x')

    if verbose: print(' - Generating DS9 region file for cosmos sources ')
    regionname = basename+'_cosmos.reg'
    sampleent  = np.where((infodat['id']>1.9e8) & (infodat['id']<2.9e8))[0]
    ras        = infodat['ra'][sampleent]
    decs       = infodat['dec'][sampleent]
    # kbs.create_DS9region(regionname,ras,decs,color='cyan',circlesize=12,textlist=None,clobber=True,point='circle')
    kbs.create_DS9region(regionname,ras,decs,color='cyan',circlesize=matchtol,textlist=None,clobber=True,point=None)

    if verbose: print(' - Generating DS9 region file for cdfs parallel sources ')
    regionname = basename+'_cdfsparallel.reg'
    sampleent  = np.where((infodat['id']>2.9e8) & (infodat['id']<4.9e8))[0]
    ras        = infodat['ra'][sampleent]
    decs       = infodat['dec'][sampleent]
    # kbs.create_DS9region(regionname,ras,decs,color='cyan',circlesize=12,textlist=None,clobber=True,point='circle')
    kbs.create_DS9region(regionname,ras,decs,color='cyan',circlesize=matchtol,textlist=None,clobber=True,point=None)

    if verbose: print(' - Generating DS9 region file for cdfs sources ')
    regionname = basename+'_cdfs.reg'
    sampleent  = np.where((infodat['id']>0.9e8) & (infodat['id']<1.9e8))[0]
    ras        = infodat['ra'][sampleent]
    decs       = infodat['dec'][sampleent]
    # kbs.create_DS9region(regionname,ras,decs,color='cyan',circlesize=12,textlist=None,clobber=True,point='circle')
    kbs.create_DS9region(regionname,ras,decs,color='cyan',circlesize=matchtol,textlist=None,clobber=True,point=None)

    if verbose: print(' - Generating DS9 region file for UDF sources ')
    regionname = basename+'_udfmosaic.reg'
    sampleent  = np.where((infodat['id']>5.9e8) & (infodat['id']<6.9e8))[0]
    ras        = infodat['ra'][sampleent]
    decs       = infodat['dec'][sampleent]
    kbs.create_DS9region(regionname,ras,decs,color='green',circlesize=20,textlist=None,clobber=True,point='diamond')

    if verbose: print(' - Generating DS9 region file for UDF10 sources ')
    regionname = basename+'_udf10.reg'
    sampleent  = np.where((infodat['id']>6.9e8) & (infodat['id']<7.9e8))[0]
    ras        = infodat['ra'][sampleent]
    decs       = infodat['dec'][sampleent]
    kbs.create_DS9region(regionname,ras,decs,color='magenta',circlesize=20,textlist=None,clobber=True,point='box')

    if verbose: print(' - Generating DS9 region file for duplicates to ignore ')
    regionname = basename+'_dupicates2ignore.reg'
    ids_ignore, duplicate_dictionary = uves.get_object_duplication_list(matchtol=matchtol)
    sampleent  = np.asarray([np.where(infodat['id'] == ig)[0][0] for ig in ids_ignore])
    ras        = infodat['ra'][sampleent]
    decs       = infodat['dec'][sampleent]
    kbs.create_DS9region(regionname,ras,decs,color='red',circlesize=20,textlist=None,clobber=True,point='cross')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_uves_FoV(figbasename,mastercat,infofile,figext='.pdf',showobjects=True,objectsAsDots=False,
                  usehighresimages=True,pointingswithnumbers=True,convolveimg=False,verbose=True):
    """
    Generata plot of CDFS and COSMOS region with UVES samples overplotted

    --- Example of use ---
    import uvEmissionlineSearch as uves
    mastercat   = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/back2backAnalysis_200213/results_master_catalog_version200213.fits'
    # infofile    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'

    # with object details
    figbasename = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVfigures/MUSE-Wide_FoV'
    uves.plot_uves_FoV(figbasename,mastercat,infofile,showobjects=True,objectsAsDots=False,pointingswithnumbers=True,verbose=True)

    # more publication friendly:
    figbasename = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVfigures/MUSE-Wide_FoV_pubfriendly'
    uves.plot_uves_FoV(figbasename,mastercat,infofile,showobjects=True,objectsAsDots=True,pointingswithnumbers=False,verbose=True)

    # Figure with pointing qualtiy summary shown:
    figbasename = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVfigures/MUSE-Wide_FoV_pubfriendly_quality'
    uves.plot_uves_FoV(figbasename,mastercat,infofile,showobjects=False,objectsAsDots=False,pointingswithnumbers=True,verbose=True)

    """
    agn, agncand = uves.get_AGN_ids()

    if verbose: print(' - Loading master cat and defining sub-samples')
    masterdat = afits.open(mastercat)[1].data
    infodat   = afits.open(infofile)[1].data

    if verbose: print('   o Objects with at least one detection')
    sel_UVdet = np.where(( ((np.abs(masterdat['ferr_CIV'])   != 99.0) & np.isfinite(masterdat['ferr_CIV']))   |
                           ((np.abs(masterdat['ferr_NV'])    != 99.0) & np.isfinite(masterdat['ferr_NV']))    |
                           ((np.abs(masterdat['ferr_HeII'])  != 99.0) & np.isfinite(masterdat['ferr_HeII']))  |
                           ((np.abs(masterdat['ferr_OIII'])  != 99.0) & np.isfinite(masterdat['ferr_OIII']))  |
                           ((np.abs(masterdat['ferr_SiIII']) != 99.0) & np.isfinite(masterdat['ferr_SiIII'])) |
                           # ((np.abs(masterdat['ferr_MgII'])  != 99.0) & np.isfinite(masterdat['ferr_MgII']))  |
                           ((np.abs(masterdat['ferr_CIII'])  != 99.0) & np.isfinite(masterdat['ferr_CIII']))   ) &
                            (masterdat['redshift'] >= 0.0) & (masterdat['duplicationID'] == 0.0) )[0]

    if verbose: print('   o LAEs ')
    # sel_LAEs = np.where( (masterdat['redshift'] >= 2.9) & (masterdat['duplicationID'] == 0.0) )[0]
    sel_LAEs = np.where( (masterdat['redshift'] >= 2.9) )[0]

    if verbose: print('   o NV emitters ')
    sel_NV   = np.where(((np.abs(masterdat['ferr_NV']) != 99.0) & np.isfinite(masterdat['ferr_NV'])) &
                        (masterdat['duplicationID'] == 0.0) )[0]

    if verbose: print('   o CIV emitters ')
    sel_CIV  = np.where(((np.abs(masterdat['ferr_CIV']) != 99.0) & np.isfinite(masterdat['ferr_CIV'])) &
                        (masterdat['duplicationID'] == 0.0) )[0]

    if verbose: print('   o HeII emitters ')
    sel_HeII = np.where(((np.abs(masterdat['ferr_HeII']) != 99.0) & np.isfinite(masterdat['ferr_HeII'])) &
                        (masterdat['duplicationID'] == 0.0) )[0]

    if verbose: print('   o OIII emitters ')
    sel_OIII = np.where(((np.abs(masterdat['ferr_OIII']) != 99.0) & np.isfinite(masterdat['ferr_OIII'])) &
                        (masterdat['duplicationID'] == 0.0) )[0]

    if verbose: print('   o SiIII emitters ')
    sel_SiIII= np.where(((np.abs(masterdat['ferr_SiIII']) != 99.0) & np.isfinite(masterdat['ferr_SiIII'])) &
                        (masterdat['duplicationID'] == 0.0) )[0]

    if verbose: print('   o CIII emitters ')
    sel_CIII = np.where(((np.abs(masterdat['ferr_CIII']) != 99.0) & np.isfinite(masterdat['ferr_CIII'])) &
                        (masterdat['duplicationID'] == 0.0) )[0]

    if verbose: print('   o MgII emitters ')
    sel_MgII = np.where(((np.abs(masterdat['ferr_MgII']) != 99.0) & np.isfinite(masterdat['ferr_MgII'])) &
                        (masterdat['duplicationID'] == 0.0) )[0]

    if verbose: print(' - Defining images and ranges')
    fields  = ['GOODS-S','COSMOS','UDF','MXDFregion']

    if usehighresimages:
        imagefiles = ['/Users/kschmidt/work/images_MAST/goodss_3dhst.v4.0.F125W_F140W_F160W_det.fits',
                      '/Users/kschmidt/work/images_MAST/cosmos_3dhst.v4.0.F125W_F140W_F160W_det.fits',
                      '/Users/kschmidt/work/images_MAST/hlsp_hlf_hst_acs-60mas_goodss_f775w_v2.0_sci.fits',
                      '/Users/kschmidt/work/images_MAST/hlsp_hlf_hst_acs-60mas_goodss_f775w_v2.0_sci.fits']
    else:
        imagefiles = ['/Users/kschmidt/work/images_MAST/goodss_3dhst.v4.0.F125W_F140W_F160W_det_lowres.fits',
                      '/Users/kschmidt/work/images_MAST/cosmos_3dhst.v4.0.F125W_F140W_F160W_det_lowres.fits',
                      '/Users/kschmidt/work/images_MAST/hlsp_hlf_hst_acs-60mas_goodss_f775w_v2.0_sci_lowres.fits',
                      '/Users/kschmidt/work/images_MAST/hlsp_hlf_hst_acs-60mas_goodss_f775w_v2.0_sci_lowres.fits']


    # imagefiles = ['/Users/kschmidt/work/images_MAST/MUSEWidePointings/wfc3_160w_candels-cdfs-15_cut_v1.0.fits',
    #               '/Users/kschmidt/work/images_MAST/cosmos_mosaic_Shrink50.fits']
    # vmin=0.0001
    # vmax=0.5

    if pointingswithnumbers:
        pointingregions = ['/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVfigures/candels_cdfs_pointings.reg',
                           '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVfigures/candels_cosmos_pointings.reg',
                           '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVfigures/candels_udf_pointings.reg',
                           '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVfigures/candels_udf_pointings.reg']
    else:
        pointingregions = ['/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVfigures/candels_cdfs_pointings_nonumbers.reg',
                           '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVfigures/candels_cosmos_pointings_nonumbers.reg',
                           '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVfigures/candels_udf_pointings_nonumbers.reg',
                           '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVfigures/candels_udf_pointings_nonumbers.reg']

    PSFcat = '/Users/kschmidt/work/MUSE/MUSEWide_PSFs/MUSEWide_PSF_Catalog_180723.fits'
    PSFdat = afits.open(PSFcat)[1].data
    QRcat  = '/Users/kschmidt/work/catalogs/MUSE_GTO/MUSE-Wide_pointings-qualtiy_info.txt'
    QRdat  = np.genfromtxt(QRcat,names=True,skip_header=2,dtype=('40a','40a','40a','40a','d','d','d','d'),comments='#')

    titletexts = ['HUDF Parallels and GOODS-S','COSMOS','The UDF mosaic and UDF10','MXDF']

    for ii, imagefile in enumerate(imagefiles):
        # if ii != 3: # only plotting MXDF
        #     continue
        plotname  = figbasename+'_'+fields[ii]+figext
        if showobjects:
            plotname  = plotname.replace(figext,'_withobj'+figext)
        if verbose: print(' - Generating '+plotname)
        pointings = pointingregions[ii]
        hdu       = afits.open(imagefile)[0]
        hud_wcs   = wcs.WCS(hdu.header)

        if 'goodss_3dhst.v4.0.F125W_F140W_F160W_det' in imagefile:
            fig = plt.figure(figsize=(5, 5))
            vmin = 0.1
            vmax = 10.
        elif 'cosmos_3dhst.v4.0.F125W_F140W_F160W_det' in imagefile:
            fig = plt.figure(figsize=(5, 5))
            vmin = 0.5
            vmax = 10.
        elif 'hlsp_xdf_hst_acswfc-60mas_hudf_f814w_v1_sci' in imagefile:
            fig = plt.figure(figsize=(5, 5))
            vmin = 0.001
            vmax = 1.
        elif 'hlsp_hlf_hst_acs-60mas_goodss_f775w_v2.0_sci' in imagefile:
            fig = plt.figure(figsize=(5, 5))
            vmin = 0.0001
            vmax = 5.0
        else:
            fig = plt.figure(figsize=(5, 5))
            vmin = 0.0001
            vmax = 0.5

        if ('COSMOS' in fields[ii]):
            fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.25, right=0.97, bottom=0.15, top=0.90)
        else:
            fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.25, right=0.97, bottom=0.15, top=0.97)

        Fsize    = 14
        lthick   = 2
        marksize = 4
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        # plt.title(titletexts[ii],fontsize=Fsize)

        ax = plt.subplot(projection=hud_wcs, label='overlays')
        ax.set_title(titletexts[ii],fontsize=Fsize)

        if convolveimg:
            if verbose: print('   starting image convolution...')
            gauss_kernel = Gaussian2DKernel(0.25)
            img2show     = convolve_fft(hdu.data,gauss_kernel)
            if verbose: print('   done; moving on.')
        else:
            img2show = hdu.data
        #ax.imshow(hdu.data, origin='lower', cmap='Greys', vmin=-0.001, vmax=0.05) #, vmin=-2.e-5, vmax=2.e-4,interpolation=None
        ax.imshow(img2show, origin='lower', cmap='Greys', vmin=vmin, vmax=vmax, norm=LogNorm()) #, vmin=-2.e-5, vmax=2.e-4,interpolation=None

        ax.coords.grid(True, color='black', ls='dotted')
        ax.coords[0].set_axislabel('Right Ascension (J2000)')
        ax.coords[1].set_axislabel('Declination (J2000)')

        if 'MXDF' in fields[ii]:
            regstr = 'circle(53.16467,-27.78537,44.0") # color=magenta width=2  font="times 10 bold roman" text={MXDF; r=44"} '
            r2 = pyregion.parse(regstr).as_imagecoord(header=hdu.header)
            patch_list, artist_list = r2.get_mpl_patches_texts()
            for pp in patch_list:
                ax.add_patch(pp)
            for tt in artist_list:
                ax.add_artist(tt)

        # overlay = ax.get_coords_overlay('fk5')
        # overlay.grid(color='white', ls='dotted')
        # overlay[0].set_axislabel('Right Ascension (J2000)')
        # overlay[1].set_axislabel('Declination (J2000)')

        ds9regs = pyregion.open(pointings).as_imagecoord(header=hdu.header)
        patch_list, artist_list = ds9regs.get_mpl_patches_texts()

        for pp in patch_list:
            ax.add_patch(pp)
        for tt in artist_list:
            field = pointings.split('_')[1]+'-'+tt.get_text()
            fieldent = np.where(PSFdat['field'] == field)[0]
            if len(fieldent) == 1:
                PSFfwhm = PSFdat['p0_sel_g'][fieldent]  # FWHM at 7050A         + PSFdat['p1_sel_g'] * (lambda - 7050)
                # pointing_quality_text = '\\\\FWHM(PSF) = '+str("%.2f" % PSFfwhm)+"''"
                pointing_quality_text = ': '+str("%.2f" % PSFfwhm)
                tt.set_text(tt.get_text()+pointing_quality_text)
                tt.set_fontsize(4)
                tt.set_color('blue')
            ax.add_artist(tt)

        # regstr = 'box(53.1243740333,-27.8516127209,0.01666667,0.01666667,340.0) # color=blue width=3 font="times 10 bold roman" text={Testing box}'
        # text 53.11950302d -27.85599899d {15} # color=blue
        # regstr = 'circle(53.1243740333,-27.8516127209,0.5") # color=magenta width=3  font="times 10 bold roman" text={115003085} '
        # r2 = pyregion.parse(regstr).as_imagecoord(header=hdu.header)
        # patch_list, artist_list = r2.get_mpl_patches_texts()

        if showobjects:
            if verbose: print('   creating object patches...')
            for oo, objid in enumerate(masterdat['id']):
                if 'goodss_3dhst.v4.0.F125W_F140W_F160W_det' in imagefile:
                    if (str(objid)[0] == '2') or (int(str(objid)[0]) >= 6):
                        continue
                elif 'cosmos_3dhst.v4.0.F125W_F140W_F160W_det' in imagefile:
                    if str(objid)[0] != '2':
                        continue
                elif 'hlsp_hlf_hst_acs-60mas_goodss_f775w_v2.0_sci' in imagefile:
                    if (int(str(objid)[0]) < 6):
                        continue
                elif 'hlsp_xdf_hst_acswfc-60mas_hudf_f814w_v1_sci' in imagefile:
                    if (int(str(objid)[0]) < 6):
                        continue

                infoent = np.where(infodat['id'] == objid)[0]

                if objectsAsDots:
                    plotmarker = 'o'
                    facecolor = 'forestgreen'
                    edgecolor = 'forestgreen'
                    symsize    = 7
                else:
                    if objid in masterdat['id'][sel_LAEs]:
                        plotmarker = 's'
                    else:
                        plotmarker = 'o'

                    if (objid in masterdat['id'][sel_CIII]) & (objid not in masterdat['id'][sel_CIV]):
                        pointcolor = 'forestgreen'
                    elif (objid in masterdat['id'][sel_CIV]) & (objid not in masterdat['id'][sel_CIII]):
                        pointcolor = 'blue'
                    elif (objid in masterdat['id'][sel_CIV]) & (objid in masterdat['id'][sel_CIII]):
                        pointcolor = 'darkorange'
                    else:
                        pointcolor = 'red'

                    if objid in masterdat['id'][sel_UVdet]:
                        facecolor  = pointcolor
                        edgecolor  = 'None'
                        symsize    = 10
                    else:
                        facecolor  = 'None'
                        edgecolor  = pointcolor
                        symsize    = 7

                    if (objid in masterdat['id'][sel_HeII]):
                        plotmarker = 'P'
                        symsize = 20

                    if (objid in masterdat['id'][sel_NV]):
                        plotmarker = 'D'
                        symsize = 20

                    if (objid in agn) or (objid in agncand):
                        plotmarker = 'X'
                        symsize = 20

                if 'MXDF' in fields[ii]:
                    symsize = symsize*2.0

                ax.scatter(infodat['ra'][infoent],infodat['dec'][infoent], transform=ax.get_transform('fk5'), marker=plotmarker,
                           s=symsize, edgecolor=edgecolor, facecolor=facecolor, alpha=0.5, zorder=1e10, lw=0.5)
            if verbose: print('   done; moving on.')

        xmin, xmax, ymin, ymax = ax.axis()
        if 'goodss_3dhst.v4.0.F125W_F140W_F160W_det' in imagefile:
            plt.xlim([xmax*0.05,xmax*0.9])
            plt.ylim([ymax*0.25,ymax*1.0])
        elif 'cosmos_3dhst.v4.0.F125W_F140W_F160W_det' in imagefile:
            plt.xlim([xmax*0.18,xmax*0.70])
            plt.ylim([ymax*0.23,ymax*0.48])
        elif 'hlsp_hlf_hst_acs-60mas_goodss_f775w_v2.0_sci' in imagefile:
            if 'MXDF' in fields[ii]:
                plt.xlim([xmax*0.38,xmax*0.445])
                plt.ylim([ymax*0.51,ymax*0.58])
            else:
                plt.xlim([xmax*0.33,xmax*0.51])
                plt.ylim([ymax*0.45,ymax*0.64])
        elif 'hlsp_xdf_hst_acswfc-60mas_hudf_f814w_v1_sci' in imagefile:
            if 'MXDF' in fields[ii]:
                plt.xlim([xmax*0.30,xmax*0.65])
                plt.ylim([ymax*0.40,ymax*0.75])
            else:
                plt.xlim([xmax*0.10,xmax*0.95])
                plt.ylim([ymax*0.15,ymax*0.95])

        # plt.xlabel('R.A.')
        # plt.ylabel('Dec.')

        # leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize/1.7},ncol=3,numpoints=1,
        #                  bbox_to_anchor=(0.5, 1.1),)  # add the legend
        # leg.get_frame().set_alpha(0.7)

        plt.savefig(plotname)
        plt.clf()
        plt.close('all')
        if verbose: print('   Succesfully stored plot to figure')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def mastercat_latextable_wrappers(mastercatalog,infofile,sortcol='id',overwrite=False,verbose=True):
    """
    Wrapper to uves.mastercat_latextable() generating a table of the emitters with either CIII or CIV detections.
    This wraper also returns a dictionary with number counts of detections and available objects for searching
    based on the content of the master catalog.

    --- Example of use ---
    import uvEmissionlineSearch as uves
    mastercat = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/back2backAnalysis_200213/results_master_catalog_version200213.fits'
    # infofile  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    detectiondic = uves.mastercat_latextable_wrappers(mastercat,infofile,verbose=True,overwrite=True,sortcol='redshift')

    """

    infodat   = afits.open(infofile)[1].data
    masterdat = afits.open(mastercatalog)[1].data

    if verbose: print(' - Defining colum sets for tables to generate')
    columnslist   = [['ra','dec','redshift','EW_0','EW_0_err',
                      'id_skelton','sep_skelton','id_rafelski','sep_rafelski','id_guo','sep_guo','id_laigle','sep_laigle'],
                     ['ra','dec','redshift','f_NV','f_CIV','f_HeII','f_OIII','f_SiIII','f_CIII','f_MgII'],
                     ['ra','dec','redshift',
                      'FR_NV1NV2','FR_CIV1CIV2','FR_OIII1OIII2','FR_SiIII1SiIII2','FR_CIII1CIII2','FR_MgII1MgII2']]
    tabstrings    = ['lyainfo','fluxes','fluxratios']

    if verbose: print(' - Collecting objects for the individual tables ')
    outname_bases = []
    goodents      = []

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # stat dictionary
    detectiondic   = collections.OrderedDict()
    areanames      = ['100fields','musewide','cdfsproper','cosmos','cdfspar1','cdfspar2','udfmosaic','udf10']
    idstart        = [1,2,3,4,6,7]
    idstartname    = areanames[2:]
    rangename      = ['allobj','allobjUVrange','CIIIorCIV','Lya', 'NV','CIV','HeII','OIII','SiIII','CIII','MgII','Non-LAE-NoCIV']
    zranges        = [[0,10.0],[0,4.9699],[1.5241,4.9699],[2.9729,6.5958],[2.8918,6.4432],[2.1114,4.9699],[1.9379,4.6411],
                      [1.8969,4.5632],[1.5514,3.9067],[1.5241,3.8548],[0.7174,2.3142],[1.5,2.1114]]
    ignoreIDlist   = [158002004,601931670,208014258,600341002]
    subsel         = ['full','LAEs']
    Nkeys_init     = 0
    for rr in rangename:
        for aa in areanames:
            for ss in subsel:
                detectiondic[aa+'_'+rr+'_'+ss] = [np.nan]
                Nkeys_init = Nkeys_init + 1

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('   o Objects with at least one detection')
    selection = np.where(( ((np.abs(masterdat['ferr_NV'])    != 99.0) & np.isfinite(masterdat['ferr_NV']))    |
                           ((np.abs(masterdat['ferr_CIV'])   != 99.0) & np.isfinite(masterdat['ferr_CIV']))   |
                           ((np.abs(masterdat['ferr_HeII'])  != 99.0) & np.isfinite(masterdat['ferr_HeII']))  |
                           ((np.abs(masterdat['ferr_OIII'])  != 99.0) & np.isfinite(masterdat['ferr_OIII']))  |
                           ((np.abs(masterdat['ferr_SiIII']) != 99.0) & np.isfinite(masterdat['ferr_SiIII'])) |
                           ((np.abs(masterdat['ferr_CIII'])  != 99.0) & np.isfinite(masterdat['ferr_CIII']))  |
                           ((np.abs(masterdat['ferr_MgII'])  != 99.0) & np.isfinite(masterdat['ferr_MgII']))   ) &
                            (masterdat['redshift'] >= 0.0) & (masterdat['redshift'] <= 6.4432) & (masterdat['duplicationID'] == 0.0) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) )[0]
    selectionTOT = np.where((masterdat['duplicationID'] == 0.0) &  (masterdat['redshift'] <= 6.4432) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) )[0]

    outname_bases.append('/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/mastercat_AnyUVLineDetection.tex')
    goodents.append(selection)
    if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
    Nwide   = str(len(masterdat['id'][selection][masterdat['id'][selection] < 5e8]))
    Nmosaic = str(len(masterdat['id'][selection][(masterdat['id'][selection] > 6e8) & (masterdat['id'][selection] < 7e8)]))
    Nudf10  = str(len(masterdat['id'][selection][masterdat['id'][selection] > 7e8]))

    Nwidetot   = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] < 5e8]))
    Nmosaictot = str(len(masterdat['id'][selectionTOT][(masterdat['id'][selectionTOT] > 6e8) & (masterdat['id'][selectionTOT] < 7e8)]))
    Nudf10tot  = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] > 7e8]))
    if verbose: print('     Of these '+Nwide+'/'+Nmosaic+'/'+Nudf10+
                      ' ('+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+
                      '%/'+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+
                      '%/'+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'%) '+
                      'are in the Wide/Mosaic/UDF10 fields ')

    if verbose:
        tabstring = ' At least one detection (all) & XXzrangeXX & '\
        +str(len(selection))+' & '+str(len(selectionTOT))+' & '+str("%5.2f" % (float(len(selection))/float(len(selectionTOT))*100.0))+'\% & '\
        +str(Nwide)+         ' & '+str(Nwidetot)+         ' & '+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+'\% & '\
        +str(Nmosaic)+       ' & '+str(Nmosaictot)+       ' & '+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+'\% & '\
        +str(Nudf10)+        ' & '+str(Nudf10tot)+        ' & '+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'\% \\\\ '
        print(tabstring)

    detectiondic['100fields_allobj_full'] = [len(selection)]
    detectiondic['musewide_allobj_full']  = [int(Nwide)]
    detectiondic['udfmosaic_allobj_full'] = [int(Nmosaic)]
    detectiondic['udf10_allobj_full']     = [int(Nudf10)]
    detectiondic['100fields_allobjUVrange_full'] = [len(selection)]
    detectiondic['musewide_allobjUVrange_full']  = [int(Nwide)]
    detectiondic['udfmosaic_allobjUVrange_full'] = [int(Nmosaic)]
    detectiondic['udf10_allobjUVrange_full']     = [int(Nudf10)]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('   o Desert obj with at least one other detection')
    selection = np.where(( ((np.abs(masterdat['ferr_NV'])    != 99.0) & np.isfinite(masterdat['ferr_NV']))    |
                           ((np.abs(masterdat['ferr_CIV'])   != 99.0) & np.isfinite(masterdat['ferr_CIV']))   |
                           ((np.abs(masterdat['ferr_HeII'])  != 99.0) & np.isfinite(masterdat['ferr_HeII']))  |
                           ((np.abs(masterdat['ferr_OIII'])  != 99.0) & np.isfinite(masterdat['ferr_OIII']))  |
                           ((np.abs(masterdat['ferr_SiIII']) != 99.0) & np.isfinite(masterdat['ferr_SiIII'])) |
                           ((np.abs(masterdat['ferr_CIII'])  != 99.0) & np.isfinite(masterdat['ferr_CIII']))  |
                           ((np.abs(masterdat['ferr_MgII'])  != 99.0) & np.isfinite(masterdat['ferr_MgII']))   ) &
                            (masterdat['redshift'] >= 1.5) & (masterdat['redshift'] < 2.9) & (masterdat['duplicationID'] == 0.0) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) )[0]

    selectionTOT = np.where((masterdat['redshift'] >= 1.5) & (masterdat['redshift'] < 2.9) &  (masterdat['duplicationID'] == 0.0) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) )[0]
    outname_bases.append('/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/mastercat_DesertWithUVLineDetection.tex')
    goodents.append(selection)
    if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
    if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
    Nwide   = str(len(masterdat['id'][selection][masterdat['id'][selection] < 5e8]))
    Nmosaic = str(len(masterdat['id'][selection][(masterdat['id'][selection] > 6e8) & (masterdat['id'][selection] < 7e8)]))
    Nudf10  = str(len(masterdat['id'][selection][masterdat['id'][selection] > 7e8]))

    Nwidetot   = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] < 5e8]))
    Nmosaictot = str(len(masterdat['id'][selectionTOT][(masterdat['id'][selectionTOT] > 6e8) & (masterdat['id'][selectionTOT] < 7e8)]))
    Nudf10tot  = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] > 7e8]))
    if verbose: print('     Of these '+Nwide+'/'+Nmosaic+'/'+Nudf10+
                      ' ('+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+
                      '%/'+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+
                      '%/'+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'%) '+
                      'are in the Wide/Mosaic/UDF10 fields ')

    if verbose:
        tabstring = ' At least one detection (Desert) & XXzrangeXX & '\
        +str(len(selection))+' & '+str(len(selectionTOT))+' & '+str("%5.2f" % (float(len(selection))/float(len(selectionTOT))*100.0))+'\% & '\
        +str(Nwide)+         ' & '+str(Nwidetot)+         ' & '+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+'\% & '\
        +str(Nmosaic)+       ' & '+str(Nmosaictot)+       ' & '+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+'\% & '\
        +str(Nudf10)+        ' & '+str(Nudf10tot)+        ' & '+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'\% \\\\ '
        print(tabstring)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('   o LAEs with at least one other detection')
    selection = np.where(( ((np.abs(masterdat['ferr_NV'])    != 99.0) & np.isfinite(masterdat['ferr_NV']))    |
                           ((np.abs(masterdat['ferr_CIV'])   != 99.0) & np.isfinite(masterdat['ferr_CIV']))   |
                           ((np.abs(masterdat['ferr_HeII'])  != 99.0) & np.isfinite(masterdat['ferr_HeII']))  |
                           ((np.abs(masterdat['ferr_OIII'])  != 99.0) & np.isfinite(masterdat['ferr_OIII']))  |
                           ((np.abs(masterdat['ferr_SiIII']) != 99.0) & np.isfinite(masterdat['ferr_SiIII'])) |
                           ((np.abs(masterdat['ferr_CIII'])  != 99.0) & np.isfinite(masterdat['ferr_CIII']))  |
                           ((np.abs(masterdat['ferr_MgII'])  != 99.0) & np.isfinite(masterdat['ferr_MgII']))   ) &
                            (masterdat['redshift'] >= 2.9) & (masterdat['redshift'] <= 6.4432) & (masterdat['duplicationID'] == 0.0) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) )[0]

    selectionTOT = np.where((masterdat['redshift'] >= 2.9) & (masterdat['redshift'] <= 6.4432) &  (masterdat['duplicationID'] == 0.0) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) )[0]
    outname_bases.append('/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/mastercat_LAEsWithUVLineDetection.tex')
    goodents.append(selection)
    if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
    if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
    Nwide   = str(len(masterdat['id'][selection][masterdat['id'][selection] < 5e8]))
    Nmosaic = str(len(masterdat['id'][selection][(masterdat['id'][selection] > 6e8) & (masterdat['id'][selection] < 7e8)]))
    Nudf10  = str(len(masterdat['id'][selection][masterdat['id'][selection] > 7e8]))

    Nwidetot   = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] < 5e8]))
    Nmosaictot = str(len(masterdat['id'][selectionTOT][(masterdat['id'][selectionTOT] > 6e8) & (masterdat['id'][selectionTOT] < 7e8)]))
    Nudf10tot  = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] > 7e8]))
    if verbose: print('     Of these '+Nwide+'/'+Nmosaic+'/'+Nudf10+
                      ' ('+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+
                      '%/'+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+
                      '%/'+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'%) '+
                      'are in the Wide/Mosaic/UDF10 fields ')

    if verbose:
        tabstring = ' At least one detection (LAEs) & XXzrangeXX & '\
        +str(len(selection))+' & '+str(len(selectionTOT))+' & '+str("%5.2f" % (float(len(selection))/float(len(selectionTOT))*100.0))+'\% & '\
        +str(Nwide)+         ' & '+str(Nwidetot)+         ' & '+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+'\% & '\
        +str(Nmosaic)+       ' & '+str(Nmosaictot)+       ' & '+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+'\% & '\
        +str(Nudf10)+        ' & '+str(Nudf10tot)+        ' & '+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'\% \\\\ '
        print(tabstring)

    detectiondic['100fields_allobj_LAEs'] = [len(selection)]
    detectiondic['musewide_allobj_LAEs']  = [int(Nwide)]
    detectiondic['udfmosaic_allobj_LAEs'] = [int(Nmosaic)]
    detectiondic['udf10_allobj_LAEs']     = [int(Nudf10)]
    detectiondic['100fields_allobjUVrange_LAEs'] = [len(selection)]
    detectiondic['musewide_allobjUVrange_LAEs']  = [int(Nwide)]
    detectiondic['udfmosaic_allobjUVrange_LAEs'] = [int(Nmosaic)]
    detectiondic['udf10_allobjUVrange_LAEs']     = [int(Nudf10)]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('   o Carbon emitters (either CIII or CIV detected) ')
    selection = np.where(( ((np.abs(masterdat['ferr_CIV'])  != 99.0) & np.isfinite(masterdat['ferr_CIV'])) |
                           ((np.abs(masterdat['ferr_CIII']) != 99.0) & np.isfinite(masterdat['ferr_CIII'])) ) &
                            (masterdat['duplicationID'] == 0.0) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) )[0]

    selectionTOT = np.where((masterdat['redshift'] >= 1.5241) & (masterdat['redshift'] <= 4.9699) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) &
                            (masterdat['duplicationID'] == 0.0) )[0]
    outname_bases.append('/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/mastercat_CarbonEmitters.tex')
    goodents.append(selection)
    if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
    if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
    Nwide   = str(len(masterdat['id'][selection][masterdat['id'][selection] < 5e8]))
    Nmosaic = str(len(masterdat['id'][selection][(masterdat['id'][selection] > 6e8) & (masterdat['id'][selection] < 7e8)]))
    Nudf10  = str(len(masterdat['id'][selection][masterdat['id'][selection] > 7e8]))

    Nwidetot   = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] < 5e8]))
    Nmosaictot = str(len(masterdat['id'][selectionTOT][(masterdat['id'][selectionTOT] > 6e8) & (masterdat['id'][selectionTOT] < 7e8)]))
    Nudf10tot  = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] > 7e8]))
    if verbose: print('     Of these '+Nwide+'/'+Nmosaic+'/'+Nudf10+
                      ' ('+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+
                      '%/'+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+
                      '%/'+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'%) '+
                      'are in the Wide/Mosaic/UDF10 fields ')

    if verbose:
        tabstring = ' Cemitter (all) & XXzrangeXX & '\
        +str(len(selection))+' & '+str(len(selectionTOT))+' & '+str("%5.2f" % (float(len(selection))/float(len(selectionTOT))*100.0))+'\% & '\
        +str(Nwide)+         ' & '+str(Nwidetot)+         ' & '+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+'\% & '\
        +str(Nmosaic)+       ' & '+str(Nmosaictot)+       ' & '+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+'\% & '\
        +str(Nudf10)+        ' & '+str(Nudf10tot)+        ' & '+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'\% \\\\ '
        print(tabstring)

    detectiondic['100fields_CIIIorCIV_full'] = [len(selection)]
    detectiondic['musewide_CIIIorCIV_full']  = [int(Nwide)]
    detectiondic['udfmosaic_CIIIorCIV_full'] = [int(Nmosaic)]
    detectiondic['udf10_CIIIorCIV_full']     = [int(Nudf10)]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('   o Deserts with Carbon detection (either CIII or CIV detected)')
    selection = np.where(( ((np.abs(masterdat['ferr_CIV'])   != 99.0) & np.isfinite(masterdat['ferr_CIV']))   |
                           ((np.abs(masterdat['ferr_CIII'])  != 99.0) & np.isfinite(masterdat['ferr_CIII']))   ) &
                            (masterdat['redshift'] >= 1.5) & (masterdat['redshift'] < 2.9) & (masterdat['duplicationID'] == 0.0) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) )[0]

    selectionTOT = np.where((masterdat['redshift'] >= 1.5) & (masterdat['redshift'] < 2.9) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) &
                            (masterdat['duplicationID'] == 0.0) )[0]

    outname_bases.append('/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/mastercat_DesertWithCarbonDetection.tex')
    goodents.append(selection)
    if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
    if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
    Nwide   = str(len(masterdat['id'][selection][masterdat['id'][selection] < 5e8]))
    Nmosaic = str(len(masterdat['id'][selection][(masterdat['id'][selection] > 6e8) & (masterdat['id'][selection] < 7e8)]))
    Nudf10  = str(len(masterdat['id'][selection][masterdat['id'][selection] > 7e8]))

    Nwidetot   = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] < 5e8]))
    Nmosaictot = str(len(masterdat['id'][selectionTOT][(masterdat['id'][selectionTOT] > 6e8) & (masterdat['id'][selectionTOT] < 7e8)]))
    Nudf10tot  = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] > 7e8]))
    if verbose: print('     Of these '+Nwide+'/'+Nmosaic+'/'+Nudf10+
                      ' ('+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+
                      '%/'+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+
                      '%/'+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'%) '+
                      'are in the Wide/Mosaic/UDF10 fields ')

    if verbose:
        tabstring = ' Cemitter (Desert) & XXzrangeXX & '\
        +str(len(selection))+' & '+str(len(selectionTOT))+' & '+str("%5.2f" % (float(len(selection))/float(len(selectionTOT))*100.0))+'\% & '\
        +str(Nwide)+         ' & '+str(Nwidetot)+         ' & '+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+'\% & '\
        +str(Nmosaic)+       ' & '+str(Nmosaictot)+       ' & '+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+'\% & '\
        +str(Nudf10)+        ' & '+str(Nudf10tot)+        ' & '+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'\% \\\\ '
        print(tabstring)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('   o LAEs with Carbon detection (either CIII or CIV detected)')
    selection = np.where(( ((np.abs(masterdat['ferr_CIV'])   != 99.0) & np.isfinite(masterdat['ferr_CIV']))   |
                           ((np.abs(masterdat['ferr_CIII'])  != 99.0) & np.isfinite(masterdat['ferr_CIII']))   ) &
                            (masterdat['redshift'] >= 2.9) & (masterdat['duplicationID'] == 0.0) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) )[0]

    selectionTOT = np.where((masterdat['redshift'] >= 2.9) & (masterdat['redshift'] <= 4.9699) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) &
                            (masterdat['duplicationID'] == 0.0) )[0]

    outname_bases.append('/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/mastercat_LAEsWithCarbonDetection.tex')
    goodents.append(selection)
    if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
    if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
    Nwide   = str(len(masterdat['id'][selection][masterdat['id'][selection] < 5e8]))
    Nmosaic = str(len(masterdat['id'][selection][(masterdat['id'][selection] > 6e8) & (masterdat['id'][selection] < 7e8)]))
    Nudf10  = str(len(masterdat['id'][selection][masterdat['id'][selection] > 7e8]))

    Nwidetot   = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] < 5e8]))
    Nmosaictot = str(len(masterdat['id'][selectionTOT][(masterdat['id'][selectionTOT] > 6e8) & (masterdat['id'][selectionTOT] < 7e8)]))
    Nudf10tot  = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] > 7e8]))
    if verbose: print('     Of these '+Nwide+'/'+Nmosaic+'/'+Nudf10+
                      ' ('+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+
                      '%/'+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+
                      '%/'+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'%) '+
                      'are in the Wide/Mosaic/UDF10 fields ')

    if verbose:
        tabstring = ' Cemitter (LAE) & XXzrangeXX & '\
        +str(len(selection))+' & '+str(len(selectionTOT))+' & '+str("%5.2f" % (float(len(selection))/float(len(selectionTOT))*100.0))+'\% & '\
        +str(Nwide)+         ' & '+str(Nwidetot)+         ' & '+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+'\% & '\
        +str(Nmosaic)+       ' & '+str(Nmosaictot)+       ' & '+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+'\% & '\
        +str(Nudf10)+        ' & '+str(Nudf10tot)+        ' & '+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'\% \\\\ '
        print(tabstring)

    detectiondic['100fields_CIIIorCIV_LAEs'] = [len(selection)]
    detectiondic['musewide_CIIIorCIV_LAEs']  = [int(Nwide)]
    detectiondic['udfmosaic_CIIIorCIV_LAEs'] = [int(Nmosaic)]
    detectiondic['udf10_CIIIorCIV_LAEs']     = [int(Nudf10)]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    f_linelist = ['ferr_'+rn for rn in rangename[4:11]]
    f_zrange   = zranges[4:11]
    for ll, linesel in enumerate(f_linelist):
        linename = linesel.split('err_')[-1]
        if verbose: print('   o '+linename+' doublet emitters ')
        selection = np.where(((np.abs(masterdat[linesel])  != 99.0) & np.isfinite(masterdat[linesel])) &
                              (masterdat['duplicationID'] == 0.0) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) )[0]

        selectionTOT = np.where((masterdat['redshift'] >= f_zrange[ll][0]) & (masterdat['redshift'] <= f_zrange[ll][1]) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) &
                                (masterdat['duplicationID'] == 0.0) )[0]

        if len(selection) > 0:
            outname_bases.append('/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/mastercat_'+
                                 linename+'Emitters.tex')
            goodents.append(selection)
        if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
        if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
        Nwide   = str(len(masterdat['id'][selection][masterdat['id'][selection] < 5e8]))
        Nmosaic = str(len(masterdat['id'][selection][(masterdat['id'][selection] > 6e8) & (masterdat['id'][selection] < 7e8)]))
        Nudf10  = str(len(masterdat['id'][selection][masterdat['id'][selection] > 7e8]))

        Nwidetot   = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] < 5e8]))
        Nmosaictot = str(len(masterdat['id'][selectionTOT][(masterdat['id'][selectionTOT] > 6e8) & (masterdat['id'][selectionTOT] < 7e8)]))
        Nudf10tot  = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] > 7e8]))
        if verbose: print('     Of these '+Nwide+'/'+Nmosaic+'/'+Nudf10+
                          ' ('+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+
                          '%/'+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+
                          '%/'+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'%) '+
                          'are in the Wide/Mosaic/UDF10 fields ')

        if verbose:
            tabstring = ' '+linesel+'(all) & '+str(f_zrange[ll])+' & '\
            +str(len(selection))+' & '+str(len(selectionTOT))+' & '+str("%5.2f" % (float(len(selection))/float(len(selectionTOT))*100.0))+'\% & '\
            +str(Nwide)+         ' & '+str(Nwidetot)+         ' & '+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+'\% & '\
            +str(Nmosaic)+       ' & '+str(Nmosaictot)+       ' & '+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+'\% & '\
            +str(Nudf10)+        ' & '+str(Nudf10tot)+        ' & '+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'\% \\\\ '
            print(tabstring)

        detectiondic['100fields_'+linesel.split('_')[-1]+'_full'] = [len(selection)]
        detectiondic['musewide_'+linesel.split('_')[-1]+'_full']  = [int(Nwide)]
        detectiondic['udfmosaic_'+linesel.split('_')[-1]+'_full'] = [int(Nmosaic)]
        detectiondic['udf10_'+linesel.split('_')[-1]+'_full']     = [int(Nudf10)]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    f_linelist = ['ferr_'+rn for rn in rangename[4:10]]
    f_zrange   = zranges[4:10]
    for ll, linesel in enumerate(f_linelist):
        linename = linesel.split('err_')[-1]
        if linename == 'NV':
            if verbose: print('   o Skipping '+linename+' doublet emitters for Desert objects ')
            continue
        if verbose: print('   o '+linename+' doublet emitters in z-desert (1.5<z<2.9)')
        selection = np.where(((np.abs(masterdat[linesel])  != 99.0) & np.isfinite(masterdat[linesel])) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) &
                                (masterdat['redshift'] >= 1.5) & (masterdat['redshift'] < 2.9) & (masterdat['duplicationID'] == 0.0) )[0]

        selectionTOT = np.where((masterdat['redshift'] >= f_zrange[ll][0]) & (masterdat['redshift'] <= f_zrange[ll][1]) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) &
                                (masterdat['redshift'] >= 1.5) & (masterdat['redshift'] < 2.9) & (masterdat['duplicationID'] == 0.0) )[0]

        if len(selection) > 0:
            outname_bases.append('/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/mastercat_DesertWith'+
                                 linename+'Detection.tex')
            goodents.append(selection)
        if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
        if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
        Nwide   = str(len(masterdat['id'][selection][masterdat['id'][selection] < 5e8]))
        Nmosaic = str(len(masterdat['id'][selection][(masterdat['id'][selection] > 6e8) & (masterdat['id'][selection] < 7e8)]))
        Nudf10  = str(len(masterdat['id'][selection][masterdat['id'][selection] > 7e8]))

        Nwidetot   = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] < 5e8]))
        Nmosaictot = str(len(masterdat['id'][selectionTOT][(masterdat['id'][selectionTOT] > 6e8) & (masterdat['id'][selectionTOT] < 7e8)]))
        Nudf10tot  = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] > 7e8]))
        if verbose: print('     Of these '+Nwide+'/'+Nmosaic+'/'+Nudf10+
                          ' ('+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+
                          '%/'+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+
                          '%/'+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'%) '+
                          'are in the Wide/Mosaic/UDF10 fields ')

        if verbose:
            tabstring = ' '+linesel+'(Desert) & '+str(f_zrange[ll])+' & '\
            +str(len(selection))+' & '+str(len(selectionTOT))+' & '+str("%5.2f" % (float(len(selection))/float(len(selectionTOT))*100.0))+'\% & '\
            +str(Nwide)+         ' & '+str(Nwidetot)+         ' & '+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+'\% & '\
            +str(Nmosaic)+       ' & '+str(Nmosaictot)+       ' & '+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+'\% & '\
            +str(Nudf10)+        ' & '+str(Nudf10tot)+        ' & '+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'\% \\\\ '
            print(tabstring)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    f_linelist = ['ferr_'+rn for rn in rangename[4:10]]
    f_zrange   = zranges[4:10]
    for ll, linesel in enumerate(f_linelist):
        linename = linesel.split('err_')[-1]
        if verbose: print('   o '+linename+' doublet emitters with Lya (z>2.9)')
        selection = np.where(((np.abs(masterdat[linesel])  != 99.0) & np.isfinite(masterdat[linesel])) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) &
                                (masterdat['redshift'] >= 2.9) & (masterdat['duplicationID'] == 0.0) )[0]

        selectionTOT = np.where((masterdat['redshift'] >= f_zrange[ll][0]) & (masterdat['redshift'] <= f_zrange[ll][1]) &
                                 (masterdat['id'] != 158002004) &
                                 (masterdat['id'] != 601931670) &
                                 (masterdat['id'] != 208014258) &
                                 (masterdat['id'] != 600341002) &
                                (masterdat['redshift'] >= 2.9) & (masterdat['duplicationID'] == 0.0) )[0]

        if len(selection) > 0:
            outname_bases.append('/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/mastercat_LAEsWith'+
                                 linename+'Detection.tex')
            goodents.append(selection)
        if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
        if verbose: print('     selection satisfied by '+str(len(selection))+' objects')
        Nwide   = str(len(masterdat['id'][selection][masterdat['id'][selection] < 5e8]))
        Nmosaic = str(len(masterdat['id'][selection][(masterdat['id'][selection] > 6e8) & (masterdat['id'][selection] < 7e8)]))
        Nudf10  = str(len(masterdat['id'][selection][masterdat['id'][selection] > 7e8]))

        Nwidetot   = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] < 5e8]))
        Nmosaictot = str(len(masterdat['id'][selectionTOT][(masterdat['id'][selectionTOT] > 6e8) & (masterdat['id'][selectionTOT] < 7e8)]))
        Nudf10tot  = str(len(masterdat['id'][selectionTOT][masterdat['id'][selectionTOT] > 7e8]))
        if verbose: print('     Of these '+Nwide+'/'+Nmosaic+'/'+Nudf10+
                          ' ('+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+
                          '%/'+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+
                          '%/'+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'%) '+
                          'are in the Wide/Mosaic/UDF10 fields ')

        if verbose:
            tabstring = ' '+linesel+'(LAE) & '+str(f_zrange[ll])+' & '\
            +str(len(selection))+' & '+str(len(selectionTOT))+' & '+str("%5.2f" % (float(len(selection))/float(len(selectionTOT))*100.0))+'\% & '\
            +str(Nwide)+         ' & '+str(Nwidetot)+         ' & '+str("%5.2f" % (float(Nwide)/float(Nwidetot)*100.))+'\% & '\
            +str(Nmosaic)+       ' & '+str(Nmosaictot)+       ' & '+str("%5.2f" % (float(Nmosaic)/float(Nmosaictot)*100.))+'\% & '\
            +str(Nudf10)+        ' & '+str(Nudf10tot)+        ' & '+str("%5.2f" % (float(Nudf10)/float(Nudf10tot)*100.))+'\% \\\\ '
            print(tabstring)

        detectiondic['100fields_'+linesel.split('_')[-1]+'_LAEs'] = [len(selection)]
        detectiondic['musewide_'+linesel.split('_')[-1]+'_LAEs']  = [int(Nwide)]
        detectiondic['udfmosaic_'+linesel.split('_')[-1]+'_LAEs'] = [int(Nmosaic)]
        detectiondic['udf10_'+linesel.split('_')[-1]+'_LAEs']     = [int(Nudf10)]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Nall = len(masterdat['redshift'][(masterdat['redshift']>0) &
    #                                  (masterdat['redshift']<4.9699) &
    #                                  (masterdat['duplicationID'] == 0) &
    #                                  (masterdat['id'] != 158002004) &
    #                                  (masterdat['id'] != 601931670) &
    #                                  (masterdat['id'] != 208014258) &
    #                                  (masterdat['id'] != 600341002)])
    # NLAE = len(masterdat['redshift'][(masterdat['redshift']>2.9) &
    #                                  (masterdat['redshift']<4.9699) &
    #                                  (masterdat['duplicationID'] == 0) &
    #                                  (masterdat['id'] != 158002004) &
    #                                  (masterdat['id'] != 601931670) &
    #                                  (masterdat['id'] != 208014258) &
    #                                  (masterdat['id'] != 600341002)])
    # print('\n All/LAEs  = '+str(Nall)+'/'+str(NLAE))

    if verbose: print('-------------------- Number counts in all 100 arcmin2 --------------------')
    for zz, zr in enumerate(zranges):
        Nall = len(masterdat['redshift'][(masterdat['redshift']>zr[0]) &
                                         (masterdat['redshift']<zr[1]) &
                                         (masterdat['duplicationID'] == 0) &
                                         (masterdat['id'] != 158002004) &
                                         (masterdat['id'] != 601931670) &
                                         (masterdat['id'] != 208014258) &
                                         (masterdat['id'] != 600341002)])
        NLAE = len(masterdat['redshift'][(masterdat['redshift']>np.max([zr[0],2.9])) &
                                         (masterdat['redshift']<zr[1]) &
                                         (masterdat['duplicationID'] == 0) &
                                         (masterdat['id'] != 158002004) &
                                         (masterdat['id'] != 601931670) &
                                         (masterdat['id'] != 208014258) &
                                         (masterdat['id'] != 600341002)])

        Nallwd = len(masterdat['redshift'][(masterdat['redshift']>zr[0]) &
                                           (masterdat['redshift']<zr[1]) &
                                           (masterdat['duplicationID'] == 0) &
                                           (masterdat['id'] != 158002004) &
                                           (masterdat['id'] != 601931670) &
                                           (masterdat['id'] != 208014258) &
                                           (masterdat['id'] != 600341002)])
        NLAEwd = len(masterdat['redshift'][(masterdat['redshift']>np.max([zr[0],2.9])) &
                                           (masterdat['redshift']<zr[1]) &
                                           (masterdat['duplicationID'] == 0) &
                                           (masterdat['id'] != 158002004) &
                                           (masterdat['id'] != 601931670) &
                                           (masterdat['id'] != 208014258) &
                                           (masterdat['id'] != 600341002)])

        if verbose: print(' [All/LAEs]_'+str("%-15s" % rangename[zz])+'  = '+str("%-4s" % Nall)+'/'+str("%-4s" % NLAE)+
                          '  (with duplication '+str("%-4s" % Nallwd)+'/'+str("%-4s" % NLAEwd)+')')

        detectiondic['100fields_'+rangename[zz]+'_full'].append(Nall)
        detectiondic['100fields_'+rangename[zz]+'_LAEs'].append(NLAE)

    if verbose: print('-------------------- Number counts in musewide (CDFS+PAR+COSMOS) --------------------')
    for zz, zr in enumerate(zranges):
        Nall = len(masterdat['redshift'][(masterdat['redshift']>zr[0]) &
                                         (masterdat['redshift']<zr[1]) &
                                         (masterdat['duplicationID'] == 0) &
                                         (np.asarray([int(str(mid)[0]) for mid in masterdat['id']]) < 5) &
                                         (masterdat['id'] != 158002004) &
                                         (masterdat['id'] != 601931670) &
                                         (masterdat['id'] != 208014258) &
                                         (masterdat['id'] != 600341002)])
        NLAE = len(masterdat['redshift'][(masterdat['redshift']>np.max([zr[0],2.9])) &
                                         (masterdat['redshift']<zr[1]) &
                                         (masterdat['duplicationID'] == 0) &
                                         (np.asarray([int(str(mid)[0]) for mid in masterdat['id']]) < 5) &
                                         (masterdat['id'] != 158002004) &
                                         (masterdat['id'] != 601931670) &
                                         (masterdat['id'] != 208014258) &
                                         (masterdat['id'] != 600341002)])

        Nallwd = len(masterdat['redshift'][(masterdat['redshift']>zr[0]) &
                                           (masterdat['redshift']<zr[1]) &
                                           (masterdat['duplicationID'] == 0) &
                                           (np.asarray([int(str(mid)[0]) for mid in masterdat['id']]) < 5) &
                                           (masterdat['id'] != 158002004) &
                                           (masterdat['id'] != 601931670) &
                                           (masterdat['id'] != 208014258) &
                                           (masterdat['id'] != 600341002)])
        NLAEwd = len(masterdat['redshift'][(masterdat['redshift']>np.max([zr[0],2.9])) &
                                           (masterdat['redshift']<zr[1]) &
                                           (masterdat['duplicationID'] == 0) &
                                           (np.asarray([int(str(mid)[0]) for mid in masterdat['id']]) < 5) &
                                           (masterdat['id'] != 158002004) &
                                           (masterdat['id'] != 601931670) &
                                           (masterdat['id'] != 208014258) &
                                           (masterdat['id'] != 600341002)])

        if verbose: print(' [All/LAEs]_'+str("%-15s" % rangename[zz])+'  = '+str("%-4s" % Nall)+'/'+str("%-4s" % NLAE)+
                          '  (with duplication '+str("%-4s" % Nallwd)+'/'+str("%-4s" % NLAEwd)+')')
        detectiondic['musewide_'+rangename[zz]+'_full'].append(Nall)
        detectiondic['musewide_'+rangename[zz]+'_LAEs'].append(NLAE)

    for ii, idst in enumerate(idstart):
        if verbose: print('-------------------- Number counts in '+idstartname[ii]+' --------------------')
        for zz, zr in enumerate(zranges):
            Nall = len(masterdat['redshift'][(masterdat['redshift']>zr[0]) &
                                             (masterdat['redshift']<zr[1]) &
                                             (masterdat['duplicationID'] == 0) &
                                             (np.asarray([int(str(mid)[0]) for mid in masterdat['id']]) == idst) &
                                             (masterdat['id'] != 158002004) &
                                             (masterdat['id'] != 601931670) &
                                             (masterdat['id'] != 208014258) &
                                             (masterdat['id'] != 600341002)])
            NLAE = len(masterdat['redshift'][(masterdat['redshift']>np.max([zr[0],2.9])) &
                                             (masterdat['redshift']<zr[1]) &
                                             (masterdat['duplicationID'] == 0) &
                                             (np.asarray([int(str(mid)[0]) for mid in masterdat['id']]) == idst) &
                                             (masterdat['id'] != 158002004) &
                                             (masterdat['id'] != 601931670) &
                                             (masterdat['id'] != 208014258) &
                                             (masterdat['id'] != 600341002)])

            Nallwd = len(masterdat['redshift'][(masterdat['redshift']>zr[0]) &
                                               (masterdat['redshift']<zr[1]) &
                                               # (masterdat['duplicationID'] == 0) &
                                               (np.asarray([int(str(mid)[0]) for mid in masterdat['id']]) == idst) &
                                               (masterdat['id'] != 158002004) &
                                               (masterdat['id'] != 601931670) &
                                               (masterdat['id'] != 208014258) &
                                               (masterdat['id'] != 600341002)])
            NLAEwd = len(masterdat['redshift'][(masterdat['redshift']>np.max([zr[0],2.9])) &
                                               (masterdat['redshift']<zr[1]) &
                                               # (masterdat['duplicationID'] == 0) &
                                               (np.asarray([int(str(mid)[0]) for mid in masterdat['id']]) == idst) &
                                               (masterdat['id'] != 158002004) &
                                               (masterdat['id'] != 601931670) &
                                               (masterdat['id'] != 208014258) &
                                               (masterdat['id'] != 600341002)])

            if verbose: print(' [All/LAEs]_'+str("%-15s" % rangename[zz])+'  = '+str("%-4s" % Nall)+'/'+str("%-4s" % NLAE)+
                              '  (with duplication '+str("%-4s" % Nallwd)+'/'+str("%-4s" % NLAEwd)+')')

            detectiondic[idstartname[ii]+'_'+rangename[zz]+'_full'].append(Nall)
            detectiondic[idstartname[ii]+'_'+rangename[zz]+'_LAEs'].append(NLAE)


    if Nkeys_init != len(detectiondic.keys()):
        print(' WARNING The intial number of keys in the detection dictionary has changed... ')
        print('         Faulty key/entry added along the way... setting trace to investigate')
        pdb.set_trace()
    else:
        print(' - - - - - - - - - - - Stats from detection dictionary - - - - - - - - - - -')
        for key in detectiondic.keys():
            if verbose:
                if detectiondic[key][1] == 0:
                    print(str("%-30s" % key)+str("%-20s" % detectiondic[key])+
                          str("%20.4f" % (float(detectiondic[key][0])))+'%')
                else:
                    print(str("%-30s" % key)+str("%-20s" % detectiondic[key])+
                          str("%20.4f" % (float(detectiondic[key][0])/float(detectiondic[key][1])*100.))+'%')

                # if detectiondic[key] == [np.nan]:
                #     print(str("%-20s" % key)+'[np.nan]')
                # else:
                #     print(str("%-20s" % key)+str(detectiondic[key]))

        uves.plot_detectiondic(detectiondic,outdir='/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/',
                               subsel='LAEs',verbose=verbose)
        uves.plot_detectiondic(detectiondic,outdir='/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/',
                               subsel='full',verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Looping over object lists and column sets to generate individual tables ')
    for gg, goodent in enumerate(goodents):
        onb = outname_bases[gg]
        for cc, columns in enumerate(columnslist):
            latexout  = onb.replace('.tex','_'+tabstrings[cc]+'.tex')
            sortindex = np.argsort(masterdat[sortcol][goodent])
            ids       = masterdat['id'][goodent][sortindex]
            uves.mastercat_latextable(latexout,masterdat,infodat,columns,ids,overwrite=overwrite,verbose=verbose)

    return detectiondic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def mastercat_latextable(latexoutput,masterdat,infodat,columns,ids,overwrite=False,verbose=True):
    """
    Generate a LaTeX table based on a master catalog restricting to certain columns and IDs

    --- Example of use ---
    see uves.mastercat_latextable_wrappers()

    """
    if os.path.isfile(latexoutput) & (overwrite == False):
        sys.exit(' The requested output ('+latexoutput+') exists and overwrite=False ')

    if verbose: print(' - Putting together table...')
    Nrow     = len(ids)

    Ncol     = 1
    colnames = 'ID & '
    for col in columns:
        if ('f_' in col):
            colnames = colnames+(uves.mastercat_latextable_colnames(col)+'   &   ')+\
                       (uves.mastercat_latextable_colnames(col.replace('_','err_'))+'   &   ')
            Ncol = Ncol + 2
        else:
            colnames = colnames+(uves.mastercat_latextable_colnames(col)+'   &   ')
            Ncol = Ncol + 1

    posstr   = 'r'*Ncol
    fout     = open(latexoutput,'w')
    fout.write("""
\\begin{table*}
\caption{\label{tab:tablelable}Table Title/Description - for extensive caption use footer below.}
\centering
\\resizebox{\\textwidth}{!}{ %% command scaling table to textwidth
\\begin{tabular}{%s}
\hline\hline
%s
\hline
""" % (posstr,colnames[:-5]+'\\\\'))

    for objid in ids:
        row_string    = str(objid)+'  & '
        objent_info   = np.where(infodat['id'] == objid)[0]
        objent_master = np.where(masterdat['id'] == objid)[0]

        for col in columns:
            if len(objent_info) == 1:
                try:
                    colval = infodat[col][objent_info]
                except:
                    colval = masterdat[col][objent_master]
            else:
                colval = masterdat[col][objent_master]

            if col in ['ra','dec']:
                row_string = row_string+str("%16.8f" % colval)+'  &  '
            elif col in ['redshift']:
                row_string = row_string+str("%12.4f" % colval)+'  &  '
            elif ('f_' in col):
                try:
                    errval = infodat[col.replace('_','err_')][objent_info]
                except:
                    errval = masterdat[col.replace('_','err_')][objent_master]

                if errval == 99.0:
                    row_string = row_string+'$<$'+str("%.2f" % colval)+'  &    &  '
                elif errval == -99.0:
                    row_string = row_string+'$>$'+str("%.2f" % colval)+'  &    &  '
                else:
                    row_string = row_string+str("%.2f" % colval)+'  &  '+str("%12.2f" % errval)+'  &  '
            elif ('name' in col) or ('reference' in col):
                row_string = row_string+str(colval[0].replace('_','\_'))+'  &  '
            else:
                row_string = row_string+str("%.2f" % colval)+'  &  '

        fout.write(row_string.replace('nan ','- ')[:-5]+' \\\\ \n')
    fout.write("""\hline
\end{tabular}
}
\\tablefoot{Main text in footer\\\\
\\tablefoottext{a}{Some information on values in the table}\\\\
\\tablefoottext{b}{Some more info on table values}\\\\
}
\end{table*}
""")

    if verbose: print(' - Wrote latex table of '+str(Ncol)+' columns and '+str(Nrow)+' rows to \n   '+latexoutput)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def mastercat_latextable_colnames(colname):
    """
    Formatting column names for LaTeX

    """
    linenames = {}
    linenames['lya']  = 'Ly$\\alpha$'
    linenames['nv']   = 'NV'
    linenames['nv1']  = 'NV$\\lambda$1239'
    linenames['nv2']  = 'NV$\\lambda$1243'
    linenames['civ']  = 'CIV'
    linenames['civ1'] = 'CIV$\\lambda$1548'
    linenames['civ2'] = 'CIV$\\lambda$1551'
    linenames['heii'] = 'HeII'
    linenames['oiii']  = 'OIII'
    linenames['oiii1'] = 'OIII$\\lambda$1661'
    linenames['oiii2'] = 'OIII$\\lambda$1666'
    linenames['siiii']  = 'SiIII'
    linenames['siiii1'] = 'SiIII$\\lambda$1883'
    linenames['siiii2'] = 'SiIII$\\lambda$1892'
    linenames['ciii']  = 'CIII'
    linenames['ciii1'] = 'CIII$\\lambda$1907'
    linenames['ciii2'] = 'CIII$\\lambda$1909'
    linenames['mgii']  = 'MgII'
    linenames['mgii1'] = 'MgII$\\lambda$2796'
    linenames['mgii2'] = 'MgII$\\lambda$2803'

    if colname.lower() == 'id':
        colname_fmt = 'ID$_\\textrm{MUSE-Wide/Deep}$'
    elif colname.lower() == 'redshift':
        colname_fmt = '$z_\\textrm{Ly$\\alpha$}$'
    elif colname.lower() == 'ra':
        colname_fmt = 'R.A.'
    elif colname.lower() == 'dec':
        colname_fmt = 'Dec.'
    elif colname.lower() == 'ew_0':
        colname_fmt = 'EW$_{0,\\textrm{Ly}\\alpha}$'
    elif colname.lower() == 'ew_0_err':
        colname_fmt = '$\\delta$EW$_{0,\\textrm{Ly}\\alpha}$'
    elif colname.lower().startswith('id_'):
        colname_fmt = 'ID$_\\textrm{'+colname.split('_')[-1][0].upper()+colname.split('_')[-1][1:].lower()+'}$'
    elif colname.lower().startswith('sep_'):
        colname_fmt = '$\\Delta_\\textrm{'+colname.split('_')[-1][0].upper()+colname.split('_')[-1][1:].lower()+'}$["]'
    elif colname.lower().startswith('f_'):
        linename = linenames[colname.split('_')[-1].lower()]
        colname_fmt = '$f_\\textrm{'+linename+'}$'
    elif colname.lower().startswith('ferr_'):
        linename = linenames[colname.split('_')[-1].lower()]
        colname_fmt = '$\\delta f_\\textrm{'+linename+'}$'
    elif colname.lower().startswith('sigma_'):
        linename = linenames[colname.split('_')[-1].lower()]
        colname_fmt = '$\sigma_\\textrm{'+linename+'}$'
    elif colname.lower().startswith('sigmaerr_'):
        linename = linenames[colname.split('_')[-1].lower()]
        colname_fmt = '$\\delta\sigma_\\textrm{'+linename+'}$'
    elif colname.lower().startswith('fr'):
        subscript = colname.split('_')[-1].lower()
        for ln in linenames.keys():
            if subscript.startswith(ln):
                if subscript.split(ln)[1].startswith('1'):
                    numerator = ln+'1'
                elif subscript.split(ln)[1].startswith('2'):
                    numerator = ln+'2'
                else:
                    numerator = ln

                denominator = subscript.split(numerator)[-1]

        ratiostring = '\\frac{f_\\textrm{'+linenames[numerator]+'}}{f_\\textrm{'+linenames[denominator]+'}}'
        if colname.lower().startswith('fr_'):
            colname_fmt = '$'+ratiostring+'$'
        else:
            colname_fmt = '$\\delta\left('+ratiostring+'\\right)$'
    else:
        colname_fmt = colname.replace('_','\_')

    return colname_fmt

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def literaturecat_latextable_wrapper(sortcol='id',overwrite=False,verbose=True):
    """
    Wrapper to uves.mastercat_latextable() generating a table of the emitters with either CIII or CIV detections.
    This wraper also returns a dictionary with number counts of detections and available objects for searching
    based on the content of the master catalog.

    --- Example of use ---
    import uvEmissionlineSearch as uves
    uves.literaturecat_latextable_wrapper(verbose=True,overwrite=True,sortcol='id')

    """
    # infofile  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    infodat   = afits.open(infofile)[1].data
    workdir   = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/'
    litcat    = workdir+'literaturecollection_emissionlinestrengths.fits'
    litdat    = afits.open(litcat)[1].data
    columns   = ['name','reference','ra','dec','redshift','f_Lya','ferr_Lya','EW0_Lya','EW0err_Lya',
                 'f_CIV','ferr_CIV','f_HeII','ferr_HeII','f_OIII','ferr_OIII','f_SiIII','ferr_SiIII','f_CIII','ferr_CIII',
                 'EW0_CIV','EW0err_CIV','EW0_HeII','EW0err_HeII','EW0_OIII','EW0err_OIII',
                 'EW0_SiIII','EW0err_SiIII','EW0_CIII','EW0err_CIII']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    latexout  = workdir+'literaturecollection_emissionlinestrengths_fluxtable.tex'
    if verbose: print(' - Putting together the literature collection latex table and storing it to\n   '+latexout)
    sortindex = np.argsort(litdat[sortcol])
    ids       = litdat['id'][sortindex]
    uves.mastercat_latextable(latexout,litdat,infodat,columns,ids,overwrite=overwrite,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # replace caption with list of references:
    refshorts = np.unique(litdat['reference'])
    refexplainer = ' '
    for refshort in refshorts:
        reflong  = lce.referencedictionary()[refshort][1]
        refexplainer = refexplainer+refshort+':'+reflong.replace('&','\&')+', '

    with open(latexout, 'r+') as fdata:
        text = fdata.read()
        text = re.sub('Main text in footer', 'References are'+refexplainer[:-2]+'.', text)
        fdata.seek(0)
        fdata.write(text)
        fdata.truncate()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def convert_FoVcoordinates():
    """
    import uvEmissionlineSearch as uves
    uves.convert_FoVcoordinates()

    """
    fieldcoord_gs = np.array([[53.062397   ,-27.80815506],[53.06840134 ,-27.82277679],[53.07440948 ,-27.83739662],[53.08042145 ,-27.85201454],[53.08643341 ,-27.86663437],[53.07892227 ,-27.80284119],[53.08493423 ,-27.81746101],[53.09094238 ,-27.83208084],[53.09695435 ,-27.84669876],[53.10297012 ,-27.86131859],[53.09545135 ,-27.79752731],[53.1014595  ,-27.81214523],[53.10747528 ,-27.82676315],[53.11348724 ,-27.84138107],[53.11950302 ,-27.85599899],[53.13603592 ,-27.8506794 ],[53.15256882 ,-27.84535599],[53.1690979  ,-27.84003258],[53.18562698 ,-27.83470535],[53.20215225 ,-27.82937813],[53.21867752 ,-27.82404709],[53.13002014 ,-27.83606148],[53.14654922 ,-27.83073997],[53.16307449 ,-27.82541656],[53.17960358 ,-27.82009125],[53.19612503 ,-27.81476212],[53.12400436 ,-27.82144356],[53.11798859 ,-27.80682755],[53.11197662 ,-27.79220963],[53.10596466 ,-27.7775898 ],[53.0999527  ,-27.76297188],[53.14052963 ,-27.81612396],[53.13451385 ,-27.80150604],[53.12849808 ,-27.78689003],[53.12248611 ,-27.77227211],[53.11647415 ,-27.75765419],[53.13900375 ,-27.76695061],[53.13299179 ,-27.75233459],[53.08943939 ,-27.78290749],[53.08343124 ,-27.76828766],[53.07291794 ,-27.78822327],[53.05038452 ,-27.77891541],[53.18253708 ,-27.73636246],[53.19457626 ,-27.76559067],[53.20059967 ,-27.78020477],[53.06690979 ,-27.77360344],[53.05638885 ,-27.79353523],[53.14950943 ,-27.74701118],[53.16602325 ,-27.74168777],[53.18855667 ,-27.75097656],[53.20662308 ,-27.79481888],[53.17203903 ,-27.75630379],[53.19010162 ,-27.80014801],[53.03385925 ,-27.78422546],[53.03986359 ,-27.7988472 ],[53.04586411 ,-27.81346703],[53.05187225 ,-27.82808876],[53.05787659 ,-27.84270859],[53.06388474 ,-27.85732841],[53.0698967  ,-27.87195015],[53.26665497 ,-27.68647957],[53.25360107 ,-27.67607117],[53.24184418 ,-27.68763161],[53.25489807 ,-27.69804001],[53.27975082 ,-27.85684013],[53.2653389  ,-27.84791946],[53.25524521 ,-27.86066055],[53.26965714 ,-27.86958313]])
    textcoord_gs  = np.array([[53.062397  , -27.80815506],[53.06840134, -27.82277679],[53.07440948, -27.83739662],[53.08042145, -27.85201454],[53.08643341, -27.86663437],[53.07892227, -27.80284119],[53.08493423, -27.81746101],[53.09094238, -27.83208084],[53.09695435, -27.84669876],[53.10297012, -27.86131859],[53.09545135, -27.79752731],[53.1014595 , -27.81214523],[53.10747528, -27.82676315],[53.11348724, -27.84138107],[53.11950302, -27.85599899],[53.13603592, -27.8506794 ],[53.15256882, -27.84535599],[53.1690979 , -27.84003258],[53.18562698, -27.83470535],[53.20215225, -27.82937813],[53.21867752, -27.82404709],[53.13002014, -27.83606148],[53.14654922, -27.83073997],[53.16307449, -27.82541656],[53.17960358, -27.82009125],[53.19612503, -27.81476212],[53.12400436, -27.82144356],[53.11798859, -27.80682755],[53.11197662, -27.79220963],[53.10596466, -27.7775898 ],[53.0999527 , -27.76297188],[53.14052963, -27.81612396],[53.13451385, -27.80150604],[53.12849808, -27.78689003],[53.12248611, -27.77227211],[53.11647415, -27.75765419],[53.13900375, -27.76695061],[53.13299179, -27.75233459],[53.08943939, -27.78290749],[53.08343124, -27.76828766],[53.07291794, -27.78822327],[53.05038452, -27.77891541],[53.18253708, -27.73636246],[53.19457626, -27.76559067],[53.20059967, -27.78020477],[53.06690979, -27.77360344],[53.05638885, -27.79353523],[53.14950943, -27.74701118],[53.16602325, -27.74168777],[53.18855667, -27.75097656],[53.20662308, -27.79481888],[53.17203903, -27.75630379],[53.19010162, -27.80014801],[53.03385925, -27.78422546],[53.03986359, -27.7988472 ],[53.04586411, -27.81346703],[53.05187225, -27.82808876],[53.05787659, -27.84270859],[53.06388474, -27.85732841],[53.0698967 , -27.87195015],[53.26665497, -27.68647957],[53.25360107, -27.67607117],[53.24184418, -27.68763161],[53.25489807, -27.69804001],[53.27975082, -27.85684013],[53.2653389 , -27.84791946],[53.25524521, -27.86066055],[53.26965714, -27.86958313],[53.15960312, -27.76493263],[53.17351151, -27.77551079],[53.18653107, -27.78676605],[53.14808273, -27.77723885],[53.16083908, -27.78763008],[53.17481613, -27.79899406],[53.13550186, -27.78861046],[53.14797974, -27.79980278],[53.16194153, -27.81017494]])

    coord_cos = np.array([[150.09121704 , 2.20111012],[150.09121704 , 2.20111012],[150.10679626 , 2.2011106 ],[150.10679626 , 2.2011106 ],[150.12236023 , 2.20111084],[150.12236023 , 2.20111084],[150.13792419 , 2.20111084],[150.13792419 , 2.20111084],[150.15348816 , 2.20111084],[150.15348816 , 2.20111084],[150.16905212 , 2.2011106 ],[150.16905212 , 2.2011106 ],[150.18463135 , 2.20111012],[150.18463135 , 2.20111012],[150.09121704 , 2.21666574],[150.09121704 , 2.21666574],[150.10679626 , 2.21666622],[150.10679626 , 2.21666622],[150.12236023 , 2.21666646],[150.12236023 , 2.21666646],[150.13792419 , 2.21666646],[150.13792419 , 2.21666646],[150.15348816 , 2.21666646],[150.15348816 , 2.21666646],[150.16905212 , 2.21666622],[150.16905212 , 2.21666622],[150.18463135 , 2.21666574],[150.18463135 , 2.21666574],[150.09121704 , 2.23222136],[150.09121704 , 2.23222136],[150.10679626 , 2.2322216 ],[150.10679626 , 2.2322216 ],[150.12236023 , 2.23222184],[150.12236023 , 2.23222184],[150.13792419 , 2.23222208],[150.13792419 , 2.23222208],[150.15348816 , 2.23222184],[150.15348816 , 2.23222184],[150.16905212 , 2.2322216 ],[150.16905212 , 2.2322216 ],[150.18463135 , 2.23222136],[150.18463135 , 2.23222136],[150.09121704 , 2.24777675],[150.09121704 , 2.24777675],[150.10679626 , 2.24777722],[150.10679626 , 2.24777722],[150.12236023 , 2.24777746],[150.12236023 , 2.24777746],[150.13792419 , 2.2477777 ],[150.13792419 , 2.2477777 ],[150.15348816 , 2.24777746],[150.15348816 , 2.24777746],[150.16905212 , 2.24777722],[150.16905212 , 2.24777722],[150.18463135 , 2.24777675],[150.18463135 , 2.24777675],[150.09121704 , 2.26333237],[150.09121704 , 2.26333237],[150.10679626 , 2.26333284],[150.10679626 , 2.26333284],[150.12236023 , 2.26333308],[150.12236023 , 2.26333308],[150.13792419 , 2.26333308],[150.13792419 , 2.26333308],[150.15348816 , 2.26333308],[150.15348816 , 2.26333308],[150.16905212 , 2.26333284],[150.16905212 , 2.26333284],[150.18463135 , 2.26333237],[150.18463135 , 2.26333237],[150.09121704 , 2.27888799],[150.09121704 , 2.27888799],[150.10678101 , 2.27888846],[150.10678101 , 2.27888846],[150.12236023 , 2.2788887 ],[150.12236023 , 2.2788887 ],[150.13792419 , 2.2788887 ],[150.13792419 , 2.2788887 ],[150.15348816 , 2.2788887 ],[150.15348816 , 2.2788887 ],[150.16906738 , 2.27888846],[150.16906738 , 2.27888846],[150.18463135 , 2.27888799],[150.18463135 , 2.27888799],[150.09121704 , 2.29444361],[150.09121704 , 2.29444361],[150.10678101 , 2.29444408],[150.10678101 , 2.29444408],[150.12236023 , 2.29444432],[150.12236023 , 2.29444432],[150.13792419 , 2.29444432],[150.13792419 , 2.29444432],[150.15348816 , 2.29444432],[150.15348816 , 2.29444432],[150.16906738 , 2.29444408],[150.16906738 , 2.29444408],[150.18463135 , 2.29444361],[150.18463135 , 2.29444361],[150.09126282 , 2.30999923],[150.09126282 , 2.30999923],[150.10681152 , 2.3099997 ],[150.10681152 , 2.3099997 ],[150.12237549 , 2.30999994],[150.12237549 , 2.30999994],[150.13792419 , 2.30999994],[150.13792419 , 2.30999994],[150.1534729 ,  2.30999994],[150.1534729 ,  2.30999994],[150.16903687 , 2.3099997 ],[150.16903687 , 2.3099997 ],[150.18458557 , 2.30999923],[150.18458557 , 2.30999923],[150.09121704 , 2.32555485],[150.09121704 , 2.32555485],[150.10678101 , 2.32555509],[150.10678101 , 2.32555509],[150.12236023 , 2.32555556],[150.12236023 , 2.32555556],[150.13792419 , 2.32555556],[150.13792419 , 2.32555556],[150.15348816 , 2.32555556],[150.15348816 , 2.32555556],[150.16906738 , 2.32555509],[150.16906738 , 2.32555509],[150.18463135 , 2.32555485],[150.18463135 , 2.32555485],[150.09121704 , 2.34111023],[150.09121704 , 2.34111023],[150.10678101 , 2.34111071],[150.10678101 , 2.34111071],[150.12236023 , 2.34111094],[150.12236023 , 2.34111094],[150.13792419 , 2.34111118],[150.13792419 , 2.34111118],[150.15348816 , 2.34111094],[150.15348816 , 2.34111094],[150.16906738 , 2.34111071],[150.16906738 , 2.34111071],[150.18463135 , 2.34111023],[150.18463135 , 2.34111023],[150.09121704 , 2.35666585],[150.09121704 , 2.35666585],[150.10678101 , 2.35666633],[150.10678101 , 2.35666633],[150.12236023 , 2.35666656],[150.12236023 , 2.35666656],[150.13792419 , 2.3566668 ],[150.13792419 , 2.3566668 ],[150.15348816 , 2.35666656],[150.15348816 , 2.35666656],[150.16906738 , 2.35666633],[150.16906738 , 2.35666633],[150.18463135 , 2.35666585],[150.18463135 , 2.35666585],[150.09121704 , 2.37222147],[150.09121704 , 2.37222147],[150.10678101 , 2.37222195],[150.10678101 , 2.37222195],[150.12236023 , 2.37222219],[150.12236023 , 2.37222219],[150.13792419 , 2.37222219],[150.13792419 , 2.37222219],[150.15348816 , 2.37222219],[150.15348816 , 2.37222219],[150.16906738 , 2.37222195],[150.16906738 , 2.37222195],[150.18463135 , 2.37222147],[150.18463135 , 2.37222147],[150.09121704 , 2.38777709],[150.09121704 , 2.38777709],[150.10678101 , 2.38777757],[150.10678101 , 2.38777757],[150.12236023 , 2.38777781],[150.12236023 , 2.38777781],[150.13792419 , 2.38777781],[150.13792419 , 2.38777781],[150.15348816 , 2.38777781],[150.15348816 , 2.38777781],[150.16906738 , 2.38777757],[150.16906738 , 2.38777757],[150.18463135 , 2.38777709],[150.18463135 , 2.38777709],[150.09121704 , 2.40333271],[150.09121704 , 2.40333271],[150.10678101 , 2.40333295],[150.10678101 , 2.40333295],[150.12236023 , 2.40333343],[150.12236023 , 2.40333343],[150.13792419 , 2.40333343],[150.13792419 , 2.40333343],[150.15348816 , 2.40333343],[150.15348816 , 2.40333343],[150.16906738 , 2.40333295],[150.16906738 , 2.40333295],[150.18463135 , 2.40333271],[150.18463135 , 2.40333271]])

    coordsets = [fieldcoord_gs,textcoord_gs,coord_cos]
    for cs in coordsets:
        uves.convert_coordinates_deg2sex(cs)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def convert_coordinates_deg2sex(coordinateset):
    """


    """
    print('------------------------')
    ravals        = acoord.Angle(coordinateset[:,0], u.degree)
    decvals       = acoord.Angle(coordinateset[:,1], u.degree)

    for cc, raval in enumerate(ravals):
        print( str(raval.to_string(unit=u.hour, sep=':'))+' , '+str(decvals[cc].to_string(sep=':')) )

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_detectiondic(detectiondic,outdir='/Users/kschmidt/Desktop/',subsel='LAEs',verbose=True):
    """
    Generate a LaTeX table based on a master catalog restricting to certain columns and IDs

    --- Example of use ---
    see uves.mastercat_latextable_wrappers()

    """


    if verbose: print(' - Setting up and generating plot')
    plotname = outdir+'detectiondicstats_'+subsel+'.pdf'
    fig = plt.figure(figsize=(7, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.97, bottom=0.40, top=0.97)
    Fsize    = 10
    lthick   = 2
    marksize = 6
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(inforstr[:-2],fontsize=Fsize)

    tnames    = []
    showkeys  = ['100fields_allobjUVrange_'+subsel,
                 'musewide_allobjUVrange_'+subsel,
                 'udfmosaic_allobjUVrange_'+subsel,
                 'udf10_allobjUVrange_'+subsel,
                 #----------,
                 '100fields_CIIIorCIV_'+subsel,
                 'musewide_CIIIorCIV_'+subsel,
                 'udfmosaic_CIIIorCIV_'+subsel,
                 'udf10_CIIIorCIV_'+subsel,
                 #----------,
                 '100fields_CIII_'+subsel,
                 'musewide_CIII_'+subsel,
                 'udfmosaic_CIII_'+subsel,
                 'udf10_CIII_'+subsel,
                 #----------,
                 '100fields_CIV_'+subsel,
                 'musewide_CIV_'+subsel,
                 'udfmosaic_CIV_'+subsel,
                 'udf10_CIV_'+subsel,
                 #----------,
                 '100fields_OIII_'+subsel,
                 'musewide_OIII_'+subsel,
                 'udfmosaic_OIII_'+subsel,
                 'udf10_OIII_'+subsel,
                 #----------,
                 '100fields_SiIII_'+subsel,
                 'musewide_SiIII_'+subsel,
                 'udfmosaic_SiIII_'+subsel,
                 'udf10_SiIII_'+subsel,
                 #----------,
                 '100fields_HeII_'+subsel,
                 'musewide_HeII_'+subsel,
                 'udfmosaic_HeII_'+subsel,
                 'udf10_HeII_'+subsel,
                 #----------,
                 '100fields_MgII_'+subsel,
                 'musewide_MgII_'+subsel,
                 'udfmosaic_MgII_'+subsel,
                 'udf10_MgII_'+subsel]

    ticknames = {'100fields_allobjUVrange_'+subsel:'100 Fields allobjUVrange ',
                 'musewide_allobjUVrange_'+subsel:'MUSE-Wide allobjUVrange ',
                 'udfmosaic_allobjUVrange_'+subsel:'UDF mosaic allobjUVrange ',
                 'udf10_allobjUVrange_'+subsel:'UDF-10 allobjUVrange ',
                 #----------,
                 '100fields_CIIIorCIV_'+subsel:'100 Fields CIII or CIV ',
                 'musewide_CIIIorCIV_'+subsel:'MUSE-Wide CIII or CIV ',
                 'udfmosaic_CIIIorCIV_'+subsel:'UDF mosaic CIII or CIV ',
                 'udf10_CIIIorCIV_'+subsel:'UDF-10 CIII or CIV ',
                 #----------,
                 '100fields_CIII_'+subsel:'100 Fields CIII ',
                 'musewide_CIII_'+subsel:'MUSE-Wide CIII ',
                 'udfmosaic_CIII_'+subsel:'UDF mosaic CIII ',
                 'udf10_CIII_'+subsel:'UDF-10 CIII ',
                 #----------,
                 '100fields_CIV_'+subsel:'100 Fields CIV ',
                 'musewide_CIV_'+subsel:'MUSE-Wide CIV ',
                 'udfmosaic_CIV_'+subsel:'UDF mosaic CIV ',
                 'udf10_CIV_'+subsel:'UDF-10 CIV ',
                 #----------,
                 '100fields_OIII_'+subsel:'100 Fields OIII ',
                 'musewide_OIII_'+subsel:'MUSE-Wide OIII ',
                 'udfmosaic_OIII_'+subsel:'UDF mosaic OIII ',
                 'udf10_OIII_'+subsel:'UDF-10 OIII ',
                 #----------,
                 '100fields_SiIII_'+subsel:'100 Fields SiIII ',
                 'musewide_SiIII_'+subsel:'MUSE-Wide SiIII ',
                 'udfmosaic_SiIII_'+subsel:'UDF mosaic SiIII ',
                 'udf10_SiIII_'+subsel:'UDF-10 SiIII ',
                 #----------,
                 '100fields_HeII_'+subsel:'100 Fields HeII ',
                 'musewide_HeII_'+subsel:'MUSE-Wide HeII ',
                 'udfmosaic_HeII_'+subsel:'UDF mosaic HeII ',
                 'udf10_HeII_'+subsel:'UDF-10 HeII ',
                 #----------,
                 '100fields_MgII_'+subsel:'100 Fields MgII ',
                 'musewide_MgII_'+subsel:'MUSE-Wide MgII ',
                 'udfmosaic_MgII_'+subsel:'UDF mosaic MgII ',
                 'udf10_MgII_'+subsel:'UDF-10 MgII '}
    Nobjdet   = []
    Nobjtot   = []
    xvalues   = []
    yvalues   = []

    for ss, dickey in enumerate(showkeys):
        tnames.append(ticknames[dickey])
        Nobjdet.append(detectiondic[dickey][0])
        Nobjtot.append(detectiondic[dickey][1])
        xvalues.append((ss+1))
        if detectiondic[dickey][1] == 0:
            yvalues.append(0.0)
        else:
            yvalues.append(float(detectiondic[dickey][0])/float(detectiondic[dickey][1])*100.)

    xerr    = [None]*len(xvalues)
    yerr    = [None]*len(xvalues)

    plt.plot(xvalues,yvalues,'o',markersize=marksize,alpha=1.0,color='darkgray')
    # plt.errorbar(xvalues,yvalues,xerr=xerr,yerr=yerr,
    #              marker='o',lw=lthick, markersize=marksize,alpha=1.0,
    #              markerfacecolor='blue',ecolor='blue',
    #              markeredgecolor='black',zorder=10)

    plt.ylabel('\% of objects with FELIS UV detection')

    plt.xticks(xvalues, tnames, rotation='vertical')

    #--------- RANGES ---------
    xrange = [0,len(xvalues)+1]

    if subsel == 'LAEs':
        yrange = [0,20]
    else:
        yrange = [0,110]
    plt.xlim(xrange)
    plt.ylim(yrange)

    numbertext = [str(Nobjdet[ii])+'/'+str(Nobjtot[ii]) for ii in np.arange(len(xvalues))]
    for tt, nt in enumerate(numbertext):
        plt.text(xvalues[tt],yrange[1]*0.85,nt, fontsize=Fsize, rotation=90,
                 horizontalalignment='center',verticalalignment='center')

    Nseperators   = 7
    linepositions = (np.arange(Nseperators)+1)*4+0.5
    for linepos in linepositions:
        plt.plot([linepos,linepos],yrange,'--',color='gray',lw=lthick)

    # if logx:
    #     plt.xscale('log')
    # if logy:
    #     plt.yscale('log')

    #--------- LEGEND ---------
    # plt.errorbar(-5000,-5000,xerr=None,yerr=1,marker='o',lw=0, markersize=marksize,alpha=1.0,
    #              markerfacecolor='k',ecolor='k',markeredgecolor='black',zorder=1,label='MW LAE (S/N\_'+key+' $>$ 3)')
    # plt.errorbar(-5000,-5000,xerr=None,yerr=1,marker='o',lw=0, markersize=marksize,alpha=1.0,
    #              markerfacecolor='None',ecolor='k',markeredgecolor='black',zorder=1,label='MW LAE (S/N\_'+key+' $<$ 3)')
    # plt.errorbar(-5000,-5000,xerr=None,yerr=1,marker=r'$\nearrow$',lw=0, markersize=marksize*2,alpha=1.0,
    #              markerfacecolor='None',ecolor='k',markeredgecolor='black',zorder=1,
    #              label='MW LAE (HST non-det. lower limit)')
    #
    #
    # plt.errorbar(-5000,-5000,xerr=None,yerr=None,marker='*',lw=0, markersize=marksize*2,alpha=1.0,
    #              markerfacecolor='None',ecolor='None',markeredgecolor='black',zorder=1,label='AGN')
    # plt.errorbar(-5000,-5000,xerr=None,yerr=None,marker='D',lw=0, markersize=marksize,alpha=1.0,
    #              markerfacecolor='None',ecolor='None',markeredgecolor='black',zorder=1,label='AGN candidate')
    #
    # leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize/1.7},ncol=3,numpoints=1,
    #                  bbox_to_anchor=(0.5, 1.1),)  # add the legend
    # leg.get_frame().set_alpha(0.7)
    #--------------------------

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_detectionFractionsInFields(plotname,verbose=True):
    """
    import uvEmissionlineSearch as uves

    plotname = '/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/table_UVdetections_detectionratios.pdf'
    uves.plot_detectionFractionsInFields(plotname)

    """
    textable = '/Users/kschmidt/work/publications/MUSE_UVemissionlines/tables/table_UVdetections.tex'
    if verbose: print(' - Loading manually (hard-coded) arrays of fractions from \n   '+textable)
    allobjarr = np.array((  [ 5.02  ,  2.18  ,  6.82 , 12.88]  ,
                            [ 4.84  ,  1.90  ,  6.59 , 14.29]  ,
                            [ 0.35  ,  0.27  ,  0.44 ,  0.47]  ,
                            [ 2.63  ,  1.58  ,  2.88 ,  7.51]  ,
                            [ 1.09  ,  0.86  ,  1.62 ,  0.65]  ,
                            [ 1.24  ,  0.25  ,  1.41 ,  5.92]  ,
                            [ 1.30  ,  0.00  ,  2.31 ,  4.07]  ,
                            [ 5.28  ,  0.96  ,  8.45 , 14.75]  ,
                            [58.82  , 25.00  , 68.42 , 54.55]  ))

    allobjarr_Nobj = np.array((  [ 2052  , 1100   ,  719  ,  233 ]  ,
                                 [ 1736  , 947    ,  607  ,  182 ]  ,
                                 [ 1997  , 1094   ,  688  ,  215 ]  ,
                                 [ 1710  , 947    ,  590  ,  173 ]  ,
                                 [ 1465  , 817    ,  495  ,  153 ]  ,
                                 [ 1451  , 803    ,  496  ,  152 ]  ,
                                 [ 1000  , 530    ,  347  ,  123 ]  ,
                                 [ 985   , 520    ,  343  ,  122 ]  ,
                                 [ 34    , 4      ,  19   ,  11  ]  ))

    laearr    = np.array((  [ 2.70  , 1.74  ,  3.05  ,  6.51]  ,
                            [ 2.62  , 1.59  ,  2.95  ,  7.32]  ,
                            [ 0.35  , 0.27  ,  0.44  ,  0.47]  ,
                            [ 2.02  , 1.38  ,  2.08  ,  5.49]  ,
                            [ 0.56  , 0.62  ,  0.63  ,  0.00]  ,
                            [ 0.42  , 0.25  ,  0.00  ,  2.84]  ,
                            [ 0.21  , 0.00  ,  0.32  ,  0.94]  ,
                            [ 1.61  , 0.78  ,  2.24  ,  3.85]  ))

    laearr_Nobj = np.array((  [ 1997  , 1094  , 688  , 215 ]  ,
                              [ 1682  , 942   , 576  , 164 ]  ,
                              [ 1997  , 1094  , 688  , 215 ]  ,
                              [ 1682  , 942   , 576  , 164 ]  ,
                              [ 1431  , 812   , 476  , 143 ]  ,
                              [ 1413  , 798   , 474  , 141 ]  ,
                              [ 947   , 525   , 316  , 106 ]  ,
                              [ 931   , 515   , 312  , 104 ]  ))

    # Any detection		& 103 & 2052 &   5.02\% & 24 & 1100 &  2.18\% & 49 & 719 &  6.82\% & 30 & 233 & 12.88\% \\
    # \ciii{} or \civ{} & 84  & 1736 &   4.84\% & 18 & 947  &  1.90\% & 40 & 607 &  6.59\% & 26 & 182 & 14.29\% \\
    # \nv			    & 7   & 1997 &   0.35\% & 3  & 1094 &  0.27\% & 3  & 688 &  0.44\% & 1  & 215 &  0.47\% \\
    # \civ				& 45  & 1710 &   2.63\% & 15 & 947  &  1.58\% & 17 & 590 &  2.88\% & 13 & 173 &  7.51\% \\
    # \heii             & 16  & 1465 &   1.09\% & 7  & 817  &  0.86\% & 8  & 495 &  1.62\% & 1  & 153 &  0.65\% \\
    # \oiii				& 18  & 1451 &   1.24\% & 2  & 803  &  0.25\% & 7  & 496 &  1.41\% & 9  & 152 &  5.92\% \\
    # \siiii			& 13  & 1000 &   1.30\% & 0  & 530  &  0.00\% & 8  & 347 &  2.31\% & 5  & 123 &  4.07\% \\
    # \ciii				& 52  & 985  &   5.28\% & 5  & 520  &  0.96\% & 29 & 343 &  8.45\% & 18 & 122 & 14.75\% \\
    # \mgii				&  20 & 34   &   58.82\% & 1 & 4    & 25.00\% & 13 & 19  & 68.42\% & 6  & 11  & 54.55\% \\
    #
    # Any detection		& 54 & 1997 &  2.70\% & 19 & 1094 &  1.74\% & 21 & 688 &  3.05\% & 14 & 215 &  6.51\% \\
    # \ciii{} or \civ{} & 44 & 1682 &  2.62\% & 15 & 942  &  1.59\% & 17 & 576 &  2.95\% & 12 & 164 &  7.32\% \\
    # \nv			    & 7  & 1997 &  0.35\% & 3  & 1094 &  0.27\% & 3  & 688 &  0.44\% & 1  & 215 &  0.47\% \\
    # \civ				& 34 & 1682 &  2.02\% & 13 & 942  &  1.38\% & 12 & 576 &  2.08\% & 9  & 164 &  5.49\% \\
    # \heii             &  8 & 1431 &  0.56\% & 5  & 812  &  0.62\% & 3  & 476 &  0.63\% & 0  & 143 &  0.00\% \\
    # \oiii				&  6 & 1413 &  0.42\% & 2  & 798  &  0.25\% & 0  & 474 &  0.00\% & 4  & 141 &  2.84\% \\
    # \siiii			&  2 & 947  &  0.21\% & 0  & 525  &  0.00\% & 1  & 316 &  0.32\% & 1  & 106 &  0.94\% \\
    # \ciii				& 15 & 931  &  1.61\% & 4  & 515  &  0.78\% & 7  & 312 &  2.24\% & 4  & 104 &  3.85\% \\

    arrayrows = {'ANY':0, 'CIIIorCIV':1, 'NV':2,'CIV':3,'HeII':4,'OIII':5,'SiIII':6,'CIII':7}#,'MgII':8}
    # ymaxvals  = {'ANY':16, 'CIIIorCIV':16, 'NV':0.7,'CIV':9.3,'HeII':2.2,'OIII':7.1,'SiIII':5.2,'CIII':16.1}#,'MgII':60}
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Setting up and generating plot:\n   '+plotname)
    fig = plt.figure(figsize=(3, 3))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.97, bottom=0.30, top=0.97)
    Fsize    = 10
    lthick   = 2
    marksize = 6
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(inforstr[:-2],fontsize=Fsize)

    linecolors = uves.linecolors(colormap='plasma')
    linenames  = uves.linenames()

    ticknames = ['All fields','1 hr\nMUSE-Wide','10 hrs\nUDF mosaic','31 hrs\nUDF-10']
    ticknames = ['All fields','MUSE-Wide: 1hr','UDF mosaic: 10hrs','UDF-10: 31hrs']

    xvalues = np.array([5,16,25,46])
    xrange  = [0,50]
    ylabel  = '\% objects with UV line detection'

    # -------------------------- Full frame --------------------------
    axfull = fig.add_subplot(111)
    axfull.set_xlabel('')
    axfull.set_ylabel(ylabel)
    # Turn off axis lines and ticks of the full window subplot
    axfull.spines['top'].set_color('none')
    axfull.spines['bottom'].set_color('none')
    axfull.spines['left'].set_color('none')
    axfull.spines['right'].set_color('none')
    axfull.tick_params(labelcolor='white', top=False, bottom=False, left=False, right=False)
    # plt.axis('off')

    # -------------------------- FULL Z RANGE --------------------------
    ax = fig.add_subplot(2, 1, 1)

    for uvline in linecolors.keys():
        if (uvline != 'MgII') & (uvline != 'CIIIorCIV'): # not showing MgII and CIIIorCIV panels
            yvalues  = allobjarr[arrayrows[uvline],:]
            linecol  = linecolors[uvline]
            plt.plot(xvalues,yvalues,'o',markersize=marksize,alpha=1.0,color=linecol,lw=lthick,zorder=10,label=linenames[uvline])
            plt.plot(xvalues[1:],yvalues[1:],'-',alpha=1.0,color=linecol,lw=lthick,zorder=10)

    ymax = 15
    plt.plot([20,20],[0,ymax],'-',color='black',alpha=1.0,zorder=1,lw=1)

    plt.xticks(xvalues, ticknames)#, rotation='vertical')
    plt.xlim([0,90])
    plt.ylim([0,ymax])

    # -------------------------- LAEs --------------------------
    ax = fig.add_subplot(2, 1, 2)

    for uvline in linecolors.keys():
       if (uvline != 'MgII') & (uvline != 'CIIIorCIV'): # not showing MgII and CIIIorCIV panels
            yvalues  = laearr[arrayrows[uvline],:]
            linecol  = linecolors[uvline]
            plt.plot(xvalues,yvalues,'o',markersize=marksize,alpha=1.0,color=linecol,lw=lthick,zorder=10,label=linenames[uvline][0])
            plt.plot(xvalues[1:],yvalues[1:],'-',alpha=1.0,color=linecol,lw=lthick,zorder=10)

    ymax = 9
    plt.plot([20,20],[0,ymax],'-',color='black',alpha=1.0,zorder=1,lw=1)

    plt.xticks(xvalues, ticknames)#, rotation='vertical')
    plt.xlim([0,90])
    plt.ylim([0,ymax])

    #--------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='lower center',prop={'size':Fsize},ncol=3,numpoints=1,
                     bbox_to_anchor=(0.5, -0.9),)  # add the legend
    leg.get_frame().set_alpha(0.7)
    #--------------------------

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Setting up and generating plot of histograms in same window')
    xoffset  = 1.0
    errtype  = 'cp' #'binomvar'
    plotname = plotname.replace('.pdf','_sharedaxis_'+errtype+'.pdf')
    fig = plt.figure(figsize=(3, 7))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.2, right=0.97, bottom=0.15, top=0.99)

    Fsize    = 10
    lthick   = 1.0
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(inforstr[:-2],fontsize=Fsize)

    axfull = fig.add_subplot(111)
    axfull.set_xlabel('')
    axfull.set_ylabel(ylabel)
    # Turn off axis lines and ticks of the full window subplot
    axfull.spines['top'].set_color('none')
    axfull.spines['bottom'].set_color('none')
    axfull.spines['left'].set_color('none')
    axfull.spines['right'].set_color('none')
    axfull.tick_params(labelcolor='white', top=False, bottom=False, left=False, right=False)

    # - - - - - - - - - - - - - -
    for ee, emline in enumerate(linecolors.keys()):
        if (emline != 'MgII') & (emline != 'CIIIorCIV'): # not showing MgII and CIIIorCIV panels
            if ee == 0:
                ax       = fig.add_subplot((len(linecolors.keys())-1)*100+11+ee)
            else:
                ax       = fig.add_subplot((len(linecolors.keys())-1)*100+11+ee, sharex=ax, sharey=None)

            linecol  = linecolors[emline]

            #- - - - - - - -
            yvalues  = laearr[arrayrows[emline],:]
            Nobjtot  = laearr_Nobj[arrayrows[emline],:]
            Nsuccess = np.round(laearr_Nobj[arrayrows[emline],:]* yvalues / 100.)
            if errtype == 'cp':
                yval_lo, yval_hi = kbs.get_clopper_pearson_confidence_interval(Nsuccess,Nobjtot,CI=95)
                yerrvals         = np.vstack((yvalues-yval_lo*100., yval_hi*100.-yvalues))
            elif errtype == 'binomvar':
                mean, var, skew, kurt = scipy.stats.binom.stats(Nobjtot, yvalues/100., moments='mvsk')
                yerrvals              = np.vstack((np.sqrt(var)/Nobjtot*100.,
                                                   np.sqrt(var)/Nobjtot*100.))
            else:
                yerrvals = np.zeros([2,4])
            ymaxval_lae = yerrvals[1]+yvalues
            ax.errorbar(xvalues+xoffset,yvalues,xerr=None,yerr=yerrvals,
                        marker='s',lw=lthick, markersize=marksize,alpha=1.0, color='None',
                        markerfacecolor=linecol,ecolor=linecol,markeredgecolor=linecol,zorder=5)

            # ax.plot(xvalues,yvalues,'s',markersize=marksize,alpha=1.0,color=linecol,lw=lthick,zorder=10)
            ax.plot(xvalues[1:]+xoffset,yvalues[1:],':',alpha=1.0,color=linecol,lw=lthick,zorder=10)
            #- - - - - - - -
            yvalues  = allobjarr[arrayrows[emline],:]
            Nobjtot  = allobjarr_Nobj[arrayrows[emline],:]
            Nsuccess = np.round(allobjarr_Nobj[arrayrows[emline],:]* yvalues / 100.)
            if errtype == 'cp':
                yval_lo, yval_hi = kbs.get_clopper_pearson_confidence_interval(Nsuccess,Nobjtot,CI=95)
                yerrvals         = np.vstack((yvalues-yval_lo*100., yval_hi*100.-yvalues))
            elif errtype == 'binomvar':
                mean, var, skew, kurt = scipy.stats.binom.stats(Nobjtot, yvalues/100., moments='mvsk')
                yerrvals              = np.vstack((np.sqrt(var)/Nobjtot*100.,
                                                   np.sqrt(var)/Nobjtot*100.))
            else:
                yerrvals = np.zeros([2,4])
            ymaxval_all = yerrvals[1]+yvalues
            ax.errorbar(xvalues-xoffset,yvalues,xerr=None,yerr=yerrvals,
                        marker='o',lw=lthick, markersize=marksize,alpha=1.0, color='None',
                        markerfacecolor=linecol,ecolor=linecol,markeredgecolor=linecol,zorder=5)

            # ax.plot(xvalues-xoffset,yvalues,'o',markersize=marksize,alpha=1.0,color=linecol,lw=lthick,zorder=10)
            ax.plot(xvalues[1:]-xoffset,yvalues[1:],'-',alpha=1.0,color=linecol,lw=lthick,zorder=10)
            #- - - - - - - -

            ax.set_xticks(xvalues)
            # allyvals = np.append(laearr[arrayrows[emline],:],allobjarr[arrayrows[emline],:])
            allyvals = np.append(ymaxval_lae,ymaxval_all)
            ymax     = np.max(allyvals)
            yrange   = [0.0-ymax*0.1,ymax+ymax*0.05]

            ax.set_ylim(xrange)
            ax.set_ylim(yrange)
            ax.plot([0,50],[0,0],'--',color='black',alpha=1.0,zorder=1,lw=1)
            ax.plot([10,10],yrange,'-',color='black',alpha=1.0,zorder=1,lw=1)
            ax.text(14,yrange[0]+(yrange[1]-yrange[0])*0.8,linenames[emline][0],fontsize=Fsize,ha='left',color=linecol)

    ax.set_xticklabels(ticknames, rotation=25, ha='right')
    # - - - - - - - - - - - - - -
    ax.plot([-10],[-10],'o',markersize=marksize,alpha=1.0,color='black',lw=lthick,zorder=10,label='Full $z$-range')
    ax.plot([-10],[-10],'s',markersize=marksize,alpha=1.0,color='black',lw=lthick,zorder=10,label='LAEs ($z>2.9$)')
    leg = plt.legend(fancybox=True, loc='lower center',prop={'size':Fsize},ncol=2,numpoints=1,
                     bbox_to_anchor=(0.4, -1.4),)  # add the legend
    leg.get_frame().set_alpha(0.7)
    # - - - - - - - - - - - - - -

    for ax in fig.get_axes():
        ax.label_outer()

    # plt.xlim([0,85])
    plt.xlim(xrange)

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def estimateGasPhaseAbundanceFromBylerFittingFunctions(linefluxfile,verbose=True):
    """
    Function estimating the gas-phase abundances (12 + log([O/H])) using the fitting formulas
    described in Byler et al. (2020) eqautions 7 and 8.

    --- Example of use ---
    import uvEmissionlineSearch as uves
    kbswork      = '/Users/kschmidt/work/'

    linefluxfile = kbswork+'catalogs/literaturecollection_emissionlinestrengths/literaturecollection_emissionlinestrengths.fits'
    id_Si3O3C3, Z_Si3O3C3, Z_Si3O3C3_err, id_He2O3C3, Z_He2O3C3, Z_He2O3C3_err  = uves.estimateGasPhaseAbundanceFromBylerFittingFunctions(linefluxfile,verbose=True)

    linefluxfile = kbswork+'MUSE/uvEmissionlineSearch/back2backAnalysis_200213/results_master_catalog_version200213.fits'
    id_Si3O3C3, Z_Si3O3C3, Z_Si3O3C3_err, id_He2O3C3, Z_He2O3C3, Z_He2O3C3_err  = uves.estimateGasPhaseAbundanceFromBylerFittingFunctions(linefluxfile,verbose=True)

    """
    fdat = afits.open(linefluxfile)[1].data

    goodent_OCC_det = np.where(np.isfinite(fdat['f_OIII2'])  & (np.abs(fdat['ferr_OIII2'])  != 99)  &
                               np.isfinite(fdat['f_CIII1'])  & (np.abs(fdat['ferr_CIII'])   != 99))[0]

    goodent_SCC_det = np.where(np.isfinite(fdat['f_SiIII1']) & (np.abs(fdat['ferr_SiIII1']) != 99)  &
                               np.isfinite(fdat['f_CIII1'])  & (np.abs(fdat['ferr_CIII'])   != 99))[0]

    goodent_HCC_det = np.where(np.isfinite(fdat['f_HeII'])   & (np.abs(fdat['ferr_HeII'])   != 99)  &
                               np.isfinite(fdat['f_CIII1'])  & (np.abs(fdat['ferr_CIII'])   != 99))[0]

    # - - - - - - - - - - - - 12 + log([O/H]) for Si3-O3-C3; Byler+20 Eq. 7 - - - - - - - - - - - -
    ent_Z_Si3O3C3  = np.intersect1d(goodent_OCC_det,goodent_SCC_det)
    if len(ent_Z_Si3O3C3) > 0:
        print(' - Calculating '+str(len(ent_Z_Si3O3C3))+' estimates of 12 + log(O/H) for Si3-O3-C3')
        id_Si3O3C3 = fdat['id'][ent_Z_Si3O3C3]
        fOIII2     = fdat['f_OIII2'][ent_Z_Si3O3C3]
        fCIII      = fdat['f_CIII'][ent_Z_Si3O3C3]
        fSiIII1    = fdat['f_SiIII1'][ent_Z_Si3O3C3]

        ferrOIII2     = fdat['ferr_OIII2'][ent_Z_Si3O3C3]
        ferrCIII      = fdat['ferr_CIII'][ent_Z_Si3O3C3]
        ferrSiIII1    = fdat['ferr_SiIII1'][ent_Z_Si3O3C3]

        OCC = np.log10( fOIII2  / fCIII )
        SCC = np.log10( fSiIII1 / fCIII )

        Z_Si3O3C3 = 3.09 + \
                    0.09  * OCC - 1.71  * OCC**2.0 - 0.73 * OCC**3.0 - \
                    16.51 * SCC - 19.84 * SCC**2.0 - 6.26 * SCC**3.0 + \
                    4.79  * OCC * SCC - 0.28 * OCC * SCC**2.0 + 1.67 * OCC**2.0 * SCC

        # Z_Si3O3C3_err = np.asarray(len(OCC)*[0.0])
        # -- error propagation see http://lectureonline.cl.msu.edu/~mmp/labs/error/e2.htm --
        # -- dR(x,y) = sqrt( (dRdx * x_err)**2 + (dRdy * y_err)**2 )
        FRerrOC  = np.abs(fOIII2/fCIII)  * np.sqrt( (ferrOIII2/fOIII2)**2   + (ferrCIII/fCIII)**2 )
        FRerrSC  = np.abs(fSiIII1/fCIII) * np.sqrt( (ferrSiIII1/fSiIII1)**2 + (ferrCIII/fCIII)**2 )

        OCCerr = np.abs( FRerrOC / (fOIII2/fCIII))   / np.log(10)
        SCCerr = np.abs( FRerrSC / (fSiIII1/fCIII))  / np.log(10)

        dRdOCC   =  0.09  - 2.0*1.71*OCC  - 3.0*0.73*OCC**2.0 + 4.79*SCC - 0.28*SCC**2.0    + 2.0*1.76*OCC*SCC
        dRdSCC   = -16.51 - 2.0*19.84*SCC - 3.0*6.26*SCC**2.0 + 4.79*OCC - 2.0*0.28*OCC*SCC + 1.67*OCC**2.0

        Z_Si3O3C3_err = np.sqrt( (dRdOCC * OCCerr)**2.0 + (dRdSCC * SCCerr)**2.0 )

        # Z_Si3O3C3_err = np.sqrt( (0.09  * OCCerr)**2.0 + (1.71  * OCCerr*OCC * 2.0)**2.0 + (0.73 * OCCerr*OCC * 3.0)**2.0 + \
        #                          (16.51 * SCCerr)**2.0 + (19.84 * SCCerr*SCC * 2.0)**2.0 + (6.26 * SCCerr*SCC * 3.0)**2.0 + \
        #                          (4.79  * OCC * SCC      * np.sqrt(((OCCerr/OCC)**2.0 + (SCCerr/SCC)**2.0 )))**2.0 +
        #                          (0.28  * OCC * SCC**2.0 * np.sqrt(((OCCerr/OCC)**2.0 + (SCCerr/SCC)**2.0 * 2.0 )))**2.0 +
        #                          (1.67  * OCC**2.0 * SCC * np.sqrt(((OCCerr/OCC)**2.0 * 2.0 + (SCCerr/SCC)**2.0 )))**2.0   )

    else:
        print(' - No estimates of 12 + log(O/H) for Si3-O3-C3 to return ')
        id_Si3O3C3    = np.nan
        Z_Si3O3C3     = np.nan
        Z_Si3O3C3_err = np.nan

    # - - - - - - - - - - - - 12 + log([O/H]) for He2-O3-C3; Byler+20 Eq. 8 - - - - - - - - - - - -
    ent_Z_He2O3C3 = np.intersect1d(goodent_OCC_det,goodent_HCC_det)
    if len(ent_Z_He2O3C3) > 0:
        print(' - Calculating '+str(len(ent_Z_He2O3C3))+' estimates of 12 + log([O/H]) for He2-O3-C3')
        id_He2O3C3 = fdat['id'][ent_Z_He2O3C3]
        fOIII2     = fdat['f_OIII2'][ent_Z_He2O3C3]
        fCIII      = fdat['f_CIII'][ent_Z_He2O3C3]
        fHeII    = fdat['f_HeII'][ent_Z_He2O3C3]

        ferrOIII2     = fdat['ferr_OIII2'][ent_Z_He2O3C3]
        ferrCIII      = fdat['ferr_CIII'][ent_Z_He2O3C3]
        ferrHeII      = fdat['ferr_HeII'][ent_Z_He2O3C3]

        OCC = np.log10( fOIII2  / fCIII )
        HCC = np.log10( fHeII   / fCIII )

        Z_He2O3C3 = 6.88 - \
                    1.13  * OCC - 0.46  * OCC**2.0 - 0.03 * OCC**3.0 - \
                    0.61  * HCC - 0.02  * HCC**2.0 - 0.04 * HCC**3.0 - \
                    0.32  * OCC * HCC + 0.03 * OCC * HCC**2.0 - 0.21 * OCC**2.0 * HCC

        Z_He2O3C3_err = np.asarray(len(OCC)*[0.0])
        # -- error propagation see http://lectureonline.cl.msu.edu/~mmp/labs/error/e2.htm --
        # -- dR(x,y) = sqrt( (dRdx * x_err)**2 + (dRdy * y_err)**2 )
        FRerrOC  = np.abs(fOIII2/fCIII)  * np.sqrt( (ferrOIII2/fOIII2)**2   + (ferrCIII/fCIII)**2 )
        FRerrHC  = np.abs(fHeII/fCIII)   * np.sqrt( (ferrHeII/fHeII)**2 + (ferrCIII/fCIII)**2 )

        OCCerr = np.abs( FRerrOC / (fOIII2/fCIII * np.log(10)) )
        HCCerr = np.abs( FRerrHC / (fHeII/fCIII  * np.log(10))  )

        dRdOCC   =  1.13 - 2.0*0.46*OCC  - 3.0*0.03*OCC**2.0 + 0.32*HCC + 0.03*HCC**2.0    - 2.0*0.21*OCC*HCC
        dRdHCC   =  -0.61 - 2.0*0.02*HCC - 3.0*0.04*HCC**2.0 + 0.32*OCC + 2.0*0.28*OCC*HCC - 0.21*OCC**2.0

        Z_He2O3C3_err = np.sqrt( (dRdOCC * OCCerr)**2.0 + (dRdHCC * HCCerr)**2.0 )

        # Z_He2O3C3_err = np.sqrt( (1.13  * OCCerr)**2.0 + (0.46  * OCCerr/OCC *2.0)**2.0 + (0.03 * OCCerr/OCC * 3.0)**2.0 + \
        #                          (0.61  * HCCerr)**2.0 + (0.02  * HCCerr/HCC *2.0)**2.0 + (0.04 * HCCerr/HCC * 3.0)**2.0 + \
        #                          (0.32  * OCC * HCC       * np.sqrt(((OCCerr/OCC)**2.0 + (HCCerr/HCC)**2.0 )))**2.0 +
        #                          (0.03  * OCC * HCC**2.0  * np.sqrt(((OCCerr/OCC)**2.0 + (HCCerr/HCC)**2.0 * 2.0 )))**2.0 +
        #                          (0.21  * OCC**2.0 * HCC  * np.sqrt(((OCCerr/OCC)**2.0 * 2.0 + (HCCerr/HCC)**2.0 )))**2.0   )
    else:
        print(' - No estimates of 12 + log(O/H) for He2-O3-C3 to return ')
        id_He2O3C3    = np.nan
        Z_He2O3C3     = np.nan
        Z_He2O3C3_err = np.nan

    # pdb.set_trace()
    return id_Si3O3C3, Z_Si3O3C3, Z_Si3O3C3_err, id_He2O3C3, Z_He2O3C3, Z_He2O3C3_err


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_GasPhaseAbundances(masterfits,outputdir,withliterature=True,overwrite=False,litsymboldot=False,verbose=True):
    """
    Plotting gas phase abundances estimate with uves.estimateGasPhaseAbundanceFromBylerFittingFunctions()

    --- Example of use ---
    import uvEmissionlineSearch as uves
    outputdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/BylerAbundanceFigures/'
    kbswork      = '/Users/kschmidt/work/'
    kbswork+'MUSE/uvEmissionlineSearch/back2backAnalysis_200213/results_master_catalog_version200213.fits'
    uves.plot_GasPhaseAbundances(masterfits,outputdir,overwrite=True,withliterature=False,verbose=True)
    uves.plot_GasPhaseAbundances(masterfits,overwrite=True,withliterature=True,verbose=True)
    """
    if withliterature:
         MUSEsymbolblackedge=True

    # - - - - - - get entries for UVES results - - - - - -
    if verbose: print('\n ----- Estimating abundances ----- ')
    if verbose: print(' - Estimate abundances for UVES detections ')
    fdat_uves         = afits.open(masterfits)[1].data

    id_Si3O3C3_uves, Z_Si3O3C3_uves, Z_Si3O3C3_err_uves, id_He2O3C3_uves, Z_He2O3C3_uves, Z_He2O3C3_err_uves  = \
        uves.estimateGasPhaseAbundanceFromBylerFittingFunctions(masterfits,verbose=True)

    bothest_id_uves             = np.intersect1d(id_Si3O3C3_uves,id_He2O3C3_uves)
    bothest_ent_uves            = []
    bothest_ent_Si3O3C3_uves    = []
    bothest_ent_He2O3C3_uves    = []
    ent_Si3O3C3_uves    = []
    ent_He2O3C3_uves    = []
    for ent, objid in enumerate(fdat_uves['id']):
        if objid in bothest_id_uves:
            bothest_ent_uves.append(ent)
        if objid in id_Si3O3C3_uves:
            ent_Si3O3C3_uves.append(ent)
            if objid in bothest_id_uves:
                bothest_ent_Si3O3C3_uves.append(ent)
        if objid in id_He2O3C3_uves:
            ent_He2O3C3_uves.append(ent)
            if objid in bothest_id_uves:
                bothest_ent_He2O3C3_uves.append(ent)

    bothest_ent_Si3O3C3_uves_Z  = []
    for ent, objid in enumerate(id_Si3O3C3_uves):
        if objid in bothest_id_uves:
            bothest_ent_Si3O3C3_uves_Z.append(ent)

    bothest_ent_He2O3C3_uves_Z  = []
    for ent, objid in enumerate(id_He2O3C3_uves):
        if objid in bothest_id_uves:
            bothest_ent_He2O3C3_uves_Z.append(ent)

    # - - - - - - get entries for literature results - - - - - -
    if verbose: print(' - Estimate abundances for literature collection ')
    kbswork           = '/Users/kschmidt/work/'
    linefluxfile_lit  = kbswork+'catalogs/literaturecollection_emissionlinestrengths/literaturecollection_emissionlinestrengths.fits'
    fdat_lit          = afits.open(linefluxfile_lit)[1].data
    id_Si3O3C3_lit, Z_Si3O3C3_lit, Z_Si3O3C3_err_lit, id_He2O3C3_lit, Z_He2O3C3_lit, Z_He2O3C3_err_lit  =\
        uves.estimateGasPhaseAbundanceFromBylerFittingFunctions(linefluxfile_lit,verbose=True)

    bothest_id_lit             = np.intersect1d(id_Si3O3C3_lit,id_He2O3C3_lit)
    bothest_ent_lit            = []
    bothest_ent_Si3O3C3_lit    = []
    bothest_ent_He2O3C3_lit    = []
    ent_Si3O3C3_lit    = []
    ent_He2O3C3_lit    = []
    for ent, objid in enumerate(fdat_lit['id']):
        if objid in bothest_id_lit:
            bothest_ent_lit.append(ent)
        if objid in id_Si3O3C3_lit:
            ent_Si3O3C3_lit.append(ent)
            if objid in bothest_id_lit:
                bothest_ent_Si3O3C3_lit.append(ent)
        if objid in id_He2O3C3_lit:
            ent_He2O3C3_lit.append(ent)
            if objid in bothest_id_lit:
                bothest_ent_He2O3C3_lit.append(ent)

    bothest_ent_Si3O3C3_lit_Z  = []
    for ent, objid in enumerate(id_Si3O3C3_lit):
        if objid in bothest_id_lit:
            bothest_ent_Si3O3C3_lit_Z.append(ent)

    bothest_ent_He2O3C3_lit_Z  = []
    for ent, objid in enumerate(id_He2O3C3_lit):
        if objid in bothest_id_lit:
            bothest_ent_He2O3C3_lit_Z.append(ent)

    # - - - - - - - - - - - - - - - - Generate Plots of abundances - - - - - - - - - - - - - - - -
    if verbose: print('\n ----- Generating plots ----- ')
    # - - - - - - Setting redshift log axes and Zrange - - - - - -
    xlog = False
    if xlog:
        xrange = [0.01,10.0]
    else:
        xrange = [-0.3,5.2]
    Zrange = [6.5,8.9]

    # - - - - - - Si3O3C3 vs He2O3C3     - - - - - -
    plotname   = outputdir+'Byler_abundanceplot_Si3O3C3vsHe2O3C3_uves.pdf'
    xvalues    = Z_Si3O3C3_uves[bothest_ent_Si3O3C3_uves_Z]
    yvalues    = Z_He2O3C3_uves[bothest_ent_He2O3C3_uves_Z]
    xerr       = Z_Si3O3C3_err_uves[bothest_ent_Si3O3C3_uves_Z]
    yerr       = Z_He2O3C3_err_uves[bothest_ent_He2O3C3_uves_Z]
    IDsALL     = bothest_id_uves
    cdatvec    = fdat_uves['redshift'][bothest_ent_uves]
    # cdatvec    = fdat_uves['EW0_CIII'][bothest_ent_uves]
    point_text = None #bothest_id_uves.astype(str)

    if withliterature:
        plotname   = plotname.replace('.pdf','AndLiterature.pdf')
        xvalues    = np.append(xvalues,    Z_Si3O3C3_lit[bothest_ent_Si3O3C3_lit_Z])
        yvalues    = np.append(yvalues,    Z_He2O3C3_lit[bothest_ent_He2O3C3_lit_Z])
        xerr       = np.append(xerr,       Z_Si3O3C3_err_lit[bothest_ent_Si3O3C3_lit_Z])
        yerr       = np.append(yerr,       Z_He2O3C3_err_lit[bothest_ent_He2O3C3_lit_Z])
        if litsymboldot:
            IDsALL     = np.append(IDsALL,bothest_id_lit*0.0+990000000000)
        else:
            IDsALL     = np.append(IDsALL,bothest_id_lit)
        cdatvec    = np.append(cdatvec,    fdat_lit['redshift'][bothest_ent_lit])
        # cdatvec    = np.append(cdatvec,    fdat_lit['EW0_CIII'][bothest_ent_lit])
        point_text = None #np.append(point_text, bothest_id_lit.astype(str))

    xlabel     = '12 + log$_{10}$(O/H) \n\Large{from SiIII1883, OIII1666 and CIII1908}'
    ylabel     = '12 + log$_{10}$(O/H) \n\Large{from HeII1640, OIII1666 and CIII1908}'
    colortype  = 'redshift'
    # colortype  = 'EW0_CIII'
    colorcode  = True

    Nbins   = 12 #np.ceil(np.sqrt(len(xvalues)))
    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,
                                                   'dummydat',linetype='onetooneWZsun',title=None,ids=IDsALL,
                                                   ylog=False,xlog=False,yrange=Zrange,xrange=Zrange,
                                                   colortype=colortype,colorcode=colorcode,cdatvec=cdatvec,
                                                   point_text=point_text,photoionizationplotparam=None,
                                                   histaxes=True,Nbins=Nbins, MUSEsymbolblackedge=MUSEsymbolblackedge,
                                                   overwrite=overwrite,verbose=verbose)

    # - - - - - - Si3O3C3 vs redshift    - - - - - -
    plotname   = outputdir+'Byler_abundanceplot_zvsSi3O3C3_uves.pdf'
    xvalues    = fdat_uves['redshift'][ent_Si3O3C3_uves]
    yvalues    = Z_Si3O3C3_uves
    xerr       = [None]*len(Z_Si3O3C3_uves)
    yerr       = Z_Si3O3C3_err_uves
    IDsALL     = id_Si3O3C3_uves
    cdatvec    = fdat_uves['EW0_CIII'][ent_Si3O3C3_uves]
    point_text = None

    if withliterature:
        plotname   = plotname.replace('.pdf','AndLiterature.pdf')
        xvalues    = np.append(xvalues,  fdat_lit['redshift'][ent_Si3O3C3_lit])
        yvalues    = np.append(yvalues,  Z_Si3O3C3_lit)
        xerr       = np.append(xerr,     [None]*len(Z_Si3O3C3_lit))
        yerr       = np.append(yerr,     Z_Si3O3C3_err_lit)
        if litsymboldot:
            IDsALL     = np.append(IDsALL,id_Si3O3C3_lit*0.0+990000000000)
        else:
            IDsALL     = np.append(IDsALL,id_Si3O3C3_lit)
        # cdatvec    = np.append(cdatvec,  fdat_lit['redshift'][ent_Si3O3C3_lit]) #fdat_lit['f_CIII'][ent_Si3O3C3_lit])
        cdatvec    = np.append(cdatvec,  fdat_lit['EW0_CIII'][ent_Si3O3C3_lit]) #fdat_lit['f_CIII'][ent_Si3O3C3_lit])
        point_text = None

    xlabel     = '$z$'
    ylabel     = '12 + log$_{10}$(O/H) \n\Large{from SiIII1883, OIII1666 and CIII1908}'
    # colortype  = 'redshift' #'f(CIII) [1e-20 erg/s/cm$^2$]'
    colortype  = 'EW0_CIII'
    colorcode  = True

    Nbins   = 12 #np.ceil(np.sqrt(len(xvalues)))
    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,
                                                   'dummydat',linetype='yZsun',title=None,ids=IDsALL,
                                                   ylog=False,xlog=xlog,yrange=Zrange,xrange=xrange,
                                                   colortype=colortype,colorcode=colorcode,cdatvec=cdatvec,
                                                   point_text=point_text,photoionizationplotparam=None,
                                                   histaxes=True,Nbins=Nbins, MUSEsymbolblackedge=MUSEsymbolblackedge,
                                                   overwrite=overwrite,verbose=verbose)


    # - - - - - - He2O3C3 vs redshift    - - - - - -
    plotname   = outputdir+'Byler_abundanceplot_zvsHe2O3C3_uves.pdf'
    xvalues    = fdat_uves['redshift'][ent_He2O3C3_uves]
    yvalues    = Z_He2O3C3_uves
    xerr       = [None]*len(Z_He2O3C3_uves)
    yerr       = Z_He2O3C3_err_uves
    IDsALL     = id_He2O3C3_uves
    cdatvec    = fdat_uves['EW0_CIII'][ent_He2O3C3_uves]
    point_text = None

    if withliterature:
        plotname   = plotname.replace('.pdf','AndLiterature.pdf')
        xvalues    = np.append(xvalues,  fdat_lit['redshift'][ent_He2O3C3_lit])
        yvalues    = np.append(yvalues,  Z_He2O3C3_lit)
        xerr       = np.append(xerr,     [None]*len(Z_He2O3C3_lit))
        yerr       = np.append(yerr,     Z_He2O3C3_err_lit)
        if litsymboldot:
            IDsALL     = np.append(IDsALL,id_He2O3C3_lit*0.0+990000000000)
        else:
            IDsALL     = np.append(IDsALL,id_He2O3C3_lit)
        # cdatvec    = np.append(cdatvec,  fdat_lit['redshift'][ent_He2O3C3_lit]) #fdat_lit['f_CIII'][ent_He2O3C3_lit])
        cdatvec    = np.append(cdatvec,  fdat_lit['EW0_CIII'][ent_He2O3C3_lit]) #fdat_lit['f_CIII'][ent_He2O3C3_lit])
        point_text = None

    xlabel     = '$z$'
    ylabel     = '12 + log$_{10}$(O/H) \n\Large{from HeII1640, OIII1666 and CIII1908}'
    # colortype  = 'redshift' # 'f(CIII) [1e-20 erg/s/cm$^2$]'
    colortype  = 'EW0_CIII'

    colorcode  = True

    Nbins   = 12 #np.ceil(np.sqrt(len(xvalues)))
    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,
                                                   'dummydat',linetype='yZsun',title=None,ids=IDsALL,
                                                   ylog=False,xlog=xlog,yrange=Zrange,xrange=xrange,
                                                   colortype=colortype,colorcode=colorcode,cdatvec=cdatvec,
                                                   point_text=point_text,photoionizationplotparam=None,
                                                   histaxes=True,Nbins=Nbins, MUSEsymbolblackedge=MUSEsymbolblackedge,
                                                   overwrite=overwrite,verbose=verbose)

    # - - - - - - Si3O3C3 vs EW(Lya)     - - - - - -
    # - - - - - - Si3O3C3 vs EW(CIII)    - - - - - -
    # - - - - - - He2O3C3 vs EW(Lya)     - - - - - -
    # - - - - - - He2O3C3 vs EW(CIII)    - - - - - -

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def perform_PyNeb_calc_main(linefluxcatalog,outputfile='./pyneb_calculations_results.txt',Nsigma=1.0,
                            skipOIIIandCIV = True,
                            curveresolution=0.1,generateExtraPlots=False,overwrite=False,verbose=True):
    """
    Main function to handle caluclations and diagostics using PyNeb

    --- Example of use ---
    import uvEmissionlineSearch as uves
    uvesdir          = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'
    outfile          = uvesdir+'PyNebCalculations/pyneb_calculations_results.txt'
    linefluxcatalog  = uvesdir+'back2backAnalysis_200213/results_master_catalog_version200213.fits'
    uves.perform_PyNeb_calc_main(linefluxcatalog,outputfile=outfile,generateExtraPlots=True,overwrite=True,verbose=True,curveresolution=0.001)

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Setting up output:\n   '+outputfile)
    if os.path.isfile(outputfile) and (overwrite == False):
        sys.exit('The output file '+outputfile+' already exists and overwrite=False ')

    fout = open(outputfile,'w')
    fout.write('# This file contains the output from uves.perform_PyNeb_calc_main() generated on '+kbs.DandTstr2()+' \n')
    fout.write('# based on the line flux catalog '+linefluxcatalog+' \n')
    fout.write('# ')
    fout.write('# id n_e n_e_min n_e_max T_e estimator \n') # n_crit
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading data in line flux catalog:\n   '+linefluxcatalog)
    fluxdat = afits.open(linefluxcatalog)[1].data
    fluxdat = fluxdat[fluxdat['duplicationID'] == 0]


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('--- Estimating n_e and T_e from set of flux ratios (running getCrossTemDen) ---')
    FR_Si3 = 'FR_SiIII1SiIII2'
    if verbose: print(' - using flux ratio '+FR_Si3)
    FRval_Si3 = fluxdat[FR_Si3]
    FRerr_Si3 = fluxdat[FR_Si3.replace('FR','FRerr')]*Nsigma

    FR_C3    = 'FR_CIII1CIII2'
    if verbose: print(' - using flux ratio '+FR_C3)
    FRval_C3 = fluxdat[FR_C3]
    FRerr_C3 = fluxdat[FR_C3.replace('FR','FRerr')]*Nsigma

    Si3 = pn.Atom('Si', 3)
    C3  = pn.Atom('C', 3)

    diags = pn.Diagnostics(addAll=True)
    diags.addDiag(label='SiIII 1883/1892', diag_tuple=('Si3', 'L(1883)/L(1892)', 'RMS([E(1883),E(1892)])') )

    # for diag in np.sort(pn.diags_dict.keys()):
    #     print('"{0}" : {1}'.format(diag, pn.diags_dict[diag]))

    goodent_Si3  = np.where(np.isfinite(FRval_Si3) & (np.abs(FRerr_Si3) != 99))[0]
    goodent_C3   = np.where(np.isfinite(FRval_C3) & (np.abs(FRerr_C3) != 99))[0]
    goodentC3Si3 = np.intersect1d(goodent_C3,goodent_Si3).astype(int)

    objids   = fluxdat['id'][goodentC3Si3]
    FR_tem   = FRval_Si3[goodentC3Si3]  # Temperature sensitive
    FR_den   = FRval_C3[goodentC3Si3]   # Densitiy sensitive

    for oo, objid in enumerate(objids):
        Te, Ne = diags.getCrossTemDen('SiIII 1883/1892', '[CIII] 1909/1907', FRval_Si3[goodentC3Si3][oo], 1./FRval_C3[goodentC3Si3][oo])
                                      # guess_tem=1e4, tol_tem = 10., tol_den = 10., max_iter = 5,
                                      # start_tem=5000, end_tem=20000, start_den=1e2, end_den=1e6)
        if verbose: print('   '+str(objid)+' (Te,Ne) = '+str((Te,Ne))+
                          '  for tem(SiIII1/SiIII2)='+str(FRval_Si3[goodentC3Si3][oo])+
                          ' and den(CIII1/CIII2)='+str(FRval_C3[goodentC3Si3][oo]))

    for oo, objid in enumerate(objids):
        Te, Ne = diags.getCrossTemDen('[CIII] 1909/1907', 'SiIII 1883/1892', 1./FRval_C3[goodentC3Si3][oo], FRval_Si3[goodentC3Si3][oo])
        #Te, Ne = diags.getCrossTemDen("[OIII] 4363/5007+", '[CIII] 1909/1907', 0.001,0.8)
        if verbose: print('   '+str(objid)+' (Te,Ne) = '+str((Te,Ne))+
                          '  for tem(CIII1/CIII2)='+str(FRval_C3[goodentC3Si3][oo])+
                          ' and den(SiIII1/SiIII2)='+str(FRval_Si3[goodentC3Si3][oo]))
    # FR1   = 'FR_CIII1CIII2'
    # FR2   = 'FR_OIII1OIII2'
    # if verbose: print(' - using flux ratios '+FR1+' and '+FR2)
    # Te, Ne = diags.getCrossTemDen('[NII] 5755/6548', '[SII] 6731/6716', 0.02, 1.0)



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('--- Estimating n_e from a single fluxratio (running getTemDen) ---')
    T_e_fix_vals = {'5k':5.e3, '10k':1.e4, '20k':2.e4}
    yvals_curve  = np.arange(0.0,2.0,curveresolution)

    if not skipOIIIandCIV:
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        FR    = 'FR_OIII1OIII2'
        FRval = fluxdat[FR]
        FRerr = fluxdat[FR.replace('FR','FRerr')]*Nsigma
        if verbose: print(' - using flux ratio '+FR)
        pn.atomicData.setDataFile('o_iii_coll_AK99.dat')
        O3 = pn.Atom('O', 3)
        # O3.printIonic(tem=10000., den=1e3, printA=True, printPop=True, printCrit=True)
        if verbose: print('   '+str(O3))

        #------------------------------------------------------------------------------
        if generateExtraPlots:
            if verbose: print(' - Setting up and generating Grotrian diagram ')
            plotname = outputfile.replace('.txt','_OIII_Grotrian.pdf')
            fig = plt.figure(figsize=(7, 5))
            fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.97, bottom=0.40, top=0.97)
            Fsize    = 10
            lthick   = 2
            marksize = 6
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif',size=Fsize)
            plt.rc('xtick', labelsize=Fsize)
            plt.rc('ytick', labelsize=Fsize)
            plt.clf()
            plt.ioff()
            #plt.title(inforstr[:-2],fontsize=Fsize)

            O3.plotGrotrian(tem=1e4, den=1e2, thresh_int=1e-3, unit = 'eV', detailed=False)

            if verbose: print('   Saving plot to '+plotname)
            plt.savefig(plotname)
            plt.clf()
            plt.close('all')
        #------------------------------------------------------------------------------

        n_e_Te3 = O3.getTemDen(yvals_curve, tem=T_e_fix_vals['5k'],  wave1=1661, wave2=1666)
        n_e_Te4 = O3.getTemDen(yvals_curve, tem=T_e_fix_vals['10k'], wave1=1661, wave2=1666)
        n_e_Te5 = O3.getTemDen(yvals_curve, tem=T_e_fix_vals['20k'], wave1=1661, wave2=1666)

        for TeKey in T_e_fix_vals:
            T_e_fix = T_e_fix_vals[TeKey]
            goodent_O3 = np.where(np.isfinite(FRval) & (np.abs(FRerr) != 99))[0]
            if len(goodent_O3) == 0:
                if verbose: print('   No good measurements found in line flux catalog for '+FR)
            else:
                if verbose: print('   Found '+str(len(goodent_O3))+' measurements in line flux catalog for '+FR)

                n_e_O3     = O3.getTemDen(FRval[goodent_O3], tem=T_e_fix, wave1=1661, wave2=1666)
                n_e_min_O3 = O3.getTemDen(FRval[goodent_O3]+Nsigma*FRerr[goodent_O3], tem=T_e_fix, wave1=1661, wave2=1666)
                n_e_max_O3 = O3.getTemDen(FRval[goodent_O3]-Nsigma*FRerr[goodent_O3], tem=T_e_fix, wave1=1661, wave2=1666)
                ylabel     = 'OIII1661/OIII1666'
                plotname   = outputfile.replace('.txt','_OIII_ne_estimates_Te'+TeKey+'.pdf')

                uves.plot_neForFR(plotname,fout,T_e_fix,FR,ylabel,FRval,FRerr,n_e_O3,n_e_min_O3,n_e_max_O3,
                                  fluxdat,goodent_O3,yvals_curve,n_e_Te3,n_e_Te4,n_e_Te5,verbose=True)


        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        FR    = 'FR_CIV1CIV2'
        FRval = fluxdat[FR]
        FRerr = fluxdat[FR.replace('FR','FRerr')]*Nsigma
        if verbose: print(' - using flux ratio '+FR)
        C4 = pn.Atom('C', 4)
        if verbose: print('   '+str(C4))

        #------------------------------------------------------------------------------
        if generateExtraPlots:
            if verbose: print(' - Setting up and generating Grotrian diagram ')
            plotname = outputfile.replace('.txt','_CIV_Grotrian.pdf')
            fig = plt.figure(figsize=(7, 5))
            fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.97, bottom=0.40, top=0.97)
            Fsize    = 10
            lthick   = 2
            marksize = 6
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif',size=Fsize)
            plt.rc('xtick', labelsize=Fsize)
            plt.rc('ytick', labelsize=Fsize)
            plt.clf()
            plt.ioff()
            #plt.title(inforstr[:-2],fontsize=Fsize)

            C4.plotGrotrian(tem=1e4, den=1e2, thresh_int=1e-3, unit = 'eV', detailed=False)

            if verbose: print('   Saving plot to '+plotname)
            plt.savefig(plotname)
            plt.clf()
            plt.close('all')
        #------------------------------------------------------------------------------
        n_e_Te3 = C4.getTemDen(yvals_curve, tem=T_e_fix_vals['5k'],  wave1=1548, wave2=1551)
        n_e_Te4 = C4.getTemDen(yvals_curve, tem=T_e_fix_vals['10k'], wave1=1548, wave2=1551)
        n_e_Te5 = C4.getTemDen(yvals_curve, tem=T_e_fix_vals['20k'], wave1=1548, wave2=1551)

        for TeKey in T_e_fix_vals:
            T_e_fix = T_e_fix_vals[TeKey]
            goodent_C4 = np.where(np.isfinite(FRval) & (np.abs(FRerr) != 99))[0]
            if len(goodent_C4) == 0:
                if verbose: print('   No good measurements found in line flux catalog for '+FR)
            else:
                if verbose: print('   Found '+str(len(goodent_C4))+' measurements in line flux catalog for '+FR)

                n_e_C4     = C4.getTemDen(FRval[goodent_C4], tem=T_e_fix, wave1=1548, wave2=1551)
                n_e_min_C4 = C4.getTemDen(FRval[goodent_C4]+Nsigma*FRerr[goodent_C4], tem=T_e_fix, wave1=1548, wave2=1551)
                n_e_max_C4 = C4.getTemDen(FRval[goodent_C4]-Nsigma*FRerr[goodent_C4], tem=T_e_fix, wave1=1548, wave2=1551)
                ylabel   = 'CIV1548/CIV1551'
                plotname = outputfile.replace('.txt','_CIV_ne_estimates_Te'+TeKey+'.pdf')

                uves.plot_neForFR(plotname,fout,T_e_fix,FR,ylabel,FRval,FRerr,n_e_C4,n_e_min_C4,n_e_max_C4,
                                  fluxdat,goodent_C4,yvals_curve,n_e_Te3,n_e_Te4,n_e_Te5,verbose=True)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    FR_Si3 = 'FR_SiIII1SiIII2'
    if verbose: print(' - using flux ratio '+FR_Si3)
    FRval_Si3 = fluxdat[FR_Si3]
    FRerr_Si3 = fluxdat[FR_Si3.replace('FR','FRerr')]*Nsigma
    if verbose: print(' - using flux ratio '+FR_Si3)
    Si3 = pn.Atom('Si', 3)
    if verbose: print('   '+str(Si3))

    #------------------------------------------------------------------------------
    if generateExtraPlots:
        if verbose: print(' - Setting up and generating Grotrian diagram ')
        plotname = outputfile.replace('.txt','_SiIII_Grotrian.pdf')
        fig = plt.figure(figsize=(7, 5))
        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.97, bottom=0.40, top=0.97)
        Fsize    = 10
        lthick   = 2
        marksize = 6
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        #plt.title(inforstr[:-2],fontsize=Fsize)

        Si3.plotGrotrian(tem=1e4, den=1e2, thresh_int=1e-3, unit = 'eV', detailed=False)

        if verbose: print('   Saving plot to '+plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')
    #------------------------------------------------------------------------------
    n_e_Te3_Si3 = Si3.getTemDen(yvals_curve, tem=T_e_fix_vals['5k'],  wave1=1883, wave2=1892)
    n_e_Te4_Si3 = Si3.getTemDen(yvals_curve, tem=T_e_fix_vals['10k'], wave1=1883, wave2=1892)
    n_e_Te5_Si3 = Si3.getTemDen(yvals_curve, tem=T_e_fix_vals['20k'], wave1=1883, wave2=1892)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    FR_C3 = 'FR_CIII1CIII2'
    FRval_C3 = fluxdat[FR_C3]
    FRerr_C3 = fluxdat[FR_C3.replace('FR','FRerr')]*Nsigma
    if verbose: print(' - using flux ratio '+FR_C3)
    C3 = pn.Atom('C', 3)
    if verbose: print('   '+str(C3))
    #------------------------------------------------------------------------------
    if generateExtraPlots:
        if verbose: print(' - Setting up and generating Grotrian diagram ')
        plotname = outputfile.replace('.txt','_CIII_Grotrian.pdf')
        fig = plt.figure(figsize=(7, 5))
        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.97, bottom=0.40, top=0.97)
        Fsize    = 10
        lthick   = 2
        marksize = 6
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        #plt.title(inforstr[:-2],fontsize=Fsize)

        C3.plotGrotrian(tem=1e4, den=1e2, thresh_int=1e-3, unit = 'eV', detailed=False)

        if verbose: print('   Saving plot to '+plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')
    #------------------------------------------------------------------------------
    n_e_Te3_C3 = C3.getTemDen(yvals_curve, tem=T_e_fix_vals['5k'],  wave1=1907, wave2=1909)
    n_e_Te4_C3 = C3.getTemDen(yvals_curve, tem=T_e_fix_vals['10k'], wave1=1907, wave2=1909)
    n_e_Te5_C3 = C3.getTemDen(yvals_curve, tem=T_e_fix_vals['20k'], wave1=1907, wave2=1909)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Calculating and plotting electron densitites')
    for TeKey in T_e_fix_vals:
        T_e_fix = T_e_fix_vals[TeKey]
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # --- SiIII ----
        goodent_Si3 = np.where(np.isfinite(FRval_Si3) & (np.abs(FRerr_Si3) != 99))[0]
        if len(goodent_Si3) == 0:
            if verbose: print('   No good measurements found in line flux catalog for '+FR_Si3)
        else:
            if verbose: print(' - Found '+str(len(goodent_Si3))+' measurements in line flux catalog for '+FR_Si3+' which are:')

            n_e_Si3      = Si3.getTemDen(FRval_Si3[goodent_Si3], tem=T_e_fix, wave1=1883, wave2=1892)
            n_e_min_Si3  = Si3.getTemDen(FRval_Si3[goodent_Si3]+Nsigma*FRerr_Si3[goodent_Si3], tem=T_e_fix, wave1=1883, wave2=1892)
            n_e_max_Si3  = Si3.getTemDen(FRval_Si3[goodent_Si3]-Nsigma*FRerr_Si3[goodent_Si3], tem=T_e_fix, wave1=1883, wave2=1892)

            s2nvalues_Si3 = fluxdat[FR_Si3.replace('FR','s2n').split('1')[0]][goodent_Si3]
            if verbose: print(' - Number of S/N(FELIS) as a double check: '+str(len(s2nvalues_Si3 )))
            if verbose: print(' - S/N(FELIS): (mean.median,min,max,std) = '+
                              str(np.mean(s2nvalues_Si3))+','+
                              str(np.median(s2nvalues_Si3))+','+
                              str(np.min(s2nvalues_Si3))+','+
                              str(np.max(s2nvalues_Si3))+','+
                              str(np.std(s2nvalues_Si3)))

            if verbose:
                for gg, gent in enumerate(goodent_Si3):
                    print('   Te='+str(TeKey)+' goodent='+str(gent)+'='+str(fluxdat['id'][gent])+
                          '  '+str(FRval_Si3[gent])+' +/- '+str(Nsigma*FRerr_Si3[gent])+
                          '   ->  n_e(SiIII) = '+str(n_e_Si3[gg])+'-'+str(n_e_min_Si3[gg])+'+'+str(n_e_max_Si3[gg]))

            ylabel   = 'SiIII1883/SiIII1892'
            plotname = outputfile.replace('.txt','_SiIII_ne_estimates_Te'+TeKey+'.pdf')
            uves.plot_neForFR(plotname,fout,T_e_fix,FR_Si3,ylabel,FRval_Si3,FRerr_Si3,n_e_Si3,n_e_min_Si3,n_e_max_Si3,
                              fluxdat,goodent_Si3,yvals_curve,n_e_Te3_Si3,n_e_Te4_Si3,n_e_Te5_Si3,verbose=True)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # --- CIII ----
        goodent_C3 = np.where(np.isfinite(FRval_C3) & (np.abs(FRerr_C3) != 99))[0]
        if len(goodent_C3) == 0:
            if verbose: print('   No good measurements found in line flux catalog for '+FR_C3)
        else:
            if verbose: print(' - Found '+str(len(goodent_C3))+' measurements in line flux catalog for '+FR_C3+' which are:')

            n_e_C3     = C3.getTemDen(FRval_C3[goodent_C3], tem=T_e_fix, wave1=1907, wave2=1909)
            n_e_min_C3 = C3.getTemDen(FRval_C3[goodent_C3]+Nsigma*FRerr_C3[goodent_C3], tem=T_e_fix, wave1=1907, wave2=1909)
            n_e_max_C3 = C3.getTemDen(FRval_C3[goodent_C3]-Nsigma*FRerr_C3[goodent_C3], tem=T_e_fix, wave1=1907, wave2=1909)

            s2nvalues_C3 = fluxdat[FR_C3.replace('FR','s2n').split('1')[0]][goodent_C3]
            if verbose: print(' - Number of S/N(FELIS) as a double check: '+str(len(s2nvalues_C3 )))
            if verbose: print(' - S/N(FELIS): (mean,median,min,max,std) = '+
                              str(np.mean(s2nvalues_C3))+','+
                              str(np.median(s2nvalues_C3))+','+
                              str(np.min(s2nvalues_C3))+','+
                              str(np.max(s2nvalues_C3))+','+
                              str(np.std(s2nvalues_C3)))

            if verbose:
                for gg, gent in enumerate(goodent_C3):
                    print('   Te='+str(TeKey)+' goodent='+str(gent)+'='+str(fluxdat['id'][gent])+
                          '  '+str(FRval_C3[gent])+' +/- '+str(Nsigma*FRerr_C3[gent])+
                          '   ->  n_e(CIII) = '+str(n_e_C3[gg])+'-'+str(n_e_min_C3[gg])+'+'+str(n_e_max_C3[gg]))

            ylabel   = 'CIII1907/CIII1909'
            plotname = outputfile.replace('.txt','_CIII_ne_estimates_Te'+TeKey+'.pdf')
            uves.plot_neForFR(plotname,fout,T_e_fix,FR_C3,ylabel,FRval_C3,FRerr_C3,n_e_C3,n_e_min_C3,n_e_max_C3,
                              fluxdat,goodent_C3,yvals_curve,n_e_Te3_C3,n_e_Te4_C3,n_e_Te5_C3,verbose=True)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # --- SiIII vs CIII ---
        plotname             = outputfile.replace('.txt','_CIIIvsSiIII_ne_estimates_Te'+TeKey+'.pdf')
        xlabel               = 'n$_\\textrm{e}$(CIII) [cm$^{-3}$]'
        ylabel               = 'n$_\\textrm{e}$(SiIII) [cm$^{-3}$]'
        goodentC3Si3         = np.intersect1d(goodent_C3,goodent_Si3).astype(int)
        n_e_goodentC3Si3_C3  = np.zeros(len(goodentC3Si3))-1
        n_e_goodentC3Si3_Si3 = np.zeros(len(goodentC3Si3))-1
        for gg, ge in enumerate(goodentC3Si3):
            n_e_goodentC3Si3_C3[gg]  = np.where(goodent_C3  == ge)[0]
            n_e_goodentC3Si3_Si3[gg] = np.where(goodent_Si3 == ge)[0]

        uves.plot_neVSne(plotname,T_e_fix,
                         n_e_C3[n_e_goodentC3Si3_C3.astype(int)],
                         n_e_min_C3[n_e_goodentC3Si3_C3.astype(int)],
                         n_e_max_C3[n_e_goodentC3Si3_C3.astype(int)],xlabel,
                         n_e_Si3[n_e_goodentC3Si3_Si3.astype(int)],
                         n_e_min_Si3[n_e_goodentC3Si3_Si3.astype(int)],
                         n_e_max_Si3[n_e_goodentC3Si3_Si3.astype(int)],ylabel,
                         fluxdat,goodentC3Si3,verbose=True,overwrite=overwrite)


        if verbose: print(' - Number of S/N(FELIS) as a double check: '+str(len(s2nvalues_Si3[n_e_goodentC3Si3_Si3.astype(int)])))
        if verbose: print(' - S/N(FELIS): (mean,median,min,max,std) = '+
                          str(np.mean(s2nvalues_Si3[n_e_goodentC3Si3_Si3.astype(int)]))+','+
                          str(np.median(s2nvalues_Si3[n_e_goodentC3Si3_Si3.astype(int)]))+','+
                          str(np.min(s2nvalues_Si3[n_e_goodentC3Si3_Si3.astype(int)]))+','+
                          str(np.max(s2nvalues_Si3[n_e_goodentC3Si3_Si3.astype(int)]))+','+
                          str(np.std(s2nvalues_Si3[n_e_goodentC3Si3_Si3.astype(int)])))

        if verbose: print(' - Number of S/N(FELIS) as a double check: '+str(len(s2nvalues_C3[n_e_goodentC3Si3_C3.astype(int)])))
        if verbose: print(' - S/N(FELIS): (mean,median,min,max,std) = '+
                          str(np.mean(s2nvalues_C3[n_e_goodentC3Si3_C3.astype(int)]))+','+
                          str(np.median(s2nvalues_C3[n_e_goodentC3Si3_C3.astype(int)]))+','+
                          str(np.min(s2nvalues_C3[n_e_goodentC3Si3_C3.astype(int)]))+','+
                          str(np.max(s2nvalues_C3[n_e_goodentC3Si3_C3.astype(int)]))+','+
                          str(np.std(s2nvalues_C3[n_e_goodentC3Si3_C3.astype(int)]))+'\n\n')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fout.close()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # EMission grids:
    # emgridcalc = outputfile.replace('.txt','_CIII_emissiongrid.pypic')
    # if os.path.isfile(emgridcalc):
    #     C3_EG = pn.EmisGrid(restore_file=emgridcalc) # Restored from a previous computation
    # else:
    #     C3_EG = pn.EmisGrid('C', 3, n_tem=30, n_den=30,
    #                         tem_min=5000., tem_max=20000.,
    #                         den_min=10., den_max=1.e8)
    #     C3_EG.save(emgridcalc)
    #
    # # O3_5007 = O3_EG.getGrid(wave=5007)
    # # O3_Te   = O3_EG.getGrid(to_eval = 'L(4363)/L(5007)')
    # C3_EG.plotImage(to_eval = 'L(1907)/L(1909)')
    # C3_EG.plotContours(to_eval = 'L(1907)/L(1909)')
    #

    # # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # # The complete list of the predefined diagnostics is stored in the pn.diags dictionary and can be listed with:
    # for diag in np.sort(pn.diags_dict.keys()):
    #     print('"{0}" : {1}'.format(diag, pn.diags_dict[diag]))
    #


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_neForFR(plotname,fout,T_e_fix,FR,ylabel,FRval,FRvalerr,n_e,n_e_min,n_e_max,
                 fluxdat,goodent,yvals_curve,n_e_Te3,n_e_Te4,n_e_Te5,verbose=True):
    """
    Function to generate plots of n_e. Used in uves.perform_PyNeb_calc_main()

    """
    #------------------------------------------------------------------------------
    if verbose: print(' - Setting up and generating plot of n_e estimates')
    fig = plt.figure(figsize=(5, 4))

    # fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.97, bottom=0.15, top=0.97)
    left, width = 0.15, 0.60
    bottom, height = 0.15, 0.60
    bottom_h = left_h = left + width + 0.01
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=left, right=0.99, bottom=bottom, top=bottom+height)
    rect_histx = [left, bottom_h, 0.705, 0.2]

    Fsize    = 12
    lthick   = 1.0
    marksize = 12
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(inforstr[:-2],fontsize=Fsize)

    plt.plot(n_e_Te5,yvals_curve,color='darkred',zorder=5,lw=lthick,label='n$_\\textrm{e}$(T$_\\textrm{e}$ = 20000K)')
    plt.plot(n_e_Te4,yvals_curve,color='indianred',zorder=5,lw=lthick,label='n$_\\textrm{e}$(T$_\\textrm{e}$ = 10000K)')
    plt.plot(n_e_Te3,yvals_curve,color='salmon',zorder=5,lw=lthick,label='n$_\\textrm{e}$(T$_\\textrm{e}$ = 5000K)')

    #----------------------------
    cmap    = plt.cm.viridis_r
    cdatvec = fluxdat['redshift'][goodent]
    clabel  = '$z$'
    cmin    = 1.5 #0.01 # 0.0, 1.4
    cmax    = 4.0 #7.5 # 10.2, 6.2
    cextend = 'neither'

    colnorm = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, num=50)
    m       = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(cmaparr)

    colvec   = []
    for ii,xval in enumerate(n_e):
        colvec.append(cmap(colnorm(cdatvec[ii])))

    colshrink = 1.0
    colaspect = 30
    colanchor = (0.0,0.5)
    cb      = plt.colorbar(m,extend=cextend,orientation='vertical',
                           pad=0.01,aspect=colaspect,shrink=colshrink,anchor=colanchor,use_gridspec=False)
    cb.set_label(clabel)
    #----------------------------
    FRvallist = np.unique(FRval[np.isfinite(FRval)])
    voffsets_indices = {}
    for FRvl in FRvallist:
        voffsets_indices[str(FRvl)] = 0

    plt.xscale('log')
    plt.ylabel(ylabel)
    plt.xlabel('n$_\\textrm{e}$ [cm$^{-3}$]')
    plt.xlim([1e2,8e6])
    plt.ylim([-0.1,2.2])

    limsizefrac = 0.05
    xvalues     = n_e
    vallim      = np.zeros(len(xvalues))
    for nn, obj_ne in enumerate(xvalues):
        mfc      = colvec[nn]#'navy'
        xlolims  = False
        xuplims  = False
        xerr_min = n_e[nn]-n_e_min[nn]
        xerr_max = n_e_max[nn]-n_e[nn]
        xerrshow = [xerr_min,xerr_max]

        xvalshow = n_e[nn]

        if ~np.isfinite(n_e_min[nn]) & ~np.isfinite(n_e_max[nn]):
            n_e[nn]     = 7e6
            n_e_min[nn] = 0
            n_e_max[nn] = 0
            xuplims     = False

            vallim[nn]  = 1
            dlog     = np.abs(np.diff(np.log10(plt.xlim()))) * limsizefrac
            xvalshow = n_e[nn]#/Nsigma * 1.0
            xerrshow = 0.0#np.abs(xvalshow - 10.**(np.log10(xvalshow)-dlog))

            mfc         = 'None'

        elif ~np.isfinite(n_e_min[nn]):
            n_e[nn]     = n_e_max[nn]
            n_e_min[nn] = +99
            n_e_max[nn] = +99

            xuplims     = True
            vallim[nn]  = 1
            dlog     = np.abs(np.diff(np.log10(plt.xlim()))) * limsizefrac
            xvalshow = n_e[nn]#/Nsigma * 1.0
            xerrshow = np.abs(xvalshow - 10.**(np.log10(xvalshow)-dlog))

        elif ~np.isfinite(n_e_max[nn]):
            n_e[nn]     = n_e_min[nn]
            n_e_min[nn] = -99
            n_e_max[nn] = -99

            xlolims     = True
            vallim[nn]  = 1
            xvalshow = n_e[nn]#/Nsigma / 1.0
            dlog     = np.abs(np.diff(np.log10(plt.xlim()))) * limsizefrac
            xerrshow = np.abs(xvalshow - 10.**(np.log10(xvalshow)+dlog))


        yvalshow = FRval[goodent][nn] #+ voffsets_indices[str(FRval[goodent][nn])]
        #voffsets_indices[str(FRval[goodent][nn])] = voffsets_indices[str(FRval[goodent][nn])] + 0.03

        # if yvalshow ==  1.5: pdb.set_trace()
        plt.errorbar(xvalshow,yvalshow,
                     xerr=[xerrshow],
                     yerr=FRvalerr[goodent][nn],
                     uplims=False,lolims=False,
                     xuplims=xuplims,xlolims=xlolims,
                     marker='.',lw=lthick, markersize=marksize,alpha=1.0,
                     markerfacecolor=mfc,ecolor=colvec[nn],
                     markeredgecolor=colvec[nn],zorder=10)

        fout.write(str(fluxdat['id'][goodent][nn])+'  '+
                   str("%15.2f" % n_e[nn])+'  '+str("%15.2f" % n_e_min[nn])+'  '+str("%15.2f" % n_e_max[nn])+'  '+
                   str("%15.2f" % T_e_fix)+'  pyneb.getTemDen('+FR+')  \n')

    leg = plt.legend(fancybox=True, loc='lower left',prop={'size':Fsize/1.3},ncol=1,numpoints=1)#,
                     # bbox_to_anchor=(0.5, 1.1),)  # add the legend
    leg.get_frame().set_alpha(0.7)

    # - - - - - - Histogram of n_e  - - - - - -
    xminsys, xmaxsys = plt.xlim() # use to get automatically expanded axes if xmin = xmax
    try:
        Nbins   = np.ceil(np.sqrt(len(xvalshow)))
    except:
        Nbins = 10
    axHistx = plt.axes(rect_histx)
    axHistx.xaxis.set_major_formatter(NullFormatter())

    binwidth_x = np.diff([xminsys,xmaxsys])/Nbins
    bindefs    = np.arange(xminsys, xmaxsys+binwidth_x, binwidth_x)
    bindefs    = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))
    axHistx.set_xscale('log')

    axHistx.hist(xvalues[np.isfinite(xvalues)], bins=bindefs,histtype='step',color='k',linestyle='-')
    axHistx.hist(xvalues[np.isfinite(xvalues) & (vallim == 0)], bins=bindefs,histtype='stepfilled',color='k',linestyle='-')
    axHistx.set_xticks([])
    axHistx.set_xlim([xminsys,xmaxsys])
    # - - - - - - - - - - - - - - - - - - - - -


    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')
    #------------------------------------------------------------------------------

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_neVSne(plotname,T_e_fix,
                n_e_1,n_e_min_1,n_e_max_1,xlabel,
                n_e_2,n_e_min_2,n_e_max_2,ylabel,
                fluxdat,goodent,verbose=True, overwrite=False):
    """
    Function to generate plots of two estimates of n_e. Used in uves.perform_PyNeb_calc_main()

    """
    Nhistbins    = 50
    # xrange       = [1e1,1e7]
    xrange       = [1e2,8e6]
    yrange       = [1e2,1e6]
    limsizefrac  = 0.05

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xvalshow_all = []
    xerrshow_all = []
    for nn, obj_ne in enumerate(n_e_1):
        xerr_min = n_e_1[nn]-n_e_min_1[nn]
        xerr_max = n_e_max_1[nn]-n_e_1[nn]
        xerrshow = [xerr_min,xerr_max]
        xvalshow = n_e_1[nn]

        if ~np.isfinite(n_e_min_1[nn]) & ~np.isfinite(n_e_max_1[nn]):
            print('WARNING - both entries for value 1 errors are not finite in uves.plot_neVSne()')
            pdb.set_trace()
        elif ~np.isfinite(n_e_min_1[nn]):
            n_e_1[nn]     = n_e_max_1[nn]
            n_e_min_1[nn] = +99
            n_e_max_1[nn] = +99

            xvalshow = n_e_1[nn]
            xerrshow = +99
        elif n_e_min_1[nn] == +99:
            xvalshow  = n_e_1[nn]
            xerrshow  = +99
        elif ~np.isfinite(n_e_max_1[nn]):
            n_e_1[nn]     = n_e_min_1[nn]
            n_e_min_1[nn] = -99
            n_e_max_1[nn] = -99

            xvalshow = n_e_1[nn]
            xerrshow = -99
        elif n_e_max_1[nn] == -99:
            xvalshow  = n_e_1[nn]
            xerrshow  = -99

        xvalshow_all.append(xvalshow)
        xerrshow_all.append(xerrshow)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    yvalshow_all = []
    yerrshow_all = []
    for nn, obj_ne in enumerate(n_e_2):
        yerr_min = n_e_2[nn]-n_e_min_2[nn]
        yerr_max = n_e_max_2[nn]-n_e_2[nn]
        yerrshow = [yerr_min,yerr_max]
        yvalshow = n_e_2[nn]

        if ~np.isfinite(n_e_min_2[nn]) & ~np.isfinite(n_e_max_2[nn]):
            print('WARNING - both entries for value 2 errors are not finite in uves.plot_neVSne()')
            pdb.set_trace()
        elif ~np.isfinite(n_e_min_2[nn]):
            n_e_2[nn]     = n_e_max_2[nn]
            n_e_min_2[nn] = +99
            n_e_max_2[nn] = +99

            yvalshow = n_e_2[nn]
            yerrshow = +99
        elif n_e_min_2[nn] == +99:
            yvalshow  = n_e_2[nn]
            yerrshow  = +99
        elif ~np.isfinite(n_e_max_2[nn]):
            n_e_2[nn]     = n_e_min_2[nn]
            n_e_min_2[nn] = -99
            n_e_max_2[nn] = -99

            yvalshow = n_e_2[nn]
            yerrshow = -99
            print(str(yerrshow)+'  '+str(nn))
        elif n_e_max_2[nn] == -99:
            yvalshow  = n_e_2[nn]
            yerrshow  = -99

        yvalshow_all.append(yvalshow)
        yerrshow_all.append(yerrshow)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xvalshow_all = np.asarray(xvalshow_all)
    yvalshow_all = np.asarray(yvalshow_all)

    if verbose:
        print(' - Objectinfo (Nobj='+str(len(xvalshow_all))+') for plot:')
        print('   IDs: '+str(fluxdat['id'][goodent]))

    Nhistbins   = np.ceil(np.sqrt(len(xvalshow_all)))
    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalshow_all,yvalshow_all,xerrshow_all,yerrshow_all,xlabel,ylabel,
                                                   'dummydat',linetype='onetoone',title=None, #'this is title',
                                                   ids=fluxdat['id'][goodent],
                                                   ylog=True,xlog=True,yrange=yrange,xrange=xrange,
                                                   colortype='zmanual',colorcode=True,cdatvec=fluxdat['redshift'][goodent],
                                                   point_text=None, #fluxdat['id'][goodent].astype(str),
                                                   photoionizationplotparam=None,
                                                   histaxes=False,Nbins=Nhistbins, showgraylimits=False,
                                                   overwrite=overwrite,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_magnitudedistributions(outputdir,infofile,masterfits, emlinelist = ['CIV','HeII','OIII','SiIII','CIII'],
                                performlinearfit=True, EWlimprint=10.0, showlimits=False, verbose=True,
                                overwrite=False, addidlabels=False):
    """
    Function to generate plots of two estimates of n_e. Used in uves.perform_PyNeb_calc_main()

    --- INPUT ---
    outputdir       The directory to store the generated plots in.
    infofile        The UVES object main infofile with object parameters.
    masterfits      The UVES master fits file with EW estimates.
    emlinelist      List of emission lines to include in plots; available: ['CIV','HeII','OIII','SiIII','CIII','MgII']
    addidlabels     To display the IDs in the plots set to True. This will also print ID lists of extreme EW objects.
    EWlimprint      The lower limit of the extreme EW emitters to plot with "addidlabels".
    showlimits      Show the upper/lower limits in the plots.
    verbose         Toggle verbosity
    overwrite       Overwrite existing figures?

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    uvesdir            = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'
    outdir_magdist     = uvesdir+'magnitudedistributions/'
    infofile           = uvesdir+'objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    masterfits         = uvesdir+'back2backAnalysis_200213/results_master_catalog_version200213.fits'

    uves.plot_magnitudedistributions(outdir_magdist,infofile,masterfits,emlinelist = ['NV','CIV','HeII','OIII','SiIII','CIII','MgII'])
    """
    masterdat = afits.open(masterfits)[1].data
    infodat   = afits.open(infofile)[1].data
    infodat   = infodat[np.where((infodat['id']<4.9e8) | (infodat['id']>5.9e8))[0]] # ignoring UDF MW mock ids

    for ii, id in enumerate(masterdat['id']): # double checking that order of objects is the same
        if infodat['id'][ii] != id:
            sys.exit('There was a mismatch in ids for entry '+str(ii))

    #------------------------------------------------------------------------------
    if verbose: print(' - Collecting indexes of magnitudes for line sets to plot ')
    magdic = collections.OrderedDict()

    if showlimits is True:
        limval = 0.0
    else:
        limval = 99

    for emline in emlinelist:

        magAB_ent   = np.where((masterdat['duplicationID'] == 0) &
                               (masterdat['contmagABerr_'+emline] != -99) & # ignore lower limits (upper limits on brightness)
                               (np.abs(masterdat['EW0err_'+emline]) != limval) & (masterdat['EW0_'+emline] != 0.0) &
                               (masterdat['EW0err_'+emline] != 99) & # don't show upper limits on EW.
                               np.isfinite(masterdat['contmagABerr_'+emline]))[0]
        magAB       = masterdat['contmagAB_'+emline][magAB_ent]
        magABerr    = masterdat['contmagABerr_'+emline][magAB_ent]
        EW          = masterdat['EW0_'+emline][magAB_ent]
        EWerr       = masterdat['EW0err_'+emline][magAB_ent]
        redshift    = masterdat['redshift'][magAB_ent]
        objids      = masterdat['id'][magAB_ent]

        magdic[emline] = {'mag':magAB, 'magerr':magABerr, 'EW':EW, 'EWerr':EWerr, 'redshift':redshift, 'id':objids}

    #------------------------------------------------------------------------------
    if verbose: print(' - Setting up and generating plots for each emission line')
    Nhistbins = 30
    xrange    = [19,31.5]
    yrange    = [0.1,400]
    colortype = 'redshift' #'zmanual' # 'z' fixes redshift range in color
    for emline in magdic.keys():
        plotname   = outputdir+'magdist_'+emline+'.pdf'
        # xlabel     = 'EW$_0$ continuum AB magnitude ('+emline+')'
        xlabel     = 'Continuum AB magnitude around '+emline
        xvalues    = magdic[emline]['mag']
        xerr       = magdic[emline]['magerr']
        ylabel     = 'EW$_0$('+emline+') [\AA]'
        yvalues    = magdic[emline]['EW']
        yerr       = magdic[emline]['EWerr']
        cdatvec    = magdic[emline]['redshift']
        if addidlabels:
            point_text  = magdic[emline]['id'].astype(str)
            ewGT100_ebt = np.where(yvalues >= EWlimprint)[0]
            print(' - The list of objects with EW_0('+emline+') >= '+str(EWlimprint)+':\n   '+str(point_text[ewGT100_ebt]))
        else:
            point_text = None

        ylog = True
        xlog = False
        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,
                                                       'dummydat',linetype=None,title=None,ids=magdic[emline]['id'],
                                                       ylog=ylog,xlog=xlog,yrange=yrange,xrange=xrange,
                                                       colortype=colortype,colorcode=True,cdatvec=cdatvec,
                                                       point_text=point_text,photoionizationplotparam=None,
                                                       histaxes=True,Nbins=Nhistbins,
                                                       overwrite=overwrite,verbose=verbose)

        if performlinearfit:
            plotname   = plotname.replace('.pdf','_ODRfit2data.pdf')

            if ylog:
                yval_fit = np.log10(yvalues)
                yerr_fit = 0.434 * yerr/yvalues  # https://faculty.washington.edu/stuve/log_error.pdf
            else:
                yval_fit = yvalues
                yerr_fit = yerr

            if xlog:
                xval_fit = np.log10(xvalues)
                xerr_fit = 0.434 * xerr/xvalues  # https://faculty.washington.edu/stuve/log_error.pdf
            else:
                xval_fit = xvalues
                xerr_fit = xerr

            fitres,rp,pp,rs,ps = kbs.fit_function_to_data_with_errors_on_both_axes(xval_fit,yval_fit,xerr_fit,yerr_fit,
                                                                                   fitfunction='linear',plotresults=plotname,
                                                                                   returnCorrelationCoeffs=True)
        yvaluesWflux     = yvalues*10**(xvalues/-2.5)
        yvaluesWflux_err = np.abs(yvaluesWflux) * np.log(10)/2.5 * xerr
        uves.plot_mocspecFELISresults_summary_plotcmds(plotname.replace('.pdf','_EWxFlux.pdf'),
                                                       xvalues,yvaluesWflux,xerr,yvaluesWflux_err,
                                                       xlabel,ylabel+' x 10\^('+xlabel+'/-2.5)',
                                                       'dummydat',linetype=None,title=None,ids=magdic[emline]['id'],
                                                       ylog=ylog,xlog=xlog,
                                                       yrange=[1e-11,1e-7],xrange=xrange,
                                                       # yrange=[yrange[0]*10**(xrange[1]/-2.5),yrange[1]*10**(xrange[0]/-2.5)],xrange=xrange,
                                                       colortype=colortype,colorcode=True,cdatvec=cdatvec,
                                                       point_text=point_text,photoionizationplotparam=None,
                                                       histaxes=True,Nbins=Nhistbins,
                                                       overwrite=overwrite,verbose=verbose)

        if performlinearfit:
            plotname   = plotname.replace('.pdf','_ODRfit2data.pdf')
            if ylog:
                yval_fit = np.log10(yvaluesWflux)
                yerr_fit = 0.434 * yvaluesWflux_err/yvaluesWflux  # https://faculty.washington.edu/stuve/log_error.pdf
            else:
                yval_fit = yvaluesWflux
                yerr_fit = yvaluesWflux_err

            if xlog:
                xval_fit = np.log10(xvalues)
                xerr_fit = 0.434 * xerr/xvalues  # https://faculty.washington.edu/stuve/log_error.pdf
            else:
                xval_fit = xvalues
                xerr_fit = xerr


            fitres,rp,pp,rs,ps = kbs.fit_function_to_data_with_errors_on_both_axes(xval_fit,yval_fit,xerr_fit,yerr_fit,
                                                                                   fitfunction='linear',plotresults=plotname,
                                                                                   returnCorrelationCoeffs=True)
    #------------------------------------------------------------------------------
    linecolors  = uves.linecolors(colormap='plasma')
    if verbose: print(' - Setting up and generating plot of histograms in same window')
    plotname = outputdir+'magdist_linescombined.pdf'

    if os.path.isfile(plotname) & (overwrite == False):
        print('\n - WARNING: the plot '+plotname+' exists and overwrite=False so moving on \n')
    else:
        fig = plt.figure(figsize=(5, 4))

        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.97, bottom=0.15, top=0.97)
        # left, width = 0.15, 0.60
        # bottom, height = 0.15, 0.60
        # bottom_h = left_h = left + width + 0.01
        # fig.subplots_adjust(wspace=0.1, hspace=0.1,left=left, right=0.99, bottom=bottom, top=bottom+height)
        # rect_histx = [left, bottom_h, 0.705, 0.2]

        Fsize    = 12
        lthick   = 1.0
        marksize = 6
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        #plt.title(inforstr[:-2],fontsize=Fsize)

        # - - - - - - Histogram of n_e  - - - - - -
        xminsys, xmaxsys = xrange
        binwidth_x = np.diff([xminsys,xmaxsys])/Nhistbins
        bindefs    = np.arange(xminsys, xmaxsys+binwidth_x, binwidth_x)
        #bindefs    = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))

        for emline in magdic.keys():
            histdata = magdic[emline]['mag']
            plt.hist(histdata, bins=bindefs,histtype='stepfilled',color=linecolors[emline],alpha=0.3,linestyle='-',label=emline)

        plt.ylabel('Number of objects')
        # plt.xlabel('EW$_0$ continuum AB magnitude')
        plt.xlabel('Continuum AB magnitude around line')
        plt.xlim(xrange)
        plt.ylim(0,11.0)
        leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)#,
                         # bbox_to_anchor=(0.5, 1.1),)  # add the legend
        leg.get_frame().set_alpha(0.7)

        if verbose: print('   Saving plot to '+plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')

    #------------------------------------------------------------------------------
    if verbose: print(' - Setting up and generating plot of histograms in same window')
    plotname = outputdir+'magdist_linescombined_sharedaxis.pdf'
    if os.path.isfile(plotname) & (overwrite == False):
        print('\n - WARNING: the plot '+plotname+' exists and overwrite=False so moving on \n')
    else:
        fig = plt.figure(figsize=(3, 5))
        fig.subplots_adjust(wspace=0.1, hspace=0.2,left=0.2, right=0.97, bottom=0.1, top=0.99)

        Fsize    = 10
        lthick   = 1.0
        marksize = 6
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        #plt.title(inforstr[:-2],fontsize=Fsize)

        # - - - - - - Histogram  - - - - - -
        xminsys, xmaxsys = xrange
        binwidth_x = np.diff([xminsys,xmaxsys])/Nhistbins
        bindefs    = np.arange(xminsys, xmaxsys+binwidth_x, binwidth_x)
        # bindefs    = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))

        axfull = fig.add_subplot(111)
        # axfull.set_xlabel('EW$_0$ continuum AB magnitude')
        axfull.set_xlabel('Continuum AB magnitude around line')
        axfull.set_ylabel('Number of objects')
        # Turn off axis lines and ticks of the full window subplot
        axfull.spines['top'].set_color('none')
        axfull.spines['bottom'].set_color('none')
        axfull.spines['left'].set_color('none')
        axfull.spines['right'].set_color('none')
        axfull.tick_params(labelcolor='white', top=False, bottom=False, left=False, right=False)

        # - - - - - - - - - - - - - -
        for ee, emline in enumerate(magdic.keys()):
            if ee == 0:
                ax       = fig.add_subplot(len(emlinelist)*100+11+ee)
            else:
                ax       = fig.add_subplot(len(emlinelist)*100+11+ee, sharex=ax, sharey=ax)
            histdata = magdic[emline]['mag']
            ax.hist(histdata, bins=bindefs,histtype='stepfilled',color=linecolors[emline],alpha=0.9,linestyle='-',label=emline,zorder=100)
            ax.grid(True,zorder=10,color='lightgrey')
            ax.set_xticks(np.arange(19,31,2))
            ax.set_yticks(np.arange(0,11,2))

            leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)#,
                             # bbox_to_anchor=(0.5, 1.1),)  # add the legend
            leg.get_frame().set_alpha(0.7)
        # - - - - - - - - - - - - - -

        for ax in fig.get_axes():
            ax.label_outer()

        plt.xlim(xrange)
        plt.ylim(0,11.0)


        if verbose: print('   Saving plot to '+plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')

    #------------------------------------------------------------------------------
    if verbose: print(' - Collecting Lya absolute mags')
    goodent  = np.where((infodat['abs_mag_UV_cont_own_median_jk100'] != 0) & (infodat['duplicationID'] == 0))[0]

    mag                 = infodat['abs_mag_UV_cont_own_median_jk100'][goodent]
    magerr              = None
    EW                  = infodat['EW_0_beta_own_median_jk100'][goodent]
    EWerr               = infodat['EW_0_beta_own_median_error_jk100'][goodent]
    if showlimits:
        EWerr[EWerr == 0.0] = -99
    else:
        EW[EWerr == 0.0]    = np.nan
        EWerr[EWerr == 0.0] = np.nan
    redshift            = infodat['redshift'][goodent]

    plotname = outputdir+'magdist_LAE_MUV.pdf'
    xlabel    = 'M$_\\textrm{UV,1500\AA}$'
    xvalues   = mag[np.isfinite(EW)]
    xerr      = None #magerr[np.isfinite(EW)]
    ylabel    = 'EW$_0$(Ly$\\alpha$) [\AA] ($\\beta = -1.97$)'
    yvalues   = EW[np.isfinite(EW)]
    yerr      = EWerr[np.isfinite(EW)]
    cdatvec   = redshift[np.isfinite(EW)]

    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,yvalues,xerr,yerr,xlabel,ylabel,
                                                   'dummydat',linetype=None,title=None,ids=infodat['id'][goodent],
                                                   ylog=True,xlog=False,yrange=[0.1,2000],xrange=[-24.9,-14.9],
                                                   colortype=colortype,colorcode=True,cdatvec=cdatvec,
                                                   point_text=None,photoionizationplotparam=None,
                                                   histaxes=True,Nbins=Nhistbins,
                                                   overwrite=overwrite,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def evaluate_velocityoffsets(linefluxcatalog,infofile,outputdir='./velocityoffsetFigures/',addDvErr=False,
                             Nlinedetections=[2,10],overwrite=False,verbose=True):
    """
    Main function to handle caluclations and diagostics using PyNeb

    --- Example of use ---
    import uvEmissionlineSearch as uves
    uvesdir          = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'
    outdir           = uvesdir+'velocityoffsetFigures/'
    linefluxcatalog  = uvesdir+'back2backAnalysis_200213/results_master_catalog_version200213.fits'
    # infofile         = uvesdir+'objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile         = uvesdir+'objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    uves.evaluate_velocityoffsets(linefluxcatalog,infofile,outputdir=outdir,Nlinedetections=[2,10],overwrite=True,addDvErr=True)

    uves.evaluate_velocityoffsets(linefluxcatalog,infofile,outputdir=outdir,Nlinedetections=[2,3],overwrite=True,addDvErr=True)
    uves.evaluate_velocityoffsets(linefluxcatalog,infofile,outputdir=outdir,Nlinedetections=[3,4],overwrite=True,addDvErr=True)
    uves.evaluate_velocityoffsets(linefluxcatalog,infofile,outputdir=outdir,Nlinedetections=[4,9],overwrite=True,addDvErr=True)

    """
    agn, agncand = uves.get_AGN_ids()
    dat_uves = afits.open(linefluxcatalog)[1].data
    infodat  = afits.open(infofile)[1].data
    infodat  = infodat[np.where((infodat['id']<4.9e8) | (infodat['id']>5.9e8))[0]] # ignoring UDF MW mock ids

    for ii, id in enumerate(dat_uves['id']): # double checking that order of objects is the same
        if infodat['id'][ii] != id:
            sys.exit('There was a mismatch in ids for entry '+str(ii))

    linenamelist = ['NV','CIV','HeII','OIII','SiIII','CIII','MgII']
    Nlines       = len(linenamelist)
    # zlineranges  = [[2.8918,6.4432],[2.1114,4.9699],[1.9379,4.6411],[1.8969,4.5632],[1.5514,3.9067],[1.5241,3.8548],[0.7174,2.3142]]
    zlineranges  = [[2.86,6.5],[2.0,5.1],[1.7,4.8],[1.7,4.7],[1.4,3.95],[1.4,3.95],[0.7174,2.5]]
    #\lya                		&  2.9729 -- 6.5958 	&

    # - - - - - - - ALL lead lines - - - - - - -
    dv_NV_ent    = np.where((dat_uves['vshift_NV'] != 99) & (dat_uves['duplicationID'] == 0) &
                          np.isfinite(dat_uves['vshift_NV']))[0]
    dv_NV        = dat_uves['vshift_NV'][dv_NV_ent]

    dv_CIV_ent   = np.where((dat_uves['vshift_CIV'] != 99) & (dat_uves['duplicationID'] == 0) &
                          np.isfinite(dat_uves['vshift_CIV']))[0]
    dv_CIV       = dat_uves['vshift_CIV'][dv_CIV_ent]

    dv_HeII_ent  = np.where((dat_uves['vshift_HeII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['vshift_HeII']))[0]
    dv_HeII      = dat_uves['vshift_HeII'][dv_HeII_ent]

    dv_OIII_ent  = np.where((dat_uves['vshift_OIII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['vshift_OIII']))[0]
    dv_OIII      = dat_uves['vshift_OIII'][dv_OIII_ent]

    dv_SiIII_ent = np.where((dat_uves['vshift_SiIII'] != 99) & (dat_uves['duplicationID'] == 0) &
                            np.isfinite(dat_uves['vshift_SiIII']))[0]
    dv_SiIII     = dat_uves['vshift_SiIII'][dv_SiIII_ent]

    dv_CIII_ent  = np.where((dat_uves['vshift_CIII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['vshift_CIII']))[0]
    dv_CIII      = dat_uves['vshift_CIII'][dv_CIII_ent]

    dv_MgII_ent  = np.where((dat_uves['vshift_MgII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['vshift_MgII']))[0]
    dv_MgII      = dat_uves['vshift_MgII'][dv_MgII_ent]

    # - - - - - - - Only Lya as lead line - - - - - - -
    dv_NV_ent_LAE    = np.where((dat_uves['vshift_NV'] != 99) & (dat_uves['duplicationID'] == 0) &
                          np.isfinite(dat_uves['vshift_NV']) & (infodat['leadline'] == 'Lya'))[0]
    dv_NV_LAE        = dat_uves['vshift_NV'][dv_NV_ent_LAE]

    dv_CIV_ent_LAE   = np.where((dat_uves['vshift_CIV'] != 99) & (dat_uves['duplicationID'] == 0) &
                          np.isfinite(dat_uves['vshift_CIV']) & (infodat['leadline'] == 'Lya'))[0]
    dv_CIV_LAE       = dat_uves['vshift_CIV'][dv_CIV_ent_LAE]

    dv_HeII_ent_LAE  = np.where((dat_uves['vshift_HeII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['vshift_HeII']) & (infodat['leadline'] == 'Lya'))[0]
    dv_HeII_LAE      = dat_uves['vshift_HeII'][dv_HeII_ent_LAE]

    dv_OIII_ent_LAE  = np.where((dat_uves['vshift_OIII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['vshift_OIII']) & (infodat['leadline'] == 'Lya'))[0]
    dv_OIII_LAE      = dat_uves['vshift_OIII'][dv_OIII_ent_LAE]

    dv_SiIII_ent_LAE = np.where((dat_uves['vshift_SiIII'] != 99) & (dat_uves['duplicationID'] == 0) &
                            np.isfinite(dat_uves['vshift_SiIII']) & (infodat['leadline'] == 'Lya'))[0]
    dv_SiIII_LAE     = dat_uves['vshift_SiIII'][dv_SiIII_ent_LAE]

    dv_CIII_ent_LAE  = np.where((dat_uves['vshift_CIII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['vshift_CIII']) & (infodat['leadline'] == 'Lya'))[0]
    dv_CIII_LAE      = dat_uves['vshift_CIII'][dv_CIII_ent_LAE]

    dv_MgII_ent_LAE  = np.where((dat_uves['vshift_MgII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['vshift_MgII']) & (infodat['leadline'] == 'Lya'))[0]
    dv_MgII_LAE      = dat_uves['vshift_MgII'][dv_MgII_ent_LAE]

    # - - - - - - - Special cases - - - - - - -
    dv_CIIICIV_ent  = np.intersect1d(dv_CIV_ent,dv_CIII_ent)
    dv_CIIICIV      = dat_uves['vshift_CIII'][dv_CIIICIV_ent] - dat_uves['vshift_CIV'][dv_CIIICIV_ent]

    print('CIV-CIII offsets:'+str(dv_CIIICIV))
    print('And corresponding IDs:'+str(dat_uves['id'][dv_CIIICIV_ent]))

    DvLyaLit_cat = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/vshift_Lya_and_MUV.txt'
    DvLyaLit_dat = np.genfromtxt(DvLyaLit_cat,skip_header=0,names=True,dtype='d,d,d,40a,d,d,d,d,d,d,d,d,d,d,d,40a')

    #------------------------------------------------------------------------------
    if verbose: print(' - Setting up and generating plot of lead line offsets intra-object comparison')
    plotname = outputdir+'evaluate_voffsets_intraobject_comparison_all.pdf'
    fig      = plt.figure(figsize=(6, 4))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.99, bottom=0.1, top=0.97)
    Fsize    = 12
    lthick   = 1.0
    marksize = 6
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(inforstr[:-2],fontsize=Fsize)

    #----------------------------
    cmap    = plt.cm.viridis_r

    clabel  = '$z$(lead line)'
    cmin    = 1.5
    cmax    = 6.5
    cextend = 'neither'

    colnorm = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, num=50)
    m       = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(cmaparr)

    # cdatvec = dat_uves['redshift']
    # colvec   = []
    # for ii,xval in enumerate(n_e):
    #     colvec.append(cmap(colnorm(cdatvec[ii])))

    colshrink = 1.0
    colaspect = 30
    colanchor = (0.0,0.5)
    cb      = plt.colorbar(m,extend=cextend,orientation='vertical',
                           pad=0.02,aspect=colaspect,shrink=colshrink,anchor=colanchor,use_gridspec=False)
    cb.set_label(clabel)
    #----------------------------

    xticklabel = linenamelist
    xvals      = np.arange(Nlines)+1.0
    dvs_ent    = [dv_NV_ent,dv_CIV_ent,dv_HeII_ent,dv_OIII_ent,dv_SiIII_ent,dv_CIII_ent,dv_MgII_ent]
    dvs        = [dv_NV,dv_CIV,dv_HeII,dv_OIII,dv_SiIII,dv_CIII,dv_MgII]

    ymin, ymax = [-1,1]
    for objent, objid in enumerate(dat_uves['id']):
        yvals = np.asarray([np.nan]*Nlines)

        detectiontracker = [0.0]*Nlines
        for ee, entlist in enumerate(dvs_ent):
            if objent in entlist:
                detectiontracker[ee] = 1.0
                obj_dvs_ent          = np.where(entlist == objent)[0]
                yvals[ee]            = dvs[ee][obj_dvs_ent][0]
                # if np.abs(yvals[ee]) > 1000: pdb.set_trace()

        if (np.sum(detectiontracker) >= Nlinedetections[0]) & (np.sum(detectiontracker) < Nlinedetections[1]):
            objcol = cmap(colnorm(dat_uves['redshift'][objent]))

            if addDvErr:
                yerr = 100.0
            else:
                yerr = None

            if objid in agn:
                objmarker = '*'
                objmarkersize = marksize+2
            elif objid in agncand:
                objmarker = 'D'
                objmarkersize = marksize
            else:
                objmarker = 'o'
                objmarkersize = marksize

            plt.errorbar(xvals,yvals,xerr=None,yerr=yerr,color=objcol,
                         marker=objmarker,lw=lthick+1, markersize=objmarkersize,alpha=1.0,
                         markerfacecolor=objcol,ecolor=objcol,
                         markeredgecolor=objcol,zorder=10)

            if ~np.isfinite(yvals[0]):
                yvals[0] = yvals[np.isfinite(yvals)][0]
            if ~np.isfinite(yvals[-1]):
                yvals[-1] = yvals[np.isfinite(yvals)][-1]
            plt.plot(xvals[np.isfinite(yvals)],yvals[np.isfinite(yvals)],':',color=objcol,#'lightgray',
                     zorder=4,lw=lthick)

            plt.text(0,yvals[0],str(objid),color=objcol,fontsize=Fsize,rotation=0,
                     horizontalalignment='center',verticalalignment='center',zorder=10)

            ymin = np.min([ymin,np.min(yvals[np.isfinite(yvals)]-110.)])
            ymax = np.max([ymax,np.max(yvals[np.isfinite(yvals)]+110.)])

    plt.ylabel('$\Delta v$(lead line)')
    plt.xlabel(' ')
    plt.xlim([-1,7.3])
    plt.ylim([ymin,ymax])

    plt.xticks(xvals, xticklabel, rotation='horizontal')

    # leg = plt.legend(fancybox=True, loc='lower left',prop={'size':Fsize/1.3},ncol=1,numpoints=1)#,
    #                  # bbox_to_anchor=(0.5, 1.1),)  # add the legend
    # leg.get_frame().set_alpha(0.7)

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

    #------------------------------------------------------------------------------
    if verbose: print(' - Setting up and generating plots of lead line offsets ')

    linename   = linenamelist+['CIV']+linenamelist
    zranges    = zlineranges+[[2.0,3.95]]+zlineranges
    leadline   = ['lead line']*Nlines + ['CIII'] + ['Lya']*Nlines
    dvs_ent    = [dv_NV_ent,dv_CIV_ent,dv_HeII_ent,dv_OIII_ent,dv_SiIII_ent,dv_CIII_ent,dv_MgII_ent,dv_CIIICIV_ent]+\
                 [dv_NV_ent_LAE,dv_CIV_ent_LAE,dv_HeII_ent_LAE,dv_OIII_ent_LAE,dv_SiIII_ent_LAE,dv_CIII_ent_LAE]#,dv_MgII_ent_LAE]
    dvs        = [dv_NV,dv_CIV,dv_HeII,dv_OIII,dv_SiIII,dv_CIII,dv_MgII,dv_CIIICIV]+\
                 [dv_NV_LAE,dv_CIV_LAE,dv_HeII_LAE,dv_OIII_LAE,dv_SiIII_LAE,dv_CIII_LAE]#,dv_MgII_LAE]

    for vv, dv in enumerate(dvs):
        if not 'CIII' in linename[vv]: continue
        objids    = dat_uves['id'][dvs_ent[vv]]
        zleadline = dat_uves['redshift'][dvs_ent[vv]]
        # ------ Correct individual velocity shifts ------
        for oo, objid in enumerate(objids):
            if zleadline[oo] > 2.9:
                if not objid in [121033078, 601381485, 720470421, 722551008, 722731033, 723311101,           115003085]:
                    znew        = infodat['z_vac_fit_jk100'][dvs_ent[vv]][oo]
                    correction  = 299792.458  * (znew-zleadline[oo])/(1+zleadline[oo])
                    print('   zdiff for '+str(objid)+': '+str(znew-zleadline[oo])+' corresponding to a correction of '+str(correction)+'km/s')
                    dv[oo]        = dv[oo] + correction
                    zleadline[oo] = znew

        # idcorrect   = 603502226
        # objent      = np.where(objids.astype(int) == idcorrect)[0]
        # correction  = infodat['peak_sep_kms_jk100'][dvs_ent[vv]][objent]
        # print('   -> Correcting '+str(idcorrect)+' from '+str(dv[objent])+' to '+str(dv[objent]+correction))
        # dv[objent]  = dv[objent] + correction
        # znewlead    = correction/299792.458*(1+zleadline[objent])+[objent]
        # print('      Corresponding to changing lead line redshift from '+str(zleadline[objent])+' to '+str(znewlead))
        # zleadline[objent] = znewlead
        #
        # idcorrect   = 604992563
        # objent      = np.where(objids.astype(int) == idcorrect)[0]
        # correction  = 299792.458  * (3.803569-zleadline[objent])/(1+zleadline[objent])
        # print('   -> Correcting '+str(idcorrect)+' from '+str(dv[objent])+' to '+str(dv[objent]+correction))
        # dv[objent]  = dv[objent] + correction
        # znewlead    = correction/299792.458*(1+zleadline[objent])+zleadline[objent]
        # print('      Corresponding to changing lead line redshift from '+str(zleadline[objent])+' to '+str(znewlead))
        # zleadline[objent] = znewlead
        # ------------------------------------------------

        histaxes  = True
        #Nhistbins = 50
        yrange    = [-990,990]

        if leadline[vv] == 'Lya':
            llstring = 'Ly$\\alpha$'
        elif leadline[vv] == 'CIII':
            llstring = 'CIII'
        else:
            llstring = leadline[vv]
        ylabel    = '$\Delta v$('+llstring+' - '+linename[vv]+') [km/s]'

        #---- vs redshift ---
        plotname = outputdir+'evaluate_voffsets_'+leadline[vv].replace(' ','')+'-'+linename[vv]+'VSredshiftLeadLine.pdf'

        xrange = zranges[vv]
        if leadline[vv] == 'Lya':
            xrange[0] = np.max([xrange[0],2.8])
            linetype  = 'horizontal'
        else:
            xrange[0] = np.max([xrange[0],1.4])
            linetype  = 'horizontalWlya'

        xlabel   = '$z$('+llstring+')'
        xerr     = None
        if addDvErr:
            yerr = [100.0] * len(dvs_ent[vv])
        else:
            yerr = None

        Nhistbins = np.ceil(np.sqrt(len(dat_uves['redshift'][dvs_ent[vv]])))
        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,dat_uves['redshift'][dvs_ent[vv]],dv,xerr,yerr,xlabel,ylabel,
                                                       'dummydat',linetype=linetype,title=None, #'this is title',
                                                       ids=dat_uves['id'][dvs_ent[vv]],
                                                       ylog=False,xlog=False,yrange=yrange,xrange=xrange,
                                                       colortype='s2nfelis',colorcode=True,
                                                       cdatvec=dat_uves['s2n_'+linename[vv]][dvs_ent[vv]],
                                                       point_text=None, #dat_uves['id'][dvs_ent[vv]].astype(str),
                                                       photoionizationplotparam=None,
                                                       histaxes=histaxes,Nbins=Nhistbins, showgraylimits=True,
                                                       overwrite=overwrite,verbose=verbose)

        Nhistbins = np.ceil(np.sqrt(len(zleadline)))
        plotname = outputdir+'evaluate_voffsets_'+leadline[vv].replace(' ','')+'-'+linename[vv]+'VSredshiftLeadLineAndLyaFit.pdf'
        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,zleadline,dv,xerr,yerr,xlabel,ylabel,
                                                       'dummydat',linetype=linetype,title=None, #'this is title',
                                                       ids=dat_uves['id'][dvs_ent[vv]],
                                                       ylog=False,xlog=False,yrange=yrange,xrange=xrange,
                                                       colortype='s2nfelis',colorcode=True,
                                                       cdatvec=dat_uves['s2n_'+linename[vv]][dvs_ent[vv]],
                                                       point_text=None, #dat_uves['id'][dvs_ent[vv]].astype(str),
                                                       photoionizationplotparam=None,
                                                       histaxes=histaxes,Nbins=Nhistbins, showgraylimits=True,
                                                       overwrite=overwrite,verbose=verbose)

        #---- vs LAE parameters ---
        if leadline[vv] == 'Lya':
            infocols = uves.get_infodat_plotcols()

            #------ Plot as a function of the estiamted Voffset from Verhamme+18 (JK estimates) ------
            if addDvErr:
                yerr = [100.0] * len(dvs_ent[vv])
            else:
                yerr = None

            xlabel   = '$\Delta v_\\textrm{Ly$\\alpha$}$ (Verhamme et al. (2018) approx.)'
            xrange   = [50,440]
            yrangemanual = [-300,710]
            plotname = outputdir+'evaluate_voffsets_'+leadline[vv].replace(' ','')+'-'+linename[vv]+'VSvoffsetVerhammeJK.pdf'

            zsys      = infodat[infocols['zsys'][0]][dvs_ent[vv]]
            zsys_err  = infodat[infocols['zsys'][1]][dvs_ent[vv]]
            xvalues   = 299792.458 * (zleadline-zsys)/(1.0+zsys) # cf. Erb+2014
            xerr      = np.abs(xvalues) * np.sqrt(2*(zsys_err/zsys)**2)

            linetype='onetoone'
            Nhistbins = np.ceil(np.sqrt(len(xvalues)))
            uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,dv,xerr,yerr,xlabel,ylabel,
                                                           'dummydat',linetype=linetype,title=None, #'this is title',
                                                           ids=dat_uves['id'][dvs_ent[vv]],
                                                           ylog=False,xlog=False,yrange=yrangemanual,xrange=xrange,
                                                           colortype='zmanual',colorcode=True,
                                                           cdatvec=dat_uves['redshift'][dvs_ent[vv]],
                                                           point_text=None, #dat_uves['id'][dvs_ent[vv]].astype(str),
                                                           photoionizationplotparam=None,
                                                           histaxes=False,Nbins=Nhistbins, showgraylimits=True,
                                                           overwrite=overwrite,verbose=verbose)
            #-----------------------------------------------------------------------------------------

            for cc, colname in enumerate(infocols.keys()):
                xlabel   = infocols[colname][2]
                xrange   = None
                plotname = outputdir+'evaluate_voffsets_'+leadline[vv].replace(' ','')+'-'+linename[vv]+'VS'+colname+'.pdf'

                if addDvErr:
                    yerr = [100.0] * len(dvs_ent[vv])
                else:
                    yerr = None

                xvalues  = infodat[infocols[colname][0]][dvs_ent[vv]]
                if len(xvalues[np.isfinite(xvalues)]) > 0:
                    if infocols[colname][1] is None:
                        xerr     = np.asarray([np.nan]*len(dvs_ent[vv]))
                    else:
                        xerr     = infodat[infocols[colname][1]][dvs_ent[vv]]

                    plot_dv      = dv
                    plot_ids     = dat_uves['id'][dvs_ent[vv]]
                    plot_cdatvec = zleadline

                    plot_dv_prelit      = plot_dv
                    plot_dv_err_prelit  = yerr

                    if colname == 'lyafwhm_kms':
                        val_lyafwhm_kms       = xvalues
                        val_lyafwhm_kms_err   = xerr
                        plotdvversion = 'zsysFWHM'
                    elif colname == 'peaksep_kms':
                        val_peaksep_kms       = xvalues
                        val_peaksep_kms_err   = xerr
                        plotdvversion = 'zsysPS'
                    else:
                        plotdvversion = None

                    histaxes  = True
                    colortype = 'redshift'
                    if 'EW' in xlabel:
                        linetype  = 'horizontal_and_nakajima18EWvsDv'
                        colortype = 'zmanual'
                        yrange    = [-400,740]
                    elif colname == 'R_e':
                        yrange    = [-200,740]
                        colortype = 'zmanual'
                        linetype  = 'horizontal'
                        histaxes  = False
                    elif colname == 'peaksep_kms':
                        linetype  = 'AV18_peaksep_wHorizontal'
                        xvalues   = 0.5 * xvalues
                        xlabel    = '1/2 $\\times$ '+xlabel
                        xrange    = [50,650]
                        yrange    = [-30,950]
                        # xerr      = np.sqrt(xerr) # errors are variances ???????????????????????????????????????

                        # appending collection of Lya velocity offsets from literature
                        goodent_lit  = np.where(np.isfinite(DvLyaLit_dat['peaksep_lya']))[0]
                        dv_lit       = DvLyaLit_dat['vshift_Lya'][goodent_lit]
                        dv_err_lit   = DvLyaLit_dat['vshifterr_Lya'][goodent_lit]
                        PS_lit       = DvLyaLit_dat['peaksep_lya'][goodent_lit]
                        PSerr_lit    = DvLyaLit_dat['peaksep_lya_err'][goodent_lit]
                        z_lit        = DvLyaLit_dat['redshift'][goodent_lit]
                        ids_lit      = np.asarray([990000000000]*len(z_lit))   # 990000000000 gives dot as symbol

                        xvalues      = np.append(xvalues,0.5 * PS_lit)
                        xerr         = np.append(xerr,PSerr_lit)
                        plot_dv      = np.append(plot_dv,dv_lit)
                        yerr         = np.append(yerr,dv_err_lit)
                        plot_cdatvec = np.append(plot_cdatvec,z_lit)
                        plot_ids     = np.append(plot_ids,ids_lit)

                        if len(xvalues) != len(plot_dv):
                            pdb.set_trace()
                        ylabel = ylabel.replace(' - CIII','')

                        #peaksep relation: 1.05 * x - 12
                        P1  = np.array([0.,-12.]) # lower  'endpoint' of correlation
                        P2  = np.array([1000.,1038.]) # higher 'endpoint' of correlation

                    elif colname == 'lyafwhm_kms':
                        linetype  = 'AV18_fwhm_wHorizontal'
                        # xerr      = np.sqrt(xerr) # errors are variances ???????????????????????????????????????
                        xrange    = [50,650]
                        yrange    = [-200,740]

                        # appending collection of Lya velocity offsets from literature
                        goodent_lit  = np.where(np.isfinite(DvLyaLit_dat['FWHM_lya']))[0]
                        dv_lit       = DvLyaLit_dat['vshift_Lya'][goodent_lit]
                        dv_err_lit   = DvLyaLit_dat['vshifterr_Lya'][goodent_lit]
                        FWHM_lit     = DvLyaLit_dat['FWHM_lya'][goodent_lit]
                        FWHMerr_lit  = DvLyaLit_dat['FWHM_lyaerr'][goodent_lit]
                        z_lit        = DvLyaLit_dat['redshift'][goodent_lit]
                        ids_lit      = np.asarray([990000000000]*len(z_lit))   # 990000000000 gives dot as symbol

                        xvalues      = np.append(xvalues,FWHM_lit)
                        xerr         = np.append(xerr,FWHMerr_lit)
                        plot_dv      = np.append(plot_dv,dv_lit)
                        yerr         = np.append(yerr,dv_err_lit)
                        plot_cdatvec = np.append(plot_cdatvec,z_lit)
                        plot_ids     = np.append(plot_ids,ids_lit)

                        if len(xvalues) != len(plot_dv):
                            pdb.set_trace()
                        ylabel = ylabel.replace(' - CIII','')

                        # fwhm relation: 0.9 * x - 36
                        P1  = np.array([0.,-36.]) # lower  'endpoint' of correlation
                        P2  = np.array([1000.,864.]) # higher 'endpoint' of correlation

                    elif 'absmagUV' in colname:
                        linetype = 'CM18_wHorizontal'
                        #xrange    = [-16.7,-24.4]
                        xrange    = [-24.0,-17.1]
                        yrange    = [-400,850]

                        # appending collection of Lya velocity offsets from literature
                        goodent_lit  = np.where(np.isfinite(DvLyaLit_dat['vshift_Lya']))[0]
                        dv_lit       = DvLyaLit_dat['vshift_Lya'][goodent_lit]
                        dv_err_lit   = DvLyaLit_dat['vshifterr_Lya'][goodent_lit]
                        MUV_lit      = DvLyaLit_dat['magabsUV'][goodent_lit]
                        MUVerr_lit   = DvLyaLit_dat['magabsUVerr'][goodent_lit]
                        z_lit        = DvLyaLit_dat['redshift'][goodent_lit]
                        #ids_lit      = np.asarray([7023430418]*len(z_lit))   # Erb+10 objects from lit collection to get square
                        ids_lit      = np.asarray([990000000000]*len(z_lit))   # 990000000000 gives dot as symbol

                        xvalues      = np.append(xvalues,MUV_lit)
                        xerr         = np.append(xerr,MUVerr_lit)
                        plot_dv      = np.append(plot_dv,dv_lit)
                        yerr         = np.append(yerr,dv_err_lit)
                        plot_cdatvec = np.append(plot_cdatvec,z_lit)
                        plot_ids     = np.append(plot_ids,ids_lit)

                        if len(xvalues) != len(plot_dv):
                            pdb.set_trace()
                        xlabel = 'M(UV)'
                        ylabel = ylabel.replace(' - CIII','')

                        removeerrbars = True
                        if removeerrbars:
                            print(' >> Removing errorbars: xerr (mean,median)='+
                                  str((np.mean(xerr[np.isfinite(xerr)]),np.median(np.median(xerr[np.isfinite(xerr)])))))
                            print(' >> Removing errorbars: yerr (mean,median)='+
                                  str((np.mean(yerr[np.isfinite(yerr)]),np.median(np.median(yerr[np.isfinite(yerr)])))))
                            xerr[np.abs(xerr) != 99] = 0.0
                            yerr[np.abs(yerr) != 99] = 0.0

                    else:
                        linetype='horizontal'

                    if verbose: print('\n - Estimating correlation coeffecients for "'+colname+'":')
                    r_pearson, pvalue_pearson   = ss.pearsonr(xvalues[np.isfinite(xvalues) & np.isfinite(xvalues)],
                                                              plot_dv[np.isfinite(xvalues) & np.isfinite(xvalues)])
                    r_spearman, pvalue_spearman = ss.spearmanr(xvalues[np.isfinite(xvalues) & np.isfinite(xvalues)],
                                                               plot_dv[np.isfinite(xvalues) & np.isfinite(xvalues)])
                    if verbose: print('   r_pearson, pvalue_pearson   = '+str("%.4f" % r_pearson)+', '+str("%.4f" % pvalue_pearson))
                    if verbose: print('   r_spearman, pvalue_spearman = '+str("%.4f" % r_spearman)+', '+str("%.4f" % pvalue_spearman)+'\n')

                    if ('lyafwhm_kms' in colname) or (colname == 'peaksep_kms'):
                        if verbose: print(' - Calculating the perpendicular euclidian distance '
                                          'to the best-fit line from Verhamme+18 for "'+colname+'":')
                        dists_perp = np.zeros(len(xvalues))
                        for dd, dpoint_x in enumerate(xvalues):
                            datapoint      = np.array([dpoint_x,plot_dv[dd]])
                            dists_perp[dd] = np.abs(np.cross(P2-P1, P1-datapoint)) / np.linalg.norm(P2-P1)
                        if verbose: print('   The distance array [km/s] looks like:\n'+str(dists_perp))
                        if verbose: print('   The stats for those distances are:')
                        if verbose: print('    mean    = '+str(np.mean(dists_perp)))
                        if verbose: print('    median  = '+str(np.median(dists_perp)))
                        if verbose: print('    std     = '+str(np.std(dists_perp)))

                    # yerr = None
                    # xerr = None
                    Nhistbins = np.ceil(np.sqrt(len(xvalues)))
                    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,plot_dv,xerr,yerr,xlabel,ylabel,
                                                                   'dummydat',linetype=linetype,title=None, #'this is title',
                                                                   ids=plot_ids,
                                                                   ylog=False,xlog=False,yrange=yrange,xrange=xrange,
                                                                   colortype=colortype,colorcode=True,
                                                                   cdatvec=plot_cdatvec,
                                                                   point_text=None, #dat_uves['id'][dvs_ent[vv]].astype(str),
                                                                   photoionizationplotparam=None,
                                                                   histaxes=histaxes,Nbins=Nhistbins, showgraylimits=True,
                                                                   overwrite=overwrite,verbose=verbose)


                    if plotdvversion is not None:
                        zsyscol      = infodat[infocols[plotdvversion][0]][dvs_ent[vv]]
                        zsyscol_err  = infodat[infocols[plotdvversion][1]][dvs_ent[vv]]
                        xvalues      = 299792.458 * (zleadline-zsyscol)/(1.0+zsyscol) # cf. Erb+2014
                        xerr         = np.abs(xvalues) * np.sqrt(2*(zsyscol_err/zsyscol)**2)

                        xlabel    = infocols[plotdvversion][2].replace('$z_\\textrm{sys}$','$\Delta v$(Ly$\\alpha$) from ')
                        Nhistbins = np.ceil(np.sqrt(len(plot_dv_prelit)))
                        uves.plot_mocspecFELISresults_summary_plotcmds(plotname.replace('.pdf','_dv.pdf'),
                                                                       xvalues,plot_dv_prelit,xerr,plot_dv_err_prelit,xlabel,ylabel,
                                                                       'dummydat',linetype='onetoone',title=None, #'this is title',
                                                                       ids=plot_ids,
                                                                       ylog=False,xlog=False,yrange=yrange,xrange=[-1000,1000],
                                                                       colortype=colortype,colorcode=True,
                                                                       cdatvec=plot_cdatvec,
                                                                       point_text=None, #dat_uves['id'][dvs_ent[vv]].astype(str),
                                                                       photoionizationplotparam=None,
                                                                       histaxes=histaxes,Nbins=Nhistbins, showgraylimits=True,
                                                                       overwrite=overwrite,verbose=verbose)
                        print('--------'+plotdvversion+'--------')
                        print(dat_uves['id'][dvs_ent[vv]].astype(str))
                        print(plot_dv_prelit)
                        print(plot_dv_err_prelit)
                        print(zsyscol)
                        print(zsys)


            plotname  = outputdir+'evaluate_voffsets_'+leadline[vv].replace(' ','')+'-'+linename[vv]+'_FWHMvsPeakSep.pdf'
            Nhistbins = np.ceil(np.sqrt(len(val_peaksep_kms)))
            uves.plot_mocspecFELISresults_summary_plotcmds(plotname,val_peaksep_kms,val_lyafwhm_kms,val_peaksep_kms_err,val_lyafwhm_kms_err,
                                                           infocols['peaksep_kms'][2],infocols['lyafwhm_kms'][2],
                                                           'dummydat',linetype='onetoone',title=None, #'this is title',
                                                           ids=plot_ids,
                                                           ylog=False,xlog=False,yrange=None,xrange=None,
                                                           colortype='vshift',colorcode=True,
                                                           cdatvec=plot_dv_prelit,
                                                           point_text=dat_uves['id'][dvs_ent[vv]].astype(str),
                                                           photoionizationplotparam=None,
                                                           histaxes=histaxes,Nbins=Nhistbins, showgraylimits=True,
                                                           overwrite=overwrite,verbose=verbose)

            plotname  = outputdir+'evaluate_voffsets_'+leadline[vv].replace(' ','')+'-'+linename[vv]+'_dvFWHMvsPeakSep.pdf'
            Nhistbins = np.ceil(np.sqrt(len(val_peaksep_kms)))
            uves.plot_mocspecFELISresults_summary_plotcmds(plotname,
                                                           1.05*val_peaksep_kms/2.-12,
                                                           0.9*val_lyafwhm_kms-34,
                                                           None,None,
                                                           'dv from '+infocols['peaksep_kms'][2],'dv from '+infocols['lyafwhm_kms'][2],
                                                           'dummydat',linetype='onetoone',title=None, #'this is title',
                                                           ids=plot_ids,
                                                           ylog=False,xlog=False,yrange=None,xrange=None,
                                                           colortype='vshift',colorcode=True,
                                                           cdatvec=plot_dv_prelit,
                                                           point_text=dat_uves['id'][dvs_ent[vv]].astype(str),
                                                           photoionizationplotparam=None,
                                                           histaxes=histaxes,Nbins=Nhistbins, showgraylimits=True,
                                                           overwrite=overwrite,verbose=verbose)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def evaluate_FELISlinewidths(linefluxcatalog,infofile,outputdir='./linewidthFigures/',overwrite=False,verbose=True):
    """
    Function to evaluate and line widths from FELIS fits. Based on evluate_velocityoffsets

    --- Example of use ---

    """
    agn, agncand = uves.get_AGN_ids()
    dat_uves = afits.open(linefluxcatalog)[1].data
    infodat  = afits.open(infofile)[1].data
    infodat  = infodat[np.where((infodat['id']<4.9e8) | (infodat['id']>5.9e8))[0]] # ignoring UDF MW mock ids

    for ii, id in enumerate(dat_uves['id']): # double checking that order of objects is the same
        if infodat['id'][ii] != id:
            sys.exit('There was a mismatch in ids for entry '+str(ii))

    linenamelist = ['NV','CIV','HeII','OIII','SiIII','CIII','MgII']
    Nlines       = len(linenamelist)
    # zlineranges  = [[2.8918,6.4432],[2.1114,4.9699],[1.9379,4.6411],[1.8969,4.5632],[1.5514,3.9067],[1.5241,3.8548],[0.7174,2.3142]]
    zlineranges  = [[2.86,6.5],[2.0,5.1],[1.7,4.8],[1.7,4.7],[1.4,3.95],[1.4,3.95],[0.7174,2.5]]
    #\lya                		&  2.9729 -- 6.5958 	&

    # - - - - - - - ALL lead lines - - - - - - -
    sig_NV_ent    = np.where((dat_uves['sigma_NV'] != 99) & (dat_uves['duplicationID'] == 0) &
                          np.isfinite(dat_uves['sigma_NV']))[0]
    sig_NV        = dat_uves['sigma_NV'][sig_NV_ent]

    sig_CIV_ent   = np.where((dat_uves['sigma_CIV'] != 99) & (dat_uves['duplicationID'] == 0) &
                          np.isfinite(dat_uves['sigma_CIV']))[0]
    sig_CIV       = dat_uves['sigma_CIV'][sig_CIV_ent]

    sig_HeII_ent  = np.where((dat_uves['sigma_HeII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['sigma_HeII']))[0]
    sig_HeII      = dat_uves['sigma_HeII'][sig_HeII_ent]

    sig_OIII_ent  = np.where((dat_uves['sigma_OIII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['sigma_OIII']))[0]
    sig_OIII      = dat_uves['sigma_OIII'][sig_OIII_ent]

    sig_SiIII_ent = np.where((dat_uves['sigma_SiIII'] != 99) & (dat_uves['duplicationID'] == 0) &
                            np.isfinite(dat_uves['sigma_SiIII']))[0]
    sig_SiIII     = dat_uves['sigma_SiIII'][sig_SiIII_ent]

    sig_CIII_ent  = np.where((dat_uves['sigma_CIII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['sigma_CIII']))[0]
    sig_CIII      = dat_uves['sigma_CIII'][sig_CIII_ent]

    sig_MgII_ent  = np.where((dat_uves['sigma_MgII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['sigma_MgII']))[0]
    sig_MgII      = dat_uves['sigma_MgII'][sig_MgII_ent]

    # - - - - - - - Only Lya as lead line - - - - - - -
    sig_NV_ent_LAE    = np.where((dat_uves['sigma_NV'] != 99) & (dat_uves['duplicationID'] == 0) &
                          np.isfinite(dat_uves['sigma_NV']) & (infodat['leadline'] == 'Lya'))[0]
    sig_NV_LAE        = dat_uves['sigma_NV'][sig_NV_ent_LAE]

    sig_CIV_ent_LAE   = np.where((dat_uves['sigma_CIV'] != 99) & (dat_uves['duplicationID'] == 0) &
                          np.isfinite(dat_uves['sigma_CIV']) & (infodat['leadline'] == 'Lya'))[0]
    sig_CIV_LAE       = dat_uves['sigma_CIV'][sig_CIV_ent_LAE]

    sig_HeII_ent_LAE  = np.where((dat_uves['sigma_HeII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['sigma_HeII']) & (infodat['leadline'] == 'Lya'))[0]
    sig_HeII_LAE      = dat_uves['sigma_HeII'][sig_HeII_ent_LAE]

    sig_OIII_ent_LAE  = np.where((dat_uves['sigma_OIII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['sigma_OIII']) & (infodat['leadline'] == 'Lya'))[0]
    sig_OIII_LAE      = dat_uves['sigma_OIII'][sig_OIII_ent_LAE]

    sig_SiIII_ent_LAE = np.where((dat_uves['sigma_SiIII'] != 99) & (dat_uves['duplicationID'] == 0) &
                            np.isfinite(dat_uves['sigma_SiIII']) & (infodat['leadline'] == 'Lya'))[0]
    sig_SiIII_LAE     = dat_uves['sigma_SiIII'][sig_SiIII_ent_LAE]

    sig_CIII_ent_LAE  = np.where((dat_uves['sigma_CIII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['sigma_CIII']) & (infodat['leadline'] == 'Lya'))[0]
    sig_CIII_LAE      = dat_uves['sigma_CIII'][sig_CIII_ent_LAE]

    sig_MgII_ent_LAE  = np.where((dat_uves['sigma_MgII'] != 99) & (dat_uves['duplicationID'] == 0) &
                           np.isfinite(dat_uves['sigma_MgII']) & (infodat['leadline'] == 'Lya'))[0]
    sig_MgII_LAE      = dat_uves['sigma_MgII'][sig_MgII_ent_LAE]

    #------------------------------------------------------------------------------
    if verbose: print(' - Setting up and generating plots of lead line widths ')

    linename   = linenamelist+linenamelist
    zranges    = zlineranges+[[2.0,3.95]]+zlineranges
    leadline   = ['lead line']*Nlines + ['Lya']*Nlines
    sigs_ent    = [sig_NV_ent,sig_CIV_ent,sig_HeII_ent,sig_OIII_ent,sig_SiIII_ent,sig_CIII_ent,sig_MgII_ent]+\
                 [sig_NV_ent_LAE,sig_CIV_ent_LAE,sig_HeII_ent_LAE,sig_OIII_ent_LAE,sig_SiIII_ent_LAE,sig_CIII_ent_LAE]#,sig_MgII_ent_LAE]
    sigs        = [sig_NV,sig_CIV,sig_HeII,sig_OIII,sig_SiIII,sig_CIII,sig_MgII]+\
                 [sig_NV_LAE,sig_CIV_LAE,sig_HeII_LAE,sig_OIII_LAE,sig_SiIII_LAE,sig_CIII_LAE]#,sig_MgII_LAE]

    for ss, sig in enumerate(sigs):
        # if not 'CIII' in linename[ss]: continue
        objids    = dat_uves['id'][sigs_ent[ss]]
        zleadline = dat_uves['redshift'][sigs_ent[ss]]

        histaxes  = False
        #Nhistbins = 50
        yrange    = [0,1.3]

        if leadline[ss] == 'Lya':
            llstring = 'Ly$\\alpha$'
        elif leadline[ss] == 'CIII':
            llstring = 'CIII'
        else:
            llstring = leadline[ss]
        ylabel    = linename[ss]+' FELIS line template width  [\AA]'

        #---- vs redshift ---
        plotname = outputdir+'evaluate_FELISlinewidths_'+leadline[ss].replace(' ','')+'-'+linename[ss]+'VSredshiftLeadLine.pdf'

        xrange = zranges[ss]
        if leadline[ss] == 'Lya':
            xrange[0] = np.max([xrange[0],2.8])
            linetype  = 'horizontal'
        else:
            xrange[0] = np.max([xrange[0],1.4])
            linetype  = 'horizontalWlya'

        xlabel   = '$z$('+llstring+')'
        xerr     = None
        yerr     = None

        Nhistbins = np.ceil(np.sqrt(len(dat_uves['redshift'][sigs_ent[ss]])))
        uves.plot_mocspecFELISresults_summary_plotcmds(plotname,dat_uves['redshift'][sigs_ent[ss]],sig,xerr,yerr,xlabel,ylabel,
                                                       'dummydat',linetype=linetype,title=None, #'this is title',
                                                       ids=dat_uves['id'][sigs_ent[ss]],
                                                       ylog=False,xlog=False,yrange=yrange,xrange=xrange,
                                                       colortype='s2nfelis',colorcode=True,
                                                       cdatvec=dat_uves['s2n_'+linename[ss]][sigs_ent[ss]],
                                                       point_text=None, #dat_uves['id'][sigs_ent[ss]].astype(str),
                                                       photoionizationplotparam=None,
                                                       histaxes=histaxes,Nbins=Nhistbins, showgraylimits=True,
                                                       overwrite=overwrite,verbose=verbose)

        #---- vs LAE parameters ---
        print(leadline[ss])
        if leadline[ss] == 'Lya':
            infocols = uves.get_infodat_plotcols()

            for cc, colname in enumerate(infocols.keys()):
                xlabel   = infocols[colname][2]
                xrange   = None
                plotname = outputdir+'evaluate_FELISlinewidths_'+leadline[ss].replace(' ','')+'-'+linename[ss]+'VS'+colname+'.pdf'

                yerr = None

                xvalues  = infodat[infocols[colname][0]][sigs_ent[ss]]
                if len(xvalues[np.isfinite(xvalues)]) > 0:
                    if infocols[colname][1] is None:
                        xerr     = np.asarray([np.nan]*len(sigs_ent[ss]))
                    else:
                        xerr     = infodat[infocols[colname][1]][sigs_ent[ss]]

                    plot_sig      = sig
                    plot_ids     = dat_uves['id'][sigs_ent[ss]]
                    plot_cdatvec = zleadline

                    plot_sig_prelit      = plot_sig
                    plot_sig_err_prelit  = yerr

                    colortype = 'redshift'
                    linetype='horizontal'

                    # yerr = None
                    # xerr = None
                    Nhistbins = np.ceil(np.sqrt(len(xvalues)))
                    uves.plot_mocspecFELISresults_summary_plotcmds(plotname,xvalues,plot_sig,xerr,yerr,xlabel,ylabel,
                                                                   'dummydat',linetype=linetype,title=None, #'this is title',
                                                                   ids=plot_ids,
                                                                   ylog=False,xlog=False,yrange=yrange,xrange=xrange,
                                                                   colortype=colortype,colorcode=True,
                                                                   cdatvec=plot_cdatvec,
                                                                   point_text=None, #dat_uves['id'][sigs_ent[ss]].astype(str),
                                                                   photoionizationplotparam=None,
                                                                   histaxes=histaxes,Nbins=Nhistbins, showgraylimits=True,
                                                                   overwrite=overwrite,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_param_for_photoionizationmodels(linefluxcatalog,outdir,Nsigma=1,infofile=None,generatePDFplots=False,
                                        addinputobj=True,addliteratureobj=False,verbose=True):
    """
    Function assembling and generating photoionization model "PDFs" of UVES objects

    --- Example of use ---
    import uvEmissionlineSearch as uves
    uvesdir          = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'
    outdir           = uvesdir+'NEOGALpdffigures/'
    linefluxcatalog  = uvesdir+'back2backAnalysis_200213/results_master_catalog_version200213.fits'
    # infofile         = uvesdir+'objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile         = uvesdir+'objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    uves.get_param_for_photoionizationmodels(linefluxcatalog,outdir,Nsigma=1,infofile=None,generatePDFplots=False,addliteratureobj=False,verbose=True)

    """
    FRdiclist = []

    coldic = {}
    coldic['NVCIV']     = 'NV1240/CIV1550'
    coldic['NVHeII']    = 'NV1240/HeII1640'
    coldic['NVOIII']    = 'NV1240/OIII1663'
    coldic['NVSiIII']   = 'NV1240/SiIII1888'
    coldic['NVCIII']    = 'NV1240/CIII1908'
    # coldic['NVMgII']    = 'NV1240/SiIII1888'  # MgII not in NEOGAL models

    coldic['CIVHeII']   = 'CIV1550/HeII1640'
    coldic['CIVOIII']   = 'CIV1550/OIII1663'
    coldic['CIVSiIII']  = 'CIV1550/SiIII1888'
    coldic['CIVCIII']   = 'CIV1550/CIII1908'
    # coldic['CIVMgII']   = 'CIV1550/'] # MgII not in NEOGAL models
    coldic['HeIIOIII']  = 'HeII1640/OIII1663'
    coldic['HeIISiIII'] = 'HeII1640/SiIII1888'
    coldic['HeIICIII']  = 'HeII1640/CIII1908'
    # coldic['HeIIMgII']  = 'HeII1640/' # not in NEOGAL models
    coldic['OIIISiIII'] = 'OIII1663/SiIII1888'
    coldic['OIIICIII']  = 'OIII1663/CIII1908'
    # coldic['OIIIMgII']  = 'OIII1663/' # not in NEOGAL models
    coldic['SiIIICIII'] = 'SiIII1888/CIII1908'
    # coldic['SiIIIMgII'] = 'SiIII1888/' # not in NEOGAL models
    # coldic['CIIIMgII']  = 'CIII1908/' # not in NEOGAL models

    coldic['HeIICIV']   = 'HeII1640/CIV1550'
    coldic['OIIICIV']   = 'OIII1663/CIV1550'
    coldic['SiIIICIV']  = 'SiIII1888/CIV1550'
    coldic['CIIICIV']   = 'CIII1908/CIV1550'
    coldic['OIIIHeII']  = 'OIII1663/HeII1640'
    coldic['SiIIIHeII'] = 'SiIII1888/HeII1640'
    coldic['CIIIHeII']  = 'CIII1908/HeII1640'
    coldic['SiIIIOIII'] = 'SiIII1888/OIII1663'
    coldic['CIIIOIII']  = 'CIII1908/OIII1663'
    coldic['CIIISiIII'] = 'CIII1908/SiIII1888'

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if addinputobj:
        obsdat    = afits.open(linefluxcatalog)[1].data
        Ndatobj   = len(obsdat['id'])
        if verbose: print(' - Defining parameter ranges for '+str(Ndatobj)+' objects in input ')
        for ii, objid in enumerate(obsdat['id']):
            if verbose:
                infostr = '   Getting info for '+str(objid)+' ('+str("%.5d" % (ii+1))+' / '+str("%.5d" % Ndatobj)+')        '
                sys.stdout.write("%s\r" % infostr)
                sys.stdout.flush()

            FRdic = {}
            FRdic['id'] = int(objid)

            for colname in obsdat.dtype.names:
                if colname in ['FR_'+key for key in coldic.keys()]:
                    frval  = obsdat[colname][ii]
                    errval = obsdat[colname.replace('FR_','FRerr_')][ii]

                    if np.isfinite(frval):
                        if errval == 99:
                            valrange = [0,frval]
                        elif errval == -99:
                            valrange = [frval,1e35]
                        else:
                            valrange = [frval-errval*Nsigma,frval+errval*Nsigma]
                    else:
                        valrange = None

                    if valrange is not None:
                        FRdic[coldic[colname.split('R_')[-1]]] = valrange

            if len(FRdic.keys()) > 1: # only add objects where there are constraints on line ratios from data
                FRdiclist.append(FRdic)
        if verbose: print('\n   done')
        Nobjdatadded = len(FRdiclist)
        if verbose: print(' - Found '+str(Nobjdatadded)+' objects with constraints from data ')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    basename = outdir+'UVESvsPhotoionizationModelParams'
    if addliteratureobj:
        basename = basename+'_wLit'
        litdir   = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/'
        litcat   = litdir+'literaturecollection_emissionlinestrengths.fits'
        litdat   = afits.open(litcat)[1].data
        Nlitobj  = len(litdat['id'])
        if verbose: print(' - Defining parameter ranges for '+str(Nlitobj)+' objects in input ')
        for ii, objid in enumerate(litdat['id']):
            if verbose:
                infostr = '   Getting info for '+str(objid)+' ('+str("%.5d" % (ii+1))+' / '+str("%.5d" % Nlitobj)+')        '
                sys.stdout.write("%s\r" % infostr)
                sys.stdout.flush()

            FRdic = {}
            FRdic['id'] = int(objid)

            for colname in litdat.dtype.names:
                if colname in ['FR_'+key for key in coldic.keys()]:
                    frval  = litdat[colname][ii]
                    errval = litdat[colname.replace('FR_','FRerr_')][ii]

                    if np.isfinite(frval):
                        if errval == 99:
                            valrange = [0,frval]
                        elif errval == -99:
                            valrange = [frval,1e35]
                        else:
                            valrange = [frval-errval*Nsigma,frval+errval*Nsigma]
                    else:
                        valrange = None

                    if valrange is not None:
                        FRdic[coldic[colname.split('R_')[-1]]] = valrange

            if len(FRdic.keys()) > 1: # only add objects where there are constraints on line ratios from data
                FRdiclist.append(FRdic)

        if verbose: print('\n   done')
        Nobjlitadded = len(FRdiclist)
        if addinputobj:
             Nobjlitadded = Nobjlitadded - Nobjdatadded
        if verbose: print(' - Found '+str(Nobjlitadded)+' objects with constraints from literature ')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # colors pulled from uves.linecolors('viridis')
    col_NEOGAL_AGN   = (0.267004, 0.004874, 0.329415, 1.0)    # purple
    col_NEOGAL_SF    = (0.190631, 0.407061, 0.556089, 1.0)    # blue
    col_BPASS_bin    = (0.20803, 0.718701, 0.472873, 1.0)     # green
    col_BPASS_sin    = (0.993248, 0.906157, 0.143936, 1.0)    # yellow

    # colors pulled from uves.linecolors('plasma')
    # col_NEOGAL_AGN   = (0.050383, 0.029803, 0.527975, 1.0)    # purple
    # col_NEOGAL_SF    = (0.610667, 0.090204, 0.619951, 1.0)    # blue
    # col_BPASS_bin    = (0.928329, 0.472975, 0.326067, 1.0)    # red
    # col_BPASS_sin    = (0.940015, 0.975158, 0.131326, 1.0)    # yellow

    paramcollections, collectionstats = pp.estimate_object_PDFs(FRdiclist, basename=basename, generatePDFplots=generatePDFplots,
                                                                col_NEOGAL_AGN=col_NEOGAL_AGN,col_NEOGAL_SF=col_NEOGAL_SF,
                                                                col_BPASS_bin=col_BPASS_bin,col_BPASS_sin=col_BPASS_sin,
                                                                maxPDFyscale=False,verbose=verbose,showemptyparamCorner=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_AGN_ids(infofile=None):
    """
    Returning list of IDs of secure and likely AGN

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves
    uvesdir    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'
    # infofile   = uvesdir+'objectinfofile_zGT1p5_3timesUDFcats_JKthesisInfo.fits'
    infofile   = uvesdir+'objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    agn, agncand = uves.get_AGN_ids(infofile=infofile)

    """
    ### From TU on 171020
    # First 24 fields bona fide high-z AGN:
    # 104014050
    # 115003085
    # 214002011
    #
    # AGN candidates:
    # 123048186 (X-ray probably too far away)
    # 123051191 (X-ray probably too far away)
    # 121033078 (because of CIV so strong, but don't know)
    agn     = [104014050,115003085,214002011]
    agncand = [123048186,123051191,121033078]

    # AGN from TUs 171020 Xray counterpart list
    # 102007068  0.338   304  0.20  2.982e-16  0.340     zSpec  AGN
    # 102008071  0.338   312  0.86  6.491e-17  0.336     zSpec  AGN
    # 102031144  0.665   290  0.20  3.498e-16  0.664     zSpec  AGN
    # 102037154  1.412   287  0.51  5.925e-17  1.413     zSpec  AGN
    # 103022086  0.670   322  0.24  1.510e-16  0.671     zSpec  AGN
    # 104014050  3.662   337  0.17  1.269e-15  3.660     zSpec  AGN
    # 106036089  0.904   344  0.82  6.638e-17  0.956     S14    AGN
    # 106048103  0.665   340  0.15  2.406e-15  0.666     zSpec  AGN
    # 108025145  0.736   407  0.49  7.387e-17  0.736     zSpec  AGN
    # 109030090  1.044   447  0.17  5.737e-16  1.043     zSpec  AGN
    # 111004005  0.604   367  0.28  4.594e-15  0.604     zSpec  AGN
    # 113001007  0.232   508  0.94  2.411e-17  0.220     S14    AGN
    # 113010038  0.577   460  0.35  2.778e-16  0.577     zSpec  AGN
    # 114024110  1.035   443  0.26  8.709e-16  1.036     zSpec  AGN
    # 114028115  1.098   509  0.30  1.525e-16  1.097     zSpec  AGN
    # 115003085  3.710   551  0.11  2.158e-15  3.700     zSpec  AGN
    # 116003060  1.364   634  0.14  3.899e-17  1.363     H14    AGN
    # 117034085  0.228   693  1.13  6.435e-17  2.302     zSpec  AGN
    # 119034073  1.015   814  0.15  4.348e-15  1.016     zSpec  AGN
    # 120023032  1.118   861  0.25  1.744e-16  1.120     zSpec  AGN
    # 123005089  0.544   640  0.51  4.219e-17  0.552     H14    AGN
    # 123051191  4.510   625  1.33  2.102e-17  2.616     H14    AGN
    agn     = agn+[102007068,102008071,102031144,102037154,103022086,104014050,106036089,106048103,
                   108025145,109030090,111004005,113001007,113010038,114024110,114028115,115003085,
                   116003060,117034085,119034073,120023032,123005089,123051191]

    # Table 5 Urrutia+18
    # 102028132
    # 105027078
    # 106036089
    # 113001007
    # 113022070
    # 116003060
    # 117034085
    # 118011046
    # 123005089
    # 123048186
    # 123051191
    # 124037072
    # 137003006
    # 139073330
    # 143041126
    # 146002220

    agn     = agn+[102028132,105027078,106036089,113001007,113022070,116003060,
                   117034085,118011046,123005089,123048186,123051191,124037072,
                   137003006,139073330,143041126,146002220]

    ### From FELIS UV search
    agn     = agn+[157001017,221004004,601381485]
    agncand = agncand+[136003119,600691153]

    # - - - - - - - - - - - - - - - - - - - - - - - - -
    agn     = np.unique(np.asarray(agn))
    agncand = np.unique(np.asarray(agncand))

    ignoreids = []#[123048186,123051191]
    for iobj in ignoreids:
        ent_agn = np.where(agn == iobj)[0]
        if len(ent_agn) > 0:
            agn[ent_agn] = -99

        ent_agncand = np.where(agncand == iobj)[0]
        if len(ent_agncand) > 0:
            agn[ent_agncand] = -99

    agn     = agn[agn != -99]
    agncand = agncand[agncand != -99]
    # - - - - - - - - - - - - - - - - - - - - - - - - -

    if infofile is not None:
        infodat    = afits.open(infofile)[1].data
        print('# id ra dec redshift duplicationID ')
        print('#------ AGN ------ ')
        for ii, objid in enumerate(agn):
            objent = np.where(infodat['id'] == objid)[0]
            if len(objent) == 1:
                print(str(objid)+' '+str(infodat['ra'][objent][0])+' '+str(infodat['dec'][objent][0])+' '+
                      str(infodat['redshift'][objent][0])+' '+str(infodat['duplicationID'][objent][0]))
            elif len(objent) == 0:
                print(str(objid)+' notUVESobj notUVESobj notUVESobj notUVESobj')
            else:
                sys.exit('Weird - more than 1 match to id'+str(objid))

        print('#------ AGN candidates ------ ')
        for ii, objid in enumerate(agncand):
            objent = np.where(infodat['id'] == objid)[0]
            if len(objent) == 1:
                print(str(objid)+' '+str(infodat['ra'][objent][0])+' '+str(infodat['dec'][objent][0])+' '+
                      str(infodat['redshift'][objent][0])+' '+str(infodat['duplicationID'][objent][0]))
            elif len(objent) == 0:
                print(str(objid)+' noinfo noinfo noinfo noinfo')
            else:
                sys.exit('Weird - more than 1 match to id'+str(objid))

    return agn, agncand
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def linecolors(colormap=None):
    """
    Return colors of lines to use in plots
    """
    linecolors  = collections.OrderedDict()

    linecolors['ANY']       = 'black'
    #linecolors['CIIIorCIV'] = 'gray'
    linecolors['NV']        = 'purple'
    linecolors['CIV']       = 'blue'
    linecolors['HeII']      = 'cyan'
    linecolors['OIII']      = 'green'
    linecolors['SiIII']     = 'orange'
    linecolors['CIII']      = 'red'
    linecolors['MgII']      = 'darkred'

    if colormap is not None:
        cmin      = 2
        cmax      = len(linecolors.keys())-1
        cmap      = matplotlib.cm.get_cmap(colormap)
        colnorm   = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
        cmaparr   = np.linspace(cmin, cmax, num=50)
        m         = plt.cm.ScalarMappable(cmap=cmap)
        m.set_array(cmaparr)

        for ll, uvline in enumerate(linecolors.keys()):
            if uvline not in ['ANY','CIIIorCIV']:
                linecolors[uvline] = cmap(colnorm(ll))

    return linecolors
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def linenames():
    """
    Return latex names of lines to use in plots
    """
    linenames   = {'ANY'       :['Any detection', 'Any UV detection'],
                   'CIIIorCIV' :['CIII or CIV',   'TBD'],
                   'NV'        :['NV',            'TBD'],
                   'CIV'       :['CIV',           'TBD'],
                   'HeII'      :['HeII',          'TBD'],
                   'OIII'      :['OIII',          'TBD'],
                   'SiIII'     :['SiIII',         'TBD'],
                   'CIII'      :['CIII',          'TBD'],
                   'MgII'      :['MgII',          'TBD']}
    return linenames

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def mapp2mabs(mag_input,redshift,mapp2Mabs=True,verbose=True):
    """
    Converting apparent magnitudes into absolute magnitudes

    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves

    redshift=[2.2838,2.3135,2.2837,2.2776,2.2928,2.2930,2.2949,2.2742,2.2935,2.3100,2.3002,2.3254,2.3064,2.2942,2.3064,2.2903,2.3139,3.0690,3.0965,3.0788,3.0691,3.0692,3.0687,3.0902,3.0919,3.0631,3.0670,3.0666,3.0975,3.0551,3.0978,3.0645,3.1013,3.0845,3.0870,3.0873]
    mag_input=[24.48,24.75,23.94,25.24,23.09,25.07,24.84,26.01,26.07,26.11,26.14,26.85,25.71,26.48,26.39,26.85,26.72,23.92,24.42,24.34,27.00,24.87,25.84,24.75,25.98,25.82,27.00,25.50,27.00,26.53,26.55,26.64,26.40,27.00,25.95,27.00]
    mabs = uves.mapp2mabs(mag_input,redshift,mapp2Mabs=True)

    """
    mag_input  = np.asarray(mag_input)
    redshift   = np.asarray(redshift)
    Nobj       = len(redshift)
    if len(mag_input) != Nobj:
        sys.exit('The number of magnitudes and redshifts provide do not match ('+str(len(mag_input))+' vs. '+str(Nobj)+')')
    # uvesdir    = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/'
    # infofile   = uvesdir+'objectinfofile_zGT1p5_3timesUDFcats_JK100fieldinfocat.fits'
    # infodat  = afits.open(infofile)[1].data
    # infodat  = infodat[np.where((infodat['id']<4.9e8) | (infodat['id']>5.9e8))[0]] # ignoring UDF MW mock ids

    cosmo    = acosmo.FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
    #cosmo    = acosmo.FlatLambdaCDM(H0=69 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.29) #Cosmology used in Erb+14
    Dlum     = cosmo.luminosity_distance(redshift).to(u.pc).value

    # if isinstance(Mapp,types.FloatType) and Av == -99: # if Av is -99, calculate it
    #     Av, Ebv = getAv(RA,DEC,band)

    Kcorrection   = (2.5 * np.log10(1.0 + redshift)) # assumes source has flat (beta = -2) SED.
                                                     #  A bluer beta will likely give you an additional
                                                     # correction of about ~0.1 mag or so.
    DM            = 5.*np.log10(Dlum)-5.             # distance modulus
    if mapp2Mabs:
        if verbose: print(' - Converting the '+str(Nobj)+' provided apparant magnitudes to absolute magnitudes (mapp2Mabs=True)')
        mag_output  = mag_input - DM + Kcorrection #- Av
    else:
        if verbose: print(' - Converting the '+str(Nobj)+' provided absolute magnitudes to apparant magnitudes (mapp2Mabs=False)')
        mag_output  = mag_input + DM - Kcorrection

    return mag_output
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plotAndFit_parameterdist(masterfits, infofile, addliterature=True, paramtype='LineOverLya', verbose=True,
                             xrange = None,
                             outdir='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/paramdistributionFigures/'):
    """
    Function to plot and fit parameter distributions

    --- INPUT ---


    --- EXAMPLE OF USE ---
    import uvEmissionlineSearch as uves


    """
    if verbose: print(' - Loading catalogs to extract data from ')
    masterdat = afits.open(masterfits)[1].data
    infodat   = afits.open(infofile)[1].data
    infodat   = infodat[np.where((infodat['id']<4.9e8) | (infodat['id']>5.9e8))[0]] # ignoring UDF MW mock ids

    for ii, id in enumerate(masterdat['id']): # double checking that order of objects is the same
        if infodat['id'][ii] != id:
            sys.exit('There was a mismatch in ids for entry '+str(ii))

    litcat  = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/' \
              'literaturecollection_emissionlinestrengths.fits'
    litdat  = afits.open(litcat)[1].data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    lines = ['NV','CIV','HeII','OIII','SiIII','CIII','MgII']
    for line in lines:
        if paramtype == 'LineOverLya':
            xlabel   = 'EW$_0$('+line+') / EW$_0$(Ly$\\alpha$) [\AA]'
            goodent  = np.where((np.abs(masterdat['ferr_'+line]) != 99) &
                                (np.abs(masterdat['EW0err_'+line]) != 99) &
                                (np.isfinite(masterdat['EW0_'+line])) &
                                (masterdat['EW0_'+line] > 0) &
                                (np.abs(infodat['EW_0_beta_own_median_error_jk100']) != 99)  &
                                (np.isfinite(infodat['EW_0_beta_own_median_jk100'])) &
                                (infodat['EW_0_beta_own_median_jk100'] > 0))[0]
            paramval = masterdat['EW0_'+line][goodent]/infodat['EW_0_beta_own_median_jk100'][goodent]

            if addliterature:
                goodent_lit  = np.where((np.abs(litdat['EW0err_'+line] != 99)) &
                                        (np.isfinite(litdat['EW0_'+line])) &
                                        (np.abs(litdat['EW0err_Lya']  != 99)) &
                                        (np.isfinite(litdat['EW0_Lya']))  )[0]
                paramval_lit = litdat['EW0_'+line][goodent_lit]/litdat['EW0_Lya'][goodent_lit]

                paramval = np.append(paramval,paramval_lit)
                Nlit     = len(paramval_lit)
                textaddition = '(incl. '+str(Nlit)+' literture values)'
            else:
                textaddition = '(no literature values incl.)'

        Nparam = len(paramval)
        if verbose: print('\n - Assembled '+str(Nparam)+' mesurements of the type '+paramtype.replace('Line',line)+' '+textaddition)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if xrange is not None:
            cutent    = np.where((paramval >= xrange[0]) & (paramval <= xrange[1]))[0]
            paramval  = paramval[cutent]
            Ninit     = Nparam
            Nparam    = len(paramval)
            if verbose: print('   xrange='+str(xrange)+' provided so limiting to the '+str(Nparam)+' values in that range')
            labeladdition = ' (of '+str(Ninit)+')'
        else:
            labeladdition = ' (all obj)'
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if len(paramval) > 2:
            plotname = outdir+'parameterdistribution_'+(paramtype.replace('Line',line))+'.pdf'
            if verbose: print(' - Setting up and generating plot')
            fig = plt.figure(figsize=(5, 5))
            fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.97, bottom=0.15, top=0.90)
            Fsize    = 15
            lthick   = 2
            marksize = 4
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif',size=Fsize)
            plt.rc('xtick', labelsize=Fsize)
            plt.rc('ytick', labelsize=Fsize)
            plt.clf()
            plt.ioff()

            xmin = np.min(paramval)
            xmax = np.max(paramval)

            Nbins   = np.ceil(np.sqrt(len(paramval)))
            binwidth_x = np.diff([xmin,xmax])/Nbins
            bindefs    = np.arange(xmin, xmax+binwidth_x, binwidth_x)

            hist = plt.hist(paramval,color="r",bins=bindefs,histtype="step",lw=1,label=r'Ntot='+str(Nparam)+labeladdition)

            plt.xlabel(xlabel)
            plt.ylabel(' Number of objects ')

            titletext = paramtype.replace('Line',line)+': mean='+str("%.4f" % np.mean(paramval))+'; std='+str("%.4f" % np.std(paramval))
            plt.title(titletext)

            leg = plt.legend(fancybox=True, loc='upper right',prop={'size':Fsize/1.0},ncol=1,numpoints=1)
            leg.get_frame().set_alpha(0.7)


            if verbose: print('   Saving plot to '+plotname)
            plt.savefig(plotname)
            plt.clf()
            plt.close('all')
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_redshiftdist(masterfits, addliterature=True, verbose=True, withancillaryhist=False, ylog=False,
                      outdir='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/paramdistributionFigures/'):
    """
    Function to plot the redshift distribution

    """
    masterdat = afits.open(masterfits)[1].data
    selection_det = np.where(( ((np.abs(masterdat['ferr_NV'])    != 99.0) & np.isfinite(masterdat['ferr_NV']))    |
                               ((np.abs(masterdat['ferr_CIV'])   != 99.0) & np.isfinite(masterdat['ferr_CIV']))   |
                               ((np.abs(masterdat['ferr_HeII'])  != 99.0) & np.isfinite(masterdat['ferr_HeII']))  |
                               ((np.abs(masterdat['ferr_OIII'])  != 99.0) & np.isfinite(masterdat['ferr_OIII']))  |
                               ((np.abs(masterdat['ferr_SiIII']) != 99.0) & np.isfinite(masterdat['ferr_SiIII'])) |
                               ((np.abs(masterdat['ferr_MgII']) != 99.0) & np.isfinite(masterdat['ferr_MgII']))   |
                               ((np.abs(masterdat['ferr_CIII'])  != 99.0) & np.isfinite(masterdat['ferr_CIII']))   ) &
                                (masterdat['redshift'] >= 0.0) & (masterdat['redshift'] <= 6.4432) & (masterdat['duplicationID'] == 0.0) &
                                (masterdat['id'] != 158002004) &
                                (masterdat['id'] != 601931670) &
                                (masterdat['id'] != 208014258) &
                                (masterdat['id'] != 600341002) )[0]

    selection_all = np.where( (masterdat['redshift'] >= 0.0) & (masterdat['redshift'] <= 6.4432) & (masterdat['duplicationID'] == 0.0) &
                              (masterdat['id'] != 158002004) &
                              (masterdat['id'] != 601931670) &
                              (masterdat['id'] != 208014258) &
                              (masterdat['id'] != 600341002) )[0]


    redshift_det = masterdat['redshift'][selection_det]
    redshift_all = masterdat['redshift'][selection_all]

    if addliterature:
        litcat  = '/Users/kschmidt/work/catalogs/literaturecollection_emissionlinestrengths/' \
                  'literaturecollection_emissionlinestrengths.fits'
        litdat  = afits.open(litcat)[1].data
        redshift_lit = litdat['redshift']

    if verbose: print(' - N_all = '+str(len(redshift_all)))
    if verbose: print(' - N_det = '+str(len(redshift_det)))
    if verbose: print(' - N_lit = '+str(len(redshift_lit)))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = outdir+'redshiftdistribution.pdf'
    fig = plt.figure(figsize=(5, 6))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.99, bottom=0.12, top=0.70)
    Fsize    = 17
    lthick   = 2.5
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    xmin = 0.0
    xmax = 8.0

    # Nbins   = np.ceil(np.sqrt(len(paramval)))
    binwidth_x = 0.28 # 0.32 #np.diff([xmin,xmax])/Nbins
    bindefs    = np.arange(xmin, xmax+binwidth_x, binwidth_x)+0.1
    bindefs[0] = 0.0
    bindefs[1] = 0.3
    bindefs[2] = 0.63

    hist_all = plt.hist(redshift_all,color="forestgreen",bins=bindefs,histtype="step",lw=lthick,zorder=5,label='$z>1.5$ MUSE-Wide and MUSE-Deep LSDCat \nobjects (this work)')
    hist_det = plt.hist(redshift_det,color="darkgreen",bins=bindefs,histtype="stepfilled",lw=lthick,zorder=8,label='$z>1.5$ MUSE objects with UV line detections \nfrom FELIS (this work)')
    hist_lit = plt.hist(redshift_lit,color="k",bins=bindefs,histtype="step",lw=lthick,zorder=10,linestyle=':',label=r'Literature comparison sample')

    if withancillaryhist:
        inami_dat = afits.open('/Users/kschmidt/work/catalogs/MUSE_GTO/MUSE_UDF_Inami17_combinedVizieR.fits')[1].data
        hist_ina = plt.hist(inami_dat['zMuse'],color="blue",bins=bindefs,histtype="step",alpha=1.0,
                            lw=0.5,label=r'MUSE-Deep Inami et al. (2017)',zorder=2)
        if verbose: print(' - N_ina = '+str(len(inami_dat['zMuse'])))

        mwdr1_dat = afits.open('/Users/kschmidt/work/catalogs/MUSE_GTO/MW_44fields_main_table_v1.0.fits')[1].data
        hist_urr = plt.hist(mwdr1_dat['Z'],color="darkorange",bins=bindefs,histtype="step",alpha=1.0,
                            lw=0.5,label=r'MUSE-Wide DR1 Urrutia et al. (2019)',zorder=2)
        if verbose: print(' - N_urr = '+str(len(mwdr1_dat['Z'])))

    plt.xlabel('Redshift')
    plt.ylabel('Number of objects')

    if ylog:
        plt.yscale('log')
        yrange = [1.0,1000]
        ydesert = 300.
        yoii    = 600.
        ylya    = 600.
        plotname = plotname.replace('.pdf','_ylog.pdf')
    else:
        yrange = [0,460]
        ydesert = 380.
        yoii    = 435.
        ylya    = 435.

    plt.ylim(yrange)
    plt.fill_between([1.5,2.9],[yrange[0],yrange[0]],[yrange[1],yrange[1]],zorder=1,color='lightgray')
    plt.text(2.2,ydesert,'MUSE\nredshift\ndesert',color='dimgray',fontsize=Fsize/1.3,zorder=1,ha='center')
    plt.text(1.35,yoii,'[OII]',color='black',fontsize=Fsize/1.3,zorder=1,ha='right')
    plt.text(3.05,ylya,'Ly$\\alpha$ $\\rightarrow$',color='black',fontsize=Fsize/1.3,zorder=1,ha='left')


    leg = plt.legend(fancybox=True, loc='lower center',prop={'size':Fsize/1.35},ncol=1,numpoints=1,
                     bbox_to_anchor=(0.5, 1.0))  # add the legend
    leg.get_frame().set_alpha(0.7)


    if verbose: print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_photometrySEDs(infofile, plotsample='desert',verbose=True,
                        outdir='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/paramdistributionFigures/'):
    """
    Function to plot the redshift distribution

    """
    if verbose: print(' - Loading photometric catalogs')
    catSkeltonGS  = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
    datSkeltonGS  = afits.open(catSkeltonGS)[1].data

    catSkeltonCOS = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat.FITS'
    datSkeltonCOS = afits.open(catSkeltonCOS)[1].data

    catRafelski   = '/Users/kschmidt/work/catalogs/rafelski/uvudf_rafelski_2015.fits'
    datRafelski   = afits.open(catRafelski)[1].data

    infofiledat   = afits.open(infofile)[1].data
    infofiledat   = infofiledat[np.where((infofiledat['id']<4.9e8) | (infofiledat['id']>5.9e8))[0]] # ignoring UDF MW mock ids

    if verbose: print(' - Getting objects satisfying sample selection "'+plotsample+'"')
    if plotsample == 'desert':
        objentries = np.where((infofiledat['redshift'] > 1.5) & (infofiledat['redshift'] < 2.9))[0]
    elif plotsample == 'laes':
        objentries = np.where((infofiledat['redshift'] > 2.9) & (infofiledat['redshift'] < 6.5))[0]
    elif plotsample == 'gt4':
        objentries = np.where((infofiledat['redshift'] > 4.0))[0]
    elif plotsample == 'lt2':
        objentries = np.where((infofiledat['redshift'] < 2.0))[0]
    else:
        sys.exit('Invalid choice of plotsampl="'+plotsample+'"')

    if verbose: print('   Found '+str(len(objentries))+' objects to plots SEDs for')
    if len(objentries) > 0:
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print(' - Plotting SEDs')
        plotname = outdir+'photometricSED_'+plotsample+'.pdf'
        fig = plt.figure(figsize=(5, 5))
        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.99, bottom=0.14, top=0.95)
        Fsize    = 17
        lthick   = 2.5
        marksize = 4
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()

        bandWeff = uves.band_waveeff('all')
        normband = 'F606W'
        normwave = bandWeff['F606W']

        for objent in objentries:
            photmatchlim=0.5
            f_out, photrefs, magABs = uves.get_SED(infofiledat,objent,photmatchlim=photmatchlim,
                                                   datSkeltonGS=datSkeltonGS, datSkeltonCOS=datSkeltonCOS, datRafelski=datRafelski)

            xvalues = []
            yvalues = []
            yerr    = []

            normflux = f_out[normband][0]

            for key in bandWeff.keys():
                xvalues.append(bandWeff[key])
                yvalues.append(f_out[key][0] / normflux)
                yerr.append(f_out[key][1])

            objcol = 'red'
            plt.errorbar(xvalues,yvalues,xerr=None,yerr=yerr,
                         marker='o',lw=0.5, markersize=3.0,alpha=0.5,
                         markerfacecolor=objcol,ecolor=objcol,color=objcol,
                         markeredgecolor='k',zorder=10)

        wavevec  = np.arange(4700,9400,10)
        betas    = [-1.2,-1.6,-1.97,-2.2,-2.4]
        for beta in betas:
            if beta == -1.97:
                linewidth = 2.0
                linesetup = 'k-'
            elif beta < -1.97:
                linewidth = 1.0
                linesetup = 'g--'
            elif beta > -1.97:
                linewidth = 1.0
                linesetup = 'b:'

            yval     = wavevec**(beta)
            norment  = np.where( np.abs(wavevec-normwave) == np.min(np.abs(wavevec-normwave)) )[0]
            plt.plot(wavevec,yval/yval[norment],linesetup,linewidth=linewidth,zorder=5,label='$\\beta$='+str(beta))

        yminsys, ymaxsys = plt.ylim()
        plt.plot([normwave,normwave],[yminsys, ymaxsys],'g:',linewidth=linewidth,zorder=5,label='$\\lambda$(Norm)')

        leg = plt.legend(fancybox=True,prop={'size':Fsize/1.35},ncol=1,numpoints=1) #, loc='upper right'

        leg.get_frame().set_alpha(0.7)

        plt.xlim([4000,11500])
        plt.ylim([0.4,5.0])

        plt.xscale('log')
        plt.yscale('log')

        plt.xlabel('$\\lambda$ [\AA]')
        plt.ylabel('flux [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]')

        if verbose: print(' - Saving plot to '+plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_nondetectionIDs(masterfits,outname='sources_extracted_using_pointsource_model.txt',verbose=True):
    """
    Function returning file with list of IDs of sources treated as object wihtout counterparts in the extractions

    """
    nondetdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_extraction_MWuves_100fields_maxdepth190808/tdose_nondetection_lists/'

    outfill_all     = nondetdir+outname
    outfill_wUVline = (nondetdir+outname).replace('.txt','_wUVlines.txt')

    nondetlists = glob.glob(nondetdir+'uves_nondetections*txt')
    if verbose: print(' - Getting the object ids of sources treated as nondetection in TDOSE extractions ')
    if verbose: print(' - Will store outputfiles ('+outname.replace('.txt','*')+') in\n   '+nondetdir)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Nnondet = 0
    fout = open(outfill_all,'w')
    fout.write('# The ids of objects treated as nondetections (no counterparts) in TDOSE extractions; i.e. with point source models\n')
    fout.write('# id\n')

    for nondetlist in nondetlists:
        if (not 'summary' in nondetlist.split('/')[-1]) & (not 'udf-0' in nondetlist.split('/')[-1]) & \
                (not 'skeltonbased' in nondetlist.split('/')[-1]):
            ids = np.atleast_1d(np.genfromtxt(nondetlist,dtype=None,comments='#'))
            if len(ids) > 0:
                for id in ids:
                    fout.write(str(id)+'\n')
                    Nnondet = Nnondet+1.0
    fout.write('# So that was '+str(Nnondet)+' objects IDs (with duplications)\n')
    fout.close()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    nondetids = np.genfromtxt(outfill_all,dtype=None,comments='#',skip_header=1)

    masterdat = afits.open(masterfits)[1].data

    selection_det = np.where(( ((np.abs(masterdat['ferr_NV'])    != 99.0) & np.isfinite(masterdat['ferr_NV']))    |
                               ((np.abs(masterdat['ferr_CIV'])   != 99.0) & np.isfinite(masterdat['ferr_CIV']))   |
                               ((np.abs(masterdat['ferr_HeII'])  != 99.0) & np.isfinite(masterdat['ferr_HeII']))  |
                               ((np.abs(masterdat['ferr_OIII'])  != 99.0) & np.isfinite(masterdat['ferr_OIII']))  |
                               ((np.abs(masterdat['ferr_SiIII']) != 99.0) & np.isfinite(masterdat['ferr_SiIII'])) |
                               ((np.abs(masterdat['ferr_MgII'])  != 99.0) & np.isfinite(masterdat['ferr_MgII']))   |
                               ((np.abs(masterdat['ferr_CIII'])  != 99.0) & np.isfinite(masterdat['ferr_CIII']))   ) &
                                (masterdat['redshift'] >= 0.0) & (masterdat['redshift'] <= 6.4432) & (masterdat['duplicationID'] == 0.0) &
                                (masterdat['id'] != 158002004) &
                                (masterdat['id'] != 601931670) &
                                (masterdat['id'] != 208014258) &
                                (masterdat['id'] != 600341002) )[0]

    UVlineids = masterdat['id'][selection_det]
    if verbose: print(' - Checking how many of these are among the '+str(len(UVlineids))+' objects with at least one UV line detection ')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Nuvnondet = 0
    fout = open(outfill_wUVline,'w')
    fout.write('# The ids of objects with at least one UV emission line treated as nondetections (no counterparts) in TDOSE extractions; i.e. with point source models\n')
    fout.write('# id\n')

    nondetids = np.genfromtxt(outfill_all,dtype=None,names=True,comments='#',skip_header=1)['id']
    for nondetid in nondetids:
        if nondetid in UVlineids:
            objent = np.where(masterdat['id'] == nondetid)[0]
            # fout.write(str(nondetid)+' # '+str(masterdat[objent])+'\n')
            fout.write(str(nondetid)+'\n')
            Nuvnondet = Nuvnondet+1.0

    fout.write('# So that was '+str(Nuvnondet)+' objects IDs (no duplications)\n')
    fout.close()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def count_EWlimits(masterfits='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/back2backAnalysis_200213/results_master_catalog_version200213.fits'):
    """

    """
    masterdat = afits.open(masterfits)[1].data

    for emline in ['NV','CIV','HeII','OIII','SiIII','CIII']:
        alldet_ent        = np.where(((np.abs(masterdat['ferr_'+emline]) != 99.0) & np.isfinite(masterdat['ferr_'+emline])) &
                                     (masterdat['redshift'] >= 0.0) & (masterdat['redshift'] <= 6.4432) & (masterdat['duplicationID'] == 0.0) &
                                     (masterdat['id'] != 158002004) &
                                     (masterdat['id'] != 601931670) &
                                     (masterdat['id'] != 208014258) &
                                     (masterdat['id'] != 600341002) )[0]

        contmag_limit_ent = np.where(((np.abs(masterdat['ferr_'+emline]) != 99.0) & np.isfinite(masterdat['ferr_'+emline])) &
                                     (masterdat['contmagABerr_'+emline] == -99) & # include lower limits (upper limits on brightness)
                                     (masterdat['redshift'] >= 0.0) & (masterdat['redshift'] <= 6.4432) & (masterdat['duplicationID'] == 0.0) &
                                     (masterdat['id'] != 158002004) &
                                     (masterdat['id'] != 601931670) &
                                     (masterdat['id'] != 208014258) &
                                     (masterdat['id'] != 600341002) )[0]

        EW0_limit_ent     = np.where(((np.abs(masterdat['ferr_'+emline]) != 99.0) & np.isfinite(masterdat['ferr_'+emline])) &
                                     ((np.abs(masterdat['EW0err_'+emline]) == 99) | (masterdat['EW0_'+emline] == 0.0)) &
                                     (masterdat['redshift'] >= 0.0) & (masterdat['redshift'] <= 6.4432) & (masterdat['duplicationID'] == 0.0) &
                                     (masterdat['id'] != 158002004) &
                                     (masterdat['id'] != 601931670) &
                                     (masterdat['id'] != 208014258) &
                                     (masterdat['id'] != 600341002) )[0]

        print(' - For '+emline+' there are the following:')
        print('   N(alldetections) = '+str(len(alldet_ent)))
        print('   N(contmaglimits) = '+str(len(contmag_limit_ent)))
        print('   N(EW0limits)     = '+str(len(EW0_limit_ent)))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def count_S2N(masterfits='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/back2backAnalysis_200213/results_master_catalog_version200213.fits',zdesert=False):
    """

    """

    if zdesert:
        zmin = 1.5
        zmax = 2.9
    else:
        zmin = 0.0
        zmax = 6.4432

    masterdat = afits.open(masterfits)[1].data

    SN5total = 0
    SN3total = 0
    for emline in ['NV','CIV','HeII','OIII','SiIII','CIII']:
        alldet_entSN5        = np.where(((np.abs(masterdat['ferr_'+emline]) != 99.0) & np.isfinite(masterdat['ferr_'+emline]) &
                                         (np.abs(masterdat['s2n_'+emline]) >= 5.0)) &
                                        (masterdat['redshift'] >= zmin) & (masterdat['redshift'] <= zmax) &
                                        (masterdat['duplicationID'] == 0.0) &
                                        (masterdat['id'] != 158002004) &
                                        (masterdat['id'] != 601931670) &
                                        (masterdat['id'] != 208014258) &
                                        (masterdat['id'] != 600341002) )[0]

        alldet_entSN3        = np.where(((np.abs(masterdat['ferr_'+emline]) != 99.0) & np.isfinite(masterdat['ferr_'+emline]) &
                                         (np.abs(masterdat['s2n_'+emline]) >= 3.0)) &
                                        (masterdat['redshift'] >= zmin) & (masterdat['redshift'] <= zmax) &
                                        (masterdat['duplicationID'] == 0.0) &
                                        (masterdat['id'] != 158002004) &
                                        (masterdat['id'] != 601931670) &
                                        (masterdat['id'] != 208014258) &
                                        (masterdat['id'] != 600341002) )[0]

        print(' - For '+emline+' there are the following:')
        print('   N(S/N>5)   = '+str(len(alldet_entSN5)))
        print('   N(S/N>3)   = '+str(len(alldet_entSN3)))
        print('   N(3<S/N<5) = '+str(len(alldet_entSN3)-len(alldet_entSN5)))


        SN5total = SN5total + int(len(alldet_entSN5))
        SN3total = SN3total + int(len(alldet_entSN3))


    selcountSN5 = np.where(( ((np.abs(masterdat['ferr_NV'])   != 99.0) & np.isfinite(masterdat['ferr_NV']))   & (np.abs(masterdat['s2n_NV'])   >=5.0) |
                             ((np.abs(masterdat['ferr_CIV'])  != 99.0) & np.isfinite(masterdat['ferr_CIV']))  & (np.abs(masterdat['s2n_CIV'])  >=5.0) |
                             ((np.abs(masterdat['ferr_HeII']) != 99.0) & np.isfinite(masterdat['ferr_HeII'])) & (np.abs(masterdat['s2n_HeII']) >=5.0) |
                             ((np.abs(masterdat['ferr_OIII']) != 99.0) & np.isfinite(masterdat['ferr_OIII'])) & (np.abs(masterdat['s2n_OIII']) >=5.0) |
                             ((np.abs(masterdat['ferr_SiIII'])!= 99.0) & np.isfinite(masterdat['ferr_SiIII']))& (np.abs(masterdat['s2n_SiIII']) >=5.0) |
                             ((np.abs(masterdat['ferr_MgII']) != 99.0) & np.isfinite(masterdat['ferr_MgII'])) & (np.abs(masterdat['s2n_MgII']) >=5.0) |
                             ((np.abs(masterdat['ferr_CIII']) != 99.0) & np.isfinite(masterdat['ferr_CIII'])) & (np.abs(masterdat['s2n_CIII']) >=5.0)  ) &
                           (masterdat['redshift'] >= zmin) & (masterdat['redshift'] <= zmax) & (masterdat['duplicationID'] == 0.0) &
                           (masterdat['id'] != 158002004) &
                           (masterdat['id'] != 601931670) &
                           (masterdat['id'] != 208014258) &
                           (masterdat['id'] != 600341002) )[0]

    selcountSN3 = np.where(( ((np.abs(masterdat['ferr_NV'])   != 99.0) & np.isfinite(masterdat['ferr_NV']))   & (np.abs(masterdat['s2n_NV'])   >=3.0) |
                             ((np.abs(masterdat['ferr_CIV'])  != 99.0) & np.isfinite(masterdat['ferr_CIV']))  & (np.abs(masterdat['s2n_CIV'])  >=3.0) |
                             ((np.abs(masterdat['ferr_HeII']) != 99.0) & np.isfinite(masterdat['ferr_HeII'])) & (np.abs(masterdat['s2n_HeII']) >=3.0) |
                             ((np.abs(masterdat['ferr_OIII']) != 99.0) & np.isfinite(masterdat['ferr_OIII'])) & (np.abs(masterdat['s2n_OIII']) >=3.0) |
                             ((np.abs(masterdat['ferr_SiIII'])!= 99.0) & np.isfinite(masterdat['ferr_SiIII']))& (np.abs(masterdat['s2n_SiIII']) >=3.0) |
                             ((np.abs(masterdat['ferr_MgII']) != 99.0) & np.isfinite(masterdat['ferr_MgII'])) & (np.abs(masterdat['s2n_MgII']) >=3.0) |
                             ((np.abs(masterdat['ferr_CIII']) != 99.0) & np.isfinite(masterdat['ferr_CIII'])) & (np.abs(masterdat['s2n_CIII']) >=3.0)  ) &
                           (masterdat['redshift'] >= zmin) & (masterdat['redshift'] <= zmax) & (masterdat['duplicationID'] == 0.0) &
                           (masterdat['id'] != 158002004) &
                           (masterdat['id'] != 601931670) &
                           (masterdat['id'] != 208014258) &
                           (masterdat['id'] != 600341002) )[0]

    print(' - In total there are the following number of total detections in the objects:')
    print('   N(S/N>5)   = '+str(SN5total))
    print('   N(S/N>3)   = '+str(SN3total))
    print('   N(3<S/N<5) = '+str(SN3total-SN5total))

    print(' - In total there are the following number of objects with at least one line:')
    print('   N(S/N>5)   = '+str(len(selcountSN5)))
    print('   N(S/N>3)   = '+str(len(selcountSN3)))
    print('   N(3<S/N<5) = '+str(len(selcountSN3)-len(selcountSN5)))


    print('-----> Counting these number of detections per objects')
    Nsingleline         = 0
    Nsingleline_withSN5 = 0
    Nmorelines          = 0
    Nmorelines_withSN5  = 0
    for ii, objid in enumerate(masterdat['id'][selcountSN3]):
        SNarr = np.array([masterdat['s2n_NV'][selcountSN3][ii], masterdat['s2n_CIV'][selcountSN3][ii],
                          masterdat['s2n_HeII'][selcountSN3][ii], masterdat['s2n_OIII'][selcountSN3][ii],
                          masterdat['s2n_SiIII'][selcountSN3][ii], masterdat['s2n_MgII'][selcountSN3][ii],
                          masterdat['s2n_CIII'][selcountSN3][ii]])

        # pdb.set_trace()
        NSN3lines = len(np.where(SNarr > 3.0)[0])
        if NSN3lines == 1:
            Nsingleline = Nsingleline + 1
            print(SNarr)
            if len(np.where(SNarr >= 5.0)[0]) == 1:
                Nsingleline_withSN5 = Nsingleline_withSN5 + 1
        elif NSN3lines > 1:
            Nmorelines  = Nmorelines + 1
            if len(np.where(SNarr >= 5.0)[0]) >= 1:
                Nmorelines_withSN5 = Nmorelines_withSN5 + 1

        NlowSN      = len(np.where((SNarr > 3.0) & (SNarr < 5.0))[0])
        NhighSN     = len(np.where(SNarr >= 5.0)[0])

    print(' - Of the objects there are the following:')
    print('   N(objects with  1 line with S/N>3)                   = '+str(Nsingleline))
    print('   of these objects with at least one S/N>=5 line are   = '+str(Nsingleline_withSN5))
    print('   N(objects with >1 line with S/N>3)                   = '+str(Nmorelines))
    print('   of these objects with at least one S/N>=5 line are   = '+str(Nmorelines_withSN5))


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
