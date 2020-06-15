# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Scripts, functions and routines to (enable) search for UV emission lines in MUSE data
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import os
import sys
import numpy as np
import MUSEWideUtilities as mu
import MXDF_exploration as mxdf
import glob
# import astropy
# import scipy
# import MiGs
import astropy.io.fits as afits
# import pyregion
# import pyfits as pyfitsOLD
# import datetime
# import shutil
# import time
# import fits2ascii as f2a
# import MUSEWidePlots as mwp
# import kbsutilities as kbs
# import tdose_utilities as tu
# from astropy import wcs
# from astropy.coordinates import SkyCoord
# import astropy.coordinates as acoord
# from astropy import units as u
# import subprocess
# import collections
# import uvEmissionlineSearch as uves
# from uncertainties import unumpy
# import matplotlib
# import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm
# from matplotlib.ticker import NullFormatter
# import NEOGALmodels as nm
# import felis_build_template as fbt
# import felis
# import literaturecollection_emissionlinestrengths as lce
# import stacking
# from itertools import combinations
# import pickle
# import re
# import pyneb as pn
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def masedaCIIIemitters_info():

    # names     = ['UDF10-22','UDF10-41','UDF10-42','UDF10-51','UDF10-64','UDF10-99','UDF10-164','UDF10-231',
    #              'UDF10-6664','UDF10-6668','UDF10-6670','UDF10-6674']
    # ras       = [53.15447,53.15288,53.15829,53.16518,53.16001,53.16308,53.16935,53.16210,53.16234,53.15342,53.16747,53.16656]
    # decs      = [-27.77144,-27.77250,-27.77745,-27.78161,-27.77100,-27.78560,-27.78498,-27.77256,-27.78444,-27.78104,-27.78183,-27.77526]
    # redshifts = [2.226,1.847,1.550,2.228,1.847,2.543,1.906,2.447,2.394,1.850,2.069,2.542]

    names     = ['UDF10-41','UDF10-42','UDF10-51','UDF10-99','UDF10-164','UDF10-231',
                 'UDF10-6664','UDF10-6668','UDF10-6670','UDF10-6674']
    ras       = [53.15288,53.15829,53.16518,53.16308,53.16935,53.16210,53.16234,53.15342,53.16747,53.16656]
    decs      = [-27.77250,-27.77745,-27.78161,-27.78560,-27.78498,-27.77256,-27.78444,-27.78104,-27.78183,-27.77526]
    redshifts = [1.847,1.550,2.228,2.543,1.906,2.447,2.394,1.850,2.069,2.542]

    return np.asarray(names), np.asarray(ras), np.asarray(decs), np.asarray(redshifts)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_narrowbandimages_masedaCIIIemitters(datacube, outputdir, clobber=False,verbose=True):
    """
    Wrapper to generate narrowband images of CIII emitters from Maseda+17

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import MXDF_exploration as mxdf
    datacube    = '/Users/kschmidt/work/MUSE/QtClassify/mxdf/mxdf_mfs-and-effvar-cube_corr_f32.fits'
    outputdir   = '/Users/kschmidt/work/MUSE/mxdf/masedaCIIIemitters/narrowbands/'

    mxdf.gen_narrowbandimages_masedaCIIIemitters(datacube,outputdir,clobber=False,verbose=True)

    """
    names, ras, decs, redshifts = mxdf.masedaCIIIemitters_info()

    wcenters = []
    dwaves   = []
    for nn, name in enumerate(names):
        wcenters.append(np.asarray([1906.68,1908.73]) * (redshifts[nn]+1))
        dwaves.append(np.asarray([1.0,1.0]) * (redshifts[nn]+1))

    mu.create_narrowband_subcube(datacube,ras,decs,5.0,5.0,wcenters,dwaves,outputdir,
                                 cube_ext=['DATA(MFS,C)','EFFVAR'],names=names,clobber=clobber)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_CIIIratioImages(overwrite=False,verbose=True):

    outdir = '/Users/kschmidt/work/MUSE/mxdf/masedaCIIIemitters/narrowbands/'

    ciii1img = outdir+'subcube_UDF10-41_5p0x5p0arcsec_narrowbandimage_cwave5428p31796dwave2p847.fits'
    ciii2img = outdir+'subcube_UDF10-41_5p0x5p0arcsec_narrowbandimage_cwave5434p15431dwave2p847.fits'
    mxdf.saveratioimg(ciii1img,ciii2img,
                      '/Users/kschmidt/work/MUSE/mxdf/masedaCIIIemitters/narrowbands/CIIIratioimag_UDF10-41.fits',overwrite=overwrite)

    ciii1img = outdir+'subcube_UDF10-42_5p0x5p0arcsec_narrowbandimage_cwave4862p034dwave2p55.fits'
    ciii2img = outdir+'subcube_UDF10-42_5p0x5p0arcsec_narrowbandimage_cwave4867p2615dwave2p55.fits'
    mxdf.saveratioimg(ciii1img,ciii2img,
                      '/Users/kschmidt/work/MUSE/mxdf/masedaCIIIemitters/narrowbands/CIIIratioimag_UDF10-42.fits',overwrite=overwrite)

    ciii1img = outdir+'subcube_UDF10-51_5p0x5p0arcsec_narrowbandimage_cwave6154p76304dwave3p228.fits'
    ciii2img = outdir+'subcube_UDF10-51_5p0x5p0arcsec_narrowbandimage_cwave6161p38044dwave3p228.fits'
    mxdf.saveratioimg(ciii1img,ciii2img,
                      '/Users/kschmidt/work/MUSE/mxdf/masedaCIIIemitters/narrowbands/CIIIratioimag_UDF10-51.fits',overwrite=overwrite)

    ciii1img = outdir+'subcube_UDF10-99_5p0x5p0arcsec_narrowbandimage_cwave6755p36724dwave3p543.fits'
    ciii2img = outdir+'subcube_UDF10-99_5p0x5p0arcsec_narrowbandimage_cwave6762p63039dwave3p543.fits'
    mxdf.saveratioimg(ciii1img,ciii2img,
                      '/Users/kschmidt/work/MUSE/mxdf/masedaCIIIemitters/narrowbands/CIIIratioimag_UDF10-99.fits',overwrite=overwrite)

    ciii1img = outdir+'subcube_UDF10-164_5p0x5p0arcsec_narrowbandimage_cwave5540p81208dwave2p906.fits'
    ciii2img = outdir+'subcube_UDF10-164_5p0x5p0arcsec_narrowbandimage_cwave5546p76938dwave2p906.fits'
    mxdf.saveratioimg(ciii1img,ciii2img,
                      '/Users/kschmidt/work/MUSE/mxdf/masedaCIIIemitters/narrowbands/CIIIratioimag_UDF10-164.fits',overwrite=overwrite)

    ciii1img = outdir+'subcube_UDF10-231_5p0x5p0arcsec_narrowbandimage_cwave6572p32596dwave3p447.fits'
    ciii2img = outdir+'subcube_UDF10-231_5p0x5p0arcsec_narrowbandimage_cwave6579p39231dwave3p447.fits'
    mxdf.saveratioimg(ciii1img,ciii2img,
                      '/Users/kschmidt/work/MUSE/mxdf/masedaCIIIemitters/narrowbands/CIIIratioimag_UDF10-231.fits',overwrite=overwrite)

    ciii1img = outdir+'subcube_UDF10-6664_5p0x5p0arcsec_narrowbandimage_cwave6471p27192dwave3p394.fits'
    ciii2img = outdir+'subcube_UDF10-6664_5p0x5p0arcsec_narrowbandimage_cwave6478p22962dwave3p394.fits'
    mxdf.saveratioimg(ciii1img,ciii2img,
                      '/Users/kschmidt/work/MUSE/mxdf/masedaCIIIemitters/narrowbands/CIIIratioimag_UDF10-6664.fits',overwrite=overwrite)

    ciii1img = outdir+'subcube_UDF10-6668_5p0x5p0arcsec_narrowbandimage_cwave5434p038dwave2p85.fits'
    ciii2img = outdir+'subcube_UDF10-6668_5p0x5p0arcsec_narrowbandimage_cwave5439p8805dwave2p85.fits'
    mxdf.saveratioimg(ciii1img,ciii2img,
                      '/Users/kschmidt/work/MUSE/mxdf/masedaCIIIemitters/narrowbands/CIIIratioimag_UDF10-6668.fits',overwrite=overwrite)

    ciii1img = outdir+'subcube_UDF10-6670_5p0x5p0arcsec_narrowbandimage_cwave5851p60092dwave3p069.fits'
    ciii2img = outdir+'subcube_UDF10-6670_5p0x5p0arcsec_narrowbandimage_cwave5857p89237dwave3p069.fits'
    mxdf.saveratioimg(ciii1img,ciii2img,
                      '/Users/kschmidt/work/MUSE/mxdf/masedaCIIIemitters/narrowbands/CIIIratioimag_UDF10-6670.fits',overwrite=overwrite)

    ciii1img = outdir+'subcube_UDF10-6674_5p0x5p0arcsec_narrowbandimage_cwave6753p46056dwave3p542.fits'
    ciii2img = outdir+'subcube_UDF10-6674_5p0x5p0arcsec_narrowbandimage_cwave6760p72166dwave3p542.fits'
    mxdf.saveratioimg(ciii1img,ciii2img,
                      '/Users/kschmidt/work/MUSE/mxdf/masedaCIIIemitters/narrowbands/CIIIratioimag_UDF10-6674.fits',overwrite=overwrite)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def saveratioimg(ciii1img,ciii2img,outname,overwrite=False):

    ciii1dat = afits.open(ciii1img)[0].data
    ciii2dat = afits.open(ciii2img)[0].data

    hduimg   = afits.PrimaryHDU(ciii1dat/ciii2dat)
    hdus     = [hduimg]
    hdulist  = afits.HDUList(hdus)                  # turn header into to hdulist
    hdulist.writeto(outname,overwrite=overwrite)  # write fits file (clobber=True overwrites excisting file)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =