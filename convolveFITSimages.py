#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# convolveFITSimages.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Convolving a list of fits images with a gaussian (PSF) kernal using the 
# astropy.nddata.convolution.convolve.convolve function
#----------------------------
#   COMMENTS
#----------------------------
# Build around the BoRG drz images where the pixelscale is not known,
# hence the keyword to provide it to build the kernel.
#----------------------------
#   INPUTS:
#----------------------------
# fitslist         : Path AND name of file containing path and names of fits images to convolve.
#                    The convolved images will be written to the same directory as the original images.
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --kernelsize     : FWHM of kernel to convolve images with (in arcsec)
#                    Default size if keyword not given is SDSS-average FWHM of ~1.4 arcsec
# --pixsize        : Size of pixel in arcsec if different from the default 0.08 arcsec of drizzled 
#                    BoRG images. Used to find dimentions of kernel 
# --sexconfig      : setting this keyword will create a default Source Extractor config file for
#                    each image. The file is placed in the same directory as the input fits fils
#                    and the outputted convolved image. This makes it easy to run SExtractor on
#                    the convolved image by simply typing
#                    > sex imagename.fits -c imagename_SEXconfig.txt 
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# convolved fits   : Will write the convolved fits images to the same directory as input images in 
#                    fitslist with the extension _convolvedFWHMXXarcsec.fits
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x convolveFITSimages.py       (only required once)
# bash> convolveFITSimages.py '/Users/kasperborelloschmidt/work/BoRG/borg_fits/fits2convolve.txt' --verbose --kernelsize 0.5 --sexconfig

#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-10-30  started by K. B. Schmidt (UCSB)
# 2013-08-14  Fixing compatibility with lists of just 1 file with np.atleast_1d. K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse                 # argument managing
import sys                      # enabling arguments to code
import os                       # enabling command line runs and executing other Python scripts with os.system('string')
import subprocess               # get output from spawned command line processes
import numpy as np              # enable opening with genfromtxt
import getopt                   # used to extract/obtain the optional input
import astropy.nddata as astro  # enabling convolution
import pyfits                   # for opening fits file
import pdb
#----------------------------
#   FUNCTIONS
#----------------------------
def pathAname(str):                         # splitting string with path and name in to
    strsplit = str.split('/')               # splitting string
    name     = strsplit[-1]                 # saving filename (last entry of catsplit)
    slash    = '/'
    path=slash.join(strsplit[0:-1])         # putting path back together
    return [path,name]
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("fitslist", type=str, help="Path and name to list of FITS files to convolve")
# ---- optional arguments ----
parser.add_argument("--kernelsize", type=float, help="FWHM of gaussian kernel to convolve with (in arcsec).")
parser.add_argument("--pixsize", type=float, help="Pixel size of input fits in arcsec. Default = 0.08 arcsec/pix.")
parser.add_argument("--sexconfig", action="store_true", help="Creates default Source Extractor config file.")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# reading list of files
files  = np.genfromtxt(args.fitslist, dtype=None, comments='#')
files  = np.atleast_1d(files) # making sure files containing only 1 fitsimage also work
Nfiles = len(files)
if args.verbose: print ':: '+sys.argv[0]+' :: Found '+str(Nfiles)+' fitsfiles to convolve'
#-------------------------------------------------------------------------------------------------------------
# creating kernel for convolution
if args.kernelsize:
    sigma = args.kernelsize/2.35 # FWHM into sigma 
else:
    sigma = 1.4/2.35             # SDSS-average FWHM seeing
if args.verbose: print ':: '+sys.argv[0]+' :: Sigma of gaussian kernel set to '+str(sigma)+' arcsec '

if args.pixsize:
    pixdim = args.pixsize
else:
    pixdim = 0.08 # arcsec - default BoRG drizzle pixel size

kernelwidth = sigma/pixdim
if args.verbose: print ':: '+sys.argv[0]+' :: Corresponding to a kernel width of '+str(kernelwidth)+' pixels with '+str(pixdim)+' arcsec/pixel'

dim        = np.ceil(kernelwidth*5)
if dim%2 == 0: dim = dim+1  # if dim is divisable by 2 (even) add 1 to get odd dimension
kerneldim  = [dim,dim]
if args.verbose: print ':: '+sys.argv[0]+' :: Kernel dimensions will be '+str(kerneldim)
kerneltype = 'gaussian' # possibilities: 'gaussian', 'boxcar', 'tophat', 'brickwall', 'airy', 'trapezoid'
kernel     = astro.make_kernel(kerneldim, kernelwidth,kerneltype) # creating kernel
#-------------------------------------------------------------------------------------------------------------
# looping over images and convolving
for ii in range(Nfiles):
    fitsIN   = files[ii]
    PNfitsIN = pathAname(fitsIN)
    hdulist  = pyfits.open(fitsIN)  # reading HDU list of fits image
    arrIN    = hdulist[0].data   # loading image into array
    imgdim   = arrIN.shape

    if dim > min(imgdim): # checking that kernel dimension is smaller than fits image dimension
        print ':: '+sys.argv[0]+' :: ERROR Dimesion of the kernel is larger than the input image:'
        print fitsIN
        print ':: '+sys.argv[0]+' :: --> ABORTING'
        sys.exit()

    arrOUT = astro.convolve(arrIN, kernel)

    # creating and writing output fits file
    replacestr = '_convolvedFWHM'+str(sigma*2.35).replace('.','p')+'arcsec.fits'
    fitsOUT    = fitsIN.replace('.fits',replacestr)  # creating name of output fits image
    hduOUT     = pyfits.PrimaryHDU(arrOUT)
    hdulistOUT = pyfits.HDUList([hduOUT])

    # --- Editing header of file (keyword max 8 characters) ---
    hdulistOUT[0].header['inputimg'] = PNfitsIN[1]
    hdulistOUT[0].header['pixsize']  = (pixdim, '[arcsec/pix] size of pixels in FITSIN')
    hdulistOUT[0].header['kFWHM']    = (sigma*2.35, '[arcsec] FWHM of kernel (ksigma*2.35)')
    hdulistOUT[0].header['ksigma']   = (sigma, '[arcsec] width (sigma) of kernel')
    hdulistOUT[0].header['kwpix']    = (kernelwidth, '[pix] width (sigma) of kernel')
    hdulistOUT[0].header['kdim']     = (dim, '[pix] dimension of kernel: KDIMxKDIM')
    hdulistOUT[0].header['ktype']    = (kerneltype, 'Kernel used by astropy.nddata.make_kernel')
    hdulistOUT[0].header['history']  = ('Convolved version of INPUTIMG created with convolveFITSimages.py')    
    # Keywords copied from input image
    hdulistOUT[0].header['RA_TARG']  = (hdulist[0].header['RA_TARG'],'right ascension of the target (deg) (J2000)')
    hdulistOUT[0].header['DEC_TARG'] = (hdulist[0].header['DEC_TARG'],'declination of the target (deg) (J2000)')
    hdulistOUT[0].header['CRPIX1']   = (hdulist[0].header['CRPIX1'],'x-coordinate of reference pixel')                
    hdulistOUT[0].header['CRPIX2']   = (hdulist[0].header['CRPIX2'],'y-coordinate of reference pixel')                
    hdulistOUT[0].header['CRVAL1']   = (hdulist[0].header['CRVAL1'],'first axis value at reference pixel')            
    hdulistOUT[0].header['CRVAL2']   = (hdulist[0].header['CRVAL2'],'second axis value at reference pixel')           
    hdulistOUT[0].header['CTYPE1']   = (hdulist[0].header['CTYPE1'],'the coordinate type for the first axis')         
    hdulistOUT[0].header['CTYPE2']   = (hdulist[0].header['CTYPE2'],'the coordinate type for the second axis')        
    hdulistOUT[0].header['CD1_1']    = (hdulist[0].header['CD1_1'],'partial of first axis coordinate w.r.t. x')      
    hdulistOUT[0].header['CD1_2']    = (hdulist[0].header['CD1_2'],'partial of first axis coordinate w.r.t. y')
    hdulistOUT[0].header['CD2_1']    = (hdulist[0].header['CD2_1'],'partial of second axis coordinate w.r.t. x')
    hdulistOUT[0].header['CD2_2']    = (hdulist[0].header['CD2_2'],'partial of second axis coordinate w.r.t. y')
    # ---------------------------------------------------------
    try:
        hdulistOUT.writeto(fitsOUT)  # write array to output if file doesn't already excists
    except IOError:
        print ':: '+sys.argv[0]+' :: The output file '+fitsOUT+' already exists. Advancing to next file...' # if it excists print note
    else:
        if args.verbose: print ':: '+sys.argv[0]+' :: Created header and wrote convolved data to: ' # if written give status
        if args.verbose: print fitsOUT

        if args.sexconfig: # creating default Source Extractor config file for convolved image.
            sexconfig = fitsOUT.replace('.fits','_SEXconfig.txt')  # creating name of output fits image
            os.system('sex -d > '+sexconfig)                       # dumping default sextractor file
            PNsexconfig = pathAname(sexconfig)
            # --- modifying default config file ---
            import fileinput
            for line in fileinput.input(sexconfig,inplace=True):
                print line.rstrip().replace('test.cat',PNsexconfig[1].replace('config.txt','output.cat'))
            for line in fileinput.input(sexconfig,inplace=True):
                print line.rstrip().replace('default.param','/Users/kasperborelloschmidt/path_executables/sextractor-2.5.0/config/default.param')
            for line in fileinput.input(sexconfig,inplace=True):
                print line.rstrip().replace('default.conv','/Users/kasperborelloschmidt/path_executables/sextractor-2.5.0/config/default.conv')
            for line in fileinput.input(sexconfig,inplace=True):
                print line.rstrip().replace('default.nnw','/Users/kasperborelloschmidt/path_executables/sextractor-2.5.0/config/default.nnw')
            for line in fileinput.input(sexconfig,inplace=True):
                print line.rstrip().replace('CHECKIMAGE_TYPE  NONE','CHECKIMAGE_TYPE  OBJECTS')
            for line in fileinput.input(sexconfig,inplace=True):
                print line.rstrip().replace('SEEING_FWHM      1.2','SEEING_FWHM      '+str(sigma*2.35))  # 0.2 for HST (BoRG)
            for line in fileinput.input(sexconfig,inplace=True):
                print line.rstrip().replace('check.fits',PNsexconfig[1].replace('config.txt','check.fits'))

            if args.verbose: print ':: '+sys.argv[0]+' :: Created the Source Extractor config file: '
            if args.verbose: print sexconfig




#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

