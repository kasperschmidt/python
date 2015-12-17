#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# extract2DspecFromMaskMosaic.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Extracting individual 2D spectra from a mask mosaic (from M Augers's reduction)
# where the idividual spectra are seperated by regions of NaNs
#----------------------------
#   COMMENTS
#----------------------------
# 
#----------------------------
#   INPUTS:
#----------------------------
# maskfits         : file containing the mask mosaic to extract spectra from
# traces           : The pixel-row location of the positive traces and an approximate width
#                    of the trace in the individual slits. These are ignored when creating the
#                    inverse variance maps. The 'nodamp' is used to determine the location
#                    of the negative traces which are also ingnored. The format expected is:
#                        [rowtrace1,widthtrace1,rowtrace2,widthtrace2,...,rowtraceN,widthtraceN]
#                    with 1 being the bottom slit and N being the top slit. 
#                    NOTE: each spectra has a slight tilt so make sure to make the cutout width large
#                          enough to account for this. Teh difference is ~5 pixels
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --nodamp        : provide the nodding amplitude used in arcsec. Default is 1.25''
# --objnames       : list of object names on the form [name1,name2,...,nameNslit] to use naming the
#                    2D spectra (from bottom to top) instead of using the default numbering.
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# fits files with 2D spectra
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x extract2DspecFromMaskMosaic.py       (only required once)
#
# bash> extract2DspecFromMaskMosaic.py /Users/kasperborelloschmidt/work/python/progs/reductioncode_MattA/testfiles/mosfire_ABsub.fits [98,20,376,20,592,20,696,20,817,20,961,20,1061,20,1259,20,1498,20,1865,20,1951,20] --verbose --objnames [testobj1,testobj2,testobj3,testobj4,testobj5,testobj6,testobj7,testobj8,testobj9,testobj10,testobj11] --nodamp 2.5
#
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-04-24  started (based on extractMOSFIRE1Dspectrum.py) by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs #
import argparse            # argument managing
import sys                 # enabling arguments to code
import getopt              # used to extract/obtain the optional input
import numpy as np         # enable opening with genfromtxt
import pyfits
import pywcs
import math
import pdb                 # for debugging
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("maskfits",  type=str, help="Name of fits file containing mask mosaic of spectra")
parser.add_argument("traces",  type=str, help="Location and width of spectral traces in slits")
# ---- optional arguments ----
parser.add_argument("--nodamp", type=float, help="Nodding amplitude (DEFAULT = 1.25 arcsec)")
parser.add_argument("--objnames", type=str, help="List of names of objects in slits on the form: [name1,name2,...,nameNslit]")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# reading fits files
mask     = pyfits.getdata(args.maskfits,ext=0)   # reading data into array   [rows,columns]
sizemask = mask.shape
hdumask  = pyfits.open(args.maskfits)            # Load the FITS hdulist using pyfits
hdrmask  = hdumask[0].header                     # get header

nodamp   = 1.25                                  # arcsec. Default nodding amplitude
if args.nodamp: nodamp = args.nodamp             # overwriting default
#-------------------------------------------------------------------------------------------------------------
# creating wavelength vector
wave     = range(sizemask[1])                    # vector to contain wavelengths
wcs      = pywcs.WCS(hdrmask)                    # Extract wcs (coordinate) information
Xmin     = 0.0
Xmax     = sizemask[1]
Ymin     = 0.0
Ymax     = sizemask[0]
pixcrd   = np.array([[Xmin,Ymin],[Xmax,Ymax]], np.float_)  # pixel coordinates to get 'sky'/cordinate values for
sky      = wcs.wcs_pix2sky(pixcrd, 1)   # convert to 'sky' values. 2nd arg. is "origin": 1 = a 1-based (Fortran-like) coordinates.
lammin   = sky[0,0]                              # wavelength of first column 
lammax   = sky[1,0]                              # wavelength of last  column 
dlam     = (lammax - lammin)/sizemask[1]
wave     = np.multiply(lammin+np.multiply(wave,dlam),1e10)  # wavelengths in angstrom
Nlam     = len(wave)
if args.verbose: print 'Created wavelength vector from input files'
#-------------------------------------------------------------------------------------------------------------
cencol   = mask[:,round(sizemask[1]/2.)] # extracting center column for defining 2Dpecs
finiteval= np.isfinite(cencol)
nanent   = np.where(finiteval == False)  # rows where values are not finite, i.e. NaN-gaps between spectra
Nspec    = len(nanent[0])/5.+1.          # the number of spectra in mosaic
if Nspec != round(Nspec):
    sys.exit('\nERROR! Number of rows with NaN is not a multiplicative of 5 as expected --> ABORTING')
if args.verbose: print 'Identified',Nspec,'spectra in the mask mosaic'
#-------------------------------------------------------------------------------------------------------------
if args.objnames:
    names = args.objnames[1:-1].split(',')
    if Nspec != len(names):
        sys.exit('\nERROR! The number of names given does not match the number of 2D spectra found --> ABORTING')
    Nspec = int(Nspec) # making sure Nspec is an interger 
else:
    Nspec = int(Nspec) # making sure Nspec is an interger 
    names = range(Nspec)
#-------------------------------------------------------------------------------------------------------------
traceval  = np.asarray(args.traces[1:-1].split(','),dtype=int)
if Nspec*2 != len(traceval):
        sys.exit('\nERROR! The number of traces given does not match the number of 2D spectra found --> ABORTING')
asperpix  = 0.1799 # ''/pix
Npixsep   = round(2*nodamp/asperpix) # seperation between positive and negative spec
#-------------------------------------------------------------------------------------------------------------
# getting min and mac row numbers of the Nspec spectra
specrows = np.zeros((Nspec, 2),dtype=int)          # list to contain min and max row number (not index) of individual 2D spectra
nanrows  = np.add(nanent,1)                        # row numbers of NaN-ros (not index)
for ii in range(Nspec):
    if ii == 0:                        # First spec
        specrows[ii,0] = 1+1
        specrows[ii,1] = nanent[0][0]
    elif ii == Nspec-1:                # Last spec
        specrows[ii,0] = specrows[ii-1,1]+5+1
        specrows[ii,1] = sizemask[0]
    else:                              # Remaining spec
        specrows[ii,0] = specrows[ii-1,1]+5+1
        specrows[ii,1] = nanent[0][5*ii]
#-------------------------------------------------------------------------------------------------------------
# write individual spectra to fits images
for jj in range(Nspec):
    spec2Dname = args.maskfits.replace('.fits','_'+str(names[jj])+'_eps.fits')

    specstart  = specrows[jj,0]-1
    specend    = specrows[jj,1]
    spec2D     = mask[specstart:specend,:]

    hduS2D     = pyfits.PrimaryHDU(spec2D)       # creating image header
    hduS2D.header.add_comment('2D spec cropped from:'+args.maskfits) # adding comment to header

    # writing hdrkeys:    '---KEY--',                   '----------------MAX LENGTH COMMENT-------------'
    hduS2D.header.append(('OBJNAME ',names[jj]         ,'Name of object in spectrum'),end=True)
    hduS2D.header.append(('ASPPIX  ',asperpix          ,'arsec/pixel'),end=True)
    hduS2D.header.append(('NODAMP  ',nodamp            ,'[arcsec] Dist. between spec = 2xNODAMP/ASPPIX'),end=True)
    hduS2D.header.append(('CTYPE1  ','LINEAR  '        ,' '),end=True)                                     
    hduS2D.header.append(('CRPIX1  ',hdrmask['CRPIX1'] ,' '),end=True)
    hduS2D.header.append(('CRVAL1  ',hdrmask['CRVAL1'] ,' '),end=True)
    hduS2D.header.append(('CD1_1   ',hdrmask['CD1_1']  ,' '),end=True)
    hduS2D.header.append(('CRPIX2  ',hdrmask['CRPIX2'] ,' '),end=True)
    hduS2D.header.append(('CRVAL2  ',hdrmask['CRVAL2'] ,' '),end=True)
    hduS2D.header.append(('CD2_2   ',hdrmask['CD2_2']  ,' '),end=True)

    hdulist = pyfits.HDUList([hduS2D])        # turn header into to hdulist
    hdulist.writeto(spec2Dname,clobber=True)     # write fits file (clobber=True overwrites excisting file)
    if args.verbose: print 'wrote 2D spectrum to '+spec2Dname

    # create and write ivar map to fits image
    ivarname = args.maskfits.replace('.fits','_'+str(names[jj])+'_ivar.fits')

    stdvec             = []
    dummy              = spec2D                                                                 # dummy array to manipulate
    traceent           = range(traceval[jj*2+1])-traceval[jj*2+1]/2+traceval[jj*2]-1-specstart  # rows with positive trace
    traceent           = np.append(traceent,traceent+int(Npixsep))                              # appending rows with negative trace
    dummy[traceent,:]  = np.nan                                                                 # filling rows of trace with NaNs
    for kk in range(spec2D.shape[1]):
        goodrows           = np.where(np.isfinite(dummy[:,kk]) == True)
        stdvec.append(np.std(spec2D[goodrows,kk]))           # standard deviation in each column KBS?? really optimal... clipping?
    ivar     = np.tile(stdvec,(spec2D.shape[0],1))           # creating array with stdvec in all rows

    hduIVAR  = pyfits.PrimaryHDU(ivar)       # creating image header
    hduIVAR.header.add_comment('Inverse variance map of '+ivarname) # adding comment to header

    # writing hdrkeys:     '---KEY--',                  '----------------MAX LENGTH COMMENT-------------'
    hduIVAR.header.append(('OBJNAME',names[jj]         ,'Name of object in spectrum'),end=True)
    hduIVAR.header.append(('ASPPIX  ',asperpix          ,'arsec/pixel'),end=True)
    hduIVAR.header.append(('NODAMP  ',nodamp            ,'[arcsec] Dist. between spec = 2xNODAMP/ASPPIX'),end=True)
    hduIVAR.header.append(('CTYPE1 ','LINEAR  '        ,' '),end=True)                                     
    hduIVAR.header.append(('CRPIX1 ',hdrmask['CRPIX1'] ,' '),end=True)
    hduIVAR.header.append(('CRVAL1 ',hdrmask['CRVAL1'] ,' '),end=True)
    hduIVAR.header.append(('CD1_1  ',hdrmask['CD1_1']  ,' '),end=True)
    hduIVAR.header.append(('CRPIX2 ',hdrmask['CRPIX2'] ,' '),end=True)
    hduIVAR.header.append(('CRVAL2 ',hdrmask['CRVAL2'] ,' '),end=True)
    hduIVAR.header.append(('CD2_2  ',hdrmask['CD2_2']  ,' '),end=True)

    hdulist = pyfits.HDUList([hduIVAR])        # turn header into to hdulist
    hdulist.writeto(ivarname,clobber=True)     # write fits file (clobber=True overwrites excisting file)
    if args.verbose: print 'wrote ivar map to     '+ivarname

if args.stop: pdb.set_trace()
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

