#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# extractMOSFIRE1Dspectrum.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Extracting the 1D spectrum from the MOSFIRE eps files based on S/N cut
# The 1D spectrum will be writting to a binary fits table.
#----------------------------
#   COMMENTS
#----------------------------
#
# NOTE!!
#
# This code is out-dated. It has been substituted by extractMOSFIRE1Dsignal2noiseSpec.py
#----------------------------
#   INPUTS:
#----------------------------
# epsfits          : eps fits file outputted by the MOSFIRE reduction pipeline
# ivarfits         : ivar fits file outputted by the MOSFIRE reduction pipeline corresponding to epsfits
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --speccenter     : Center of spectrum. Default value is max of S/N map median-collapsed along dispersion direction.
# --snthresh       : S/N threshold for extracting 1D spec. Rows where S/N in median-collapsed profile drop below this
#                    value are included in the extracted spectrum. Default value is 1.0.
# --aperwidth      : To use an extraction aperture of a fixed widht used this keyword to set the width.
#                    Expect odd integer so width of extraction aperture becomes: speccenter +/- (aperwidth-1)/2.
#                    FWHM is ~2.92 pixels for YJ spec and ~2.87 for HK spec according to McLean et al 2012.
#                    The percentage of light captured (assuming normal distribution) is 
#                    scipy.special.erf(FWHM*2.35/2/np.sqrt(2)). Hence, using aperwidth = 7 pix = 2.4 FWHM => 99.5% light
# --noddamp        : The nodding amplitude in arcsec. This is used to extract the negative residual on each side of 
#                    the positive spectrum which also contain science information. This is done by applying the same 
#                    extraction aperture, multiplying the extraction by -1 and adding it to the 1D spectrum.
#                    The nodd amplitude is in the YOFFSET header keyword (pixel scale is in PSCALE)
#                    NOTE: Only effective if --aperwidth is also set, i.e., not with default S/N extractions
# --verbose        : set -verbose to get info/messages printed to the screen
# --plot           : save pdf diagnostic plots
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# fitstable        : File (named as epsfits with extension _1Dspec) containing a binary fits table with:
#                          WAVELENGTH        wavelength vector in Angstroms
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x extractMOSFIRE1Dspectrum.py       (only required once)
# bash> extractMOSFIRE1Dspectrum.py borg_0751+2917_Y_borg_0751+2917_0175_eps.fits borg_0751+2917_Y_borg_0751+2917_0175_ivar.fits --verbose --aperwidth 7 --plot --speccenter 63  --noddamp 1.25 ; --snthresh 0.4

# bash> extractMOSFIRE1Dspectrum.py borg_0751+2917_Y_borg_0751+2917_0675_eps.fits borg_0751+2917_Y_borg_0751+2917_0675_ivar.fits --verbose --aperwidth 7 --plot --speccenter 150 ; --snthresh 0.4

# bash> extractMOSFIRE1Dspectrum.py eps_Y_m130101_0290-0324_S14.fits ivar_Y_m130101_0290-0324_S14.fits --verbose --aperwidth 7 --plot  --speccenter 29
#
#
# --- April 2013 run ---
# bash> extractMOSFIRE1Dspectrum.py eps_Y_m130424_0075-0124_S16.fits ivar_Y_m130424_0075-0124_S16.fits --aperwidth 7 --noddamp 1.25 --speccenter 25 --plot --verbose

# bash> extractMOSFIRE1Dspectrum.py borg_0951+3304_Y_borg_0951+3304_0281_eps.fits  borg_0951+3304_Y_borg_0951+3304_0281_ivar.fits --aperwidth 7 --noddamp 1.25 --speccenter 60 --plot --verbose

# bash> extractMOSFIRE1Dspectrum.py borg_1437+5043_Y_borg_1437+5043_0070_r2_T12e_eps.fits borg_1437+5043_Y_borg_1437+5043_0070_r2_T12e_ivar.fits --aperwidth 7 --speccenter 29 --plot --verbose

# bash> extractMOSFIRE1Dspectrum.py borg_1510+1115_J_borg_1510+1115_0745_eps.fits borg_1510+1115_J_borg_1510+1115_0745_ivar.fits --aperwidth 7 --speccenter 26 --plot --verbose

#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-01-14  started by K. B. Schmidt (UCSB)
# 2013-01-23  added noddamp by K. B. Schmidt (UCSB)
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
parser.add_argument("epsfits",  type=str, help="Name of eps  fits file outouted by MOSFIRE pipeline")
parser.add_argument("ivarfits", type=str, help="Name of ivar fits file outouted by MOSFIRE pipeline")
# ---- optional arguments ----
parser.add_argument("--speccenter", type=float, help="Center (row) of spectrum. If not provided S/N peak used")
parser.add_argument("--snthresh", type=float, help="S/N threshold of median S/N/pix within extraction aperture")
parser.add_argument("--noddamp", type=float, help="Nodding amplitude in arcsec")
parser.add_argument("--aperwidth", type=int, help="Width of fixed extraction aperture to use (instead of S/N cut)")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--plot", action="store_true", help="Save diagnostic pdf plots")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
sys.exit('Program stopped as it has been replaced by extractMOSFIRE1Dsignal2noiseSpec.py')
#-------------------------------------------------------------------------------------------------------------
# reading fits files
eps      = pyfits.getdata(args.epsfits,ext=0)   # reading data into array   [rows,columns]
sizeeps  = eps.shape
hdueps   = pyfits.open(args.epsfits)            # Load the FITS hdulist using pyfits
hdreps   = hdueps[0].header                     # get header

ivar     = pyfits.getdata(args.ivarfits,ext=0)  # reading data into array   [rows,columns]
sizeivar = ivar.shape
hduivar  = pyfits.open(args.ivarfits)           # Load the FITS hdulist using pyfits
hdrivar  = hduivar[0].header                    # get header

if sizeeps != sizeivar: sys.exit('  ERROR: Size of data arrays loaded from eps '+str(sizeeps)+' and ivar '+str(sizeivar)+' fits are not the same --> ABORTING')
ivar[np.where(ivar <= 0)] = np.nan               # making sure there are no negative values in ivar array
SNR      = np.multiply(eps,np.sqrt(ivar))        # array with S/N
medSNR   = np.median(SNR,1)                      # S/N profile in spatial direction
#-------------------------------------------------------------------------------------------------------------
# if nodding position is given get pixel val
if args.noddamp:
    pixscale = hdreps['PSCALE']
    noddamppix = np.round(args.noddamp/pixscale) # the nodd amplitude in pixels
#-------------------------------------------------------------------------------------------------------------
# Creating plot
if args.plot:
     import matplotlib.pyplot as plt
     figname = args.epsfits.replace('.fits','_medSNR.pdf')
     fig = plt.figure()
     plt.plot(range(sizeeps[0]),medSNR)
     fig.savefig(figname)
     if args.verbose: print 'Created the figure '+figname
#-------------------------------------------------------------------------------------------------------------
# creating wavelength vector
wave     = range(sizeeps[1])                       # vector to contain wavelengths
wcs      = pywcs.WCS(hdreps)                    # Extract wcs (coordinate) information

Xmin     = 0.0
Xmax     = sizeeps[1]
Ymin     = 0.0
Ymax     = sizeeps[0]
pixcrd   = np.array([[Xmin,Ymin],[Xmax,Ymax]], np.float_)  # pixel coordinates to get 'sky'/cordinate values for
sky      = wcs.wcs_pix2sky(pixcrd, 1)  # convert to 'sky' values. 2nd arg. is "origin": 1 = a 1-based (Fortran-like) coordinates.
lammin   = sky[0,0]                             # wavelength of first column 
lammax   = sky[1,0]                             # wavelength of last  column 
dlam     = (lammax - lammin)/sizeeps[1]
wave     = np.multiply(lammin+np.multiply(wave,dlam),1e10)  # wavelengths in angstrom
#-------------------------------------------------------------------------------------------------------------
# define center of spectrum
if args.speccenter:                 # if speccenter is provided us it
    speccenter = args.speccenter
else:                               # or else estimate it from median SN profile
    maxmed     = max(medSNR)        # maximum value of median-S/N slit-profile
    speccenter = np.where(medSNR == maxmed)
    speccenter = float(speccenter[0])
#-------------------------------------------------------------------------------------------------------------
# getting rows of extraction aperture
SNthresh = 1                               # signal-to-noise threshold that S/N/pix is above when included in extraction aperture 
if args.snthresh: SNthresh = args.snthresh # overwriting default if given on commandline

if not args.aperwidth:  # if a fixed width of aperture is not given use S/N
    goodrows = np.where(medSNR > SNthresh)[0]
else:
    if not args.noddamp:
        goodrows = np.arange(args.aperwidth)-(args.aperwidth-1)/2.0+speccenter
        goodrows = goodrows.astype(int) # convert array to integers

        apereps   = eps[goodrows,:]                       # extraction aperture applied to 2D spec
        aperivar  = ivar[goodrows,:]                      # extraction aperture applied to inverse variance map
    else:
        goodrowsTOP = np.arange(args.aperwidth)-(args.aperwidth-1)/2.0+(speccenter+noddamppix) # rows of top negative spec
        goodrowsCEN = np.arange(args.aperwidth)-(args.aperwidth-1)/2.0+speccenter              # rows of central positive spec
        goodrowsBOT = np.arange(args.aperwidth)-(args.aperwidth-1)/2.0+(speccenter-noddamppix) # rows of bottom negative spec
        goodrows = np.concatenate((goodrowsBOT,goodrowsCEN,goodrowsTOP),axis=0) # combine goodrows
        goodrows = goodrows.astype(int) # convert array to integers

        posmask  = [1]*args.aperwidth
        negmask  = [-1]*args.aperwidth
        rowsmask = np.concatenate((negmask,posmask,negmask),axis=0) # combine goodrows
        rowsmask = np.tile(rowsmask,(Xmax,1)).T                     # dublicating rows and transposing them to have right dim.

        apereps   = np.multiply(eps[goodrows,:],rowsmask)                       # extraction aperture applied to 2D spec
        aperivar  = np.multiply(ivar[goodrows,:],rowsmask)                      # extraction aperture applied to inverse variance map
#-------------------------------------------------------------------------------------------------------------
# create 1D spectrum: SIMPLE SUM
specsum     = sum(apereps,1)                        # summing spectrum (no normalisation with inverse variance)
specerr     = np.std(apereps,0)                     # standard deviation of values in each column
#-------------------------------------------------------------------------------------------------------------
# create 1D spectrum: NORMALIZED BY INVERSE VARIANCE
sumivar     = sum(aperivar,1)                       # summing variances in extraction aperture
apersigma2  = np.divide(np.multiply(aperivar,0.0)+1,aperivar) # variance map in extraction aperture
sumsigma2   = sum(apersigma2,1)                     # summing variances in extraction aperture

specivarerr = np.mean(apersigma2,0)                 # average variance in each column

epsivarprod = sum(np.multiply(apereps,aperivar),1)  # the sum of the signal normalized by variance
specivar    = np.divide(epsivarprod,sumivar)        # summing spectrum normalising with inverse variance
#-------------------------------------------------------------------------------------------------------------
# create 1D spectra: MEAN AND NORM OF GOOD PIXELS (1/0 INVERSE VARIANCE MASK)
pixNaNs           = np.isnan(ivar)
ivarmask          = np.multiply(ivar,0.0)+1.0 # inverse variance mask indicating good (1) and bad (0) pixels
ivarmask[pixNaNs] = 0.0
ivarmaskaper      = ivarmask[goodrows,:]      # apply extraction aperture to inverse variance mask

specsummask       = specsum*0.0
specerrmask       = specsum*0.0
specivarmask      = specsum*0.0
specivarerrmask   = specsum*0.0
noise             = specsum*0.0
for ii in range(sizeeps[1]):                            # looping over columns and filling spectrum
    goodpix = np.where(ivarmaskaper[:,ii] == 1.0)
    if goodpix == []:                                   # if all pixels in column are bad set values to 0.0
        specsummask[ii]  = -99.0
        specivarmask[ii] = -99.0
        specsummask[ii]  = 0.0
        specivarmask[ii] = 0.0
    else:
        specsummask[ii]  = np.sum(apereps[goodpix,ii]) # the sum of good pixels stored as spectral value
        specerrmask[ii]  = np.std(apereps[goodpix,ii])  # standard deviation of pixels

        epsivarprodmask  = np.sum(np.multiply(apereps[goodpix,ii],aperivar[goodpix,ii]))  # sum of signal normalized by variance
        specivarmask[ii] = np.divide(epsivarprodmask,np.sum(aperivar[goodpix,ii]))                # normalising by inverse variance sum
        specivarerrmask[ii] = np.mean(apersigma2[goodpix,ii]) # mean variance in good pix
        noise[ii]        = np.sqrt(np.sum(apersigma2[goodpix,ii]))  # the total noise in the good pixels. sqrt(sum of variances)
                                                                    # can be used to estimate signal to noise: S/N = specsummask/noise
#-------------------------------------------------------------------------------------------------------------
# writing output to fits table
outname  = args.epsfits.replace('.fits','_1Dspec.fits')
if args.verbose: print 'writing binary table with 1D spectrum to '+outname

col1     = pyfits.Column(name='WAVELENGTH', format='F', array=wave)
col2     = pyfits.Column(name='SPEC_SUM', format='F', array=specsum)
col3     = pyfits.Column(name='SPECERR', format='F', array=specerr)
col4     = pyfits.Column(name='SPEC_NORM', format='F', array=specivar)
col5     = pyfits.Column(name='SPECERR_NORM', format='F', array=specivarerr)
col6     = pyfits.Column(name='SPEC_SUM_IVARMASK', format='F', array=specsummask)
col7     = pyfits.Column(name='SPECERR_IVARMASK', format='F', array=specerrmask)
col8     = pyfits.Column(name='SPEC_NORM_IVARMASK', format='F', array=specivarmask)
col9     = pyfits.Column(name='SPECERR_NORM_IVARMASK', format='F', array=specivarerrmask)
col10    = pyfits.Column(name='NOISE', format='F', array=noise)
cols     = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])

tbhdu    = pyfits.new_table(cols)          # creating table header

# writing hdrkeys:   '---KEY--',                  '----------------MAX LENGTH COMMENT-------------'
tbhdu.header.append(('SPECCEN' ,speccenter       ,'The central row spectrum extracted around'),end=True)
tbhdu.header.append(('NROWSUM' ,len(goodrows)    ,'The number of rows combined into spec'),end=True)
if not args.aperwidth: tbhdu.header.append(('SNTHRESH',SNthresh         ,'The S/N threshold of extracted spec rows'),end=True)
if args.aperwidth: tbhdu.header.append(('APWIDTH',args.aperwidth         ,'Width of the extraction aperture'),end=True)

hdu      = pyfits.PrimaryHDU()             # creating primary (minimal) header
thdulist = pyfits.HDUList([hdu, tbhdu])    # combine primary and table header to hdulist
thdulist.writeto(outname,clobber=True)     # write fits file (clobber=True overwrites excisting file)
#-------------------------------------------------------------------------------------------------------------
# write signal-to-noise ratios per pixel to image
SNRfits = args.epsfits.replace('.fits','_SNRmap.fits')
if args.verbose: print 'writing signal-to-noise map to '+SNRfits
hduSNR  = pyfits.PrimaryHDU(SNR)          # creating image header
hduSNR.header.add_comment('A signal-to-noise map of '+args.epsfits) # adding comment to header
hdulist = pyfits.HDUList([hduSNR])        # turn header into to hdulist
hdulist.writeto(SNRfits,clobber=True)     # write fits file (clobber=True overwrites excisting file)
#-------------------------------------------------------------------------------------------------------------
# write ivar mask to image
ivarmaskfits = args.epsfits.replace('.fits','_ivarmask.fits')
if args.verbose: print 'writing ivar mask to '+ivarmaskfits
hduIVM  = pyfits.PrimaryHDU(ivarmask)          # creating image header
hduIVM.header.add_comment('Inverse variance mask (1=good, 0=bad) for '+args.epsfits) # adding comment to header
hdulist = pyfits.HDUList([hduIVM])             # turn header into to hdulist
hdulist.writeto(ivarmaskfits,clobber=True)     # write fits file (clobber=True overwrites excisting file)
#-------------------------------------------------------------------------------------------------------------
# Creating plot
if args.plot:
     import matplotlib.pyplot as plt
     figname = args.epsfits.replace('.fits','_SPEC.pdf')
     fig = plt.figure()

     if args.aperwidth: labeltext = 'MeanSpecMask: CenPix = '+str(speccenter)+', AperWidth = '+str(args.aperwidth)
     if not args.aperwidth: labeltext = 'MeanSpecMask: CenPix = '+str(speccenter)+', Nrows = '+str(len(goodrows))+', S/N thresh = '+str(SNthresh)

     plt.errorbar(wave,specsummask, yerr=specivarerrmask, fmt='.',label='Mean variance in ivar mask: $<\sigma^2>$')
     plt.plot(wave,specsummask,label=labeltext)

     plt.plot(wave,specivarmask,label='Inverse variance weighted spec')

     #if args.aperwidth: labeltext = 'SumSpec: CenPix = '+str(speccenter)+', AperWidth = '+str(args.aperwidth)
     #if not args.aperwidth: labeltext ='SumSpec: CenPix = '+str(speccenter)+', Nrows = '+str(len(goodrows))+', S/N thresh = '+str(SNthresh)
     #plt.errorbar(wave,specsum, yerr=specerr, fmt='.',label='StdDev on values')
     #plt.plot(wave,specsum,label=labeltext)

     leg = plt.legend(loc='lower right',prop={'size':10})
     fig.savefig(figname)
     if args.verbose: print 'Created the figure '+figname
#-------------------------------------------------------------------------------------------------------------
# Creating plot
if args.plot:
    import pyfits
    import matplotlib.pyplot as plt
    figname = args.epsfits.replace('.fits','_SPEC_SN.pdf')
    data = pyfits.getdata(outname,ext=0)
    fig = plt.figure()

    labeltext = 'S/N spectrum (i.e. $\sum S_i /\sqrt{\sum \sigma^2_i}$ where $i$ is over rows in 2D mask)'

    plt.plot(wave, data.SPEC_SUM_IVARMASK/data.NOISE,label=labeltext)

    leg = plt.legend(fancybox=True, loc='lower right',prop={'size':10})
    leg.get_frame().set_alpha(0.7)

    fig.savefig(figname)
    if args.verbose: print 'Created the figure '+figname
    plt.show()
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

