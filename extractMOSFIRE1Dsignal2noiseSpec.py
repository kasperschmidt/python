#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# extractMOSFIRE1Dsignal2noiseSpec.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Extracting a 1D spectrum and S/N spectrom from the MOSFIRE eps and ivar files
# The 1D spectrum will be writting to a binary fits table.
#----------------------------
#   COMMENTS
#----------------------------
# This routine is based on extractMOSFIRE1Dspectrum.py and is essentially a cleaned
# version of that code with further add-ons, improvements and fixes.
#
# As you can see the noddamp keyword has disappeared compared to extractMOSFIRE1Dspectrum.py
# This happened after realizing that the *eps*.fits are shifted versions of the A-B diff images
# such that all the flux falls in the center row of the spectrum. See illutration in notes from
# 130505.
#----------------------------
#   INPUTS:
#----------------------------
# epsfits          : eps MOSFIRE reduction pipeline output file
# ivarfits         : ivar (inverse variance) MOSFIRE reduction pipeline output file
#                    If no ivar excists provide the string 'noivar' and a dummy ivar
#                    of 1s will be created.
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --matt           : set this keyword if the 2D spec is extracted with extract2DspecFromMaskMosaic.py
#                    from a mask mosaic reduced with Matt Augers reduction pipeline.
# --speccenter     : Center of spectrum. Default value is max of S/N map median-collapsed along dispersion direction.
#                    -- DEFAULT = max(dispersion-collapsed 2D spec) --
# --aperwidth      : Set constant width of extraction aperture. Expects an odd integer so width of extraction aperture
#                    becomes: speccenter +/- (aperwidth-1)/2. FWHM is ~2.92 pixels for YJ spec and ~2.87 for HK spec 
#                    according to McLean et al 2012. The percentage of light captured (assuming a normal distribution)
#                    is scipy.special.erf(FWHM*2.35/2/np.sqrt(2)). Hence, using aperwidth = 7 pix = 2.4 FWHM => 99.5% light
#                    -- DEFAULT = 7 pixels --
# --redshift       : This keyword expects a float giving the assumed redshift of the object. This will overplot a list
#                    of lines on the spectra for reference.
# --ivarscale      : For reduction post 130813 the inverse variance maps are too small. The correction factor appears to
#                    be sqrt(texp) where texp is the exposure time per frame. Use this keyword to correct for it.
# --verbose        : set -verbose to get info/messages printed to the screen
# --noplot         : set --noplot to skip timeconsuming plotting
# --eps            : saving created plots as eps files
# --png            : saving created plots as png files
# --show           : showning plots on screen
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# fitstable        : File (named as epsfits with extension _1DSNspec) containing a binary fits table
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x extractMOSFIRE1Dsignal2noiseSpec.py       (only required once)
#
# --- April 2013 run ---
# bash> extractMOSFIRE1Dsignal2noiseSpec.py borg_1510+1115_J_borg_1510+1115_0745_eps.fits borg_1510+1115_J_borg_1510+1115_0745_ivar.fits --aperwidth 7 --speccenter 26 --verbose
#
# bash> extractMOSFIRE1Dsignal2noiseSpec.py borg_1510+1115_Y_borg_1510+1115_0668_eps.fits borg_1510+1115_Y_borg_1510+1115_0668_ivar.fits --aperwidth 15 --speccenter 19 --verbose --redshift 0.644

# bash> extractMOSFIRE1Dsignal2noiseSpec.py borg_1510+1115_Y_borg_1510+1115_0745_eps.fits borg_1510+1115_Y_borg_1510+1115_0745_ivar.fits --aperwidth 15 --speccenter 26 --verbose
# BAD centerpix = Nrows/2.0 + offset/0.1799 + noddamp/0.1799 = 26.11 for Nrows = 34, offset = 0.39 arcsec, and noddamp 1.25 arcsec
# GOOD centerpix = Nrowseps/2.0 + offset/0.1799

# bash> extractMOSFIRE1Dsignal2noiseSpec.py borg_1510+1115_Y_borg_1510+1115_1524_eps.fits borg_1510+1115_Y_borg_1510+1115_1524_ivar.fits --aperwidth 15 --speccenter 17 --verbose
# BAD centerpix = Nrows/2.0 + offset/0.1799 + noddamp/0.1799 = 25.84 for Nrows = 33, offset = 0.43 arcsec, and noddamp 1.25 arcsec
# GOOD centerpix = Nrowseps/2.0 + offset/0.1799

# bash> extractMOSFIRE1Dsignal2noiseSpec.py rectified_HIP68767_HIP68767_Y_A-B_KBSmod_eps.fits noivar --aperwidth 15 --speccenter 1032 --verbose

# --- Extraction from raw data ---
# borg_1510+1115_J e.g. 'm130426_0308.fits', 'm130426_0310.fits', 'm130426_0312.fits'
# bash> extractMOSFIRE1Dsignal2noiseSpec.py m130426_0308.fits noivar --aperwidth 7 --speccenter 1636 --verbose

# borg_1510+1115_Y e.g.  'm130425_0217.fits', 'm130425_0219.fits', 'm130425_0221.fits'
# bash> extractMOSFIRE1Dsignal2noiseSpec.py m130425_0217.fits noivar --aperwidth 7 --speccenter 1636 --verbose

# borg_0951+3304_Y e.g. 'm130426_0176.fits', 'm130426_0178.fits', 'm130426_0180.fits'
# bash> extractMOSFIRE1Dsignal2noiseSpec.py m130426_0176.fits noivar --aperwidth 7 --speccenter 1636 --verbose

# borg_1437+5043_Y e.g.  'm130425_0174.fits', 'm130425_0176.fits', 'm130425_0178.fits'
# bash> extractMOSFIRE1Dsignal2noiseSpec.py m130425_0174.fits noivar --aperwidth 7 --speccenter 1636 --verbose

# --- Extraction spec from Matt-reduction ---
# bash> extractMOSFIRE1Dsignal2noiseSpec.py mosfire_ABsub_testobj1_eps.fits mosfire_ABsub_testobj1_ivar.fits --aperwidth 20 --speccenter 97 --verbose --matt 

# --- Marusa's "lines" 130625 ---
# bash> extractMOSFIRE1Dsignal2noiseSpec.py miki11M_Y_1_eps.fits miki11M_Y_1_ivar.fits --aperwidth 7 --speccenter 57 --verbose # 57,63

# bash> extractMOSFIRE1Dsignal2noiseSpec.py miki11M_J_7_eps.fits miki11M_J_7_ivar.fits --aperwidth 7 --speccenter 37 --verbose # 28,37

# bash> extractMOSFIRE1Dsignal2noiseSpec.py miki17M_Y_1_eps.fits miki17M_Y_1_ivar.fits --aperwidth 9 --speccenter 23 --verbose # 23,18
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-04-24  started (based on extractMOSFIRE1Dspectrum.py) by K. B. Schmidt (UCSB)
# 2013-05-15  added matt keyword. K. B. Schmidt (UCSB)
# 2013-06-05  added 'noivar' option for ivarfits. K. B. Schmidt (UCSB)
# 2013-08-15  added ivarscale keyword. K. B. Schmidt (UCSB)
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
parser.add_argument("--matt", action="store_true", help="Set if input is from Matt Augers reduction pipeline")
parser.add_argument("--speccenter", type=float, help="Center (row) of spectrum. If not provided S/N peak used")
parser.add_argument("--aperwidth", type=int, help="Width of fixed extraction aperture to use (instead of S/N cut)")
parser.add_argument("--redshift", type=float, help="Redshift of lines list to overplot spectra")
parser.add_argument("--ivarscale", action="store_true", help="Scale inverse variance maps with sqrt(texp) (reduction 130813)")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--noplot", action="store_true", help="Turn plotting off")
parser.add_argument("--eps", action="store_true", help="Turn plots into eps files")
parser.add_argument("--png", action="store_true", help="Turn plots into png files")
parser.add_argument("--show", action="store_true", help="Showing plots on screen")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# reading fits files
eps      = pyfits.getdata(args.epsfits,ext=0)   # reading data into array   [rows,columns]
sizeeps  = eps.shape
hdueps   = pyfits.open(args.epsfits)            # Load the FITS hdulist using pyfits
hdreps   = hdueps[0].header                     # get header

if args.ivarfits == 'noivar':
    ivar     = np.ones(eps.shape)
    sizeivar = ivar.shape
    hdrivar  = hdreps
else:
    ivar     = pyfits.getdata(args.ivarfits,ext=0)  # reading data into array   [rows,columns]
    sizeivar = ivar.shape
    hduivar  = pyfits.open(args.ivarfits)           # Load the FITS hdulist using pyfits
    hdrivar  = hduivar[0].header                    # get header
    if args.ivarscale: # for 130813 reductions the ivar maps are not scaled properly - correct for that.
         texp = hdrivar['TRUITIME']
         ivar = ivar*np.sqrt(texp)
if sizeeps != sizeivar: sys.exit('  ERROR: Size of data arrays loaded from eps '+str(sizeeps)+' and ivar '+str(sizeivar)+' fits are not the same --> ABORTING')

ivar[np.where(ivar <= 0)] = np.nan               # making sure there are no negative values in ivar array
#ivarmask = ivar*0.0                              # mask of all zeros
#ivarmask[np.isfinite(ivar)] = 1.0                # setting values in mask to 1 if entries in ivar are neither inf nor nan
SNR      = np.multiply(eps,np.sqrt(ivar))        # array with S/N
medSNR   = np.median(SNR,1)                      # S/N profile in spatial direction
#-------------------------------------------------------------------------------------------------------------
# Setting extraction aperture width
awidth = 7 # default value
if args.aperwidth: awidth = args.aperwidth
if args.verbose: print 'Width of aperture set to ',awidth,' pixels'
#-------------------------------------------------------------------------------------------------------------
# creating wavelength vector
wave     = range(sizeeps[1])                       # vector to contain wavelengths
wcs      = pywcs.WCS(hdreps)                       # Extract wcs (coordinate) information
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
Nlam     = len(wave)
if args.verbose: print 'Created wavelength vector from input files'
#-------------------------------------------------------------------------------------------------------------
# Load emission line list if redshift is given. 
if args.redshift: 
    ELs = np.genfromtxt('/Users/kasperborelloschmidt/work/emissionlinelist.txt', dtype=None, comments='#')
    Nobj = len(ELs)
    ELwave  = ELs['f0']
    ELname  = ELs['f1']
#-------------------------------------------------------------------------------------------------------------
# Setting center of spectrum
speccenter = float(np.where(medSNR == max(medSNR))[0])
if args.speccenter: speccenter = args.speccenter
if args.verbose: print 'Center pixel of spec (spatial direction) set to pixel (row-1) ',speccenter
#-------------------------------------------------------------------------------------------------------------
# getting rows of extraction aperture
goodrows     = np.arange(int(awidth))-(int(awidth)-1)/2+int(speccenter)        # rows of central positive spec
if args.verbose: print 'Defined extraction apertures'
#-------------------------------------------------------------------------------------------------------------
# creating aperture eps and ivar arrays for MOSFIRE DRP
apereps      = eps[goodrows,:]                                  # extraction aperture applied to 2D spec
aperivar     = ivar[goodrows,:]                                 # extraction aperture applied to inverse variance map

if args.matt: # appending negative spectrum if fits files are from Matt Auger's reduction
    specsep      = int(np.round(2*hdreps['NODAMP']/hdreps['ASPPIX']))               # seperation between spec
    goodrows_neg = np.arange(int(awidth))-(int(awidth)-1)/2+int(speccenter)+specsep # the rows of the negative spectrum

    apereps      = np.append(apereps,np.multiply(eps[goodrows_neg,:],-1.),0)        # adding rows for negative spec
    aperivar     = np.append(aperivar,ivar[goodrows_neg,:],0)                       # adding rows for negative spec
    aperivar[np.where(ivar <= 0)] = np.nan                                          # setting negative values in ivar array to NaNs
#-------------------------------------------------------------------------------------------------------------
# Creating the actual 1D spectra, error vectors, and noise spectrum 
specsum            = []
specmean           = []
specnormivar       = []
noise1D            = []
specerr_stdsignal  = []
specerr_meansigma2 = []
for ii in range(Nlam): # looping over wavelengths to make sure no ivar NaN-pixels are include in spectra
    entmask = np.isfinite(aperivar[:,ii])
    if len(np.where(entmask == True)[0]) == 0:
        specsum.append(0.0)
        specmean.append(0.0)
        specnormivar.append(0.0)
        noise1D.append(999.0)
        specerr_stdsignal.append(-99.0)
        specerr_meansigma2.append(-99.0)        
    else:
        specsum.append(np.sum(apereps[entmask,ii]))     # 1D spectrum: SIMPLE SUM
        specmean.append(np.mean(apereps[entmask,ii]))   # 1D spectrum: SIMPLE MEAN
        epsivarprod  = np.sum(np.multiply(apereps[entmask,ii],aperivar[entmask,ii]),0) # sum of signal normalized by variance
        sumivar      = np.sum(aperivar[entmask,ii])                                    # summing variances in extraction aperture
        if sumivar == 0: pdb.set_trace()
        specnormivar.append(epsivarprod/sumivar)        # 1D spectrum: INVERSE VARIANCE WEIGHTED, i.e., sum(dat*ivar)/sum(ivar)
        specerr_stdsignal.append(np.std(apereps[entmask,ii])) # standard deviation of values in each column
        apersigma2 = np.divide(aperivar[entmask,ii]*0.0+1.0,aperivar[entmask,ii])  # variance map in extraction aperture
        specerr_meansigma2.append(np.mean(apersigma2))                             # average variance in each column

        ones        = np.multiply(aperivar[entmask,ii],0.0)+1.0
        variancearr = np.divide(ones,aperivar[entmask,ii])
        noise1D.append(np.sqrt(np.sum(variancearr,0)))
#-------------------------------------------------------------------------------------------------------------
# checking spectrum for NaN pixels
pixcheck = np.isfinite(specsum)
goodpix  = np.where(pixcheck == True)[0]
Nbadpix  = len(np.where(pixcheck == False)[0])
if Nbadpix >= 1:
    print 'NB! Found '+str(Nbadpix)+' pixel with False result of np.isfinite. These are not written to fits file'
#-------------------------------------------------------------------------------------------------------------
# writing output to fits table
outname  = args.epsfits.replace('.fits','_1Dspec.fits')
if args.verbose: print 'writing binary table with 1D spectra to '+outname

col1     = pyfits.Column(name='WAVELENGTH', format='F', array= np.asarray(wave)[goodpix] )
col2     = pyfits.Column(name='SPEC1D_SUM', format='F', array=np.asarray(specsum)[goodpix] )
col3     = pyfits.Column(name='SPEC1D_MEAN', format='F', array=np.asarray(specmean)[goodpix] )
col4     = pyfits.Column(name='SPEC1D_NORMIVAR', format='F', array=np.asarray(specnormivar)[goodpix] )
col5     = pyfits.Column(name='SPEC1DERR_STDSIGNAL', format='F', array=np.asarray(specerr_stdsignal)[goodpix] )
col6     = pyfits.Column(name='SPEC1DERR_MEANSIG2', format='F', array=np.asarray(specerr_meansigma2)[goodpix] )
col7     = pyfits.Column(name='NOISE', format='F', array=np.asarray(noise1D)[goodpix] )
cols     = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7])

tbhdu    = pyfits.new_table(cols)          # creating table header

# writing hdrkeys:   '---KEY--',                  '----------------MAX LENGTH COMMENT-------------'
tbhdu.header.append(('SPECCEN ' ,speccenter       ,'The central row spectrum extracted around'),end=True)
tbhdu.header.append(('APWIDTH ' ,awidth           ,'Width of the extraction aperture'),end=True)
tbhdu.header.append(('TARGRA  ' ,hdreps['TARGRA'] ,'RA  of object'),end=True)
tbhdu.header.append(('TARGDEC ' ,hdreps['TARGDEC'],'DEC of object'),end=True)
tbhdu.header.append(('RA      ' ,hdreps['TARGRA'] ,'RA  of pointing'),end=True)
tbhdu.header.append(('DEC     ' ,hdreps['TARGDEC'],'DEC of pointing'),end=True)


hdu      = pyfits.PrimaryHDU()             # creating primary (minimal) header
thdulist = pyfits.HDUList([hdu, tbhdu])    # combine primary and table header to hdulist
thdulist.writeto(outname,clobber=True)     # write fits file (clobber=True overwrites excisting file)
#-------------------------------------------------------------------------------------------------------------
# write signal-to-noise ratios per pixel to image
SNRfits = args.epsfits.replace('.fits','_SNRmap.fits')
if args.verbose: print 'writing signal-to-noise map to '+SNRfits
hduSNR  = pyfits.PrimaryHDU(SNR)          # creating image header
hduSNR.header.add_comment('A signal-to-noise map [e/s * sqrt(ivar)] of '+args.epsfits) # adding comment to header
hdulist = pyfits.HDUList([hduSNR])        # turn header into to hdulist
hdulist.writeto(SNRfits,clobber=True)     # write fits file (clobber=True overwrites excisting file)
#-------------------------------------------------------------------------------------------------------------
if not args.noplot:
    if args.verbose: print ' '
    #-------------------------------------------------------------------------------------------------------------
    # Loading created fitsfile for plotting
    specdat   = pyfits.open(outname)
    specdatTB = specdat[1].data     # assuming first extension is a table and loading it
    goodent   = np.where(np.asarray(noise1D) != 999.0)
    if Nbadpix == 0:
        specdatTB = specdatTB[goodent]  # only using values where at least 1 pixel didn't have ivar = nan or inf (goodent)
    else:
        specdatTB = specdatTB

    if args.verbose: print 'Columns read from '+outname+' are:'
    if args.verbose: print specdatTB.columns  # printing the columns of the data table
    #-------------------------------------------------------------------------------------------------------------
    # PLOTTING
    import matplotlib.pyplot as plt 
    Fsize = 10
    plt.rc('text', usetex=True)                    # enabling LaTex rendering of text
    plt.rc('font', family='serif',size=Fsize)      # setting text font
    plt.rc('xtick', labelsize=Fsize) 
    plt.rc('ytick', labelsize=Fsize) 

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    plotname = args.epsfits.replace('.fits','_spec1D_sum')
    if args.verbose: print 'Creating figure '+plotname
    fig = plt.figure()  # create a figure object
    fig.clf()                                        # clearing figure
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure: subplot(numRows, numCols, plotNum)
    ax.grid(True,linestyle='-',color='0.75')

    Xmin = np.min(specdatTB.WAVELENGTH)
    Xmax = np.max(specdatTB.WAVELENGTH)
    Ymin = np.min(specdatTB.SPEC1D_SUM)
    Ymax = np.max(specdatTB.SPEC1D_SUM)
    #============================
    if args.redshift:
        ELwaveshift = ELwave*(1+args.redshift)
        goodwave = np.where((ELwaveshift > Xmin) & (ELwaveshift < Xmax))
        plt.vlines(ELwaveshift[goodwave],Ymin,Ymax, colors='0.5', linestyles='solid')
        for jj in range(len(goodwave[0])):
            ax.text(ELwaveshift[goodwave[0][jj]],Ymin,ELname[goodwave[0][jj]],horizontalalignment='center',verticalalignment='top',rotation=90)
    #============================

    labelstr = '$<\sigma^2>$, i.e., mean variance in extraction aperture'
    plt.errorbar(specdatTB.WAVELENGTH,specdatTB.SPEC1D_SUM,yerr=specdatTB.SPEC1DERR_MEANSIG2,label=labelstr, fmt='.')

    labelstr = '$\sum$S, i.e., summed signal in extraction aperture'
    plt.plot(specdatTB.WAVELENGTH,specdatTB.SPEC1D_SUM,label=labelstr)

    ax.set_xlabel('$\lambda$ [\AA]', fontsize=Fsize)
    ax.set_ylabel('electrons / second', fontsize=Fsize)
    ax.set_title(plotname.replace('_','\_'))

    plt.ylim(0,800)
    
    leg = plt.legend(fancybox=True, loc='upper left')  # add the legend in the middle of the plot
    leg.get_frame().set_alpha(0.7)  
  
    fig.savefig(plotname+'.pdf')
    if args.eps: fig.savefig(plotname+'.eps')
    if args.png: fig.savefig(plotname+'.png')
    if args.show: plt.show()  # draw plot on screen

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    plotname = args.epsfits.replace('.fits','_spec1D_normivar')
    if args.verbose: print 'Creating figure '+plotname
    fig = plt.figure()  # create a figure object
    fig.clf()                                        # clearing figure
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure: subplot(numRows, numCols, plotNum)
    ax.grid(True,linestyle='-',color='0.75')

    Xmin = np.min(specdatTB.WAVELENGTH)
    Xmax = np.max(specdatTB.WAVELENGTH)
    Ymin = np.min(specdatTB.SPEC1D_NORMIVAR)
    Ymax = np.max(specdatTB.SPEC1D_NORMIVAR)
    #============================
    if args.redshift:
        ELwaveshift = ELwave*(1+args.redshift)
        goodwave = np.where((ELwaveshift > Xmin) & (ELwaveshift < Xmax))
        plt.vlines(ELwaveshift[goodwave],Ymin,Ymax, colors='0.5', linestyles='solid')
        for jj in range(len(goodwave[0])):
            ax.text(ELwaveshift[goodwave[0][jj]],Ymin,ELname[goodwave[0][jj]],horizontalalignment='center',verticalalignment='top',rotation=90)
    #============================

    labelstr = '$<\sigma^2>$, i.e., mean variance in extraction aperture'
    plt.errorbar(specdatTB.WAVELENGTH,specdatTB.SPEC1D_NORMIVAR,yerr=specdatTB.SPEC1DERR_MEANSIG2,label=labelstr, fmt='.')

    labelstr = '$\sum$S * ivar / $\sum$ivar, i.e., inverse variance wigthed spectrum'
    plt.plot(specdatTB.WAVELENGTH,specdatTB.SPEC1D_NORMIVAR,label=labelstr)

    ax.set_xlabel('$\lambda$ [\AA]', fontsize=Fsize)
    ax.set_ylabel('electrons / second', fontsize=Fsize)
    ax.set_title(plotname.replace('_','\_'))

    leg = plt.legend(fancybox=True, loc='upper left')  # add the legend in the middle of the plot
    leg.get_frame().set_alpha(0.7)  
  
    fig.savefig(plotname+'.pdf')
    if args.eps: fig.savefig(plotname+'.eps')
    if args.png: fig.savefig(plotname+'.png')
    if args.show: plt.show()  # draw plot on screen

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    plotname = args.epsfits.replace('.fits','_spec1D_sum_SNR')
    if args.verbose: print 'Creating figure '+plotname
    fig = plt.figure()  # create a figure object
    fig.clf()                                        # clearing figure
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure: subplot(numRows, numCols, plotNum)
    ax.grid(True,linestyle='-',color='0.75')

    Xmin = np.min(specdatTB.WAVELENGTH)
    Xmax = np.max(specdatTB.WAVELENGTH)
    Ymin = np.min(specdatTB.SPEC1D_SUM/specdatTB.NOISE)
    Ymax = np.max(specdatTB.SPEC1D_SUM/specdatTB.NOISE)
    #============================
    if args.redshift:
        ELwaveshift = ELwave*(1+args.redshift)
        goodwave = np.where((ELwaveshift > Xmin) & (ELwaveshift < Xmax))
        plt.vlines(ELwaveshift[goodwave],Ymin,Ymax, colors='0.5', linestyles='solid')
        for jj in range(len(goodwave[0])):
            ax.text(ELwaveshift[goodwave[0][jj]],Ymin,ELname[goodwave[0][jj]],horizontalalignment='center',verticalalignment='top',rotation=90)
    #============================

    labelstr = '$<\sigma^2>$, i.e., mean variance in extraction aperture'
    plt.errorbar(specdatTB.WAVELENGTH,specdatTB.SPEC1D_SUM/specdatTB.NOISE,yerr=specdatTB.SPEC1DERR_MEANSIG2,label=labelstr, fmt='.')

    labelstr = '$\sum$S / N'
    plt.plot(specdatTB.WAVELENGTH,specdatTB.SPEC1D_SUM/specdatTB.NOISE,label=labelstr)

    ax.set_xlabel('$\lambda$ [\AA]', fontsize=Fsize)
    ax.set_ylabel('S / N', fontsize=Fsize)
    ax.set_title(plotname.replace('_','\_'))

    leg = plt.legend(fancybox=True, loc='upper left')  # add the legend in the middle of the plot
    leg.get_frame().set_alpha(0.7)  
  
    fig.savefig(plotname+'.pdf')
    if args.eps: fig.savefig(plotname+'.eps')
    if args.png: fig.savefig(plotname+'.png')
    if args.show: plt.show()  # draw plot on screen

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    plotname = args.epsfits.replace('.fits','_spec1D_normivar_SNR')
    if args.verbose: print 'Creating figure '+plotname
    fig = plt.figure()  # create a figure object
    fig.clf()                                        # clearing figure
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure: subplot(numRows, numCols, plotNum)
    ax.grid(True,linestyle='-',color='0.75')

    Xmin = np.min(specdatTB.WAVELENGTH)
    Xmax = np.max(specdatTB.WAVELENGTH)
    Ymin = np.min(specdatTB.SPEC1D_NORMIVAR/specdatTB.NOISE)
    Ymax = np.max(specdatTB.SPEC1D_NORMIVAR/specdatTB.NOISE)
    #============================
    if args.redshift:
        ELwaveshift = ELwave*(1+args.redshift)
        goodwave = np.where((ELwaveshift > Xmin) & (ELwaveshift < Xmax))
        plt.vlines(ELwaveshift[goodwave],Ymin,Ymax, colors='0.5', linestyles='solid')
        for jj in range(len(goodwave[0])):
            ax.text(ELwaveshift[goodwave[0][jj]],Ymin,ELname[goodwave[0][jj]],horizontalalignment='center',verticalalignment='top',rotation=90)
    #============================

    labelstr = '$<\sigma^2>$, i.e., mean variance in extraction aperture'
    plt.errorbar(specdatTB.WAVELENGTH,specdatTB.SPEC1D_NORMIVAR/specdatTB.NOISE,yerr=specdatTB.SPEC1DERR_MEANSIG2,label=labelstr, fmt='.')

    labelstr = '$\sum$S * ivar / $\sum$ivar / N'
    plt.plot(specdatTB.WAVELENGTH,specdatTB.SPEC1D_NORMIVAR/specdatTB.NOISE,label=labelstr)

    ax.set_xlabel('$\lambda$ [\AA]', fontsize=Fsize)
    ax.set_ylabel('S / N', fontsize=Fsize)
    ax.set_title(plotname.replace('_','\_'))

    leg = plt.legend(fancybox=True, loc='upper left')  # add the legend in the middle of the plot
    leg.get_frame().set_alpha(0.7)  
  
    fig.savefig(plotname+'.pdf')
    if args.eps: fig.savefig(plotname+'.eps')
    if args.png: fig.savefig(plotname+'.png')
    if args.show: plt.show()  # draw plot on screen


#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

