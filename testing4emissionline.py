#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# testing4emissionline.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Testing if an emission line in a MOSFIRE spectrum is real or not
#----------------------------
#   COMMENTS
#----------------------------
# Uses the tests from emissionlinetests.py
#----------------------------
#   INPUTS:
#----------------------------
# epsfits          : eps MOSFIRE reduction pipeline output file
# ivarfits         : ivar (inverse variance) MOSFIRE reduction pipeline output file
# fits1D           : fits file with 1D extraction (created with extractMOSFIRE1Dsignal2noiseSpec.py)
# linecoord        : central coordinates of line(s) to test. Expects list on the form:
#                        [xpix1,ypix1,xpix2,ypix2,...,xpixN,ypixN]
#                    where N are the number of lines to test in the spectrum
#                    NOTE the actual center of the line can be obtained by using the offset given on the
#                    *xml.html page produced by MAGMA when creating the masks, knowing that there are
#                    0.1799 ''/pix, and how many pixels the slit span (in the *_pair*.fits file). 
#                    The center will then be given by
#                         centerpix = Nrows/2.0 + offset/0.1799
#                    To turn that into the position in the *eps.fits add Npix = Noddamp/0.1799''/pix
#                    NB! a more straight forward (correct) way of doing it is (see notes 130602):
#                        (Npix width of eps file)/2 + (offset from *xml.html)/(0.1799''/pix) 
# datainfo         : vetor giving size of slit, and estimated FWHM PSF of observations.
#                        [slitlength,slitwidth,psf_fwhm]
#                    all given in arc seconds
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --nocont         : if no clear continuum is present provide the expected pixel position of the object
#                    in the spatial direction. Used to estimate whether the line is offset wrt the object
#                    or not. Default is to use the max value of the median-collapsed slit profile
# --plot           : set --plot to display diagnostic plots (needs to close plots as they appear)
# --verbose        : set --verbose to get info/messages printed to the screen
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# Directoryfile    : Output file containing the name of the directories the output has been
#                    put into. Directories will be moved to directory of catlist
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x testing4emissionline.py       (only required once)
# bash> testing4emissionline.py borg_1510+1115_Y_borg_1510+1115_0668_eps.fits borg_1510+1115_Y_borg_1510+1115_0668_ivar.fits borg_1510+1115_Y_borg_1510+1115_0668_eps_1Dspec.fits [956,22,988,22,934,22] [7.0,0.7,0.5] --verbose --plot

# -- Marusa's "line" 130625 ---
# bash> testing4emissionline.py miki11M_Y_1_eps.fits miki11M_Y_1_ivar.fits miki11M_Y_1_eps_1Dspec.fits [832,57] [14.99,0.70,0.8] --verbose --plot --nocont 63
#offset = 98/2.+ 2.43/0.1799 ~ 63

# bash> testing4emissionline.py miki11M_J_7_eps.fits miki11M_J_7_ivar.fits miki11M_J_7_eps_1Dspec.fits [247,37] [7.01,0.70,0.8] --verbose --plot --nocont 28
#offset = 51/2 + 0.51/0.1799 ~ 28

# bash> testing4emissionline.py miki17M_Y_1_eps.fits miki17M_Y_1_ivar.fits miki17M_Y_1_eps_1Dspec.fits [963,18] [7.01,0.70,0.8] --verbose --plot --nocont 23
#offset = 52/2 - 0.52/0.1799 ~ 23


#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-05-03  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import emissionlinetests as elt  # The tests to perform
import argparse                   # argument managing
import sys                        # enabling arguments to code
import os                         # enabling command line runs and executing other Python scripts with os.system('string')
import numpy as np                # enable opening with genfromtxt
import pyfits
import pdb                        # for debugging with pdb.set_trace()
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("epsfits"  , type=str, help="Name of eps  fits file outouted by MOSFIRE pipeline")
parser.add_argument("ivarfits" , type=str, help="Name of ivar fits file outouted by MOSFIRE pipeline")
parser.add_argument("fits1D"   , type=str, help="Name of *_1Dspec.fits file created with extractMOSFIRE1Dsignal2noiseSpec.py")
parser.add_argument("linecoord", type=str, help="The center coordinates of line(s) to test")
parser.add_argument("datainfo" , type=str, help="Info about the data, i.e., [slitlength,slitwidth,psf_fwhm] ")
# ---- optional arguments ----
parser.add_argument("--nocont", type=int, help="Pixel position in spatial direction of object if no continuum present")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging ")
parser.add_argument("--plot", action="store_true", help="Set this flag to plot diagnostic plots")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# Loading data for spectrum
eps      = pyfits.getdata(args.epsfits,ext=0)   # reading data into array   [rows,columns]
sizeeps  = eps.shape
hdueps   = pyfits.open(args.epsfits)            # Load the FITS hdulist using pyfits
hdreps   = hdueps[0].header                     # get header

ivar     = pyfits.getdata(args.ivarfits,ext=0)  # reading data into array   [rows,columns]
sizeivar = ivar.shape
hduivar  = pyfits.open(args.ivarfits)           # Load the FITS hdulist using pyfits
hdrivar  = hduivar[0].header                    # get header

spec1D   = pyfits.open(args.fits1D)
spec1DTB = spec1D[1].data                       # assuming first extension is a table and loading it

Nwave    = len(spec1DTB.WAVELENGTH)
dlam     = (np.max(spec1DTB.WAVELENGTH)-np.min(spec1DTB.WAVELENGTH))/Nwave

linecen  = args.linecoord[1:-1].split(',')
Nlines   = len(linecen)/2.0

datainfo = args.datainfo[1:-1].split(',')
datainfo = [float(datainfo[0]),float(datainfo[1]),float(datainfo[2])]
pixscale = hdreps['PSCALE']
noddamp  = hdreps['YOFFSET']
#obj_pos  = hdreps['CRVAL2']  # Not working but it should be: https://groups.google.com/forum/#!topic/mosfire-drp/7CGGhwe12rI

if Nlines != round(Nlines): sys.exit("Didn't find 2 coordinates (xpix,ypix) for each of the lines to test --> ABORTING")
if args.verbose: print 'Was provided coordinates for ',Nlines,' emission lines to test in spectrum '
if args.verbose: print 'If all answers are no the line is probably real (maybe means not enabled yet)'
#-------------------------------------------------------------------------------------------------------------
# testing lines

plotval = 0
if args.plot: plotval = 1

if args.verbose: ' '
for ii in range(int(Nlines)): # looping over number of lines to test
    pix   = [float(linecen[0::2][ii]),float(linecen[1::2][ii])] # center pixel values
    cwave = np.min(spec1DTB.WAVELENGTH) + dlam*pix[0]           # central wavelength

    if args.verbose: print '-------------------------------------- LINE NO.',ii+1,'-------------------------------------- '
    if args.verbose: print '      [x,y] =',pix,' ; central wavelength = ',"%.2f" % cwave,' A'
    #---------------------------------------------------------------------------------------------------------
    # EL assymetric?
    skewlim = 1.0 # an (arbitratry) limit on weather the line is skewed
    coeffgauss,coeffskew  = elt.assymetry(spec1DTB.SPEC1D_NORMIVAR,pix[0],plot=plotval)

    if coeffskew[3] > skewlim : result = 'no'
    if coeffskew[3] <= skewlim : result = 'yes'
    if args.verbose: 
        print '    - Emission line symetric?                                                    ',result
        print '      The gaussian coefficients of the fit were [A,mu,sigma]:'
        print '          ',coeffgauss
        print '      A skewed gaussian fits has "skewness"    : ',"%.2f" % coeffskew[3]
    #---------------------------------------------------------------------------------------------------------
    # EL size > PSF?
    NsigSIE   = 4.0
    Xmin      = round(coeffgauss[1]-NsigSIE*coeffgauss[2])
    Xmax      = round(coeffgauss[1]+NsigSIE*coeffgauss[2])
    Ymin      = round(pix[1]-NsigSIE*coeffgauss[2])
    Ymax      = round(pix[1]+NsigSIE*coeffgauss[2])
    if Xmin < 0         : Xmin = 0
    if Xmax > sizeeps[1]: Xmax = sizeeps[1]
    if Ymin < 0         : Ymin = 0
    if Ymax > sizeeps[0]: Ymax = sizeeps[0]
    epscutout = eps[Ymin:Ymax,Xmin:Xmax] # cutting out Nsigma region around linecenter 
    LCcutout  = [round((Xmax-Xmin)/2.),round((Xmax-Xmin)/2.)] # the approximate linecenter in cutout
    gaus2Dfit_param, gaus2Dfit_cov  = elt.size2D(epscutout,LCcutout,plot=plotval)

    Xsize     = 2.3548*gaus2Dfit_param[3] # FWHM in dispersion direction 
    Ysize     = 2.3548*gaus2Dfit_param[4] # FWHM in spatial direction
    PSFsize   = datainfo[2]               # FWHM of PSF given by user

    resultx   = 'yes'
    resulty   = 'yes'
    if (Xsize > PSFsize): resultx = 'no'
    if (Ysize > PSFsize): resulty = 'no'

    if args.verbose: 
        print '    - Emssion line smaller than PSF in x?                                        ',resultx
        print '    - Emssion line smaller than PSF in y?                                        ',resulty
    #---------------------------------------------------------------------------------------------------------
    # EL on sky line residual?
    NsigSKY    =  1.0 # number of sigmas around center to check
    lamrange   =  [cwave - NsigSKY*coeffgauss[2]*dlam,cwave + NsigSKY*coeffgauss[2]*dlam]
    linematch  =  elt.skyline(lamrange)
    Slimit     =  10 # limit on the line strenght. Lines below this limit are ingnored. Estimated by looking at MOSFIRE reduction
    Nstrong    =  len(np.where(linematch[:,1] > Slimit)[0])
    if Nstrong >  0 : result = 'yes'
    if Nstrong == 0 : result = 'no'

    if args.verbose: 
        print '    - Emission line on OH sky line residual?                                     ',result
        print '      Looked in wavelength [A] range [',"%.2f" % lamrange[0],',',"%.2f" % lamrange[1],']'
        print '      for OH lines with relative strengths >',Slimit
    #---------------------------------------------------------------------------------------------------------
    # EL offset in slit?
    line_ypos          = pix[1]                                   # Given y pixel position of line
    yvalcen            = eps[:,round(coeffgauss[1])]              # values in y-direction at center of gauss fit
    maxy_pix           = np.where(yvalcen == np.max(yvalcen))[0]  # line position assuming max val at center
    slitprofile_med    = np.median(eps,1)                         # median-collapsed slit profile
    obj_pos            = np.where(slitprofile_med == np.max(slitprofile_med))[0] # object-position = max of slit profile
    if args.nocont: obj_pos[0] = args.nocont
    diffpos            = abs((line_ypos - obj_pos))               # offset in pixels
    diffpos_as         = diffpos[0]*pixscale                      # offset in arc seconds
    offsettol          = 0.50 # Tolerance on offset in arcsec
    if diffpos_as >  offsettol: result = 'yes'
    if diffpos_as <= offsettol: result = 'no' 
    if args.verbose: print '    - Emission line offset >',str(offsettol),"'' in slit wrt. object position?                ",result
    if args.verbose: print '      Offset between obj and EL =',"%.2f" % diffpos_as,"'' with pixscale =","%.2f" % pixscale,"''"
    #---------------------------------------------------------------------------------------------------------
    # EL = Noise peak?
    
    NMCMC  = -99.99
    result = 'maybe'
    if args.verbose: print '    - Emission line similar to noise peak?                                       ',result
    if args.verbose: print '      The number of MCMC noise realizations with similar strenght were',NMCMC
    #---------------------------------------------------------------------------------------------------------
    # Shadow image(s) detected?
    shift           = round(2*noddamp/pixscale)  # expected shift of shadow image in pixels
    peakratio       = 0.5 # The shadow peak should be ~0.5 of the A+B emission line
    rattol          = 0.1 # tolerance on difference in 2D gauss fit amplitude
    # --- SHADOW ABOVE ---
    shadowcen_above = [pix[0],pix[1]+shift]
    Xmin            = round(coeffgauss[1]-NsigSIE*coeffgauss[2])
    Xmax            = round(coeffgauss[1]+NsigSIE*coeffgauss[2])
    Ymin            = round(shadowcen_above[1]-NsigSIE*coeffgauss[2])
    Ymax            = round(shadowcen_above[1]+NsigSIE*coeffgauss[2])
    if Xmin < 0         : Xmin = 0
    if Xmax > sizeeps[1]: Xmax = sizeeps[1]
    if Ymin < 0         : Ymin = 0
    if Ymax > sizeeps[0]: Ymax = sizeeps[0]
    epscutout       = eps[Ymin:Ymax,Xmin:Xmax] # cutting out Nsigma region around linecenter 
    LCcutout        = [round((Xmax-Xmin)/2.),shadowcen_above[1]-Ymin] # the approximate linecenter in cutout
    SHabove_param, SHabove_cov  = elt.size2D(-1*epscutout,LCcutout,plot=plotval)

    diff_SHabove    = (SHabove_param[0]/gaus2Dfit_param[0]-peakratio)
    resultabove     = 'yes'
    if abs(diff_SHabove) < rattol: resultabove = 'no'

    # --- SHADOW ABOVE ---
    shadowcen_below = [pix[0],pix[1]-shift]
    Xmin            = round(coeffgauss[1]-NsigSIE*coeffgauss[2])
    Xmax            = round(coeffgauss[1]+NsigSIE*coeffgauss[2])
    Ymin            = round(shadowcen_below[1]-NsigSIE*coeffgauss[2])
    Ymax            = round(shadowcen_below[1]+NsigSIE*coeffgauss[2])
    if Xmin < 0         : Xmin = 0
    if Xmax > sizeeps[1]: Xmax = sizeeps[1]
    if Ymin < 0         : Ymin = 0
    if Ymax > sizeeps[0]: Ymax = sizeeps[0]
    epscutout       = eps[Ymin:Ymax,Xmin:Xmax] # cutting out Nsigma region around linecenter 
    LCcutout        = [round((Xmax-Xmin)/2.),shadowcen_below[1]-Ymin] # the approximate linecenter in cutout
    SHbelow_param, SHbelow_cov  = elt.size2D(-1*epscutout,LCcutout,plot=plotval)

    diff_SHbelow    = (SHbelow_param[0]/gaus2Dfit_param[0]-peakratio)
    resultbelow     = 'yes'
    if abs(diff_SHbelow) < rattol: resultbelow = 'no'

    result = 'maybe'
    if args.verbose: 
        print '    - Is a shadow missing above the line?                                        ',resultabove
        print '      A_2Dgaussfit_shadow/A_2Dgaussfit_line - ',"%.2f" % peakratio,' = ',"%.2f" % diff_SHabove
        print '    - Is a shadow missing below the line?                                        ',resultbelow
        print '      A_2Dgaussfit_shadow/A_2Dgaussfit_line - ',"%.2f" % peakratio,' = ',"%.2f" % diff_SHbelow
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

