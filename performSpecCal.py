#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# performSpecCal.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# script using calibrate1Dspectrum.py to calibrate 1D spectra
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# paramfile        : file containing the needed input. A template can be found 
#                    in performSpecCal_input.param but it can be created easily
#                    with performSpecCal_createParamFile.py
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --plot           : set --plot to display diagnostic plots (needs to close plots as they appear)
# --verbose        : set --verbose to get info/messages printed to the screen
# --stop           : stoppping program before end for de-bugging
# --show           : showning plots on screen for manipulation and saving
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# 
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x performSpecCal.py       (only required once)
# bash> performSpecCal.py /Users/kasperborelloschmidt/work/python/myprogs/performSpecCal_input.param --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-05-24  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs        # 
import calibrate1Dspectrum as cal # The calibration steps
import datetime
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
parser.add_argument("paramfile"  , type=str  , help="Path and name to file containing input parameters.")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging ")
parser.add_argument("--plot", action="store_true", help="Set this flag to plot diagnostic plots")
parser.add_argument("--show", action="store_true", help="Showing plots on screen")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
if args.plot: 
    plotvalue = 1
else:
    plotvalue = 0 # default = no plots
#-------------------------------------------------------------------------------------------------------------
#                                            LOAD INPUT PARAMETERS
#-------------------------------------------------------------------------------------------------------------
parameters        = np.genfromtxt(args.paramfile,dtype="S",comments='#',delimiter='\n')
tpfile            = parameters[0]
disp              = float(parameters[1])
outputfile        = parameters[2]
spec1D_tel        = parameters[3]
spec1Dstandard    = parameters[4]
magV_tel          = float(parameters[5])
nightdate_tel     = parameters[6]
utcstart_tel      = parameters[7]
utcstop_tel       = parameters[8]
spec1D_bri        = parameters[9]
magCAT_bri        = float(parameters[10])
nightdate_bri     = parameters[11]
utcstart_bri      = parameters[12]
utcstop_bri       = parameters[13]
spec1D_sci        = parameters[14]
nightdate_sci     = parameters[15]
utcstart_sci      = parameters[16]
utcstop_sci       = parameters[17]
throughput        = np.genfromtxt(tpfile) # reading througput into array
#-------------------------------------------------------------------------------------------------------------
#                           A:                      TELLURIC STAR
#-------------------------------------------------------------------------------------------------------------
# reading wavelengths, spectrum and ra dec from file 
wave_tel, spec1D_tel_eps, radec_tel = cal.readfitsspec1D(spec1D_tel,spec='SPEC1D_SUM')

# get E(B-V) for ra and dec
AV, extval_tel = kbs.getAv(radec_tel[0],radec_tel[1],'F098M') # filter only important when using AV

# interpolating throughput to telluric wavelengths
tpinterp_tel       = kbs.interpn(throughput[:,0]*10**4,throughput[:,1],wave_tel) 

# loading expected initial 'template' spectrum of standard star (in flux units: [erg/s/cm2/A]
standat            = np.genfromtxt(spec1Dstandard) 
wave_stan          = standat[:,0]
flux_stan          = standat[:,1]
spec_stan_flux     = kbs.interpn(wave_stan,flux_stan,wave_tel)

# rescale 'template' spec to match telluric's magnitude
magABstan          = cal.magFromSpec(wave_tel,spec_stan_flux,tpinterp_tel) # magnitude of standard in band
magABtel           = magV_tel+magABstan  # assuming telluric A0V (Vega) star => magABtel_band - magABtel_V = magABstan_band
spec_stanscale     = cal.scalespec(wave_tel,spec_stan_flux,tpinterp_tel,magABtel) # scaling standard spectrum to match telluric
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# model spectrum of a telluric standard (i.e. manipulating standard) for testing
telmodel           = 0
modelstr           = ' '
if telmodel == 1: 
    modelstr       = '(model)'
    noise          = np.random.randn(len(wave_tel))*0.01 + 0.0
    model_flux     = spec_stanscale * 0.1 + spec_stanscale * noise
    model_scale    = cal.scalespec(wave_tel,model_flux,tpinterp_tel,magABtel)
    model_eps      = cal.convertunits(wave_tel,model_scale,76000,tpinterp_tel,1.0,1.0855,conv='flux2eps')
    spec1D_tel_eps = model_eps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# correct for airmass (in [e/s])
spec1D_tel_AMCOR  = cal.correct_airmass(wave_tel,spec1D_tel_eps,radec_tel,nightdate_tel,utcstart_tel,utcstop_tel,plot=0)

# correct for galactic extinction (in [e/s])
spec1D_tel_EXTCOR = cal.correct_galacticext(wave_tel,spec1D_tel_AMCOR,extval_tel,extlaw='Cardelli')

# correct for telluric absorption and telescope sensitivity
telcorrection     = cal.correct_telluric(wave_tel,spec1D_tel_EXTCOR,spec_stanscale)   # correction factors [erg/s/cm2/A]/[e/s/pix]
spec1D_tel_TELCOR = spec1D_tel_EXTCOR * telcorrection                                 # correcting spectrum and coverting to flux
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# diagnostic plots
if args.plot:
    import pylab as plt
    from matplotlib import rc
    rc('text', usetex=True)
    rc('font', family='serif')
    plt.clf()
    f, (ax1,ax2) = plt.subplots(2) # , sharex=True
    ax1.set_title('SPEC OF TELLURIC STANDARD')
    ax1.plot(wave_tel,spec1D_tel_eps   ,'r-' ,label='Observed telluric [e/s] '+modelstr)
    ax1.plot(wave_tel,spec1D_tel_AMCOR ,'g-' ,label='Airmass corrected obs. telluric')
    ax1.plot(wave_tel,spec1D_tel_EXTCOR,'m-' ,label='Gal. Ext. + AM corrected obs. telluric')
    ax1.set_ylabel("[e/s]")
    ax1.set_xlabel("$\lambda$ / [\AA]")
    #ax1.set_yscale('log')
    ax1.legend(fancybox=True, loc='lower left')  # add the legend in the middle of the plot

    ax2.plot(wave_tel,spec_stanscale/1.0e-17   ,'y-' ,label='Scaled standard w. magAB='+str(magABtel))
    ax2.plot(wave_tel,spec1D_tel_TELCOR/1.0e-17,'b--' ,label='Telluric + AM + GE corrected obs. telluric')
    ax2.set_ylabel("[$10^{-17}$ erg/s/cm$^2$/ \AA]")
    ax2.set_xlabel("$\lambda$ / [\AA]")
    #ax2.set_yscale('log')
    ax2.legend(fancybox=True, loc='upper right')  # add the legend in the middle of the plot

    fname = outputfile.replace('.fits','_PLOTtelluric.pdf')
    print ' Wrote to '+fname
    f.savefig(fname)
    if args.show: plot.show()  # draw plot on screen
#-------------------------------------------------------------------------------------------------------------
#                           B:                       BRIGHT OBJ
#-------------------------------------------------------------------------------------------------------------
# reading wavelengths, spectrum and ra dec from file 
wave_bri, spec1D_bri_eps, radec_bri = cal.readfitsspec1D(spec1D_bri,spec='SPEC1D_SUM')

# make sure only parts of spectrum with telluric correction is included
telcorrpix     = np.where((wave_bri >= min(wave_tel)) & (wave_bri <= max(wave_tel)))[0]
wave_bri       = wave_bri[telcorrpix]
spec1D_bri_eps = spec1D_bri_eps[telcorrpix]

# get E(B-V) for ra and dec
AV, extval_bri = kbs.getAv(radec_bri[0],radec_bri[1],'F098M') # filter only important when using AV

# interpolating throughput to bright object's wavelengths
tpinterp_bri      = kbs.interpn(throughput[:,0]*10**4,throughput[:,1],wave_bri) 

# correct for airmass (in [e/s])
spec1D_bri_AMCOR  = cal.correct_airmass(wave_bri,spec1D_bri_eps,radec_bri,nightdate_bri,utcstart_bri,utcstop_bri,plot=0)

# correct for galactic extinction (in [e/s])
spec1D_bri_EXTCOR = cal.correct_galacticext(wave_bri,spec1D_bri_AMCOR,extval_bri,extlaw='Cardelli')

# correct for telluric absorption and telescope sensitivity (and go to flux units)
spec1D_bri_TELCOR = spec1D_bri_EXTCOR * telcorrection

# observed magnitude of object in band after corrections
magABbri          = cal.magFromSpec(wave_bri,spec1D_bri_TELCOR,tpinterp_bri)

# flux calibrating, i.e. scaling spectrum to match catalog value
spec1D_bri_SCALE  = cal.scalespec(wave_bri,spec1D_bri_TELCOR,tpinterp_bri,magCAT_bri)

# obtaining the fluxcalbration factors
fluxcalib         = spec1D_bri_SCALE/spec1D_bri_TELCOR

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# diagnostic plots
if args.plot:
    import pylab as plt
    from matplotlib import rc
    rc('text', usetex=True)
    rc('font', family='serif')
    plt.clf()
    f, (ax1,ax2) = plt.subplots(2) # , sharex=True
    ax1.set_title('SPEC OF FLUX CALIBRATOR')
    ax1.plot(wave_tel,spec1D_bri_eps   ,'r-' ,label='Observed bright object [e/s]')
    ax1.plot(wave_tel,spec1D_bri_AMCOR ,'g-' ,label='Airmass corrected obs. bright')
    ax1.plot(wave_tel,spec1D_bri_EXTCOR,'m-' ,label='Gal. Ext. + AM corrected obs. bright')
    ax1.set_ylabel("[e/s]")
    ax1.set_xlabel("$\lambda$ / [\AA]")
    #ax1.set_yscale('log')
    ax1.legend(fancybox=True, loc='upper right')  # add the legend in the middle of the plot

    ax2.plot(wave_tel,spec1D_bri_TELCOR/1.0e-17,'y-' ,label='Telluric + AM + GE corrected (magAB='+str(magABbri)+')')
    ax2.plot(wave_tel,spec1D_bri_SCALE/1.0e-17 ,'k-' ,label='T + AM + GE + flux calib. to magAB='+str(magCAT_bri))
    ax2.set_ylabel("[$10^{-17}$ erg/s/cm$^2$/ \AA]")
    ax2.set_xlabel("$\lambda$ / [\AA]")
    #ax2.set_yscale('log')
    ax2.legend(fancybox=True, loc='upper right')  # add the legend in the middle of the plot

    fname = outputfile.replace('.fits','_PLOTfluxcalibrator.pdf')
    print ' Wrote to '+fname
    f.savefig(fname)
    if args.show: plot.show()  # draw plot on screen
#-------------------------------------------------------------------------------------------------------------
#                           C.1:                 SCIENCE OBJ (FLUX)
#-------------------------------------------------------------------------------------------------------------
# reading wavelengths, spectrum and ra dec from file 
wave_sci, spec1D_sci_eps, radec_sci = cal.readfitsspec1D(spec1D_sci,spec='SPEC1D_SUM')

# make sure only parts of spectrum with telluric correction is included
telcorrpix     = np.where((wave_sci >= min(wave_tel)) & (wave_sci <= max(wave_tel)))[0]
wave_sci       = wave_sci[telcorrpix]
spec1D_sci_eps = spec1D_sci_eps[telcorrpix]

# get E(B-V) for ra and dec
AV, extval_sci = kbs.getAv(radec_sci[0],radec_sci[1],'F098M') # filter only important when using AV

# interpolating throughput to science object's wavelengths
tpinterp_sci      = kbs.interpn(throughput[:,0]*10**4,throughput[:,1],wave_sci) 

# correct for airmass
spec1D_sci_AMCOR  = cal.correct_airmass(wave_sci,spec1D_sci_eps,radec_sci,nightdate_sci,utcstart_sci,utcstop_sci,plot=0)

# correct for galactic extinction
spec1D_sci_EXTCOR = cal.correct_galacticext(wave_sci,spec1D_sci_AMCOR,extval_sci,extlaw='Cardelli')

# correct for telluric absorption and telescope losses
spec1D_sci_TELCOR = spec1D_sci_EXTCOR * telcorrection

# flux calibrate corrected spectrum
spec1D_sci_FCAL   = spec1D_sci_TELCOR * fluxcalib

# observed magnitude of object in band after corrections
magABsci_TELCOR   = cal.magFromSpec(wave_sci,spec1D_sci_EXTCOR,tpinterp_sci)
magABsci_FCAL     = cal.magFromSpec(wave_sci,spec1D_sci_FCAL  ,tpinterp_sci)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# diagnostic plots
if args.plot:
    import pylab as plt
    from matplotlib import rc
    rc('text', usetex=True)
    rc('font', family='serif')
    plt.clf()
    f, (ax1,ax2) = plt.subplots(2) # , sharex=True
    ax1.set_title('SPEC OF SCIENCE OBJECT')
    ax1.plot(wave_tel,spec1D_sci_eps   ,'r-' ,label='Observed science object [e/s]')
    ax1.plot(wave_tel,spec1D_sci_AMCOR ,'g-' ,label='Airmass corrected science object')
    ax1.plot(wave_tel,spec1D_sci_EXTCOR,'m-' ,label='Gal. Ext. + AM corrected science object')
    ax1.set_ylabel("[e/s]")
    ax1.set_xlabel("$\lambda$ / [\AA]")
    #ax1.set_yscale('log')
    ax1.legend(fancybox=True, loc='upper right')  # add the legend in the middle of the plot

    ax2.plot(wave_tel,spec1D_sci_TELCOR/1.0e-17,'y-' ,label='Telluric + AM + GE corrected science obj.')
    ax2.plot(wave_tel,spec1D_sci_FCAL/1.0e-17  ,'k-' ,label='T + AM + GE + flux calib. science obj.')
    ax2.set_ylabel("[$10^{-17}$ erg/s/cm$^2$/ \AA]")
    ax2.set_xlabel("$\lambda$ / [\AA]")
    #ax2.set_yscale('log')
    ax2.legend(fancybox=True, loc='upper right')  # add the legend in the middle of the plot

    fname = outputfile.replace('.fits','_PLOTscience.pdf')
    print ' Wrote to '+fname
    f.savefig(fname)
    if args.show: plot.show()  # draw plot on screen
#-------------------------------------------------------------------------------------------------------------
#                           C.2:                 SCIENCE OBJ (NOISE)
#-------------------------------------------------------------------------------------------------------------
# reading wavelengths, spectrum and ra dec from file 
wave_noi, spec1D_noi_eps, radec_noi = cal.readfitsspec1D(spec1D_sci,spec='NOISE')

# make sure only parts of spectrum with telluric correction is included
telcorrpix     = np.where((wave_noi >= min(wave_tel)) & (wave_noi <= max(wave_tel)))[0]
wave_noi       = wave_noi[telcorrpix]
spec1D_noi_eps = spec1D_noi_eps[telcorrpix]

# interpolating throughput to science object's wavelengths
tpinterp_noi      = kbs.interpn(throughput[:,0]*10**4,throughput[:,1],wave_noi) 

# correct for airmass
spec1D_noi_AMCOR  = cal.correct_airmass(wave_noi,spec1D_noi_eps,radec_noi,nightdate_sci,utcstart_sci,utcstop_sci,plot=0)

# correct for galactic exstinction
spec1D_noi_EXTCOR = cal.correct_galacticext(wave_noi,spec1D_noi_AMCOR,extval_sci,extlaw='Cardelli')

# correct for telluric absorption and telescope losses
spec1D_noi_TELCOR = spec1D_noi_EXTCOR * telcorrection

# flux calibrate corrected spectrum
spec1D_noi_FCAL   = spec1D_noi_TELCOR * fluxcalib

# spectrum containing variance, i.e., noise**2
spec1D_var_FCAL   = spec1D_noi_FCAL*spec1D_noi_FCAL
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# diagnostic plots
if args.plot:
    import pylab as plt
    from matplotlib import rc
    rc('text', usetex=True)
    rc('font', family='serif')
    plt.clf()
    f, (ax1,ax2) = plt.subplots(2) # , sharex=True
    ax1.set_title('NOISE OF SCIENCE OBJECT')
    ax1.plot(wave_tel,spec1D_noi_eps   ,'r-' ,label='Noise of science object [e/s]')
    ax1.plot(wave_tel,spec1D_noi_AMCOR ,'g-' ,label='Airmass corrected noise')
    ax1.plot(wave_tel,spec1D_noi_EXTCOR,'m-' ,label='Gal. Ext. + AM corrected noise')
    ax1.set_ylabel("[e/s]")
    ax1.set_xlabel("$\lambda$ / [\AA]")
    #ax1.set_yscale('log')
    ax1.legend(fancybox=True, loc='upper right')  # add the legend in the middle of the plot

    ax2.plot(wave_tel,spec1D_noi_TELCOR/1.0e-17,'y-' ,label='Telluric + AM + GE corrected noise')
    ax2.plot(wave_tel,spec1D_noi_FCAL/1.0e-17  ,'k-' ,label='T + AM + GE + flux calib. noise')
    ax2.set_ylabel("[$10^{-17}$ erg/s/cm$^2$/ \AA]")
    ax2.set_xlabel("$\lambda$ / [\AA]")
    #ax2.set_yscale('log')
    ax2.legend(fancybox=True, loc='upper right')  # add the legend in the middle of the plot

    fname = outputfile.replace('.fits','_PLOTsciencenoise.pdf')
    print ' Wrote to '+fname
    f.savefig(fname)
    if args.show: plot.show()  # draw plot on screen
#-------------------------------------------------------------------------------------------------------------
#                                               SAVE OUTPUT
#-------------------------------------------------------------------------------------------------------------
# writing output to fits table
if args.verbose: print 'Writing calibrated 1D spectrum to '+outputfile

col1     = pyfits.Column(name='WAVELENGTH', format='F', array=wave_sci)
col2     = pyfits.Column(name='SCI_SPEC'  , format='F', array=spec1D_sci_FCAL)
col3     = pyfits.Column(name='SCI_NOISE' , format='F', array=spec1D_noi_FCAL)
col4     = pyfits.Column(name='SCI_VAR'   , format='F', array=spec1D_var_FCAL)
cols     = pyfits.ColDefs([col1, col2, col3, col4])

tbhdu    = pyfits.new_table(cols)          # creating table header

# writing hdrkeys:   '---KEY--',                  '----------------MAX LENGTH COMMENT-------------'
tbhdu.header.append(('PARFILE ' ,args.paramfile   ,'Paramter file used with performSpecCal.py'),end=True)
tbhdu.header.append(('RA      ' ,radec_sci[0]     ,'RA  of science object (mask)'),end=True)
tbhdu.header.append(('DEC     ' ,radec_sci[1]     ,'DEC of science object (mask)'),end=True)
tbhdu.header.append(('RA_TEL  ' ,radec_tel[0]     ,'RA  of telluric standard (mask)'),end=True)
tbhdu.header.append(('DEC_TEL ' ,radec_tel[1]     ,'DEC of telluric standard (mask)'),end=True)
tbhdu.header.append(('RA_BRI  ' ,radec_bri[0]     ,'RA  of flux calibrator (bright object) (mask)'),end=True)
tbhdu.header.append(('DEC_BRI ' ,radec_bri[1]     ,'DEC of flux calibrator (bright object) (mask)'),end=True)

hdu      = pyfits.PrimaryHDU()             # creating primary (minimal) header
thdulist = pyfits.HDUList([hdu, tbhdu])    # combine primary and table header to hdulist
thdulist.writeto(outputfile,clobber=True)  # write fits file (clobber=True overwrites excisting file)




#-------------------------------------------------------------------------------------------------------------
# calculate magnitude(s)
testmag = 0
if testmag == 1:
    magABspec_corrEXT   =  cal.magFromSpec(wave_sci,spec1D_corrEXT,tpinterp_sci)
    magABspec           =  cal.magFromSpec(wave_sci,spec1D_sci_flux,tpinterp_sci)
    magABspec_stanscale =  cal.magFromSpec(wave_sci,spec_stanscale,tpinterp_sci)
    print 'magABspec_corrEXT,magABspec,magABstan',magABspec_corrEXT,magABspec,magABstan
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

