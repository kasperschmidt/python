#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# performSpecCal_createParamFile.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Script for creating a template performSpecCal_*.param file needed to run
# performSpecCal.py for a single object
#----------------------------
#   COMMENTS
#----------------------------
# For multiple runs see e.g. script
# /Users/kasperborelloschmidt/work/observing/130423_MOSFIRE_BoRG/data/130606reduction/paramcommands.sh
#----------------------------
#   INPUTS:
#----------------------------
# band             : band of observations: 'Y' or 'J'
# scienceobject    : fits file containing 1D spectrum of science object from 
#                    extractMOSFIRE1Dsignal2noiseSpec.py
# fluxcalibrator   : fits file containing 1D spectrum of bright flux calibrator from 
#                    extractMOSFIRE1Dsignal2noiseSpec.py
# magFC            : Extinction corrected catalog magnitude of object in 'band' 
#                    (corrected for slitlosses e.g. via estimateslitloss.py)
# telluric         : fits file containing 1D spectrum of telluric standard from 
#                    extractMOSFIRE1Dsignal2noiseSpec.py
# magVtel          : V band magnitude of telluric. See e.g. 
#                    http://www1.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/astrocat/hipparcos
# utcstartsci      : Star of observations in UT for Science object (and fluxcalibrator). Format: 2013-04-26 09:00
# utcstopsci       : End of observations in UT for Science object (and fluxcalibrator). Format: 2013-04-26 13:00
# utcstarttel      : Star of observations in UT for tellutic. Format: 2013-04-26 09:00
# utcstoptel       : End of observations in UT for telluric. Format: 2013-04-26 13:00
# paramfile        : name to save parameter file to
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# paramfile to be run with performSpecCal.py
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x performSpecCal_createParamFile.py       (only required once)
# bash> performSpecCal_createParamFile.py Y borg_1437+5043_130424_Y_borg_1437+5043_0560_r2_T12c_eps.fits borg_1437+5043_130424_Y_305_r1_eps.fits 99.9 rectified_HIP53735_130424_HIP53735_Y_A-B_KBSmod_eps.fits 9.9 '2013-04-25 07:00'   '2013-04-25 08:00'       '2013-04-25 08:10'   '2013-04-25 08:20' performSpecCal_this_is_output.param --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-06-06  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import numpy as np         # enable opening with genfromtxt
import pdb                 # for debugging with pdb.set_trace()
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("band", type=str, help="Band of observations")
parser.add_argument("scienceobject", type=str, help="Science object filename")
parser.add_argument("fluxcalibrator", type=str, help="Flux calibrator filename")
parser.add_argument("magFC", type=float, help="Magnitude of flux calibrator in 'band'")
parser.add_argument("telluric", type=str, help="Telluric standard filename")
parser.add_argument("magVtel", type=float, help="Magnitude in V band of Telluric")
parser.add_argument("utcstartsci", type=str, help="UT start of science observations")
parser.add_argument("utcstopsci", type=str, help="UT end of science observations")
parser.add_argument("utcstarttel", type=str, help="UT start of telluric observations")
parser.add_argument("utcstoptel", type=str, help="UT end of telluric observations")
parser.add_argument("paramfile", type=str, help="Name of outputfile (*.param)")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

disp          = {'J':1.3028, 'Y':1.0855}
outputcalspec = args.scienceobject.replace('.fits','_CALIBRATED.fits')
cattel        = '/Users/kasperborelloschmidt/work/observing/calibrationstars/alpha_lyr_stis_005.ascii'
nightdate0    = args.scienceobject.split('_')[2]
nightdate     = '20'+nightdate0[0:2]+'/'+nightdate0[2:4]+'/'+nightdate0[4:6]

inputlist = ('/Users/kasperborelloschmidt/work/MOSFIRE/XTcalc_dir/mosfire/'+args.band+'_tp_tot.txt',
            disp[args.band],
            outputcalspec,
            args.telluric,
            cattel,
            args.magVtel,
            nightdate,
            args.utcstarttel,
            args.utcstoptel,
            args.fluxcalibrator,
            args.magFC,
            nightdate,
            args.utcstartsci,
            args.utcstopsci,
            args.scienceobject,
            nightdate,
            args.utcstartsci,
            args.utcstopsci)

#-------------------------------------------------------------------------------------------------------------
output = open(args.paramfile,"w")

output.write('''
# Input parameters and spectra to run performSpecCal.py which flux calibrates
# and corrects 1D (MOSFIRE) spectra from telluric absorption, telescope losses,
# and galactic extinction.
# The following input are needed (in that order!)
#
# ---- GENERAL INFO ----
# tpfile            : File containing the throughput of filter spectrum is taken in
#                     Expects file with columns: wavelength[1e-6m], througput
%s

# Dipersion of filter [A/pix]. See e.g. McLean et al. 2012 (J = 1.3028, Y = 1.0855
%s

# name of output fits file to contain calibrated spectrum of object(s) etc.
%s

# ---- INFO FOR TELLURIC STANDARD ----
# spec1D_tel       : name of file containing the 1D spectrum of the telluric standard
#                    File can be created with extractMOSFIRE1Dsignal2noiseSpec.py
%s

# spec1Dstandard   : name of file containing the 'template' spectrum of telluric type, e.g. Vega (A0V) spec
#                    Expects file with two columns: wave[A] , flux[erg/s/cm2/A]
%s

# magV_tel         : catalog V band magnitude of telluric
#                    for HIP stars see http://www1.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/astrocat/hipparcos
%s

# nightdate_tel    : The date of observation night to calculate airmass correction for. String of type '2013/04/25'
%s

# utcstart_tel     : The start of the observations in UTC. String of type 'YYYY-mm-dd HH:MM'. 
#                    UT times are in headers (TIME-OBS and TIME-END keywords) of raw fits file
#                    To obtain this from a local time and date use 
#                       bash> skycalc_java 
#                    or do: 
#                       >>>   time = datetime.datetime.strptime('2013-04-25 23:00', '$Y-$m-$d $H:$M') (if summertime subtract 1 hour!)
#                       >>>   time = time.replace(tzinfo=pytz.timezone('US/Aleutian'))
#                       >>>   time.astimezone(pytz.timezone('utc'))
%s

# utcstop_tel      : Similar to utcstart but for end of the observations in UTC. String of type 'YYYY-mm-dd HH:MM'. 
%s

# ---- INFO FOR BRIGHT OBJECT (FLUX CALIBRATOR) ----
# spec1D_bright    : name of file containing the 1D spectrum to flux calibrate with
#                    File can be created with extractMOSFIRE1Dsignal2noiseSpec.py
%s

# mag_bright       : expected (catalog) AB magnitude of object in spectral band, i.e. Y_AB for Y-band spectroscopy
%s

# nightdate_bright : The date of observation night to calculate airmass correction for. String of type '2013/04/25'
%s

# utcstart_bright  : The start of the observations in UTC. String of type 'YYYY-mm-dd HH:MM'. 
%s

# utcstop_bright   : Similar to utcstart but for end of the observations in UTC. String of type 'YYYY-mm-dd HH:MM'. 
%s

# ---- INFO FOR (FAINT) OBJECT OF INTEREST ----
# spec1D_sci       : name of file containing the 1D science spectrum to calibrate (including ivar map)
#                    File can be created with extractMOSFIRE1Dsignal2noiseSpec.py
%s

# nightdate_sci    : The date of observation night to calculate airmass correction for. String of type '2013/04/25'
%s

# utcstart_sci     : The start of the observations in UTC. String of type 'YYYY-mm-dd HH:MM'. 
%s

# utcstop_sci      : Similar to utcstart but for end of the observations in UTC. String of type 'YYYY-mm-dd HH:MM'. 
%s
''' % inputlist)

output.close()
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------


