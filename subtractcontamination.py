#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# subtractcontamination.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# loading 2D.fits image extracted via the 3D-HST pipeline and a
# fits image with the contamination subtracted SCI extension
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# twodfits         : the 2D fits file containing the grism to subtract contamination from
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
#
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> subtractcontamination.py 'MACS0717.5+3745-020-full140121_00846-G102.2D.fits' --verbose
# bash> subtractcontamination.py 'MACS0717.5+3745-020-full140121_02068-G102.2D.fits' --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2014-01-24  started by K. B. Schmidt (STScI)
#----------------------------
#   MODULES
#----------------------------
import argparse            # argument managing
import sys                 # enabling arguments to code
import os
import pdb                 # for debugging with pdb.set_trace()
import pyfits
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("twodfits", type=str, help="2D fits file to subtract contamination from")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print ' -  Loading ',args.twodfits
hduimg  = pyfits.open(args.twodfits) # Load the FITS hdulist
hdrsci  = hduimg[4].header    # extracting science header
sci     = hduimg[4].data
contam  = hduimg[7].data

filename, fileext =  os.path.splitext(args.twodfits)
output = filename+'_ContamSub'+fileext
pyfits.writeto(output, sci-contam, hdrsci, clobber=False)
if args.verbose: print ' -  Wrote ',output
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------