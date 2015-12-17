#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# calculatearea4fitsimage.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Using the header innformation of a fits image to calculate the area on the sky
#----------------------------
#   COMMENTS
#----------------------------
# 
#----------------------------
#   INPUTS:
#----------------------------
# fitsiimage       : Name of fits image to calculate area of.
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --exclude        : exclude pixels with the provided value (e.g. 'nan', '0', '2.5', '>0' or '<0')
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
# bash> calculatearea4fitsimage.py /Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p0/borg_0440-5244/borg_0440-5244_f125w_wfc3ir_drz.fits --verbose --exclude nan
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-08-17  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs #
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import getopt              # used to extract/obtain the optional input
import numpy as np         # enable opening with genfromtxt
import pdb                 # for debugging with pdb.set_trace()
import pyfits
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("fitsimage", type=str, help="Image to calculate area of.")
# ---- optional arguments ----
parser.add_argument("--exclude", type=str, help="Pixel values to exclude from estimate (e.g. 'nan', '0', '2.5', '>0' or '<0')")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
# reading fits iamge
if args.verbose: print ' - Loading ',args.fitsimage
imgHDU = pyfits.open(args.fitsimage)
img    = imgHDU[0].data
imghdr = imgHDU[0].header  
simg   = img.shape
if args.verbose: print ' - shape of fits image is       ',simg
#-------------------------------------------------------------------------------------------------------------
cd11      = imghdr['CD1_1']
cd21      = imghdr['CD2_1']
pixscaleX = np.sqrt(np.power(cd11,2)+np.power(cd21,2)) * 3600. 

cd22      = imghdr['CD2_2']
cd12      = imghdr['CD1_2']
pixscaleY = np.sqrt(np.power(cd22,2)+np.power(cd12,2)) * 3600. 

apixel   = pixscaleX*pixscaleY
if args.verbose: print ' - size (area) of each pixel is ',apixel,'square arc seconds'

Npix = float(simg[0]*simg[1]) # Number of pixels in image

if args.exclude == 'nan':
    Npix = len(np.where(img != np.nan)[0])
elif args.exclude == '<0':
    Npix = len(np.where(img < 0)[0])
elif args.exclude == '>0':
    Npix = len(np.where(img > 0)[0])
elif args.exclude:
    Npix = len(np.where(img != float(args.exclude))[0])
if args.verbose: print ' - Pixels to calculate area for ',Npix,' (out of max',simg[0]*simg[1],')' 

size = Npix*apixel # total size in square arc seconds
if args.verbose: 
    print ' - Hence the area of ',args.fitsimage.split('/')[-1],' is ',str("%.4f" % (size/3600.)),' square arc minutes'
    print '                                ',str("%.4f" % (size)),' square arc seconds'
    print '                                ',str("%.8f" % (size/3600./3600.)),' square degrees'

#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------

