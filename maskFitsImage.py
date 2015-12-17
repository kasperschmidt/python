#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# maskFitsImage.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Masking a fits image by setting values to 0 from a mask and a given
# value in that mask. E.g. a map of exposure times can be used to create
# a new fits image containing pixels with that exposure time. Suitable
# for mosaics of images.
#----------------------------
#   COMMENTS
#---------------------------- 
#
#----------------------------
#   INPUTS:
#----------------------------
# imagefits        : Path AND name of image to create new masked image from
# maskfits         : Path AND name of (exposure time) mask create us for new image
# maskvalue        : list of values or range of values to include from mask when
#                    creating new image. Accepted formats are:
#                        '(minvval,maxval)'       # include all values between minval and maxval (both included)
#                        '[val1,val2,...,valN]'   # include individual values in list
#                    Note that the parenthesis/brackets are used to determine whether 
#                    the input is a range or a list of values.
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --verbose        : set --verbose to get info/messages printed to the screen
# --show           : showning plots on screen
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# outputimage      : New fits image where all pixels not satisfying mask criteria are set to 0
#                    The size of the image is the same as fitsimage and maskimage
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x maskFitsImage.py       (only required once)
# bash> maskFitsImage.py borg_1437+5043_f125w_wfc3ir_drz.fits borg_1437+5043_f125w_wfc3ir_drz_wht_MASK_smooth.fits '(42000,70000)' --verbose 

# bash> maskFitsImage.py borg_1437+5043_f125w_wfc3ir_drz.fits ../../BoRG_1437+5043_byhandexposuremask.fits '[1]' --verbose 
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-04-16  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs #
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import numpy as np         # enable opening with genfromtxt
import pdb                 # for debugging with pdb.set_trace()
import pyfits
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("imagefits", type=str, help="Path and name of file containing fist image to mask")
parser.add_argument("maskfits", type=str, help="Path and name of mask to use")
parser.add_argument("maskvalue", type=str, help="Values in mask to extract")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--show", action="store_true", help="Showing plots on screen")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# getting mask values 
mval      = map(float, args.maskvalue[1:-1].split(','))
Nmval     = len(mval)
if args.maskvalue[0] == '(':
    valtype = 'range'
    if Nmval != 2: sys.exit('ERROR: Range of mask values "(...)" must contain 2 values only. Contains '+str(Nmval)+' --> ABORTING')
if args.maskvalue[0] == '[':
    valtype = 'list'

if args.verbose: print 'The mask values to extract data for was a '+valtype+' of values.'
#-------------------------------------------------------------------------------------------------------------
# reading fits files
img       = pyfits.getdata(args.imagefits,ext=0)  # reading data into array   [rows,columns]
sizeimg   = img.shape
hduimg    = pyfits.open(args.imagefits)           # Load the FITS hdulist using pyfits
hdrimg    = hduimg[0].header                      # get header

mask      = pyfits.getdata(args.maskfits,ext=0)   # reading data into array   [rows,columns]
sizemask  = mask.shape
hdumask   = pyfits.open(args.maskfits)            # Load the FITS hdulist using pyfits
hdrmask   = hdumask[0].header                     # get header

if sizeimg != sizemask:
    sys.exit('ERROR: The size of the input image '+str(sizeimg)+' and the size of the mask '+str(sizemask)+' do no agree --> ABORTING')

#-------------------------------------------------------------------------------------------------------------
# masking image if a range is given
if valtype == 'range':
    newimage = img
    if '_rms' in imagefits.split('/')[-1]:
        newimage[np.where(mask < mval[0])] = np.nan
        newimage[np.where(mask > mval[1])] = np.nan
    else:
        newimage[np.where(mask < mval[0])] = 0.0
        newimage[np.where(mask > mval[1])] = 0.0
#-------------------------------------------------------------------------------------------------------------
# masking image if a list is given
if valtype == 'list':
    if '_rms' in args.imagefits.split('/')[-1]: # makin nan array if image is rms map
        newimage = img*np.nan
    else:
        newimage = img*0.0

    for ii in range(Nmval):
        newimage[np.where(mask == mval[ii])] = img[np.where(mask == mval[ii])]
#-------------------------------------------------------------------------------------------------------------
# write new image to file
vallist  = args.maskvalue[1:-1].replace('.','p').replace(',','_')
filename = args.imagefits.replace('.fit','_MASK'+valtype+'_'+vallist+'.fit')
if args.verbose: print 'New image will be written to '+filename
pyfits.writeto(filename, newimage, hdrimg, clobber=True)
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

