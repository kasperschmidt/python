#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# addsources2fits.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Adds sources to a 'detection image' to be used with sextractor forcing it to perform
# photometry at the specified location.
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# sources              : list of sources to add. Expects the format:
#                            'x1,y1,type1,x2,y2,type2,...,xN,yN,typeN'
#                        where x and y are the coordinates and type can be either 
#                        'pixel' (not enabled as of 130909) or 'wcs' [in degrees]
# fitsimage            : fits image to add sources to
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --verbose            : set -verbose to get info/messages printed to the screen
# --stop               : stoppping program before end for de-bugging
# --help               : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
#
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> addsources2fits.py '219.22405,50.72597,wcs,219.22027,50.71563,wcs' /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/borg_1437+5043/borg_1437+5043_f125w_wfc3ir_drz.fits  --verbose 
#
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-09-09  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import commands            # get output from spawned command line processes
import getopt              # used to extract/obtain the optional input
import numpy as np         # enable opening with genfromtxt
import pyfits              #
import pywcs               #
import pdb                 # for debugging with pdb.set_trace()
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("sources", type=str, help="Sources to add to fits image")
parser.add_argument("fitsimage", type=str, help="Fits image to add sources to")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
# getting source info
strsplit = np.asarray(args.sources.split(','))
Nsources = len(strsplit)/3.0
if Nsources != int(len(strsplit)/3.0):
    sys.exit('The string of sources provided does not seem to have the right number of arguments --> ABORTING')
xcoord   = strsplit[list(np.arange(0,Nsources*3,3))].astype(np.float)
ycoord   = strsplit[list(np.arange(1,Nsources*3,3))].astype(np.float)
types    = strsplit[list(np.arange(2,Nsources*3,3))]
#-------------------------------------------------------------------------------------------------------------
# reading fits image
imgHDU = pyfits.open(args.fitsimage)
img    = imgHDU[0].data
imghdr = imgHDU[0].header  
simg   = img.shape
if args.verbose: print ' - shape of fits image is       ',simg
#-------------------------------------------------------------------------------------------------------------
# Getting pixel position for given ra and dec
xpix, ypix = [],[]
for ss in xrange(int(Nsources)):
    if types[ss] == 'wcs':
        wcs   = pywcs.WCS(imghdr)  # Extract wcs (coordinate) information
        radec = np.array([[xcoord[ss],ycoord[ss]]], np.float_)
        pix   = wcs.wcs_sky2pix(radec, 1)  # convert ra and dec to pixel
        xpix.append(pix[0][0])
        ypix.append(pix[0][1])      
    elif types[ss] == 'pixel':
        xpix.append(xcoord[ss])
        ypix.append(ycoord[ss])
    else:
        sys.exit('Invalid type of coordiante: '+str(types[ss])+' Valid types are "wcs" and "pixel" --> ABORTING')
#-------------------------------------------------------------------------------------------------------------
# adding sources to new image
imgsourceadd = img.copy()
#imgsourceadd = np.array(simg) # empty array 

#KBS turn the following values into keywords so they are easy changable
n_dx, n_dy = 15.0, 15.0
intensity  = 0.5

for ii in xrange(int(Nsources)):
   x_int = np.round(xpix[ii])
   y_int = np.round(ypix[ii])
   for dx in np.arange(-n_dx,n_dx):
        for dy in np.arange(-n_dy,n_dy):
            dist2 = (xpix[ii]-(x_int+dx))**2.0 + (ypix[ii]-(y_int+dy))**2.0
            flux = intensity/(2.0+dist2)
            imgsourceadd[y_int+dy,x_int+dx] = flux + imgsourceadd[y_int+dy,x_int+dx] 
#-------------------------------------------------------------------------------------------------------------
outname = args.fitsimage.replace('.fits','_sourcesadded.fits')
hdu = pyfits.PrimaryHDU(imgsourceadd,header=imghdr)
hdu.writeto(outname,clobber=False)
if args.verbose: print ' - Wrote image with added sources to '+outname
#-------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
