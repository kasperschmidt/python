#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# createBoRGcutouts.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Creating cutouts/postage stamps of borg objects
#----------------------------
#   COMMENTS
#----------------------------
# To run for multiple objects one can creating a script similar to 
# /Users/kasperborelloschmidt/work/BoRG/borgdata/createcutout_BoRG_higzcandidates.sh 
#----------------------------
#   INPUTS:
#----------------------------
# imagedir         : Path to the directory containing the drz fits images
# objname          : Name (of object) used when naming/saving cutouts.
# coords           : The coordinates of the object to create cutout for and 
#                    the coordinate type. Either coordinates are given in
#                    pixel number or in ra dec (degrees). Hence, 3 values
#                    are expected. For example:
#                        1721.28 703.124 pixels
#                        157.7353387 38.0673682 wcs
# pswidth          : The width of the postage stamps in degrees. 
#                    Final size of cutout is 2*postagewidth X 2*postagewidth
# outputdir        : directory to put postage stamps in (already excists)
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
# bash> chmod +x createBoRGcutouts.py       (only required once)
# bash> createBoRGcutouts.py version_2p0/borg_0440-5244/ borg_0440-5244_682_pix 1488.73 1113.19 pixels 0.00138 BoRG_higzcandidates_ALL_pstamps130206 --verbose

# bash> createBoRGcutouts.py version_2p0/borg_0440-5244/ borg_0440-5244_682_wcs 69.9455842724 -52.7320161788 wcs 0.00138 BoRG_higzcandidates_ALL_pstamps130206 --verbose


# bash> createBoRGcutouts.py ../BoRG_1437+5043_followup/cat_byhandmasking130418/region3/borg_1437+5043/ borg_1437+5043_region3_0447 219.18983 50.73406 wcs 0.00138 /Users/kasperborelloschmidt/Desktop/cutouts/ --verbose


# Hint... to ease running multiple objects create a shell script like BoRG/borgdata/createcutout_BoRG_higzcandidates.sh
# 
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-02-06  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs #
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import commands            # get output from spawned command line processes
import getopt              # used to extract/obtain the optional input
import numpy as np         # enable opening with genfromtxt
import pdb                 # for debugging
import cutout              # 
import pyfits              #
import pywcs               #
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("imagedir", type=str, help="Directory containing drz fits images to cutout from")
parser.add_argument("objname", type=str, help="Name (of objects) to use when creating cutouts")
parser.add_argument("coords", type=str, nargs=3, help="The coordinates to use:  xval yval pixels  OR  ra dec wcs")
parser.add_argument("pswidth", type=float, help="Size of postage stamps in degrees")
parser.add_argument("outputdir", type=str, help="Directory to save images in (already created)")
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
# finding images
drzimages = commands.getoutput('ls '+args.imagedir+'/*drz*fits').split('\n')
Nimg = len(drzimages)
if args.verbose: print 'Found ',Nimg,' *drz*.fits images in '+args.imagedir+' to create cutouts from'
#-------------------------------------------------------------------------------------------------------------
# creating cutouts
for ii in range(Nimg): # looping over images
    image   = pyfits.getdata(drzimages[ii],ext=0)  # reading image into array   [rows,columns]
    hdulist = pyfits.open(drzimages[ii])     # Load the FITS hdulist using pyfits
    hdr     = hdulist[0].header                  # get header

    if args.coords[2] == 'pixels':
        wcs = pywcs.WCS(hdr)  # Extract wcs (coordinate) information
        pix = np.array([[args.coords[0],args.coords[1]]], np.float_)
        sky = wcs.wcs_pix2sky(pix, 1)  # convert to 'sky' values. 2nd arg. is "origin": 1 = a 1-based (Fortran-like) coordinates.
        ra  = sky[0][0]
        dec = sky[0][1]
    elif args.coords[2] == 'wcs':
        ra  = float(args.coords[0])
        dec = float(args.coords[1])
    else:
        sys.exit('Coordinate units '+args.coords[2]+' not allowed --> ABORTING')

    outname = args.outputdir+'/'+args.objname+'_cutoutFROM_'+drzimages[ii].split('/')[-1]
#    pdb.set_trace()
    cutout.cutout(drzimages[ii],ra,dec,xw=args.pswidth,yw=args.pswidth,units='wcs',outfile=outname)
    if args.verbose: print 'Created the cutout '+outname
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

