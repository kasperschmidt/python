#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# getPixvals4BoRGobject.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# small script returning pixel values in the *_rms.fits, *_drz.fits and *_segm.fits files for a given BoRG id
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# rmsfits          : Path AND basename of fits files, i.e., file name without *_rms.fits, *_drz.fits and *_segm.fits
# drzfits          :
# segmfits         :
# objectid         : ID of object to print pixel values for
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
#
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x getPixvals4BoRGobject.py       (only required once)
# bash> getPixvals4BoRGobject.py '/Users/kasperborelloschmidt/Desktop/borg_0637-7518/borg_0637-7518_f125w_wfc3ir_rms_norm.fits' '/Users/kasperborelloschmidt/Desktop/borg_0637-7518/borg_0637-7518_f606w_wfc3uvis_drz.fits' '/Users/kasperborelloschmidt/Desktop/borg_0637-7518/borg_0637-7518_final_sources_F125_segm.fits' 8 --verbose

# bash> getPixvals4BoRGobject.py '/Users/kasperborelloschmidt/Desktop/borg_0637-7518/borg_0637-7518_f606w_wfc3uvis_rms.fits' '/Users/kasperborelloschmidt/Desktop/borg_0637-7518/borg_0637-7518_f606w_wfc3uvis_drz.fits' '/Users/kasperborelloschmidt/Desktop/borg_0637-7518/borg_0637-7518_final_sources_F125_segm.fits' 8 --verbose

#
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-11-05  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import subprocess # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
import numpy as np# enable opening with genfromtxt
import pyfits
#----------------------------
#   FUNCTIONS
#----------------------------
def pathAname(str):                         # splitting string with path and name in to
    strsplit = str.split('/')               # splitting string
    name     = strsplit[-1]                 # saving filename (last entry of catsplit)
    slash    = '/'
    path=slash.join(strsplit[0:-1])+'/'     # putting path back together
    return [path,name]
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("rmsfits",  type=str, help="RMS fits")
parser.add_argument("drzfits",  type=str, help="DRZ fits")
parser.add_argument("segmfits", type=str, help="SEGM fits")
parser.add_argument("objectid", type=int, help="Object ID to look up in segmentation map")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")
parser.add_argument("-k", "--keywords", type=str, help="Provide list of keywords to run BPZ with in string")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# reading fits files
RMS     = pyfits.open(args.rmsfits)
RMSimg  = RMS[0].data

DRZ     = pyfits.open(args.drzfits)
DRZimg  = DRZ[0].data

SEGM    = pyfits.open(args.segmfits)
SEGMimg = SEGM[0].data

#print SEGMimg.shape, args.objectid
indices = np.where(SEGMimg==args.objectid)
xvals   = indices[1]
yvals   = indices[0]


np.set_printoptions(threshold=np.nan)  # to enable printing all elements for sub-array
print ':: '+sys.argv[0]+' :: Found '+str(SEGMimg[indices].size)+' pixels in the segmentation map for objects '+str(args.objectid)
#print RMSimg[indices].size
#print DRZimg[indices].size
#print SEGMimg[indices].size
print ' '
print ':: '+sys.argv[0]+' :: the RMS values:'
print '   MIN and MAX values are: '+str(min(RMSimg[indices]))+' and '+str(max((RMSimg[indices])))
ent = np.where(RMSimg[indices]==max(RMSimg[indices]))
ent[0]
xent = xvals[ent[0]]
yent = yvals[ent[0]]
print '   the pixel coordinates of max values: ('+str(xent[0]+1.0)+','+str(yent[0]+1.0)+')'

#print ' '
#print ':: '+sys.argv[0]+' :: the DRZ values:'
#print DRZimg[indices]

#print ' '
#print ':: '+sys.argv[0]+' :: the SEGM values:'
#print SEGMimg[indices]
#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

