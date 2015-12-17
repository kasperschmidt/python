#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# createBoRG1437mask.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Creating a mask of BoRG1437 by hand representing the different
# exposure times of different regions of the composite image. 
#----------------------------
#   COMMENTS
#----------------------------
# 
#----------------------------
#   INPUTS:
#----------------------------
# 
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --verbose        : set -verbose to get info/messages printed to the screen
# --show           : showning plots on screen
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# BoRG1437mask    : The created fits image mask
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x createBoRG1437mask.py       (only required once)
# bash> createBoRG1437mask.py --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-04-18  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs #
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import numpy as np         # enable opening with genfromtxt
import pdb                 # for debugging with pdb.set_trace()
import pyfits              #
#-------------------------------------------------------------------------------------------------------------
def acoef(x1,x2,y1,y2): # slope of line
    a = (y2-y1) / (x2-x1)
    return a
def bcoef(x1,x2,y1,y2): # intersection of line with y axis
    b = y1 - (y2-y1) / (x2-x1) * x1
    return b
def yvalue(x,list): # intersection of line with y axis
    y = list[1]*x+list[2]
    return y
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
#parser.add_argument("catlist", type=str, help="Provide path and name of file containing list of catalog names")
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
# The coordinates of the 10 points defining 'intersections' or corners of the exposures
xcoord = [17.,808.,1465.,687.,1207.,2268.,1788.,1389.,2896.,2156.]
ycoord = [1223.,2525.,362.,1986.,2286.,1628.,852.,638.,1473.,2805.]
#-------------------------------------------------------------------------------------------------------------
# The coefficients of the 8 lines defining the FOVs
lines = [['lineA0',acoef(xcoord[0],xcoord[2],ycoord[0],ycoord[2]),bcoef(xcoord[0],xcoord[2],ycoord[0],ycoord[2])],
         ['lineB1',acoef(xcoord[0],xcoord[1],ycoord[0],ycoord[1]),bcoef(xcoord[0],xcoord[1],ycoord[0],ycoord[1])],
         ['lineC2',acoef(xcoord[1],xcoord[5],ycoord[1],ycoord[5]),bcoef(xcoord[1],xcoord[5],ycoord[1],ycoord[5])],
         ['lineD3',acoef(xcoord[2],xcoord[5],ycoord[2],ycoord[5]),bcoef(xcoord[2],xcoord[5],ycoord[2],ycoord[5])],
         ['lineE4',acoef(xcoord[3],xcoord[7],ycoord[3],ycoord[7]),bcoef(xcoord[3],xcoord[7],ycoord[3],ycoord[7])],
         ['lineF5',acoef(xcoord[3],xcoord[9],ycoord[3],ycoord[9]),bcoef(xcoord[3],xcoord[9],ycoord[3],ycoord[9])],
         ['lineG6',acoef(xcoord[7],xcoord[8],ycoord[7],ycoord[8]),bcoef(xcoord[7],xcoord[8],ycoord[7],ycoord[8])],
         ['lineH7',acoef(xcoord[8],xcoord[9],ycoord[8],ycoord[9]),bcoef(xcoord[8],xcoord[9],ycoord[8],ycoord[9])]]
#-------------------------------------------------------------------------------------------------------------
# Creating mask array to fill with values
dimension = (3338, 3302)
mask      = np.zeros(dimension)
#-------------------------------------------------------------------------------------------------------------
# Defining BoRG_1437+5043 region 1 (see notes from 130421)
for ii in range(dimension[0]):     # looping over y-values
    for jj in range(dimension[1]): # looping over x-values
        x = jj+1
        y = ii+1
        if (y > yvalue(x,lines[0])) and (y > yvalue(x,lines[3])) and (y < yvalue(x,lines[1])) and (y < yvalue(x,lines[4])):
            mask[ii,jj] = 1.0
        if (y > yvalue(x,lines[4])) and (y < yvalue(x,lines[1])) and (y < yvalue(x,lines[2])) and (y > yvalue(x,lines[5])):
            mask[ii,jj] = 1.0
        if (y > yvalue(x,lines[4])) and (y < yvalue(x,lines[6])) and (y > yvalue(x,lines[3])):
            mask[ii,jj] = 1.0
#-------------------------------------------------------------------------------------------------------------
# Defining BoRG_1437+5043 region 2 (see notes from 130421)
for ii in range(dimension[0]):     # looping over y-values
    for jj in range(dimension[1]): # looping over x-values
        x = jj+1
        y = ii+1
        if (y > yvalue(x,lines[6])) and (y > yvalue(x,lines[3])) and (y > yvalue(x,lines[4])) and (y < yvalue(x,lines[5])) and (y < yvalue(x,lines[2])):
            mask[ii,jj] = 2.0
#-------------------------------------------------------------------------------------------------------------
# Defining BoRG_1437+5043 region 3 (see notes from 130421)
for ii in range(dimension[0]):     # looping over y-values
    for jj in range(dimension[1]): # looping over x-values
        x = jj+1
        y = ii+1
        if (y > yvalue(x,lines[6])) and (y < yvalue(x,lines[3])) and (y < yvalue(x,lines[2])):
            mask[ii,jj] = 3.0
        if (y < yvalue(x,lines[5])) and (y < yvalue(x,lines[7])) and (y > yvalue(x,lines[2])) and (y > yvalue(x,lines[6])):
            mask[ii,jj] = 3.0
#-------------------------------------------------------------------------------------------------------------
# write new image to file
filename = 'BoRG_1437+5043_byhandexposuremask.fits'
if args.verbose: print 'Mask of BoRG+1437+5043 will be written to '+filename
pyfits.writeto(filename, mask, clobber=True)

#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

