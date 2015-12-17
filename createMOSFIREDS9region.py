#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# createMOSFIREDS9region.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Creating a DS9 region file resembling the mosfire layout around a
# given position. Convenient for finding charts.
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# RADEC            : RA and DEC to drac region file around (center of MOSFIRE field)
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --PA             : postion angle (N-->E) in degrees to rotate regions with. 
#                    Default value is 0.0
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# ds9 region file  : Region file with the name MOSFIREds9_RAxxxDECyyyPAzzz.reg put in working directory
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x createMOSFIREDS9region.py       (only required once)
# bash> createMOSFIREDS9region.py 117.71050  29.28142 --PA 59.0 --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-12-13  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import commands   # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
import numpy as np# enable opening with genfromtxt
#----------------------------
#   FUNCTIONS
#----------------------------
#
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("RADEC", type=float, nargs=2, help="RA DEC [deg] to draw MOSFIRE around")
# ---- optional arguments ----
parser.add_argument("--PA", type=float, help="Position angle N-->E [deg] Default value is 0.0")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
if args.PA:  # poistion angle to rotate with
    PAdeg = args.PA
else:
    PAdeg = 0.0
PAstr = str(PAdeg)
#-------------------------------------------------------------------------------------------------------------
# writing output target list
RADECname = 'ra'+str(args.RADEC[0])+'dec'+str(args.RADEC[1])+'PA'+str(PAdeg).replace('.','p')
RADECname = RADECname.replace('.','p')
regfile   = 'MOSFIREds9_'+RADECname+'.reg'
outcat    = open(regfile,"w")

outcat.write('# Region file format: DS9 version 4.1 \n')
outcat.write('# Created with: '+sys.argv[0]+' \n')
outcat.write('# Regions show the MOSFIRE layout. \n')
outcat.write('# \n')
outcat.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1  \n')
outcat.write('fk5   \n')

# MOSFIRE FOV
RADECstr  = str(args.RADEC[0])+','+str(args.RADEC[1])
outcat.write('box('+RADECstr+',368.4",368.4",'+PAstr+') # color=magenta width=3 font="times 10 bold" text={MOSFIRE FOV}   \n')
outcat.write('circle('+RADECstr+',216") # color=magenta width=3   \n')
outcat.write('box('+RADECstr+',240",360",'+PAstr+') # color=magenta width=3   \n')

# GUIDE CAM
PArad = PAdeg*(2.0*np.pi)/360.0 # PA in radians
RAGC  = args.RADEC[0] + np.sin(PArad)*6.6/60.0
DECGC = args.RADEC[1] + np.cos(PArad)*6.6/60.0

# ---- KBS 121218 -----
#PArad = np.deg2rad(PAdeg) #*(2.0*np.pi)/360.0 # PA in radians
#d1    = np.deg2rad(90.0-args.RADEC[1])
#d2    = np.deg2rad(90.0-args.RADEC[1]+6.6)
#ra1   = np.deg2rad(args.RADEC[0])
#ra2   = ra1 - np.arccos( (np.cos(PArad)-np.sin(d1)*np.sin(d2)) / (np.cos(d1)*np.cos(d2)) )

#print np.rad2deg(np.arccos( (np.cos(PArad)-np.sin(d1)*np.sin(d2)) / (np.cos(d1)*np.cos(d2)) ))
#print np.sin(d1)*np.sin(d2) - np.cos(PArad), np.sin(d1)*np.sin(d2), np.cos(PArad), np.cos(d1)*np.cos(d2)
#print '---',ra1,d1,ra2,d2
#print args.RADEC,PAdeg
#RAGC  = np.rad2deg(ra2)
#DECGC = np.rad2deg(d2)
#print RAGC,DECGC
#=-----------------------------------

#RADECstrGC = str(RAGC)+','+str(DECGC)
#outcat.write('box('+RADECstrGC+',174",174",'+PAstr+') # color=magenta width=3 font="times 10 bold" text={Guide Cam}  \n')

# COMPASS
#outcat.write('# compass('+RADECstrGC+',2") compass=fk5 {N} {E} 1 1 color=magenta width=2 font="times 10 bold"   \n')

outcat.close()
if args.verbose: print 'Wrote mosfire region file for (RA,Dec) = ('+str(args.RADEC[0])+','+str(args.RADEC[1])+') to ./'+regfile


#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

