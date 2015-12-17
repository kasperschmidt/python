#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# calculatePAofpair.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Calculating the position angle (degrees East of North) between two objects
#----------------------------
#   COMMENTS
#----------------------------
# 
#----------------------------
#   INPUTS:
#----------------------------
# ra1              : The RA in degrees of main object. Will try to split the string using space.
#                    if len(split) = 1 [RA] assumed to be fractional degrees 
#                    if len(split) = 3 [RA] assumed to be sexagesimal on the format hh mm ss.s
# dec1             : The Dec in degrees of main object. Will try to split the string using space.
#                    if len(split) = 1 [Dec] assumed to be fractional degrees 
#                    if len(split) = 3 [Dec] assumed to be sexagesimal on the format dd mm ss.s
# ra2              : The RA in degrees of second object. Will try to slit the string using space.
#                    if len(split) = 1 [RA] assumed to be fractional degrees 
#                    if len(split) = 3 [RA] assumed to be sexagesimal on the format hh mm ss.s
# dec2             : The Dec in degrees of second object. Will try to slit the string using space.
#                    if len(split) = 1 [Dec] assumed to be fractional degrees 
#                    if len(split) = 3 [Dec] assumed to be sexagesimal on the format dd mm ss.s
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --dist           : the distance between the two objects in arcsec. If given the offset in east-west
#                    and north-south will be estimated
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
# bash> calculatePAofpair.py '20 57 52.560' '+00 06 35.28' '20 57 50.681' '+00 06 06.65' --verbose

# bash> calculatePAofpair.py 17.55091 -2.38913 17.55494 -2.38513 --verbose # PA~45 degrees 
# bash> calculatePAofpair.py '01 10 12.218' '-02 23 20.87' '01 10 13.186' '-02 23 06.47' --verbose # PA~45 degrees 
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-09-17  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse            # argument managing
import sys                 # enabling arguments to code
import numpy as np         # enable opening with genfromtxt
import commands
import pdb                 # for debugging with pdb.set_trace()
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("ra1", type=str, help="Ra of priamry obj")
parser.add_argument("dec1", type=str, help="Dec of priamry obj")
parser.add_argument("ra2", type=str, help="Ra of secondary obj")
parser.add_argument("dec2", type=str, help="Dec of secondary obj")
# ---- optional arguments ----
parser.add_argument("--dist", type=float, help="Distance between objects in arcsec")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
#if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
rasplit  = args.ra1.split(' ')
decsplit = args.dec1.split(' ')
if (len(rasplit) == 3) & (len(decsplit) == 3):
    cmd     = 'skycoor -d '+args.ra1.replace(' ',':')+' '+args.dec1.replace(' ',':')
    radec1  = commands.getoutput(cmd).split()
    raval1  = float(radec1[0])/15.0
    decval1 = float(radec1[1])
elif (len(rasplit) == 1) & (len(decsplit) == 1):
    raval1  = float(args.ra1)/15.0
    decval1 = args.dec1
else:
    sys.exit('Invalid format of ra1/dec1 --> ABORTING')

rasplit  = args.ra2.split(' ')
decsplit = args.dec2.split(' ')
if (len(rasplit) == 3) & (len(decsplit) == 3):
    cmd     = 'skycoor -d '+args.ra2.replace(' ',':')+' '+args.dec2.replace(' ',':')
    radec2  = commands.getoutput(cmd).split()
    raval2  = float(radec2[0])/15.0
    decval2 = float(radec2[1])
elif (len(rasplit) == 1) & (len(decsplit) == 1):
    raval2  = float(args.ra2)/15.0
    decval2 = args.dec2
else:
    sys.exit('Invalid format of ra2/dec2 --> ABORTING')
    
#-------------------------------------------------------------------------------------------------------------
# Convert ra and dec to PA via IDL 'posang'
posangstr = 'posang,1, '+str(raval1)+','+str(decval1)+','+str(raval2)+','+str(decval2)+',PAang & print, PAang'
cmd       = "echo '"+posangstr+"' > getPAtemp.pro"
echoout   = commands.getoutput(cmd)
angout    = commands.getoutput('idl < getPAtemp.pro').split('\n')[-1]
rmout     = commands.getoutput('rm getPAtemp.pro')

angpos    = float(angout)
if angpos < 0.0: angpos = 360.0+angpos # making sure the angle is positive

if args.verbose: print ' - The position angle (east of north) was estimated to be ',angpos
#-------------------------------------------------------------------------------------------------------------
# calculate offset in east-west and north-south direction
if args.dist:
    NSoffset = args.dist*np.cos(angpos*np.pi/180.)
    EWoffset = args.dist*np.sin(angpos*np.pi/180.)
    if args.verbose: 
        print ' - Using distance ',args.dist,'arcsec the offset from primary (1st ra-dec set) to secondary (2nd ra-dec set) is:'
        print '   [W -> E, S -> N] =',[EWoffset,NSoffset]
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
#if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------

