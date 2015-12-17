#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# mosfirecand2borgcat.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Script taking one of Michele Trenti's catalogs of candidates for
# MOSFIRE follow-up (final*.mosfire) and extracting the corresponding BoRG data.
# This makes it much easier to create DS9 refion files (BoRGcat2DS9region.py) 
# and mosfire target lists (BoRGcat2MOSFIREtarglist.py) for these objects.
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# mosfirecand      : Path and name of final*.mosfire catalog to extract data for (final*.mosfire)
# borgcat          : Path and name of BoRG catalog to extract data from (*multiband.cat)
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# catalog          : Catalog containing BoRG data for objects in mosfirecand catalog.
#                    Will be put in same directory as mosfirecand
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x mosfirecand2borgcat.py       (only required once)
# bash> mosfirecand2borgcat.py '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/final_catalog_borg_1033+5051.mosfire' '/Users/kasperborelloschmidt/work/BoRG/borgdata/catalogs_borg1033+5051/borg_1033+5051_multiband.cat' --verbose 

# bash> mosfirecand2borgcat.py '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/final_catalog_borg_0751+2917.mosfire' '/Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p0/multiband/borg_0751+2917_multiband.cat' --verbose 
#
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-12-17  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import commands   # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
import numpy as np# enable opening with genfromtxt
import matplotlib as plt
#----------------------------
#   FUNCTIONS
#----------------------------
#
def arr2str(arr):                           # taking line of np.array or np.void and turning into string for easy passing to file
    strout = str(arr).replace(",",' ').replace("'",' ').replace('[',' ').replace(']',' ').replace('(',' ').replace(')',' ')+'  \n'
    return strout
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("mosfirecand", type=str, help="Path and name of file with list of objects to extract data for")
parser.add_argument("borgcat", type=str, help="Path and name of BoRG catalog to extract data from")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# Reading catalogs
mos      = np.genfromtxt(args.mosfirecand, dtype=None, comments='#')
Nmos     = len(mos)
mosID    = mos['f0']

borg     = np.genfromtxt(args.borgcat, dtype=None, comments='#')
Nborg    = len(borg)
borgID   = borg['f0']

outcat = args.mosfirecand.replace('.mosfire','_mosfire.cat')
output = open(outcat,"w")
output.write('# BoRG catalog data extracted with '+sys.argv[0]+' \n')
output.write('# Objects from '+args.mosfirecand+' \n')
output.write('# BoRG data taken from '+args.borgcat+' \n')
output.write('# \n')

for ii in range(Nmos):
    goodent    = np.where(borgID == mosID[ii])
    if len(goodent) == 1:
        #outstr = arr2str(mos[ii])
        #output.write(outstr)
        outstr = arr2str(borg[goodent])
        output.write(outstr)
    else:
        print ' '
        print ':: '+sys.argv[0]+' :: ERROR :: More than 1 or no match to object '+str(mosID[ii]+1)
        print ' '
print ':: '+sys.argv[0]+' :: Wrote output to: '
print outcat
output.close()

#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

# 1 ID
# 2 Z_B
# 3 Z_B_MIN
# 4 Z_B_MAX
# 5 T_B
# 6 ODDS    
# 7 Z_ML
# 8 T_ML
# 9 CHI-SQUARED
# 10 M_0
