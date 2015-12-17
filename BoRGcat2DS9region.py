#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# BoRGcat2DS9region.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Turning a BoRG catalog (on the *_multiband.cat format for instance created with
# createBoRGmultibandcat.py) into a DS9 region file that can be overplotted on the fits images
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# catalog          : Path and name of BoRG catalog to turn into a DS9 region file
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --color          : to change color from (default = red) provide name of color to use
# --size           : to change size of circles (default = 2'') provide value in arcsec
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# regionfile       : Created region file. Will be place in the same directory as the catalog
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x BoRGcat2DS9region.py       (only required once)
# bash> BoRGcat2DS9region.py '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121108BPZ/121108223229_BPZrun_borg_0751+2917_multiband_modified_BPZinput/borg_0751+2917_multiband_modified_BPZinput_postBPZmodified.cat' --verbose
#
# bash> BoRGcat2DS9region.py '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121108BPZ/12110823279_BPZrun_borg_1033+5051_multiband_modified_BPZinput/borg_1033+5051_multiband_modified_BPZinput_postBPZmodified.cat' --verbose
#
# bash> BoRGcat2DS9region.py '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/modifiedcatalogs131301_brightobj/borg_1033+5051_multiband_modified.cat' --verbose
#
#
#
# --- STARS ---
# bash> BoRGcat2DS9region.py '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/modifiedcatalogs121121starsMAGcut/borg_0751+2917_multiband_modifiedSTARS.cat' --verbose --color 'yellow' --size 1.0
#
# bash> BoRGcat2DS9region.py '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/modifiedcatalogs121121starsMAGcut/borg_1033+5051_multiband_modifiedSTARS.cat' --verbose --color 'yellow' --size 1.0
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-11-21  started by K. B. Schmidt (UCSB)
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
parser.add_argument("catalog", type=str, help="Provide path and name of file containing list of catalog names")
# ---- optional arguments ----
parser.add_argument("--color", type=str, help="Color to use (default = 'red')")
parser.add_argument("--size", type=float, help="Size of circles (default = 2.0)")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# reading catalog
cat  = np.genfromtxt(args.catalog, dtype=None, comments='#')
Nobj = len(cat)

#-------------------------------------------------------------------------------------------------------------
# writing output target list
regfile  = args.catalog.replace('.cat','.reg')
outcat   = open(regfile,"w")

outcat.write('# Region file format: DS9 version 4.1 \n')
outcat.write('# Create with: '+sys.argv[0]+' \n')
outcat.write('# File marks the '+str(Nobj)+' objects from: \n')
outcat.write('# '+args.catalog+' \n')
outcat.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1   \n')
outcat.write('fk5  \n')

RAall   = cat['f4']      # object RA
DECall  = cat['f5']      # object DEC
IDall   = cat['f0']      # object ID
if args.color: 
    COLOR   = args.color            # radius of circle in arcsec
else:
    COLOR   = 'red'
if args.size: 
    SIZE    = args.size            # radius of circle in arcsec
else:
    SIZE    = 2.0

for ii in range(Nobj):
    RA   = RAall[ii]      # object RA
    DEC  = DECall[ii]     # object DEC
    ID   = IDall[ii]      # object ID
    outcat.write('circle('+str(RA)+','+str(DEC)+','+str(SIZE)+'") # color='+COLOR+' width=2 font="times 10 bold" text={ID: '+str(ID)+'}  \n')

outcat.close()
if args.verbose: print ':: '+sys.argv[0]+' :: turned input catalog into the region file: '+regfile

#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

