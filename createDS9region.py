#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# createDS9region.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Read an ascii or fits catalog of objects and create a
# DS9 region file (marking objects with cirles)
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# catalog          : Path and name of catalog to turn into a DS9 region file
# IDRADECcol       : the column number (for ascii input) or column names (for fits input)
#                    for the ID, RA, and Dec columns to use for DS9 region
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --color          : to change color from (default = red) provide name of color to use
# --size           : to change size of circles (default = 2'') provide value in arcsec
# --regionmame     : give path and name of output region file to change from the default
#                    naming scheme: catalog.replace('extension','_ds9region+extension')
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# regionfile       : DS9 region file created.
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> createDS9region.py ~/work/GitHub/GLASS/GrismReduction/candidatesOfInterest/MACS0717_v2p0.txt 0 1 2 --color green --size 1 --verbose

# bash> createDS9region.py hlsp_clash_hst_acs-ir_macs0717_cat.txt 0 1 2 --color cyan --size 0.8 --verbose

# bash> createDS9region.py /Users/kasperborelloschmidt/work/GLASS/MACS0717test/hlsp_clash_hst_acs-ir_macs0717_cat_reformat_flux_AB25.fits id ra dec --color magenta --size 1 --verbose

# bash> createDS9region.py /Users/kasperborelloschmidt/work/GLASS/MACS0717test/hlsp_clash_hst_acs-ir_macs0717_cat_reformat_flux_AB25_noF160Wcut.fits  id ra dec --color red --size 1.2 --verbose

# bash> createDS9region.py hlsp_clash_hst_ir_macs0717_cat.txt 0 1 2 --color yellow --size 0.9 --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-12-16  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import commands   # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
import numpy as np# enable opening with genfromtxt
import pdb
import pyfits
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("catalog", type=str, help="Provide path and name of file containing list of catalog names")
parser.add_argument("IDRADECcol", type=str, nargs=3, help="column number (ascii input) or column names (fits input) for ID, RA, and Dec")
# ---- optional arguments ----
parser.add_argument("--color", type=str, help="Color to use (default = 'red')")
parser.add_argument("--size", type=float, help="Size of circles (default = 2.0)")
parser.add_argument("--regionname", type=str, help="Name of output region file in case default should not be used")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
# reading catalog
if args.catalog[-5:-1] == '.fit': # if input file is fits table
    data     = pyfits.open(args.catalog)[1].data
    Nobj     = len(data[args.IDRADECcol[0]])
    IDall    = data[args.IDRADECcol[0]]
    RAall    = data[args.IDRADECcol[1]]
    Decall   = data[args.IDRADECcol[2]]

else: # if not fits table, assume it's an ascii table
    data     = np.genfromtxt(args.catalog, dtype=None, comments='#')
    Nobj     = len(data['f'+str(int(args.IDRADECcol[0]))])
    IDall    = data['f'+str(int(args.IDRADECcol[0]))]
    RAall    = data['f'+str(int(args.IDRADECcol[1]))]
    Decall   = data['f'+str(int(args.IDRADECcol[2]))]
#-------------------------------------------------------------------------------------------------------------
# writing output target list
if args.regionname:
    regfile =  args.regionname
else:
    regfile  = args.catalog.split('.')[0]+'_ds9region.reg'
outcat   = open(regfile,"w")

outcat.write('fk5\n')

if args.color:
    COLOR   = args.color            # radius of circle in arcsec
else:
    COLOR   = 'red'

if args.size:
    SIZE    = args.size            # radius of circle in arcsec
else:
    SIZE    = 2.0

for ii in range(Nobj):
    ID   = IDall[ii]      # object ID
    RA   = RAall[ii]      # object RA
    Dec  = Decall[ii]     # object DEC

    outcat.write('circle('+str(RA)+','+str(Dec)+','+str(SIZE)+'") # color='+COLOR+' width=2 font="helvetica 14 normal roman" text={'+str(ID)+'}  \n')

outcat.close()
if args.verbose: print ' - Turned input catalog into the region file: \n   '+regfile
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------

