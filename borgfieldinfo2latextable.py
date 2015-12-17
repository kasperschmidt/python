#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# borgfieldinfo2latextable.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Turning borg field info file into latex table rows
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# infofile         : infofile to turn into latex table rows
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
# bash> borgfieldinfo2latextable.py /Users/kasperborelloschmidt/work/BoRG/borgdata/130724toPascal/fieldinfo.txt --verbose 
# bash> borgfieldinfo2latextable.py /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_field_info_ALL.txt  --verbose 
# bash> borgfieldinfo2latextable.py /Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p0/borg_field_info.txt  --verbose 
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-07-24  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import commands            # get output from spawned command line processes
import getopt              # used to extract/obtain the optional input
import numpy as np         # enable opening with genfromtxt
import pdb                 # for debugging with pdb.set_trace()
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("infofile", type=str, help="BoRG field info file")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
# Reading ascii
dat    = np.genfromtxt(args.infofile, dtype=None, comments='#')
field  = dat['f0']
filter = dat['f1']
texp   = dat['f2']
area   = dat['f3']
mlim   = dat['f4']
ra     = dat['f5']
dec    = dat['f6']
gl     = dat['f7']
gb     = dat['f8']
ebv    = dat['f9']
nomzp  = dat['f10']
av     = dat['f11']
corzp  = dat['f12']

ufield, uindex = np.unique(field,return_index=True)
ura    = ra[uindex]
udec   = dec[uindex]
uebv   = ebv[uindex]
uarea  = area[uindex]
Nfield = len(ufield)
#-------------------------------------------------------------------------------------------------------------
#columns of table:
#field ra dec texp606 mlim606 texp098 mlim098 texp105 mlim105 texp125 mlim125 texp160 mlim160 area125 ebv 

for ii in range(Nfield):
    row = ' '+ufield[ii]+' & '+str("%.4f" % ura[ii])+' & '+str("%.4f" % udec[ii])

    ent606 = np.where((field == ufield[ii]) & (filter == 'F606W'))[0]
    if len(ent606) == 1:
        row = row+' & '+str("%.0f" % texp[ent606[0]])+' & '+str("%.2f" % mlim[ent606[0]])
    elif len(ent606) == 0:
        row = row+' & $\dots$ & $\dots$'
        
    ent098 = np.where((field == ufield[ii]) & (filter == 'F098M'))[0]
    row = row+' & '+str("%.0f" % texp[ent098[0]])+' & '+str("%.2f" % mlim[ent098[0]])

    ent105 = np.where((field == ufield[ii]) & (filter == 'F105W'))[0]
    if len(ent105) == 1:
        row = row+' & '+str("%.0f" % texp[ent105[0]])+' & '+str("%.2f" % mlim[ent105[0]])
    elif len(ent105) == 0:
        row = row+' & $\dots$ & $\dots$'

    ent125 = np.where((field == ufield[ii]) & (filter == 'F125W'))[0]
    row = row+' & '+str("%.0f" % texp[ent125[0]])+' & '+str("%.2f" % mlim[ent125[0]])

    ent160 = np.where((field == ufield[ii]) & (filter == 'F160W'))[0]
    row = row+' & '+str("%.0f" % texp[ent160[0]])+' & '+str("%.2f" % mlim[ent160[0]])

    row = row+' & '+str("%.2f" % uarea[ii])+' & '+str("%.3f" % uebv[ii])+' \\\\'

    if args.verbose: print row

if args.verbose: print '\n - \n - \n' 
#-------------------------------------------------------------------------------------------------------------
# turn high z candidates into latex table rows
candfile = '/Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_higzcandidates_ALL.txt'
cand     = np.genfromtxt(candfile, dtype=None, comments='#')
Nlines   = len(cand['f0'])
catID     = cand['f0']    
SN_J      = cand['f1']
Jmag      = cand['f2']
Jmagerr   = cand['f3']
YJcol     = cand['f4']
YJcolerr  = cand['f5']
JHcol     = cand['f6']
JHcolerr  = cand['f7']
xpix      = cand['f8']
ypix      = cand['f9']
catflag   = cand['f10']
stellarity= cand['f11']
SN_V      = cand['f12']
SN_Y      = cand['f13']
SN_H      = cand['f14']
value1    = cand['f15']


for jj in range(Nlines):
    row = ''
    row = row+str("%.0f" % catID[jj])
    row = row+' & '+str("%.2f" % Jmag[jj])
    row = row+' $\pm$ '+str("%.2f" % Jmagerr[jj])
    row = row+' & '+str("%.1f" % YJcol[jj])
    row = row+' $\pm$ '+str("%.1f" % YJcolerr[jj])
    row = row+' & '+str("%.1f" % JHcol[jj])
    row = row+' $\pm$ '+str("%.1f" % JHcolerr[jj])
    #row = row+' & '+str("%.1f" % xpix[jj])
    #row = row+' & '+str("%.1f" % ypix[jj])
    #row = row+' & '+str("%.1f" % catflag[jj])
    #row = row+' & '+str("%.1f" % stellarity[jj])
    row = row+' & '+str("%.1f" % SN_V[jj])
    row = row+' & '+str("%.1f" % SN_Y[jj])
    row = row+' & '+str("%.1f" % SN_J[jj])
    row = row+' & '+str("%.1f" % SN_H[jj])+' \\\\'
    #row = row+' & '+str("%.0f" % value1[jj])
    if args.verbose: print row
    
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
