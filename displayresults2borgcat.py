#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# displayresults2borgcat.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Script taking the output from a display_v2 inspection and extracting
# the correspondng info from the BoRG catalogs. 
# It also makes for an easy way to format the objects into MOSFIRE target 
# list lines by using the printMOSFIREtargets and copypasting the output into
# the targetlists.
#----------------------------
#   COMMENTS
#----------------------------
# 
#----------------------------
#   INPUTS:
#----------------------------
# displayresults   : file with the results from running display_v2
# borgcat          : Path and name of directory containing the BoRG catalogs 
#                    to extract data from (*multiband.cat)
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --printMOSFIREtargets : set this keywrod to get a line printed for each object
#                         in the MOSFIRE-targetlist format (with priority from 
#                         displayresilts file).
#                         Note: this can also be done with BoRGcat2MOSFIREtarglist.py
#                               on the files create with this script but then the
#                               results from display_v2 (the priorities) has to be added
#                               by hand.
# --verbose             : set -verbose to get info/messages printed to the screen
# --show                : showning plots on screen
# --stop                : stoppping program before end for de-bugging
# --help                : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# Directoryfile    : Output file containing the name of the directories the output has been
#                    put into. Directories will be moved to directory of catlist
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x displayresults2borgcat.py       (only required once)
# bash> displayresults2borgcat.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/all_final_dropouts_cats_DisplayInspection130122.txt /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/multiband/ --verbose --printMOSFIREtargets
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-04-13  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs #
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import commands            # get output from spawned command line processes
import numpy as np
import pdb                 # for debugging with pdb.set_trace()
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("displayresults", type=str, help="Path and name of file containing results from running display_v2")
parser.add_argument("borgcat", type=str, help="Directory containing the catalogs (*_multiband.cat) to extract info from")
# ---- optional arguments ----
parser.add_argument("--printMOSFIREtargets", action="store_true", help="Plot MOSFIRE targetlist line for each object")
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
# Reading catalogs
dispres     = np.genfromtxt(args.displayresults, dtype=None, comments='#')
Nobj        = len(dispres)
names       = dispres['f0']
priority    = dispres['f1']
fields      = []
for ii in range(Nobj): fields.append(names[ii][0:14])
ufields     = np.unique(fields)

cmdstr      = 'ls '+args.borgcat+'*_multiband.cat | cat'
multicats   = commands.getoutput(cmdstr).split('\n')
multifields = []
for ii in range(len(multicats)): multifields.append(multicats[ii].split('/')[-1].split('_multiband')[0])

print ':: '+sys.argv[0]+' :: Wrote output to: '
for jj in range(len(ufields)):
    ent     = np.where(np.asarray(multifields) == ufields[jj])
    if len(ent) == 1:
        ent = ent[0]
    else:
        sys.exit('ERROR: More than 1 match to '+ufields[jj]+"; there shouldn't be --> ABORTING")

    borgdat = np.genfromtxt(multicats[ent], dtype=None, comments='#')
    borgID  = borgdat['f1']
    outcat  = multicats[ent].replace('.cat','_objectextract.cat')
    output  = open(outcat,"w")
    output.write('# BoRG catalog data extracted with '+sys.argv[0]+' \n')
    output.write('# Objects from '+args.displayresults+' \n')
    output.write('# BoRG data taken from '+multicats[ent]+' \n')
    output.write('# \n')

    pvalues   = priority[np.where(np.asarray(fields) == ufields[jj])]
    objfield  = names[np.where(np.asarray(fields) == ufields[jj])]
    Nobjfield = len(objfield)
    for kk in range(Nobjfield):
        goodent    = np.where(borgID == objfield[kk])
        outstr = kbs.arr2str(borgdat[goodent])
        output.write(outstr)

        if args.printMOSFIREtargets:
            skycorrOUT  = commands.getoutput('skycoor '+str(borgdat[goodent]['f4'][0])+' '+str(borgdat[goodent]['f5'][0]))
            radec       = skycorrOUT.split(' ')
            mag         = borgdat[goodent]['f45'][0] # F125W (J) automag
            MOSstr   = objfield[kk]+' '+str(pvalues[kk])+' '+str(mag)+' '+radec[0]+' '+radec[1]+' 2000.0 2000.0'
            print MOSstr

    print '  --> '+outcat
    output.close()



#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

