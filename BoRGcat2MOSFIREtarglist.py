#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# BoRGcat2MOSFIREtarglist.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Procedure turning a BoRG catalog into a template MOSFIRE target list which
# can be loaded into MAGMA and be used to create MOSFIRE masks.
#----------------------------
#   COMMENTS
#----------------------------
# A key component of the target list is the priority as the MOSFIRE slit mask
# optimisation aims at finding the mask with the highes total priority, i.e., 
# maximising Sigma_i(prirotiy_i) where i counts the objecst chosen for the mask.
#
# When modfying the values in the target list, e.g., by giving higher priority
# to a set of particular objects it is stronly suggested to rename the file (for
# instance by removing the RAW extension in the file name) as the script will
# overwrite the file in case a new targetlist is created from the same input catalog.
#
# If the --alignmentstars keyword is not used, no allignment stars are added to the
# target list and these have to be added manually before the list is fed to MAGMA.
#----------------------------
#   INPUTS:
#----------------------------
# borgcat          : Path and name of BoRG catalog
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --astrooffset    : The astrometric offset (expects RAoffset DECoffset in arcsec.) 
#                    relative to absolute astrometry to apply: 
#                       RA_absolute  = RA_borgcat  + RAoffset[0]
#                       Dec_absolute = Dec_borgcat + DECoffset[1]
#                    This could for instance be obtained from matching stars in the field 
#                    to the 2MASS catlogs at http://irsa.ipac.caltech.edu
# --alignmentstars : Path and name to a catalog (on the BoRG format) containing the
#                    alignment stars MAGMA can use when optimizing.
#                    If this keyword is not set no alignment stars (priority = -1) are 
#                    put in the outputted target list and these stars then have to be 
#                    manually before the target list is fed to MAGMA.
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# targetlist       : Output file containing the MOSFIRE target list with the following columns:
#                        1) Target Name
#                        2) Priority (-1 for alignment stars)
#                        3) Magnitude/Brightness: F125W (J) automag
#                        4) Right Ascension
#                        5) Declination
#                        6) Epoch
#                        7) Equinox
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x BoRGcat2MOSFIREtarglist.py       (only required once)
# bash> BoRGcat2MOSFIREtarglist.py '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121025/borg_1437+5043_multiband_modified.cat' --verbose --alignmentstars '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121026/borg_1437+5043_multiband_modifiedSTARS.cat'
#
# BoRGcat2MOSFIREtarglist.py '/Users/kasperborelloschmidt/work/BoRG/borg_cats/borg_1437+5043_multiband.cat' --verbose 
#
# OBSERVATION DEC2012
#
# bash> BoRGcat2MOSFIREtarglist.py '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121108BPZ/121108223229_BPZrun_borg_0751+2917_multiband_modified_BPZinput/borg_0751+2917_multiband_modified_BPZinput_postBPZmodified.cat' --verbose --alignmentstars '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/modifiedcatalogs121214stel0p6magLT20/borg_0751+2917_multiband_modifiedSTARS.cat' --astrooffset 0.16479492 -0.16136169

# bash> BoRGcat2MOSFIREtarglist.py '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/modifiedcatalogs131301_brightobj/borg_1033+5051_multiband_modified.cat' --verbose

#
# bash> BoRGcat2MOSFIREtarglist.py '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121108BPZ/12110823279_BPZrun_borg_1033+5051_multiband_modified_BPZinput/borg_1033+5051_multiband_modified_BPZinput_postBPZmodified.cat' --verbose --alignmentstars '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/modifiedcatalogs121214stel0p6magLT20/borg_1033+5051_multiband_modifiedSTARS.cat' --astrooffset 1.0711670 -0.21286011
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-11-25  started by K. B. Schmidt (UCSB)
# 2012-12-07  chagned default priority to depend on brightness and size. K. B. Schmidt (UCSB)
# 2012-12-10  added astrooffset keyword. K.B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import commands   # enable easy capture of commandline output
import numpy as np
import pdb
#----------------------------
#   FUNCTIONS
#----------------------------
# Uses 'skycoor' to convert R.A. and Dec. format
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("borgcat", type=str, help="Path and name of borg catalog file")
# ---- optional arguments ----
parser.add_argument("--astrooffset", type=float, nargs=2, help="Astrometric offset to apply: RAoffset Decoffset [arcsec]")
parser.add_argument("--alignmentstars", type=str, help="Path and name of file containing alignment stars (priority -1)")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# reading input catalog
borgcat  = np.genfromtxt(args.borgcat, dtype=None, comments='#')
Nobj     = len(borgcat['f1'])

#-------------------------------------------------------------------------------------------------------------
# Creating default priority (indicated with '-')
Pvec = [40] * Nobj   # default priority vector of length Nobj
Aminorax       = borgcat['f7']
Bmajorax       = borgcat['f8']
Areaellipse    = np.pi*Aminorax*Bmajorax
F125Wautoflux  = borgcat['f40']

Pvec = Pvec + F125Wautoflux/2.0 + 600/Areaellipse
Pvec = np.ceil(Pvec)*(-1)  # putting minus on priority to indicate default value

#-------------------------------------------------------------------------------------------------------------
# Applying atrometric offset to RA and Dec
if args.astrooffset:
    aoff = [args.astrooffset[0]/3600.0,args.astrooffset[1]/3600.0]  # turning input into degrees
else:
    aoff = [0.0,0.0]

RAvec  = borgcat['f4'] + np.array([aoff[0]]*Nobj)
Decvec = borgcat['f5'] + np.array([aoff[1]]*Nobj)
#-------------------------------------------------------------------------------------------------------------
# writing output target list
targlist = args.borgcat.replace('.cat','_MOSFIREtargetlistRAW.txt')
outcat   = open(targlist,"w")
outcat.write('# MOSFIRE target list to be used with MAGMA  \n')
outcat.write('# Created with '+sys.argv[0]+' \n')
outcat.write('# from '+args.borgcat+' \n')
outcat.write('# File contains '+str(Nobj)+' objects \n')
outcat.write('#  \n')
if args.alignmentstars: 
    outcat.write('# Appending '+args.alignmentstars+' alignment stars catalog \n')
else:
    outcat.write('# Appending no alignment stars catalog \n')
outcat.write('#  \n')
outcat.write('# The columns are: \n')
outcat.write('#    1 Target Name  \n')
outcat.write('#    2 Priority (-1 for alignment stars) \n')
outcat.write('#    3 Magnitude/Brightness: F125W (J) automag \n')
outcat.write('#    4 Right Ascension \n')
outcat.write('#    5 Declination  \n')
outcat.write('#    6 Epoch \n')
outcat.write('#    7 Equinox \n')
outcat.write('#  \n')

objID    = borgcat['f0']     # used id name is not a string (happens if the file has been modified for BPZ)
name     = borgcat['f1']
priority = Pvec              # default priority vector
mag      = borgcat['f45']
ra       = RAvec  
dec      = Decvec 
epoch    = [2000.0] * Nobj   # creating vector of length Nobj 
equinox  = [2000.0] * Nobj   # creating vector of length Nobj 

for ii in range(Nobj):       # looping over objects and wirting to target list
    skycorrOUT  = commands.getoutput('skycoor '+str(ra[ii])+' '+str(dec[ii]))
    radec       = skycorrOUT.split(' ')
    if isinstance(name[ii],str): # checking that the name is a string
        objstring   = name[ii]+' '+str(priority[ii])+' '+str(mag[ii])+' '+radec[0]+' '+radec[1]+' '+str(epoch[ii])+' '+str(equinox[ii])+' \n'
    else: 
        objstring   = str(objID[ii])+' '+str(priority[ii])+' '+str(mag[ii])+' '+radec[0]+' '+radec[1]+' '+str(epoch[ii])+' '+str(equinox[ii])+' \n'
    outcat.write(objstring)

#-------------------------------------------------------------------------------------------------------------
# adding alignment stars to catalog if provided
if args.alignmentstars:
    starcat  = np.genfromtxt(args.alignmentstars, dtype=None, comments='#')
    Nstars   = starcat['f1'].size

    # Applying atrometric offset to RA and Dec
    RAvecstar  = starcat['f4'] + np.array([aoff[0]]*Nstars)
    Decvecstar = starcat['f5'] + np.array([aoff[1]]*Nstars)

    name     = starcat['f1']
    priority = [-1] * Nstars       # default priority vector of length Nstars
    mag      = starcat['f45']
    ra       = RAvecstar
    dec      = Decvecstar
    epoch    = [2000.0] * Nstars   # creating vector of length Nstars 
    equinox  = [2000.0] * Nstars   # creating vector of length Nstars 

    outcat.write('# \n')
    outcat.write('# Alignment stars \n')

    if Nstars > 1:
        for jj in range(Nstars):       # looping over objects and wirting to target list
            skycorrOUT  = commands.getoutput('skycoor '+str(ra[jj])+' '+str(dec[jj]))
            radec       = skycorrOUT.split(' ')
            objstring   = name[jj]+' '+str(priority[jj])+' '+str(mag[jj])+' '+radec[0]+' '+radec[1]+' '+str(epoch[jj])+' '+str(equinox[jj])+' \n'
            outcat.write(objstring)
    else:
        objstring = str(name)+' '+str(priority[0])+' '+str(mag)+' '+radec[0]+' '+radec[1]+' '+str(epoch[0])+' '+str(equinox[0])+' \n'
        outcat.write(objstring)

    if args.verbose: print ':: '+sys.argv[0]+' :: Adding the '+str(Nstars)+' alignment stars to the target list'

#-------------------------------------------------------------------------------------------------------------
outcat.close()
if args.verbose: print ':: '+sys.argv[0]+' :: wrote target list to '+targlist
#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

