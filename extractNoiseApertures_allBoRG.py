#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# extractNoiseApertures_allBoRG.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# running extractNoiseApertures_allBoRG.py for all BoRG fields and saving the results
# to dictionary
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# outputdic        : name of dictionary (npz file) to save results to
# outputdir        : directory to put dictionary and fits images from extractNoiseApertures in
# Raper            : Radius of apertures to extract given in arcsec (Braadley et al. 2012 uses 0.32 for 5 sigma limits)
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --Naper          : The number of apertures to extract (will be used with extractNoiseApertures.py)
#                    default is using extractNoiseApertures_ongrid.py.
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
# --- On Grid ---
# bash> extractNoiseApertures_allBoRG.py extractNoiseApertures_allBoRG_output130830_ongrid_0p32.npz emptyAperturesBoRG130830_ongrid_0p32  0.32 --verbose

# bash> extractNoiseApertures_allBoRG.py extractNoiseApertures_allBoRG_output130918_ongrid_0p32.npz emptyAperturesBoRG130918_ongrid_0p32 0.32 --verbose

# --- 2000 apertures ---
# bash> extractNoiseApertures_allBoRG.py extractNoiseApertures_allBoRG_output130726_2k_0p32.npz emptyAperturesBoRG130830_2k_0p32 2000 0.32 --verbose

# bash> extractNoiseApertures_allBoRG.py extractNoiseApertures_allBoRG_output130726_2k_0p40.npz emptyAperturesBoRG130830_2k_0p40 2000 0.40 --verbose

#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-07-19  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import getopt              # used to extract/obtain the optional input
import numpy as np         # enable opening with genfromtxt
import pdb                 # for debugging with pdb.set_trace()
import commands
import glob
from time import localtime, strftime
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("outputdic", type=str, help="Dictionary to save results to")
parser.add_argument("outputdir", type=str, help="Directory to contain output")
parser.add_argument("Raper", type=float, help="Radius of apertures in arcsec")
# ---- optional arguments ----
parser.add_argument("--Naper", type=int, help="Number of apertures to extract")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
# get borg field names and files to run on
datapath = '/Users/kasperborelloschmidt/work/BoRG/borgdata/'

# --- release v2.0 ---
v2img = glob.glob(datapath+'version_2p0/borg_*/*f125*drz.fits')
v2rms = glob.glob(datapath+'version_2p0/borg_*/*f125*rms.fits')
v2seg = glob.glob(datapath+'version_2p0/segm/borg_*F125*segm.fits')

ent = np.where(np.asarray([v2img[ii].split('/')[-1][0:14] for ii in range(len(v2img))]) == 'borg_1510+1115')[0]
del v2img[ent]
del v2rms[ent]
del v2seg[ent]

# --- cycle 19 pt 1 ---
c19p1img = glob.glob(datapath+'cycle19_new_nov14_cat/borg_*/*f125*drz.fits')
c19p1rms = glob.glob(datapath+'cycle19_new_nov14_cat/borg_*/*f125*rms.fits')
c19p1seg = glob.glob(datapath+'cycle19_new_nov14_cat/borg_*/final*125*segm.fits')

# --- cycle 19 pt 2 ---
c19p2img = glob.glob(datapath+'cycle19_all_130617/cat_130625/borg_*/*f125*drz.fits')
c19p2rms = glob.glob(datapath+'cycle19_all_130617/cat_130625/borg_*/*f125*rms.fits')
c19p2seg = glob.glob(datapath+'cycle19_all_130617/cat_130625/borg_*/final*125*segm.fits')

# --- COMBINING ---
totimg     = v2img + c19p1img + c19p2img
totrms     = v2rms + c19p1rms + c19p2rms
totseg     = v2seg + c19p1seg + c19p2seg
fieldnames = [totimg[ii].split('/')[-1][0:14] for ii in range(len(totimg))]
Nfield     = len(fieldnames)
#-------------------------------------------------------------------------------------------------------------
# create empty apertures fill output dictionary 
outdir  = args.outputdir
if outdir[-1] != '/': outdir = outdir+'/'
    
dict = {}
for ii in xrange(Nfield):
    if args.verbose: print ' - Getting the apertures (r='+str(args.Raper)+'") from --> ',fieldnames[ii],' <-- on '+strftime("%a, %d %b %Y %H:%M:%S", localtime())
    imgfile = totimg[ii]
    rmsfile = totrms[ii]
    segfile = totseg[ii]
    if args.Naper:
        cmdstr ='extractNoiseApertures.py '+imgfile+' '+rmsfile+' '+segfile+' '+str(args.Naper)+' '+str(args.Raper)+' --outputdir '+outdir+' --verbose'
    else:
        cmdstr ='extractNoiseApertures_ongrid.py '+imgfile+' '+rmsfile+' '+segfile+' '+str(args.Raper)+' --outputdir '+outdir+' --verbose'
    os.system(cmdstr)
    dict[fieldnames[ii]+'_hdr'] = commands.getoutput('less '+glob.glob(outdir+fieldnames[ii]+'*_aperturearray.txt')[0]+' | head -n 4')
    dict[fieldnames[ii]] = np.loadtxt(glob.glob(outdir+fieldnames[ii]+'*_aperturearray.txt')[0])
#-------------------------------------------------------------------------------------------------------------
# Save dictionary to binary
np.savez(outdir+args.outputdic,**dict) # save array as binary file
if args.verbose: print 'Saved output dictionary to ',args.outputdic
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------

