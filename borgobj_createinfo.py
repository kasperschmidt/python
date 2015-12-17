#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# addsources2fits.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Adds sources to a 'detection image' to be used with sextractor forcing it to perform
# photometry at the specified location.
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# sources              : list of sources to add. Expects the format:
#                            'x1,y1,type1,x2,y2,type2,...,xN,yN,typeN'
#                        where x and y are the coordinates and type can be either 
#                        'pixel' (not enabled as of 130909) or 'radec'
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --verbose            : set -verbose to get info/messages printed to the screen
# --stop               : stoppping program before end for de-bugging
# --help               : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
#
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> borgobj_extractinfo.py borg_1437+5043_0637 /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2/multiband/borg_1437+5043_multiband_region2.cat --dropoutcat /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2/borg_1437+5043/final_dropouts.cat  --verbose 
#
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-08-09  started by K. B. Schmidt (UCSB)
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
parser.add_argument("object", type=str, help="Object to extract info for")
parser.add_argument("multibandcat", type=str, help="multiband catalog containing object")
# ---- optional arguments ----
parser.add_argument("--dropoutcat", type=str, help="multiband catalog containing object")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
# turning name into ID
trimID = int(args.object.split('_')[-1])
if args.verbose: print ' - Extracted the ID '+str(trimID)+' from the object name ',args.object
#-------------------------------------------------------------------------------------------------------------
# Reading ascii
dat    = np.genfromtxt(args.multibandcat, dtype=None, comments='#')

#   1 id                    = dat['f0']
cid                   = dat['f1']
#   3 x                     = dat['f2']
#   4 y                     = dat['f3']
ra                    = dat['f4']
dec                   = dat['f5']
#   7 kron_radius           = dat['f6']
#   8 a_image               = dat['f7']
#   9 b_image               = dat['f8']
#  10 theta_image           = dat['f9']
#  11 fwhm_image            = dat['f10']
#  12 fwhm_world            = dat['f11']
#  13 stellarity            = dat['f12']
#  14 seflags               = dat['f0']
#  15 f606w_bkgrd           = dat['f0']
#  16 f606w_isoflux         = dat['f0']
#  17 f606w_isofluxerr      = dat['f0']
#  18 f606w_isosnr          = dat['f0']
#  19 f606w_autoflux        = dat['f0']
#  20 f606w_autofluxerr     = dat['f0']
#  21 f606w_autosnr         = dat['f0']
#  22 f606w_isomag          = dat['f0']
#  23 f606w_isomagerr       = dat['f0']
f606w_automag         = dat['f23']
f606w_automagerr      = dat['f24']
#  26 f098m_bkgrd           = dat['f0']
#  27 f098m_isoflux         = dat['f0']
#  28 f098m_isofluxerr      = dat['f0']
#  29 f098m_isosnr          = dat['f0']
#f098m_autoflux        = dat['f29']
#f098m_autofluxerr     = dat['f30']
#  32 f098m_autosnr         = dat['f0']
#  33 f098m_isomag          = dat['f0']
#  34 f098m_isomagerr       = dat['f0']
f098m_automag         = dat['f34']
f098m_automagerr      = dat['f35']
#  37 f125w_bkgrd           = dat['f0']
#  38 f125w_isoflux         = dat['f0']
#  39 f125w_isofluxerr      = dat['f0']
#  40 f125w_isosnr          = dat['f0']
#  41 f125w_autoflux        = dat['f0']
#  42 f125w_autofluxerr     = dat['f0']
#  43 f125w_autosnr         = dat['f0']
#  44 f125w_isomag          = dat['f0']
#  45 f125w_isomagerr       = dat['f0']
f125w_automag         = dat['f45']
f125w_automagerr      = dat['f46']
#  48 f160w_bkgrd           = dat['f0']
#  49 f160w_isoflux         = dat['f0']
#  50 f160w_isofluxerr      = dat['f0']
#  51 f160w_isosnr          = dat['f0']
#  52 f160w_autoflux        = dat['f0']
#  53 f160w_autofluxerr     = dat['f0']
#  54 f160w_autosnr         = dat['f0']
#  55 f160w_isomag          = dat['f0']
#  56 f160w_isomagerr       = dat['f0']
f160w_automag         = dat['f56']
f160w_automagerr      = dat['f57']
#  59 f105w_bkgrd           = dat['f0']
#  60 f105w_isoflux         = dat['f0']
#  61 f105w_isofluxerr      = dat['f0']
#  62 f105w_isosnr          = dat['f0']
#  63 f105w_autoflux        = dat['f0']
#  64 f105w_autofluxerr     = dat['f0']
#  65 f105w_autosnr         = dat['f0']
#  66 f105w_isomag          = dat['f0']
#  67 f105w_isomagerr       = dat['f0']
f105w_automag         = dat['f67']
f105w_automagerr      = dat['f68']

#-------------------------------------------------------------------------------------------------------------
if args.verbose: print ' - Info from multiband catalog: '
ent = np.where(cid == args.object)[0]

row = ''
row = row+cid[ent][0]
row = row+' & '+str("%.2f" % f606w_automag[ent])
row = row+' $\pm$ '+str("%.2f" % f606w_automagerr[ent])
row = row+' & '+str("%.2f" % f098m_automag[ent])
row = row+' $\pm$ '+str("%.2f" % f098m_automagerr[ent])
row = row+' & '+str("%.2f" % f105w_automag[ent])
row = row+' $\pm$ '+str("%.2f" % f105w_automagerr[ent])
row = row+' & '+str("%.2f" % f125w_automag[ent])
row = row+' $\pm$ '+str("%.2f" % f125w_automagerr[ent])
row = row+' & '+str("%.2f" % f160w_automag[ent])
row = row+' $\pm$ '+str("%.2f" % f160w_automagerr[ent])+' \\\\'

if args.verbose: print row
#-------------------------------------------------------------------------------------------------------------
# turn high z candidates into latex table rows
if args.dropoutcat:
    if args.verbose: print ' - Info from dropout catalog: '
    cand     = np.genfromtxt(args.dropoutcat, dtype=None, comments='#')
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

    ent = np.where(catID == trimID)[0]
    row = ''
    row = row+args.object
    row = row+' & '+str("%.2f" % Jmag[ent])
    row = row+' $\pm$ '+str("%.2f" % Jmagerr[ent])
    row = row+' & '+str("%.1f" % YJcol[ent])
    row = row+' $\pm$ '+str("%.1f" % YJcolerr[ent])
    row = row+' & '+str("%.1f" % JHcol[ent])
    row = row+' $\pm$ '+str("%.1f" % JHcolerr[ent])
    #row = row+' & '+str("%.1f" % xpix[ent])
    #row = row+' & '+str("%.1f" % ypix[ent])
    #row = row+' & '+str("%.1f" % catflag[ent])
    #row = row+' & '+str("%.1f" % stellarity[ent])
    row = row+' & '+str("%.1f" % SN_V[ent])
    row = row+' & '+str("%.1f" % SN_Y[ent])
    row = row+' & '+str("%.1f" % SN_J[ent])
    row = row+' & '+str("%.1f" % SN_H[ent])+' \\\\'
    #row = row+' & '+str("%.0f" % value1[ent])
    
    if args.verbose: print row
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
