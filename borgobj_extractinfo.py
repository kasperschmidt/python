#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# borgobj_extractinfo.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# This script makes it easy to extract the available info from a given BoRG object
# (and turns it into latex table rows).
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# object               : The name of the object to extract info for
# multibandcat         : File containing the multiband cant of object
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --dropoutcat         : If a dropout catalog excists this can be provided and info will 
#                        be extracted for the object from here as well.
# --sexcatextcol       : list of columns to extract in the original sextractor catalogs.
#                        Expects the format 'path,col1,col2,...,colN' where colX are the column 
#                        numbers to extract starting from 0 and path is the path to the directory
#                        containing the catalogs (end path with a /).
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

# -- getting aperture flux info in r=0.32'' (8 pixel diameter) apertures --
#   6 FLUX_APER       Flux vector within fixed circular aperture(s)   [count] 2,4,6,8,10 aperture diameter(s) in pixels
#  11 FLUXERR_APER    RMS error vector for aperture flux(es)          [count]

# bash> borgobj_extractinfo.py borg_1437+5043_0638 /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/multiband/borg_1437+5043_multiband.cat --sexcatextcol '/Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/borg_1437+5043/,8,13'   --verbose 

# bash> borgobj_extractinfo.py borg_1437+5043_0597 /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/multiband/borg_1437+5043_multiband.cat --sexcatextcol '/Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/borg_1437+5043/,8,13'   --verbose 

# bash> borgobj_extractinfo.py borg_1437+5043_0559 /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/multiband/borg_1437+5043_multiband.cat --sexcatextcol '/Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/borg_1437+5043/,8,13'   --verbose 

# bash> borgobj_extractinfo.py borg_1437+5043_0239 /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/multiband/borg_1437+5043_multiband.cat --sexcatextcol '/Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/borg_1437+5043/,8,13'   --verbose 

# bash> borgobj_extractinfo.py borg_1437+5043_0070 /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/multiband/borg_1437+5043_multiband.cat --sexcatextcol '/Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/borg_1437+5043/,8,13'   --verbose 

#----------------------------
#   BUGS
#----------------------------
# Only works for BoRG1437+5043 where there is F105W data. Need to add keyword to ignore F105W data.
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-09-09  started by K. B. Schmidt (UCSB)
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
parser.add_argument("--sexcatextcol", type=str, help="Sextractor catalog's path and columns to extract")
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
def getcatcolval(catalog,object,column):
    '''
    Function returning the value of a column for a given object in a provided catalog
    Note that columns start at 0.
    '''
    catdat     = np.genfromtxt(catalog, dtype=None, comments='#')    
    ent        = np.where(catdat['f0'] == object)[0]
    colval     = catdat['f'+str(int(column))][ent]
    return colval
#-------------------------------------------------------------------------------------------------------------
if args.sexcatextcol:
    if args.verbose: print ' - Info from original sextractor catalogs: '        
    path   = args.sexcatextcol.split(',')[0]
    cols   = args.sexcatextcol.split(',')[1:]
    Ncol   = int(len(cols))
    Nbands = 5
    str606 = 'final_sources_F606.cat'
    str098 = 'final_sources_F098.cat'
    str105 = 'final_sources_F105.cat'
    str125 = 'final_sources_F125.cat'
    str160 = 'final_sources_F160.cat'

    catfiles = [path+str606,path+str098,path+str105,path+str125,path+str160]

    row = ' '
    row = row+args.object
    for bb in xrange(Nbands):
        for cc in xrange(Ncol):
            val = getcatcolval(catfiles[bb],trimID,cols[cc])
            row = row+' & '+str("%.5f" % val)
            
    row = row+' \\\\'
    if args.verbose: print row
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
