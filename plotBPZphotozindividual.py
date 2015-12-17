#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# plotBPZphotozindividual.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Given a BPZ output file individual figures of p(z) are plotted and 
# the data is potentially save to fits file (can be converted to ascii 
# with fits2ascii.py)
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# bpzfile          : Path AND name of BPZ output file
# objid            : The object(s) to plot and extract data for
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --savedat        : If set the p(z) data for the individual objects is saved to
#                    a fits table
# --ascii          : If savadat is set use this keyword to enabling saving of ascii 
#                    version of fits table
# --verbose        : set -verbose to get info/messages printed to the screen
# --eps            : saving created plots as eps files
# --png            : saving created plots as png files
# --show           : showning plots on screen for manipulation and saving
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# 
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> plotBPZphotozindividual.py bpzfile 560 345 --verbose

# bash> plotBPZphotozindividual.py /Users/kasperborelloschmidt/work/BoRG/BPZfullrundir/hdfnprior130122/130122193253_BPZrun_borg_0440-5244_multiband_BPZinput/borg_0440-5244_multiband_BPZinput.probs [647,650] --verbose

# MOSFIRE SPEC
# -- HDFN PRIOR ---
# bash> plotBPZphotozindividual.py /Users/kasperborelloschmidt/work/BoRG/BPZfullrun130313/hdfnprior130313/130313194859_BPZrun_borg_0951+3304_multiband_BPZinput/borg_0951+3304_multiband_BPZinput.probs [180,277] --verbose --savedat --ascii

# bash> plotBPZphotozindividual.py /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2/BPZ130419/hdfnprior130419/13041920549_BPZrun_borg_1437+5043_multiband_BPZinput/borg_1437+5043_multiband_BPZinput.probs [637]  --verbose --savedat --ascii

# bash> plotBPZphotozindividual.py /Users/kasperborelloschmidt/work/BoRG/BPZfullrun130313/hdfnprior130313/13031319516_BPZrun_borg_1510+1115_multiband_BPZinput/borg_1510+1115_multiband_BPZinput.probs [354,1218,1487,1524,1705]  --verbose --savedat --ascii

# -- FLAT PRIOR ---
# bash> plotBPZphotozindividual.py /Users/kasperborelloschmidt/work/BoRG/BPZfullrun130313/flatprior130313/130313211844_BPZrun_borg_0951+3304_multiband_BPZinput/borg_0951+3304_multiband_BPZinput.probs [180,277] --verbose --savedat --ascii

# bash> plotBPZphotozindividual.py /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2/BPZ130419/flatprior130419/130419205621_BPZrun_borg_1437+5043_multiband_BPZinput/borg_1437+5043_multiband_BPZinput.probs [637] --verbose --savedat --ascii
  
# bash> plotBPZphotozindividual.py /Users/kasperborelloschmidt/work/BoRG/BPZfullrun130313/flatprior130313/130313211936_BPZrun_borg_1510+1115_multiband_BPZinput/borg_1510+1115_multiband_BPZinput.probs [354,1218,1487,1524,1705]  --verbose --savedat --ascii
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-07-11  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import getopt              # used to extract/obtain the optional input
import numpy as np         # enable opening with genfromtxt
import matplotlib.pyplot as plt   # importing plotting packages
import pdb                 # for debugging with pdb.set_trace()
import pyfits
import commands            #
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("bpzfile", type=str, help="BPZ output to extract data from")
parser.add_argument("objid", type=str, help="list of IDs to plot (and extract data for): [ID1,ID2,...,IDN]")
# ---- optional arguments ----
parser.add_argument("--savedat", action="store_true", help="Set keyword to save data to fits table")
parser.add_argument("--ascii", action="store_true", help="Save ascii version of fits table as well")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--eps", action="store_true", help="Turn plots into eps files")
parser.add_argument("--png", action="store_true", help="Turn plots into png files")
parser.add_argument("--show", action="store_true", help="Showing plots on screen")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
probfile = args.bpzfile
pdat     = np.genfromtxt(probfile, dtype=None, comments='#')
objidstr = args.objid[1:-1].split(',')
Nobj     = len(objidstr)
objid    = [int(objidstr[jj]) for jj in xrange(Nobj)]

dz       = 0.0100
zmin     = 0.0100
zmax     = 10.0100
zval     = np.arange(zmin,zmax,dz)

for ii in xrange(Nobj):
    #----------------- PLOT NORMALIZED P(Z) -----------------
    plotname = args.bpzfile.split('multiband_BPZinput.pr')[0]+str(objid[ii])+'_photozprob.pdf'
    if args.eps: plotname.replace('.pdf','.eps')
    if args.png: plotname.replace('.pdf','.png')
    if args.verbose: print '\n - Plotting photo-z data to '+plotname
    objdat   = list(pdat[objid[ii]-1])
    #dat650 = 
    #dat647 = list(pdat[646])

    lwidth = 4
    plt.rcParams['axes.linewidth'] = lwidth/2.
    plt.rcParams['xtick.major.pad']='12'  # distance of tick labels from axis
    plt.rcParams['ytick.major.pad']='12'  # distance of tick labels from axis
    fig = plt.figure()
    fig.clf()
    plt.rcParams.update({'font.size': 40})
    ax = fig.add_subplot(1, 1, 1)
    ax.set_position([0.2,0.2,0.75,0.75])
    #ax.grid(True,linestyle='-',color='0.75')#,lw=lwidth
    ax.tick_params('both', length=7, width=lwidth/2., which='major')

    plt.plot(zval,np.divide(objdat[1:],dz),'k-',linewidth=lwidth)

    ax.set_xlabel('$z$')
    ax.set_ylabel('p$(z)$')

    fig.savefig(plotname)
    if args.show: plot.show()  # draw plot on screen

    # numerical integration
    #import scipy.integrate as si
    #si.simps(np.divide(dat650[1:],dz),np.arange(0.0100,10.0100,dz))
        
    if args.savedat:
        #----------------- SAVING DATA TO FITS TABLE -----------------
        fitsname = args.bpzfile.split('multiband_BPZinput.pr')[0]+str(objid[ii])+'_photozprob.fits'
        if args.verbose: print ' - Writing photo-z data to '+fitsname

        col1  = pyfits.Column(name='Z' , format='F', array=np.asarray(zval))
        col2  = pyfits.Column(name='PZ', format='F', array=np.asarray(objdat[1:]))
        cols  = pyfits.ColDefs([col1, col2])
        tbhdu = pyfits.new_table(cols)          # creating table header

        # writing hdrkeys:   '---KEY--',                  '----------------MAX LENGTH COMMENT-------------'
        tbhdu.header.append(('OBJID   ' ,args.objid[ii]   ,'Paramter file used with performSpecCal.py'),end=True)
        tbhdu.header.append(('DZ      ' ,dz               ,'Redshift intervals'),end=True)
        tbhdu.header.append(('ZMIN    ' ,zmin             ,'Minimum redshift'),end=True)
        tbhdu.header.append(('ZMAX    ' ,zmax             ,'Maximum redshift'),end=True)

        hdu      = pyfits.PrimaryHDU()             # creating primary (minimal) header
        thdulist = pyfits.HDUList([hdu, tbhdu])    # combine primary and table header to hdulist
        thdulist.writeto(fitsname,clobber=True)  # write fits file (clobber=True overwrites excisting file)  

        if args.ascii: # turning fits table into ascii if requested
            if args.verbose: print ' - Turning fits table into ascii'
            cmdstr = 'fits2ascii.py '+fitsname
            cmdout = commands.getoutput(cmdstr)
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
