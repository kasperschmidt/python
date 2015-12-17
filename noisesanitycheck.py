#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# noisesanitycheck.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# short script checking S/N distribution as described in section 18.5
# (notesLT130616) of notes.
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# searchstring     : string used with glob to get files to check
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
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
# bash> noisesanitycheck.py '/Users/kasperborelloschmidt/work/observing/130423_MOSFIRE_BoRG/data/130813reduction/*1Dspec_CALIBRATED.fits' --verbose

# bash> noisesanitycheck.py '/Users/kasperborelloschmidt/work/observing/130423_MOSFIRE_BoRG/data/130813reduction_woIVARSCALE/*1Dspec_CALIBRATED.fits' --verbose

# bash> noisesanitycheck.py '/Users/kasperborelloschmidt/work/observing/130423_MOSFIRE_BoRG/data/130820reduction/*1Dspec_CALIBRATED.fits' --verbose

#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-08-15  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs #
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import commands            # get output from spawned command line processes
import getopt              # used to extract/obtain the optional input
import numpy as np         # enable opening with genfromtxt
import matplotlib.pyplot as plt   # importing plotting packages
import pdb                 # for debugging with pdb.set_trace()
import pyfits
import glob
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("searchstring", type=str, help="String to get filenames with glob")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--eps", action="store_true", help="Turn plots into eps files")
parser.add_argument("--png", action="store_true", help="Turn plots into png files")
parser.add_argument("--show", action="store_true", help="Showing plots on screen")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
files = glob.glob(args.searchstring)
Nfiles = len(files)
path   = '/'.join(args.searchstring.split('/')[:-1])+'/'

stdevtot  = []
meantot   = []
mediantot = []

if args.verbose: print ' - Plotting histograms to '+path+'*_noisesanitycheck_histogram.pdf'
for ff in xrange(Nfiles):
    Fsize = 13
    plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
    plt.rc('font', family='serif',size=Fsize)           # setting text font
    plt.rc('xtick', labelsize=Fsize) 
    plt.rc('ytick', labelsize=Fsize) 
    name  = files[ff].split('/')[-1]
    plotname = path+name+'_noisesanitycheck_histogram.pdf'
    if args.eps: plotname = plotname.replace('.pdf','.eps')
    if args.png: plotname = plotname.replace('.pdf','.png')
    plt.clf()

    dat   = pyfits.open(files[ff]) # reading data
    datTB = dat[1].data  # assuming first extension is a table and loading it
    SN = datTB['SCI_SPEC']/datTB['SCI_NOISE']
    hist = plt.hist(SN[~np.isnan(SN)],color="k",bins=100,histtype="step",lw=2,label=name.replace('_','\_')) # Hist of non-nan values
    hist = plt.hist(np.random.normal(0.0, 1.0, len(SN)),color="r",bins=100,histtype="step",lw=1,label=r'normal w $\mu = 0$ and $\sigma = 1$') # Hist of non-nan values
    
    print name,'stddev(SN) =',np.std(SN[~np.isnan(SN)]),'  mean(SN) =',np.mean(SN[~np.isnan(SN)])
    stdevtot.append(np.std(SN[~np.isnan(SN)]))
    meantot.append(np.mean(SN[~np.isnan(SN)]))
    mediantot.append(np.median(SN[~np.isnan(SN)]))

    leg = plt.legend(fancybox=True, loc='upper left',numpoints=1,prop={'size':12})
    leg.get_frame().set_alpha(0.6)

    plt.xlabel(r'S/N')
    plt.ylabel(r'Count')

    #plt.xlim(-5,5)
    #plt.ylim(0,10)

    plt.savefig(plotname)
    if args.show: plot.show()  # draw plot on screen   

if args.verbose: 
    print '\n - in summary the average values are :'
    print '<stdev>  = ',np.mean(np.asarray(stdevtot))
    print '<mean>   = ',np.mean(np.asarray(meantot))
    print '<median> = ',np.mean(np.asarray(mediantot))
            
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
#                                                  END
#-------------------------------------------------------------------------------------------------------------
