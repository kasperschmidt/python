#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# plotMOSFIREspec.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Script making simple plots of mosfire 1Dspecs and ivar files.
# (The eps files can be plotted with plotMOSFIREepsfiles.py)
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
#
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --search_1Dspec  : string used with glob to get 1Dspec.fits files to plot
# --search_ivar    : string used with glob to get *ivar.fits  files to plot
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
# bash> plotMOSFIREspec.py --search_1Dspec '*1Dspec.fits' --search_ivar '*ivar.fits' --verbose
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
#parser.add_argument("searchstring", type=str, help="String to get filenames with glob")
# ---- optional arguments ----
parser.add_argument("--search_1Dspec", type=str, help="String to get 1Dspec filenames with glob")
parser.add_argument("--search_ivar", type=str, help="String to get ivar filenames with glob")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--eps", action="store_true", help="Turn plots into eps files")
parser.add_argument("--png", action="store_true", help="Turn plots into png files")
parser.add_argument("--show", action="store_true", help="Showing plots on screen")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------
#                                             1Dspec
#-------------------------------------------------------------------------------------------------------------
if args.search_1Dspec:
    files_1D  = glob.glob(args.search_1Dspec)
    Nfiles_1D = len(files_1D)
    if len(args.search_1Dspec.split('/')) > 1:
        path_1D   = '/'.join(args.search_1Dspec.split('/')[:-1])+'/'
    else:
        path_1D   = './'

    if args.verbose: print ' - Plotting S/N histogram of 1Dspec to '+path_1D+'*_plotMOSFIREspec_SNhist1Dspec.pdf'
    for ff in xrange(Nfiles_1D):
        Fsize = 13
        plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
        plt.rc('font', family='serif',size=Fsize)           # setting text font
        plt.rc('xtick', labelsize=Fsize) 
        plt.rc('ytick', labelsize=Fsize) 
        name  = files_1D[ff].split('/')[-1]
        plotname = path_1D+name.replace('.fits','_plotMOSFIREspec_SNhist1Dspec.pdf')
        if args.eps: plotname = plotname.replace('.pdf','.eps')
        if args.png: plotname = plotname.replace('.pdf','.png')
        plt.clf()

        dat   = pyfits.open(files_1D[ff]) # reading data
        datTB = dat[1].data  # assuming first extension is a table and loading it
        SN = datTB['SPEC1D_SUM']/datTB['NOISE']
        SNstdev = np.std(SN[~np.isnan(SN)])
        SNmean  = np.mean(SN[~np.isnan(SN)])
        SNmedian  = np.median(SN[~np.isnan(SN)])
    
        hist = plt.hist(SN[~np.isnan(SN)],color="k",bins=100,histtype="step",lw=2,label='stddev(SN) = '+str("%.2f" % SNstdev)+' ; mean(SN) = '+str("%.2f" % SNmean)+' ; median(SN) = '+str("%.2f" % SNmedian)) # Hist of non-nan values
    
        print name,'stddev(SN) =',SNstdev,'  mean(SN) =',SNmean
    
        leg = plt.legend(fancybox=True, loc='upper left',numpoints=1,prop={'size':12})
        leg.get_frame().set_alpha(0.6)

        plt.xlabel(r'S/N')
        plt.ylabel(r'Count')

        #plt.xlim(-5,5)
        #plt.ylim(0,10)

        plt.savefig(plotname)
        if args.show: plot.show()  # draw plot on screen   
#-------------------------------------------------------------------------------------------------------------
#                                              ivar
#-------------------------------------------------------------------------------------------------------------
if args.search_ivar:
    files_ivar  = glob.glob(args.search_ivar)
    Nfiles_ivar = len(files_ivar)
    if len(args.search_ivar.split('/')) > 1:
        path_ivar   = '/'.join(args.search_ivar.split('/')[:-1])+'/'
    else:
        path_ivar   = './'

    if args.verbose: print ' - Plotting histogram of ivar pix to '+path_ivar+'*_plotMOSFIREspec_histivar.pdf'
    for ff in xrange(Nfiles_ivar):
        Fsize = 13
        plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
        plt.rc('font', family='serif',size=Fsize)           # setting text font
        plt.rc('xtick', labelsize=Fsize) 
        plt.rc('ytick', labelsize=Fsize) 
        name  = files_ivar[ff].split('/')[-1]
        plotname = path_ivar+name.replace('.fits','_plotMOSFIREspec_histivar.pdf')
        if args.eps: plotname = plotname.replace('.pdf','.eps')
        if args.png: plotname = plotname.replace('.pdf','.png')
        plt.clf()

        ivar     = pyfits.getdata(files_ivar[ff],ext=0)  # reading data into array   [rows,columns]
        sizeivar = ivar.shape
        ivarval  = ivar[round(sizeivar[0]/2.),:] # extracting values at central row
        
        ivarstdev = np.std(ivarval[~np.isnan(ivarval)])
        ivarmean  = np.mean(ivarval[~np.isnan(ivarval)])
        ivarmedeian  = np.median(ivarval[~np.isnan(ivarval)])
    
        hist = plt.hist(ivarval[~np.isnan(ivarval)],color="k",bins=100,histtype="step",lw=2,label='stddev(ivar\_center) = '+str("%.2f" % ivarstdev)+' ; mean(ivar\_center) = '+str("%.2f" % ivarmean)+' ; '+str("%.2f" % ivarmedeian)) # Hist of non-nan values
    
        print name,'stddev(ivar) =',ivarstdev,'  mean(ivar) =',ivarmean
    
        leg = plt.legend(fancybox=True, loc='upper left',numpoints=1,prop={'size':12})
        leg.get_frame().set_alpha(0.6)

        plt.xlabel(r'ivar value')
        plt.ylabel(r'Count')

        #plt.xlim(-5,5)
        #plt.ylim(0,10)

        plt.savefig(plotname)
        if args.show: plot.show()  # draw plot on screen   
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
#                                                  END
#-------------------------------------------------------------------------------------------------------------
