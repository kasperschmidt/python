#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# plotMOSFIREepsfiles.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# reading and plotting *eps*.fits files outputted from MOSFIRE reduction
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# fitslist         : list of *eps*.fits images to plot, e.g., 
#                    image_Y_eps.fits image_J_eps.fits image_H_eps.fits
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --namebase       : name-base of output plot(s). Default is using name of first fits file 
# --sumrows        : number of rows to sum spectrum over around peak in collapse spec (max +/- sumrows) 
# --verbose        : set -verbose to get info/messages printed to the screen
# --eps            : saving created plots as eps files
# --png            : saving created plots as png files
# --stop           : stoppping program at before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# plots            : Plots of the spectra. Saved in same directory as fits files
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x plotMOSFIREepsfiles.py       (only required once)
# bash> plotMOSFIREepsfiles.py borg_0751+2917_Y_borg_0751+2917_0175_eps.fits borg_0751+2917_2_J_borg_0751+2917_0175_eps.fits --namebase borg_0751+2917_0175 --sumrows 10 --verbose

# plotMOSFIREepsfiles.py LONGSLIT-25x0.7_J_eps.fits --namebase NGC1227 --sumrows 35 --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-01-04  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs #
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import commands            # get output from spawned command line processes
import getopt              # used to extract/obtain the optional input
import pyfits
import pywcs
import numpy as np         # enable opening with genfromtxt
import matplotlib.pyplot as plt 
import pdb                 # for debugging
from scipy import signal
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("fitslist", nargs='+', help="list of *eps*.fits images to plot")
# ---- optional arguments ----
parser.add_argument("--namebase", type=str, help="name-base for plot(s) to use instead of fits name")
parser.add_argument("--sumrows", type=int, help="number of rows to sum spectrum over around peak in cloppsed spec (max +/- sumrows)")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--eps", action="store_true", help="Turn plots into eps files")
parser.add_argument("--png", action="store_true", help="Turn plots into png files")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
Nimg = len(args.fitslist)
if args.verbose: print 'Number of spectra to be plotted: ',Nimg

if args.sumrows: # number of rows to sum over
    Nrowsum = args.sumrows
else:
    Nrowsum = 10 # default value

for ii in range(Nimg): # Creating array of data for fits files
#    print pyfits.info(args.fitslist[ii])
    dat    = pyfits.getdata(args.fitslist[ii],ext=0)  # reading data into array   [rows,columns]
    size   = dat.shape
    wave   = range(size[1])                      # vector to contain wavelengths

    hdulist = pyfits.open(args.fitslist[ii])     # Load the FITS hdulist using pyfits
    hdr     = hdulist[0].header                  # get header
    wcs = pywcs.WCS(hdr)  # Extract wcs (coordinate) information

    Xmin = 0.0
    Xmax = size[1]
    Ymin = 0.0
    Ymax = size[0]
    pixcrd = np.array([[Xmin,Ymin],[Xmax,Ymax]], np.float_)  # pixel coordinates to get 'sky'/cordinate values for
    sky = wcs.wcs_pix2sky(pixcrd, 1)  # convert to 'sky' values. 2nd arg. is "origin": 1 = a 1-based (Fortran-like) coordinates.

    lammin = sky[0,0]  # wavelength of first column 
    lammax = sky[1,0]  # wavelength of last  column 

    dlam   = (lammax - lammin)/size[1]
    wave   = np.multiply(lammin+np.multiply(wave,dlam),1e10)  # wavelengths in angstrom

    pixcrd2 = wcs.wcs_sky2pix(sky, 1)     # Convert the same coordinates back to pixel coordinates.
    assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6 # Sanity check; should be the same as original pixcrd, modulo floating-point error.

    rowsum = np.sum(dat,axis=1) # collapse columns of spectrum (slit profile)
    maxrow = np.where(rowsum == max(rowsum))  # number of row containing center of spectrum 
    spec1D = np.sum(dat[maxrow[0]-Nrowsum:maxrow[0]+Nrowsum,:],axis=0) # collapse specified rows (max slit profile +/- Nrowsum)

    for hh in range(size[1]):
        maxcol = np.max(dat[:,hh])             # max value in column hh
        maxent = np.where(dat[:,hh] == maxcol) # row number (entry) of maximum value
        if hh == 0: 
            maxcolent = maxent
        else:
            maxcolent = np.append(maxcolent,maxent)

    if ii == 0:  # appending to excisting arrays
        spec1Dtotal = spec1D
        wavetotal   = wave
        rowall      = size[0]
        colall      = size[1]
        rowtot      = rowsum
        speccenter  = maxrow
        colmaxall   = maxcolent
    else:
        spec1Dtotal = np.append(spec1Dtotal,spec1D)
        wavetotal   = np.append(wavetotal,wave) 
        rowall      = np.append(rowall,size[0])
        colall      = np.append(colall,size[1])
        rowtot      = np.append(rowtot,rowsum)
        speccenter  = np.append(speccenter,sum(rowall[0:-1])+maxrow)
        colmaxall   = np.append(colmaxall,maxcolent)


#-------------------------------------------------------------------------------------------------------------
# PLOTTING
plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=13)           # setting text font
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if args.namebase:
    base = args.namebase
else:
    base = args.fitslist[0].replace('.fits','')
plotname = base+'_spec1D'

if args.verbose: print 'Creating figure '+plotname
fig = plt.figure()  # create a figure object
fig.clf()                                        # clearing figure
ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
plt.title('MOSFIRE Spectrum: '+base.replace('_','\_'))

#sc = ax.scatter(wavetotal,spec1Dtotal,marker = 'o');
plt.plot(wavetotal,spec1Dtotal, 'k-',label='Spectrum: collapsed '+str(2*Nrowsum)+' rows')

ksize = 9 # size of kernel for median fileter (has to be odd)
medfiltspec = signal.medfilt(spec1Dtotal,kernel_size=ksize) # median filtering spectrum
plt.plot(wavetotal,medfiltspec, 'm-',label='Median spectrum. Kernel='+str(ksize)+' pix')

ax.grid(True,linestyle='-',color='0.75')
ax.set_xlabel('$\lambda$ [Angstrom]')
ax.set_ylabel('e/s')

# plt.vlines(12800, min(spec1Dtotal),max(spec1Dtotal), color='r', lw=4, linestyles='solid',label='Pa$\\beta$ at 12820 Angstrom') # Paschen beta line at 12800 angstrom

leg = plt.legend(fancybox=True, loc='lower right')  # add the legend in the middle of the plot
leg.get_frame().set_alpha(0.7)                      # set the alpha value of the legend: it will be translucent

plt.xlim([min(wavetotal),max(wavetotal)])
plt.ylim([min(spec1Dtotal),max(spec1Dtotal)])
fig.savefig(plotname+'.pdf')
if args.eps: fig.savefig(plotname+'.eps')
if args.png: fig.savefig(plotname+'.png')
# show()  # draw plot on screen

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if args.namebase:
    base = args.namebase
else:
    base = args.fitslist[0].replace('.fits','')
plotname = base+'_rows1D'

if args.verbose: print 'Creating figure '+plotname
fig = plt.figure()  # create a figure object
fig.clf()                                        # clearing figure
ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure

xval = range(np.sum(rowall))
#sc = ax.scatter(xval,rowtot,marker = 'o');
plt.plot(xval,rowtot, 'k-o')

ax.grid(True,linestyle='-',color='0.75')
ax.set_xlabel('rows (accumulative)')
ax.set_ylabel('e/s')

# drawing vertical lines to indicate spectrum seperation
plt.vlines(0, min(rowtot), max(rowtot), color='k', linestyles='solid')
if Nimg > 1:
    for jj in range(len(rowall)): # loop over row numbers of spectra
        plt.vlines(sum(rowall[0:jj]), min(rowtot), max(rowtot), color='k', linestyles='solid')
        plt.vlines(speccenter[jj]-Nrowsum,min(rowtot), max(rowtot), color='k', linestyles='--')
        plt.vlines(speccenter[jj]+Nrowsum,min(rowtot), max(rowtot), color='k', linestyles='--')
else:
    plt.vlines(speccenter[0]-Nrowsum,min(rowtot), max(rowtot), color='k', linestyles='--')
    plt.vlines(speccenter[0]+Nrowsum,min(rowtot), max(rowtot), color='k', linestyles='--')

plt.vlines(max(xval), min(rowtot), max(rowtot), color='k', linestyles='solid')

plt.xlim([min(xval)-5,max(xval)+5])
fig.savefig(plotname+'.pdf')
if args.eps: fig.savefig(plotname+'.eps')
if args.png: fig.savefig(plotname+'.png')
# show()  # draw plot on screen

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if args.namebase:
    base = args.namebase
else:
    base = args.fitslist[0].replace('.fits','')
plotname = base+'_cols1D'

if args.verbose: print 'Creating figure '+plotname
fig = plt.figure()  # create a figure object
fig.clf()                                        # clearing figure
ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure

xval = range(np.sum(colall))
plt.scatter(xval,colmaxall,marker='o',s=1)

ax.grid(True,linestyle='-',color='0.75')
ax.set_xlabel('columns (concattenated)')
ax.set_ylabel('e/s')

# drawing vertical lines to indicate spectrum seperation
plt.vlines(0, min(colmaxall), max(colmaxall), color='k', linestyles='solid')
if Nimg > 1:
    for jj in range(len(colall)): # loop over row numbers of spectra
        plt.vlines(sum(colall[0:jj]), min(colmaxall), max(colmaxall), color='k', linestyles='solid')

plt.vlines(max(xval), min(colmaxall), max(colmaxall), color='k', linestyles='solid')

plt.xlim([min(xval)-5,max(xval)+5])
plt.ylim([np.median(colmaxall)-Nrowsum*2,np.median(colmaxall)+Nrowsum*2])
fig.savefig(plotname+'.pdf')
if args.eps: fig.savefig(plotname+'.eps')
if args.png: fig.savefig(plotname+'.png')
# show()  # draw plot on screen

#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

