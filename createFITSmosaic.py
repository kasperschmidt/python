#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# createFITSmosaic.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Script creating a mosaic of fits image postage stamps.
#----------------------------
#   COMMENTS
#----------------------------
# The fits images (postage stamps) can be created with createBoRGcutouts.py
# Creating multiple mosaics can be done with a script like
# /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_higzcandidates_ALL_pstamps130321_3arcsec/createfitsmosaics_BoRG_higzcandidates.sh
#
#Also take a look at the scripts in borgdata/BoRG13_postagestamps_130818/
#
# Note that the stretch is fixed (need to turn Vmin and Vmax into optional keywords)
#----------------------------
#   INPUTS:
#----------------------------
# fitslist         : comma separated list of fits images to turn into mosaic 
# namebase         : namebase/id to use in output mosaics/plots
# mosaicshape      : shape of mosaic (not enabled - default just creating 1xNimg mosaic)
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --bandnames      : comma seprarated list of names of band plotted (written on images)
# --circle         : to draw circle in each frame (around object) provide
#                       ra[deg] dec[deg] radius[deg]
# --colormap       : provide the colormap to use (DEFAULT is 'gray')
# --stretchval     : provide the stretch to use (DEFAULT is 'power'). Cuts are hardcoded
# --verbose        : set -verbose to get info/messages printed to the screen
# --eps            : saving created plots as eps files
# --png            : saving created plots as png files
# --show           : showning plots on screen
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
#
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x createFITSmosaic.py       (only required once)
# bash> createFITSmosaic.py borg_1510+1115_1705_c19_cutoutFROM_borg_1510+1115_f606w_wfc3uvis_drz.fits,borg_1510+1115_1705_c19_cutoutFROM_borg_1510+1115_f098m_wfc3ir_drz.fits,borg_1510+1115_1705_c19_cutoutFROM_borg_1510+1115_f125w_wfc3ir_drz.fits,borg_1510+1115_1705_c19_cutoutFROM_borg_1510+1115_f160w_wfc3ir_drz.fits BoRG_1510+1115_1705 --bandnames F606W,F098M,F125W,F160W --verbose --circle 227.54000 11.25111 0.000138888
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-03-20  started by K. B. Schmidt (UCSB)
# 2013-08-18  added colormap keyword. K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs #
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import numpy as np         # enable opening with genfromtxt
import matplotlib.pyplot as plt   # importing plotting packages
import aplpy               # for plotting fits images
import pdb                 # for debugging with pdb.set_trace()
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("fitslist", type=str, help="Comma seprated list of fits files to turn into mosaic")
parser.add_argument("namebase", type=str, help="Prefix to use in output names")
# ---- optional arguments ----
parser.add_argument("--bandnames", type=str, help="Comma seprated list of names of bands given in fits images")
parser.add_argument("--circle", type=float, nargs=3, help="Mark object with circle. Provide RA[deg] Dec[deg] Radius[deg]")
parser.add_argument("--colormap", type=str, help="Color map to use (default is 'gray')")
parser.add_argument("--stretchval", type=str, help="Color map stretch to use (default is 'power')")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--eps", action="store_true", help="Turn plots into eps files")
parser.add_argument("--png", action="store_true", help="Turn plots into png files")
parser.add_argument("--show", action="store_true", help="Showing plots on screen")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
imglist = args.fitslist.split(',')
Nimg    = len(imglist)

Nrow  = 1
Ncol  = Nimg

if args.bandnames: band = args.bandnames.split(',')

#-------------------------------------------------------------------------------------------------------------
# PLOTTING
plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=16)           # setting text font
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
imgsize  = 4.

plotname = './'+args.namebase+'_mosaic'
if args.verbose: print ':: '+sys.argv[0]+' :: Creating figure '+plotname
fig = plt.figure(figsize=(imgsize*Nimg,imgsize))  # create a figure object
fig.clf()                                        # clearing figure

for ii in range(Nimg):
    Nplot = ii

    xmin = 1./Nimg*ii+0.01
    ymin = 0.05
    dx   = 1./Nimg-2*0.01
    dy   = 0.90

    #print xmin,ymin,dx,dy
    img = aplpy.FITSFigure(imglist[ii],figure=fig,subplot=[xmin,ymin,dx,dy])
    if args.colormap:
        clrmap = args.colormap # colormaps: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
        txtcol = 'k'
    else:
        clrmap = 'gray'
        txtcol = 'r'

    if args.stretchval == 'power':
        stretchval = args.stretchval
        vminval    = -2e-2
        vmaxval    = 2e-2
        vmidval    = None
    elif args.stretchval == 'log':
        stretchval = args.stretchval
        vminval    = -0.02
        vmaxval    = 0.04
        vmidval    = -10.
    else: #default
        stretchval = 'power'
        vminval    = -2e-2
        vmaxval    = 2e-2
        vmidval    = None
        
    img.show_colorscale(cmap=clrmap,stretch=stretchval,vmin=vminval,vmax=vmaxval,vmid=vmidval) 
    img.show_grid()
    img.set_grid_alpha(0.2)

    #img.set_tick_labels_font(size='x-small')
    #img.set_axis_labels_font(size='small')
    img.tick_labels.hide_x()
    img.axis_labels.hide_x()
    img.tick_labels.hide_y()
    img.axis_labels.hide_y()

    if args.circle: img.show_circles(args.circle[0],args.circle[1],args.circle[2],edgecolor=txtcol,lw=2)
    if ii == 0: img.add_label(0.5,0.90,args.namebase.replace('_',' '),color=txtcol,relative=True,horizontalalignment='center',size=22)
    if args.bandnames: img.add_label(0.5,0.10,band[ii],color=txtcol,relative=True,horizontalalignment='center',size=30)

fig.savefig(plotname+'.pdf')
if args.eps: fig.savefig(plotname+'.eps')
if args.png: fig.savefig(plotname+'.png')
if args.show: plot.show()  # draw plot on screen

#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

