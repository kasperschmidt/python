#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# plotBPZresults.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Plotting results from running runBPZmultiplecats.py
#----------------------------
#   COMMENTS
#----------------------------
# In the case of BoRG catalogs the files can be prepared using 
# borgcat2bpzcat.py
#----------------------------
#   INPUTS:
#----------------------------
# bpzcat           : Path AND name of file containing the output from the BPZ run to visualize 
#                    File ends on *_bpz.cat
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --savecounts     : save the histogram counts to text file for plotting
# --Nzbin          : --Nzbin int gives the number of redshifts bins to use when creating the
#                    animated gif from the pngs, hence keyword only effective with --png keyword
# --zrange         : Range of redshifts to bin [zmin,zmax]. Default is [0,10]"
# --verbose        : set -verbose to get info/messages printed to the screen
# --eps            : saving created plots as eps files
# --png            : saving created plots as png files
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# Directoryfile    : Output file containing the name of the directories the output has been
#                    put into. Directories will be moved to directory of catlist
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x plotBPZresults.py       (only required once)
# bash> plotBPZresults.py '/Users/kasperborelloschmidt/work/BoRG/BPZmultitest121015/121018xxxxxx_BPZrun_borg_0110-0224_multiband_BPZinputTESTfull/borg_0110-0224_multiband_BPZinputTESTfull_bpz.cat' --verbose --png --Nzbin 2 --savecounts
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-10-18 started by K. B. Schmidt (UCSB)
# 2014-03-13 added zrange keyword. K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import subprocess # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
import numpy as np
#----------------------------
#   FUNCTIONS
#----------------------------
def pathAname(str):                         # splitting string with path and name in to
    strsplit = str.split('/')               # splitting string
    name     = strsplit[-1]                 # saving filename (last entry of catsplit)
    slash    = '/'
    path=slash.join(strsplit[0:-1])+'/'     # putting path back together
    return [path,name]
def DandTstr():                             # Creating a string with date and time on the format yymmddhhmmss
    import datetime
    now    = datetime.datetime.now()
    now    = str(now)
    now0   = now.split(' ')
    date   = now0[0].split('-')
    date1  = date[0].split('20')
    time   = now0[1].split(':')
    HHMM   = ''.join(time[0:2])
    SS     = time[2].split('.')
    HHMMSS = HHMM+str(SS[0])
    DandT = str(date1[1])+''.join(date[1:3])+HHMMSS
    return str(DandT)
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("bpzcat", type=str, help="Path and name of bpz output catalog, i.e., *_bpz.cat")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--eps", action="store_true", help="Turn plots into eps files")
parser.add_argument("--png", action="store_true", help="Turn plots into png files and stich them togehter to aniated gif")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")
parser.add_argument("--Nzbin", type=int, help="number of redshift bins to use")
parser.add_argument("--zrange", type=float, nargs=2, help="Range of redshifts to bin [zmin,zmax]. Default is [0,10]")
parser.add_argument("--savecounts", action="store_true", help="Saving the histogram counts to text file")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
pnBPZ    = pathAname(args.bpzcat)
namebase = pnBPZ[1].split('_bpz.')
namebase = namebase[0]
namebase_tex = namebase.replace('_', '\_')

# reading BPZ results
cat  = np.genfromtxt(args.bpzcat, dtype=None, comments='#')
Nobj = len(cat)
catload = True

# reading the correspodning photometric catalog if it exists
photcat = pnBPZ[0]+'/'+namebase+'.cat'
if os.path.exists(photcat): # checking if the directory excists
    if args.verbose: print ':: '+sys.argv[0]+' :: Reading the photometric catalog '+photcat
    photcat  = np.genfromtxt(photcat, dtype=None, comments='#')
    Nobjphot = len(photcat)
    photload = True
    if Nobj != Nobjphot:
        sys.exit('ERROR: The number of objects in the two catalogs are different ('+Nobj+'!='+Nobjphot+') --> ABORTING')
else:                                  # if not move it
    if args.verbose: print ':: '+sys.argv[0]+' :: No photometric catalog present to be read'
    if args.verbose: print ':: '+sys.argv[0]+' :: Was expecting to find: '+photcat
#-------------------------------------------------------------------------------------------------------------
# PLOTTING
import matplotlib.pyplot as plt 
from matplotlib  import cm
from astropy import cosmology 

plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=13)           # setting text font
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 
plt.ioff() # turn off interactive mode, i.e. prevent plotting windows from popping up

plotdir = pnBPZ[0]+'/plots'
if not os.path.exists(plotdir): # checking if a directory for plots excists
    os.mkdir(plotdir)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if photload: # only creating this figure if a photometric catalog was loaded
    plotname = plotdir+'/'+namebase+'_raANDdec'
    if args.verbose: print ':: '+sys.argv[0]+' :: Creating figure '+plotname
    fig = plt.figure()  # create a figure object
    fig.clf()                                        # clearing figure
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    sc = ax.scatter(photcat['f4'],photcat['f5'],c=cat['f1'],marker = 'o', cmap = cm.jet );
    cb = fig.colorbar(sc)
    ax.grid(True,linestyle='-',color='0.75')
    cb.set_label('$z_\mathrm{phot}$')
    ax.set_xlabel('R.A.')
    ax.set_ylabel('Dec.')
    ax.set_title(namebase_tex)

    fig.savefig(plotname+'.pdf')
    if args.eps: fig.savefig(plotname+'.eps')
    if args.png: fig.savefig(plotname+'.png')
    # show()  # draw plot on screen

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if photload: # only creating this figure if a photometric catalog was loaded
    if args.Nzbin:  # checking if number of redshift bins is specified 
        Ndeltaz = args.Nzbin
    else:
        Ndeltaz = 20   # default number of redshift bins

    if args.zrange:
        minz = args.zrange[0]
        maxz = args.zrange[1]
    else: # Default redshift range
        minz    = 0.0   # min(cat['f1'])
        maxz    = 10.0  # max(cat['f1'])
    zbinval = [minz]
    deltaz  = (maxz-minz)/Ndeltaz
    if args.png: 
        pngnames   = []
        bincenters = []
    for ii in range(Ndeltaz):
        xval = photcat['f4']
        yval = photcat['f5']
        zval = cat['f1']
        lowlim = minz+ii*deltaz
        uplim  = minz+(ii+1)*deltaz
        zbinval.append(uplim)
        ent  = ((zval>=lowlim) & (zval<=uplim)).nonzero()   # entries where redshifts are in ii'th deltaz
        if len(ent[0]) > 0:
            plotname = plotdir+'/'+namebase+'_deltaz'+str(ii)
            if args.verbose: print ':: '+sys.argv[0]+' :: Creating figure '+plotname
            fig = plt.figure()                               # create a figure object
            fig.clf()                                        # clearing figure
            ax = fig.add_subplot(1, 1, 1)                    # create an axes object in the figure
            ax.scatter(xval,yval,marker='o', facecolors='none', edgecolors='0.75',zorder=1)
            sc = ax.scatter(xval[ent],yval[ent],c=zval[ent],marker='o',cmap=cm.jet,zorder=2)
            cb = fig.colorbar(sc)
            cb.set_label('$z_\mathrm{phot}$')
            ax.set_xlabel('R.A.')
            ax.set_ylabel('Dec.')
            ax.set_title(namebase_tex+'\n'+str("%.2f" % lowlim)+' $< z_\mathrm{phot} <$ '+str("%.2f" % uplim)+'  :  '+str("%.2f" % cosmology.lookback_time(lowlim).value)+' $< t_\mathrm{lookback}/[Gyr] <$ '+str("%.2f" % cosmology.lookback_time(uplim).value))

            string = str(len(ent[0]))+' objects shown'
            ax.text(0.95, 0.05, string, horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)

            fig.savefig(plotname+'.pdf')
            if args.eps: fig.savefig(plotname+'.eps')
            if args.png: 
                fig.savefig(plotname+'.png')
                pngnames.append(plotname+'.png')      # save filename for movie creation
                bincenters.append(lowlim+deltaz/2.)    # storing bin centers for histogram gif (below)
    if args.png:
        if len(pngnames) > 0:  # cheking that images have actually been created
            if args.verbose: print ':: '+sys.argv[0]+' :: Turning images into animated gif'
            import Image
            import numpy
            import images2gif

            imagelist  = []                        # list of images to turn into gif
            durlist    = []                        # list of duration of each images
            for image in pngnames:
                im      = Image.open(image)         # opening image
                imarr   = numpy.asarray(im)         # putting image into numpy arry
                imagelist.append(imarr)             # adding image to list
                durlist.append(0.5)

            images2gif.writeGif(plotdir+'/'+namebase+'_deltaz.gif',imagelist,duration=durlist)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if catload:
    if args.Nzbin:  # checking if number of redshift bins is specified 
        Nbin = args.Nzbin
    else:
        Nbin = 60   # default number of redshift bins
    plotname = plotdir+'/'+namebase+'_zHIST'
    if args.verbose: print ':: '+sys.argv[0]+' :: Creating figure '+plotname
    fig = plt.figure()  # create a figure object
    fig.clf()                                        # clearing figure
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    rmin = minz
    rmax = maxz
    histvals, binvals, patches = plt.hist(cat['f1'],bins=Nbin,range=[rmin,rmax],facecolor='g',alpha=0.75,zorder=2)
    ax.set_xlabel('$z_\mathrm{phot}$')
    ax.set_ylabel('\#')
    ax.set_title(namebase_tex)

    if photload and not args.Nzbin: # plotting vertical lines at redshift slices if Nzbin not given
        for ii in range(len(zbinval)):
            ax.axvline(x=zbinval[ii],color='black',linestyle='--',zorder=1)

    fig.savefig(plotname+'.pdf')
    if args.eps: fig.savefig(plotname+'.eps')
    if args.png: fig.savefig(plotname+'.png')

    if args.savecounts: # if histogram counts are to be saved do that
        countsfile = args.bpzcat.replace('.cat','_DzCOUNTS.txt')
        cfile      = open(countsfile, 'w')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        cfile.write('# File containing the histogram counts, i.e., number counts in redshift bins \n')
        cfile.write('# for the objects in: \n')
        cfile.write('# '+args.bpzcat+' \n')
        cfile.write('# Created with plotBPZresults.py \n')
        cfile.write('# Histogram is running from '+str(rmin)+' to '+str(rmax)+' \n')
        cfile.write('# The columns are: \n')
        cfile.write('#  \n')
        cfile.write('# binmin   binmax   count   \n')
        for ii in range(len(histvals)):
            cfile.write(str(binvals[ii])+'   '+str(binvals[ii+1])+'   '+str(histvals[ii])+'  \n')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        cfile.close() 

    if args.png and photload:  # preparing gif movie
        if len(bincenters) > 0:  # cheking that images have actually been created
            if args.verbose: print ':: '+sys.argv[0]+' :: Turning histogram into animated gif'
            import Image
            import numpy
            import images2gif

            histpngs = []
            for jj in range(len(bincenters)):
                fig.clf() # clear figure to redraw with only one axvline
                # --- drawing commands from above ---
                ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
                histvals,binvals,patches=plt.hist(cat['f1'],bins=Nbin,range=[minz,maxz],facecolor='g',alpha=0.75,zorder=2)
                ax.set_xlabel('$z_\mathrm{phot}$')
                ax.set_ylabel('\#')
                ax.set_title(namebase_tex)
                #------------------------------------
                ax.axvline(x=bincenters[jj],color='black',linestyle='--',zorder=0)
                pngname = plotname+str(jj)+'.png'
                fig.savefig(pngname)
                histpngs.append(pngname)

            imagelist  = []                        # list of images to turn into gif
            durlist    = []                        # list of duration of each images
            for image in histpngs:
                im      = Image.open(image)         # opening image
                imarr   = numpy.asarray(im)         # putting image into numpy arry
                imagelist.append(imarr)             # adding image to list
                durlist.append(0.5)

            images2gif.writeGif(plotdir+'/'+namebase+'_zHIST.gif',imagelist,duration=durlist)

#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

