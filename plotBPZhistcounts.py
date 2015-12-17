#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# plotBPZhistcounts.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Plotting multiple *_DzCOUNTS.txt files created with plotBPZresults.py
#----------------------------
#   COMMENTS
#----------------------------
# Code must be run in directory where the dirlist-directories are located
#----------------------------
#   INPUTS:
#----------------------------
# dirlist          : Path AND name of file containing names of directories containing BPZ output.
#                    If the results have been created with runBPZmultiplecats.py the directories
#                    contain names including a time stamp yymmddhhmmss, e.g., 121019120545_
#                    NOTE: same input as web4BPZgifs.py and modifyBoRGcats_postBPZ.py use
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --area           : --area str provides the path and name of the borg_field_info.txt file (or a similar
#                    file) which contains information about the input catalogs. If nothing is given 
#                    it is assumed that the fields are of the same size and the objects counts can be
#                    directly compared.
# --binsize        : --binsite integer enables the user to provide the size of the histogram bins.
#                    The default is binsize 100
# --nbins          : --nbins integer defines the number of bins to use in each hostgram as an alternative
#                    to defining the binsize (i.e., nbins overwrites binsize input)
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
# bash> chmod +x plotBPZhistcounts.py       (only required once)
# bash> plotBPZhistcounts.py '/Users/kasperborelloschmidt/work/BoRG/BPZmultitest121015/resultdirectorylist.txt' --verbose --png
# bash> plotBPZhistcounts.py '/Users/kasperborelloschmidt/work/BoRG/BoRGfullrun/resultdirectorylist121019run.txt' --verbose --nbins 50 --area '/Users/kasperborelloschmidt/work/BoRG/borg_field_info.txt'
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-10-19 started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import subprocess # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
import numpy as np
import commands
#----------------------------
#   FUNCTIONS
#----------------------------
def pathAname(str):                         # splitting string with path and name in to
    strsplit = str.split('/')               # splitting string
    name     = strsplit[-1]                 # saving filename (last entry of catsplit)
    slash    = '/'
    path=slash.join(strsplit[0:-1])         # putting path back together
    return [path,name]
def gauss(x,A,mu,sigma):                           # Define gauss function:
    import numpy
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("dirlist", type=str, help="List of directories to search for *_DzCOUNTS.txt files")
# ---- optional arguments ----
parser.add_argument("--area", type=str, help="field_info.txt file with areas of fields in arcmin**2 (default assumes fields have same size)")
parser.add_argument("--binsize", type=int, help="The size of the histogram bins (default = 100)")
parser.add_argument("--nbins", type=int, help="Number of bins to use for histgrams (overwrites --binsize input)")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--eps", action="store_true", help="Turn plots into eps files")
parser.add_argument("--png", action="store_true", help="Turn plots into png files and stich them togehter to aniated gif")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# READING DATA
buffer    = open(args.dirlist).readlines()  # read file into buffer
Nfiles    = 0                               # resetting counter 
PNdirlist = pathAname(args.dirlist)
if args.area:                               # reading infofile if provided
    infodat  = np.genfromtxt(args.area, dtype=None, comments='#')   # read file into buffer
    fields, indices = np.unique(infodat['f0'], return_index=True)   # returning unique fields and indices of list
    areasALL = infodat['f3']
    areas    = areasALL[indices]                                    # getting areas corresponding to fields

for l in buffer:                            # loop over objects and read counts
    dir = l.split('\n')                     # removing returns
    dir = dir[0]         
    cfpath   = dir+'/*_DzCOUNTS.txt'        # search string
    countfile= commands.getoutput('ls '+cfpath)
    if os.path.isfile(countfile):           # if file excists read it
        cat  = np.genfromtxt(countfile, dtype=None, comments='#')
        PNcountfile = pathAname(countfile)
        catname = PNcountfile[-1].split('_multiband')

        try:
            binminALL = np.vstack([binminALL,cat['f0']])            # if array is degined stack new values
        except NameError:
            binminALL = cat['f0']                                   # creating array if not defined

        try:
            binmaxALL = np.vstack([binmaxALL,cat['f1']])            # if array is degined stack new values
        except NameError:
            binmaxALL = cat['f1']                                   # creating array if not defined

        if args.area:
            fieldarea  = areas[np.where(fields == catname[0])]
            if catname[0]=='borg_1510+1115': fieldarea=[3.4754773]  # info missin in input file... easy fix
        else:
            fieldarea = 1.0

        try:
            countsALL = np.vstack([countsALL,cat['f2']/fieldarea])  # if array is defined stack new values
        except NameError:
            countsALL = cat['f2']/fieldarea                         # creating array if not defined

        Nfiles = Nfiles + 1
    else:
        print ':: '+sys.argv[0]+' :: ERROR :: The '+countfile+' was not found in '+dir[0]

#print 'shapes: ',countsALL.shape,binminALL.shape,binmaxALL.shape
Nzbin     = countsALL[0,:].size
Nobjfield = countsALL.sum(axis=1)  # sum of each row     (# objects in each field)
Nobjbin   = countsALL.sum(axis=0)  # sum of each column  (# objects in each redshift bin)
minval    = countsALL.min(0)       # minimum value of each column
maxval    = countsALL.max(0)       # maximum value of each column

if args.verbose: print ':: '+sys.argv[0]+' :: Each catlog contains '+str(Nzbin)+' redshift bins'
#-------------------------------------------------------------------------------------------------------------
# PLOTTING
import matplotlib.pyplot as plt 
from matplotlib  import cm
from astropy import cosmology 

plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=13)           # setting text font
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 

plotdir = PNdirlist[0]+'/HISTplots'
if not os.path.exists(plotdir): # checking if a directory for plots excists
    os.mkdir(plotdir)
    if args.verbose: print ':: '+sys.argv[0]+' :: The directory '+str(plotdir)+' did not exist so created it'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
from scipy.stats import poisson as pos  # importing a poisson distribution
from scipy.optimize import curve_fit    # package to fit curves to data
import matplotlib
import matplotlib.mlab as mlab          # 

for jj in range(Nzbin):   # looping over redshift bins
#for jj in range(2):   # looping over redshift bins
    if Nobjbin[jj] >= 1:                              # only create plot if there are actually objects in redshift bin
        rmin = minval[jj]+0.0
        rmax = maxval[jj]+0.0
        if args.binsize: 
            binsize = args.binsize
        else:
            binsize = 100.00                      # setting default bin size if not given
        if args.nbins:                            # if Nbins is given instead overwrite binsize
            Nbin = args.nbins
        else:
            Nbin = np.ceil((rmax-rmin)/binsize)   # calculating number of bins (rounded up)
        # ------------------------------------------------------
        # creating histogram field-by-field to color code individual 'bricks'
        binsize    = (rmax-rmin)/Nbin
        # loading color map (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps) for coloring bricks
        colorlist  = plt.get_cmap('gist_rainbow')
        #colorlist  = plt.get_cmap('Set1')

        fig = plt.figure()                    # create a figure object
        fig.clf()                             # clearing figure
        ax = fig.add_subplot(1, 1, 1)         # create an axes object in the figure
        plotname = plotdir+'/'+PNdirlist[1].replace('.txt','_zHIST'+str(jj))
        if args.verbose: print ':: '+sys.argv[0]+' :: Creating figure '+plotname

        yoff = [0.0]*Nbin
        for kk in  xrange(Nfiles):
            bincount = countsALL[kk,jj]   # the number of objects in the jj'th redshift bin for the kk'th file
            histvec  = [0.0]*Nbin
            binlow   = []
            binhigh  = []
            for ff in range(Nbin):
                binlow.append(ff*binsize+rmin)
                binhigh.append((ff+1)*binsize+rmin) 
                if (bincount-binlow[ff]>10e-10) & (bincount-binhigh[ff]<=10e-10): histvec[ff] = 1.0
                if (ff==0) & (bincount-binlow[ff]>=0.0) & (bincount-binhigh[ff]<=0.0): histvec[ff] = 1.0 # first bin

            coluse = colorlist(kk/(Nfiles+0.0))
            plt.bar(binlow,histvec,binsize,bottom=yoff,color=coluse)
            yoff = map(sum, zip(yoff,histvec))

            string = 'field \#'+str(kk)
            xcol   = np.floor(kk/16.)
            ylines = (kk/16.-xcol)*10.0
            ax.text(0.40+xcol*0.15,0.97-ylines*0.05, string, horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,color=coluse)


        if args.area: 
            ax.set_xlabel('$N_\mathrm{objects}/\mathrm{arcmin}^2$')
        else:
            ax.set_xlabel('$N_\mathrm{objects}')
        ax.set_ylabel('$N_\mathrm{fields}$')
        ax.set_title(str(binminALL[0,jj])+' $< z <$ '+str(binmaxALL[0,jj]))

        plt.grid(True)
        plt.legend()
        fig.savefig(plotname+'.pdf')
        # ------------------------------------------------------
        # over-plotting histogram to simply fitting the gaussian to the distribution
        histinput = countsALL[:,jj]           # extracting row (number counts for in redshift for all fields)
        histvals, binvals, patches = plt.hist(histinput,bins=Nbin,range=[rmin,rmax],facecolor='none',alpha=0.75,zorder=2)#,linewidth=2.0,linestyle='dotted'
        bin_centres = (binvals[:-1] + binvals[1:])/2

        # fitting hist to poisson and gauss distribution
        mu    = np.mean(histinput)        # the mean of histogram input
        sigma = np.std(histinput)         # the standard deviation of histogram input
        var   = sigma**2                  # the variance of histogram input

        # if the variance is high (say var>1000), the normal distribution with mean=var and 
        # variance=var is a good approximation. For large k, k! also becomes crazy to calculate
        #coeff, var_matrix = curve_fit(poisson, bin_centres, histvals, p0=mu) # fitting gauss to data
        #poissonfit = poisson(bin_centres,*coeff)
        #plt.plot(bin_centres,poissonfit,'k.',linewidth=2)                    # overplotting poisson fit

        G0 = [1, mu, sigma]                                                 # initial guess (A, mu and sigma)
        try:
            coeff, var_matrix = curve_fit(gauss, bin_centres, histvals, p0=G0)   # fitting gauss to data
            gaussfit = gauss(bin_centres,*coeff)
            labelstr = 'Gauss fit ($\mu=$'+str("%.2f" % coeff[1])+', $\sigma=$'+str("%.2f" % coeff[2])+')'
            plt.plot(bin_centres,gaussfit,'k--',linewidth=2,label=labelstr)                     # overplotting gauss fit
        except RuntimeError:
            labelstr = 'NB: SciPy Gauss fit failed...'
            ax.text(0.95, 0.95,labelstr,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
        # ------------------------------------------------------
        fig.savefig(plotname+'.pdf')
        if args.eps: fig.savefig(plotname+'.eps')
        if args.png: fig.savefig(plotname+'.png')

#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

