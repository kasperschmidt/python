#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# plotBorgMultibandCatalogs.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Plotting the content of a list of BoRG *_multiband.cat catalogs
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# catalognames     : Path AND name of file containing path and names of files to plot.
#                    Code will loop over the individual files and assign different colors to
#                    each catalog's points.
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --verbose        : set -verbose to get info/messages printed to the screen
# --eps            : saving created plots as eps files
# --png            : saving created plots as png files
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# plots will be saved in the same directory as the catalognames file in a sub-directory called plots
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x plotBorgMultibandCatalogs.py       (only required once)
# bash> plotBorgMultibandCatalogs.py '/Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p2/multiband/BoRGcatalognames.txt' --verbose 
#
#
# bash> plotBorgMultibandCatalogs.py '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121108/modifiedcatnames.txt' --verbose 
#
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-05-13  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import commands   # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
import numpy as np# enable opening with genfromtxt
import matplotlib.pyplot as plt
#----------------------------
#   FUNCTIONS
#----------------------------
def pathAname(str):                         # splitting string with path and name in to
    strsplit = str.split('/')               # splitting string
    name     = strsplit[-1]                 # saving filename (last entry of catsplit)
    slash    = '/'
    path=slash.join(strsplit[0:-1])         # putting path back together
    return [path,name]
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("catalognames", type=str, help="Path and name of file containing list of catalog names")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--eps", action="store_true", help="Turn plots into eps files")
parser.add_argument("--png", action="store_true", help="Turn plots into png files")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# Reading ascii
pnCATS = pathAname(args.catalognames)
cats   = np.genfromtxt(args.catalognames,dtype=None,comments='#')
try: 
    len(cats)
    Ncat   = len(cats)
except:
    Ncat   = 1
    cats   = np.append(cats, cats) # creating array with two elements to be able to index

if args.verbose: print ':: '+sys.argv[0]+' :: Found '+str(Ncat)+' catalogs to plot'
#-------------------------------------------------------------------------------------------------------------
# PLOTTING
plotdir = pnCATS[0]+'/plots'
if not os.path.exists(plotdir): # checking if a directory for plots excists
    os.mkdir(plotdir)
if args.verbose: print ':: '+sys.argv[0]+' :: Plots will be save to '+str(plotdir)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
plotname = plotdir+'/colorlegend'
if args.verbose: print ':: '+sys.argv[0]+' :: Creating figure '+plotname
plt.rc('text', usetex=False)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=10)           # setting text font
plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10) 
fig = plt.figure()  # create a figure object
fig.clf()                                        # clearing figure
ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
colorlist  = plt.get_cmap('gist_rainbow') # loading color map (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
plt.xlim((0,1))
plt.ylim((0,1))

for ii in range(Ncat):  # looping over catalogs
    pndat  = pathAname(cats[ii])
    coluse = colorlist(ii/(Ncat+0.0))
    string = pndat[1]
    xcol   = np.floor(ii/25)
    ylines = (ii/25.-xcol)*10.0
    ax.text(0.0+xcol*0.33,0.97-ylines*0.10, string, horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,color=coluse)

fig.savefig(plotname+'.pdf')
if args.eps: fig.savefig(plotname+'.eps')
if args.png: fig.savefig(plotname+'.png')


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
plotname = plotdir+'/stellarityVSf125isomag'
if args.verbose: print ':: '+sys.argv[0]+' :: Creating figure '+plotname
plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=13)           # setting text font
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 
fig = plt.figure()  # create a figure object
fig.clf()                                        # clearing figure
ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
colorlist  = plt.get_cmap('gist_rainbow') # loading color map (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
plt.xlim((30,10))
plt.ylim((-0.1,1.1))
pointsize = 5

for ii in range(Ncat):  # looping over catalogs
    dat    = np.genfromtxt(cats[ii],dtype=None,comments='#')  # loading data
    xval   = dat['f43']
    yval   = dat['f12']
    coluse = colorlist(ii/(Ncat+0.0))

    ax.scatter(xval,yval,marker='o',facecolors=coluse,edgecolors='none',s=pointsize,zorder=1)

ax.grid(True,linestyle='-',color='0.75')
ax.set_xlabel('ISOMAG$_\mathrm{F125W}$')
ax.set_ylabel('Stellarity')
#ax.set_title(namebase_tex)
# show()  # draw plot on screen

fig.savefig(plotname+'.pdf')
if args.eps: fig.savefig(plotname+'.eps')
if args.png: fig.savefig(plotname+'.png')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
plotname = plotdir+'/stellarityVSf125automag'
if args.verbose: print ':: '+sys.argv[0]+' :: Creating figure '+plotname
plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=13)           # setting text font
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 
fig = plt.figure()  # create a figure object
fig.clf()                                        # clearing figure
ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
colorlist  = plt.get_cmap('gist_rainbow') # loading color map (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
plt.xlim((30,10))
plt.ylim((-0.1,1.1))
pointsize = 5

for ii in range(Ncat):  # looping over catalogs
    dat    = np.genfromtxt(cats[ii],dtype=None,comments='#')  # loading data
    xval   = dat['f45']
    yval   = dat['f12']
    coluse = colorlist(ii/(Ncat+0.0))

    ax.scatter(xval,yval,marker='o',facecolors=coluse,edgecolors='none',s=pointsize,zorder=1)

ax.grid(True,linestyle='-',color='0.75')
ax.set_xlabel('automag$_\mathrm{F125W}$')
ax.set_ylabel('Stellarity')
#ax.set_title(namebase_tex)
# show()  # draw plot on screen

fig.savefig(plotname+'.pdf')
if args.eps: fig.savefig(plotname+'.eps')
if args.png: fig.savefig(plotname+'.png')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
plotname = plotdir+'/YmJvsJmH_isomag'
if args.verbose: print ':: '+sys.argv[0]+' :: Creating figure '+plotname
plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=13)           # setting text font
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 
fig = plt.figure()  # create a figure object
fig.clf()                                        # clearing figure
ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
colorlist  = plt.get_cmap('gist_rainbow') # loading color map (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
Xrange = [-1,1]
Yrange = [1,3.5]
plt.xlim((Xrange[0],Xrange[1]))
plt.ylim((Yrange[0],Yrange[1]))

for ii in range(Ncat):  # looping over catalogs
    dat    = np.genfromtxt(cats[ii],dtype=None,comments='#')  # loading data
    xval   = dat['f43']-dat['f54']
    yval   = dat['f32']-dat['f43']
    pointsize = dat['f39']
    coluse = colorlist(ii/(Ncat+0.0))

    ax.scatter(xval,yval,marker='o',facecolors=coluse,edgecolors='none',s=pointsize,zorder=1)

plt.plot([Xrange[0],0.02],[1.75,1.75],color='k',linestyle='-')
plt.plot([0.02,0.32],[1.75,3.75],color='k',linestyle='-')

plt.plot([Xrange[0],0.5],[1.75,1.75],color='k',linestyle='--')
plt.plot([0.5,0.5],[1.75,Yrange[1]],color='k',linestyle='--')

# legend by hand
ax.text(0.92,0.95,'(S/N)$_\mathrm{J}$ = 5',horizontalalignment='right',verticalalignment='center',transform=ax.transAxes,color='k')
ax.scatter([Xrange[0]+(Xrange[1]-Xrange[0])*0.95]*2,[Yrange[0]+(Yrange[1]-Yrange[0])*0.95]*2,facecolors='k',marker='o',s=5)

ax.text(0.92,0.90,'(S/N)$_\mathrm{J}$ = 8',horizontalalignment='right',verticalalignment='center',transform=ax.transAxes,color='k')
ax.scatter([Xrange[0]+(Xrange[1]-Xrange[0])*0.95]*2,[Yrange[0]+(Yrange[1]-Yrange[0])*0.90]*2,facecolors='k',marker='o',s=8)


ax.grid(True,linestyle='-',color='0.75')
ax.set_xlabel('J - H (ISOMAG)')
ax.set_ylabel('Y - J (ISOMAG)')
#ax.set_title(namebase_tex)
# show()  # draw plot on screen

fig.savefig(plotname+'.pdf')
if args.eps: fig.savefig(plotname+'.eps')
if args.png: fig.savefig(plotname+'.png')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
plotname = plotdir+'/YmJvsJmH_automag'
if args.verbose: print ':: '+sys.argv[0]+' :: Creating figure '+plotname
plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=13)           # setting text font
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 
fig = plt.figure()  # create a figure object
fig.clf()                                        # clearing figure
ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
colorlist  = plt.get_cmap('gist_rainbow') # loading color map (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
Xrange = [-1,1]
Yrange = [1,3.5]
plt.xlim((Xrange[0],Xrange[1]))
plt.ylim((Yrange[0],Yrange[1]))

for ii in range(Ncat):  # looping over catalogs
    dat    = np.genfromtxt(cats[ii],dtype=None,comments='#')  # loading data
    xval   = dat['f45']-dat['f56']
    yval   = dat['f34']-dat['f45']
    pointsize = dat['f42']
    coluse = colorlist(ii/(Ncat+0.0))

    ax.scatter(xval,yval,marker='o',facecolors=coluse,edgecolors='none',s=pointsize,zorder=1)

plt.plot([Xrange[0],0.02],[1.75,1.75],color='k',linestyle='-')
plt.plot([0.02,0.32],[1.75,3.75],color='k',linestyle='-')

plt.plot([Xrange[0],0.5],[1.75,1.75],color='k',linestyle='--')
plt.plot([0.5,0.5],[1.75,Yrange[1]],color='k',linestyle='--')

# legend by hand
ax.text(0.92,0.95,'(S/N)$_\mathrm{J auto}$ = 5',horizontalalignment='right',verticalalignment='center',transform=ax.transAxes,color='k')
ax.scatter([Xrange[0]+(Xrange[1]-Xrange[0])*0.95]*2,[Yrange[0]+(Yrange[1]-Yrange[0])*0.95]*2,facecolors='k',marker='o',s=5)

ax.text(0.92,0.90,'(S/N)$_\mathrm{J auto}$ = 8',horizontalalignment='right',verticalalignment='center',transform=ax.transAxes,color='k')
ax.scatter([Xrange[0]+(Xrange[1]-Xrange[0])*0.95]*2,[Yrange[0]+(Yrange[1]-Yrange[0])*0.90]*2,facecolors='k',marker='o',s=8)


ax.grid(True,linestyle='-',color='0.75')
ax.set_xlabel('J - H (automag)')
ax.set_ylabel('Y - J (automag)')
#ax.set_title(namebase_tex)
# show()  # draw plot on screen

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

