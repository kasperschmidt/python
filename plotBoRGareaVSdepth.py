#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# plotBoRGareaVSdepth.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Plotting an updated version of Figure 1 in Bradley et al. 2012
# And figure of literature alpha values vs recshift
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
# bash> plotBoRGareaVSdepth.py --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-08-25  started by K. B. Schmidt (UCSB)
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
import datetime            # print currecnt date and time with datetime.datetime.now()
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
#parser.add_argument("catlist", type=str, help="Provide path and name of file containing list of catalog names")
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

#Hippies:
xHIPPIES = [27.05]
yHIPPIES = [np.log10(122.8)]

#Candels:
#xCANDELS = [27.1,27.9,27.85,27.3]
#yCANDELS = [2.4,2.0,1.78,1.44]
xCANDELS = [27.1,27.9]
yCANDELS = [np.log10(260),np.log10(120)]

#HUDF09
xUDF = [29.0]
yUDF = [np.log10(14)]

#HUDF12
xUDF12 = [29.5]
yUDF12 = [np.log10(4.5)]

#ERS
xERS = [27.5]
yERS = [np.log10(40)]

#[BoRG09,BoRG12,BoRG13,BoRG09+12+13]
B09lim = [25.9,26.6,26.8,26.4,26.6,26.6,27.1,26.6,26.6,26.5,26.1,27.0,26.7,26.8,26.5,26.2,26.7,26.3,27.1,26.7,27.2,26.0,26.9,26.7,26.6,27.2,27.0,26.5,26.6]
B12lim = [26.9,26.8,26.7,25.1,27.3,27.2,27.0,26.7,27.1,26.8,26.2,27.1,26.8,27.0,27.2,27.0,26.8,26.8,26.1,26.9,27.1,27.0,27.0,21.6,26.4,26.4,26.7,26.8,26.8,27.3]
B13lim = [26.74,26.52,26.77,26.83,27.17,26.83,27.32,26.85,27.00,27.04,27.74,27.54,27.04,27.67,26.23,27.11]
Balllim = [25.9,26.6,26.8,26.4,26.6,26.6,27.1,26.6,26.6,26.5,26.1,27.0,26.7,26.8,26.5,26.2,26.7,26.3,27.1,26.7,27.2,26.0,26.9,26.7,26.6,27.2,27.0,26.5,26.6,26.9,26.8,26.7,25.1,27.3,27.2,27.0,26.7,27.1,26.8,26.2,27.1,26.8,27.0,27.2,27.0,26.8,26.8,26.1,26.9,27.1,27.0,27.0,21.6,26.4,26.4,26.7,26.8,26.8,27.3,26.74,26.52,26.77,26.83,27.17,26.83,27.32,26.85,27.00,27.04,27.74,27.54,27.04,27.67,26.23,27.11]

xBoRG = [np.median(B09lim),np.median(B12lim),np.median(B13lim),np.median(Balllim)]
yBoRG = [np.log10(135),np.log10(139),np.log10(70),np.log10(350)]

labLT7 = [       'Oesch+ 10',                'Reddy \& Steidel 09',   'Bouwens+ 07',  'Bouwens+ 12',  'Bouwens+ 11', 'McLure+ 13', 'Schencker+ 13']
zLT7   = [0.75,1.25,1.75,1.5,1.9,2.5,              2.3,3.05,            3.8,5.0,5.9,     5.0,6.0,      7.0 , 7.0 ,7.0]
aLT7   = [-1.52,-1.84,-1.60,-1.46,-1.60,-1.75,    -1.73,-1.73,       -1.76,-1.69,-1.77, -1.79,-1.73,  -2.01, -1.90, -1.87]
dapLT7 = [0.25,0.15,0.21,0.54,0.51,0.11,           0.07,0.13,          0.05,0.09,0.16,   0.12,0.20,    0.21,0.14,0.18]
damLT7 = [0.25,0.15,0.21,0.54,0.51,0.11,           0.07,0.13,          0.05,0.09,0.16,   0.12,0.20,    0.21,0.15,0.17]

lab8 = ['Bradley+12','Oesch+12','Bouwens+11','Schenker+13','McLure+13','BoRG13 \n(This work)','This work 8$\sigma$']
z8   = [7.75,7.85,8.15,8.25,8.35,         7.95,8.05 ] # old 2   7.95,8.05 ] # old  7.95,8.05 ]
a8   = [-1.98,-2.06,-1.91,-1.94,-2.02,   -1.87,-2.08] # old 2  -1.71,-2.07] # old -2.02,-2.15]
dap8 = [0.23,0.45,0.32,0.21,0.24,         0.26, 0.30] # old 2   0.26,0.29 ] # old 0.28,0.29 ]
dam8 = [0.22,0.37,0.32,0.24,0.23,         0.26, 0.29] # old 2   0.25,0.28 ] # old 0.28,0.28 ]

#-------------------------------------------------------------------------------------------------------------
# PLOTTING
import matplotlib.pyplot as plt 
Fsize = 18
plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=Fsize)           # setting text font
plt.rc('xtick', labelsize=Fsize) 
plt.rc('ytick', labelsize=Fsize) 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
plotname = 'BoRGareaVSdepth.pdf'
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting area VS depth figure to',plotname
fig = plt.figure()  # create a figure object
fig.clf()                                        # clearing figure
ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure: subplot(numRows, numCols, plotNum)
#ax.grid(True,linestyle='-',color='0.75',zorder=-1)

#----------------------------------------------
#plt.plot(xHIPPIES[0],yHIPPIES[0],'mo',zorder=5,markersize=10)
#plt.text(xHIPPIES[0]+0.1,yHIPPIES[0],'HIPPIES',color='magenta',
#         fontsize=Fsize,zorder=5,verticalalignment='top', horizontalalignment='left')

#----------------------------------------------
plt.plot(xCANDELS[0],yCANDELS[0],'bo',zorder=5,markersize=10)
plt.text(xCANDELS[0]+0.1,yCANDELS[0],'CANDELS wide',color='blue',
         fontsize=Fsize,zorder=5,verticalalignment='center', horizontalalignment='left')
#----------------------------------------------
plt.plot(xCANDELS[1],yCANDELS[1],'ro',zorder=5,markersize=10)
plt.text(xCANDELS[1]+0.1,yCANDELS[1],'CANDELS deep',color='red',
         fontsize=Fsize,zorder=5,verticalalignment='center', horizontalalignment='left')
#----------------------------------------------
#plt.plot(xCANDELS[2],yCANDELS[2],'ro',zorder=5,markersize=10)
#plt.text(xCANDELS[2]+0.1,yCANDELS[2],'CANDELS2',color='red',
#         fontsize=Fsize,zorder=5,verticalalignment='center', horizontalalignment='left')
#----------------------------------------------
#plt.plot(xCANDELS[3],yCANDELS[3],'ro',zorder=5,markersize=10)
#plt.text(xCANDELS[3]+0.1,yCANDELS[3],'CANDELS3',color='red',
#         fontsize=Fsize,zorder=5,verticalalignment='center', horizontalalignment='left')

#----------------------------------------------
plt.plot(xUDF[0],yUDF[0],'go',zorder=5,markersize=10)
plt.text(xUDF[0]-0.1,yUDF[0],'HUDF09',color='green',
         fontsize=Fsize,zorder=5,verticalalignment='center', horizontalalignment='right')

#----------------------------------------------
plt.plot(xERS[0],yERS[0],color='0.5',marker='o',zorder=5,markersize=10)
plt.text(xERS[0]-0.1,yERS[0],'ERS',color='0.5',
         fontsize=Fsize,zorder=5,verticalalignment='center', horizontalalignment='right')

#----------------------------------------------
plt.plot(xBoRG[0],yBoRG[0],'ko',zorder=5,markersize=10)
plt.text(xBoRG[0]-0.1,yBoRG[0],'BoRG09',color='black',
         fontsize=Fsize,zorder=5,verticalalignment='center', horizontalalignment='right')
#----------------------------------------------
plt.plot(xBoRG[1],yBoRG[1],'ko',zorder=5,markersize=10)
plt.text(xBoRG[1]+0.1,yBoRG[1],'BoRG12',color='black',
         fontsize=Fsize,zorder=5,verticalalignment='bottom', horizontalalignment='left')
#----------------------------------------------
plt.plot(xBoRG[2],yBoRG[2],'ko',zorder=5,markersize=10)
plt.text(xBoRG[2]-0.1,yBoRG[2],'BoRG13',color='black',
         fontsize=Fsize,zorder=5,verticalalignment='center', horizontalalignment='right')
#----------------------------------------------
plt.plot(xBoRG[3],yBoRG[3],'ko',zorder=5,markersize=10)
plt.text(xBoRG[3]+0.1,yBoRG[3],'BoRG09+12+13',color='black',
         fontsize=Fsize,zorder=5,verticalalignment='center', horizontalalignment='left')

ax.set_xlim(26,29.1)
ax.set_ylim(1,2.6)

ax.set_xlabel('F125W 5$\sigma$ limiting mag (AB)', fontsize=Fsize)
ax.set_ylabel('log(Area/arcmin$^2$)', fontsize=Fsize)
#ax.set_title('title of plot')

#leg = plt.legend(fancybox=True, loc='upper left')  # add the legend in the middle of the plot
#leg.get_frame().set_alpha(0.7)    

fig.savefig(plotname)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
plotname = 'alphaVSz.pdf'
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting alpha VS redshift figure to',plotname
fig = plt.figure()  # create a figure object
fig.clf()                                        # clearing figure

ylimits = [-2.5,-1.3]
# --------------------- x < 7 ---------------------
fig.subplots_adjust(wspace=.02,top=0.95)
ax = fig.add_subplot(1, 2, 1)  # create an axes object in the figure: subplot(numRows, numCols, plotNum)
#ax.grid(True,linestyle='-',color='0.75',zorder=-1)

plt.errorbar(zLT7[0:6],aLT7[0:6],yerr=[damLT7[0:6],dapLT7[0:6]]        ,fmt='bo',zorder=5,markersize=7)
plt.errorbar(zLT7[6:8],aLT7[6:8],yerr=[damLT7[6:8],dapLT7[6:8]]        ,fmt='ro',zorder=5,markersize=7)
plt.errorbar(zLT7[8:11],aLT7[8:11],yerr=[damLT7[8:11],dapLT7[8:11]]    ,fmt='mo',zorder=5,markersize=7)
plt.errorbar(zLT7[11:13],aLT7[11:13],yerr=[damLT7[11:13],dapLT7[11:13]],fmt='go',zorder=5,markersize=7)
plt.errorbar(zLT7[13:],aLT7[13:],yerr=[damLT7[13:],dapLT7[13:]]        ,fmt='ko',zorder=5,markersize=7)

plt.text(1.2,-2.47,labLT7[0],color='blue',fontsize=Fsize,verticalalignment='bottom',horizontalalignment='center',rotation='vertical')
plt.text(2.9,-2.47,labLT7[1],color='red',fontsize=Fsize,verticalalignment='bottom',horizontalalignment='center',rotation='vertical')
plt.text(4.9,-2.47,labLT7[2],color='magenta',fontsize=Fsize,verticalalignment='bottom',horizontalalignment='center',rotation='vertical')
plt.text(5.5,-2.47,labLT7[3],color='green',fontsize=Fsize,verticalalignment='bottom',horizontalalignment='center',rotation='vertical')
plt.text(7.0,-2.47,labLT7[4],color='black',fontsize=Fsize,verticalalignment='bottom',horizontalalignment='center',rotation='vertical')



ax.set_xlim(0,7.5)
ax.set_ylim(ylimits)

ax.set_xlabel('$z$', fontsize=Fsize)
ax.set_ylabel('$\\alpha$', fontsize=Fsize)

# --------------------- x = 8 ---------------------
ax = fig.add_subplot(1, 2, 2)  # create an axes object in the figure: subplot(numRows, numCols, plotNum)
#ax.grid(True,linestyle='-',color='0.75',zorder=-1)

plt.errorbar(z8[0:5],a8[0:5], yerr=[dam8[0:5],dap8[0:5]],fmt='ko',zorder=5,markersize=7)
plt.errorbar(z8[5:],a8[5:], yerr=[dam8[5:],dap8[5:]],fmt='ks',zorder=5,markersize=10)

for ii in xrange(len(lab8)):
    plt.text(z8[ii],-1.57,lab8[ii],color='black',fontsize=Fsize,zorder=5,verticalalignment='bottom', 
             horizontalalignment='center',rotation='vertical')

ax.set_yticklabels('0')
ax.set_xticklabels('8.0')
ax.set_xticks([8.0])

ax.set_xlim(7.7,8.4)
ax.set_ylim(ylimits)

#ax.set_title('title of plot')

#leg = plt.legend(fancybox=True, loc='upper left')  # add the legend in the middle of the plot
#leg.get_frame().set_alpha(0.7)    

fig.savefig(plotname)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
plotname = 'alphaVSz_wtheory.pdf'
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting alpha VS redshift with theory figure to',plotname
fig = plt.figure()  # create a figure object
fig.clf()                                        # clearing figure

ylimits = [-2.32,-0.9]
ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure: subplot(numRows, numCols, plotNum)
#ax.grid(True,linestyle='-',color='0.75',zorder=-1)

theory_alpha  = np.array([-1.29,-1.63,-1.60,-1.68,-1.73,-1.77,-1.76,-1.84,-1.92,-2.18])
theory_alphap = np.array([0.05,0.04,0.04,0.05,0.07,0.11,0.14,0.12,0.11,0.25])
theory_alpham = np.array([0.05,0.02,0.06,0.07,0.05,0.05,0.12,0.17,0.15,0.02])
theory_z      = np.array([0.3,1,2,3,4,5,6,7,8,10])

plt.fill_between(theory_z,theory_alpha-theory_alpham,theory_alpha+theory_alphap,alpha=0.4,color='k')
plt.plot(np.zeros(len(a8))+8.0,a8,'ko',markersize=12,mfc='None',markeredgewidth=1.5)

plt.errorbar(zLT7[6:8],aLT7[6:8],yerr=[damLT7[6:8],dapLT7[6:8]]        ,fmt='ro',zorder=5,markersize=12,markeredgewidth=1.5)
plt.errorbar(zLT7[0:6],aLT7[0:6],yerr=[damLT7[0:6],dapLT7[0:6]]        ,fmt='bo',zorder=5,markersize=12,markeredgewidth=1.5)
plt.errorbar(zLT7[8:11],aLT7[8:11],yerr=[damLT7[8:11],dapLT7[8:11]]    ,fmt='mo',zorder=5,markersize=12,markeredgewidth=1.5)
plt.errorbar(zLT7[11:13],aLT7[11:13],yerr=[damLT7[11:13],dapLT7[11:13]],fmt='go',zorder=5,markersize=12,markeredgewidth=1.5)
plt.errorbar(zLT7[13],aLT7[13],yerr=np.max([damLT7[13],dapLT7[13]])    ,fmt='co',zorder=5,markersize=12,markeredgewidth=1.5)
plt.errorbar(zLT7[14],aLT7[14],yerr=np.max([damLT7[14],dapLT7[14]])    ,
             fmt='o',zorder=5,markersize=12,color='#663399',markeredgewidth=1.5)
plt.errorbar(zLT7[15],aLT7[15],yerr=np.max([damLT7[15],dapLT7[15]])    ,
             fmt='o',zorder=5,markersize=12,color='#FF6600',markeredgewidth=1.5)

plt.text(1.2,-0.92,labLT7[0],color='blue',fontsize=Fsize,verticalalignment='top',horizontalalignment='center',rotation='vertical')
plt.text(2.9,-0.92,labLT7[1],color='red',fontsize=Fsize,verticalalignment='top',horizontalalignment='center',rotation='vertical')
plt.text(4.9,-0.92,labLT7[2],color='magenta',fontsize=Fsize,verticalalignment='top',horizontalalignment='center',rotation='vertical')
plt.text(5.5,-0.92,labLT7[3],color='green',fontsize=Fsize,verticalalignment='top',horizontalalignment='center',rotation='vertical')
plt.text(6.6,-0.92,labLT7[4],color='cyan',fontsize=Fsize,verticalalignment='top',horizontalalignment='center',rotation='vertical')
plt.text(7.0,-0.92,labLT7[5],color='#663399',fontsize=Fsize,verticalalignment='top',horizontalalignment='center',rotation='vertical')
plt.text(7.4,-0.92,labLT7[6],color='#FF6600',fontsize=Fsize,verticalalignment='top',horizontalalignment='center',rotation='vertical')

# - - - - - This work point - - - - -
plt.text(8.3,-0.92,lab8[5],color='black',fontsize=Fsize,zorder=5,verticalalignment='top',horizontalalignment='center',rotation='vertical')
plt.errorbar(8.0,a8[5], yerr=dap8[5],fmt='ks',zorder=5,markersize=15,markeredgewidth=1.5)
# - - - - - - - - - - - - - - - - - -

plt.text(3.0,-2.2,'Tacchella+ 13 model',color='black',alpha=0.6,fontsize=Fsize,zorder=5,verticalalignment='top',horizontalalignment='center')

ax.set_xlim(0,10.5)
ax.set_ylim(ylimits)

ax.set_xlabel('$z$', fontsize=Fsize)
ax.set_ylabel('$\\alpha$', fontsize=Fsize)

fig.savefig(plotname)

#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
