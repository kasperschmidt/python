#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# plotTrenti2011overdensity.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Plotting the Trenti et al. 2011 overdensity objects
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
# bash> plotTrenti2011overdensity.py --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-09-10  started by K. B. Schmidt (UCSB)
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
# conversiong from cps [e/s] to Flampda [erg/cm2/A/s] from fits headers
PHOTFLAM_098 = 6.05015669999999E-20 #[erg/cm2/A/e]
PHOTFLAM_105 =        3.0385782E-20 #[erg/cm2/A/e]
PHOTFLAM_125 =        2.2483504E-20 #[erg/cm2/A/e]
PHOTFLAM_160 =        1.9275707E-20 #[erg/cm2/A/e]
PHOTFLAM_606 =        1.1706338E-19 #[erg/cm2/A/e]

# aperture values for object T12a,T12b,T12c,T12d,T12e
plotsymb          = ['ko-','ro-','go-','mo-','bo-']
plotlabel         = ['T12a','T12b','T12c','T12d','T12e']
apervalBoRG13_606 = np.array([0.01964,0.00909,0.14869,0.03027,0.03643])
apervalBoRG13_098 = np.array([0.02520,-0.00206,0.08568,0.00027,0.06101])
apervalBoRG13_105 = np.array([0.75906,-0.15390,0.19096,-0.06253,0.14556])
apervalBoRG13_125 = np.array([1.16491,0.00367,0.23813,0.01042,0.43950])
apervalBoRG13_160 = np.array([1.03015,0.04636,0.22203,-0.06992,0.34926])

ZP606, ZP098, ZP105, ZP125, ZP160 = 26.0814041237, 25.68102452, 26.27, 26.2473632068, 25.9558992372 # from borg_fldinfo.py
#ZP606, ZP098, ZP105, ZP125, ZP160 = 26.0440, 25.6650, 26.2554, 26.2363, 25.9480 # from dirlist_deredzeropoints.txt 

# getting 1sigma fluxes in bands from 5sigma limiting mags
Fonesig606 = np.exp( (27.7006676315 - 2.5*np.log10(1./5.) - ZP606)*np.log(10)/-2.5 )
Fonesig098 = np.exp( (27.6543699199 - 2.5*np.log10(1./5.) - ZP098)*np.log(10)/-2.5 )
#onesig105 =  - 2.5*np.log10(1./5.)
Fonesig125 = np.exp( (27.7391672948 - 2.5*np.log10(1./5.) - ZP125)*np.log(10)/-2.5 )
Fonesig160 = np.exp( (27.4884069583 - 2.5*np.log10(1./5.) - ZP160)*np.log(10)/-2.5 )

aper606onesig = apervalBoRG13_606.copy()
aper098onesig = apervalBoRG13_098.copy()
#aper105onesig = apervalBoRG13_105.copy()
aper125onesig = apervalBoRG13_125.copy()
aper160onesig = apervalBoRG13_160.copy()

# Setting fluxes to 1sigma values if smaller
aper606onesig[np.where( apervalBoRG13_606 < Fonesig606)] = Fonesig606
aper098onesig[np.where( apervalBoRG13_098 < Fonesig098)] = Fonesig098
#aper105onesig = 
aper125onesig[np.where( apervalBoRG13_125 < Fonesig125)] = Fonesig125
aper160onesig[np.where( apervalBoRG13_160 < Fonesig160)] = Fonesig160

YJcolBoRG13 = -2.5 * np.log10( (aper098onesig*PHOTFLAM_098) / (aper125onesig*PHOTFLAM_125) ) + ZP098 - ZP125
JHcolBoRG13 = -2.5 * np.log10( (aper125onesig*PHOTFLAM_125) / (aper160onesig*PHOTFLAM_160) ) + ZP125 - ZP160

YJcolBoRG12 = np.array([2.7,1.8,1.8,1.8,2.0])
JHcolBoRG12 = np.array([0.0,-0.3,-0.7,0.0,-0.4])

V_606_BoRG09 = np.array([28.6,29.0,28.3,28.5,29.0]) # all lower limits 
Y_098_BoRG09 = np.array([28.8,29.4,29.3,29.2,29.3]) # all lower limits 
J_125_BoRG09 = np.array([25.9,27.5,27.3,27.4,27.2])
H_160_BoRG09 = np.array([26.0,27.8,27.6,27.7,27.6])

YJcolBoRG09 = Y_098_BoRG09-J_125_BoRG09
JHcolBoRG09 = J_125_BoRG09-H_160_BoRG09

#-------------------------------------------------------------------------------------------------------------
# PLOTTING
import matplotlib.pyplot as plt 
Fsize = 15
plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=Fsize)           # setting text font
plt.rc('xtick', labelsize=Fsize) 
plt.rc('ytick', labelsize=Fsize) 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
plotname = 'Trenti2011overdensity_YJJHspace.pdf'
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting the objects in YJ-JH color space to ',plotname
fig = plt.figure()  # create a figure object
fig.clf()                                        # clearing figure
ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure: subplot(numRows, numCols, plotNum)
#ax.grid(True,linestyle='-',color='0.75',zorder=-1)

xmin, xmax, ymin, ymax = -1.0, 1, -1.5, 3.4

for ii in xrange(5): # looping over objects
    plt.plot([JHcolBoRG13[ii],JHcolBoRG12[ii],JHcolBoRG09[ii]],
             [YJcolBoRG13[ii],YJcolBoRG12[ii],YJcolBoRG09[ii]],plotsymb[ii],label=plotlabel[ii])


plt.plot([xmin,0.02],[1.75,1.75],'k--',label=r'color selection')
plt.plot([0.02,0.02+0.15*(ymax-1.75)],[1.75,ymax],'k--')

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(r'$J_{125} - H_{160}$', fontsize=Fsize)
ax.set_ylabel(r'$Y_{098} - J_{125}$', fontsize=Fsize)
#ax.set_title('title of plot')

leg = plt.legend(fancybox=True, loc='upper right')  # add the legend in the middle of the plot
leg.get_frame().set_alpha(0.7)    

fig.savefig(plotname)
if args.show: plot.show()  # draw plot on screen   

#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
