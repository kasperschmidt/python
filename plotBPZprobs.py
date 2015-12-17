#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# plotBPZprobs.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Plotting the probability distributions from BPZ 
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# probfile         : File containing probability distributions outputted by the BPZ run
# outputdir        : directory to save the plots in (assuming it has already been created) - end with '/'
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --skip           : provide number of objects to skip, i.e. from 0 to N
#                    NB! note that using this keyword has the consequence that ./html/probplots/probplotALL
#                        only contains the N to Ntotal objects, i.e., not the skipped objects.
# --verbose        : set -verbose to get info/messages printed to the screen
# --eps            : saving created plots as eps files
# --pdf            : saving created plots as pdf files
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
#
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x plotBPZprobs.py       (only required once)
# bash> plotBPZprobs.py 'BPZdir/121108120242_BPZrun_borg_2351-4332_multiband_byhand_BPZinput/borg_2351-4332_multiband_byhand_BPZinput.probs' 'BPZdir/121108120242_BPZrun_borg_2351-4332_multiband_byhand_BPZinput/' --verbose
#
# bash> plotBPZprobs.py 'BPZdir/121108133927_BPZrun_borg_2351-4332_multiband_byhand_BPZinput/borg_2351-4332_multiband_byhand_BPZinput.probs' 'BPZdir/121108133927_BPZrun_borg_2351-4332_multiband_byhand_BPZinput/probplotsNEW/' --verbose
#
# --- in BoRG/modifiedBoRGcats121108BPZ/ directory ---
# bash> plotBPZprobs.py '121109005728_BPZrun_borg_2351-4332_multiband_modified_BPZinput/borg_2351-4332_multiband_modified_BPZinput.probs' '121109005728_BPZrun_borg_2351-4332_multiband_modified_BPZinput/html/probplots/' --verbose
#
# bash> plotBPZprobs.py '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121108BPZ/121108223229_BPZrun_borg_0751+2917_multiband_modified_BPZinput/borg_0751+2917_multiband_modified_BPZinput.probs' '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121108BPZ/121108223229_BPZrun_borg_0751+2917_multiband_modified_BPZinput/html/probplots/' --verbose
#
# bash> plotBPZprobs.py '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121108BPZ/12110823279_BPZrun_borg_1033+5051_multiband_modified_BPZinput/borg_1033+5051_multiband_modified_BPZinput.probs' '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121108BPZ/12110823279_BPZrun_borg_1033+5051_multiband_modified_BPZinput/html/probplots/' --verbose
#
#
#
# for i in 1211*; do plotBPZprobs.py $i/*.probs $i/html/probplots/ ; done 
#----------------------------
#   BUGS
#----------------------------
#
# Only runs ~768 objects before running out of memory (can't find leak)
# In that case restart using the --skip keyword. This memory leak should
# have been fixed by adding a fig.close() statement on 130227 by KBS
#    --> Fixed by KBS on 130702 by moving fig = plt.figure() outside loop
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-11-08  started by K. B. Schmidt (UCSB)
# 2013-01-06  added --skip. K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import commands   # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
import numpy as np# enable opening with genfromtxt
import matplotlib as plt
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
parser.add_argument("probfile", type=str, help="Output file from BPZ containing the probability distributions")
parser.add_argument("outputdir", type=str, help="Directory to save plots in")
# ---- optional arguments ----
parser.add_argument("--skip", type=int, help="Provide integer indicating no. of objects to skip, i.e. skip first N objects")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--eps", action="store_true", help="Turn plots into eps files")
parser.add_argument("--pdf", action="store_true", help="Turn plots into pdf files")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# reading probability distributions file
buffer    = open(args.probfile).readlines()  # read file into buffer
Nobj = len(buffer)-1
if args.verbose: print ':: '+sys.argv[0]+' :: Found '+str(Nobj)+' to create probplots for'
PNinput      = pathAname(args.probfile)
fieldnameTEX = PNinput[1].split('_')
fieldnameTEX = fieldnameTEX[0]+'\_'+fieldnameTEX[1]
#-------------------------------------------------------------------------------------------------------------
hdr   = buffer[0].split('(')
hdr   = hdr[2].split(')')
nos   = hdr[0].split(',')
xvals = np.arange(float(nos[0]),float(nos[1]),float(nos[2]))
Nbins = len(xvals)
#-------------------------------------------------------------------------------------------------------------
proball = [0.0]*Nbins
#-------------------------------------------------------------------------------------------------------------
# PLOTTING
import matplotlib.pyplot as plt 
from matplotlib  import cm
from astropy import cosmology 

plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=13)           # setting text font
plt.rc('xtick', labelsize=13) 
plt.rc('ytick', labelsize=13) 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
fig = plt.figure()      # create a figure object
for ii in xrange(Nobj): # looping over objects
    if args.skip and ii <= args.skip: continue
    line   = buffer[ii+1].split()
    yvals = line[1:Nbins+1]
    proball = np.add(proball,np.array(map(float, yvals)))  # adding all probabilities together
    plotname = args.outputdir+'probplot'+line[0]
    if args.verbose: print ':: '+sys.argv[0]+' :: Creating figure '+plotname
    fig.clf()           # clearing figure
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    plt.plot(xvals,yvals)
    ax.grid(True,linestyle='-',color='0.75')
    ax.set_xlabel('z')
    ax.set_ylabel('P$(z)$')
    ax.set_title(fieldnameTEX)
    fig.savefig(plotname+'.png')
    if args.eps: fig.savefig(plotname+'.eps')
    if args.pdf: fig.savefig(plotname+'.pdf')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if not args.skip: # only create probplot all if skip keyword is not given
    plotname = args.outputdir+'probplotALL'
    if args.verbose: print ':: '+sys.argv[0]+' :: Creating figure '+plotname
    fig = plt.figure()  # create a figure object
    fig.clf()                                        # clearing figure
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    plt.plot(xvals,proball)
    ax.grid(True,linestyle='-',color='0.75')
    ax.set_xlabel('z')
    ax.set_ylabel('P$(z)_\mathrm{total}$ ('+str(Nobj)+' objects)')
    ax.set_title(fieldnameTEX)
    fig.savefig(plotname+'.png')
    if args.eps: fig.savefig(plotname+'.eps')
    if args.pdf: fig.savefig(plotname+'.pdf')
    #show()  # draw plot on screen
#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

