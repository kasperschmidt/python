#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# estimateNfromLF.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Estimating the number of objects down to Lmin provided a set of 
# Luminoisity function parameters and a volume curve, V(L)
#----------------------------
#   COMMENTS
#----------------------------
# 
#----------------------------
#   INPUTS:
#----------------------------
# lfparam         : the 3 luminoisty function paramters alpha lstar phistar 
#                   in that order. Expect units [#] [1e44 erg/s] [Mpc-3]
# volumedict      : Dictionary containing the total volume curve in ['total_5sig']
#                   File can be crearted with mpd_pymc_estimateVolumes_BoRG13volumes.npz
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --range          : The range to integrate over. Accepts either 1=L, 2=mapp or 3=Mabs min
#                    and max range given as e.g. type min max = 2 24.5 30 
# --verbose        : set -verbose to get info/messages printed to the screen
# --show           : showning plot full volume curve and part integrated over
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# Directoryfile    : Output file containing the name of the directories the output has been
#                    put into. Directories will be moved to directory of catlist
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x runBPZmultiplecats.py       (only required once)
# bash> estimateNfromLF.py -1.98 0.76 4.3e-4 mpd_pymc_estimateVolumes_BoRG13andUDFERSvolumes_130817.npz --range 3 -22.5 -19.5 --verbose
# bash> estimateNfromLF.py -1.98 0.76 4.3e-4 mpd_pymc_estimateVolumes_BoRG13volumes_130817.npz --range 3 -22.5 -19.5 --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-01-28  started by K. B. Schmidt (UCSB)
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
import scipy.integrate
import pdb                 # for debugging with pdb.set_trace()
from cosmocalc import cosmocalc
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("lfparam", type=float, nargs=3, help="Luminosity function paramters")
parser.add_argument("volumedict", type=str, help="Dictionary containing volumes")
# ---- optional arguments ----
parser.add_argument("--range", type=float, nargs=3, help="The range to integrate over (kind min max). Kind: 1=L, 2=mapp, 3=Mabs")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--eps", action="store_true", help="Turn plots into eps files")
parser.add_argument("--png", action="store_true", help="Turn plots into png files")
parser.add_argument("--show", action="store_true", help="Showing plots on screen")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
# load file for plotting
voldict = np.load(args.volumedict) # loading dictionary with volumes created above
#-------------------------------------------------------------------------------------------------------------
vdict = voldict['total_5sig']
Ldict = voldict['Lval']
ent   = np.where(Ldict > 0) # all entries as default

if args.range and (args.range[0] == 1): # Lcut
    ent = np.where((Ldict > args.range[1]) & (Ldict < args.range[2]))
    Lmin, Lmax = args.range[1], args.range[2]
elif args.range and (args.range[0] == 2): # mapp cut
    Mabs   = kbs.L2Mabs(Ldict,MUVsun=5.48)
    mapp   = kbs.magabs2app(Mabs,'dummy','dummy','dummy',Av='dummy',band='Jbradley2012',cos='WMAP7BAOH0')
    ent    = np.where((mapp > args.range[1]) & (mapp < args.range[2]))
    Lmin   = kbs.Mabs2L(kbs.magapp2abs(args.range[2],'dummy','dummy','dummy',Av='dummy',band='Jbradley2012',cos='WMAP7BAOH0'),MUVsun=5.48)
    Lmax   = kbs.Mabs2L(kbs.magapp2abs(args.range[1],'dummy','dummy','dummy',Av='dummy',band='Jbradley2012',cos='WMAP7BAOH0'),MUVsun=5.48)
elif args.range and (args.range[0] == 3): # Mabs cut
    Mabs   = kbs.L2Mabs(Ldict,MUVsun=5.48)
    ent = np.where((Mabs > args.range[1]) & (Mabs < args.range[2]))
    Lmin   = kbs.Mabs2L(args.range[2],MUVsun=5.48)
    Lmax   = kbs.Mabs2L(args.range[1],MUVsun=5.48)
elif args.range:
    print ' ERROR wrong choice of range type. Chose between 1=L, 2=mapp and 3=Mabs --> ABORTING \n'
    pdb.set_trace()

vfct = vdict[ent]
Lval = Ldict[ent]
Lvalmin = np.min(Lval)
Lvalmax = np.max(Lval)

alpha   = args.lfparam[0]
Lstar   = args.lfparam[1]
phistar = args.lfparam[2]
#-------------------------------------------------------------------------------------------------------------
def schechter(L,alpha,Lstar,phistar):
    '''
    Calculating the Schechter function value. 
        Phi(L)dL = phi* x (L/L*) x exp(L/L*) d(L/L*)
    Units returned are:
        [phistar] = Mpc^-3
    Which corresponds to
    '''
    power   = (L/Lstar)**alpha
    expon   = np.exp(-(L/Lstar))
    schval  = phistar / Lstar * power * expon
    return schval
#-------------------------------------------------------------------------------------------------------------
def volume(L,Lval,volbins):
    '''
    Function linearly interpolating between volume bins to get a volume value
    '''
    if L < np.min(Lval):
        vol = 0.0
    elif L > np.max(Lval):
        sys.exit(':: ERROR :: Input Volume L in volume function is larger than max(Lvector) --> ABORTING') 
    else:
        newL = np.sort(np.append(Lval,L))
        newV = kbs.interpn(Lval,volbins,newL)
        ent  = np.where(newL == L)
        vol  = newV[ent]
    return vol
#-------------------------------------------------------------------------------------------------------------
test = 0
if test == 1:
    vfct  = np.zeros(len(Lval))+1e5 # areas all set to 1e5 for sanity check agains Micheles results.
    vdict = np.zeros(len(Ldict))+1e5
    #Michele's result: nsources[0.00043, -20.26, -1.98, -22.5, -19.5, 100000.0] = 28.3102

if args.show:
    import pylab as plt
    plt.plot(Ldict,vdict,'ko',ms=7)
    plt.plot(Lval,vfct,'co',ms=4)
    plt.xscale('log')
    plt.show()

Nobjrange = scipy.integrate.trapz(vfct*schechter(Lval,alpha,Lstar,phistar),Lval)

Nobjrange2, error = scipy.integrate.quad(lambda ll: volume(ll,Ldict,vdict)*schechter(ll,alpha,Lstar,phistar),Lvalmin,Lvalmax)
Nobj,       error = scipy.integrate.quad(lambda ll: volume(ll,Ldict,vdict)*schechter(ll,alpha,Lstar,phistar),Lmin,Lmax)


#print Nobjrange,Nobjrange2,Nobj
if args.verbose: print ' - Estimated number of objects in  L/[1e44erg/s] range ['+str("%.2e" % Lmin)+','+str("%.2e" % Lmax)+'] is '+str("%.2e" % Nobj)+' Mpc^-3'
    
#zobj    = 8
#Uvolume = cosmocalc(zobj,H0=70.,WM=0.3)['VCM_Gpc3']*1e9 # in volume of Universe in Mpc^3

#-------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
