#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# extractNoiseApertures_allBoRG_evaluate.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Plot and evaluate an output dictionary from extractNoiseApertures_allBoRG.py 
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# npzfile          : npz file containing the dictionary to evaluate
# raper            : radius of apertures in arcsec
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --sigmacut       : the N sigma to count values outside (DEFAULT is sigmacut = 5.0)
# --drawfile       : to prevent drawing samples every time a filename can be provided
#                    containing a dictionary of the draws and outliers. If keyword is
#                    not set a dictionary of this type is automatically created and written to
#                    'npzfile'_drawnsamplesANDoutliers.npz
#                    note if drawfile is dummy then the draws are ignored and only the input
#                    npzfile is evaluated.
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
# bash> extractNoiseApertures_allBoRG_evaluate.py extractNoiseApertures_allBoRG_output130723_2k_0p40.npz 0.40  --verbose
# bash> extractNoiseApertures_allBoRG_evaluate.py extractNoiseApertures_allBoRG_output13'723_2k_0p32.npz 0.32 --verbose

# bash> extractNoiseApertures_allBoRG_evaluate.py extractNoiseApertures_allBoRG_output130726_2k_0p40.npz 0.40 --verbose
# bash> extractNoiseApertures_allBoRG_evaluate.py extractNoiseApertures_allBoRG_output130726_2k_0p32.npz 0.32 --verbose --nodraw


# bash> extractNoiseApertures_allBoRG_evaluate.py extractNoiseApertures_allBoRG_output130726_2k_0p32.npz 0.32 --verbose --drawfile extractNoiseApertures_allBoRG_output130726_2k_0p32_drawnsamplesANDoutliers_4Ddraw130829.npz
# --- extractNoiseApertures_allBoRG_output130726_2k_0p32_drawnsamplesANDoutliers.npz
# --- extractNoiseApertures_allBoRG_output130726_2k_0p32_drawnsamplesANDoutliers_test130829.npz

# bash> extractNoiseApertures_allBoRG_evaluate.py extractNoiseApertures_allBoRG_output130726_2k_0p32.npz 0.32 --drawfile extractNoiseApertures_allBoRG_output130726_2k_0p32_drawnsamplesANDoutliers_130816allfields.npz --verbose

# --- Apertures on a grid ---
# bash> extractNoiseApertures_allBoRG_evaluate.py extractNoiseApertures_allBoRG_output130830_ongrid_0p32.npz 0.32 --drawfile dummy --verbose

# bash> extractNoiseApertures_allBoRG_evaluate.py extractNoiseApertures_allBoRG_output130918_ongrid_0p32.npz 0.32 --drawfile dummy --verbose

#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-07-19  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import getopt              # used to extract/obtain the optional input
import numpy as np         # enable opening with genfromtxt
import pdb                 # for debugging with pdb.set_trace()
import commands
import glob
from scipy.optimize import curve_fit
import drawsamplefromdistribution as dfd
import drawsamplefromdistribution_multiD as dfdMD
from time import localtime, strftime
#-------------------------------------------------------------------------------------------------------------
def gauss(x, *p):
    """
    Gaussian/normal distribution
    """
    A, mu, sigma = p
    t = (x-mu)/sigma
    return A*np.exp(-t**2/2.)
#-------------------------------------------------------------------------------------------------------------
def getaperval(fields,dict,datacol,Nsig=5.0,verbose=1):
    """
    Getting the aperture values for all fields, std and the number of N sigma outliers

    input:
       fields    list of fields
       dict      dictionary containing data
       datacol   list of columns containing aperture values to use for each field - starting from 0
       Nsig      the nummber of sigmas to count objects outside (default is 5 sigma)
    """
    Nfield        = len(fields)
    aperval       = np.array([])
    apervalnonorm = np.array([])
    stddev        = np.array([])
    NsigavTOT     = np.array([])

    Napertot      = 0
    for ii in xrange(Nfields):
        av      = dict[fields[ii]][:,datacol[ii]]
        sd      = np.std(av)                  # the standard deviation of the aperture values
        mean    = np.mean(av)
    
        Nav     = len(av)
        Nsigav  = len(np.where( abs(av-mean)/sd > Nsig)[0])
        #if verbose == 1: print ' - '+str(Nsig)+' sigma outliers in ',fields[ii],' out of ',Nav,' apertures : ',Nsigav
        
        NsigavTOT     = np.append(NsigavTOT,Nsigav)
        stddev        = np.append(stddev, sd )
        aperval       = np.append(aperval, (av-mean)/sd )  # combining aperture values normalized by standard dev (1 sigma)
        apervalnonorm = np.append(apervalnonorm, av )
        apervalnonorm_1sig = apervalnonorm.copy()
        apervalnonorm_1sig[np.where(np.abs(apervalnonorm) < sd)] = sd # setting fluxes < 1sigma to 1sigma value

        Napertot = Napertot + Nav # counting the aperturealues
        if verbose == 1: 
            print ' - ',fields[ii],' contains ',Nav,' apertures (total apertures: ',Napertot,')'
            if fields[ii] == 'borg_2203+1851xx':
                print 'Aper no. 1:',dict[fields[ii]][29425,:] 
                print 'Aper no. 2:',dict[fields[ii]][29429,:] 
                print 'Aper no. 3:',dict[fields[ii]][29430,:] 
            if fields[ii] == 'borg_1459+7146xx':
                print 'Aper no. 4:',dict[fields[ii]][25834,:] 
                
    if verbose == 1: print '\n - in Total that is  : ',np.sum(NsigavTOT)
    Naper   = len(aperval)
        
    return stddev, aperval, NsigavTOT, Naper, apervalnonorm, apervalnonorm_1sig
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("npzfile", type=str, help="Dictionary to evaluate (npz file from extractNoiseApertures_allBoRG_evaluate.py)")
parser.add_argument("raper", type=float, help="Radius of aperture in arcsec")
# ---- optional arguments ----
parser.add_argument("--sigmacut", type=float, help="Sigma value to count apertures outside")
parser.add_argument("--drawfile", type=str, help="File with dictionary containing draws and outlier counts")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--eps", action="store_true", help="Turn plots into eps files")
parser.add_argument("--png", action="store_true", help="Turn plots into png files")
parser.add_argument("--show", action="store_true", help="Showing plots on screen")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
# load dictionary
dict = np.load(args.npzfile)
keys    = dict.keys()
keylen  = [len(keys[ii]) for ii in range(len(keys))]
fields  = np.asarray(keys)[np.where(np.asarray(keylen) == 14)]
Nfields = len(fields)
if args.verbose: print ' - Found '+str(Nfields)+' fields in dictionary '+args.npzfile

sigcut  = 5.0 # default sigma value
if args.sigmacut: sigcut = args.sigmacut

Vcol = []
Ycol = []
Jcol = []
Hcol = []

for hh in xrange(Nfields):
    hdr = np.append(dict[fields[hh]+'_hdr'],np.asarray('10'))[0]
    colnames = hdr.split('#')[-1].split(' ')[1:]

    vent = [colnames[ii].find('f606w') for ii in range(len(colnames))]
    if len(np.where(np.asarray(vent) != -1)[0]) == 0: vent = [colnames[ii].find('f600lp') for ii in range(len(colnames))]
        
    Vcol.append(np.where(np.asarray(vent) != -1)[0][0])
    Ycol.append(np.where(np.asarray([colnames[ii].find('f098m') for ii in range(len(colnames))]) != -1)[0][0])
    Jcol.append(np.where(np.asarray([colnames[ii].find('f125w') for ii in range(len(colnames))]) != -1)[0][0])
    Hcol.append(np.where(np.asarray([colnames[ii].find('f160w') for ii in range(len(colnames))]) != -1)[0][0])
    
stddev_V, aperval_V, N5avTOT_V, Naper_V, apervalnonorm_V, apervalnonorm_1sig_V = getaperval(fields,dict,Vcol,Nsig=sigcut,verbose=0)
stddev_Y, aperval_Y, N5avTOT_Y, Naper_Y, apervalnonorm_Y, apervalnonorm_1sig_Y = getaperval(fields,dict,Ycol,Nsig=sigcut,verbose=0)
stddev_J, aperval_J, N5avTOT_J, Naper_J, apervalnonorm_J, apervalnonorm_1sig_J = getaperval(fields,dict,Jcol,Nsig=sigcut,verbose=1)
stddev_H, aperval_H, N5avTOT_H, Naper_H, apervalnonorm_H, apervalnonorm_1sig_H = getaperval(fields,dict,Hcol,Nsig=sigcut,verbose=0)

dict.close()
#-------------------------------------------------------------------------------------------------------------
# fitting gaussian to J band distribution
Nbins     = 200 # Note that cutting in 4D needs larger bins
histparam = np.histogram(aperval_J,bins=Nbins)
binsize   = (np.max(aperval_J)-np.min(aperval_J))/Nbins
counts    = histparam[0]
binedge   = histparam[1]
bincen    = binedge[0:-1]+binsize/2.
aperstd   = np.std(aperval_J)

p0                 = [np.max(counts), 0, aperstd]                       # guess for the fitting coefficients (A, mu and sigma)
gaussguess         = gauss(bincen, *p0)
cgauss, var_matrix = curve_fit(gauss, bincen, counts, p0=p0)  # fitting gauss to data
gaussfit           = gauss(bincen, *cgauss)                   # Get the fitted curve
#-------------------------------------------------------------------------------------------------------------
# counting outliers...
Ncut       = len(np.where(abs(aperval_J) > sigcut)[0])
NcutsigL   = len(np.where( aperval_J < cgauss[1]-5*cgauss[2] )[0])
NcutsigH   = len(np.where( aperval_J > cgauss[1]+5*cgauss[2] )[0])
Ncutsig    = NcutsigL+NcutsigH
fraccut    = (0.0+Ncut)/(0.0+Naper_J)
fraccutsig = (0.0+Ncutsig)/(0.0+Naper_J)

if args.verbose: 
    print '\n - Of the '+str(Naper_J)+' apertures in the dictionary...'
    print '   Apertures with |counts| > 5      :   '+str(Ncut)+' ~ '+str("%.2e" % fraccut)
    print '   Apertures with |counts| > 5sigma :   '+str(Ncutsig)+' ~ '+str("%.2e" % fraccutsig)
#-------------------------------------------------------------------------------------------------------------
# The outliers
#outent = np.where((abs(aperval_J) > 5) & (abs(aperval_H) > 2.5) & (abs(aperval_V) < 1.5))[0] # entries of outliers
outent = np.where((abs(aperval_J) > 8) & (abs(aperval_H) > 2.5) & (abs(aperval_V) < 1.5))[0] # entries of outliers
Nout   = len(outent)
Vout, Yout, Jout, Hout = aperval_V[outent], aperval_Y[outent], aperval_J[outent], aperval_H[outent]
if args.verbose: print ' - Requiring Jsig > 5 and Hsig > 2.5 and Vsig < 1.5 results in '+str(Nout)+' apertures\n'
#-------------------------------------------------------------------------------------------------------------
# Getting colors
# conversiong from cps [e/s] to Flampda [erg/cm2/A/s] from fits headers
PHOTFLAM_098 = 6.05015670000000E-20 #[erg/cm2/A/e]
PHOTFLAM_125 =        2.2483504E-20 #[erg/cm2/A/e]
PHOTFLAM_160 =        1.9275707E-20 #[erg/cm2/A/e]
FlamY, FlamJ, FlamH = apervalnonorm_Y*PHOTFLAM_098, apervalnonorm_J*PHOTFLAM_125, apervalnonorm_H*PHOTFLAM_160
#FlamY, FlamJ, FlamH = apervalnonorm_1sig_Y*PHOTFLAM_098, apervalnonorm_J*PHOTFLAM_125, apervalnonorm_H*PHOTFLAM_160

PHOTFNU_098 =        1.9636467E-07 # [Jy*sec/electron]
PHOTFNU_125 =        1.1692158E-07 # [Jy*sec/electron]
PHOTFNU_160 =        1.5187624E-07 # [Jy*sec/electron]
FnuY, FnuJ, FnuH = apervalnonorm_Y*PHOTFNU_098, apervalnonorm_J*PHOTFNU_125, apervalnonorm_H*PHOTFNU_160

ZP098, ZP125, ZP160 = 25.68, 26.25, 25.96

YJcol = -2.5 * np.log10( FlamY / FlamJ ) + ZP098 - ZP125
JHcol = -2.5 * np.log10( FlamJ / FlamH ) + ZP125 - ZP160

#YJcolnu = -2.5 * np.log10( FnuY / FnuJ ) + ZP098 - ZP125
#JHcolnu = -2.5 * np.log10( FnuJ / FnuH ) + ZP125 - ZP160

#-------------------------------------------------------------------------------------------------------------
# The outliers -- all cuts w/o JH
outentwoJH = np.where((abs(aperval_J) > 5) & 
                     (abs(aperval_H) > 2.5) & 
                     (abs(aperval_V) < 1.5) &
                     (YJcol > 1.75))[0]
NoutwoJH   = len(outentwoJH)
VoutwoJH, YoutwoJH, JoutwoJH, HoutwoJH = aperval_V[outentwoJH], aperval_Y[outentwoJH], aperval_J[outentwoJH], aperval_H[outentwoJH]
if args.verbose: print ' - Apertures satisfying the S/N_J, S/N_H, S/N_V, YJcol: '+str(NoutwoJH)+' apertures\n'
#-------------------------------------------------------------------------------------------------------------
# The outliers -- all cuts
outentall = np.where((abs(aperval_J) > 5) & 
                     (abs(aperval_H) > 2.5) & 
                     (abs(aperval_V) < 1.5) &
                     (YJcol > 1.75) &
                     ((JHcol)<0.02+0.15*(YJcol-1.75)))[0]
Noutall   = len(outentall)
Voutall, Youtall, Joutall, Houtall = aperval_V[outentall], aperval_Y[outentall], aperval_J[outentall], aperval_H[outentall]
if args.verbose: print ' - Apertures satisfying all cuts (S/N_J, S/N_H, S/N_V, YJcol, JHcol+YJcol): '+str(Noutall)+' apertures\n'
#-------------------------------------------------------------------------------------------------------------
# Cut aperture data on V and H 
VHcut   = np.where((aperval_V < 1.5) & (aperval_H > 2.5))[0]
VvalVHcut, YvalVHcut, JvalVHcut, HvalVHcut    = aperval_V[VHcut], aperval_Y[VHcut], aperval_J[VHcut], aperval_H[VHcut]
#-------------------------------------------------------------------------------------------------------------
# Cut aprture data on V
Vnondet = np.where(aperval_V < 1.5)[0] # entries of apertures with S/N_V < 1.5
Vval, Yval, Jval, Hval    = aperval_V[Vnondet], aperval_Y[Vnondet], aperval_J[Vnondet], aperval_H[Vnondet]
#-------------------------------------------------------------------------------------------------------------
oneANDtwo = 0
if oneANDtwo == 1:
    #-------------------------------------------------------------------------------------------------------------
    # Apertures prior to J and H S/N cuts
    VYJ   = np.where((aperval_V < 1.5) & (YJcol > 1.75))[0]
    #-------------------------------------------------------------------------------------------------------------

    VHcut   = np.where((aperval_V < 1.5) & (aperval_H > 2.5))[0]
    VvalVHcut, YvalVHcut, JvalVHcut, HvalVHcut    = aperval_V[VHcut], aperval_Y[VHcut], aperval_J[VHcut], aperval_H[VHcut]

    Vnondet = np.where(aperval_V < 1.5)[0] # entries of apertures with S/N_V < 1.5
    Vval, Yval, Jval, Hval    = aperval_V[Vnondet], aperval_Y[Vnondet], aperval_J[Vnondet], aperval_H[Vnondet]

    if args.verbose: print ' - Drawing from cumulative distribution (after V and H S/N cuts)'

    Awfc3         = 123.0*136.0 # http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c03_optimum_instr4.html#345239
    Aaper         = np.pi * args.raper**2.
    NdrawsVHcut   = int(round(Awfc3/Aaper))           # The number of apertures in a wfc3 field.
    drawvalVHcut  = np.zeros((NdrawsVHcut,2,Nfields)) # Array to contain values of draws
    NoutVHcutdraw = np.zeros(Nfields)                 # Arrau to contain coutn of outliers

    #for ff in xrange(Nfields):
    for ff in xrange(1):
        if args.verbose: print '   Pseudo -',fields[ff],' on '+strftime("%a, %d %b %Y %H:%M:%S", localtime())
        #draws = dfd.drawsamplefromdistribution(JvalVHcut,Nbins,points2D=None,bins2D=None,verbose=True,size=NdrawsVHcut)
        draws = dfdMD.drawsample(JvalVHcut,Nbins,verbose=True,size=NdrawsVHcut)
        drawvalVHcut[:,:,ff] = draws
        NoutVHcutdraw[ff] = len(np.where(np.abs(drawvalVHcut[:,0,ff]) > 5.0)[0])  # counting outliers
    #-------------------------------------------------------------------------------------------------------------
    if args.verbose: print ' - Drawing from 2D distribution of J and H (after V S/N cut)'
    NdrawsJH      = NdrawsVHcut
    drawJHval     = np.zeros((NdrawsJH,4,Nfields)) # Array to contain values of draws
    NoutdrawJH    = np.zeros(Nfields)              # Array to contain coutn of outliers
    distJH        = np.asarray([Jval,Hval])        # array with distribution to draw from
    
    #for ff in xrange(Nfields):
    for ff in xrange(1):
        if args.verbose: print '   Pseudo -',fields[ff],' on '+strftime("%a, %d %b %Y %H:%M:%S", localtime())
        #drawval2D  = dfd.drawsamplefromdistribution(Jval,Nbins,points2D=Hval,bins2D=Nbins,verbose=True,size=NdrawsJH,plot=False)
        drawval2D = dfdMD.drawsample(distJH,Nbins,verbose=True,size=NdrawsJH,Ncut=3)
        dfdMD.plot2D(distJH[0,:],distJH[1,:],drawval2D[:,0:2],drawval2D[:,2:4],
                   Nbins,save='drawsamplefromdistributionMD_draw2D.pdf')
        drawJHval[:,:,ff] = drawval2D
        NoutdrawJH[ff] = len(np.where((np.abs(drawJHval[:,0,ff]) > 5.0) & (np.abs(drawJHval[:,2,ff]) > 2.5))[0]) # count outliers    
    pdb.set_trace()

#-------------------------------------------------------------------------------------------------------------
if (not args.drawfile):
    if args.verbose: print ' - Drawing from 4D distribution prior to any cuts'
    Awfc3         = 123.0*136.0 # http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c03_optimum_instr4.html#345239
    Aaper         = np.pi * args.raper**2.
    Ndraws4D      = int(round(Awfc3/Aaper))        # The number of apertures in a wfc3 field.
    draw4Dval     = np.zeros((Ndraws4D,8,Nfields)) # Array to contain values of draws
    Noutdraw4D    = np.zeros(Nfields)              # Array to contain coutn of outliers    
    dist4D        = np.asarray([aperval_V,aperval_Y,aperval_J,aperval_H]) # array with distribution to draw from

    #for ff in xrange(Nfields):
    for ff in xrange(1):
        if args.verbose: print '   Pseudo -',fields[ff],' on '+strftime("%a, %d %b %Y %H:%M:%S", localtime())
        draws4D = dfdMD.drawsample(dist4D,Nbins,verbose=False,size=Ndraws4D)

        V4Ddraw = draws4D[:,0]
        Y4Ddraw = draws4D[:,2]
        J4Ddraw = draws4D[:,4]
        H4Ddraw = draws4D[:,6]

        FlamY   = Y4Ddraw*np.std(apervalnonorm_Y) * PHOTFLAM_098 # multiplying with 1sigma to get e/s and then converting
        FlamJ   = J4Ddraw*np.std(apervalnonorm_J) * PHOTFLAM_125 # multiplying with 1sigma to get e/s and then converting
        FlamH   = H4Ddraw*np.std(apervalnonorm_H) * PHOTFLAM_160 # multiplying with 1sigma to get e/s and then converting
        YJcol4D = -2.5 * np.log10( FlamY / FlamJ ) + ZP098 - ZP125
        JHcol4D = -2.5 * np.log10( FlamJ / FlamH ) + ZP125 - ZP160
        
        Noutdraw4D[ff] = len(np.where((np.abs(J4Ddraw) > 5.0) & 
                                      (np.abs(H4Ddraw) > 2.5) & 
                                      (np.abs(V4Ddraw) < 1.5) &
                                      (np.abs(YJcol4D) > 1.75))[0]) # count outliers

        draw4Dval[:,:,ff] = draws4D # filling output array
    #-------------------------------------------------------------------------------------------------------------
    # Write draws and outlier counts to file
    filename=args.npzfile.replace('.npz','_drawnsamplesANDoutliers.npz')
    data = {}
    #data['NdrawsVHcut']   = NdrawsVHcut
    #data['drawvalVHcut']  = drawvalVHcut
    #data['NoutVHcutdraw'] = NoutVHcutdraw

    #data['NdrawsJH']      = NdrawsJH
    #data['drawJHval']     = drawJHval
    #data['NoutdrawJH']    = NoutdrawJH

    data['Ndraws4D']      = Ndraws4D
    data['draw4Dval']     = draw4Dval
    data['Noutdraw4D']    = Noutdraw4D

    np.savez(filename,**data) # save dictionary as binary file
    if args.verbose: print ' - Saved drawn samples and outlier counts to ',filename

    if args.verbose: print ' --- Done drawing... stoppping. '
    pdb.set_trace()
elif args.drawfile == 'dummy':
    # creating dummy arrays
    Ndraws4D      = 100
    draw4Dval     = np.zeros((Ndraws4D,8,Nfields)) # Array to contain values of draws
    Noutdraw4D    = np.zeros(Nfields)              # Array to contain coutn of outliers    
    dist4D        = np.asarray([aperval_V,aperval_Y,aperval_J,aperval_H]) # array with distribution to draw from
else:
    # Read draw file containing dictionary
    data = np.load(args.drawfile)
    #NdrawsVHcut   = data['NdrawsVHcut']
    #drawvalVHcut  = data['drawvalVHcut']
    #NoutVHcutdraw = data['NoutVHcutdraw']

    #NdrawsJH      = data['NdrawsJH']
    #drawJHval     = data['drawJHval']
    #NoutdrawJH    = data['NoutdrawJH']

    Ndraws4D      = data['Ndraws4D']
    draw4Dval     = data['draw4Dval']
    Noutdraw4D    = data['Noutdraw4D']
    dist4D        = np.asarray([aperval_V,aperval_Y,aperval_J,aperval_H])
#-------------------------------------------------------------------------------------------------------------
Noutdraw4Dtot = np.sum(Noutdraw4D)
# Getting colors for draws

YJcol4Dall, JHcol4Dall = np.zeros((Nfields,Ndraws4D)), np.zeros((Nfields,Ndraws4D))
    
for ff in xrange(Nfields):
    V4Ddraw = draw4Dval[:,0,ff]
    Y4Ddraw = draw4Dval[:,2,ff]
    J4Ddraw = draw4Dval[:,4,ff]
    H4Ddraw = draw4Dval[:,6,ff]

    FlamY   = Y4Ddraw*np.std(apervalnonorm_Y) * PHOTFLAM_098 # multiplying with 1sigma to get e/s and then converting
    FlamJ   = J4Ddraw*np.std(apervalnonorm_J) * PHOTFLAM_125 # multiplying with 1sigma to get e/s and then converting
    FlamH   = H4Ddraw*np.std(apervalnonorm_H) * PHOTFLAM_160 # multiplying with 1sigma to get e/s and then converting
    YJcol4Dall[ff,:] = -2.5 * np.log10( FlamY / FlamJ ) + ZP098 - ZP125
    JHcol4Dall[ff,:] = -2.5 * np.log10( FlamJ / FlamH ) + ZP125 - ZP160    
#-------------------------------------------------------------------------------------------------------------
#                                       PLOTTING
#-------------------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt

Fsize = 18
plt.rc('text', usetex=True)                      # enabling LaTex rendering of text
plt.rc('font', family='serif',size=Fsize)           # setting text font
plt.rc('xtick', labelsize=Fsize) 
plt.rc('ytick', labelsize=Fsize) 

# creating plot of k-logL* plane with samples
plotname1 = args.npzfile.replace('.npz','_aperturehist.pdf')
if args.eps: plotname1 = plotname1.replace('.pdf','.eps')
if args.png: plotname1 = plotname1.replace('.pdf','.png')
if args.verbose: print '\n - Plotting histogram of apertures to',plotname1
plt.clf()

#plt.hist(aperval_V,bins=Nbins,color="b",histtype="step",lw=1)
#plt.hist(aperval_Y,bins=Nbins,color="m",histtype="step",lw=1)
#plt.hist(aperval_H,bins=Nbins,color="c",histtype="step",lw=1)

plt.hist(aperval_J,bins=Nbins,color="k",histtype="step",lw=3)

#plt.plot(bincen,gaussguess,'b:',label=r'Gaussian guess [A, $\mu$, $\sigma$] = ['+str("%.2f" % p0[0])+','+str("%.2f" % p0[1])+','+str("%.2f" % p0[2])+']')
plt.plot(bincen,gaussfit,'k--',label=r'Gaussian fit to J band dist.',lw=2)
#   [A, $\mu$, $\sigma$] = ['+str("%.2f" % cgauss[0])+','+str("%.2f" % cgauss[1])+','+str("%.2f" % cgauss[2])+']'

#plt.plot(bincen,counts,'go')

# filling sigma regions:
xmin = cgauss[1]-7*cgauss[2]
xmax = cgauss[1]+7*cgauss[2]
ymin = 20.0 #0.0
ymax = np.max(counts)+0.2*np.max(counts)
plt.fill_between([cgauss[1]-1*cgauss[2],cgauss[1]+1*cgauss[2]],[ymin,ymin],[ymax,ymax],alpha=0.1,color='k')
plt.fill_between([cgauss[1]-2*cgauss[2],cgauss[1]+2*cgauss[2]],[ymin,ymin],[ymax,ymax],alpha=0.1,color='k')
plt.fill_between([cgauss[1]-3*cgauss[2],cgauss[1]+3*cgauss[2]],[ymin,ymin],[ymax,ymax],alpha=0.1,color='k')
plt.fill_between([cgauss[1]-4*cgauss[2],cgauss[1]+4*cgauss[2]],[ymin,ymin],[ymax,ymax],alpha=0.1,color='k')
plt.fill_between([cgauss[1]-5*cgauss[2],cgauss[1]+5*cgauss[2]],[ymin,ymin],[ymax,ymax],alpha=0.1,color='k')
plt.hist(np.linspace(-1010,-1000,10),bins=10,color='k',alpha=0.50,label=r'1$\sigma$-5$\sigma$ of Gaussian fit')

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
#plt.xlabel(r'(count - mean(count)) / std(count) $\sim$ $\sigma$ ')
plt.xlabel(r'S/N in J')
plt.ylabel('Count')
plt.yscale('log')

leg = plt.legend(fancybox=True, loc='lower center',numpoints=1,prop={'size':12})
leg.get_frame().set_alpha(0.6)

plt.savefig(plotname1)
if args.show: plot.show()  # draw plot on screen

#-------------------------------------------------------------------------------------------------------------
from matplotlib.colors import LogNorm

plotname = args.npzfile.replace('.npz','_2DhistJH.pdf')
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting 2D histogram of apertures to',plotname
plt.clf()

heatmap, xedges, yedges = np.histogram2d(Hval,Jval,bins=(Nbins,Nbins))
extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]

plt.imshow(heatmap, extent=extent, interpolation='nearest',norm=LogNorm(vmin=1, vmax=1000))
#plt.plot(aperval_J,aperval_H,'r.')

xmin, xmax, ymin, ymax = -10, 10, -10, 10

#plt.plot([-5,-5],[ymin,-2.5],'k--',label=r'$\sigma$ cuts')
#plt.plot([-5,-5],[ymin,-2.5],'k--',label=r'S/N cuts')
#plt.plot([-5,-5],[2.5,ymax],'k--')

#plt.plot([5,5],[ymin,-2.5],'k--')
plt.plot([5,5],[2.5,ymax],'k-',lw=3)
plt.plot([8,8],[2.5,ymax],'k-',lw=3)

#plt.plot([xmin,-5],[2.5,2.5],'k--')
plt.plot([5,xmax],[2.5,2.5],'k-',lw=3)

#plt.plot([xmin,-5],[-2.5,-2.5],'k--')
#plt.plot([5,xmax],[-2.5,-2.5],'k--')

#plt.plot(Joutall,Houtall,'ro',label=r'S/N cut outliers: N='+str(Noutall))
#plt.plot(Joutall,Houtall,'ro',label=r'S/N + color cuts: N='+str(Noutall))


plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
#plt.xlabel(r'(count - mean(count)) / std(count) $\sim$ $\sigma$ in J')
plt.xlabel(r'S/N in J aperture')
#plt.ylabel(r'(count - mean(count)) / std(count) $\sim$ $\sigma$ in H')
plt.ylabel(r'S/N in H aperture')

cbar = plt.colorbar()
cbar.set_label(r'2D Histogram of apertures with S/N in V $<$ 1.5')

#leg = plt.legend(fancybox=True, loc='lower center',numpoints=1,prop={'size':12})
#leg.get_frame().set_alpha(0.6)

plt.savefig(plotname)
if args.show: plot.show()  # draw plot on screen
#-------------------------------------------------------------------------------------------------------------
plotname = args.npzfile.replace('.npz','_outlierYJHcol.pdf')
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting histogram of outliers after field-draws to',plotname
plt.clf()

xmin, xmax, ymin, ymax = -1.5, 1, -2, 3.5

plt.plot(JHcol[outent],YJcol[outent],'ko',label=r'outliers: N='+str(Nout))

for ff in xrange(Nfields):
    outent = np.where((abs(draw4Dval[:,4,ff]) > 5.0) & 
                      (abs(draw4Dval[:,6,ff]) > 2.5) & 
                      (abs(draw4Dval[:,0,ff]) < 1.5))[0] # entries of outliers

    if ff == 1 and args.drawfile != 'dummy':
        plt.plot(JHcol4Dall[ff,outent],YJcol4Dall[ff,outent],'r.',label=r'drawn outliers: N='+str(Noutdraw4Dtot))
    elif args.drawfile == 'dummy':
        dum = 'my'
    else:
        plt.plot(JHcol4Dall[ff,outent],YJcol4Dall[ff,outent],'r.',)

plt.plot([xmin,0.02],[1.75,1.75],'k--',label=r'color selection')
plt.plot([0.02,0.02+0.15*(ymax-1.75)],[1.75,ymax],'k--')

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xlabel(r'$J_{125} - H_{160}$')
plt.ylabel(r'$Y_{098} - J_{125}$')

leg = plt.legend(fancybox=True, loc='upper left',numpoints=1,prop={'size':12})
leg.get_frame().set_alpha(0.6)

plt.savefig(plotname)
if args.show: plot.show()  # draw plot on screen       
#-------------------------------------------------------------------------------------------------------------
plotname = args.npzfile.replace('.npz','_histdrawoutliers.pdf')
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting histogram of outliers after field-draws to',plotname
plt.clf()

hist = plt.hist(Noutdraw4D,bins=10,color="k",histtype="step",lw=3)

#leg = plt.legend(fancybox=True, loc='upper left',numpoints=1,prop={'size':12})
#leg.get_frame().set_alpha(0.6)

plt.xlabel(r'Expected outliers in '+str(Ndraws4D)+' random empty apertures per field')
plt.ylabel(r'Count')

plt.xlim(0,np.max(hist[1])+1)
plt.ylim(0,np.max(hist[0])+5)

plt.savefig(plotname)
if args.show: plot.show()  # draw plot on screen   

#-------------------------------------------------------------------------------------------------------------
plotname = args.npzfile.replace('.npz','_histdrawoutliers_norm2k.pdf')
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting histogram of outliers after field-draws to',plotname
plt.clf()

hist = plt.hist(Noutdraw4D*2000.0/Ndraws4D,bins=10,color="k",histtype="step",lw=3)

#leg = plt.legend(fancybox=True, loc='upper left',numpoints=1,prop={'size':12})
#leg.get_frame().set_alpha(0.6)

plt.xlim(0,np.max(hist[1])+1)
plt.ylim(0,np.max(hist[0])+5)

plt.xlabel(r'Expected outliers in 2000 random empty apertures per field')
plt.ylabel(r'Count')

plt.savefig(plotname)
if args.show: plot.show()  # draw plot on screen   
#-------------------------------------------------------------------------------------------------------------
#What draw to plot distributions for
ff = 0
#-------------------------------------------------------------------------------------------------------------
plotname = args.npzfile.replace('.npz','_2DhistDraw_VY.pdf')
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting 2D hist of YJ draws to ',plotname
dfdMD.plot2D(dist4D[0,:],dist4D[1,:],draw4Dval[:,0:2,ff],draw4Dval[:,2:4,ff],Nbins,xlab='S/N V',ylab='S/N Y',save=plotname)
#-------------------------------------------------------------------------------------------------------------
plotname = args.npzfile.replace('.npz','_2DhistDraw_YJ.pdf')
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting 2D hist of YJ draws to ',plotname
dfdMD.plot2D(dist4D[1,:],dist4D[2,:],draw4Dval[:,2:4,ff],draw4Dval[:,4:6,ff],Nbins,xlab='S/N Y',ylab='S/N J',save=plotname)
#-------------------------------------------------------------------------------------------------------------
plotname = args.npzfile.replace('.npz','_2DhistDraw_JH.pdf')
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting 2D hist of YJ draws to ',plotname
dfdMD.plot2D(dist4D[2,:],dist4D[3,:],draw4Dval[:,4:6,ff],draw4Dval[:,6:8,ff],Nbins,xlab='S/N J',ylab='S/N H',save=plotname)
#-------------------------------------------------------------------------------------------------------------
plotname = args.npzfile.replace('.npz','_1DhistDraw_V.pdf')
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting 1D display of V draws to ',plotname
dfdMD.plot1D(dist4D[0,:],draw4Dval[:,0:2,ff],Nbins,xlab='S/N V',save=plotname)
#-------------------------------------------------------------------------------------------------------------
plotname = args.npzfile.replace('.npz','_1DhistDraw_Y.pdf')
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting 1D display of Y draws to ',plotname
dfdMD.plot1D(dist4D[1,:],draw4Dval[:,2:4,ff],Nbins,xlab='S/N Y',save=plotname)
#-------------------------------------------------------------------------------------------------------------
plotname = args.npzfile.replace('.npz','_1DhistDraw_J.pdf')
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting 1D display of J draws to ',plotname
dfdMD.plot1D(dist4D[2,:],draw4Dval[:,4:6,ff],Nbins,xlab='S/N J',save=plotname)
#-------------------------------------------------------------------------------------------------------------
plotname = args.npzfile.replace('.npz','_1DhistDraw_H.pdf')
if args.eps: plotname = plotname.replace('.pdf','.eps')
if args.png: plotname = plotname.replace('.pdf','.png')
if args.verbose: print ' - Plotting 1D display of H draws to ',plotname
dfdMD.plot1D(dist4D[3,:],draw4Dval[:,6:8,ff],Nbins,xlab='S/N H',save=plotname)    
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------

