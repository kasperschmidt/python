"""
#----------------------------
#   NAME
#----------------------------
# kbsutilities.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Collection of subroutines and scripts convenient for coding
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# >>> import kbsutilities as kbs
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-01-04  started by K. B. Schmidt (UCSB)
#----------------------------
"""
#-------------------------------------------------------------------------------------------------------------
# importing modules
import kbsutilities as kbs
import fileinput
import sys
import math
import datetime
import numpy as np
import types
import pdb
import MiGs
import fits2ascii as f2a
import commands
import time
import pyfits
import os
import glob
#from wand.image import Image as wImage
#from wand.color import Color as wColor
from PyPDF2 import PdfFileReader, PdfFileMerger
import collections
from reproject import reproject_interp

### MATPLOTLIB ###
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle

### SCIPY ###
from scipy import odr
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter

### ASTROPYSICS ###
import astropysics
#import astropysics.obstools     #KBS170728 produce : TypeError: 'float' object cannot be interpreted as an index
#from astropysics import coords  #KBS170728 produce : TypeError: 'float' object cannot be interpreted as an index
from astropysics.constants import choose_cosmology

### ASTROPY ###
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.table as atab
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.io.fits as afits
from astropy import wcs

### IMAGE REGISTRATION ###
from image_registration import chi2_shift
from image_registration.fft_tools import shift


#-------------------------------------------------------------------------------------------------------------
def test4float(str): 
    """testing whehter a string contains letters (i.e. whether it can be converted to a float)"""
    try:
        float(str)
        return True
    except ValueError:
        return False

#-------------------------------------------------------------------------------------------------------------
def replaceAll(file,searchExp,replaceExp): 
    """search and replace in file"""
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

#-------------------------------------------------------------------------------------------------------------
def arr2str(arr): 
    """Taking line of np.array or np.void and turning into string for easy passing to file"""
    strout = str(arr).replace(",",' ').replace("'",' ').replace('[',' ').replace(']',' ').replace('(',' ').replace(')',' ')+'  \n'
    return strout

#-------------------------------------------------------------------------------------------------------------
def pathAname(str): 
    """splitting string with path and name in two"""
    strsplit = str.split('/')               # splitting string
    name     = strsplit[-1]                 # saving filename (last entry of catsplit)
    slash    = '/'
    path=slash.join(strsplit[0:-1])+'/'     # putting path back together
    return [path,name]

#-------------------------------------------------------------------------------------------------------------
def DandTstr(): 
    """Creating a string with date and time on the format yymmddhhmmss"""

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

#-------------------------------------------------------------------------------------------------------------
def DandTstr2():
    """Creating a string with date and time on the format %Y-%m-%d %H:%M:%S"""
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#-------------------------------------------------------------------------------------------------------------
def divideWith0(numerator,denominator,badval): 
    """Replacing division with 0 in ratios with value instead of crashing"""
    return np.array([numeratorij/denominatorij if denominatorij != 0.0 else badval for numeratorij, denominatorij in zip(numerator,denominator)])

#-------------------------------------------------------------------------------------------------------------
def fluxlimit(mABobj,SNobj,SNcut=5,lam=10480,Dlam=1,Npix=1,Verbose=0): 
    """
    Calculating the limiting flux given aparent AB magnitude of an object and the obtained (average) S/N
    Default is a 5sigma limiting magnitude for the MOSFIRE Y band at 10480 Angstrom 
    (http://www2.keck.hawaii.edu/inst/mosfire/filters.html)
    To get limiting _emissionline_ flux the wavelength width and number-of-pixels values can be used. 

    ---- INPUT ----
    mABobj  The AB magnitude of the object
    SNobj   The SN of object at the selected wavelength

    ---- OPTIONAL INPUT ----
    SNcut   The SN limiting magnitude to calculate
             DEFAULT = 5 sigma
    lam     The wavelength at which the limiting magnitude flux is estimated
    Dlam    The widht of a line. Used to turn limiting flux into emission line flux
             DEFAULT = 1 (i.e. not a line)
    Npix    The number of pixels the lines (Dlam) spans over
             DEFAULT = 1 (i.e. not a line)
    Verbose Set to 1 for verbosity

    ---- EXAMPLES ----
    import kbs.utilities as kbs
    Flim     = kbs.fluxlimit(19.29,4,SNcut=3,lam=13000,Dlam=1,Npix=1,Verbose=1)
    Flimline = kbs.fluxlimit(20.28,8,SNcut=3,lam=17993,Dlam=14,Npix=1,Verbose=1)
    """

    zpAB  = 48.60
    mlim  = float(mABobj)-2.5*math.log10(float(SNcut)/float(SNobj))
    Fnu   = 10**((mlim+zpAB)/-2.5)                   # in erg/s/cm2/Hz
    Flam  =  2.99792458e+18 * Fnu / lam**2           # in erg/s/cm2/A

    Flam  = Flam * Dlam / math.sqrt(Npix)            # potentially changing Flam to limiting line flux

    if Verbose == 1: 
        print('The ',SNcut,' sigma limiting flux was estimated to be ',Flam,' erg/s/cm^2/A')
        print('corresponding to a limiting magnitude of ',mlim)
    return Flam

#-------------------------------------------------------------------------------------------------------------
def getAv(RA,DEC,filter):
    """
    returns redening Av (and E(B-V)) for ra and dec (degrees; scalar or numpy arrays) for a given
    HST filter usning E(B-V) from the Schlegel dust maps (also returned)

       Avval, EBVval = kbs.getAv(51,219,'F125W')
    The extinction corrected apparent mag is then
       magband_corr = magband - Av_band

    Could also correct for extinction using:
        extlaw    = astropysics.obstools.CardelliExtinction(EBmV=extval, Rv=Rvval)
        magcorr   = extlaw.correctPhotometry(mag,bandwavelength)

    """
    if isinstance(RA,types.FloatType): 
        Nvals = 1
    elif isinstance(RA,types.IntType): 
        Nvals = 1
    else:
        Nvals = range(len(RA))

    if Nvals > 1:
        gall        = []
        galb        = []
        for ii in Nvals: # looping over RA and Decs and converting to galactic coordiantes
            gcoords = coords.ICRSCoordinates(RA[ii],DEC[ii]).convert(coords.GalacticCoordinates)
            gall.append(gcoords.l.degrees)
            galb.append(gcoords.b.degrees)
    else:
        gcoords = coords.ICRSCoordinates(RA,DEC).convert(coords.GalacticCoordinates)
        gall = gcoords.l.degrees
        galb = gcoords.b.degrees

    #dustmaps    = '/Users/kasperborelloschmidt/work/python/progs/dust_getvalV0p1/fits/SFD_dust_4096_%s.fits'
    dustmaps    = '/Users/kschmidt/work/dustmaps/SFD_dust_4096_%s.fits'
    Ebv         = astropysics.obstools.get_SFD_dust(gall,galb,dustmaps,interpolate=True) # redening from Schlegel maps

    av_ebv = {} # ebv2Av values for HST filters from Larry; CCM reddening curve with R_V = 3.1
    av_ebv['F300X']  = 6.78362003559
    av_ebv['F475X']  = 3.79441819047
    av_ebv['F475W']  = 3.82839055809
    av_ebv['F606W']  = 3.01882984135
    av_ebv['F600LP'] = 2.24159324026
    av_ebv['F098M']  = 1.29502816006
    av_ebv['F105W']  = 1.18148250758
    av_ebv['F125W']  = 0.893036743585
    av_ebv['F160W']  = 0.633710427959

    try:
        av_ebv[filter] 
    except KeyError:
        sys.exit(':: kbs.getAv :: The filter '+filter+' is not accepted as input --> ABORTING')

    Av          = av_ebv[filter] * Ebv

    return Av,Ebv
#-------------------------------------------------------------------------------------------------------------
def getAv_area(RAcenter,DECcenter,radius,filter='F125W',valreturn=None,stepsize=[0.5,0.5],verbose=True):
    """
    use kbs.getAv to get the dust extinction in a circular aperture of a given radius

    stepsize     Step size in degrees when generating grid to estimate extinction over
    valreturn    Determines the output. Chose between 'median', 'mean', ...
                 If 'None' all values will be returned

    A, EBV, grid = kbs.getAv_area(322.35977,-7.6908861,1.7/60.,filter='F125W',valreturn='median',stepsize=[0.5/60.,0.5/60.])

    """

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    center   = SkyCoord(ra=RAcenter*u.degree, dec=DECcenter*u.degree, frame='fk5')
    ra_grid  = np.arange(RAcenter,RAcenter+2*radius,stepsize[0])-radius
    dec_grid = np.arange(DECcenter,DECcenter+2*radius,stepsize[0])-radius
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Setting up grid to estimate extinction on')
    grid_coords = []
    for rag in ra_grid:
        for decg in dec_grid:
            c_check = SkyCoord(ra=rag*u.degree, dec=decg*u.degree, frame='fk5')
            sep     = center.separation(c_check)  # Differing frames handled correctly

            if sep < radius*u.deg:
                grid_coords.append([rag,decg])


    Ngridpoints = len(grid_coords)
    if Ngridpoints == 0:
        sys.exit(' Grid has 0 grid points. Check your radius ('+str(radius)+') and stepsize ('+str(stepsize)+')')
    coords_arr  = np.asarray(grid_coords)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ds9gridfile = './gridregions.reg'
    if verbose: print(' - Saving grid to DS9 region file in '+ds9gridfile)
    kbs.create_DS9region(ds9gridfile,coords_arr[:,0],coords_arr[:,1],clobber=True)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Getting extinction for '+str(Ngridpoints)+' gridpoints ')
    Avvals  = []
    EBVvals = []
    for coord in coords_arr:
        cc = SkyCoord(ra=coord[0]*u.degree, dec=coord[1]*u.degree, frame='fk5')
        if verbose: print '   Getting A and E(B-V) at (ra, dec) = (',str("%.8f" % cc.ra.deg),',',\
            str("%.8f" % cc.dec.deg),')  ',
        Avval, EBVval = kbs.getAv(cc.ra.deg,cc.dec.deg,filter)
        Avvals.append(Avval)
        EBVvals.append(EBVval)
        if verbose: print '  A, E(B-V) = ',str("%.2f" % Avval),str("%.2f" % EBVval)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if valreturn == 'median':
        if verbose: print ' - Returning median value on calculated grid points'
        Av_final   = np.median(np.asarray(Avvals))
        EBV_final  = np.median(np.asarray(EBVvals))
        gridreturn = Ngridpoints
    elif valreturn == 'mean':
        if verbose: print ' - Returning mean value on calculated grid points'
        Av_final   = np.mean(np.asarray(Avvals))
        EBV_final  = np.mean(np.asarray(EBVvals))
        gridreturn = Ngridpoints
    else:
        if verbose: print ' - Returning all values calculated'
        Av_final   = np.asarray(Avvals)
        EBV_final  = np.asarray(EBVvals)
        gridreturn = coords_arr

    return Av_final, EBV_final, gridreturn

#-------------------------------------------------------------------------------------------------------------
def create_DS9region(outputfile,ralist,declist,color='red',circlesize=0.5,textlist=None,clobber=False):
    """
    Generate a basic DS9 region file with circles around a list of coordinates

    ralist
    declist
    color
    size
    text

    """
    if type(color) is str:
        color = [color]*len(ralist)

    if not clobber:
        if os.path.isfile(outputfile):
            sys.exit('File already exists and clobber = False --> ABORTING')
    fout = open(outputfile,'w')

    fout.write("# Region file format: DS9 version 4.1 \nfk5\n")

    for rr, ra in enumerate(ralist):
        string = 'circle('+str(ra)+','+str(declist[rr])+','+str(circlesize)+'") # color='+color[rr]+' width=3 '

        if textlist is not None:
            string = string+' font="times 10 bold roman" text={'+textlist[rr]+'}'

        fout.write(string+' \n')

    fout.close()

#-------------------------------------------------------------------------------------------------------------
def magapp2abs(Mapp,zobj,RA,DEC,Av=-99,band='Jbradley2012',cos='WMAP7BAOH0'):
    """
    Converting apparent magnitude(s) into absolut magnitude(s)

    Av    : The extinction. If not given it's estimated from the Schlegel maps (time consuming)
            Note that RA and DEC is only used if Av is not given; otherwise they are 'dummys'
    
    band  : the band to do the calculations for. The default is to use the J band
            conversion used in Bradley et al. 2012. In this case the (extinction correted) 
            J-band magnitude is expected and MUV = MJ125 - 47.14 is returned. This 
            corresponds to 
                  Mabs = mobs - 5.0 * (np.log10(lumdist) - 1.0) + (2.5 * np.log10(1.0 + zobj))
            With k-correction (last term) assuming the source has a flat (beta = -2) SED using
            a 0.3 0.7 0.7 cosmology  
            NB! for band='Jbradley2012' zobj, RA, DEC and Av are all dummy values
    cos   : the cosmology to use, e.g.
            'WMAP7BAOH0' (Default) from http://lambda.gsfc.nasa.gov/product/map/dr4/params/lcdm_sz_lens_wmap7_bao_h0.cfm
            'WMAP7'                from http://lambda.gsfc.nasa.gov/product/map/dr4/params/lcdm_sz_lens_wmap7.cfm
    """    
    if band == 'Jbradley2012':
        Mabs          = np.array([Mapp - 47.14])
    else:
        cosmo = choose_cosmology(cos)
        Dlum          = coords.funcs.cosmo_z_to_dist(zobj, zerr=None, disttype='luminosity')*1e6 # luminosity distance in pc
        Kcorrection   = (2.5 * np.log10(1.0 + zobj)) # assumes source has flat (beta = -2) SED.  
                                                     # A bluer beta will likely give you an additional
                                                     # correction of about ~0.1 mag or so.
        if isinstance(Mapp,types.FloatType) and Av == -99: # if Av is -99, calculate it
            Av, Ebv = getAv(RA,DEC,band) 
            Mabs    = Mapp - 5*np.log10(Dlum)+5 + Kcorrection - Av # corrected absolut magnitude of objects
        else:
            Mabs    = None
    return Mabs
#-------------------------------------------------------------------------------------------------------------
def magabs2app(Mabs,zobj,RA,DEC,Av=-99,band=None,cos='WMAP7BAOH0'):
    """
    Converting absolute magnitude(s) into apparent magnitude(s)

    Av    : The extinction. If not given it's estimated from the Schlegel maps (time consuming)
            Note that RA and DEC is only used if Av is not given; otherwise they are 'dummys'
    
    band  : the band to do the calculations for. The default is to use the J band
            conversion used in Bradley et al. 2012. In this case the (extinction correted) 
            J-band magnitude is expected and MJ125 = MUV + 47.14 is returned. This 
            corresponds to inverting 
                  Mabs = mobs - 5.0 * (np.log10(lumdist) - 1.0) + (2.5 * np.log10(1.0 + zobj))
            With k-correction (last term) assuming the source has a flat (beta = -2) SED using
            a 0.3 0.7 0.7 cosmology.
            NB! for band='Jbradley2012' zobj, RA, DEC and Av are all dummy values
    cos   : the cosmology to use, e.g.
            'WMAP7BAOH0' (Default) from http://lambda.gsfc.nasa.gov/product/map/dr4/params/lcdm_sz_lens_wmap7_bao_h0.cfm
            'WMAP7'                from http://lambda.gsfc.nasa.gov/product/map/dr4/params/lcdm_sz_lens_wmap7.cfm
    """
    if band == 'Jbradley2012':
        Mapp          = np.array([Mabs + 47.14])
    else:
        cosmo = choose_cosmology(cos) 
        Dlum          = coords.funcs.cosmo_z_to_dist(zobj, zerr=None, disttype='luminosity')*1e6 # luminosity distance in pc
        Kcorrection   = (2.5 * np.log10(1.0 + zobj)) # assumes source has flat (beta = -2) SED.  
                                                     # A bluer beta will likely give you an additional
                                                     # correction of about ~0.1 mag or so.
        if isinstance(Mabs,types.FloatType) and Av == -99: # if Av is -99, calculate it
            Av, Ebv = getAv(RA,DEC,band) 
            Mapp    = Mabs + 5*np.log10(Dlum) - 5 - Kcorrection + Av # corrected absolut magnitude of objects
        else:
            Mapp    = None
    return Mapp
#-------------------------------------------------------------------------------------------------------------
def Mabs2L(Mabs,MUVsun=5.5):
    """
    Converting absolute magnitude(s) to luminosity in erg/s
    Using a default absolut magnitude of the sun (in UV) of 5.5 from http://www.ucolick.org/~cnaw/sun.html
    """
    Lsun        = 3.839e-11 # 1e44 erg/s
    Lobj        = 10**((MUVsun-Mabs)/2.5)*Lsun  # Luminosity in erg/s
    return Lobj
#-------------------------------------------------------------------------------------------------------------
def L2Mabs(Lobj,MUVsun=5.5):
    """
    Converting luminsoity 10^44 erg/s into absolute magnitude(s)
    Using a default absolut magnitude of the sun (in UV) of 5.5 from http://www.ucolick.org/~cnaw/sun.html
    """
    Lsun        = 3.839e-11 # 1e44 erg/s
    Mabs        = MUVsun - 2.5*np.log10(Lobj/Lsun)
    return Mabs
#-------------------------------------------------------------------------------------------------------------
def interpn(*args, **kw):
    """Interpolation on N-Dimensions 

    ai = interpn(x, y, z, ..., a, xi, yi, zi, ...)
    where the arrays x, y, z, ... define a rectangular grid
    and a.shape == (len(x), len(y), len(z), ...)

    KBS:
    Taken from http://projects.scipy.org/scipy/ticket/1727#comment:3
    An alternative is to use scipy.interpolate.LinearNDInterpolator
    but slow according to http://stackoverflow.com/questions/14119892/python-4d-linear-interpolation-on-a-rectangular-grid
    and had problems getting it to work

    -- OPTIONAL INPUT --
    method     the interpolation method to use. Options are 'linear','nearest', 'zero', 'slinear', 'quadratic', 'cubic'

    -- EAXMPLE --
    newy = kbs.interpn(oldx,oldy,newx)

    """
    method = kw.pop('method', 'linear')
    if kw:
        raise ValueError("Unknown arguments: " % kw.keys())
    nd = (len(args)-1)//2
    if len(args) != 2*nd+1:
        raise ValueError("Wrong number of arguments")
    q = args[:nd]
    qi = args[nd+1:]
    a = args[nd]
    for j in range(nd):
        a = interp1d(q[j], a, axis=j, kind=method)(qi[j])
    return a
#-------------------------------------------------------------------------------------------------------------
def simulate_schechter_distribution(alpha, L_star, L_min, N,trunmax=10):
    """ 
    Generate N samples from a Schechter distribution, which is like a gamma distribution 
    but with a negative alpha parameter and cut off on the left somewhere above zero so that
    it converges.
        
    If you pass in stupid enough parameters then it will get stuck in a loop forever, and it
    will be all your own fault.
        
    Based on algorithm in http://www.math.leidenuniv.nl/~gill/teaching/astro/stanSchechter.pdf

    KBS:-------------------------------------------------------------------------------------
          Code taken from https://gist.github.com/joezuntz/5056136 and modified.
          Schechter distribution with -1 < alpha+1 (k) < -0

          trunmax : To prevent an invinite loop trunmax gives the maximum allowed run time [s].
                    If this time is surpased any found entries are retured or an array of 0s
        -------------------------------------------------------------------------------------
    """
    output = []
    n      = 0
    Nvals  = N
    t0     = time.time()    
    while n<N:
        t1   = time.time()
        Lgam = np.random.gamma(scale=L_star, shape=alpha+2, size=N)  # drawing values from gamma dist with k+1
        Lcut = Lgam[Lgam>L_min]                                      # removing L values from gamma dist > L_min
        ucut = np.random.uniform(size=Lcut.size)                     # random values [0:1]
        Lval = Lcut[ucut<L_min/Lcut]                                 # only keeping L values where ucut < L_min/L
        output.append(Lval)                                          # append thes to output array
        n+=Lval.size                                                 # increase counter

        if (t1-t0) > trunmax:                                        # check that runtime is not too long
            Nvals = n                                                # set Nvals to found values
            if Nvals < 2.: 
                output.append(np.zeros(N))                           # if not even 2 values were found return array of 0s
                Nvals  = N                                           # updating Nvals
            n += N-n                                                 # make sure loop ends
    values = np.concatenate(output)[:Nvals]                          # generate output by reformatting
    return values
#-------------------------------------------------------------------------------------------------------------
def sex2deg(rasex,decsex):
    """ 
    converting ra and dec strings with sexagesimal values to float of degrees using skycoor

    """
    skycoorout = commands.getoutput('skycoor -d '+rasex+' '+decsex)    
    outsplit   = skycoorout.split()
    radeg      = float(outsplit[0])
    decdeg     = float(outsplit[1])
    return radeg, decdeg
#-------------------------------------------------------------------------------------------------------------
def correctmag4extinction(mag,ebv,band,Rvval=3.1,verbose=1):
    """
    Correcting a numpy array of magnitudes with given E(B-V) for galactic extionction
    Using the Cardelli et al. 1989 extintion law

    Example of usage:
       import kbsutilities as kbs
       mag  = np.array([26.75,26.25,25.88,26.36,27.31,26.97,27.32,25.77,27.07,26.91,27.64,26.67,27.04,27.16])
       ebv  = np.array([0.038,0.013,0.013,0.028,0.08,0.07,0.083,0.013,0.046,0.046,0.046,0.046,0.046,0.026])
       band = np.asarray(['F125W']*len(mag))
       magcorr = kbs.correctmag4extinction(mag,ebv,band,Rvval=3.1,verbose=1)
    """
    extval  = ebv
    extlaw  = astropysics.obstools.CardelliExtinction(EBmV=extval, Rv=Rvval)
    Nmag    = len(mag)
    wave    = np.zeros(Nmag)

    if verbose == 1: print 'Correcting the '+str(Nmag)+' magnitudes using Cardelli ext. law with Rv=',Rvval
    
    for ii in xrange(Nmag):
        bandobj = band[ii]

        if   bandobj == 'F606W': 
            wave[ii] = 5887.40
        elif bandobj == 'F098M':
            wave[ii] = 9864.10            
        elif bandobj == 'F125W':
            wave[ii] = 12486.00
        elif bandobj == 'F160W':
            wave[ii] = 15369.00
        else:
            sys.exit('The band '+str(bandobj)+' for magnitude number '+str(ii+1)+' is not a valid choice --> ABORTING')

    magcorr = extlaw.correctPhotometry(mag,wave)    
    if verbose == 1: 
        for jj in xrange(Nmag):
            print '   Corrected '+str(mag[jj])+' in '+band[jj]+' to    '+str(magcorr[jj])
    return magcorr
#-------------------------------------------------------------------------------------------------------------
def drawnbinom(n,p,size=1):
    """
    Drawing value from binomial distribution.
    Draw done using the prescription in Equation (C2) of Kelly et al. (2008)

    Note that the returned value is a float as opposed to the values returned by
        nbinomdraw  = scipy.stats.nbinom.rvs(N, p, size=1)
        nbinomdraw  = np.random.negative_binomial(N, p, size=1)
    This enables draws from distributions with very small probabilities where the
    long integers fail (returns -2147483648) in the above methods.
    
    Parameters
    ----------
        n : Number of objects in sample
        p : Probability of drawing an object from parent sample

    Returns
    -------
        A numpy array of size 'size' containing draws
    """
    t0   = time.time()    
    m    = np.zeros(size) # array to contain draws
    nlim = 1.0e6 # the limit deciding what approach to take

    if n <= nlim: # for very large n the u-vector becomes too large
        for ii in xrange(size):
            u = np.random.uniform(low=0.0, high=1.0, size=n)
            m[ii] = int(np.sum(np.log10(u)/np.log10(1-p)))    
    else:       # so in that case slit up in sub samples instead
        #print 'n > nlim = ',nlim
        nusub = np.floor(n/nlim) # sub arrays of size Nlim to draw
        nrest = n - nlim*nusub   # size of final subarray

        for ii in xrange(size):
            madd = 0.0
#            for jj in xrange(int(nusub)):
            for jj in xrange(1): # only doing it once and then multiplying by number of sub arrays
                usub  = np.random.uniform(low=0.0, high=1.0, size=nlim)
                msub  = int(np.sum(np.log10(usub)/np.log10(1-p)))    
                madd  = madd + msub * nusub
            if nrest > 0:
                urest = np.random.uniform(low=0.0, high=1.0, size=nrest)
                mrest = int(np.sum(np.log10(urest)/np.log10(1-p)))    
                madd  = madd + mrest
            
            m[ii]    = int(madd)

# looping is way too slow ...
#        sumval = 0.0
#        for ii in xrange(size):
#            for jj in xrange(int(n)):
#                u = np.random.uniform(low=0.0, high=1.0, size=1)
#                sumval = sumval + np.log10(u)/np.log10(1-p)
#            m[ii] = int(sumval)
        
    N = m + n
    #print 'Took ',time.time()-t0,' sec.'
    return N
#-------------------------------------------------------------------------------------------------------------
def appendfitstable(tab1,tab2,newtab='kbs_appendfitstable_results.fits',clobber=False):
    """
    Appending 1 fits table to another.
    It is assumed that the two tables contain the same columns.
    see http://pythonhosted.org/pyfits/users_guide/users_table.html#appending-tables

    Note that columns with object IDs are also added, hence, the be aware of duplicate ids

    Parameters
    ----------
        tab1     primariy fits table
        tab2     fits table to append to tab1
                 (should contain the same columns as tab1)
        clobber  overwrite output file if it already exists

    Returns
    -------
        the name 'newtab' of the created table

    Example
    -------
    import kbsutilities as kbs
    tab1   = 'simulatedsamples/dataarraySim_pdistschechter_Ntot1000_k-0p5_Lstar0p5_LJlim0p1_Nobj17.fits'
    tab2   = 'simulatedsamples/dataarraySim_pdistschechter_Ntot2000_k-0p5_Lstar0p5_LJlim0p1_Nobj25.fits'
    newtab = 'simulatedsamples/testname.fits'
    output = kbs.appendfitstable(tab1,tab2,newtab=newtab)

    """
    t1     = pyfits.open(tab1)
    t2     = pyfits.open(tab2)

    nrows1 = t1[1].data.shape[0] # counting rows in t1
    nrows2 = t2[1].data.shape[0] # counting rows in t2

    nrows  = nrows1 + nrows2 # total number of rows in the table to be generated
    hdu    = pyfits.new_table(t1[1].columns, nrows=nrows)

    for name in t1[1].columns.names:
        hdu.data.field(name)[nrows1:]=t2[1].data.field(name)

    hdu.writeto(newtab,clobber=clobber)

    return newtab
#-------------------------------------------------------------------------------------------------------------
def confcontours(xpoints,ypoints,binx=200,biny=200):
    """
    Function estimating confidence contours for a given 2D distribution of points.

    @return: gridsigma, extent

    which can be plotted with for instance
    plt.contour(gridsigma.transpose(),[1,2,3],extent=extent,origin='lower',colors=['r','r','r'],label='contours',zorder=5)
    """
    from fast_kde import fast_kde # used to create confidence curves for contours
    xmin        = np.min(xpoints)
    xmax        = np.max(xpoints)
    ymin        = np.min(ypoints)
    ymax        = np.max(ypoints)
    extent      = [xmax,xmin,ymin,ymax]

    Nval        = binx*biny

    kde_grid    = fast_kde(ypoints,xpoints, gridsize=(binx,biny), weights=None,extents=[ymin,ymax,xmin,xmax])

    binarea     = (xmax-xmin)/binx * (ymax-ymin)/biny
    kde_int     = kde_grid * binarea # ~integrated value in grid
    kde_flat    = np.ravel(kde_int)
    sortindex   = np.argsort(kde_int,axis=None)[::-1]
    gridsigma   = np.zeros((binx,biny))

    sum = 0.0
    for ss in xrange(Nval):
        xx  = np.where(kde_int == kde_flat[sortindex[ss]])
        sum = sum + np.sum(kde_int[xx])
        if (sum < 0.68): gridsigma[xx] = 1.0
        if (sum > 0.68) and (sum < 0.95): gridsigma[xx] = 2.0
        if (sum > 0.95) and (sum < 0.99): gridsigma[xx] = 3.0

    return gridsigma, extent

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def createLaTeXtab(fitstable,verbose=False):
    """
    Turn binary fits table into LaTeX table

    --- INPUT ---
    fitstable   : path and name to fits catalog to turn into table

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    cat = 'dropoutselection131216_fullIR/hlsp_clash_hst_ir_macs0717_cat_reformat_flux_AB25_noF160Wcut_Bdrops_V09.fits'
    table = kbs.createLaTeXtab(cat,verbose=True)

    """
    if os.path.exists(fitstable):
        if verbose: print ' - Loading fitstable in ',fitstable
        dat   = pyfits.open(fitstable)
        datTB = dat[1].data
        cols  = datTB.columns.names
        Ncol  = len(cols)
        Nobj  = len(datTB[cols[0]])
        if verbose: print ' - Found ',Nobj,' objects to create table rows for'
        if verbose: print '   and ',Ncol,' columns of data'
    else:
        sys.exit('The provided catalog does not exist --> ABORTING')

    rows = '\n'
    for oo in xrange(Nobj):
        rowlist = [str(np.asarray(datTB[oo])[ii]) for ii in xrange(Ncol)]
        rows = rows + ' & '.join(rowlist)+'\\\\ \n'

    colheads = '\colhead{' + '} & \colhead{'.join(cols)+'}'
    colheads = colheads.replace('_','\\_')

    tablestr  = """
    \\tabletypesize{\\tiny}
    \\begin{deluxetable*}{%s}
    \\tablecolumns{%s}
    \\tablewidth{0pt}
    \\tablecaption{Title of table...}
    \\tablehead{%s}
    \startdata
    %s
    \enddata
    \\tablecomments{Comments...\\\\
    \\tnm{a}{table note...}\\\\
    }
    \label{tab:label}
    \end{deluxetable*}
    """ % (''.join(['c']*Ncol),Ncol,colheads,rows)

    return tablestr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def createLaTeXtab_CLASH(fitsCLASHcat,IDs,verbose=False,calcstuff=False):
    """
    Turn binary fits table into LaTeX table

    --- INPUT ---
    fitsCLASHcat   : path and name to fits CLASH catalog to use to turn into table
    IDs            : list of object ids to include in table

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    cat   = '/Users/kasperborelloschmidt/work/GLASS/MACS0717test/hlsp_clash_hst_ir_macs0717_cat.txt.FITS'
    table = kbs.createLaTeXtab_CLASH(cat,[859,1730],verbose=True)

    """
    if os.path.exists(fitsCLASHcat):
        if verbose: print ' - Loading CLASH catalog in ',fitsCLASHcat
        dat   = pyfits.open(fitsCLASHcat)
        datTB = dat[1].data

        Nobj    = len(IDs)
        entries = np.zeros(Nobj).astype(dtype=int)-99
        if verbose: print ' - Grabbing data for requested objects'
        for ii in xrange(Nobj):
            objent = np.where(datTB['id'] == IDs[ii])[0]
            if len(objent) == 0:
                if verbose: print '   No match to ID ',IDs[ii]
            elif len(objent) > 1:
                if verbose: print '   More than one match to ID ',IDs[ii],' --> skipping'
            else:
                entries[ii] = int(objent)

        if Nobj != len(np.where(entries != -99)[0]):
            sys.exit('Did not manage to extract data for all the '+str(Nobj)+' objects requested --> ABORTING')

        datTB = datTB[:][entries]
        cols  = datTB.columns.names
        Ncol  = len(cols)
        Nobj  = len(datTB[cols[0]])
        if verbose: print ' - Extracted ',Ncol,' colukmns for the requested ',Nobj,' objects.'
    else:
        sys.exit('The provided catalog does not exist --> ABORTING')

    if verbose: print ' - Putting data rows together.'
    headlist = ['ID','$\\alpha_\\textrm{J2000}$','$\\delta_\\textrm{J2000}$','F225W mag', 'F275W mag', 'F336W mag', 'F390W mag', 'F435W mag', 'F475W mag', 'F555W mag', 'F606W mag', 'F625W mag', 'F775W mag', 'F814W mag', 'F850lp mag', 'F105W mag', 'F110W mag', 'F125W mag', 'F140W mag', 'F160W mag', '$z_\\textrm{phot}$']

    Ncolfinal = len(headlist)
    colheads = '\colhead{' + '} & \colhead{'.join(headlist)+'}'

    fmt  = '%.2f' # format for numeric values
    fmtz = '%.1f' # format for numeric values
    rows = '\n'
    for oo in xrange(Nobj):
        rows = rows+str(int(datTB['id'][oo]))+' & '+ \
        str('%.6f' % datTB['ra'][oo])+' & '+ \
        str('%.6f' % datTB['dec'][oo])+' & ' \
        '$ '+str(fmt % datTB['f225w_mag'][oo])+' \\pm '+str(fmt % datTB['f225w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f275w_mag'][oo])+' \\pm '+str(fmt % datTB['f275w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f336w_mag'][oo])+' \\pm '+str(fmt % datTB['f336w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f390w_mag'][oo])+' \\pm '+str(fmt % datTB['f390w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f435w_mag'][oo])+' \\pm '+str(fmt % datTB['f435w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f475w_mag'][oo])+' \\pm '+str(fmt % datTB['f475w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f555w_mag'][oo])+' \\pm '+str(fmt % datTB['f555w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f606w_mag'][oo])+' \\pm '+str(fmt % datTB['f606w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f625w_mag'][oo])+' \\pm '+str(fmt % datTB['f625w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f775w_mag'][oo])+' \\pm '+str(fmt % datTB['f775w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f814w_mag'][oo])+' \\pm '+str(fmt % datTB['f814w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f850lp_mag'][oo])+' \\pm '+str(fmt % datTB['f850lp_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f105w_mag'][oo])+' \\pm '+str(fmt % datTB['f105w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f110w_mag'][oo])+' \\pm '+str(fmt % datTB['f110w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f125w_mag'][oo])+' \\pm '+str(fmt % datTB['f125w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f140w_mag'][oo])+' \\pm '+str(fmt % datTB['f140w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f160w_mag'][oo])+' \\pm '+str(fmt % datTB['f160w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmtz % datTB['zb'][oo]) + \
        ' ^{+'+str(fmtz % (datTB['zbmax'][oo]-datTB['zb'][oo])) + \
        '}_{-'+str(fmtz % (datTB['zb'][oo]-datTB['zbmin'][oo]))+' }$ \\\\ \n'

    tablestr  = """
    \\tabletypesize{\\tiny}
    \\begin{deluxetable*}{%s}
    \\tablecolumns{%s}
    \\tablewidth{0pt}
    \\tablecaption{Title of table...}
    \\tablehead{%s}
    \startdata
    %s
    \enddata
    \\tablecomments{Comments...\\\\
    \\tnm{a}{table note...}\\\\
    }
    \label{tab:label}
    \end{deluxetable*}
    """ % (''.join(['c']*Ncolfinal),Ncolfinal,colheads,rows)

    # ===== Calculate stuff for objects =====
    if calcstuff==True:
        if verbose: print ' - Calculating stuff for objects:'
        for oo in xrange(Nobj):
            mags  = np.array([datTB['f125w_mag'][oo],datTB['f140w_mag'][oo],datTB['f160w_mag'][oo]])
            errs  = np.array([datTB['f125w_magerr'][oo],datTB['f140w_magerr'][oo],datTB['f160w_magerr'][oo]])

            magav = np.average(mags[mags != -99],weights=1/errs[mags != -99]**2) # wighted average of magnitudes
            #if verbose: print '   obj,magav,zphot,mags[mags!=-99] = ',IDs[oo],magav,datTB['zb'][oo],mags[mags != -99]

            #if verbose: print IDs[oo],datTB['ra'][oo],datTB['dec'][oo]
            #if verbose: print IDs[oo],'   ',str(fmt % datTB['f105w_mag'][oo]),',',str(fmt % datTB['f105w_magerr'][oo])
            if verbose: print IDs[oo],'   ',str(fmt % datTB['f140w_mag'][oo]),',',str(fmt % datTB['f140w_magerr'][oo])


    return tablestr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def createSCImCONTAM(spec2Dfits,verbose=False,clobber=False):
    """
    Read a 2D fits spectrum from the 3D-HST/GLASS reduction and produce a
    new fits image containing the contamination subtracted 2D spectrum.

    --- INPUT ---
    spec2Dfits   : path and anme to 2D spectrum (ssumes it contains the extensions SCI and CONTAM)

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    outputimg = kbs.createSCImCONTAM(spec2Dfits,verbose=True)

    """
    if os.path.isfile(spec2Dfits):
        if verbose: print ' - Loading fits file'
        twod       = pyfits.open(spec2Dfits)
        SCI        = twod['SCI'].data
        CONTAM     = twod['CONTAM'].data
        SCImCONTAM = SCI - CONTAM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def void2file(npvoid,outputfile,verbose=False,clobber=False):
    """
    Wirting a numpy void (from np.genfromtxt) to a text file
    """

    if os.path.isfile(outputfile) and (clobber == False):
        sys.exit(' clobber == false so not creating '+outputfile+' as it already exists')
    else:
        f = open(outputfile,'w')
        colnames = npvoid.dtype.names
        hdr = '# '+str(colnames).replace(',',' ').replace(')',' ').replace('(',' ').replace("'",' ').replace("'",' ').replace(']',' ').replace('[',' ')+'\n'
        f.write(hdr)
        for ii in xrange(len(npvoid[colnames[0]])):
            outstr = str(npvoid[ii].tolist()).replace(',',' ').replace(')',' ').replace('(',' ').replace("'",' ').replace("'",' ').replace(']',' ').replace('[',' ')
            f.write("%s\n" % outstr)

        f.close()
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def crossmatch(ralist,declist,
               catalog='skelton_goodss',idcol='id',racol='ra',deccol='dec',catext=1,verbose=True):
    """
    Crossmatch a list of coordinates to a catalog and return IDs, ra, dec, and r_match[arcsec]

    --- INPUT ---
    ralist        List of RA to crossmatch to catalog
    declist       List of Dec to crossmatch to catalog
    catalog       Fits catalog to find matches in and return ID, coordinates and r_match
    idcol         ID column name in fits catalog
    racol         RA column name in fits catalog
    deccol        Dec column name in fits catalog
    catext        Fits extension containing catalog data in fits catalog
    verbose       Toggle verbosity

    --- EXAMPLE OF USE ---

    id, ra, dec, rmatch = kbs.crossmatch([53.090976],[-27.957849]) # object 12 in skelton catalog

    """
    if len(ralist) != len(declist):
        sys.exit('The provided "ralist" and "declist" have different lengths --> ABORTING')
    else:
        Nobj = len(ralist)
        if verbose: print ' - Will return best match to the '+str(Nobj)+' coordinate sets provide '
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading catalog to get IDs'
    if catalog == 'skelton_goodss':
        catpath   = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/'
        catmatch  = catpath+'goodss_3dhst.v4.1.cat.FITS'
    else:
        catmatch  = catalog

    catdat = pyfits.open(catmatch)[catext].data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Finding IDs for provided objects '
    objids        = np.zeros(Nobj)-99.0
    ra_objids     = np.zeros(Nobj)-99.0
    dec_objids    = np.zeros(Nobj)-99.0
    rmatch_objids = np.zeros(Nobj)-99.0
    for rr in xrange(Nobj):
        objra              = ralist[rr]
        objdec             = declist[rr]
        if verbose:
            infostr = '   Finding ID for (ra,dec) = ('+str("%12.8f" % objra)+', '+str("%12.8f" % objdec)+')'
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        rmatch             = np.sqrt( (np.cos(np.deg2rad(objdec))*(catdat[racol]-objra))**2.0 +
                                      (catdat[deccol]-objdec)**2.0 )
        objent             = np.where(rmatch == np.min(rmatch))[0]
        objids[rr]         = catdat[idcol][objent]
        ra_objids[rr]      = catdat[racol][objent]
        dec_objids[rr]     = catdat[deccol][objent]
        rmatch_objids[rr]  = rmatch[objent]*3600. # rmatch in arcsec

    if verbose: print '\n   ... done '
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Returning ID, RA, Dec, r_match[arcsec] '
    return objids, ra_objids, dec_objids, rmatch_objids
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def crossmatch2cat(radeccat='/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_candels-cdfs_v0.2.fits',
                   idcol='ID',racol='RA',deccol='DEC',catext=1,
                   matchcat='skelton_goodss',
                   m_idcol='id',m_racol='ra',m_deccol='dec',m_catext=1,IDstrings=False,
                   writetofile='./kbscrossmatch_defaultname',clobber=False,verbose=True):
    """

    Generate catlogs of crossmatches to a given fits catalog

    --- INPUT ---

    radeccat      Fits catalog with RA and Dec to crossmatch to matchcat
    idcol         ID column name in fits radeccat
    racol         RA column name in fits radeccat
    deccol        Dec column name in fits radeccat
    catext        Fits extension containing catalog data in fits radeccat
    rmatchcat     Fits catalog to extract crossmatches and IDs from
    m_idcol       ID column name in fits rmatchcat
    m_racol       RA column name in fits rmatchcat
    m_deccol      Dec column name in fits rmatchcat
    m_catext      Fits extension containing catalog data in fits rmatchcat
    IDstrings     Save IDs as strings as opposed to the default integers
    writetofile   Generate ascii and fits output of crossmatches
                  If writetofile='None' nothing will be written, and the crossmatach will just be returned
    clobber       Overwrite existing output files?
    verbose       Toggle verbosity

    --- EXAMPLE OF USE ---

    radeccat = '/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_candels-cdfs_v0.2.fits'
    matchcat = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
    id, ra, dec, rmacth = kbs.crossmatch2cat(radeccat,matchcat=matchcat,writetofile='./MUSEcdfs_cm2_3DHSTgoods',clobber=False)

    """
    catdat    = pyfits.open(radeccat)[catext].data
    objs_id   = catdat[idcol]
    objs_ra   = catdat[racol]
    objs_dec  = catdat[deccol]

    id_match, ra_match, dec_match, r_match = kbs.crossmatch(objs_ra,objs_dec,catalog=matchcat,
                                                            idcol=m_idcol,racol=m_racol,deccol=m_deccol,
                                                            catext=m_catext,verbose=verbose)

    if writetofile != 'None':
        asciifile = writetofile+'.txt'

        if os.path.isfile(asciifile) and (clobber == False):
            if verbose: print ' - WARNING: '+asciifile+' exists and clobber=False so no file generated'
        else:
            if verbose: print ' - Generating '+asciifile
            fout = open(asciifile,'w')
            fout.write('# crossmatch using kbsutilities.crossmatch2cat() and '
                       'kbsutilities.crossmatch() on '+kbs.DandTstr2()+' \n')
            fout.write('# crossmatched objects in '+radeccat+' to '+matchcat+'\n')
            fout.write('# ID ra dec   ID_match ra_match dec_match r_match_arcsec \n')
            for ii, objid in enumerate(objs_id):
                if IDstrings:
                    IDoutstring       = str("%20s"   % objid)
                    IDoutstring_match = str("%20s"   % id_match[ii])
                else:
                    IDoutstring       = str("%20i"   % int(objid))
                    IDoutstring_match = str("%20i"   % int(id_match[ii]))

                objstr  = IDoutstring                    +'  '+\
                          str("%16.8f" % objs_ra[ii])    +'  '+\
                          str("%16.8f" % objs_dec[ii])   +'  '+\
                          IDoutstring_match              +'  '+\
                          str("%16.8f" % ra_match[ii])   +'  '+\
                          str("%16.8f" % dec_match[ii])  +'  '+\
                          str("%16.8f" % r_match[ii])
                fout.write(objstr+'\n')
            fout.close()

        fitsfile  = writetofile+'.fits'
        if os.path.isfile(fitsfile) and (clobber == False):
            if verbose: print ' - WARNING: '+fitsfile+' exists and clobber=False so no file generated'
        else:
            if verbose: print ' - Generating '+fitsfile
            fitspath   = kbs.pathAname(asciifile)[0]

            if IDstrings:
                ffmt = ['A20','D','D','A20','D','D','D']
            else:
                ffmt = ['J','D','D','J','D','D','D']

            outputfile = f2a.ascii2fits(asciifile,asciinames=True,skip_header=2,outpath=fitspath,
                                        verbose=verbose,fitsformat=ffmt)

    return id_match, ra_match, dec_match, r_match

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def generate_mergedPDF(searchstring='./*.pdf',outputdocument='./allpdfs.pdf',verbose=True):
    """
    Genrate PDF containing all individual PDF plots found by globbing on provided search string

    --- INPUT ---
    searchstring    String used with glob to find files to merger
    outputdocument  Name of document to store mergerd file in
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---

    searchs = './spectraplots160325/*.pdf'
    outdoc  = './spectraplots160325_allpdfs.pdf'
    kbs.generate_mergedPDF(searchstring=searchs,outputdocument=outdoc)

    """
    pdf_files = glob.glob(searchstring)

    if len(pdf_files) > 1:
        if verbose: print ' - Found '+str(len(pdf_files))+' to merge ',
        merger    = PdfFileMerger()

        if verbose: print ' ... merging ',
        for fname in pdf_files:
            merger.append(PdfFileReader(fname, "rb"))

        merger.write(outputdocument)
        if verbose: print '   Worte output to '+outputdocument
    else:
        if verbose: print ' - Found less than 2 pdf files globbing for '+searchstring+' --> No output generated'
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def convert_pdf2png(pdffile,res=300,clobber=False,quality='75',resize='100%',verbose=True):
    """
    Converting a pdf image to png using command line conversion

    --- EXAMPLE OF USE ---

    kbs.convert_pdf2png('./spectraplots160325/10210079_3dhstSpecAndInfo.pdf')
    """
    output = pdffile.replace('.pdf','.png')
    if os.path.isfile(output) and (clobber == False):
        print ' WARNING: '+output+' exists and clobber=False so skipping '
    else:
        convertcmd = 'convert '+pdffile+' -quality '+quality+' '+output
        cmdout = commands.getoutput(convertcmd)
        if verbose and cmdout != '':
            print '--> convert output:'
            print cmdout
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def convert_wavelength(lambdainput,version='air2vac',verbose=True):
    """
    Converting array of wavelengths either:
     - from vacuum to air (dry air at 1 atm pressure and 15C with 0.045% CO2 by volume following Morton et al. 2000)
     - from air to vacuum (following the VALD3 tools http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
                           solution derived by N. Piskunov)

    input wavelength should be given in Angstrom.

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    import numpy as np
    lam_vac  = np.array([1215.670,5008.208])
    lam_air  = kbs.convert_wavelength(lam_vac,version='vac2air')
    lam_vac2 = kbs.convert_wavelength(lam_air,version='air2vac')

    """
    if version == 'air2vac': # expression from  N. Piskunov (http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion)
        s = 10**4 / lambdainput
        n = 1.0 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + \
            0.0001599740894897 / (38.92568793293 - s**2)
        lambdaout = lambdainput * n

    elif version == 'vac2air': # expression from Morton et al. (2000) ApJS, 130:403
        s = 10**4 / lambdainput
        n = 1.0 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)
        lambdaout = lambdainput / n

    else:
        sys.exit('Invalid "version" provided; choose between "air2vac" and "vac2air"')

    return lambdaout
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def velocityoffset2dwave(redshift,voffset,lam_rest,verbose=True):
    """
    Converting a velcoity offset into a wavelength difference for a given rest-frame wavelength
    E.g., what does a 340 km/s velocity offset between Lya and CIII] for a redshift z_Lya=7.733
    source (Stark et al. 2016) mean for the wavelength shift of CIII] wrt. the wavelength predicted
    based on the Lya redshift?

    lam_obs, lam_offset, dlam = kbs.velocityoffset2dwave(7.733,340.0,1909.0)

    """

    cc         = 299792.458 # km/s
    lam_obs    = (redshift + 1.0) * lam_rest

    z_offset   = (voffset * (redshift + 1.0) / cc ) # systemic redshift if redshift = z_lya
    lam_offset = (redshift - z_offset + 1.0 ) * lam_rest

    dlam       = lam_obs - lam_offset

    if verbose:
        print ' - For a line at '+str(lam_rest)+'A from a redshift '+str(redshift)+' object '
        print '   the predicted observed wavelength of the line would be             : ',lam_obs,' A '
        print ' - Accounting for a velcoity offset of '+str(voffset)+'km/s '
        print '   the predicted observed wavelength of the line becomes              : ',lam_offset,' A '
        print '   which corresponds to an expected wavelength shift of the line of   : ',dlam,' A '

    return lam_obs, lam_offset, dlam
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def matplotlib_colornames():
    """
    Returning dictionary containg the names and hex values for defined matplotlib colors.
    Also see: http://matplotlib.org/examples/color/named_colors.html
    """
    cnames = {
    'aliceblue':            '#F0F8FF',
    'antiquewhite':         '#FAEBD7',
    'aqua':                 '#00FFFF',
    'aquamarine':           '#7FFFD4',
    'azure':                '#F0FFFF',
    'beige':                '#F5F5DC',
    'bisque':               '#FFE4C4',
    'black':                '#000000',
    'blanchedalmond':       '#FFEBCD',
    'blue':                 '#0000FF',
    'blueviolet':           '#8A2BE2',
    'brown':                '#A52A2A',
    'burlywood':            '#DEB887',
    'cadetblue':            '#5F9EA0',
    'chartreuse':           '#7FFF00',
    'chocolate':            '#D2691E',
    'coral':                '#FF7F50',
    'cornflowerblue':       '#6495ED',
    'cornsilk':             '#FFF8DC',
    'crimson':              '#DC143C',
    'cyan':                 '#00FFFF',
    'darkblue':             '#00008B',
    'darkcyan':             '#008B8B',
    'darkgoldenrod':        '#B8860B',
    'darkgray':             '#A9A9A9',
    'darkgreen':            '#006400',
    'darkkhaki':            '#BDB76B',
    'darkmagenta':          '#8B008B',
    'darkolivegreen':       '#556B2F',
    'darkorange':           '#FF8C00',
    'darkorchid':           '#9932CC',
    'darkred':              '#8B0000',
    'darksalmon':           '#E9967A',
    'darkseagreen':         '#8FBC8F',
    'darkslateblue':        '#483D8B',
    'darkslategray':        '#2F4F4F',
    'darkturquoise':        '#00CED1',
    'darkviolet':           '#9400D3',
    'deeppink':             '#FF1493',
    'deepskyblue':          '#00BFFF',
    'dimgray':              '#696969',
    'dodgerblue':           '#1E90FF',
    'firebrick':            '#B22222',
    'floralwhite':          '#FFFAF0',
    'forestgreen':          '#228B22',
    'fuchsia':              '#FF00FF',
    'gainsboro':            '#DCDCDC',
    'ghostwhite':           '#F8F8FF',
    'gold':                 '#FFD700',
    'goldenrod':            '#DAA520',
    'gray':                 '#808080',
    'green':                '#008000',
    'greenyellow':          '#ADFF2F',
    'honeydew':             '#F0FFF0',
    'hotpink':              '#FF69B4',
    'indianred':            '#CD5C5C',
    'indigo':               '#4B0082',
    'ivory':                '#FFFFF0',
    'khaki':                '#F0E68C',
    'lavender':             '#E6E6FA',
    'lavenderblush':        '#FFF0F5',
    'lawngreen':            '#7CFC00',
    'lemonchiffon':         '#FFFACD',
    'lightblue':            '#ADD8E6',
    'lightcoral':           '#F08080',
    'lightcyan':            '#E0FFFF',
    'lightgoldenrodyellow': '#FAFAD2',
    'lightgreen':           '#90EE90',
    'lightgray':            '#D3D3D3',
    'lightpink':            '#FFB6C1',
    'lightsalmon':          '#FFA07A',
    'lightseagreen':        '#20B2AA',
    'lightskyblue':         '#87CEFA',
    'lightslategray':       '#778899',
    'lightsteelblue':       '#B0C4DE',
    'lightyellow':          '#FFFFE0',
    'lime':                 '#00FF00',
    'limegreen':            '#32CD32',
    'linen':                '#FAF0E6',
    'magenta':              '#FF00FF',
    'maroon':               '#800000',
    'mediumaquamarine':     '#66CDAA',
    'mediumblue':           '#0000CD',
    'mediumorchid':         '#BA55D3',
    'mediumpurple':         '#9370DB',
    'mediumseagreen':       '#3CB371',
    'mediumslateblue':      '#7B68EE',
    'mediumspringgreen':    '#00FA9A',
    'mediumturquoise':      '#48D1CC',
    'mediumvioletred':      '#C71585',
    'midnightblue':         '#191970',
    'mintcream':            '#F5FFFA',
    'mistyrose':            '#FFE4E1',
    'moccasin':             '#FFE4B5',
    'navajowhite':          '#FFDEAD',
    'navy':                 '#000080',
    'oldlace':              '#FDF5E6',
    'olive':                '#808000',
    'olivedrab':            '#6B8E23',
    'orange':               '#FFA500',
    'orangered':            '#FF4500',
    'orchid':               '#DA70D6',
    'palegoldenrod':        '#EEE8AA',
    'palegreen':            '#98FB98',
    'paleturquoise':        '#AFEEEE',
    'palevioletred':        '#DB7093',
    'papayawhip':           '#FFEFD5',
    'peachpuff':            '#FFDAB9',
    'peru':                 '#CD853F',
    'pink':                 '#FFC0CB',
    'plum':                 '#DDA0DD',
    'powderblue':           '#B0E0E6',
    'purple':               '#800080',
    'red':                  '#FF0000',
    'rosybrown':            '#BC8F8F',
    'royalblue':            '#4169E1',
    'saddlebrown':          '#8B4513',
    'salmon':               '#FA8072',
    'sandybrown':           '#FAA460',
    'seagreen':             '#2E8B57',
    'seashell':             '#FFF5EE',
    'sienna':               '#A0522D',
    'silver':               '#C0C0C0',
    'skyblue':              '#87CEEB',
    'slateblue':            '#6A5ACD',
    'slategray':            '#708090',
    'snow':                 '#FFFAFA',
    'springgreen':          '#00FF7F',
    'steelblue':            '#4682B4',
    'tan':                  '#D2B48C',
    'teal':                 '#008080',
    'thistle':              '#D8BFD8',
    'tomato':               '#FF6347',
    'turquoise':            '#40E0D0',
    'violet':               '#EE82EE',
    'wheat':                '#F5DEB3',
    'white':                '#FFFFFF',
    'whitesmoke':           '#F5F5F5',
    'yellow':               '#FFFF00',
    'yellowgreen':          '#9ACD32'}

    return cnames

#-------------------------------------------------------------------------------------------------------------
def plot_linecoverage(lines,filters,outname='testfigure__RENAME__.pdf', figuresize_x=8, Fsize=10,
                      redshiftrange=[0.,17.],zspacing=1.0,markfilter=None,verticalmarker=None,verbose=True):
    """

    Plotting an overview of coverage of emission lines for given filters as a function of redshift

    --- INPUT ---
    lines           list of lines in MiGs.linelistdic to show in plot
    filters         list of filter names to illustrate line location for in plot
    outname         name of output plot to generate
    figuresize_x    Size of figure in x direction
    Fsize           Text font size
    redshiftrange   Range of redshift to include in plot
    zspacing        Spacing between redshift axis ticks
    markfilter      List of filters to outline in black
    verticalmarker  Add vertical lines at given redshifts providing [[redshift1,name1,color1],[redshift2,name2,color2],...]
    verbose         toggle verbosity

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    lines   = ['lya','civ2','ciii2','oii2','oiii2','ha','pah6p2','h29p7','neii12p8']
    filters = ['muse','wfc3_g102','wfc3_g141','niriss_f115w','niriss_f150w','niriss_f200w','nirspec_g140h/F100lp','nirspec_g235h/F170lp','nirspec_g395h/F290lp','nirspec_prismclear','nircam_f277w','nircam_f356w','nircam_f444w','miri_lrs','miri_mrs']
    outname = '/Users/kschmidt/work/linecoverageplot.pdf'
    kbs.plot_linecoverage(lines,filters,outname=outname,verbose=True)

    lines   = ['lya','oii2','hb','oiii2','ha','paa','brg']
    filters = ['muse','wfc3_g102','wfc3_g141','niriss_f115w','niriss_f150w','niriss_f200w','nircam_f277w','nircam_f356w','nircam_f444w']
    markfilter = ['nircam_f277w','nircam_f356w','nircam_f444w']
    kbs.plot_linecoverage(lines,filters,outname=outname,verbose=True,verticalmarker=[[0.308,'$z$(A2744) = 0.308','black'],[4.5,'$z$(MUSE LAE) = 4.5','gray'],[8.38,'$z$(YD4) = 8.38','lightgray']],markfilter=markfilter,redshiftrange=[0.0,13.0],figuresize_x=12,Fsize=15)

    """
    filterinfo = kbs.get_filterinfo()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading line list dictionary from MiGs.linelistdic'
    linedic      = MiGs.linelistdic(listversion='full') # loading line list for plots

    import matplotlib
    cmap = matplotlib.cm.get_cmap('rainbow')
    norm = matplotlib.colors.Normalize(vmin=0, vmax=len(filters))

    filtercolors = [ cmap(norm(filternumber)) for filternumber in np.arange(len(filters))]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Plotting figure '
    figuresize_y = len(lines) #6
    fig, ax      = plt.subplots(figsize=(figuresize_x,figuresize_y))
    plt.rc('text', usetex=True)                         # enabling LaTex rendering of text
    plt.rc('font', family='serif',size=Fsize)           # setting text font
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    left   = 0.12   # the left side of the subplots of the figure
    right  = 0.72   # the right side of the subplots of the figure
    bottom = 1.0/figuresize_y       # the bottom of the subplots of the figure
    top    = 1.0-1.0/figuresize_y   # the top of the subplots of the figure
    wspace = 0.20   # the amount of width reserved for blank space between subplots
    hspace = 0.20   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    plt.gca().xaxis.grid(b=True, which='major',linestyle='--',zorder=0,color='lightgray')
    # plt.minorticks_on()
    # plt.gca().xaxis.grid(b=True, which='minor',linestyle=':',zorder=0,color='gray')
    plt.gca().yaxis.grid(False)

    xlabel  = 'Redshift'
    ylabel  = ''
    yrange  = [0.5,len(lines)+0.5]

    linelabels = []
    dy         = 1.0/len(filters)
    for ll, line in enumerate(lines):
        linelabels.append(linedic[line][0])
        for ff,filt in enumerate(filters):
            if filt not in filterinfo.keys():
                sys.exit(' Filter info for "'+filt+'" was not found in filterinfo dictionary from kbs.get_filterinfo')
            label  = filterinfo[filt]['label']
            zrange = np.asarray(filterinfo[filt]['waverange'])/linedic[line][1]-1.0
            yhigh   = [ll+1.5-ff*dy]*2
            ylow  = [ll+1.5-(ff+1)*dy]*2

            if markfilter is not None:
                if filt in markfilter:
                    ecolor="black"
                else:
                    ecolor=filtercolors[ff]
            else:
                ecolor=filtercolors[ff]

            if ll == 0:
                plt.fill_between(zrange,yhigh,ylow,alpha=0.7,facecolor=filtercolors[ff],label=label,zorder=20,edgecolor=ecolor)
            else:
                plt.fill_between(zrange,yhigh,ylow,alpha=0.7,facecolor=filtercolors[ff],zorder=20,edgecolor=ecolor)

        if ll < len(lines)-1:
            plt.plot(redshiftrange,[ll+1.5]*2,'k-',linewidth=1.0,zorder=30)

    ylabelposition = np.arange(len(lines))+1
    plt.yticks(ylabelposition, linelabels)

    if verticalmarker is not None:
        for vm in verticalmarker:
            plt.plot([vm[0],vm[0]],yrange,'-',color=vm[2],label=vm[1],zorder=100)

    redshiftticks        = np.arange(redshiftrange[0],redshiftrange[1]+1, zspacing)
    if zspacing == np.round(zspacing):
        redshiftlabels   = redshiftticks.astype(int).astype(str)
        plt.xticks(redshiftticks,redshiftlabels)
    else:
        redshiftlabels   = redshiftticks.astype(float).astype(str)
        plt.xticks(redshiftticks,redshiftlabels,rotation=90)
    # redshiftlabels[1::2] = ''

    #plt.yticks('')
    #plt.ylabels(linelabels)

    #--------- LEGEND ---------
    anchorpos = (1.0, 0.99)
    leg = plt.legend(fancybox=True,numpoints=1, loc='upper left',prop={'size':Fsize},ncol=1,
                     bbox_to_anchor=anchorpos)  # add the legend
    leg.get_frame().set_alpha(0.7)
    #--------------------------

    plt.xlim(redshiftrange)
    plt.ylim(yrange)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    ax2 = ax.twiny()  # ax2 is responsible for "top" axis and "right" axis
    ax2.set_xticks(redshiftticks)
    if zspacing == np.round(zspacing):
        ax2.set_xticklabels(redshiftlabels)
    else:
        ax2.set_xticklabels(redshiftlabels,rotation=90)
    # ax2.axis["right"].major_ticklabels.set_visible(False)
    # ax2.axis["top"].major_ticklabels.set_visible(True)
    ax2.set_xlabel(xlabel)
    ax2.set_xlim(redshiftrange)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.savefig(outname)
    plt.clf()
    plt.close('all')
    if verbose: print ' - Saved figure to ',outname


#-------------------------------------------------------------------------------------------------------------
def get_filterinfo():
    """
    band widths taken from http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=HST&gname2=ACS_WFC

    """
    filterinfo = collections.OrderedDict()

    filterinfo['bluemuse']        = {'label':'Blue-MUSE',        'waverange':[3700,6000]}
    filterinfo['muse']            = {'label':'MUSE',             'waverange':[4800,9300]}

    filterinfo['mosfire']         = {'label':'MOSFIRE',          'waverange':[9000,25000]}

    filterinfo['wfc3_g102']       = {'label':'HST WFC3 G102',        'waverange':[8000,11719]}
    filterinfo['wfc3_g141']       = {'label':'HST WFC3 G141',        'waverange':[10404,17747]}

    # NIRISS
    filterinfo['niriss_f090w']    = {'label':'NIRISS F090W',     'waverange':[9973,13061]}
    filterinfo['niriss_f115w']    = {'label':'NIRISS F115W',     'waverange':[9973,13061]}
    filterinfo['niriss_f150w']    = {'label':'NIRISS F150W',     'waverange':[13020,16968]}
    filterinfo['niriss_f200w']    = {'label':'NIRISS F200W',     'waverange':[17201,22595]}

    # NIRCAM
    filterinfo['nircam_f277w']    = {'label':'NIRCAM F277W',     'waverange':[23657,32166]}
    filterinfo['nircam_f356w']    = {'label':'NIRCAM F356W',     'waverange':[30696,40787]}
    filterinfo['nircam_f444w']    = {'label':'NIRCAM F444W',     'waverange':[38013,50996]}

    # NIRSpec info from: https://jwst-docs.stsci.edu/display/JTI/NIRSpec+Dispersers+and+Filters
    filterinfo['nirspec_g140mf070lp']  = {'label':'NIRSpec G140M/F070LP',    'waverange':[7000,12700]}
    filterinfo['nirspec_g140m/F100lp'] = {'label':'NIRSpec G140M/F100LP',    'waverange':[9700,18900]}
    filterinfo['nirspec_g235m/F170lp'] = {'label':'NIRSpec G235M/F170LP',    'waverange':[16600,31700]}
    filterinfo['nirspec_g395m/F290lp'] = {'label':'NIRSpec G395M/F290LP',    'waverange':[28700,52700]}
    filterinfo['nirspec_g140h/F070lp'] = {'label':'NIRSpec G140H/F070LP',    'waverange':[7000,12700]}
    filterinfo['nirspec_g140h/F100lp'] = {'label':'NIRSpec G140H/F100LP',    'waverange':[9700,18900]}
    filterinfo['nirspec_g235h/F170lp'] = {'label':'NIRSpec G235H/F170LP',    'waverange':[16600,31700]}
    filterinfo['nirspec_g395h/F290lp'] = {'label':'NIRSpec G395H/F290LP',    'waverange':[28700,52700]}
    filterinfo['nirspec_prismclear']   = {'label':'NIRSpec PRISM/CLEAR',     'waverange':[6000,53000]}

    # MIRI
    filterinfo['miri_F560W']   = {'label':'MIRI F560W',     'waverange':[48847,64243]}
    filterinfo['miri_F770W']   = {'label':'MIRI F770W',     'waverange':[64740,88353]}
    filterinfo['miri_F1000W']  = {'label':'MIRI F1000W',    'waverange':[87628,111018]}
    filterinfo['miri_F1130W']  = {'label':'MIRI F1130W',    'waverange':[106390,119837]}
    filterinfo['miri_F1280W']  = {'label':'MIRI F1280W',    'waverange':[112516,143369]}
    filterinfo['miri_F1500W']  = {'label':'MIRI F1500W',    'waverange':[131189,171411]}
    filterinfo['miri_F1800W']  = {'label':'MIRI F1800W',    'waverange':[160279,202851]}
    filterinfo['miri_F2100W']  = {'label':'MIRI F2100W',    'waverange':[178759,244051]}
    filterinfo['miri_F2550W']  = {'label':'MIRI F2550W',    'waverange':[223391,299940]}

    # MIRI LRS https://jwst-docs.stsci.edu/display/JTI/MIRI+Low-Resolution+Spectroscopy
    filterinfo['miri_lrs']   = {'label':'MIRI LRS',     'waverange':[50000,120000]}

    # MIRI MRS https://jwst-docs.stsci.edu/display/JTI/MIRI+Medium-Resolution+Spectroscopy
    filterinfo['miri_mrs']   = {'label':'MIRI MRS IFU',       'waverange':[49000,288000]}

    return filterinfo
#-------------------------------------------------------------------------------------------------------------
def combine_fits_tables(fitstables,outputtable,genDS9reg=True,racol='RA',deccol='DEC',idcol=None):
    """
    combine a set of fits tables into a new fits table.
    Assumes that input table contain the same columns.

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    import glob

    tabdir    = '/Users/kschmidt/work/MUSE/MWDR1galfitmodeling/galfit_wrapper_results_final/'
    magtables = glob.glob(tabdir+'galfit_wrapper*cdfs*775w*.fits')
    outtable  = tabdir+'galfit_wrapper_output_cat_acs_775w_DR1objWOskeltonID.fits'
    kbs.combine_fits_tables(magtables,outtable,genDS9reg=True,idcol='UNIQUE_ID')

    """
    if os.path.isfile(outputtable):
        sys.exit(' The output already exists - not overwriting the file \n'+outputtable)
    else:
        print(' - Will attempt to combine (append) the '+str(len(fitstables))+' fits table provided ')
    newtable = atab.Table.read(fitstables[0], format='fits')
    for fitstab in fitstables[1:]:
        addtab   = atab.Table.read(fitstab, format='fits')
        newtable = atab.vstack([newtable, addtab])

    newtable.write(outputtable)
    print(' - Wrote new table to \n   '+outputtable)

    if genDS9reg:
        if idcol is None:
            textlist = None
        else:
            textlist = newtable[idcol]
        kbs.create_DS9region(outputtable.replace('.fits','.reg'),newtable[racol],newtable[deccol],
                             color='red',circlesize=0.5,textlist=textlist,clobber=False)
#-------------------------------------------------------------------------------------------------------------
def odr_fitfunction_linear(p, x):
    a, b = p
    return a*x + b
def odr_fitfunction_quad(p, x):
    a, b, c = p
    return a * x*x + b*x + c
def fit_function_to_data_with_errors_on_both_axes(xval,yval,xerr,yerr,fitfunction='linear',plotresults='./ODRfit2data.pdf'):
    """
    Use scipy's Orthogonal Distance Regression (ODR) to fit a function to data with errors on both axes.
    Essential a 2D leas squares, where the diagonal convariance matrices of the data is accounted for.

    Example taken from https://micropore.wordpress.com/2017/02/07/python-fit-with-error-on-both-axis/

    --- INPUT ---


    --- EXAMPLE OF USE ---


    """
    print(' - Setting up model to fit ('+fitfunction+')')
    if fitfunction == 'linear':
        quad_model = odr.Model(kbs.odr_fitfunction_linear)
        initguess  = [0., 1.]
    elif fitfunction == 'quadratic':
        quad_model = odr.Model(kbs.odr_fitfunction_quad)
        initguess  = [0., 1., 1.]
    else:
        sys.exit(' - fitfunction='+str(fitfunction)+' is not a valid choice')

    print(' - Making sure all inputs are float arrays ')
    xval = np.asarray(xval).astype(float)
    yval = np.asarray(yval).astype(float)
    xerr = np.asarray(xerr).astype(float)
    yerr = np.asarray(yerr).astype(float)

    print(' - Create a ReadData object from the provided data ')
    data = odr.RealData(xval, yval, sx=xerr, sy=yerr)

    print(' - Set up ODR with the mode and the data ')
    odrclass = odr.ODR(data, quad_model, beta0=initguess)

    print(' - Run the Orthogonal Distance Regression on the data ')
    fitresults = odrclass.run()

    print(' - The 1sigma parameter estimates of the '+fitfunction+' from the ODR fit are ')
    popt = fitresults.beta
    perr = fitresults.sd_beta
    for i in range(len(popt)):
        print('   '+str(popt[i])+' +/- '+str(perr[i]))

    if plotresults:
        nstd = 5. # to draw 5-sigma intervals
        popt_up = popt + nstd * perr
        popt_dw = popt - nstd * perr

        x_fit = np.linspace(min(xval), max(xval), 100)
        if fitfunction == 'linear':
            fit = kbs.odr_fitfunction_linear(popt, x_fit)
            fit_up = kbs.odr_fitfunction_linear(popt_up, x_fit)
            fit_dw= kbs.odr_fitfunction_linear(popt_dw, x_fit)
        elif fitfunction == 'quadratic':
            fit = kbs.odr_fitfunction_quad(popt, x_fit)
            fit_up = kbs.odr_fitfunction_quad(popt_up, x_fit)
            fit_dw= kbs.odr_fitfunction_quad(popt_dw, x_fit)

        fig = plt.figure(figsize=(5,5))
        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
        lthick = 1
        Fsize  = 12
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        plt.title('Orthogonal Distance Regression fit to data w. error on both axes', fontsize=Fsize)

        plt.errorbar(xval, yval, yerr=yerr, xerr=xerr, hold=True, ecolor='k', fmt='none', label='data')
        plt.fill_between(x_fit, fit_up, fit_dw, alpha=0.25, label='5-sigma interval',color='red')
        plt.plot(x_fit, fit, 'r', lw=2, label='Best fit curve')

        plt.xlabel('x values', fontsize=Fsize)
        plt.ylabel('y values', fontsize=Fsize)

        leg = plt.legend(fancybox=True, loc='lower right',prop={'size':Fsize},ncol=1,numpoints=1)
        leg.get_frame().set_alpha(0.7)

        plt.savefig(plotresults)
    return fitresults
#-------------------------------------------------------------------------------------------------------------
def UniverseEpochsInFractionsOfHumaLifetime(humanlifetime=80,zevaluate=["MACS J1149",0.54],H0=70, Omega_m0=0.3, T_CMB0=2.725):
    """
    Function printing important times and epochs in the evolution of the universe in terms of a human lifetime.

    --- INPUT ___
    humanlifetime    Assumed human lifetime in years
    universeage      Total age of universe in Gyr

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    kbs.UniverseEpochsInFractionsOfHumaLifetime(humanlifetime=80)

    """
    print(' - Setting up cosmological calculator for input H0='+str(H0)+', Omega_m0='+str(Omega_m0)+
          ' and T_CMB0='+str(T_CMB0))
    cosmo               = FlatLambdaCDM(H0=H0, Om0=Omega_m0, Tcmb0=T_CMB0)
    TotalAgeUniverse    = cosmo.age(0.0)
    print('   This results in an Universe age of '+str(TotalAgeUniverse))

    print(' - Defining parameters')
    secINmin            = 60.0              # sec
    secINhour           = secINmin * 60.0   # min/hour
    secINday            = secINhour * 24.0  # hour/day
    secINmonth          = secINday * 30.0   # days/month
    secINyear           = secINday * 365.25 # day/year
    secINhumanlifetime  = humanlifetime * secINyear
    yearINuniverse      = TotalAgeUniverse.value * 10e9

    conversionfactor    =  yearINuniverse/secINhumanlifetime # years of Universe / sec in human lifetime

    print(' - Defining epochs (redshifts) to evaluate:')
    epochs = collections.OrderedDict()
    epochs['Last Scattering (CMB)']     = 1100.0
    epochs['Dark Ages']                 = 100.0
    epochs['EoR start']                 = 15.0
    epochs['EoR "mid"']                 = 10.0
    epochs['EoR end']                   = 6.0
    epochs['SF heyday']                 = 2.0
    epochs[zevaluate[0]]                = zevaluate[1]

    for key in epochs.keys():
        redshift = epochs[key]
        Uage     = cosmo.age(redshift)
        humansec = (Uage.value*1e9)/conversionfactor

        print('\n    - '+str("%.30s" % key)+' @ z = '+str("%10.4f" % redshift)+' ~ Age = '+str(Uage))
        print('       For a human lifetime of '+str(humanlifetime)+' years this corresponds to:')
        print('       human minutes = '+str("%10.4f" % (humansec/secINmin)))
        print('       human hours   = '+str("%10.4f" % (humansec/secINhour)))
        print('       human days    = '+str("%10.4f" % (humansec/secINday)))
        print('       human months  = '+str("%10.4f" % (humansec/secINmonth)))
        print('       human years   = '+str("%10.4f" % (humansec/secINyear)))

#-------------------------------------------------------------------------------------------------------------
def reproject_fitsimage(fitsimage,fitsoutput,fitsimageext=0,radec=None,naxis=None,overwrite=False,verbose=True):
    """
    Reprojecting, i.e., rotating fits image around given position so North is up and East is left.

    --- INPUT ---
    fitsimage         FITS file containing image to reporject in extension "fitsimageext".
    fitsoutput        Path an name of file to store reprojects (rotated) image to.
    fitsimageext      FITS extension of image to reproject.
    radec             R.A. and Dec. [deg] of point in image to use as reference pixel in new image.
                      If None the central pixel of the image will be used.
    naxis             The output image size [x,y] to store from reprojected image. By default the output will
                      be a square image of size max(NAXIS1,NAXIS2) of fitsimage.
    overwrite         Overwrite existing output image?
    verbose           Toggle verbosity.

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    imgpath     = '/Users/kschmidt/work/MUSE/MUSEGalaxyGroups/CGR84_acs/acs_2.0_cutouts/'
    fitsimage   = imgpath+'0001_150.05075000_2.59641000_acs_I_100006+0235_unrot_sci_20.fits'

    fitsoutput  = fitsimage.replace('unrot','rot')
    kbs.reproject_fitsimage(fitsimage,fitsoutput,fitsimageext=0,overwrite=False)

    fitsoutput  = fitsimage.replace('unrot','rotandcut')
    radec = [150.04876,2.5994931]
    naxis = [1000,1500]
    kbs.reproject_fitsimage(fitsimage,fitsoutput,radec=radec,naxis=naxis,overwrite=False)

    """
    if verbose:
        print(' - Will reproject (rotate) fits image:')
        print('   '+fitsimage+' (extension = '+str(fitsimageext)+')')

    imghdu     = afits.open(fitsimage)[fitsimageext]
    hdrin      = imghdu.header
    wcsin      = wcs.WCS(hdrin)
    pixscale   = wcs.utils.proj_plane_pixel_scales(wcsin)
    hdrout     = hdrin.copy()

    if radec is None:
        pixpos     = np.asarray(imghdu.data.shape)/2
        radec      = wcsin.wcs_pix2world(pixpos[1],pixpos[0],1)
        if verbose: print(' - Using central coordinate set (ra,dec) = ('+str(radec[0])+','+str(radec[1])+
                          ') as reference pixel in new image ')
    else:
        if verbose: print(' - Using input coordinate set (ra,dec) = ('+str(radec[0])+','+str(radec[1])+
                          ') as reference pixel in new image ')
        pixpos     = wcsin.wcs_world2pix(radec[0],radec[1],1)

    hdrout['NAXIS1'] = np.max([hdrin['NAXIS1'],hdrin['NAXIS2']])
    hdrout['NAXIS2'] = np.max([hdrin['NAXIS1'],hdrin['NAXIS2']])
    hdrout['CRPIX1'] = float(pixpos[0])
    hdrout['CRPIX2'] = float(pixpos[1])
    hdrout['CRVAL1'] = float(radec[0])
    hdrout['CRVAL2'] = float(radec[1])
    hdrout['CD1_1']  = -1*np.abs(pixscale[0])
    hdrout['CD1_2']  = 0.0
    hdrout['CD2_2']  = np.abs(pixscale[0])
    hdrout['CD2_1']  = 0.0

    array, footprint = reproject_interp(imghdu, hdrout)

    if naxis is not None:
        if verbose: print(' - Cutting out projected image array to store around reference coordinate')
        afits.writeto(fitsoutput, array, hdrout, overwrite=overwrite)
        imghdu     = afits.open(fitsoutput)[0]
        hdrin      = imghdu.header
        # array = array[int(hdrout['CRPIX2']-hdrout['NAXIS2']/2):int(hdrout['CRPIX2']+hdrout['NAXIS2']/2),
        #               int(hdrout['CRPIX1']-hdrout['NAXIS1']/2):int(hdrout['CRPIX1']+hdrout['NAXIS1']/2)]

        hdrout = hdrin.copy()
        hdrout['NAXIS1'] = naxis[0]
        hdrout['NAXIS2'] = naxis[1]
        hdrout['CRPIX1'] = naxis[0]/2.
        hdrout['CRPIX2'] = naxis[1]/2.
        array, footprint = reproject_interp(imghdu, hdrout)

    afits.writeto(fitsoutput, array, hdrout, overwrite=overwrite)
#-------------------------------------------------------------------------------------------------------------
def plot_GALFITmodel(galfit_imgblock,colormap='viridis',showcomponentnumbers=False,
                     figsize=(10,3),logcolor=True,vscale=0.99,verbose=True,addcircles=None):
    """
    Plotting the FITS extensions of a GALFIT imgblock
    Plotting based on plot_fitsimage, but using sub-plots to get overview.

    --- INPUT ___
    galfit_imgblock       The GALFIT imgblock to plot
    colormap              Colormap to use for plotting
    showcomponentnumbers  To add component numbers on model panel, set to True
    figsize               Size of figure. Labels will stay fixed so small size = large labels.
    logcolor              If true, color map will be in log
    vscale                The scale to apply when estimating vmin and vmax. If type(vscale) is not float,
                          it is assumed that list of vmin and vmax is provided instead.
    verbose               Toggle verbosity
    addcircles            To add circles to individual frames provide list of [x,y,r,color]


    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    path       = '/Users/kschmidt/work/MUSE/MUSEGalaxyGroups/GALFIT/galfit_wrapper_results_final/imgblocks/'

    galfitimg  = path+'imgblock_CGr32_289_acs_814w_289.0_190307.fits'
    kbs.plot_GALFITmodel(galfitimg,colormap='nipy_spectral',vscale=[1e-4,1.0],logcolor=True)
    """
    fileext    = galfit_imgblock.split('.')[-1]
    outputfile = galfit_imgblock.replace('.'+fileext,'_overview.pdf')
    if verbose: print(' - Will plot GALFIT model summary to:\n   '+outputfile)
    if verbose: print(' - Using    logartihmic colormap : '+str(logcolor))
    if verbose: print('            vscale               : '+str(vscale))

    dat_refimg   = afits.open(galfit_imgblock)[1].data
    dat_model    = afits.open(galfit_imgblock)[2].data
    dat_residual = afits.open(galfit_imgblock)[3].data

    if type(vscale) is not float:
        vmin = vscale[0]
        vmax = vscale[1]
    else:
        vmin,vmax  = kbs.get_vminvmax(dat_refimg, logscale=logcolor, scale=vscale, verbose=False)
    if verbose: print('            i.e. [vmin,vmax]     : ['+str(vmin)+','+str(vmax)+']')

    if logcolor:
        normclass = LogNorm()
    else:
        normclass = None

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.4, hspace=0.1,left=0.05, right=0.95, bottom=0.05, top=0.98)
    Fsize  = 12
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    Nrow = 1
    Ncol = 3
    for ii, imgdat in enumerate([dat_refimg,dat_model,dat_residual]):
        ax = plt.subplot2grid((Nrow,Ncol), (0, ii), colspan=1, rowspan=1)
        im = plt.imshow(imgdat, cmap=colormap, norm=normclass, origin='lower',vmin=vmin, vmax=vmax)
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        ax.grid(False)

        if (ii == 1) & showcomponentnumbers: # Adding component numbers to model panel
            hdrinfo = afits.open(galfit_imgblock)[2].header

            Ncomp   = 0
            skycomp = False
            for key in hdrinfo.keys():
                if 'COMP_' in key:
                    if 'sky' not in hdrinfo[key]:
                        Ncomp = Ncomp+1
                    else:
                        skycomp = True
            infostr = 'GALFIT components = '+str(Ncomp)
            if skycomp:
                infostr = infostr+' (+sky)'
            ax.text(0,0,infostr,fontsize=Fsize-3, color='black', ha='left', va='bottom')

            for cc in np.arange(Ncomp):
                XCinfo = hdrinfo[str(cc+1)+'_XC'].replace('[','').replace(']','')
                YCinfo = hdrinfo[str(cc+1)+'_YC'].replace('[','').replace(']','')
                xpix = float(XCinfo.split()[0])
                ypix = float(YCinfo.split()[0])
                ax.text(xpix,ypix,str(cc+1),fontsize=Fsize, color='black', ha='center', va='center')

        if addcircles is not None:
            for cirsetup in addcircles:
                circ = Circle((cirsetup[0],cirsetup[1]),cirsetup[2],color=cirsetup[3],
                              alpha=0.8,linestyle='-',linewidth=3, fill=False, label=None)
                ax.add_patch(circ)


    plt.savefig(outputfile)
    fig.clf()

#-------------------------------------------------------------------------------------------------------------
def plot_fitsimage(fitsfile,outputfile,fitsext=0,colormap='viridis',figsize=(4,3),addcircles=None,
                   logcolor=True,vscale=0.99,verbose=True):
    """

    --- INPUT
    fitsfile       FITS file containing image to plot
    outputfile     Name of image to generate
    fitsext        Extension of FITS file containing image to plot
    colormap       Colormap to use for plotting
    figsize        Size of figure. Labels will stay fixed so small size = large labels.
    addcircles     To add circles to image, provide list of [x,y,r,color] for the circles to add.
    logcolor       If true, color map will be in log
    vscale         The scale to apply when estimating vmin and vmax. If type(vscale) is not float,
                   it is assumed that list of vmin and vmax is provided instead.
    verbose        Toggle verbosity

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    path       = '/Users/kschmidt/work/MUSE/MUSEGalaxyGroups/GALFIT/galfit_wrapper_results_final/imgblocks/'
    fitsfile   = path+'imgblock_CGr32_289_acs_814w_289.0_190307.fits'
    outputfile = fitsfile.replace('.fits','_model.pdf')
    kbs.plot_fitsimage(fitsfile,outputfile,fitsext=2,colormap='nipy_spectral',vscale=0.99,logcolor=True)

    """
    if verbose: print(' - Will plot image to:\n   '+outputfile)
    if verbose: print(' - Using    logartihmic colormap : '+str(logcolor))
    if verbose: print('            vscale               : '+str(vscale))
    imgdat = afits.open(fitsfile)[fitsext].data

    if type(vscale) is not float:
        vmin = vscale[0]
        vmax = vscale[1]
    else:
        vmin,vmax  = kbs.get_vminvmax(imgdat, logscale=logcolor, scale=vscale, verbose=False)
    if verbose: print('            i.e. [vmin,vmax]     : ['+str(vmin)+','+str(vmax)+']')

    if logcolor:
        normclass = LogNorm()
    else:
        normclass = None

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.90, bottom=0.1, top=0.95)
    Fsize  = 12
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    ax = plt.gca()
    im = plt.imshow(imgdat, cmap=colormap, norm=normclass, origin='lower',vmin=vmin, vmax=vmax)

    if addcircles is not None:
        for cirsetup in addcircles:
            circ = Circle((cirsetup[0],cirsetup[1]),cirsetup[2],color=cirsetup[3],
                          alpha=0.8,linestyle='-',linewidth=3, fill=False, label=None)
            ax.add_patch(circ)

    divider = make_axes_locatable(ax)
    cax     = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.grid(False)
    plt.savefig(outputfile)
    fig.clf()

#-------------------------------------------------------------------------------------------------------------
def get_vminvmax(dataarray, scale=0.95, logscale=False, verbose=True):
    """
    Define vmin and vmax for an input image

    --- INPUT ---
    dataarray      Data array to estimate vmin and vmax for
    scale          The scale factor ]0,1] to use for the vmin and vmax values
    logscale       If logscal is True, negative calues will be ignored in array when setting vmin and vmax
    verbose        Toggle verbosity

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    vmin,vmax = kbs.get_vminvmax(data, scale=0.99)

    """
    if scale > 1 or scale <= 0:
        sys.exit(' Scale should be contained in ]0,1]; it is not as it is scale='+str(scale))
    else:
        if verbose: print(' - Will estimate vmin/vmax using a scale of '+str(scale))

    includemask = np.isfinite(dataarray)

    if logscale:
        positivemask = (dataarray > 0)
        includemask  = (includemask & positivemask)

    flatdata  = dataarray[includemask].flatten()

    if len(flatdata) > 1:
        sortdata   = np.sort(flatdata)
        length     = len(sortdata)
        vmax_index = int(np.floor(scale*length))-1
        vmin_index = int(np.ceil((1-scale)*length))
        vmin       = sortdata[vmin_index]
        vmax       = sortdata[vmax_index]
    else:
        if verbose: print(' WARNING: No finite pixel values in data array')
        vmin = np.min(dataarray)
        vmax = np.max(dataarray)

    if verbose: print(' - Returning vmin,vmax = '+str(vmin)+','+str(vmax))
    return vmin, vmax

#-------------------------------------------------------------------------------------------------------------
def align_arrays(fix_image_arr, reg_image_arr, reg_image_arr_noise=None, verbose=True):
    """
    Aligning two (image) arrays by (re-)registering the array content using the
    image_registration package

    --- INPUT ---
    fix_img_arr         The image array to keep fixed registering (aligning) the reg_img_arr to.
    reg_img_arr         The image array to register to fix_img_arr
    reg_img_arr_noise   Per pixel noise/uncertainty in the reg_img_arr
    verbose             Toggle verbosity

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    arr_shifts, arr_aligned = kbs.align_arrays(fix_image_arr, reg_image_arr)

    The output contains:    arr_shifts = xoff, yoff, xoff_err, yoff_err  ; the shifts to apply for alignment
                            arr_aligned                                  ; shifts applied to reg_image_arr

    """
    if verbose: print(' - Determining shifts for optimal (minimum chi2) array alignment.')
    xoff, yoff, xoff_err, yoff_err   = chi2_shift(fix_image_arr, reg_image_arr, reg_image_arr_noise,
                                                  return_error=True, upsample_factor='auto')
    if verbose: print('   Shifts were determined to be (xoff, yoff, xoff_err, yoff_err) = '+
                      str([xoff, yoff, xoff_err, yoff_err]))

    if verbose: print(' - Applying shifts to reg_img_arr to produced re-registerd (aligned) array')
    arr_aligned  = shift.shiftnd(reg_image_arr, (-yoff, -xoff))

    if verbose: print(' - Returning shifts and shifted array')
    arr_shifts = xoff, yoff, xoff_err, yoff_err
    return arr_shifts, arr_aligned # arr_shifts contains xoff, yoff, xoff_err, yoff_err

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def calc_concentration(fitsimage,verbose=True):
    """
    Calculating the concentration (C in CAS) of an image as defined in Conselice (2013) section 4.1.

    --- INPUT ---
    fitsimage      image to calculcate clumpiness off

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu

    """
    sys.exit(' Not written yet')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def calc_asymmetry(fitsimage,verbose=True):
    """
    Calculating the asymmetry (A in CAS) of an image as defined in Conselice (2013) section 4.2.

    --- INPUT ---
    fitsimage      image to calculcate clumpiness off

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu

    """
    sys.exit(' Not written yet')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def calc_clumpiness(fitsimage,sigmasmooth,imgext=0,bckimg=None,bckimgext=0,negative2zero=True, mask=None,
                    r_center=None, plotimages=None, verbose=True):
    """
    Calculating the clumpiness (S in CAS) of an image as defined in Conselice (2013) section 4.3 - well, actually
    using definition in Eq. 5 of Christopher J. Conselice, Cui Yang and Asa F. L. Bluck (2009).

    --- INPUT ---
    fitsimage      image to calculcate clumpiness off
    sigmasmooth    the width of the Gaussian smoothing kernel used
    imgext         FITS extecions of image
    bckimg         If a background region is defined, provide that as an image
    bckimgext      FITS extension of background image
    negative2zero  Set all negative valus in fits image to 0? Done by Conselic in original paper.
    segmap         To apply a mask to the fitsimage (post smoothing) and the skyimg supply it here.
                   This mask could be resulting from a SExtractor segmentation map.
    r_center       Provide a radius in pixels to mark central region to include in estimate. The remaining
                   pixels will be set to 0, hence, corresponding to a central circular mask applied to image.
    showimage      Provde a path and file name to plot the image used to estimate clumpiness after masking etc.
    verbose        Toggle verbosity

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs

    imgdir       = '/Volumes/DATABCKUP1/TDOSEextractions/190325_MWDR1_OIIemitters_apertureExt/tdose_cutouts/'
    fitsimage    = imgdir+'acs_814w_candels-cdfs-46_cut_v1.0_id146074360_cutout4p0x4p0arcsec.fits'
    sigmasmooth  = 5.0 #pixels
    Sval         = kbs.calc_clumpiness(fitsimage,sigmasmooth,negative2zero=True,r_center=False,plotimages=imgdir+'clumpinessimages.pdf')
    """
    if verbose: print(' - Estimating clumpiness as defined by '
                      'C. J. Conselice, C. Yang and A. F. L. Bluck (2009).')
    dataimg   = afits.open(fitsimage)[imgext].data
    if negative2zero:
        dataimg[dataimg < 0] = 0

    if mask is not None:
        dataimg[mask] = 0.0

    dataimg_smooth = gaussian_filter(dataimg, sigma=sigmasmooth)

    if r_center: # generate circular mask and apply to image
        imgxsize, imgysize = dataimg.shape
        y,x  = np.ogrid[-imgysize/2.0:imgysize/2.0, -imgxsize/2.0:imgxsize/2.0]
        circularmask = x*x + y*y <= r_center*r_center
        dataimg[~circularmask] = 0.0 # setting pixels outside circle to 0
        dataimg_smooth[~circularmask] = 0.0

    if bckimg is not None:
        Barr = afits.open(bckimg)[bckimgext].data
        Barr[dataimg == 0] = 0.0
        Barr_smooth    = gaussian_filter(Barr, sigma=sigmasmooth)
    else:
        Barr        = dataimg * 0.0
        Barr_smooth = dataimg * 0.0

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plotimages is not None:
        plotlog = True
        vmin,vmax  = kbs.get_vminvmax(dataimg, logscale=plotlog, scale=0.95, verbose=False)
        if verbose: print(' - Will plot using [vmin,vmax]     : ['+str(vmin)+','+str(vmax)+']')
        if plotlog:
            normclass = LogNorm()
        else:
            normclass = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        fig = plt.figure(figsize=(5,5))
        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.90, bottom=0.1, top=0.95)
        Fsize  = 12
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        ax = plt.gca()
        im = plt.imshow(dataimg, cmap='viridis', norm=normclass, origin='lower',vmin=vmin, vmax=vmax)
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        ax.grid(False)
        plt.savefig(plotimages)
        fig.clf()
        plt.close()
        if verbose: print(' - saved plot to '+plotimages)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        fig = plt.figure(figsize=(5,5))
        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.90, bottom=0.1, top=0.95)
        Fsize  = 12
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        ax = plt.gca()
        im = plt.imshow(dataimg_smooth, cmap='viridis', norm=normclass, origin='lower',vmin=vmin, vmax=vmax)
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        ax.grid(False)
        plt.savefig(plotimages.replace('.pdf','_smoothed.pdf'))
        fig.clf()
        plt.close()
        if verbose: print(' - saved plot to '+plotimages.replace('.pdf','_smoothed.pdf'))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    imgpart = np.sum(dataimg-dataimg_smooth) / np.sum(dataimg)
    bckpart = np.sum(Barr-Barr_smooth)       / np.sum(dataimg)
    Sval    = 10.0 * (imgpart - bckpart)
    if verbose: print('   Clumpiness values estimated to be '+str(Sval))
    return Sval
#-------------------------------------------------------------------------------------------------------------
#                                                  END
#-------------------------------------------------------------------------------------------------------------
