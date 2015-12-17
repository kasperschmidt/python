"""
#----------------------------
#   NAME
#----------------------------
# calibrate1Dspectrum.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# A set of subroutines and function used when calibrating 1D spectra.
# Performs correction for telluric absorption, galactic ectinction, airmass
# etc.
#----------------------------
#   COMMENTS
#----------------------------
# The fits table expected by .readfitsspec1D should contain the WAVELENGTH 
# keywords and is similar to what is created with extractMOSFIRE1Dsignal2noiseSpec.py
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# >>> import calibrate1Dspectrum as cal
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-05-24  started by K. B. Schmidt (UCSB)
#----------------------------
"""
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs    # 
import numpy as np            # enable opening with genfromtxt
import astropysics.obstools   #.CardelliExtinction  #.Site.nightTable
import astropysics
import pyfits
import commands               # get output from spawned command line processes
import scipy                  # for integration etc.
import observer               # used to obtain almanacs and airmass values 
import pywcs
import pytz                   # getting timezone information
import datetime               # used to convert utc
import ephem                  # time manipulations
import sys
import pdb                    # for debugging with pdb.set_trace()
#-------------------------------------------------------------------------------------------------------------
__version__ = 1.0 
__author__ = "K. B. Schmidt (UCSB)"
#-------------------------------------------------------------------------------------------------------------
#                                            MAIN ROUTINES
#-------------------------------------------------------------------------------------------------------------
def convertunits(wave,spec1D,area,trans,epsilon,disp,conv='eps2flux'): 
    """
    ---- PURPOSE ----
    Convertin the units of a 1D spectrum in [e/s/pix] to [erg/s/cm2/A]

    ---- INPUT ----
    wave            wavelengths of spec1D in [A]
    spec1D          1D spectrum in [e/s/pix]
    area            Collecting area of telescope in [cm2]. Keck has area=76000 cm
    trans           Transmission (or rather total throughput) curve of pass-band
                    interpolated to the wavelength of spec1D
    epsilon         Fraction of light of object falling on the N pixels in 
                    spatial direction of slit prior to collapsing
    disp            dispersion in [A/pix]
    conv            Conversion to perform. Choices are:
                      eps2flux   e/s/pix      ->  erg/s/cm2/A  DEFAULT
                      flux2eps   erg/s/cm2/A  ->  e/s/pix

    ---- OUTPUT ----
    spec1D          

    ---- EXAMPLE OF USAGE ----
    # covert units from e/s to erg/s/cm2/A
    import calibrate1Dspectrum as cal
    spec1D_tel_flux  = cal.convertunits(wave,spec1D_eps,76000,throughput_interp,epsilon,1.0855,conv='eps2flux')

    """
    h        = 6.6261e-27      # cm2*g/s
    c        = 2.9979e+10      # cm/s
    lhc      = wave*10**-8/c/h # photon energy in 1/erg
    convfact = area * lhc * trans * epsilon * disp

    if   conv == 'eps2flux':
        spec1Dconv = spec1D/convfact
    elif conv == 'flux2eps':
        spec1Dconv = spec1D*convfact
    else:
        sys.exit(':: convertunits :: ERROR: Chosen conversion not available --> ABORTING')

    return spec1Dconv
#-------------------------------------------------------------------------------------------------------------
def correct_galacticext(wave,spec1D,extval,extlaw='Cardelli',plot=0): 
    """
    ---- PURPOSE ----
    Correcting 1D spectrum for galactic extinction using the Cardelli et al. 1989
    extinction law

    ---- INPUT ----
    wave            wavelengths of spec1D in [A]
    spec1D          1D spectrum in [e/s] or [erg/s/cm2/A] or ...
    extval          redenning used for extinction correction. 
                    For 'Cardelli' the value if E(B-V) is expected
                    For 'Calzetti' A0 is expected
    extlaw          Extinction law to use. Options are:
                      'Cardelli' from Cardelli et al. 1989 (DEFAULT)
                      'Calzetti' from Calzetti et al. 1994

    ---- OUTPUT ----
    spec_corr       The extinction corrected spectrum

    ---- EXAMPLE OF USAGE ----

    """
    if extlaw == 'Cardelli':
        extlaw = astropysics.obstools.CardelliExtinction(EBmV=extval, Rv=3.1)
        #extlaw = CardelliExtinction_self(wave,extval, Rv=3.1)
    elif extlaw == 'Calzetti':
        extlaw = astropysics.obstools.CalzettiExtinction(A0=extval)

    spec_in   = astropysics.spec.Spectrum(wave, spec1D, err=None, ivar=None, unit='wl', name='spec1D', copy=True, sort=True)
    spec_corr = extlaw.correctSpectrum(spec_in,newspec=True)

    return spec_corr.flux
#-------------------------------------------------------------------------------------------------------------
    def correctSpectrum(self,spec,newspec=True):
        """
        Uses the supplied extinction law to correct a spectrum for extinction.
        
        if newspec is True, a copy of the supplied spectrum will have the 
        extinction correction applied
        
        returns the corrected spectrum
        """
    
        if newspec:
            spec = spec.copy()
            
        oldunit = spec.unit
        spec.unit = 'wavelength-angstrom'
        corr = 10**(self(spec.x)/2.5)
        spec.flux *= corr
        spec.err *= corr
        
        spec.unit = oldunit
        return spec

#-------------------------------------------------------------------------------------------------------------
def correct_airmass(wave,spec,radec,date,start,stop,plot=0): 
    """
    ---- PURPOSE ----
    Correcting 1D spectrum for atmospheric distortion caused by airmass

    ---- INPUT ----
    wave            wavelengths of spec1D in [A]
    spec            1D spectrum in [erg/s/cm2/A] (or [e/s])
    radec           ra and dec of objects on sky, i.e., [ra,dec]
    date            The evening date of obsnight to calculate airmass correction for. String of type '2013/04/25'
    start           The start of the observations in UTC. String of type 'YYYY-mm-dd HH:MM'. 
                    Used to obtain the (average) airmass of the spectrum.
    stop            The end of the observations in UTC. String of type 'YYYY-mm-dd HH:MM'. 
                    Used to obtain the (average) airmass of the spectrum.
    ---- OUTPUT ----
    speccorr        the spectrum corrected for airmass

    ---- EXAMPLE OF USAGE ----
    >>> correct_airmass(wave,spec,[227.53731361,11.26085244],'2013/04/25','2013-04-26 10:00','2013-04-26 15:00',plot=1)

    """
    localtimeBAD, utctime, AMvec = airmass(radec[0],radec[1],date,plot=plot,location='keck')

    t_data     = [str(t.datetime())[0:16] for t in utctime]  # turning UTC times in to date and time ('rounding' to 5 minute intervals)
    utcstart   = datetime.datetime.strptime(start, '%Y-%m-%d %H:%M') # converting string to datetime instance
    utcstop    = datetime.datetime.strptime(stop, '%Y-%m-%d %H:%M')  # converting string to datetime instance

    utcstart_s = ephem.Date(utcstart.replace(tzinfo=pytz.timezone('utc'))) # pytz.timezone('US/Hawaii') # utcstart in seconds
    utcstop_s  = ephem.Date(utcstop.replace(tzinfo=pytz.timezone('utc')))  # pytz.timezone('US/Hawaii') # utcstop in seconds

    if (utcstart_s > utctime[-1]) or (utcstart_s < utctime[0]):
        sys.exit('ERROR: UTC Start date is not in current night (outside airmass data range) --> ABORTING')
    if (utcstop_s  > utctime[-1]) or (utcstop_s  < utctime[0]):
        sys.exit('ERROR: UTC Stop date is not in current night (outside airmass data range) --> ABORTING')

    utcstart_diff = abs(np.asarray(utctime)-float(utcstart_s))
    utcstop_diff  = abs(np.asarray(utctime)-float(utcstop_s))

    startent   = np.where(utcstart_diff == np.min(utcstart_diff))[0]
    stopent    = np.where(utcstop_diff == np.min(utcstop_diff))[0] 

    speccorr   = spec*0.0
    for ii in range(len(wave)):
        AM           = np.mean(AMvec[startent[0]:stopent[0]+1])                  # the airmass for wavelength ii
        kappaval     = kappa(wave[ii],site='maunakea',intmethod='linear',plot=0) # the extinction coefficient for wavelength ii
        correction   = 10**(kappaval*AM/2.5)                                     # correction value to apply at each wavelength
        speccorr[ii] = spec[ii]*correction

    if plot == 1: 
        import pylab as plt
        plt.clf()
        plt.plot(wave,spec,'r-',label='input spectrum')
        plt.plot(wave,speccorr,'b--',label='airmass corrected spectrum')
        plt.legend(fancybox=True, loc='upper right')  # add the legend in the middle of the plot
        plt.show()

    return speccorr
#-------------------------------------------------------------------------------------------------------------
def correct_telluric(wave,spec,specStandard,plot=0): 
    """
    ---- PURPOSE ----
    returns the correction factor for each wavelength obtained from comparing observations
    of a telluric standard with a catalog spectrum of the spectral type (e.g. Vega for A0V) 

    ---- INPUT ----
    wave            wavelengths of spec1D in [A]
    spec            spectrum in [erg/s/cm2/A] to compare to standard 
                    (or in [e/s] so flux conversion is included in the 'correction')
    specStandard    intrinsic spectrum of the standard star to compare with in [erg/s/cm2/A]
                    (intepolated to wave and rescaled to matct expected magnitude)

    ---- OUTPUT ----


    ---- EXAMPLE OF USAGE ----
    telcorrection  = cal.correct_telluric(wave,spec_telluric,spec_vega_rescaled)
    sepc_telluric_corrected = spec_telluric * telcorrection
    """
    correction = specStandard/spec

    return correction

#-------------------------------------------------------------------------------------------------------------
#                                               UTILITIES
#-------------------------------------------------------------------------------------------------------------
def readfitsspec1D(fitstable,spec='SPEC',wave='WAVELENGTH'):
    """
    ---- PURPOSE ----
    Reading a binary fits table to obtain the 1D spectrum and the
    corresponding wavelengths.

    ---- INPUT ----
    fitstable       path an name of fits table containing data
    spec            name of the column in the fits table containing the 1D spectrum
                    DEFAULT = 'SPEC'
    wave            name of the column in the fits table containing the wavelengths
                    DEFAULT = 'WAVELENGTH'

    ---- OUTPUT ----
    spec1D          numpy array with 1D spectrum
    wave            numpy array with wavelengts
    coords          [ra,dec] in degrees of object

    ---- EXAMPLE OF USAGE ----
    >>> import calibrate1Dspectrum as cal
    >>> wave, spec1D = cal.readfitsspec1D('borg_1510+1115_J_borg_1510+1115_0745_eps_1Dspec.fits',spec='SPEC1D_SUM')

    """
    specdat   = pyfits.open(fitstable)
    specdatTB = specdat[1].data
    spec1D    = specdatTB[spec]
    wave      = specdatTB[wave]
    ra        = specdat[1].header['RA']
    dec       = specdat[1].header['DEC']

    return wave, spec1D, [ra,dec]
#-------------------------------------------------------------------------------------------------------------
def scalespec(wave,spec,trans,mag):
    """
    ---- PURPOSE ----
    Rescales spectrum to provided mag

    ---- INPUT ----
    wave            wavelengths of spec1D in [A]
    spec            spectrum in [erg/s/cm2/A] to scale
    trans           the transmission curve used to calucate mag
                    (interpolated to same wavelengths as spec)
    mag             The desired integrated magnitude of the spectrum after scaling.
                    This will be determined by integrating spec over wave with magFromSpec

    ---- OUTPUT ----
    spec_scale      spectrum with scaled flux values

    ---- EXAMPLE OF USAGE ----

    """
    magspec    = magFromSpec(wave,spec,trans)          # calculate magnitude of input spectrum
    scale      = 10**((mag-magspec)/2.5)               # scale factor between mags
    spec_scale = spec / scale                          # scaling spectrum
    magscale   = magFromSpec(wave,spec_scale,trans)    # testing that the right magnitude is obtained
    lim        = 0.01                                  # accepted difference in scaled magnitude and desired magnitude

    if np.abs(magscale-mag) > lim: 
        sys.exit('ERROR: Scaled spectrum has mag '+str(magscale)+' != '+str(lim)+' of the desired mag '+str(mag)+' --> ABORTING')

    return spec_scale
#-------------------------------------------------------------------------------------------------------------
def magFromSpec(wave,spec,trans):
    """
    ---- PURPOSE ----
    Calculating the AB magnitude of a spectrum given in flux units [erg/s/cm2/A]

    ---- INPUT ----
    wave            wavelengths of spec1D in [A]
    spec            spectrum in [erg/s/cm2/A]
    trans           the transmission curve of the filter to integrate over
                    (interpolated to same wavelengths as spec)

    ---- OUTPUT ----
    magAB           The AB magnitude from integrating spectrum

    ---- EXAMPLE OF USAGE ----

    """
    #lammin = np.min(wave)
    #lammax = np.max(wave) 
    #intSpecTransQuad,errST = scipy.integrate.quad(spectrans,lammin,lammax,args=(wave,spec,trans))

    cval   = 2.99792458e18 # A/s
    STval  = np.multiply(spec,trans)
    Tval   = np.divide(cval*trans,(wave)**2)

    intSpecTrans  = scipy.integrate.trapz(STval, x=wave)
    intTrans      = scipy.integrate.trapz(Tval, x=wave)

    magAB = -2.5 * ( np.log10(intSpecTrans) - np.log10(intTrans) ) - 48.6

    return magAB
#-------------------------------------------------------------------------------------------------------------
def spectrans(waveval,wavevec,specvec,transvec):
    """
    Retunring the (interpolated) value of the spectrum for a given wavelength
    """
    STvec = specvec*transvec
    ent   = np.where(wavevec == waveval)

    if len(ent[0]) == 1 :
        STval = STvec[ent]
    else:
        wavenew    = np.sort(np.append(wavevec,waveval))
        STvecnew   = kbs.interpn(wavevec,STvec,wavenew)
        ent        = np.where(wavenew == waveval)
        STval      = STvecnew[ent]

    return STval

#-------------------------------------------------------------------------------------------------------------
def kappa(waveval,site='maunakea',intmethod='linear',plot=0):
    """
    ---- PURPOSE ----
    Returning the (interpolated) value of the extinction coefficient
    needed to correct magnitudes for atmospheric extinction (air mass)

    ---- INPUT ----
    wave            wavelength in [A]
    site            name of site (data) to use when interpolating
    intmethod       the interpolation to perform

    ---- OUTPUT ----
    kappa(lambda)   The atmospheric instinction coefficient for lambda[A]

    ---- EXAMPLE OF USAGE ----

    """

    if site == 'maunakea':
        #data taken from http://www.gemini.edu/?q=node/10790#Mauna%20Kea
        wavedat  = np.array([0.310,0.320,0.340,0.360,0.380,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.800,0.900,1.25,1.65,2.20])*1e4
        kappadat = np.array([1.37,0.82,0.51,0.37,0.30,0.25,0.17,0.13,0.12,0.11,0.11,0.10,0.07,0.05,0.015,0.015,0.033])
    else:
        sys.exit(':: kappa :: ERROR: Chosen site not available --> ABORTING') 

    wavenew    = np.sort(np.append(wavedat,waveval))

    if waveval < np.min(wavedat):
        print 'NB! Selected wavelength is shorter than shortest wavelength in data. '
        print '    Setting output to kappa of shortest wavelength in data.' 
        kappa      = kappadat[0]
        kappanew   = np.insert(kappadat,0,kappa)
    elif waveval > np.max(wavedat):
        print 'NB! Selected wavelength is longer than longest wavelength in data. '
        print '    Setting output to kappa of longest wavelength in data.' 
        kappa      = kappadat[-1]
        kappanew   = np.append(kappadat,kappa)
    else:
        kappanew   = kbs.interpn(wavedat,kappadat,wavenew,method=intmethod)
        kappa      = kappanew[np.where(wavenew == waveval)]

    if plot == 1:
        import pylab as plt
        plt.clf()
        plt.plot(wavedat,kappadat,'ro',label='data for '+site)
        plt.plot(wavenew,kappanew,'r--',label='interpolation of data')

        #kappanew   = kbs.interpn(wavedat,kappadat,wavenew,method='cubic')
        #plt.plot(wavenew,kappanew,'g--',label='cubic')
        #kappanew   = kbs.interpn(wavedat,kappadat,wavenew,method='nearest')
        #plt.plot(wavenew,kappanew,'b--',label='nearest')
        #kappanew   = kbs.interpn(wavedat,kappadat,wavenew,method='slinear')
        #plt.plot(wavenew,kappanew,'y--',label='slinear')
        #kappanew   = kbs.interpn(wavedat,kappadat,wavenew,method='quadratic')
        #plt.plot(wavenew,kappanew,'m--',label='quadratic')
        #kappanew   = kbs.interpn(wavedat,kappadat,wavenew,method='zero')
        #plt.plot(wavenew,kappanew,'r:',label='zero')

        plt.plot(waveval,kappa,'go',label='obtained value')
        plt.legend(fancybox=True, loc='upper right')  # add the legend in the middle of the plot
        plt.show()

    return kappa


#-------------------------------------------------------------------------------------------------------------
def airmass(ra,dec,date,plot=0,location='keck'):
    """
    ---- PURPOSE ----
    Obtain the airmass for a given position and date
    Uses the observer.py scripts from http://www.ucolick.org/~magee/observer/
    Converts ra and dec using skycoor from the commandline
    ---- INPUT ----
    ra              right ascension of object
    dec             declination of object
    data            string with date of almanac to use. Format: 'yyyy/mm/dd'

    ---- OPTIONAL INPUT ----
    plot            set to 1 to plot the airmass curve
    location        site for which the airmass is calulcated.
                    Default = 'keck' but see observer.site_list() for options

    ---- OUTPUT ----
    local           list of local time in seconds  NB! don't use this as it seems to be off by 1 our... use UTC!!
                    Can be converted using: t_data = [t.datetime() for t in local]
    utc             list of utc time in seconds
                    Can be converted using: t_data = [t.datetime() for t in utc]
    airmass         list of airmas values

    ---- EXAMPLE OF USAGE ----

    """
    obs       = observer.Observer(location)

    radecsex  = commands.getoutput('skycoor '+str(ra)+' '+str(dec))
    rasex     = radecsex.split(' ')[0].replace(':',' ')
    decsex    = radecsex.split(' ')[1].replace(':',' ')
    target    = obs.target('target', rasex, decsex)

    obs.almanac(date)
    #obs.almanac_data # printing almanac data to screen

    obs.airmass(target)
    #obs.airmass_data # printing airmass data to screen

    local     = obs.airmass_data[0].local # NB!! note that for BoRG objects observed on 130425 this was off by 1 hour
    utc       = obs.airmass_data[0].utc
    airmass   = obs.airmass_data[0].airmass

    if plot == 1:
        observer.plots.plot_airmass(obs, date.replace('/','')+'_airmasplot.png',telescope='keck1')

    return local, utc, airmass
#-------------------------------------------------------------------------------------------------------------
def CardelliExtinction_self(wave,EBmV,Rv=3.1):
    """
    ---- PURPOSE ----
    Returning the absolute extinction A(lambda)/A(V) at wavelength lambda for the 
    milky way extinction law from Cardelli et al. 1989

    ---- INPUT ----
    wave            the wavelength to determine correction factors for given in [A]
    EBmV            The reddening e.g. obtained from the Schlegel maps with kbs.getAv(ra,decs)
    Rv              A(V)/E(B-V) Default value is set to 3.1

    ---- OPTIONAL INPUT ----

    ---- OUTPUT ----
    extcorr         numpy array with extinction corrections to apply to spectrum

    """
    x        = 1e4/wave    #converting lambda to 1/microns as used in Cardelli et al. 1989
    a        = np.ndarray(x.shape,x.dtype)
    b        = np.ndarray(x.shape,x.dtype)
        
    if any((x<0.3)|(10<x)): # checking that all wavelengths are in proper range
        raise ValueError('Some wavelengths outside the Cardelli et al. 1989 extinction curve range')

    # defining entries of x in the various spectral renages
    irs      = (0.3 <= x) & (x <= 1.1)
    opts     = (1.1 <= x) & (x <= 3.3)
    nuv1s    = (3.3 <= x) & (x <= 5.9)
    nuv2s    = (5.9 <= x) & (x <= 8)
    fuvs     = (8 <= x) & (x <= 10)

    #Cardelli et al. 1989 Infrared
    a[irs]   = .574*x[irs]**1.61
    b[irs]   = -0.527*x[irs]**1.61

    #Cardelli et al. 1989 NIR/optical
    a[opts]  = np.polyval((.32999,-.7753,.01979,.72085,-.02427,-.50447,.17699,1),x[opts]-1.82)
    b[opts]  = np.polyval((-2.09002,5.3026,-.62251,-5.38434,1.07233,2.28305,1.41338,0),x[opts]-1.82)

    #Cardelli et al. 1989 NUV
    a[nuv1s] = 1.752-.316*x[nuv1s]-0.104/((x[nuv1s]-4.67)**2+.341)
    b[nuv1s] = -3.09+1.825*x[nuv1s]+1.206/((x[nuv1s]-4.62)**2+.263)
    y        = x[nuv2s]-5.9
    Fa       = -.04473*y**2-.009779*y**3
    Fb       = -.2130*y**2-.1207*y**3
    a[nuv2s] = 1.752-.316*x[nuv2s]-0.104/((x[nuv2s]-4.67)**2+.341)+Fa
    b[nuv2s] = -3.09+1.825*x[nuv2s]+1.206/((x[nuv2s]-4.62)**2+.263)+Fb

    #Cardelli et al. 1989 FUV
    a[fuvs]  = np.polyval((-.070,.137,-.628,-1.073),x[fuvs]-8)
    b[fuvs]  = np.polyval((.374,-.42,4.257,13.67),x[fuvs]-8)
        
    AloAv    = a + b/Rv
    return AloAv

#-------------------------------------------------------------------------------------------------------------
#                                                  END
#-------------------------------------------------------------------------------------------------------------


