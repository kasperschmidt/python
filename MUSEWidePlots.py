# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import os
import glob
import pyfits
import kbsutilities as kbs
import numpy as np
import MUSEWidePlots as mwp
import pdb
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def Spec1DCombinationImage_e24ANDe36(date='1701XX',verbose=True):
    """
    Wrapper for gen_Spec1DCombinationImage()

    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp
    mwp.Spec1DCombinationImage_e24ANDe36(date='170130')

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading catalogs and defining ID and redshift lists'
    photcat_e24 = '/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits'
    e24_id      = pyfits.open(photcat_e24)[1].data['UNIQUE_ID']
    e24_z       = pyfits.open(photcat_e24)[1].data['REDSHIFT']

    photcat_e36 = '/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_e36_v0.2.fits'
    e36_id      = pyfits.open(photcat_e36)[1].data['ID']
    e36_z       = pyfits.open(photcat_e36)[1].data['REDSHIFT']

    IDlist      = np.append(e24_id,e36_id)
    redshifts   = np.append(e24_z, e36_z)
    specstring  = 'spectrum_*IIIII*.fits'
    wavegrid    = [4800,9300,3]
    clobber     = True
    fluxext     = 'FLUX'
    waveext     = 'WAVE_AIR'
    NrowPerSpec = 1
    specdir     = '/Users/kschmidt/work/MUSE/spectra1D/Arche170127/spectra/'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print '\n - Generating image for non-median-subtracted 1D spectra'
    outimage    = specdir+'Spec1DCOmbinationImage_e24ANDe36_'+date+'.fits'
    imgarray    = mwp.gen_Spec1DCombinationImage(IDlist,redshifts,specdir,specstring=specstring,NrowPerSpec=NrowPerSpec,
                                                 wavelengthgrid=wavegrid,fluxext=fluxext,waveext=waveext,
                                                 imagename=outimage,clobber=clobber,verbose=verbose,
                                                 scale='None',subtract='None')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print '\n - Generating image for median-subtracted 1D spectra'
    outimage    = specdir+'Spec1DCOmbinationImage_e24ANDe36_'+date+'_mediansubtracted.fits'
    imgarray    = mwp.gen_Spec1DCombinationImage(IDlist,redshifts,specdir,specstring=specstring,NrowPerSpec=NrowPerSpec,
                                                 wavelengthgrid=wavegrid,fluxext=fluxext,waveext=waveext,
                                                 imagename=outimage,clobber=clobber,verbose=verbose,
                                                 scale='None',subtract='median')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print '\n - Generating image for max-scaled and median-subtracted 1D spectra'
    outimage    = specdir+'Spec1DCOmbinationImage_e24ANDe36_'+date+'_maxscale_mediansubtracted.fits'
    imgarray    = mwp.gen_Spec1DCombinationImage(IDlist,redshifts,specdir,specstring=specstring,NrowPerSpec=NrowPerSpec,
                                                 wavelengthgrid=wavegrid,fluxext=fluxext,waveext=waveext,
                                                 imagename=outimage,clobber=clobber,verbose=verbose,
                                                 scale='max',subtract='median')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_Spec1DCombinationImage(IDlist,redshifts,spec1Ddir,specstring='spectrum_*IIIII*.fits',NrowPerSpec=1,
                               wavelengthgrid=[4800,9300,5],fluxext='FLUX',waveext='WAVE_AIR',
                               imagename='Spec1DCOmbinationImage_RENAME.fits',subtract='None',scale='None',
                               clobber=False,verbose=True):
    """

    Generate Fits file/image containing ordered (sub)samples of the extracted 1D spectra.

    --- INPUT ---
    IDlist          List of IDs to use in the seach for spectra (inserting into "specstring"
    redshift        Redshift for the IDs in "IDlist"
    spec1Ddir       Directory containing the 1D spectra to contain into a fits image.
    specstring      String used to serach for spectra. "IIIII" will be replaced with the ID from
                    idlist and the resulting string will be fed to glob searching in the "spec1Ddir"
    NrowPerSpec     The number of rows to use for each spectrum
    wavelengthgrid  The wavelength grid to interpolate spectra onto (i.e., defining the wavelength
                    spacing and the number of columns in output image)
    fluxext         Name of fits extention containing fluxes
    waveext         Name of fits extention containing wavelengths
    imagename       Name of fits file to stroe array to
    subtract        Subtract a given flux value from spectrum. Provide float, list of values or use keywords:
                       'median'   Subtract the median flux value from each spectrum
                       'mean'     Subtract the mean flux value from each spectrum
    scale           Scale each spectrum by a given value. Provide float, list of values or use keywords:
                       'max'      Scale each spectrum by maxmimum flux value

    clobber         If true, fits image will be overwritten.
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---
    See Spec1DCombinationImage_e24ANDe36()

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Setting up output image dimensions'

    Nspec     = len(IDlist)
    Nrows     = Nspec*NrowPerSpec
    zsort_ent = np.argsort(np.asarray(redshifts))
    zsort     = np.asarray(redshifts)[zsort_ent]
    IDsort    = np.asarray(IDlist)[zsort_ent]

    wavegrid  = np.arange(wavelengthgrid[0],wavelengthgrid[1],wavelengthgrid[2])
    Ncols     = len(wavegrid)+2
    imgarray  = np.zeros([Nrows,Ncols])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Filling output image with fluxes interpolated to requested wavelength grid'
    for ss in xrange(Nspec):
        objid   = str(IDsort[ss])
        globstr = spec1Ddir+specstring.replace('IIIII',objid)
        spec    = glob.glob(globstr)

        if len(spec) == 0:
            if verbose: print ' - WARNING: No spectrum found for '+objid
        else:
            if len(spec) > 1:
                if verbose: print ' - WARNING: More than two spectra found when globbing for:\n            '+globstr
                if verbose: print ' -          Will plot the first spectrum found'

            flux       = pyfits.open(spec[0])[1].data[fluxext]
            wave       = pyfits.open(spec[0])[1].data[waveext]
            fluxinterp = kbs.interpn(wave,flux,wavegrid)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if scale != 'None':
                if type(scale) == str:
                    if scale == 'max':
                        fluxinterp = fluxinterp/np.max(fluxinterp)
                    else:
                        if verbose: print ' WARNING "'+scale+'" invalid value of scale keyword; not scaling performed'
                else:
                    if len(scale) > 1:
                        fluxinterp = fluxinterp/scale[ss]
                    else:
                        fluxinterp = fluxinterp/scale
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if subtract != 'None':
                if type(subtract) == str:
                    if subtract == 'median':
                        fluxinterp = fluxinterp - np.median(fluxinterp)
                    elif subtract == 'mean':
                        fluxinterp = fluxinterp - np.mean(fluxinterp)
                    else:
                        if verbose: print ' WARNING "'+scale+'" invalid value of scale keyword; not scaling performed'
                else:
                    if len(subtract) > 1:
                        fluxinterp = fluxinterp/subtract[ss]
                    else:
                        fluxinterp = fluxinterp/subtract
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            imgarray[ss*NrowPerSpec:ss*NrowPerSpec+NrowPerSpec,0]  = zsort[ss]
            imgarray[ss*NrowPerSpec:ss*NrowPerSpec+NrowPerSpec,1]  = int(objid)
            imgarray[ss*NrowPerSpec:ss*NrowPerSpec+NrowPerSpec,2:] = fluxinterp
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Storing image to fits file '+imagename
    pyfits.writeto(imagename,imgarray,clobber=clobber)

    return imgarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =