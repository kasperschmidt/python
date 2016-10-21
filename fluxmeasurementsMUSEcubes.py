# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# routines and wrappers to perform (force) flux measurements on the MUSE data cubes using LSDCat
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import sys
import pyfits
import commands
import numpy as np
import datetime
import astropy.wcs as wcs
import fluxmeasurementsMUSEcubes as fmm
import matplotlib.pyplot as plt
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def measure_fluxes_multiplefields(linecatsearchstr='linecat*.fits',fields=None,filed_ids=None,
                                  combineoutput=True,verbose=True):
    """
    Wrapper around fmm.measure_fluxes() to mesure fluxes on cubes from multiple fields.
    Setting combineoutput=True the flux catalogs will combined to one master catalog.

    """

    print '--- not enabled yet ---'
    return
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def measure_fluxes(linecatalog, field='cdfs', field_id=15, SNthreshinit=1.0, SNthreshanal=1.0, cubeversion='_v1.0',
                   fhdu='MEDFILTERED_DATA', ferrhdu='EFF_STAT', ffhdu='FILTERED_DATA',fferhdu='FILTERED_STAT',
                   ffsnhdu='SIGNALTONOISE',rmin=3, rmax=6, dataparentpath='/Volumes/DATABCKUP3/MUSE/', clobber=False,
                   verbose=True):
    """
    Retrieve the flux measurements for an input LSDCat catalog (e.g., created with
    fmm.save_LSDCatFriendlyFitsFile) and return ...

    --- INPUT ---
    linecatalog      LSDCat catalog of (line) positions at which fluxes should be measured.
                     A mock catalog (without running LSDCat) can be generated with
                     fmm.save_LSDCatFriendlyFitsFile() to perform forced flux measurements.

    field='cdfs'
    field_id=15
    SNthreshinit=1.0
    SNthreshanal=1.0
    cubeversion      The version string used in the flux cube filenames
    fhdu             Fits extension containing flux values
    ferrhdu          Fits extension containing flux error
    ffhdu            Fits extension containing the fileterd flux
    fferhdu          Fits extension containing the fileterd flux error
    ffsnhdu          Fits extension containing the signal to noise
    rmin=3
    rmax=6
    clobber          Set to True to overwrite output if it already exists
    verbose          Toggle verbosity


    --- EXAMPLE OF USE ---
    import fluxmeasurementsMUSEcubes as fmm
    linecatalog = './linecat4LSDCat_measure_161018.fits'
    fluxcatalog = fmm.measure_fluxes(linecatalog, field='cdfs', field_id=15, verbose=True)

    """
    # paths and filenames
    fieldname     = 'candels-'+field+'-'+str(field_id)
    #linecat_base=`basename ${linecatalog} .fits`

    field_path   = dataparentpath+'/'+fieldname+'/'
    fluxcube     = field_path+'median_filtered_DATACUBE_'+fieldname+cubeversion+'.fits_effnoised.fits'
    SNcube       = field_path+'s2n_'+fieldname+'.fits'
    filteredcube = field_path+'spec_cced_spat_cced_median_filtered_DATACUBE_'+fieldname+cubeversion+'_effnoised.fits'

    if verbose:
        print ' - Will measure fluxes using the following setup:'
        print """
    Field                 : %s
    Field ID              : %s
    Catalog               : %s
    fluxcube              : %s
    filtered cube         : %s
    S/N cube              : %s
    S/N thresh (init,anal): [%s,%s]
    Apertire radius range : [%s,%s]
    """ % (field,field_id,linecatalog,fluxcube,filteredcube,fluxcube,SNthreshinit,SNthreshanal,rmin,rmax)

    if verbose: print ' - Putting together lsd_cat_measure.py from input '
    measure_cmd  = 'lsd_cat_measure.py ' \
                   ' --inputcat '+linecatalog+\
                   ' --thresh '+str(SNthreshinit)+\
                   ' --threshana '+str(SNthreshanal)+\
                   ' --fluxcube '+fluxcube+\
                   ' --fhdu '+fhdu+\
                   ' --ferrhdu '+ferrhdu+\
                   ' --filteredfluxcube '+filteredcube+\
                   ' --ffhdu '+ffhdu+\
                   ' --fferhdu '+fferhdu+\
                   ' --sncube '+SNcube+\
                   ' --ffsnhdu '+ffsnhdu+\
                   ' --rmin '+str(rmin)+\
                   ' --rmax '+str(rmax)

    if clobber:
        if verbose: print ' - NB: Clobber=True so adding "--clobber" to lsd_cat_measure.py command'
        measure_cmd = measure_cmd+' --clobber '

    if verbose: print ' - The command to spawn  :\n    '+measure_cmd+'\n'
    nowstr  = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    if verbose: print ' - Measuring fluxes (might take a while)'
    if verbose: print '   Started on '+nowstr
    if verbose: print '   --------------------------------- lsd_cat_measure.py output -----------------------------------'
    lsdcout = commands.getoutput(measure_cmd)
    print lsdcout
    if verbose: print '   -----------------------------------------------------------------------------------------------'
    nowstr  = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    if verbose: print '   Finished on '+nowstr

    fluxcatalog = lsdcout.split('DONE!!! Wrote FITS catalog ')[-1].split(' ')[0]
    return fluxcatalog
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_cubepixelpos(cube,ra,dec,wavelength,cubeextension='MEDFILTERED_DATA',verbose=True):
    """
    Obtain the pixel position in a

    --- EXAMPLE OF USE ---

    import fluxmeasurementsMUSEcubes as fmm

    cube = '/Volumes/DATABCKUP3/MUSE//candels-cdfs-15/median_filtered_DATACUBE_candels-cdfs-15_v1.0.fits_effnoised.fits'
    ra_pix, dec_pix, lam_pix = fmm.get_cubepixelpos(cube,53.113493,-27.858925,6468.3)  ; lam_pix=1375
    ra_pix, dec_pix, lam_pix = fmm.get_cubepixelpos(cube,53.113493,-27.858925,8506.61) ; lam_pix=3006

    """
    cubehdr         = pyfits.open(cube)[cubeextension].header
    lam_pix         = fmm.get_cubepixelpos_wavelength(cubehdr,wavelength)
    ra_pix, dec_pix = fmm.get_cubepixelpos_spatial(cubehdr,ra,dec)

    if verbose:
        print '                            ra            dec           wavelength '
        print ' - Coordinates    :'+str("%15.8f" % ra)+' '+str("%15.8f" % dec)+' '+str("%15.4f" % wavelength)
        print ' - Pixel position :'+str("%15.4f" % ra_pix)+' '+str("%15.4f" % dec_pix)+' '+str("%15.4f" % lam_pix)

    return ra_pix, dec_pix, lam_pix
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_cubepixelpos_spatial(cubeheader,ra,dec,verbose=True):
    """
    Get pixel position in MUSE cube in spatial direction

    """
    del cubeheader['COMMENT'] # wcs doesn't like the comment "These data have been ZAPped!" so deleting it
    wcsobj          = wcs.WCS(cubeheader)
    lamdummy        = -99
    ra_pix, dec_pix = wcsobj.wcs_world2pix(np.array([[ra,dec,lamdummy]]),1)[0][0:2] # wcsobj.wcs_pix2world()

    return ra_pix, dec_pix
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_cubepixelpos_wavelength(cubeheader,wavelength,indexed=0,verbose=True):
    """
    Get pixel position in MUSE cube in wavelenght direction

    """
    resolution  = cubeheader['CD3_3'] # wavelength resolution in Angstrom
    waveref_val = cubeheader['CRVAL3']
    Nwaves      = cubeheader['NAXIS3']

    wavevector  = np.arange(waveref_val,waveref_val+Nwaves*resolution,resolution)

    wavediff    = np.abs(wavevector-wavelength)
    lam_pix     = np.where(wavediff == np.min(wavediff))

    if len(lam_pix) > 1:
        if verbose: print ' - WARNING: more than 1 wavelength slice matches minimum; using the lowest wavelength'
    lam_pix     = lam_pix[0]

    if indexed == 1:
        lam_pix = lam_pix + 1.0
    return lam_pix
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def save_LSDCatFriendlyFitsFile(outputname,lineIDs,objIDs,x_pix,y_pix,lam_pix,clobber=False,
                                radecwave=False,coordcube=None,verbose=True):
    """
    Save LSDCat friendly fits file which can be used as input for lsd_cat_measure.py to make
    forced flux measurements of lines or at expected line locations in the MUSE data cubes

    --- INPUT ---
    outputname  File name of output fits file to generate
    lineIDs     Unique IDs of the lines found in the objects. These IDs can potentially
                be 'appended' to the objIDs to make it easier to tie them to the parent object.
    objIDs      Unique IDs of objects to store in file. These IDs tie togehter lines in
                lineIDs list coming from the same objects.
    ra_pix      Pixel postion of line (SN peak) along the x direction in the cube
    dec_pix     Pixel postion of line (SN peak) along the y direction in the cube
    lam_pix     Pixel postion of line (SN peak) along the wavelength direction
    clobber     Set to true to overwrite existing output catalog
    radecwave   If True x_pix, y_pix and lam_pix are expected to be ra [deg], dec [deg] observed wavelength [A]
                These will then be converted to pixel positions in the cube given in coordcube
    coordcube   If x_pix, y_pix and lam_pix are given as ra, dec and wavelength (radecwave=True) provide
                the cube to be used for the coordinate (wcs) transformation from coordinates to pixel postions.
    verbose     Toggle verbosity

    --- EXAMPLE OF USE ---
    import fluxmeasurementsMUSEcubes as fmm

    # >>> Example for 3 objects with 4,5 and 1 lines each <<<
    outputname = '/Users/kschmidt/work/MUSE/ciii_candidates/linecat4LSDCat_measure_161018.fits'
    lineIDs    = ['1150308501','1150308502','1150308503','1150388804']+['88801','88802','88803','88804','88805']+['99901']
    objIDs     = ['11503085']*4+[888]*5+[999]*1
    x_pix      = [120]*4+[222]*5+[333]*1
    y_pix      = [277]*4+[222]*5+[333]*1
    lam_pix    = [80,103,778,875]+[60,800,900,1200,1250]+[1111]

    fmm.save_LSDCatFriendlyFitsFile(outputname,lineIDs,objIDs,x_pix,y_pix,lam_pix,clobber=False)

    # >>> Example for multiple lines for a single cdfs 15 cube ra and dec. <<<
    linecat  = '/Users/kschmidt/work/MUSE/ciii_candidates/linecat4LSDCat_measure_cdfs15oiiemitter_161021.fits'
    obj_z    = 6474.0/3729.0-1.0
    obj_id   = '11503085'
    obj_ra   = 53.113493
    obj_dec  = -27.858925
    wave_obs = np.array([3726,3729,4340.47,9999,4861.33,4959,5007,6562.8])*(obj_z+1) # NB: Halpha outside MUSE range
    objIDs   = np.asarray( [obj_id]*len(wave_obs) )
    lineIDs  = np.asarray( [obj_id+str("%.3d" % (nn+1)) for nn in xrange(len(wave_obs))] )
    ras      = np.asarray( [obj_ra]*len(wave_obs) )
    decs     = np.asarray( [obj_dec]*len(wave_obs) )
    cube     = '/Volumes/DATABCKUP3/MUSE//candels-cdfs-15/median_filtered_DATACUBE_candels-cdfs-15_v1.0.fits_effnoised.fits'

    fmm.save_LSDCatFriendlyFitsFile(linecat,lineIDs,objIDs,ras,decs,wave_obs,radecwave=True,coordcube=cube,clobber=False)

    fluxcatalog = fmm.measure_fluxes(linecat, field='cdfs', field_id=15, verbose=True, clobber=False)

    """
    if radecwave:
        if verbose: print ' - RA, Dec and wavelength provided; converting coordinates to pixel positions in cube:' \
                          '\n   '+coordcube
        for ii in xrange(len(x_pix)):
            if (lam_pix[ii] > 4775.0) & (lam_pix[ii] < 9325.0):
                ra_pix, dec_pix, wave_pix = fmm.get_cubepixelpos(coordcube,x_pix[ii],y_pix[ii],lam_pix[ii],verbose=False)
                x_pix[ii]   = ra_pix
                y_pix[ii]   = dec_pix
                lam_pix[ii] = wave_pix
            else:
                if verbose: print '   WARNING skipping line at '+str(lam_pix[ii])+'A as it is outside MUSE wavelength range'
                lam_pix[ii] = -99

        x_pix   = np.delete(x_pix,   np.where(lam_pix == -99))
        y_pix   = np.delete(y_pix,   np.where(lam_pix == -99))
        lineIDs = np.delete(lineIDs, np.where(lam_pix == -99))
        objIDs  = np.delete(objIDs,  np.where(lam_pix == -99))
        lam_pix = np.delete(lam_pix, np.where(lam_pix == -99))

    datadic = {}
    datadic['I']          = np.asarray(lineIDs)
    datadic['ID']         = np.asarray(objIDs)
    datadic['X_PEAK_SN']  = np.asarray(x_pix)
    datadic['Y_PEAK_SN']  = np.asarray(y_pix)
    datadic['Z_PEAK_SN']  = np.asarray(lam_pix)

    Nlines                = len(datadic['I'])
    if verbose: print ' - Checking dimensions of input data (they should all have Nlines='+str(Nlines)+' entries)'
    lenarr = [len(datadic['ID']),len(datadic['X_PEAK_SN']),len(datadic['Y_PEAK_SN']),len(datadic['Z_PEAK_SN'])]

    if lenarr != [len(datadic['I'])]*4:
        sys.exit(' Dimensions of input data are not as expected ')
    else:
        if verbose: '   Input data looks good'
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Setting up array for fits output \n   '+outputname
    columndefs = []
    fitsformat = ['A20','A20','D','D','D']
    for kk, key in enumerate(['I','ID','X_PEAK_SN','Y_PEAK_SN','Z_PEAK_SN']):
        try:
            #dtype=[('I', 'S20'), ('ID', 'S20'), ('X_PEAK_SN', '>f8'), ('Y_PEAK_SN', '>f8'), ('Z_PEAK_SN', '>f8')])
            columndefs.append(pyfits.Column(name=key  , format=fitsformat[kk], array=datadic[key]))
        except:
            print ' ----ERROR---- in defining columns for fits file --> stopping with pdb.set_trace() to investigate'
            pdb.set_trace()

    cols     = pyfits.ColDefs(columndefs)
    tbhdu    = pyfits.new_table(cols)          # creating table header
    hdu      = pyfits.PrimaryHDU()             # creating primary (minimal) header
    thdulist = pyfits.HDUList([hdu, tbhdu])    # combine primary and table header to hdulist
    thdulist.writeto(outputname,clobber=clobber)  # write fits file (clobber=True overwrites excisting file)
    if verbose: print '   Wrote array to output file; done.'
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =