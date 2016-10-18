# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# routines and wrappers to perform (force) flux measurements on the MUSE data cubes using LSDCat
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import sys
import pyfits
import commands
import numpy as np
import datetime
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
def get_MUSEcubePixelPos(outputname,IDs,ra_pix,dec_pix,lam_pix,clobber=False,verbose=True):
    """
    Save LSDCat friend fits file which can be used as input for lsd_cat_measure.py to
    make forced flux measurements of lines or at expected line locations in the data cubes

    """
    print '--- not enabled yet ---'
    return

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def save_LSDCatFriendlyFitsFile(outputname,lineIDs,objIDs,x_pix,y_pix,lam_pix,clobber=False,verbose=True):
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

    --- EXAMPLE OF USE ---
    import fluxmeasurementsMUSEcubes as fmm

    # Example for 3 objects with 4,5 and 1 lines each
    outputname = '/Users/kschmidt/work/MUSE/ciii_candidates/linecat4LSDCat_measure_161018.fits'
    lineIDs    = ['1150308501','1150308502','1150308503','1150388804']+['88801','88802','88803','88804','88805']+['99901']
    objIDs     = ['11503085']*4+[888]*5+[999]*1
    x_pix      = [120]*4+[222]*5+[333]*1
    y_pix      = [277]*4+[222]*5+[333]*1
    lam_pix    = [80,103,778,875]+[60,800,900,1200,1250]+[1111]

    fmm.save_LSDCatFriendlyFitsFile(outputname,lineIDs,objIDs,x_pix,y_pix,lam_pix,clobber=False)

    """
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