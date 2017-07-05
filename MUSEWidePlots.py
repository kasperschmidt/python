# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import os
import glob
import pyfits
import kbsutilities as kbs
import numpy as np
import MUSEWidePlots as mwp
import sys
import matplotlib.pyplot as plt
import pdb
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def loadcatalogs(verbose=True):
    """
    Loading catalogs and returning combined MUSE-Wide catalog IDs, ra, dec, redshifts etc.

    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp
    IDlist, ra, dec, redshifts = mwp.loadcatalogs()

    """
    if verbose: print ' - Loading catalogs and defining ID, RA, Dec and redshift lists'
    photcat_e24 = '/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits'
    e24_id      = pyfits.open(photcat_e24)[1].data['UNIQUE_ID']
    e24_ra      = pyfits.open(photcat_e24)[1].data['RA']
    e24_dec     = pyfits.open(photcat_e24)[1].data['DEC']
    e24_z       = pyfits.open(photcat_e24)[1].data['REDSHIFT']

    #photcat_e36 = '/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_e36_v0.2.fits'
    photcat_e36 = '/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_e36_v1.0.fits'
    e36_id      = pyfits.open(photcat_e36)[1].data['ID']
    e36_ra      = pyfits.open(photcat_e36)[1].data['RA']
    e36_dec     =  pyfits.open(photcat_e36)[1].data['DEC']
    e36_z       = pyfits.open(photcat_e36)[1].data['REDSHIFT']

    IDlist      = np.append(e24_id,e36_id)
    ra          = np.append(e24_ra,e36_ra)
    dec         = np.append(e24_dec,e36_dec)
    redshifts   = np.append(e24_z, e36_z)

    return IDlist, ra, dec, redshifts
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def Spec1DCombinationImage_e24ANDe36(date='1701XX',verbose=True):
    """
    Wrapper for gen_Spec1DCombinationImage()

    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp
    imagearray = mwp.Spec1DCombinationImage_e24ANDe36(date='170130')

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IDlist, ra, dec, redshifts = mwp.loadcatalogs(verbose=verbose)

    specstring  = 'spectrum_*IIIII*.fits'
    wavegrid    = [4800,9300,3]
    clobber     = True
    fluxext     = 'FLUX'
    waveext     = 'WAVE_AIR'
    NrowPerSpec = 1
    specdir     = '/Users/kschmidt/work/MUSE/spectra1D/Arche170127/spectra/'

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #if verbose: print '\n - Makeing sub-selection of objects'

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print '\n - Generating image for non-median-subtracted 1D spectra'
    outimage    = specdir+'Spec1DCOmbinationImage_e24ANDe36_'+date+'.fits'
    imgarray    = mwp.gen_Spec1DCombinationImage(IDlist,redshifts,specdir,specstring=specstring,NrowPerSpec=NrowPerSpec,
                                                 wavelengthgrid=wavegrid,fluxext=fluxext,waveext=waveext,
                                                 imagename=outimage,clobber=clobber,verbose=verbose,
                                                 scale='None',subtract='None')
    returnname  = outimage
    returnarray = imgarray
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

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print '\n - Returning resulting image array from '+returnname
    return returnarray
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
    Ninfopixs = 2
    Ncols     = len(wavegrid)+Ninfopixs
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
    hdu = pyfits.PrimaryHDU(imgarray)
    # writing hdrkeys:'---KEY--',                             '----------------MAX LENGTH COMMENT-------------'
    hdu.header.append(('EXTNAME ', 'SPECIMG'                  ,'Nem of extension'),end=True)
    hdu.header.append(('CRPIX1  ', 1.0+Ninfopixs              ,''),end=True)
    hdu.header.append(('CDELT1  ', np.unique(np.diff(wavegrid))[0]  ,''),end=True)
    hdu.header.append(('CRVAL1  ', wavegrid[0]                ,''),end=True)
    hdu.header.append(('CUNIT1  ', 'Angstrom'                 ,''),end=True)
    hdu.header.append(('CTYPE1  ', 'WAVE    '                 ,''),end=True)
    hdu.header.append(('CRPIX2  ', 1.0                        ,''),end=True)
    hdu.header.append(('CRVAL2  ', 1.0                        ,''),end=True)
    hdu.header.append(('CDELT2  ', 1.0                        ,''),end=True)
    hdu.header.append(('CUNIT2  ','Nspec  '                   ,''),end=True)
    hdu.header.append(('CTYPE2  ','INT  '                     ,''),end=True)

    hdu.writeto(imagename,clobber=clobber)

    return imgarray

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_DS9regionfile(outputfile_base,circlesize=0.5,width=3,addIDlabel=True,fontsize=12,clobber=False,verbose=True):
    """
    Generate DS9 region files of objects in MUSE-Wide catalogs

    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp
    outbase  = '/Users/kschmidt/work/catalogs/MUSE_GTO//MUSE-Wide_objects'
    file_cdfs, file_cosmos = mwp.gen_DS9regionfile(outbase,circlesize=1.0,width=3,fontsize=15,addIDlabel=True,clobber=True)


    """
    if verbose: print ' - Generating MUSE-Wide region files. Saving them to :'
    outputfile_CDFS   = outputfile_base+'_cdfs.reg'
    outputfile_COSMOS = outputfile_base+'_cosmos.reg'
    if verbose: print '   '+outputfile_CDFS
    if verbose: print '   '+outputfile_COSMOS

    IDlist, ra, dec, redshifts = mwp.loadcatalogs(verbose=verbose)

    for outputfile in [outputfile_CDFS,outputfile_COSMOS]:
        if not clobber:
            if os.path.isfile(outputfile):
                sys.exit('The file '+outputfile+' already exists and clobber = False ')
        fout = open(outputfile,'w')

        fout.write("# Region file format: DS9 version 4.1 \nfk5\n")
        if 'cosmos' in outputfile:
            fieldent = np.where(ra > 100)[0]
        else:
            fieldent = np.where(ra < 100)[0]

        RAs  = ra[fieldent]
        DECs = dec[fieldent]
        IDs  = IDlist[fieldent]
        zs   = redshifts[fieldent]

        for rr, raval in enumerate(RAs):
            color = 'white'
            if zs[rr] > 2.9: color = 'red'

            string = 'circle('+str(raval)+','+str(DECs[rr])+','+str(circlesize)+'") # color='+color+' width='+str(width)+' '

            if addIDlabel:
                string = string+' font="times '+str(fontsize)+' bold roman" text={'+str(IDs[rr])+'}'

            fout.write(string+' \n')

        fout.close()

    return outputfile_CDFS, outputfile_COSMOS

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_redshifthistograms(plotname='/Users/kschmidt/work/MUSE/MUSEWide_sourcehistogram.pdf',
                            zranges=[[0,2.5],[2.5,7]],colors=['black','red'],labels=None,fill=False,
                            bins=None,xrange=None,ylabel=' ',verbose=True):
    """

    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp
    IDs, zs = mwp.plot_redshifthistograms(labels=['Other (H$\\alpha$, [OIII], [OII])','Ly$\\alpha$'],fill=True)

    IDs, zs = mwp.plot_redshifthistograms(zranges=[[2.6,7]],colors=['red'],labels=['Ly$\\alpha$'],fill=True,xrange=[2.6,7])


    """
    if verbose: print ' - Loading source catalogs'
    IDlist, ra, dec, redshifts = mwp.loadcatalogs(verbose=verbose)

    if verbose: print ' - Setting up and generating histogram of MUSE-Wide sources in \n   '+plotname
    fig = plt.figure(figsize=(10, 3))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.98, bottom=0.2, top=0.95)
    Fsize    = 15
    lthick   = 2
    marksize = 3
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title('M^\star',fontsize=Fsize)

    if xrange is None:
        xrange = [np.min(redshifts),np.max(redshifts)]

    if bins is None:
        bin_dz = 0.1
        bins   = np.arange(np.min(redshifts),np.max(redshifts)+bin_dz,bin_dz)

    outID = []
    outz  = []
    for zz, zrange in enumerate(zranges):
        goodent   = np.where((redshifts > zrange[0]) &(redshifts <= zrange[1]))[0]
        goodIDs   = IDlist[goodent]
        goodz     = redshifts[goodent]
        goodcolor = colors[zz]
        zmin      = np.min(goodz)
        zmax      = np.max(goodz)
        outID.append(goodIDs)
        outz.append(goodz)

        infostr   = '   Histinfo:'
        if labels is not None:
            goodlabel   = labels[zz]
            infostr     = infostr+'  label = '+goodlabel

        infostr   = infostr+'  zrange = '+str(zrange)+'  Nobj = '+str(len(goodent))+'  color = '+goodcolor+\
                    '  zmin = '+str("%.5f" % zmin)+'  zmax = '+str("%.5f" % zmax)
        if verbose: print infostr

        hist = plt.hist(goodz,color=goodcolor,bins=bins,histtype="step",lw=lthick,label=goodlabel,
                        fill=fill,fc=goodcolor)

    plt.xlim(xrange)
    # xticksstrings = ['A','B','C','D','E','F']
    # plt.xticks(xvec, xticksstrings)
    plt.xlabel(r'redshift', fontsize=Fsize)

    #plt.ylim(yrange)
    plt.ylabel(ylabel, fontsize=Fsize)

    if labels is not None:
        #--------- LEGEND ---------
        anchorpos = (0.5, 1.2)
        leg = plt.legend(fancybox=True,numpoints=1, loc='upper center',prop={'size':Fsize-3},ncol=len(labels))#,
                         #bbox_to_anchor=anchorpos)  # add the legend
        leg.get_frame().set_alpha(0.7)
        #--------------------------

    if verbose: print '   Saving plot... '
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

    return outID, outz
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def compare_TDOSEspecs(tdosespecssearch='/Volumes/DATABCKUP2/TDOSEextractions/spec_comparisons/tdose_cdfs-04_170630/*fits',
                       caruanaspecssearch='/Volumes/DATABCKUP2/TDOSEextractions/spec_comparisons/caruana_cdfs-04/*',
                       urrutiaspecssearch='/Users/kschmidt/work/MUSE/spectra1D/Arche170127/spectra/spectrum_104*fits',
                       TDOSE2MWidmatch='/Users/kschmidt/work/catalogs/MUSE_GTO/MW_1-24_main_table_v3.2.fits',
                       colors=['black','red','green'],labels=['TDOSE Gauss','Caruana','Urrutia'],
                       outdir=None,skipobj=False,dwave=70,verbose=True):
    """
    Generate plots and print comparison statistics for TDOSE, Caruana and Utturia 1D spectra

    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp
    mwp.compare_TDOSEspecs()

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loaking for spectra by globbing using the input path-strings'
    tdosespecs   = glob.glob(tdosespecssearch)
    caruanaspecs = glob.glob(caruanaspecssearch)
    urrutiaspecs = glob.glob(urrutiaspecssearch)

    Ntdose =  len(tdosespecs)
    if Ntdose == 0:
        sys.exit('No tdose spectra found in '+tdosespecssearch)
    else:
        if verbose: print ' - Will generate comparison plots for the '+str(Ntdose)+' TDOSE spectra found '

    if len(caruanaspecs) == 0:
        if verbose: print(' - WARNING No Caruana spectra found in \n   '+caruanaspecssearch)
    if len(urrutiaspecs) == 0:
        if verbose: print(' - WARNING No Urrutia spectra found in \n   '+urrutiaspecssearch)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Putting together lists of IDs'
    if verbose: print '   Generate ID list for TDOSE spectra'
    ids_tdose = [float(tt.split('/')[-1].split('_')[-1].split('.')[0]) for tt in tdosespecs]
    ids_tdose = np.asarray(ids_tdose)

    if verbose: print '   Generate ID list for Caruana spectra'
    ids_car = [float(tt.split('/')[-1].split('_')[-1].split('.')[0]) for tt in caruanaspecs]
    ids_car = np.asarray(ids_car)

    if len(urrutiaspecs) > 0:
        if verbose: print '   Generate ID list for Urrutia spectra'
        ids_urr_MW = [float(tt.split('/')[-1].split('_')[-1].split('.')[0]) for tt in urrutiaspecs]
        ids_urr_MW = np.asarray(ids_urr_MW)

        idmatch    = pyfits.open(TDOSE2MWidmatch)[1].data
        ids_urr    = []
        rmatch_urr = []
        for MWid in ids_urr_MW:
            if float(str(MWid)[:3]) <= 124: # add 0 to id for CDFS 1-24 fields
                MWid_9digit = str(int(MWid))[:3]+'0'+str(int(MWid))[3:]
            else:
                MWid_9digit = str(int(MWid))
            matchent = np.where(idmatch['UNIQUE_ID'] == MWid_9digit)[0]
            ids_urr.append(idmatch['Guo_ID'][matchent])
            rmatch_urr.append(idmatch['Guo_sep'][matchent])

        ids_urr    = np.asarray(ids_urr)
        rmatch_urr = np.asarray(rmatch_urr)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if outdir is None:
        outdir = '/'.join(os.path.abspath(tdosespecs[0]).split('/')[:-1])+'/'
    if verbose: print ' - Will store output figures in \n   '+outdir

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Looping over objects and generating plots'

    for ii, tspec in enumerate(tdosespecs):
        filename = tspec.split('/')[-1]
        objidstr = filename.split('_')[-1].split('.fit')[0]
        plotname = outdir+filename.replace('.fits','_TCUcomparison.pdf')
        skipthisobj = False
        if os.path.isfile(plotname) & skipobj:
            skipthisobj = True

        if verbose:
            infostr = ' - Generate plot for ID_TDOSE = '+objidstr+' ('+str(ii+1)+'/'+str(Ntdose)+') '
            if skipthisobj:
                infostr = infostr+'... plot exists --> skipobj'
            else:
                infostr = infostr+'                           '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        if skipthisobj: continue
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        dat_tdose    = pyfits.open(tspec)[1].data
        t_flux    = dat_tdose['flux']
        t_fluxerr = dat_tdose['fluxerror']
        t_wave    = dat_tdose['wave']
        t_s2n     = dat_tdose['s2n']

        cent   = np.where(ids_car == ids_tdose[ii])[0]
        if len(cent) == 1:
            c_flux    = pyfits.open(caruanaspecs[cent[0]])[0].data
            c_fluxerr = pyfits.open(caruanaspecs[cent[0]])[1].data
            c_wave    = t_wave

        else:
            c_flux    = None
            c_fluxerr = None
            c_wave    = None

        if len(urrutiaspecs) > 0:
            uent   = np.where(ids_urr == str(int(ids_tdose[ii])))[0]
            if (len(uent) == 1) & (uent != 0):
                dat_urr   = pyfits.open(urrutiaspecs[uent[0]])[1].data
                u_flux    = dat_urr['FLUX']
                u_fluxerr = dat_urr['FLUXERR']
                u_wave    = dat_urr['WAVE_AIR']
                u_rmatch  = rmatch_urr[uent]
            else:
                u_flux    = None
                u_fluxerr = None
                u_wave    = None
                u_rmatch  = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        maxs2n    = np.max(t_s2n[ np.isfinite(t_s2n) ])
        entfix0   = np.where(t_s2n == maxs2n)[0][0]
        wavezoom0 = t_wave[entfix0]

        wavezoom1    = 5500
        diff1         = np.abs(t_wave-wavezoom1)
        entfix1       = np.where(diff1 == np.min(diff1))[0][0]

        wavezoom2    = 8500
        diff2         = np.abs(t_wave-wavezoom2)
        entfix2       = np.where(diff2 == np.min(diff2))[0][0]

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Ncol         = 3
        Nrow         = 5
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        fig = plt.figure(figsize=(Ncol*2., Nrow))
        fig.subplots_adjust(wspace=0.4, hspace=0.6,left=0.1, right=0.98, bottom=0.07, top=0.96)
        Fsize  = 5
        lthick = 1
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        #plt.title(plotname.split('/')[-1].replace('_','\_'),fontsize=Fsize)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        drawboxes = [ [wavezoom0-dwave,wavezoom0+dwave],
                      [wavezoom1-dwave,wavezoom1+dwave],
                      [wavezoom2-dwave,wavezoom2+dwave] ]


        #---------------------------- Text... ----------------------------
        rowval  = 0
        rowspan = 1
        colval  = 0
        colspan = 3
        ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

        mwp.gen_compare_text(ax,filename.replace('_','\_'),
                             t_wave,t_flux,t_fluxerr,c_wave,c_flux,c_fluxerr,u_wave,u_flux,u_fluxerr,
                             u_rmatch,ids_urr_MW[uent],colors,fontsize=Fsize)
        #---------------------------- Full Spec Flux  ----------------------------
        spectitle = 'Full Flux Spectrum'
        rowval  = 1
        rowspan = 1
        colval  = 0
        colspan = 3
        ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

        mwp.gen_compare_spec(ax,spectitle,t_wave,t_flux,t_fluxerr,c_wave,c_flux,c_fluxerr,u_wave,u_flux,u_fluxerr,
                             xrange=[4800,9300],markwave=False,drawboxes=drawboxes,fontsize=Fsize,colors=colors,
                             labels=colors,plotSNcurve=False,shownoise=True,lthick=lthick,fillalpha=0.30,verbose=True)


        #---------------------------- Full Spec S/N  ----------------------------
        spectitle = 'Full S/N Spectrum'
        rowval  = 2
        rowspan = 1
        colval  = 0
        colspan = 3
        ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

        mwp.gen_compare_spec(ax,spectitle,t_wave,t_flux,t_fluxerr,c_wave,c_flux,c_fluxerr,u_wave,u_flux,u_fluxerr,
                             xrange=[4800,9300],markwave=False,drawboxes=drawboxes,fontsize=Fsize,colors=colors,
                             labels=colors,plotSNcurve=True,shownoise=True,lthick=lthick,fillalpha=0.30,verbose=True)

        #---------------------------- Zoom Flux 0 ----------------------------
        spectitle = '$\\lambda$(Max(S/N))$\\pm$'+str(dwave)+'\AA'
        rowval  = 3
        rowspan = 1
        colval  = 0
        colspan = 1
        ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

        mwp.gen_compare_spec(ax,spectitle,t_wave,t_flux,t_fluxerr,c_wave,c_flux,c_fluxerr,u_wave,u_flux,u_fluxerr,
                             xrange=drawboxes[0],markwave=wavezoom0,drawboxes=False,fontsize=Fsize,colors=colors,
                             labels=colors,plotSNcurve=False,shownoise=True,lthick=lthick,fillalpha=0.30,verbose=True)

        #---------------------------- Zoom Flux 1 ----------------------------
        spectitle = str(wavezoom1)+'$\\pm$'+str(dwave)+'\AA'
        rowval  = 3
        rowspan = 1
        colval  = 1
        colspan = 1
        ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

        mwp.gen_compare_spec(ax,spectitle,t_wave,t_flux,t_fluxerr,c_wave,c_flux,c_fluxerr,u_wave,u_flux,u_fluxerr,
                             xrange=drawboxes[1],markwave=wavezoom1,drawboxes=False,fontsize=Fsize,colors=colors,
                             labels=colors,plotSNcurve=False,shownoise=True,lthick=lthick,fillalpha=0.30,verbose=True)

        #---------------------------- Zoom Flux 2 ----------------------------
        spectitle = str(wavezoom2)+'$\\pm$'+str(dwave)+'\AA'
        rowval  = 3
        rowspan = 1
        colval  = 2
        colspan = 1
        ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

        mwp.gen_compare_spec(ax,spectitle,t_wave,t_flux,t_fluxerr,c_wave,c_flux,c_fluxerr,u_wave,u_flux,u_fluxerr,
                             xrange=drawboxes[2],markwave=wavezoom2,drawboxes=False,fontsize=Fsize,colors=colors,
                             labels=colors,plotSNcurve=False,shownoise=True,lthick=lthick,fillalpha=0.30,verbose=True)

        #---------------------------- Zoom S/N 0 ----------------------------
        spectitle = '$\\lambda$(Max(S/N))$\\pm$'+str(dwave)+'\AA'
        rowval  = 4
        rowspan = 1
        colval  = 0
        colspan = 1
        ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

        mwp.gen_compare_spec(ax,spectitle,t_wave,t_flux,t_fluxerr,c_wave,c_flux,c_fluxerr,u_wave,u_flux,u_fluxerr,
                             xrange=drawboxes[0],markwave=wavezoom0,drawboxes=False,fontsize=Fsize,colors=colors,
                             labels=colors,plotSNcurve=True,shownoise=True,lthick=lthick,fillalpha=0.30,verbose=True)

        #---------------------------- Zoom S/N 1 ----------------------------
        spectitle = str(wavezoom1)+'$\\pm$'+str(dwave)+'\AA'
        rowval  = 4
        rowspan = 1
        colval  = 1
        colspan = 1
        ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

        mwp.gen_compare_spec(ax,spectitle,t_wave,t_flux,t_fluxerr,c_wave,c_flux,c_fluxerr,u_wave,u_flux,u_fluxerr,
                             xrange=drawboxes[1],markwave=wavezoom1,drawboxes=False,fontsize=Fsize,colors=colors,
                             labels=colors,plotSNcurve=True,shownoise=True,lthick=lthick,fillalpha=0.30,verbose=True)

        #---------------------------- Zoom S/N 2 ----------------------------
        spectitle = str(wavezoom2)+'$\\pm$'+str(dwave)+'\AA'
        rowval  = 4
        rowspan = 1
        colval  = 2
        colspan = 1
        ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

        mwp.gen_compare_spec(ax,spectitle,t_wave,t_flux,t_fluxerr,c_wave,c_flux,c_fluxerr,u_wave,u_flux,u_fluxerr,
                             xrange=drawboxes[2],markwave=wavezoom2,drawboxes=False,fontsize=Fsize,colors=colors,
                             labels=colors,plotSNcurve=True,shownoise=True,lthick=lthick,fillalpha=0.30,verbose=True)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #if verbose: print ' - Saving plot to',plotname
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')

    if verbose: print '\n ... done '
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_compare_text(ax,title,t_wave,t_flux,t_fluxerr,c_wave,c_flux,c_fluxerr,u_wave,u_flux,u_fluxerr,u_rmatch,MWid,
                     colors,fontsize=5):
    """
    plot of text stats
    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def produce_strings(wave,flux,fluxerr):
        mea_f          = np.mean(flux[np.isfinite(flux)])
        med_f          = np.median(flux[np.isfinite(flux)])
        max_f          = np.max(flux[np.isfinite(flux)])
        max_f_wave     = wave[np.where(flux == max_f)[0]][0]

        s2n            = flux/fluxerr
        mea_s2n        = np.mean(s2n[np.isfinite(s2n)])
        med_s2n        = np.median(s2n[np.isfinite(s2n)])
        max_s2n        = np.max(s2n[np.isfinite(s2n)])
        max_s2n_wave   = wave[np.where(s2n == max_s2n)[0]][0]

        numberstring   = '\\verb+'+str("%12.2f" % mea_f)+'+'+\
                         '\\verb+'+str("%18.2f" % med_f)+'+'+\
                         '\\verb+'+str("%25.2f" % max_f)+'+'+\
                         '\\verb+'+str("%20.2f"  % max_f_wave)+'+'\
                         '\\verb+'+str("%10.2f" % mea_s2n)+'+'+\
                         '\\verb+'+str("%12.2f" % med_s2n)+'+'+\
                         '\\verb+'+str("%15.2f" % max_s2n)+'+'+\
                         '\\verb+'+str("%17.2f"  % max_s2n_wave)+'+'

        return numberstring
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ax.text(45.0,110.0,title,horizontalalignment='center',verticalalignment='center',fontsize=fontsize)

    colstring = '$<$f/[1e-20cgs]$>$ \\verb+  + f\_med/[1e-20cgs] \\verb+  + ' \
                'f\_max/[1e-20cgs] \\verb+  + $\\lambda$(f\_max) \\verb+  +'+\
                '$<$S/N$>$ \\verb+  + S/N\_med \\verb+  + S/N\_max \\verb+  + $\\lambda$(S/N\_max)'

    ax.text(8,85.0,colstring,horizontalalignment='left',verticalalignment='center',fontsize=fontsize)
    # --------------------------- TDOSE ---------------------------
    ax.text(-8.0,60.0,'('+colors[0]+')',horizontalalignment='left',verticalalignment='center',fontsize=fontsize,
            color=colors[0])

    numberstring = produce_strings(t_wave,t_flux,t_fluxerr)
    ax.text(0.0,60.0,'TDOSE:\\verb+  +'+numberstring,horizontalalignment='left',verticalalignment='center',fontsize=fontsize)

    # --------------------------- CARUANA ---------------------------
    ax.text(-8.0,35.0,'('+colors[1]+')',horizontalalignment='left',verticalalignment='center',fontsize=fontsize,
            color=colors[1])

    numberstring = produce_strings(c_wave,c_flux,c_fluxerr)
    ax.text(0.0,35.0,'Caruana:\\verb+ +'+numberstring,horizontalalignment='left',verticalalignment='center',fontsize=fontsize)

    # --------------------------- URRUTIA ---------------------------
    ax.text(-8.0,10.0,'('+colors[2]+')',horizontalalignment='left',verticalalignment='center',fontsize=fontsize,
            color=colors[2])
    if u_flux is not None:
        numberstring = produce_strings(u_wave,u_flux,u_fluxerr)
        ax.text(0.0,10.0,'Urrutia:\\verb+  +'+numberstring,horizontalalignment='left',verticalalignment='center',fontsize=fontsize)
        ax.text(0.0,-5.0,'MUSE-Wide ID = '+str(int(MWid[0]))+' and Guo Rmatch = '+str("%5.2f" % u_rmatch)+'"',
                horizontalalignment='left',verticalalignment='center',fontsize=fontsize)
    else:
        str_urr = 'Urrutia:\\verb+         +No match between MUSE-Wide and Guo catalog'
        ax.text(0.0,10.0,str_urr,horizontalalignment='left',verticalalignment='center',fontsize=fontsize)

    ax.set_ylim([0,100])
    ax.set_xlim([0,100])

    ax.axis('off')
    ax.set_xlabel(' ')
    ax.set_ylabel(' ')
    ax.set_xticks([])
    ax.set_yticks([])

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_compare_spec(ax,spectitle,
                     t_wave,t_flux,t_fluxerr,c_wave,c_flux,c_fluxerr,u_wave,u_flux,u_fluxerr,
                     xrange=[4800,9300],markwave=5500,drawboxes=False,fontsize=6,colors=['red','blue','green'],
                     labels=['t','c','u'],plotSNcurve=False,shownoise=True,lthick=2,fillalpha=0.30,verbose=True):
    """
    Plotting commands for compare_TDOSEspecs() plotting
    """
    ax.set_title(spectitle,fontsize=fontsize, loc='left',y=0.9)

    if markwave is not None:
        diff    = np.abs(np.asarray(t_wave)-markwave)
        markent = np.where(diff == np.min(diff))[0]
        if len(markent) == 0:
            if verbose: print ' WARNING the "markwave" is outside provided plotting xrange '
            markwave = False

        if plotSNcurve:
            t_s2n = t_flux/t_fluxerr
            ax.plot(t_wave,t_s2n,color=colors[0],lw=lthick)
            yrangeplot = ax.get_ylim()

            if c_flux is not None:
                c_s2n = c_flux/c_fluxerr
                ax.plot(c_wave,c_s2n,color=colors[1],lw=lthick)

            if u_flux is not None:
                u_s2n = u_flux/u_fluxerr
                ax.plot(u_wave,u_s2n,color=colors[2],lw=lthick)

            if markwave:
                ax.plot([t_wave[markent],t_wave[markent]],ax.get_ylim(),color='magenta',linestyle='--',lw=lthick)
            ylabel = 'S/N'
        else:
            ax.plot(t_wave,t_flux,color=colors[0],lw=lthick)
            yrangeplot = ax.get_ylim()
            if c_flux is not None:
                ax.plot(c_wave,c_flux,color=colors[1],lw=lthick)
            if u_flux is not None:
                ax.plot(u_wave,u_flux,color=colors[2],lw=lthick)

            if shownoise:
                plt.fill_between(t_wave,t_flux-t_fluxerr,t_flux+t_fluxerr,alpha=fillalpha,color=colors[0])
                if c_flux is not None:
                    plt.fill_between(c_wave,c_flux-c_fluxerr,c_flux+c_fluxerr,alpha=fillalpha,color=colors[1])
                if u_flux is not None:
                    plt.fill_between(u_wave,u_flux-u_fluxerr,u_flux+u_fluxerr,alpha=fillalpha,color=colors[2])

            if markwave:
                ax.plot([t_wave[markent],t_wave[markent]],ax.get_ylim(),color='magenta',linestyle='--',lw=lthick)
            ylabel = 'Flux [1e-20cgs]'


        ax.set_xlabel('Wavelength [\AA]', fontsize=fontsize)
        ax.set_ylabel(ylabel,fontsize=fontsize)

        ax.set_xlim(xrange)
        ax.set_ylim(yrangeplot)

        ax.plot(xrange,[0,0],'-',color='gray',lw=0.5)

        if drawboxes:
            for drawbox in drawboxes:
                plt.fill_between(drawbox,[yrangeplot[0],yrangeplot[0]],[yrangeplot[1],yrangeplot[1]],
                         color='magenta',alpha=fillalpha)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =