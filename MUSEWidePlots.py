# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import os
import glob
import pyfits
import kbsutilities as kbs
import numpy as np
import MUSEWidePlots as mwp
import MUSEWideUtilities as mwu
import sys
import matplotlib.pyplot as plt
import collections
import MiGs
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
    mwp.compare_TDOSEspecs(skipobj=True)

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
        statdic = mwu.get_specstat(wave,flux,fluxerr)
        numberstring   = '\\verb+'+str("%12.2f" % statdic['mea_f'])+'+'+\
                         '\\verb+'+str("%18.2f" % statdic['med_f'])+'+'+\
                         '\\verb+'+str("%25.2f" % statdic['max_f'])+'+'+\
                         '\\verb+'+str("%20.2f"  % statdic['max_f_wave'])+'+'\
                         '\\verb+'+str("%10.2f" % statdic['mea_s2n'])+'+'+\
                         '\\verb+'+str("%12.2f" % statdic['med_s2n'])+'+'+\
                         '\\verb+'+str("%15.2f" % statdic['max_s2n'])+'+'+\
                         '\\verb+'+str("%17.2f"  % statdic['max_s2n_wave'])+'+'

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

            if markwave:
                ax.plot([t_wave[markent],t_wave[markent]],ax.get_ylim(),color='magenta',linestyle='--',lw=lthick)

            if u_flux is not None:
                u_s2n = u_flux/u_fluxerr
                ax.plot(u_wave,u_s2n,color=colors[2],lw=lthick)

            if c_flux is not None:
                c_s2n = c_flux/c_fluxerr
                ax.plot(c_wave,c_s2n,color=colors[1],lw=lthick)

            ax.plot(t_wave,t_s2n,color=colors[0],lw=lthick)

            ylabel = 'S/N'
        else:
            ax.plot(t_wave,t_flux,color=colors[0],lw=lthick) # plot to get yrange
            yrangeplot = ax.get_ylim()

            if markwave:
                ax.plot([t_wave[markent],t_wave[markent]],ax.get_ylim(),color='magenta',linestyle='--',lw=lthick)

            if shownoise:
                if u_flux is not None:
                    plt.fill_between(u_wave,u_flux-u_fluxerr,u_flux+u_fluxerr,alpha=fillalpha,color=colors[2])
                if c_flux is not None:
                    plt.fill_between(c_wave,c_flux-c_fluxerr,c_flux+c_fluxerr,alpha=fillalpha,color=colors[1])
                plt.fill_between(t_wave,t_flux-t_fluxerr,t_flux+t_fluxerr,alpha=fillalpha,color=colors[0])

            if u_flux is not None:
                ax.plot(u_wave,u_flux,color=colors[2],lw=lthick)
            if c_flux is not None:
                ax.plot(c_wave,c_flux,color=colors[1],lw=lthick)
            ax.plot(t_wave,t_flux,color=colors[0],lw=lthick) # plot in front

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
def plot_comparison_stats(tdosespecssearch='/Volumes/DATABCKUP2/TDOSEextractions/spec_comparisons/tdose_cdfs-04_170630/*fits',
                          caruanaspecssearch='/Volumes/DATABCKUP2/TDOSEextractions/spec_comparisons/caruana_cdfs-04/*',
                          colors=['black','red'],labels=['TDOSE Gauss','Caruana'],wavetolerance=3,
                          outdir='./',verbose=True):
    """
    Plotting commands for compare_TDOSEspecs() plotting

    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp
    mwp.plot_comparison_stats(wavetolerance=3,outdir='/Volumes/DATABCKUP2/TDOSEextractions/spec_comparisons/tdose_cdfs-04_170630/')
    """

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loaking for spectra by globbing using the input path-strings'
    tdosespecs   = glob.glob(tdosespecssearch)
    caruanaspecs = glob.glob(caruanaspecssearch)

    Ntdose =  len(tdosespecs)
    if Ntdose == 0:
        sys.exit('No tdose spectra found in '+tdosespecssearch)
    else:
        if verbose: print ' - Will collect and plot stats for the '+str(Ntdose)+' TDOSE spectra found '

    if len(caruanaspecs) == 0:
        sys.exit(' No Caruana spectra found in '+caruanaspecssearch)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Putting together lists of IDs'
    if verbose: print '   Generate ID list for TDOSE spectra'
    ids_tdose = [float(tt.split('/')[-1].split('_')[-1].split('.')[0]) for tt in tdosespecs]
    ids_tdose = np.asarray(ids_tdose)

    if verbose: print '   Generate ID list for Caruana spectra'
    ids_car = [float(tt.split('/')[-1].split('_')[-1].split('.')[0]) for tt in caruanaspecs]
    ids_car = np.asarray(ids_car)


    if len(ids_tdose) != len(ids_car):
        sys.exit(' The number of spectra (IDs) do not match between the TDOSE ('+
                 str(len(ids_tdose))+') and the Caruana ('+str(len(ids_car))+') spectra')

    t_mea_f      = []
    t_med_f      = []
    t_max_f      = []
    t_std_f      = []
    t_max_f_wave = []
    t_mea_s2n      = []
    t_med_s2n      = []
    t_max_s2n      = []
    t_std_s2n      = []
    t_max_s2n_wave = []

    c_mea_f      = []
    c_med_f      = []
    c_max_f      = []
    c_std_f      = []
    c_max_f_wave = []
    c_mea_s2n      = []
    c_med_s2n      = []
    c_max_s2n      = []
    c_std_s2n      = []
    c_max_s2n_wave = []

    plotname = outdir+'spec_stat_histograms.pdf'
    for ii, tspec in enumerate(tdosespecs):
        dat_tdose    = pyfits.open(tspec)[1].data
        t_flux    = dat_tdose['flux']
        t_fluxerr = dat_tdose['fluxerror']
        t_wave    = dat_tdose['wave']
        t_statdic = mwu.get_specstat(t_wave,t_flux,t_fluxerr)

        t_mea_f.append(t_statdic['mea_f'])
        t_med_f.append(t_statdic['med_f'])
        t_std_f.append(t_statdic['std_f'])
        t_mea_s2n.append(t_statdic['mea_s2n'])
        t_med_s2n.append(t_statdic['med_s2n'])
        t_std_s2n.append(t_statdic['std_s2n'])

        cent   = np.where(ids_car == ids_tdose[ii])[0]
        if len(cent) == 1:
            c_flux    = pyfits.open(caruanaspecs[cent[0]])[0].data
            c_fluxerr = pyfits.open(caruanaspecs[cent[0]])[1].data
            c_wave    = t_wave
            c_statdic = mwu.get_specstat(c_wave,c_flux,c_fluxerr)

            c_mea_f.append(c_statdic['mea_f'])
            c_med_f.append(c_statdic['med_f'])
            c_std_f.append(c_statdic['std_f'])
            c_mea_s2n.append(c_statdic['mea_s2n'])
            c_med_s2n.append(c_statdic['med_s2n'])
            c_std_s2n.append(c_statdic['std_s2n'])

            if np.abs(c_statdic['max_f_wave'] - t_statdic['max_f_wave']) < wavetolerance:
                t_max_f.append(t_statdic['max_f'])
                t_max_f_wave.append(t_statdic['max_f_wave'])

                c_max_f.append(c_statdic['max_f'])
                c_max_f_wave.append(c_statdic['max_f_wave'])

            if np.abs(c_statdic['max_s2n_wave'] - t_statdic['max_s2n_wave']) < wavetolerance:
                t_max_s2n.append(t_statdic['max_s2n'])
                t_max_s2n_wave.append(t_statdic['max_s2n_wave'])

                c_max_s2n.append(c_statdic['max_s2n'])
                c_max_s2n_wave.append(c_statdic['max_s2n_wave'])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Ncol         = 4
    Nrow         = 2
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=(Ncol, Nrow))
    fig.subplots_adjust(wspace=0.5, hspace=0.5,left=0.09, right=0.98, bottom=0.15, top=0.96)
    Fsize  = 4
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('/')[-1].replace('_','\_'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #---------------------------- f mean hist ----------------------------
    rowval  = 0
    rowspan = 1
    colval  = 0
    colspan = 1
    ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

    bins = np.arange(-5,50,2)
    mwp.plot_comp_hist(ax,'Mean Flux',t_mea_f,c_mea_f,colors,Fsize,lthick,bins)

    #---------------------------- f median hist ----------------------------
    rowval  = 0
    rowspan = 1
    colval  = 1
    colspan = 1
    ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

    mwp.plot_comp_hist(ax,'Median Flux',t_med_f,c_med_f,colors,Fsize,lthick,bins)

    #---------------------------- f std hist ----------------------------
    rowval  = 0
    rowspan = 1
    colval  = 2
    colspan = 1
    ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

    mwp.plot_comp_hist(ax,'Standard Deviation Flux',t_std_f,c_std_f,colors,Fsize,lthick,bins)

    #---------------------------- f max hist ----------------------------
    rowval  = 0
    rowspan = 1
    colval  = 3
    colspan = 1
    ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

    bins = np.arange(100,600,20)
    mwp.plot_comp_hist(ax,'Max Flux',t_max_f,c_max_f,colors,Fsize,lthick,bins)

    #---------------------------- S/N mean hist ----------------------------
    rowval  = 1
    rowspan = 1
    colval  = 0
    colspan = 1
    ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

    bins = np.arange(0,2,0.05)
    mwp.plot_comp_hist(ax,'Mean S/N',t_mea_s2n,c_mea_s2n,colors,Fsize,lthick,bins)

    #---------------------------- S/N median hist ----------------------------
    rowval  = 1
    rowspan = 1
    colval  = 1
    colspan = 1
    ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

    mwp.plot_comp_hist(ax,'Median S/N',t_med_s2n,c_med_s2n,colors,Fsize,lthick,bins)

    #---------------------------- f std hist ----------------------------
    rowval  = 1
    rowspan = 1
    colval  = 2
    colspan = 1
    ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

    mwp.plot_comp_hist(ax,'Standard Deviation S/N',t_std_s2n,c_std_s2n,colors,Fsize,lthick,bins)

    #---------------------------- S/N max hist ----------------------------
    rowval  = 1
    rowspan = 1
    colval  = 3
    colspan = 1
    ax = plt.subplot2grid((Nrow,Ncol), (rowval, colval), colspan=colspan, rowspan=rowspan)

    bins = np.arange(-5,50,2)
    mwp.plot_comp_hist(ax,'Max S/N',t_max_s2n,c_max_s2n,colors,Fsize,lthick,bins)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #if verbose: print ' - Saving plot to',plotname
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_comp_hist(ax,xlabel,hist1,hist2,colors,fontsize,lthick,bins):
    """
    Plotting for plot_comparison_stats()
    """
    #ax.set_title(title,fontsize=fontsize, loc='left',y=0.9)

    ax.hist(hist1,bins=bins,color=colors[0],lw=lthick,histtype='step')
    ax.hist(hist2,bins=bins,color=colors[1],lw=lthick,histtype='step')

    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel('\#'  , fontsize=fontsize)

    ax.plot(ax.get_xlim(),[0,0],'-',color='gray',lw=0.5)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview(spectra, labels, wavecols, fluxcols, fluxerrcols, redshift, voffset=0,
                        outputfigure='default', plotSN=False, skyspectra=None, wavecols_sky=None, fluxcols_sky=None,
                        yrangefull=None, verbose=True):
    """

    Plotting overview with zoom-ins of 1D spectrum.

    Genralization of the more restrictive plotting script
    ciiiEmitterCandidates.plot_MUSElya_forsample()

    This can be used to plot multiple spectra (e.g., PSF weighted 1D MUSE, TDOSE 1D MUSE and 3D-HST 1D spectra) of
    individual objects useful for line searches and identifications.

    --- INPUT ---
    spectra       List of spectra to include in plot
    labels        Labels to use for spectra
    wavecols      Column names of entries in "spectra" files containing wavelengths     [in angstrom]
    fluxcols      Column names of entries in "spectra" files containing fluxes          [1e-20 cgs]
    fluxerrcols   Column names of entries in "spectra" files containing flux errors     [1e-20 cgs]
    redshift      Redshift to position emission line markes at
    voffset       Velocity offset [km/s] wrt. to the emission line markers to mark. Defaults to offset of Lya line
                  assuming the redshift provided is the Lya redshift.
    outputfigure  Path and name of figure to generate
    plotSN        If true a S/N version (flux/fluxerr instead of flux) of the figure will be generated
    skyspectra    To show the estimated sky provide "sky spectra" for each of the "spectra" provided.
                  Use None if no sky spectrum should be plotted
    wavecols_sky  Columns containing wavelengts in sky spectra
    fluxcols_sky  Columns containing flux in sky spectra

    yrangefull    To fix the y-range for the overview panel including all spectra set it with this keyword
    verbose       Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp

    datapath     = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/TDOSEoutput171019_106003018/'
    spectra      = [datapath+'tdose_spectrum_candels-cdfs-06_modelimg_0106003018-0106003018.fits',datapath+'spectrum_10603018.fits',datapath+'spectrum_10603018_goodss-02-G141_22684_MiG1Dreformat.fits']
    labels       = ['TDOSE','PSFext','3D-HST']
    wavecols     = ['wave','WAVE_AIR','WAVE_AIR']
    fluxcols     = ['flux','FLUX','FLUX']
    fluxerrcols  = ['fluxerror','FLUXERR','FLUXERR']
    zLya         = 2.97756004333
    voffset      = 235
    skyspectra   = ['/Users/kschmidt/work/MUSE/skyspectra/SKY_SPECTRUM_candels-cdfs-06_av.fits',None,None]
    wavecols_sky = ['lambda',None,None]
    fluxcols_sky = ['data',None,None]
    yrangefull   = [-1000,2000]

    mwp.plot_1DspecOverview(spectra, labels, wavecols, fluxcols, fluxerrcols, zLya, voffset=voffset, skyspectra=skyspectra, wavecols_sky=wavecols_sky, fluxcols_sky=fluxcols_sky, outputfigure='default', yrangefull=yrangefull, plotSN=False)

    """
    voffset = -1.0 * voffset
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if outputfigure == 'default':
        plotpath = '/'.join(spectra[0].split('/')[:-1])
        specfigure = plotpath+'/overview_1DspecWzooms_z'+str(redshift).replace('.','p')+\
                     '_voffset'+str(voffset).replace('.','p')+'.pdf'
    else:
        specfigure = outputfigure

    if plotSN:
        specfigure = specfigure.replace('.pdf','_SN.pdf')
    if verbose: print ' - 1D overview figure will be saved to:\n   '+specfigure
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    redshiftplot = redshift + (voffset*(redshift+1.0) / 299792.458)
    if verbose: print ' - Will plot emission line markers using at redshift '+str("%.6f" % redshift)+\
                      ' (z~'+str("%.6f" % redshiftplot)+' including (Lya) velocity offset of '+str(voffset)+'km/s)'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Nspec = len(spectra)
    if verbose: print ' - Loading the '+str(Nspec)+' spectra provided for plotting'

    datadic = collections.OrderedDict()

    for ss, specname in enumerate(spectra):
        spec_dat      = pyfits.open(specname)[1].data
        spec_wave     = spec_dat[wavecols[ss]]
        spec_flux     = spec_dat[fluxcols[ss]]
        spec_ferr     = spec_dat[fluxerrcols[ss]]
        spec_S2N      = spec_flux/spec_ferr
        spec_filllow  = spec_flux-spec_ferr
        spec_fillhigh = spec_flux+spec_ferr
        spec_wavecov  = [np.min(spec_wave),np.max(spec_wave)]

        if skyspectra[ss] is not None:
            spec_dat      = pyfits.open(skyspectra[ss])[1].data
            spec_wave_sky = spec_dat[wavecols_sky[ss]]
            spec_flux_sky = spec_dat[fluxcols_sky[ss]]
        else:
            spec_wave_sky = None
            spec_flux_sky = None

        datadic[specname] = {'spec_wave':spec_wave,
                             'spec_flux':spec_flux,
                             'spec_ferr':spec_ferr,
                             'spec_S2N':spec_S2N,
                             'spec_filllow':spec_filllow,
                             'spec_fillhigh':spec_fillhigh,
                             'spec_wavecov':spec_wavecov,
                             'spec_wave_sky':spec_wave_sky,
                             'spec_flux_sky':spec_flux_sky}

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Defining quanteties to plot; checking wavelength coverage of lines given redshift'
    llistdic      = MiGs.linelistdic(listversion='full') # loading line list for plots

    plot_Lyg      = False; checkwave = 973.*(1+redshift)

    plot_Lyb      = False; checkwave = 1026.*(1+redshift)
    for specname in spectra:
        wavecov = datadic[specname]['spec_wavecov']
        if (checkwave > wavecov[0]) & (checkwave < wavecov[1]):
            plot_Lyb = True

    plot_Lya      = True

    plot_CII      = False; checkwave = 1335.*(1+redshift)
    for specname in spectra:
        wavecov = datadic[specname]['spec_wavecov']
        if (checkwave > wavecov[0]) & (checkwave < wavecov[1]):
            plot_CII      = True

    plot_SiIVOIV  = False; checkwave = 1400.*(1+redshift)
    for specname in spectra:
        wavecov = datadic[specname]['spec_wavecov']
        if (checkwave > wavecov[0]) & (checkwave < wavecov[1]):
            plot_SiIVOIV  = True

    plot_CIV      = False; checkwave = 1549.*(1+redshift)
    for specname in spectra:
        wavecov = datadic[specname]['spec_wavecov']
        if (checkwave > wavecov[0]) & (checkwave < wavecov[1]):
            plot_CIV      = True

    plot_HeII     = False; checkwave = 1640.*(1+redshift)
    for specname in spectra:
        wavecov = datadic[specname]['spec_wavecov']
        if (checkwave > wavecov[0]) & (checkwave < wavecov[1]):
            plot_HeII      = True

    plot_CIII     = False; checkwave = 1908.*(1+redshift)
    for specname in spectra:
        wavecov = datadic[specname]['spec_wavecov']
        if (checkwave > wavecov[0]) & (checkwave < wavecov[1]):
            plot_CIII     = True

    plot_CIIb     = False; checkwave = 2326.*(1+redshift)

    plot_MgII     = False; checkwave = 2795.*(1+redshift)
    for specname in spectra:
        wavecov = datadic[specname]['spec_wavecov']
        if (checkwave > wavecov[0]) & (checkwave < wavecov[1]):
            plot_MgII     = True

    plot_OII      = False; checkwave = 3727.5*(1+redshift)
    for specname in spectra:
        wavecov = datadic[specname]['spec_wavecov']
        if (checkwave > wavecov[0]) & (checkwave < wavecov[1]):
            plot_OII      = True

    plot_Hd       = False; checkwave = 4101.74*(1+redshift)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Plotting figure '
    figuresize_x = 13
    figuresize_y = 10
    fig          = plt.figure(figsize=(figuresize_x,figuresize_y))
    Fsize        = 10
    LW           = 2
    plt.rc('text', usetex=True)                         # enabling LaTex rendering of text
    plt.rc('font', family='serif',size=Fsize)           # setting text font
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    left   = 0.06   # the left side of the subplots of the figure
    right  = 0.98   # the right side of the subplots of the figure
    bottom = 0.05   # the bottom of the subplots of the figure
    top    = 0.98   # the top of the subplots of the figure
    wspace = 0.20   # the amount of width reserved for blank space between subplots
    hspace = 0.20   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    speccols          = ['blue','green','red','magenta','cyan','orange']
    #speccol           = 'blue'
    xlabel            = '$\lambda$ / [\AA]'
    ylabel            = '$f_\lambda / [10^{-20}$erg/s/cm$^2$/\\AA]'
    if plotSN:
        ylabel        = 'S/N'
    col_linemarker    = 'gray'
    wavescale         = 1.0

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_Lyb:
        plt.subplot(4, 3, 1) # Lybeta+OIV
        windowcenter = 1030.0
        windowwidth  = 20.0
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        yrange = mwp.plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,
                                                   windowcenter,redshift,xrange,plotSN=plotSN,labels=None)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        mwp.plot_1DspecOverview_plotlines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,
                                          xrange,yrange,redshift,LW,wavetype='vac')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_Lyb = xrange
        yrange_Lyb = yrange
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.subplot(4, 3, 1)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nLy$\\beta$ + OVI doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_Lya:
        plt.subplot(4, 3, 2) # Lyalpha+NV
        windowcenter = 1229.0
        windowwidth  = 28.0
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        yrange = mwp.plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,
                                                   windowcenter,redshift,xrange,plotSN=plotSN,labels=None)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        mwp.plot_1DspecOverview_plotlines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,
                                          xrange,yrange,redshift,LW,wavetype='vac')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_Lya = xrange
        yrange_Lya = yrange
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.subplot(4, 3, 2)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nLy$\\alpha$ + NV doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_CII:
        plt.subplot(4, 3, 3) # CII
        windowcenter = 1335.
        windowwidth  = 10.
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        yrange = mwp.plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,
                                                   windowcenter,redshift,xrange,plotSN=plotSN,labels=None)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        mwp.plot_1DspecOverview_plotlines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,
                                          xrange,yrange,redshift,LW,wavetype='vac')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_CII = xrange
        yrange_CII = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.subplot(4, 3, 3)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nCII$\\lambda$1336',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_SiIVOIV:
        plt.subplot(4, 3, 4) # SiIV+OIV
        windowcenter = 1397.0
        windowwidth  = 13
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        yrange = mwp.plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,
                                                   windowcenter,redshift,xrange,plotSN=plotSN,labels=None)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        mwp.plot_1DspecOverview_plotlines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,
                                          xrange,yrange,redshift,LW,wavetype='vac')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_SiIVOIV = xrange
        yrange_SiIVOIV = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.subplot(4, 3, 4)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nSiIV and OIV] doublets',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_CIV:
        plt.subplot(4, 3, 5) # CIV
        windowcenter = 1549.
        windowwidth  = 10.
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        yrange = mwp.plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,
                                                   windowcenter,redshift,xrange,plotSN=plotSN,labels=None)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        mwp.plot_1DspecOverview_plotlines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,
                                          xrange,yrange,redshift,LW,wavetype='vac')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_CIV = xrange
        yrange_CIV = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.subplot(4, 3, 5)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nCIV doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_HeII:
        plt.subplot(4, 3, 6) # HeII
        windowcenter = 1653.0
        windowwidth  = 26
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        yrange = mwp.plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,
                                                   windowcenter,redshift,xrange,plotSN=plotSN,labels=None)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        mwp.plot_1DspecOverview_plotlines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,
                                          xrange,yrange,redshift,LW,wavetype='vac')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_HeII = xrange
        yrange_HeII = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.subplot(4, 3, 6)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nHeII and OIII] doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_CIII:
        plt.subplot(4, 3, 7) # CIII]
        windowcenter = 1908.0
        windowwidth  = 10
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        yrange = mwp.plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,
                                                   windowcenter,redshift,xrange,plotSN=plotSN,labels=None)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        mwp.plot_1DspecOverview_plotlines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,
                                          xrange,yrange,redshift,LW,wavetype='vac')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_CIII = xrange
        yrange_CIII = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.subplot(4, 3, 7)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nCIII] doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_MgII:
        plt.subplot(4, 3, 8) # MgII
        windowcenter = 2795.
        windowwidth  = 10
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        yrange = mwp.plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,
                                                   windowcenter,redshift,xrange,plotSN=plotSN,labels=None)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        mwp.plot_1DspecOverview_plotlines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,
                                          xrange,yrange,redshift,LW,wavetype='vac')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_MgII = xrange
        yrange_MgII = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.subplot(4, 3, 8)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nMgII doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_OII:
        plt.subplot(4, 3, 9) # [OII]
        windowcenter = 3727.5
        windowwidth  = 10
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        yrange = mwp.plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,
                                                   windowcenter,redshift,xrange,plotSN=plotSN,labels=None)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        mwp.plot_1DspecOverview_plotlines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,
                                          xrange,yrange,redshift,LW,wavetype='vac')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_OII = xrange
        yrange_OII = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.subplot(4, 3, 9)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\n[OII] doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.subplot(4, 3, (10,12)) # Full specs
    windowcenter = 7500.0
    xrangefull   = [4500,18000]
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    yrange = mwp.plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,
                                               windowcenter,redshift,xrangefull,plotSN=plotSN,labels=labels)
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    if not plotSN:
        if yrangefull is None:
            yrangefull  = yrange
    else:
        yrangefull  = yrange
    Dyrangefull = yrangefull[1]-yrangefull[0]
    # --- "ZOOM BOXES" ---
    if plot_Lyb:
        plt.plot(xrange_Lyb,np.zeros(2)+yrange_Lyb[0],'-',color='black',lw=LW)
        plt.plot(xrange_Lyb,np.zeros(2)+yrange_Lyb[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_Lyb[0],yrange_Lyb,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_Lyb[1],yrange_Lyb,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_Lyb)),yrange_Lyb[1]+0.03*Dyrangefull,'Ly$\\beta$+OVI',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_Lya:
        plt.plot(xrange_Lya,np.zeros(2)+yrange_Lya[0],'-',color='black',lw=LW)
        plt.plot(xrange_Lya,np.zeros(2)+yrange_Lya[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_Lya[0],yrange_Lya,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_Lya[1],yrange_Lya,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_Lya)),yrange_Lya[1]+0.03*Dyrangefull,'Ly$\\alpha$+NV',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_CII:
        plt.plot(xrange_CII,np.zeros(2)+yrange_CII[0],'-',color='black',lw=LW)
        plt.plot(xrange_CII,np.zeros(2)+yrange_CII[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CII[0],yrange_CII,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CII[1],yrange_CII,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_CII)),yrange_CII[1]+0.03*Dyrangefull,'CII',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_SiIVOIV:
        plt.plot(xrange_SiIVOIV,np.zeros(2)+yrange_SiIVOIV[0],'-',color='black',lw=LW)
        plt.plot(xrange_SiIVOIV,np.zeros(2)+yrange_SiIVOIV[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_SiIVOIV[0],yrange_SiIVOIV,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_SiIVOIV[1],yrange_SiIVOIV,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_SiIVOIV)),yrange_SiIVOIV[1]+0.03*Dyrangefull,'SiIV+OIV]',rotation='vertical',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_CIV:
        plt.plot(xrange_CIV,np.zeros(2)+yrange_CIV[0],'-',color='black',lw=LW)
        plt.plot(xrange_CIV,np.zeros(2)+yrange_CIV[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CIV[0],yrange_CIV,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CIV[1],yrange_CIV,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_CIV)),yrange_CIV[1]+0.03*Dyrangefull,'CIV',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_HeII:
        plt.plot(xrange_HeII,np.zeros(2)+yrange_HeII[0],'-',color='black',lw=LW)
        plt.plot(xrange_HeII,np.zeros(2)+yrange_HeII[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_HeII[0],yrange_HeII,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_HeII[1],yrange_HeII,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_HeII)),yrange_HeII[1]+0.03*Dyrangefull,'HeII+OIII]',rotation='vertical',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_CIII:
        plt.plot(xrange_CIII,np.zeros(2)+yrange_CIII[0],'-',color='black',lw=LW)
        plt.plot(xrange_CIII,np.zeros(2)+yrange_CIII[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CIII[0],yrange_CIII,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CIII[1],yrange_CIII,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_CIII)),yrange_CIII[1]+0.03*Dyrangefull,'CIII]',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    # if plot_CIIb:
    #     plt.plot(xrange_CIIb,np.zeros(2)+yrange_CIIb[0],'-',color='black',lw=LW)
    #     plt.plot(xrange_CIIb,np.zeros(2)+yrange_CIIb[1],'-',color='black',lw=LW)
    #     plt.plot(np.zeros(2)+xrange_CIIb[0],yrange_CIIb,'-',color='black',lw=LW)
    #     plt.plot(np.zeros(2)+xrange_CIIb[1],yrange_CIIb,'-',color='black',lw=LW)
    #     plt.text(np.mean(np.asarray(xrange_CIIb)),yrange_CIIb[1]+0.03*Dyrangefull,'CII]',
    #              color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_MgII:
        plt.plot(xrange_MgII,np.zeros(2)+yrange_MgII[0],'-',color='black',lw=LW)
        plt.plot(xrange_MgII,np.zeros(2)+yrange_MgII[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_MgII[0],yrange_MgII,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_MgII[1],yrange_MgII,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_MgII)),yrange_MgII[1]+0.03*Dyrangefull,'MgII',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_OII:
        plt.plot(xrange_OII,np.zeros(2)+yrange_OII[0],'-',color='black',lw=LW)
        plt.plot(xrange_OII,np.zeros(2)+yrange_OII[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_OII[0],yrange_OII,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_OII[1],yrange_OII,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_OII)),yrange_OII[1]+0.03*Dyrangefull,'[OII]',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    # if plot_Hd:
    #     plt.plot(xrange_Hd,np.zeros(2)+yrange_Hd[0],'-',color='black',lw=LW)
    #     plt.plot(xrange_Hd,np.zeros(2)+yrange_Hd[1],'-',color='black',lw=LW)
    #     plt.plot(np.zeros(2)+xrange_Hd[0],yrange_Hd,'-',color='black',lw=LW)
    #     plt.plot(np.zeros(2)+xrange_Hd[1],yrange_Hd,'-',color='black',lw=LW)
    #     plt.text(np.mean(np.asarray(xrange_Hd)),yrange_Hd[1]+0.03*Dyrangefull,'H$\delta$',
    #              color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    #--------- LEGEND ---------
    anchorpos = (0.5, 1.2)
    leg = plt.legend(fancybox=True,numpoints=1, loc='upper center',prop={'size':Fsize},ncol=len(labels))#,
                     #bbox_to_anchor=anchorpos)  # add the legend
    leg.get_frame().set_alpha(0.7)
    #--------------------------


    plt.xlim(xrangefull)
    plt.ylim(yrangefull)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.savefig(specfigure, dpi=300) # dpi = dot per inch for rasterized points
    plt.clf()
    plt.close('all')
    if verbose: print ' - Saved figure to ',specfigure
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,colors,windowcenter,redshift,xrange,plotSN=False,labels=None):
    """

    """
    if labels is None:
        labels = [None]*len(spectra)
    yrangecomb = [0,10]

    for ss, specname in enumerate(spectra):
        waveent = (datadic[specname]['spec_wave'] > xrange[0]) & (datadic[specname]['spec_wave'] < xrange[1])

        if plotSN:
            plt.plot(datadic[specname]['spec_wave'], datadic[specname]['spec_S2N'], '-',
                     alpha=0.8,color=colors[ss],label=labels[ss])

            try:
                fluxmin = np.min(np.asarray([0,  np.min(datadic[specname]['spec_S2N'][waveent]) ]))
                fluxmax = np.max(np.asarray([10, np.max(datadic[specname]['spec_S2N'][waveent]) ]))
                yrange  = [fluxmin,1.1*fluxmax]
            except:
                yrange = [0,10]

        else:
            plt.plot(datadic[specname]['spec_wave'], datadic[specname]['spec_flux'], '-',
                     alpha=0.8,color=colors[ss],label=labels[ss])

            plt.fill_between(datadic[specname]['spec_wave'],datadic[specname]['spec_filllow'],datadic[specname]['spec_fillhigh'],
                             alpha=0.20,color=colors[ss])

            try:
                fluxmin = np.min(np.asarray([0,  np.min(datadic[specname]['spec_flux'][waveent]) ]))
                fluxmax = np.max(np.asarray([10, np.max(datadic[specname]['spec_flux'][waveent]) ]))
                yrange  = [fluxmin,1.1*fluxmax]
            except:
                yrange = [0,500]

        if skyspectra[ss] is not None:
            sky_w   = datadic[specname][wavecols_sky[ss]]
            sky_f   = datadic[specname][fluxcols_sky[ss]]
            skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]))[0]
            skywave = sky_w[skyent]
            skylow  = np.zeros(len(skywave))
            skyflux = sky_f[skyent]
            skyhigh = skyflux
            #skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])

            plt.fill_between(skywave,skylow+yrange[0],skyhigh+yrange[0],alpha=0.3,color='black')

        if yrange[0] < yrangecomb[0]:
            yrangecomb[0] = yrange[0]

        if yrange[1] > yrangecomb[1]:
            yrangecomb[1] = yrange[1]

    return yrangecomb
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview_plotlines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,LW,wavetype='vac'):
    """

    """
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    for ll in llistdic.keys():
        linedat      = llistdic[ll]
        linename     = linedat[0]
        if wavetype == 'vac':
            linewave = linedat[1]
        else:
            linewave = kbs.convert_wavelength(linedat[1],version='vac2air')
            pdb.set_trace()
        horalign     = linedat[2]
        lineposition = linewave*(redshift+1.0)/wavescale

        if (lineposition > xrange[0]) & (lineposition < xrange[1]):
            plt.plot(np.zeros(2)+lineposition,yrange,color=col_linemarker,alpha=0.7,linestyle='-',linewidth=LW)
            if horalign == 'right':
                xpos = lineposition-0.2*windowwidth
            elif horalign == 'left':
                xpos = lineposition+0.2*windowwidth
            else:
                xpos = lineposition

            if ll in ['oiv1','oiv2','ovi1','ovi2']:
                ypos = yrange[1]*0.85
            else:
                ypos = yrange[1]*0.95

            plt.text(xpos,ypos,linename,color=col_linemarker,size=Fsize,
                     rotation='horizontal',horizontalalignment=horalign,verticalalignment='top')

            if voffset != 0.0:
                zoffset  = voffset*(redshift+1.0) / 299792.458
                range    = np.sort(np.asarray([((redshift+zoffset)+1)*linewave/wavescale, (redshift +1)*linewave/wavescale]))
                lineymin = yrange[0]
                lineymax = yrange[1]
                plt.fill_between(range,np.zeros(2)+lineymin,np.zeros(2)+lineymax,alpha=0.5,color=col_linemarker)
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =