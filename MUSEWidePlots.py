# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import os
import glob
import astropy.io.fits as afits
import kbsutilities as kbs
import numpy as np
import MUSEWidePlots as mwp
import MUSEWideUtilities as mwu
import tdose_extract_spectra as tes
import sys
import matplotlib.pyplot as plt
import collections
import MiGs
import aplpy
import tdose_utilities as tu
import felis
import scipy
import commands
import pdb
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def loadcatalogs(verbose=True):
    """
    Loading catalogs and returning combined MUSE-Wide catalog IDs, ra, dec, redshifts etc.

    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp
    IDlist, ra, dec, redshifts = mwp.loadcatalogs()

    """
    if verbose: print(' - Loading catalogs and defining ID, RA, Dec and redshift lists')
    photcat_e24 = '/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits'
    e24_id      = afits.open(photcat_e24)[1].data['UNIQUE_ID']
    e24_ra      = afits.open(photcat_e24)[1].data['RA']
    e24_dec     = afits.open(photcat_e24)[1].data['DEC']
    e24_z       = afits.open(photcat_e24)[1].data['REDSHIFT']

    #photcat_e36 = '/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_e36_v0.2.fits'
    photcat_e36 = '/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_e36_v1.0.fits'
    e36_id      = afits.open(photcat_e36)[1].data['ID']
    e36_ra      = afits.open(photcat_e36)[1].data['RA']
    e36_dec     = afits.open(photcat_e36)[1].data['DEC']
    e36_z       = afits.open(photcat_e36)[1].data['REDSHIFT']

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
    #if verbose: print('\n - Makeing sub-selection of objects')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Generating image for non-median-subtracted 1D spectra')
    outimage    = specdir+'Spec1DCOmbinationImage_e24ANDe36_'+date+'.fits'
    imgarray    = mwp.gen_Spec1DCombinationImage(IDlist,redshifts,specdir,specstring=specstring,NrowPerSpec=NrowPerSpec,
                                                 wavelengthgrid=wavegrid,fluxext=fluxext,waveext=waveext,
                                                 imagename=outimage,clobber=clobber,verbose=verbose,
                                                 scale='None',subtract='None')
    returnname  = outimage
    returnarray = imgarray
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Generating image for median-subtracted 1D spectra')
    outimage    = specdir+'Spec1DCOmbinationImage_e24ANDe36_'+date+'_mediansubtracted.fits'
    imgarray    = mwp.gen_Spec1DCombinationImage(IDlist,redshifts,specdir,specstring=specstring,NrowPerSpec=NrowPerSpec,
                                                 wavelengthgrid=wavegrid,fluxext=fluxext,waveext=waveext,
                                                 imagename=outimage,clobber=clobber,verbose=verbose,
                                                 scale='None',subtract='median')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Generating image for max-scaled and median-subtracted 1D spectra')
    outimage    = specdir+'Spec1DCOmbinationImage_e24ANDe36_'+date+'_maxscale_mediansubtracted.fits'
    imgarray    = mwp.gen_Spec1DCombinationImage(IDlist,redshifts,specdir,specstring=specstring,NrowPerSpec=NrowPerSpec,
                                                 wavelengthgrid=wavegrid,fluxext=fluxext,waveext=waveext,
                                                 imagename=outimage,clobber=clobber,verbose=verbose,
                                                 scale='max',subtract='median')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Returning resulting image array from '+returnname)
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
    if verbose: print(' - Setting up output image dimensions')

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
    if verbose: print(' - Filling output image with fluxes interpolated to requested wavelength grid')
    for ss in xrange(Nspec):
        objid   = str(IDsort[ss])
        globstr = spec1Ddir+specstring.replace('IIIII',objid)
        spec    = glob.glob(globstr)

        if len(spec) == 0:
            if verbose: print(' - WARNING: No spectrum found for '+objid)
        else:
            if len(spec) > 1:
                if verbose: print(' - WARNING: More than two spectra found when globbing for:\n            '+globstr)
                if verbose: print(' -          Will plot the first spectrum found')

            flux       = afits.open(spec[0])[1].data[fluxext]
            wave       = afits.open(spec[0])[1].data[waveext]
            fluxinterp = kbs.interpn(wave,flux,wavegrid)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if scale != 'None':
                if type(scale) == str:
                    if scale == 'max':
                        fluxinterp = fluxinterp/np.max(fluxinterp)
                    else:
                        if verbose: print(' WARNING "'+scale+'" invalid value of scale keyword; not scaling performed')
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
                        if verbose: print(' WARNING "'+scale+'" invalid value of scale keyword; not scaling performed')
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
    if verbose: print(' - Storing image to fits file '+imagename)
    hdu = afits.PrimaryHDU(imgarray)
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
    if verbose: print(' - Generating MUSE-Wide region files. Saving them to :')
    outputfile_CDFS   = outputfile_base+'_cdfs.reg'
    outputfile_COSMOS = outputfile_base+'_cosmos.reg'
    if verbose: print('   '+outputfile_CDFS)
    if verbose: print('   '+outputfile_COSMOS)

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
    if verbose: print(' - Loading source catalogs')
    IDlist, ra, dec, redshifts = mwp.loadcatalogs(verbose=verbose)

    if verbose: print(' - Setting up and generating histogram of MUSE-Wide sources in \n   '+plotname)
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
        if verbose: print(infostr)

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

    if verbose: print('   Saving plot... ')
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
    if verbose: print(' - Loaking for spectra by globbing using the input path-strings')
    tdosespecs   = glob.glob(tdosespecssearch)
    caruanaspecs = glob.glob(caruanaspecssearch)
    urrutiaspecs = glob.glob(urrutiaspecssearch)

    Ntdose =  len(tdosespecs)
    if Ntdose == 0:
        sys.exit('No tdose spectra found in '+tdosespecssearch)
    else:
        if verbose: print(' - Will generate comparison plots for the '+str(Ntdose)+' TDOSE spectra found ')

    if len(caruanaspecs) == 0:
        if verbose: print(' - WARNING No Caruana spectra found in \n   '+caruanaspecssearch)
    if len(urrutiaspecs) == 0:
        if verbose: print(' - WARNING No Urrutia spectra found in \n   '+urrutiaspecssearch)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Putting together lists of IDs')
    if verbose: print('   Generate ID list for TDOSE spectra')
    ids_tdose = [float(tt.split('/')[-1].split('_')[-1].split('.')[0]) for tt in tdosespecs]
    ids_tdose = np.asarray(ids_tdose)

    if verbose: print('   Generate ID list for Caruana spectra')
    ids_car = [float(tt.split('/')[-1].split('_')[-1].split('.')[0]) for tt in caruanaspecs]
    ids_car = np.asarray(ids_car)

    if len(urrutiaspecs) > 0:
        if verbose: print('   Generate ID list for Urrutia spectra')
        ids_urr_MW = [float(tt.split('/')[-1].split('_')[-1].split('.')[0]) for tt in urrutiaspecs]
        ids_urr_MW = np.asarray(ids_urr_MW)

        idmatch    = afits.open(TDOSE2MWidmatch)[1].data
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
    if verbose: print(' - Will store output figures in \n   '+outdir)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Looping over objects and generating plots')

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
        dat_tdose    = afits.open(tspec)[1].data
        t_flux    = dat_tdose['flux']
        t_fluxerr = dat_tdose['fluxerror']
        t_wave    = dat_tdose['wave']
        t_s2n     = dat_tdose['s2n']

        cent   = np.where(ids_car == ids_tdose[ii])[0]
        if len(cent) == 1:
            c_flux    = afits.open(caruanaspecs[cent[0]])[0].data
            c_fluxerr = afits.open(caruanaspecs[cent[0]])[1].data
            c_wave    = t_wave

        else:
            c_flux    = None
            c_fluxerr = None
            c_wave    = None

        if len(urrutiaspecs) > 0:
            uent   = np.where(ids_urr == str(int(ids_tdose[ii])))[0]
            if (len(uent) == 1) & (uent != 0):
                dat_urr   = afits.open(urrutiaspecs[uent[0]])[1].data
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
        #if verbose: print(' - Saving plot to'+plotname)
        plt.savefig(plotname)
        plt.clf()
        plt.close('all')

    if verbose: print('\n ... done ')
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
            if verbose: print(' WARNING the "markwave" is outside provided plotting xrange ')
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
    if verbose: print(' - Loaking for spectra by globbing using the input path-strings')
    tdosespecs   = glob.glob(tdosespecssearch)
    caruanaspecs = glob.glob(caruanaspecssearch)

    Ntdose =  len(tdosespecs)
    if Ntdose == 0:
        sys.exit('No tdose spectra found in '+tdosespecssearch)
    else:
        if verbose: print(' - Will collect and plot stats for the '+str(Ntdose)+' TDOSE spectra found ')

    if len(caruanaspecs) == 0:
        sys.exit(' No Caruana spectra found in '+caruanaspecssearch)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Putting together lists of IDs')
    if verbose: print('   Generate ID list for TDOSE spectra')
    ids_tdose = [float(tt.split('/')[-1].split('_')[-1].split('.')[0]) for tt in tdosespecs]
    ids_tdose = np.asarray(ids_tdose)

    if verbose: print('   Generate ID list for Caruana spectra')
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
        dat_tdose    = afits.open(tspec)[1].data
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
            c_flux    = afits.open(caruanaspecs[cent[0]])[0].data
            c_fluxerr = afits.open(caruanaspecs[cent[0]])[1].data
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
    #if verbose: print(' - Saving plot to',plotname)
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
def plot_1DspecOverview(spectra, labels, wavecols, fluxcols, fluxerrcols, redshift, voffset=0, show_error=True,
                        outputfigure='default', plotSN=False, skyspectra=None, wavecols_sky=None, fluxcols_sky=None,
                        yrangefull=None, xrangefull=[4500,18000], speccols=None, linenames=None,
                        col_matrix=False, col_matrix_title='The Color Matrix', col_matrix_text='Text',
                        col_matrix_labels=['xparam','yparam'], col_matrix_ranges=[[0,3],[0,3]],
                        col_matrix_binranges=None,col_matrix_p1dat=None,col_matrix_p2dat=None,
                        showFELISresults=False,FELISvetting=None,
                        verbose=True):
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
    linenames     List of (9) line names to generate zoom-ins around. If None hardcoded setups will be used

    col_matrix    These set of commands enables illustrating the binning of componsite (stacked) spectra
                  showin in the overview. This will add a matrix highlighting the binning location.
                  See uves.stack_composite_plotNxNspecs() for use.

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
    skyspectra   = ['/Users/kschmidt/work/MUSE/skyspectra/SKY_SPECTRUM_candels-cdfs-06_av.fits',None,'/Users/kschmidt/work/MUSE/skytable.fits']
    wavecols_sky = ['lambda',None,'lam']
    fluxcols_sky = ['data',None,'flux']
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
        specfigure = outputfigure.replace('.pdf','_z'+str(redshift).replace('.','p')+\
                     '_voffset'+str(voffset).replace('.','p')+'.pdf')

    if plotSN:
        specfigure = specfigure.replace('.pdf','_SN.pdf')
    if verbose: print(' - 1D overview figure will be saved to:\n   '+specfigure)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    redshiftplot = redshift + (voffset*(redshift+1.0) / 299792.458)
    if verbose: print(' - Will plot emission line markers at redshift '+str("%.6f" % redshift)+\
                      ' (z~'+str("%.6f" % redshiftplot)+' incl. (Lya) velocity offset = '+str(voffset)+'km/s)')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Nspec = len(spectra)
    if verbose: print(' - Loading the '+str(Nspec)+' spectra provided for plotting')

    datadic = collections.OrderedDict()

    for ss, specname in enumerate(spectra):
        spec_dat      = afits.open(specname)[1].data
        spec_wave     = spec_dat[wavecols[ss]]
        spec_flux     = spec_dat[fluxcols[ss]]
        spec_ferr     = spec_dat[fluxerrcols[ss]]
        spec_S2N      = spec_flux/spec_ferr
        spec_filllow  = spec_flux-spec_ferr
        spec_fillhigh = spec_flux+spec_ferr
        spec_wavecov  = [np.min(spec_wave),np.max(spec_wave)]

        if skyspectra[ss] is not None:
            spec_dat      = afits.open(skyspectra[ss])[1].data
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
    if verbose: print(' - Defining quanteties to plot; checking wavelength coverage of lines given redshift')
    llistdic      = MiGs.linelistdic(listversion='full') # loading line list for plots

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Plotting figure ')
    figuresize_x = 13
    figuresize_y = 8
    fig          = plt.figure(figsize=(figuresize_x,figuresize_y))
    Fsize        = 12
    LW           = 2
    plt.rc('text', usetex=True)                         # enabling LaTex rendering of text
    plt.rc('font', family='serif',size=Fsize)           # setting text font
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    left   = 0.06   # the left side of the subplots of the figure
    right  = 0.98   # the right side of the subplots of the figure
    bottom = 0.07   # the bottom of the subplots of the figure
    top    = 0.98   # the top of the subplots of the figure
    wspace = 0.25   # the amount of width reserved for blank space between subplots
    hspace = 0.32   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    if speccols is None:
        speccols          = ['blue','green','red','magenta','cyan','orange','purple','yellow','skyblue']

    xlabel            = '$\lambda$ / [\AA]'
    ylabel            = '$f_\lambda / [10^{-20}$erg/s/cm$^2$/\\AA]'
    if plotSN:
        ylabel        = 'S/N'
    col_linemarker    = 'gray'
    wavescale         = 1.0
    nrows             = 4
    ncols             = 3

    subplotinfo  = mwp.plot_1DspecOverview_subplotinfo()
    xrangedic    = {}
    yrangedic    = {}

    if linenames is None:
        if redshift > 2.5:
            linenames    = ['Lyb','Lya','NV','CII','SiIVOIV','CIV','HeII','CIII','MgII']
        elif redshift == 0.0:
            linenames    = ['Lyb','Lya+NV','CII','SiIVOIV','CIV','HeII','CIII','MgII','OII']
        elif (redshift > 1.5) & (redshift <= 2.5):
            linenames    = ['CIV','HeII','CIII','MgII','OII','NeIII','Hg','Hb','OIII']
        else:
            linenames    = ['MgII','OII','NeIII','Hg','Hb','OIII','HeI','Ha','SII']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for ll, linename in enumerate(linenames):
        linewave      = subplotinfo[linename][0]
        windowcenter  = subplotinfo[linename][1]
        windowwidth   = subplotinfo[linename][2]
        linelatex     = subplotinfo[linename][3]
        checkwave     = linewave*(1+redshift)
        speccoverage  = mwp.plot_1DspecOverview_checkcoverage(spectra,datadic,checkwave)

        if windowcenter*(redshift+1.0) > 10000.:
            windowwidth  = 70

        plotindex = ll+1
        xrangedic[linename], yrangedic[linename] = \
            mwp.plot_1DspecOverview_subplot(nrows,ncols,plotindex,linename,datadic,spectra,skyspectra,wavecols_sky,
                                            fluxcols_sky,speccols,plotSN,voffset,llistdic,Fsize,col_linemarker,LW,
                                            xlabel,ylabel,show_error,linelatex,windowcenter,windowwidth,wavescale,redshift,
                                            showFELISresults=showFELISresults,FELISvetting=FELISvetting,specinrange=speccoverage)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if col_matrix:
        plt.subplot(4, 3, (10,11)) # Full specs
    else:
        plt.subplot(4, 3, (10,12)) # Full specs
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    yrange = mwp.plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,
                                               xrangefull,show_error=show_error,plotSN=plotSN,labels=labels)
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    if not plotSN:
        if yrangefull is None:
            yrangefull  = yrange
    else:
        yrangefull  = yrange
    Dyrangefull = yrangefull[1]-yrangefull[0]

    # --- "ZOOM BOXES" ---
    boxzorder = 150

    for ll, linename in enumerate(linenames):
        linewave      = subplotinfo[linename][0]
        checkwave     = linewave*(1+redshift)
        speccoverage  = mwp.plot_1DspecOverview_checkcoverage(spectra,datadic,checkwave)
        if speccoverage:
            mwp.plot_1DspecOverview_genbox(subplotinfo[linename][3],xrangedic[linename],yrangedic[linename],
                                           LW,boxzorder,Dyrangefull,col_linemarker,Fsize)

    #--------- LEGEND ---------
    anchorpos = (0.5, 1.2)
    if len(labels) > 4:
        ncols = 4
    else:
        ncols = len(labels)
    leg = plt.legend(fancybox=True,numpoints=1, loc='upper center',prop={'size':Fsize},ncol=ncols)#,
                     #bbox_to_anchor=anchorpos)  # add the legend
    leg.get_frame().set_alpha(0.7)
    #--------------------------


    plt.xlim(xrangefull)
    plt.ylim(yrangefull)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if col_matrix:
        plt.subplot(4, 3, 12) # Color matrix
        plt.title(col_matrix_title)

        NxNmatrix = False
        if NxNmatrix:
            cmatrix      = np.zeros([3,3])
            pixval       = 0
            for xx in np.arange(np.sqrt(len(spectra))):
                for yy in np.arange(np.sqrt(len(spectra)))[::-1]:
                    cmatrix[int(yy),int(xx)] = pixval
                    pixtext = col_matrix_text[pixval]
                    #pixtext = 'pix'+str(int(xx))+str(int(yy))+'\\\\(val='+str(pixval)+')'
                    plt.text(xx,yy,pixtext,zorder=30,color='gray',#speccols[int(pixval)],
                             size=Fsize*2,horizontalalignment='center',verticalalignment='center')
                    pixval    = pixval + 1

            cmap         = plt.cm.gray
            cmap_norm    = plt.Normalize(cmatrix.min(), cmatrix.max())
            cmatrix_cmap = cmap(cmap_norm(cmatrix))

            pixval = 0
            for xx in np.arange(np.sqrt(len(spectra))):
                for yy in np.arange(np.sqrt(len(spectra)))[::-1]:
                    cmatrix_cmap[int(yy),int(xx),:3] = speccols[pixval]
                    pixval = pixval + 1

            plt.imshow(cmatrix_cmap,zorder=10)
            plt.xticks([])
            plt.yticks([])
            plt.xlabel(str(col_matrix_ranges[0][0])+' $<$ '+col_matrix_labels[0]+' $<$ '+str(col_matrix_ranges[0][1]))
            plt.ylabel(str(col_matrix_ranges[1][0])+' $<$ '+col_matrix_labels[1]+' $<$ '+str(col_matrix_ranges[1][1]))
        else:
            xr_full       = [col_matrix_binranges[0][0][0],col_matrix_binranges[0][-1][1]]
            yr_full       = [col_matrix_binranges[1][0][0],col_matrix_binranges[1][-1][1]]

            if (col_matrix_p1dat is not None) & (col_matrix_p2dat is not None):
                Nloops = int(np.sqrt(len(spectra)))
                plt.plot(col_matrix_p1dat,col_matrix_p2dat,'ok',zorder=20,markersize=1.5)
            if (col_matrix_p1dat is not None) & (col_matrix_p2dat is None):
                Nloops = int(len(spectra))
                plt.hist(col_matrix_p1dat,color='black',bins=20,histtype="step",lw=LW,zorder=20,
                         alpha=0.6,fill=True,fc='black')

            boxval = 0
            for xx in np.arange(Nloops):
                xr = col_matrix_binranges[0][xx]
                for yy in np.arange(Nloops):
                    yr = col_matrix_binranges[1][yy]
                    if col_matrix_p2dat is not None:
                        boxcol  = speccols[boxval]
                        boxtext = col_matrix_text[boxval]
                    else:
                        boxcol  = speccols[xx]
                        boxtext = col_matrix_text[xx]

                    plt.fill_between(xr,[yr[0],yr[0]],[yr[1],yr[1]],alpha=1.00,color=boxcol,zorder=10)

                    if col_matrix_p2dat is not None:
                        text_yoffset = np.abs(np.diff(yr)/2.)
                    else:
                        text_yoffset = np.abs(np.diff(yr)*3./4.)
                    plt.text(xr[0]+np.abs(np.diff(xr)/2.),yr[0]+text_yoffset,
                             boxtext,zorder=30,color='black',#backgroundcolor = 'gray',#speccols[int(pixval)],
                             size=Fsize*1.5,horizontalalignment='center',verticalalignment='center',
                             bbox=dict(boxstyle="round",edgecolor='none',facecolor='white',alpha=0.6))

                    boxval  = boxval + 1


            plt.xlim(xr_full)
            plt.ylim(yr_full)
            plt.xlabel(col_matrix_labels[0])
            plt.ylabel(col_matrix_labels[1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.savefig(specfigure, dpi=300) # dpi = dot per inch for rasterized points
    plt.clf()
    plt.close('all')
    if verbose: print(' - Saved figure to ',specfigure)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview_checkcoverage(spectra,datadic,checkwave):
    """
    Function checking if plot should be generated, i.e., if there is spectral coverage of line
    """
    coverage  = False

    for specname in spectra:
        wavecov = datadic[specname]['spec_wavecov']
        if (checkwave > wavecov[0]) & (checkwave < wavecov[1]):
            coverage = True

    return coverage


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview_subplotinfo():
    """
    Retruning dictionary with subplot in for generating subplots of different lines
    """
    subplotinfodic = {}

    #                             linewave    wavecenter,    wavewidth,    latexname
    subplotinfodic['Lyg']     =  [ 973.0,      973.0,        20.0,         'Ly$\gamma$']
    subplotinfodic['Lyb+OVI'] =  [1026.0,     1030.0,        20.0,         'Ly$\\beta$ + OVI']
    subplotinfodic['Lya+NV']  =  [1216.0,     1229.0,        30.0,         'Ly$\\alpha$ + NV']
    subplotinfodic['Lya']     =  [1216.0,     1216.0,        10.0,         'Ly$\\alpha$']
    subplotinfodic['NV']      =  [1216.0,     1241.0,        10.0,         'NV']
    subplotinfodic['CII']     =  [1335.0,     1335.0,        10.0,         'CII$\\lambda$1336']
    subplotinfodic['SiIVOIV'] =  [1400.0,     1400.0,        13.0,         'SiIV + OIV]']
    subplotinfodic['CIV']     =  [1549.0,     1549.0,        10.0,         'CIV']
    subplotinfodic['HeII']    =  [1640.0,     1640.0,        10.0,         'HeII']
    subplotinfodic['OIII1663']=  [1663.0,     1663.0,        10.0,         'OIII]']
    subplotinfodic['SiIII']   =  [1888.0,     1886.0,        10.0,         'SiIII']
    subplotinfodic['CIII']    =  [1909.0,     1908.0,        10.0,         'CIII]']
    subplotinfodic['CIIb']    =  [2326.0,     2326.0,        26.0,         'CIIb']
    subplotinfodic['MgII']    =  [2795.0,     2798.0,        10.0,         'MgII']
    subplotinfodic['OII']     =  [3727.0,     3727.5,        10.0,         '[OII]']
    subplotinfodic['NeIII']   =  [3869.0,     3869.0,        30.0,         '[NeIII]$\\lambda$3869']
    subplotinfodic['Hd']      =  [4101.0,     4101.0,        10.0,         'H\delta']
    subplotinfodic['Hg']      =  [4340.0,     4340.0,        30.0,         'H$\gamma$']
    subplotinfodic['Hb']      =  [4862.0,     4862.0,        30.0,         'H$\\beta$']
    subplotinfodic['OIII']    =  [5007.0,     4985.0,        60.0,         '[OIII]']
    subplotinfodic['HeI']     =  [5877.0,     5877.0,        30.0,         'HeI$\\lambda$5877']
    subplotinfodic['Ha']      =  [6563.0,     6563.0,        30.0,         'H$\\alpha$']
    subplotinfodic['SII']     =  [6717.0,     6723.0,        30.0,         '[SII]']

    return subplotinfodic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview_subplot(nrows,ncols,plotindex,
                                line,datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,plotSN,voffset,
                                llistdic,Fsize,col_linemarker,LW,xlabel,ylabel,show_error,
                                linestring,windowcenter,windowwidth,wavescale,redshift,
                                showFELISresults,FELISvetting,specinrange=True):
    """
    Function generating sub plots for plot_1DspecOverview()
    """
    plt.subplot(nrows,ncols,plotindex) # CII
    if specinrange:
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        yrange = mwp.plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,speccols,
                                                   xrange,show_error=show_error,plotSN=plotSN,labels=None)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        mwp.plot_1DspecOverview_plotlines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,
                                          xrange,yrange,redshift,LW,wavetype='vac')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        if not plotSN:
            objid             = int(spectra[0].split('_')[-1].split('.fit')[0])
            S2Nmin            = 3.0
            vshiftmax         = 1000.0

            if showFELISresults:
                line_felis = line.split('166')[0]
                FELISsummaryfile  = glob.glob(showFELISresults.replace('LLLL',line_felis))

                if len(FELISsummaryfile) == 1:
                    pickledir         = showFELISresults.split('CCresults')[0]
                    outkeystr         = FELISsummaryfile[0].split('template'+line_felis)[-1].split('.tx')[0]

                    fmt              = '12a,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,300a,300a'
                    summarydat       = np.genfromtxt(FELISsummaryfile[0],skip_header=25,dtype=fmt,comments='#',names=True)
                    objsummary_ent   = np.where(summarydat['id'].astype(int) == objid)[0]
                    zobj             = summarydat['z_spec'][objsummary_ent]

                    if len(objsummary_ent) != 0:
                        FELISpicklefile = glob.glob(pickledir+'*'+str(objid)+'*template'+line_felis+'*'+outkeystr+'*.pkl')

                        maxS2N = summarydat['FELIS_S2Nmax'][objsummary_ent]
                        vshift = summarydat['vshift_CCmatch'][objsummary_ent]
                        if (maxS2N > S2Nmin) & (np.abs(vshift) < vshiftmax):

                            if FELISvetting:
                                dat_vetfelis  = np.genfromtxt(FELISvetting,  dtype=None,comments='#',names=True,skip_header=12)
                                Nobj_fvet     = len(np.unique(dat_vetfelis['id']))

                                vetobjent = np.where(dat_vetfelis['id'] == objid)[0]
                                vetresult = dat_vetfelis['trust'+line_felis][vetobjent][0]
                            else:
                                vetresult = None

                            # if 'OIII' in line: pdb.set_trace()
                            mwp.plot_1DspecOverview_addFELIStemplate(FELISpicklefile[0],spectra[0],xrange,
                                                                     FELISvettingResult=vetresult,FELIScolor='red')

        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No coverage of\n'+linestring,
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        xrange = [0,1]
        yrange = [0,1]

    plt.xlim(xrange)
    plt.ylim(yrange)
    return list(xrange), list(yrange)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview_genbox(boxstring,xrange,yrange,LW,boxzorder,Dyrangefull,col_linemarker,Fsize):
    """
    Function generating boxes for full-spec view in plot_1DspecOverview()
    """
    plt.plot(xrange,np.zeros(2)+yrange[0],'-',color='black',lw=LW,zorder=boxzorder)
    plt.plot(xrange,np.zeros(2)+yrange[1],'-',color='black',lw=LW,zorder=boxzorder)
    plt.plot(np.zeros(2)+xrange[0],yrange,'-',color='black',lw=LW,zorder=boxzorder)
    plt.plot(np.zeros(2)+xrange[1],yrange,'-',color='black',lw=LW,zorder=boxzorder)
    plt.text(np.mean(np.asarray(xrange)),yrange[1]+0.03*Dyrangefull,boxstring,
             color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom',zorder=boxzorder)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview_plotspecs(datadic,spectra,skyspectra,wavecols_sky,fluxcols_sky,colors,xrange,
                                  show_error=True,plotSN=False,labels=None):
    """

    """
    if labels is None:
        labels = [None]*len(spectra)
    yrangecomb = [0,10]

    for ss, specname in enumerate(spectra):
        waveent = (datadic[specname]['spec_wave'] > xrange[0]) & (datadic[specname]['spec_wave'] < xrange[1])

        if plotSN:
            plt.step(datadic[specname]['spec_wave'], datadic[specname]['spec_S2N'], '-',
                     alpha=1.0,color=colors[ss],label=labels[ss],zorder=100,where='mid')

            plt.plot(datadic[specname]['spec_wave'], datadic[specname]['spec_S2N']*0.0+3.0, '--',
                     alpha=1.0,color='gray',zorder=100,lw=1)

            try:
                specvals = datadic[specname]['spec_S2N'][waveent][np.isfinite(datadic[specname]['spec_S2N'][waveent])]
                fluxmin  = np.min(np.asarray([0,  np.min(specvals) ]))
                fluxmax  = np.max(np.asarray([10, np.max(specvals) ]))
                yrange   = [fluxmin,1.1*fluxmax]
            except:
                yrange = [0,10]

        else:
            plt.step(datadic[specname]['spec_wave'], datadic[specname]['spec_flux'], '-',
                     alpha=1.0,color=colors[ss],label=labels[ss],zorder=100,where='mid')

            if show_error:
                plt.fill_between(datadic[specname]['spec_wave'],datadic[specname]['spec_filllow'],datadic[specname]['spec_fillhigh'],
                                 alpha=0.20,color=colors[ss],zorder=100,step='mid')

            try:
                specvals = datadic[specname]['spec_flux'][waveent][np.isfinite(datadic[specname]['spec_flux'][waveent])]
                fluxmin  = np.min(np.asarray([0,  np.min(specvals) ]))
                fluxmax  = np.max(np.asarray([10, np.max(specvals) ]))
                yrange   = [fluxmin,1.1*fluxmax]
            except:
                yrange = [0,100]

        if yrange[0] < yrangecomb[0]:
            yrangecomb[0] = yrange[0]

        if yrange[1] > yrangecomb[1]:
            yrangecomb[1] = yrange[1]

    for ss, specname in enumerate(spectra):
        if plotSN:
            waveent = (datadic[specname]['spec_wave'] > xrange[0]) & (datadic[specname]['spec_wave'] < xrange[1])
            meanerr = np.mean(datadic[specname]['spec_ferr'][np.isfinite(datadic[specname]['spec_ferr'])])
        else:
            meanerr = 1.0

        if skyspectra[ss] is not None:
            sky_w   = datadic[specname]['spec_wave_sky']
            sky_f   = datadic[specname]['spec_flux_sky']

            if skyspectra[ss] == '/Users/kschmidt/work/MUSE/skytable.fits':
                sky_w = sky_w * 1e4
                if plotSN:
                    sky_f = sky_f / 1e5
                else:
                    sky_f = sky_f / 1e3

            if '/Users/kschmidt/work/MUSE/spectra_noise/median_eff_noise_spectrum' in skyspectra[ss]:
                sky_f   = 2*sky_f

            sky_f   = sky_f / meanerr # scale to fit in S/N windows

            skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]) &
                               (sky_w > np.min(datadic[specname]['spec_wave'])) &
                               (sky_w < np.max(datadic[specname]['spec_wave'])))[0]
            skywave = sky_w[skyent]
            skylow  = np.zeros(len(skywave))
            skyflux = sky_f[skyent]
            skyhigh = skyflux
            #skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])

            plt.fill_between(skywave,skylow+yrangecomb[0],skyhigh+yrangecomb[0],alpha=1.0,color='black',zorder=1)

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
            plt.plot(np.zeros(2)+lineposition,yrange,color=col_linemarker,alpha=0.7,linestyle='-',linewidth=LW,zorder=20)
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
                     rotation='horizontal',horizontalalignment=horalign,verticalalignment='top',zorder=20)

            if np.abs(voffset) != 0.0:
                zoffset  = voffset*(redshift+1.0) / 299792.458
                range    = np.sort(np.asarray([((redshift+zoffset)+1)*linewave/wavescale, (redshift +1)*linewave/wavescale]))
                lineymin = yrange[0]
                lineymax = yrange[1]
                try:
                    plt.fill_between(range,np.zeros(2)+lineymin,np.zeros(2)+lineymax,alpha=0.5,color=col_linemarker,zorder=10)
                except:
                    pdb.set_trace()
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview_addFELIStemplate(FELISpicklefile,specname,xrange,FELISvettingResult=None,FELIScolor='red'):
    """
    Overplotting the best-fit template from FELIS

    """
    loaddic  = felis.load_picklefile(FELISpicklefile)

    spec_wave, spec_flux, spec_df, spec_s2n = felis.load_spectrum(specname,verbose=False)
    waveent = (spec_wave > xrange[0]) & (spec_wave < xrange[1])
    spec_wave, spec_flux, spec_df, spec_s2n = spec_wave[waveent], spec_flux[waveent], spec_df[waveent], spec_s2n[waveent]

    spec_dic  = loaddic[specname]
    z_spec    = spec_dic['zspec']

    # Getting the entry of the template with maximum S/N in the CC results
    max_S2N = np.max(spec_dic['S2NCCmaxvec'])
    besttemplate_ent = np.where(spec_dic['ccresultsarr_S2N'] == max_S2N)[0][0] # Template of max S/N

    template         = spec_dic['templatevec'][besttemplate_ent]
    max_z            = spec_dic['zCCmaxvec'][besttemplate_ent]

    t_wave_init, t_flux_init, t_df_init, t_s2n_init = felis.load_spectrum(template,verbose=False)
    func       = scipy.interpolate.interp1d(t_wave_init,t_flux_init,kind='linear',fill_value="extrapolate")
    t_flux     = func(spec_wave/(1.0+z_spec))

    # Getting the entry in the CC flux scalings vector for the given template where S/N is max
    max_S2N_ent      = np.where(spec_dic['ccresultsarr_S2N'][besttemplate_ent,:] == max_S2N)[0][0]

    t_wave_obs = t_wave_init*(1.0+z_spec)
    Npix = len(spec_flux)
    template_triplelength            = np.zeros(3*Npix)
    template_triplelength[0:Npix]    = t_flux
    template_shift_S2Nmax            = np.roll(template_triplelength, int(max_S2N_ent+np.floor(Npix/2.)))[Npix:-Npix]
    template_shift_S2Nmax_normalized = template_shift_S2Nmax/np.trapz(template_shift_S2Nmax,spec_wave/(1.0+z_spec))

    flux_scale_S2Nmax     = spec_dic['ccresultsarray_flux'][besttemplate_ent,max_S2N_ent]
    max_wave              = spec_dic['wavelengths'][max_S2N_ent]

    # moving template to observed frame and scaling it
    template_obs_flux = template_shift_S2Nmax_normalized * (flux_scale_S2Nmax / (1.0+z_spec))

    #scale by continuum offset
    template_obs_flux = template_obs_flux + spec_dic['continuumlevel']

    skipfelistemplate = False
    if FELISvettingResult is None:
        linestyle  = '-'
        FELIScolor = 'orange'
    elif FELISvettingResult == 1:
        linestyle = '-'
    elif FELISvettingResult == 9:
        linestyle = '--'
    elif FELISvettingResult == 0:
        #skipfelistemplate = True
        linestyle = '.'

    if not skipfelistemplate:
        plt.plot(spec_wave, template_obs_flux, linestyle, alpha=0.8 , color=FELIScolor,label=None, lw=5)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_FoVoverview(ras,decs,names,sizeFoV,outputdir='./',pointings=None,showregions=True,
                     pmin=1,pmax=99,vmin=-0.003,vmax=0.003,stretch='linear',invert=False,
                     cleanimage=False,fontsize=30,clobber=False,verbose=True):
    """
    Plot overview of MUSE-Wide object (ra dec) field-of-view


    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp
    import afits

    objdata   = afits.open('/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits')[1].data
    names     = [str(id) for id in objdata['id']][-5:]
    ras       = objdata['ra'][-5:]
    decs      = objdata['dec'][-5:]
    pointings = objdata['pointing'][-5:]
    sizeFoV   = [5,5]
    outputdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVoverviews/'

    mwp.plot_FoVoverview(ras,decs,names,sizeFoV,pointings=pointings,outputdir=outputdir,cleanimage=False,clobber=False,showregions=True)

    """
    outnames = ['FoVoverview_'+str(name)+'.png' for name in names]

    if (len(ras) != len(decs)) or (len(ras) != len(outnames)):
        sys.exit('The number of ras ('+str(len(ras))+') and decs ('+str(len(ras))+') provided do not match')

    Nplots = len(ras)
    if verbose: print(' - Lopping over '+str(Nplots)+' coordinate pares to generate plots for ')
    for rr,ra in enumerate(ras):
        plotname = outputdir+outnames[rr]
        infostr  = ' - Generating '+plotname+' ('+str(rr+1)+'/'+str(Nplots)+') '
        if os.path.isfile(plotname) & (clobber == False):
            infostr = infostr+' ... plot exists --> skipobj'
        else:
            infostr = infostr+'                           '
        if verbose: print('\n'+infostr)
        if os.path.isfile(plotname) & (clobber == False): continue

        if ra > 100:  # <------------------ COSMOS
            regfile = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVoverviews/MUSE-Wide_objects_cosmos.reg'
            #regfile = '/Users/kschmidt/work/catalogs/MUSE_GTO/MUSE-Wide_objects_cosmos.reg'
            imgfile = '/Users/kschmidt/work/images_MAST/hlsp_candels_hst_wfc3_cos-tot_f125w_v1.0_drz.fits'
        else:         # <------------------ CDFS
            regfile = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/FoVoverviews/MUSE-Wide_objects_cdfs.reg'
            #regfile = '/Users/kschmidt/work/catalogs/MUSE_GTO/MUSE-Wide_objects_cdfs.reg'
            imgfile = '/Users/kschmidt/work/images_MAST/hlsp_candels_hst_wfc3_gs-tot_f125w_v1.0_drz.fits'

        if pointings is not None: # change image file if poitnings are provided
            imgfile = '/Users/kschmidt/work/images_MAST/MUSEWidePointings/acs_814w_'+pointings[rr]+'_cut_v1.0.fits'

        xcen        = ras[rr]
        ycen        = decs[rr]
        xwindow     = sizeFoV[0]/3600. # arcsec -> degrees
        ywindow     = sizeFoV[1]/3600.
        regionfiles = [regfile]
        cutname     = './FoVoverview_tempcutout.fits'
        cutout      = tu.extract_subimage(imgfile,xcen,ycen,[xwindow*3600.,ywindow*3600.],outname=cutname,
                                          clobber=clobber,verbose=verbose)

        # ------------------------ PLOTTING FOV ------------------------
        fsize = fontsize
        plt.rc('text', usetex=True)                         # enabling LaTex rendering of text
        plt.rc('font', family='serif',size=fsize)           # enable serif rendering of text
        fig = aplpy.FITSFigure(cutname)
        #fig.show_rgb(rgbclean)
        fig.show_grayscale(vmin=vmin,vmax=vmax,pmin=pmin, pmax=pmax,stretch=stretch, invert=invert)

        # - - - - ZOOM - - - -
        #fig.recenter(xcen, ycen, width=xwindow, height=ywindow)  # degrees

        # - - - - TICKS - - - -
        fig.tick_labels.set_font(size=fsize,family='serif')
        # fig.ticks.set_xspacing(1.0/3600.)  # degrees
        # fig.ticks.set_yspacing(1.0/3600.)  # degrees
        fig.tick_labels.set_xformat('hh:mm:ss.s')
        fig.tick_labels.set_yformat('dd:mm:ss.s')
        fig.ticks.set_linewidth(4)

        fig.axis_labels.set_font(size=fsize,family='serif')

        # - - - - GRID - - - -
        fig.add_grid()
        fig.grid.set_color('white')
        fig.grid.set_alpha(0.3)

        if cleanimage:
            fig.ticks.hide()
            fig.axis_labels.hide()
            fig.tick_labels.hide()
            fig.grid.hide()

        # - - - - COLOR BAR - - - -
        fig.add_colorbar()
        fig.colorbar.set_width(0.3)
        fig.colorbar.set_location('right')

        # - - - - CONTOURS - - - -
        # if contourimg:
        #     fig.show_contour(contourimg, levels=Clevel, colors=Ccolor,alpha=0.5,linewidths=3)

        # - - - - DS9 REGION - - - -
        if showregions:
            for region in regionfiles:
                if os.path.isfile(region):
                    fig.show_regions(region)
                else:
                    if verbose: print('WARNING: The region file '+region+' does not exist, hence nothing shown')

        # - - - - LABELS - - - -
        # fig.add_label(34.455, 54.112, 'My favorite galaxy')
        fig.add_label(0.5,1.05,'Object: '+names[rr].replace('_','\_'),color='black',relative=True,family='serif',size=fsize)

        # - - - - SAVE - - - -
        fig.save(plotname,dpi=300)
        #fig.save(plotname)
        plt.clf()
        plt.close('all')

        os.remove(cutname)



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_MgIIemitterUDF884spec(datestamp='190515',smoothsigma=0,verbose=True):
    """
    Function plotting the spectra for MgII emitter 884 in UDF03

    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp
    mwp.plot_MgIIemitterUDF884spec(smoothsigma=0)

    """
    figuredir       = '/Users/kschmidt/work/MUSE/MgIIemittersUDF3/tdose_spectra/'+datestamp+'/'

    specdir         = figuredir
    spec_modelimg   = specdir+'tdose_spectrum_modelimg_0000000884-0000000884.fits'
    spec_gauss      = specdir+'tdose_spectrum_gauss_0000000884-0000000884.fits'

    objz            = 0.737
    #-------------------------------------------------------------------------------------------------------
    if verbose: print(' - Defining emission line lists ')
    linelistdic  = MiGs.linelistdic()

    for kk, key in enumerate(linelistdic.keys()):
        if kk == 0:
            linelist_all = np.array([linelistdic[key][1]*(1.0+objz),linelistdic[key][0]])
        else:
            linelist_all = np.vstack((linelist_all,[linelistdic[key][1]*(1.0+objz),linelistdic[key][0]]))

    #-------------------------------------------------------------------------------------------------------
    if verbose: print(' - Setting up plot ranges and line lists ')
    specs        = [spec_modelimg,spec_gauss]

    filelist     = [specs[0]]
    col          = ['black']
    labels       = ['Multi-component Sersic model']

    compspec     = [specs[1]]
    comp_colors  = ['red']
    comp_labels  = ['Single-component Gauss model']

    xranges      = [[4800,9300],[4800,9300]]
    ylogval      = False

    if datestamp == '190510': # MUSE-Wide version
        yranges_full = [[-100,2500],[-1,50]]
        yranges_zoom = [[100,350],[-1,10]]
    elif datestamp == '190515':
        yranges_full = [[-100,2500],[-1,250]]
        yranges_zoom = [[100,400],[-1,75]]
    else:
        sys.exit('Invalid datestamp provided')

    plotnames    = [figuredir+'/tdose_1Dspectra_singleVSmulticomp_UDF884_full_flux.pdf',
                    figuredir+'/tdose_1Dspectra_singleVSmulticomp_UDF884_zoom_flux.pdf']

    linesetup = {}
    linesetup[plotnames[0]] = [linelist_all], ['black']
    linesetup[plotnames[1]] = [linelist_all], ['black']

    #-------------------------------------------------------------------------------------------------------
    if verbose: print(' - Plotting spectra ')
    showfluxnoise = True
    for pp, pname in enumerate(plotnames):
        linelist, linecols = linesetup[pname]

        # - - - - - - - - - - FLUX AND S/N PLOTS - - - - - - - - - -
        plotname  = pname
        xrange    = xranges[pp]
        if pp == 0:
            yrange    = yranges_full[0]
        else:
            yrange    = yranges_zoom[0]
        tes.plot_1Dspecs(filelist,plotname=plotname,colors=col,labels=labels,plotSNcurve=False,
                         comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                         comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',
                         xrange=xrange,yrange=yrange,showspecs=False,shownoise=showfluxnoise,verbose=True,pubversion=True,
                         showlinelists=linelist,linelistcolors=linecols,smooth=smoothsigma,ylog=ylogval)

        plotname  = pname.replace('flux.pdf','s2n.pdf')
        if pp == 0:
            yrange    = yranges_full[1]
        else:
            yrange    = yranges_zoom[1]
        tes.plot_1Dspecs(filelist,plotname=plotname,colors=col,labels=labels,plotSNcurve=True,
                         comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                         comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',
                         xrange=xrange,yrange=yrange,showspecs=False,shownoise=showfluxnoise,verbose=True,pubversion=True,
                         showlinelists=linelist,linelistcolors=linecols,smooth=smoothsigma,ylog=ylogval)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_MgIIemitterUDF939spec(datestamp='190515',smoothsigma=0,verbose=True):
    """
    Function plotting the spectra for MgII emitter 884 in UDF03

    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp
    mwp.plot_MgIIemitterUDF939spec(datestamp='190515',smoothsigma=0)

    """
    figuredir       = '/Users/kschmidt/work/MUSE/MgIIemittersUDF3/tdose_spectra/'+datestamp+'/'

    specdir         = figuredir
    spec_modelimg   = specdir+'tdose_spectrum_modelimg_0000000939-0000000939.fits'
    spec_gauss      = specdir+'tdose_spectrum_gauss_0000000939-0000000939.fits'

    objz            = 1.296
    #-------------------------------------------------------------------------------------------------------
    if verbose: print(' - Defining emission line lists ')
    linelistdic  = MiGs.linelistdic()

    for kk, key in enumerate(linelistdic.keys()):
        if kk == 0:
            linelist_all = np.array([linelistdic[key][1]*(1.0+objz),linelistdic[key][0]])
        else:
            linelist_all = np.vstack((linelist_all,[linelistdic[key][1]*(1.0+objz),linelistdic[key][0]]))

    #-------------------------------------------------------------------------------------------------------
    if verbose: print(' - Setting up plot ranges and line lists ')
    specs        = [spec_modelimg,spec_gauss]

    filelist     = [specs[0]]
    col          = ['black']
    labels       = ['Multi-component Sersic model']

    compspec     = [specs[1]]
    comp_colors  = ['red']
    comp_labels  = ['Single-component Gauss model']

    xranges      = [[4800,9300],[4800,9300]]
    ylogval      = False

    if datestamp == '190510': # MUSE-Wide version
        yranges_full = [[-50,800],[-1,15]]
        yranges_zoom = [[-10,150],[-1,3]]
    elif datestamp == '190515':
        yranges_full = [[-50,1000],[-1,80]]
        yranges_zoom = [[-10,200],[-1,10]]
    else:
        sys.exit('Invalid datestamp provided; valid values are:')

    plotnames    = [figuredir+'/tdose_1Dspectra_singleVSmulticomp_UDF939_full_flux.pdf',
                    figuredir+'/tdose_1Dspectra_singleVSmulticomp_UDF939_zoom_flux.pdf']

    linesetup = {}
    linesetup[plotnames[0]] = [linelist_all], ['black']
    linesetup[plotnames[1]] = [linelist_all], ['black']

    #-------------------------------------------------------------------------------------------------------
    if verbose: print(' - Plotting spectra ')
    showfluxnoise = True
    for pp, pname in enumerate(plotnames):
        linelist, linecols = linesetup[pname]

        # - - - - - - - - - - FLUX AND S/N PLOTS - - - - - - - - - -
        plotname  = pname
        xrange    = xranges[pp]
        if pp == 0:
            yrange    = yranges_full[0]
        else:
            yrange    = yranges_zoom[0]
        tes.plot_1Dspecs(filelist,plotname=plotname,colors=col,labels=labels,plotSNcurve=False,
                         comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                         comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',
                         xrange=xrange,yrange=yrange,showspecs=False,shownoise=showfluxnoise,verbose=True,pubversion=True,
                         showlinelists=linelist,linelistcolors=linecols,smooth=smoothsigma,ylog=ylogval)

        plotname  = pname.replace('flux.pdf','s2n.pdf')
        if pp == 0:
            yrange    = yranges_full[1]
        else:
            yrange    = yranges_zoom[1]
        tes.plot_1Dspecs(filelist,plotname=plotname,colors=col,labels=labels,plotSNcurve=True,
                         comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                         comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',
                         xrange=xrange,yrange=yrange,showspecs=False,shownoise=showfluxnoise,verbose=True,pubversion=True,
                         showlinelists=linelist,linelistcolors=linecols,smooth=smoothsigma,ylog=ylogval)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def subcatsMWfootprint_diagnostics(catname='Skelton',plotdir='/Users/kschmidt/work/MUSE/MWv2_analysis/continuum_source_selection/',
                                   skeltonwhitakermag='814',xrange=None,bins=None,verbose=True):
    """
    Diagnostics of the sub-catalogs in the MUSE-Wide footprint generated with mwu.create_subcatalogs_inMWfootprint()

    --- EXAMPLE OF USE ---
    import MUSEWidePlots as mwp
    mwp.subcatsMWfootprint_diagnostics(catname='Skelton',skeltonwhitakermag='814')

    """

    ids      = np.array([])
    mags     = np.array([])
    magnames = np.array([])

    if (catname.lower() == 'skelton') or (catname.lower() == 'skelton_goodss') or (catname.lower() == 'all'):
        photcat_goodss = '/Users/kschmidt/work/catalogs/MUSE_GTO/goodss_3dhst.v4.1_inMUSEWideFootprint.fits'
        photdat_goodss = afits.open(photcat_goodss)[1].data
        ids            = np.append(ids,photdat_goodss['id']+1100000000)
        if skeltonwhitakermag in ['775','606']:
            magcol         = 'f_F'+skeltonwhitakermag+'W'
        else:
            magcol         = 'f_F'+skeltonwhitakermag+'Wcand'
        mags           = np.append(mags,25.0-2.5*np.log10(photdat_goodss[magcol]))
        magnames       = np.append(magnames,magcol.replace('_','\_'))

    if (catname.lower() == 'skelton') or (catname.lower() == 'skelton_cosmos') or (catname.lower() == 'all'):
        photcat_cosmos = '/Users/kschmidt/work/catalogs/MUSE_GTO/cosmos_3dhst.v4.1_inMUSEWideFootprint.fits'
        photdat_cosmos = afits.open(photcat_cosmos)[1].data
        ids            = np.append(ids,photdat_cosmos['id']+2100000000)
        magcol         = 'f_F'+skeltonwhitakermag+'W'
        mags           = np.append(mags,25.0-2.5*np.log10(photdat_cosmos[magcol]))
        magnames       = np.append(magnames,magcol.replace('_','\_'))

    if (catname.lower() == 'whitaker') or (catname.lower() == 'all'):
        photcat  = '/Users/kschmidt/work/catalogs/MUSE_GTO/hlsp_hlf_hst_60mas_goodss_v2.0_catalog_inMUSEWideFootprint.fits'
        photdat  = afits.open(photcat)[1].data
        ids      = np.append(ids,photdat['id']+1200000000)
        magcol   = 'f_f'+skeltonwhitakermag+'w'
        mags     = np.append(mags,25.0-2.5*np.log10(photdat[magcol]))
        magnames = np.append(magnames,magcol.replace('_','\_'))

    if (catname.lower() == 'laigle') or (catname.lower() == 'all'):
        photcat  = '/Users/kschmidt/work/catalogs/MUSE_GTO/cosmos2015_laigle_v1.1_candelsregion_inMUSEWideFootprint.fits'
        photdat  = afits.open(photcat)[1].data
        ids      = np.append(ids,photdat['NUMBER']+2200000000)
        magcol   = 'V_MAG_ISO'
        mags     = np.append(mags,photdat[magcol])
        magnames = np.append(magnames,magcol.replace('_','\_'))

    if len(ids) == 0:
        sys.exit('No IDs available for "catname='+str(catname)+'"')

    goodent    = np.where((mags < 40) & (mags > 5) & np.isfinite(mags))[0]
    mags_good  = mags[goodent]
    ids_good   = ids[goodent]

    Nbad       = len(ids) - len(ids_good)
    Ncosmos    = len(np.where(ids_good > 1.9e9)[0])
    Ngoodss    = len(np.where(ids_good < 1.9e9)[0])
    Ntotal     = Ngoodss+Ncosmos

    if verbose: print(' - Read the catalog selection "'+catname+'" finding the following number of sources:')
    if verbose: print('   (discarding '+str(Nbad)+' sources for not being finite or having poor mags)')
    if verbose: print('   Total   :   '+str(Ntotal))
    if verbose: print('   GOODS-S :   '+str(Ngoodss))
    if verbose: print('   COSMOS  :   '+str(Ncosmos))

    # - - - - - - - - - - - - - - - - - - - - PLOTTING - - - - - - - - - - - - - - - - - - - -
    if catname.lower() == 'all':
        magext = 'm'+skeltonwhitakermag
    else:
        magext = magcol
    plotname = plotdir+'mag_histogram_'+catname.lower()+'_'+magext+'.pdf'
    if verbose: print(' - Setting up and generating histogram of MUSE-Wide sources in \n   '+plotname)
    fig = plt.figure(figsize=(5, 4))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.2, right=0.95, bottom=0.2, top=0.95)
    Fsize    = 14
    lthick   = 1.5
    marksize = 3
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title('M^\star',fontsize=Fsize)

    if xrange is None:
        xrange = [np.min(mags_good),np.max(mags_good)]

    if (bins is None):
        bin_dz = 0.1
        bins   = np.arange(np.min(mags_good),np.max(mags_good)+bin_dz,bin_dz)
        if xrange is not None:
            bins   = np.arange(np.min(xrange),np.max(xrange)+bin_dz,bin_dz)


    magranges = [[0,24],[24,25],[25,26],[26,99],[0,26]]
    colors    = ['blue','green','orange','red','black']

    for mm, magrange in enumerate(magranges):
        goodent   = np.where((mags_good > magrange[0]) & (mags_good <= magrange[1]))[0]
        Ngood     = len(goodent)

        if Ngood>1:
            goodIDs   = ids[goodent]
            goodmag   = mags_good[goodent]
            goodcolor = colors[mm]
            magmin    = np.min(goodmag)
            magmax    = np.max(goodmag)

            infostr   = '   Histinfo:'

            percent  = float(Ngood)/float(Ntotal)*100.
            label    = str(magrange[0])+'$<$mag$<=$'+str(magrange[1])+' \n('+str(Ngood)+' obj; '+str('%.2f' % percent)+'\%)'

            if mm < len(magranges)-1:
                fillval = True
                linest  = '-'
            else:
                fillval = False
                linest  = ':'
            hist     = plt.hist(goodmag,color=goodcolor,bins=bins,histtype="step",lw=lthick,label=label,ls=linest,
                                fill=fillval,fc=goodcolor,alpha=0.8)

    plt.xlim(xrange)
    plt.xlabel('AB magnitude \n('+', '.join(list(magnames))+')', fontsize=Fsize)

    #plt.ylim(yrange)
    plt.ylabel(catname.replace('_','\_')+' catalog objects\nover MUSE-Wide 100 field footprint', fontsize=Fsize)

    #--------- LEGEND ---------
    anchorpos = (0.5, 1.2)
    leg = plt.legend(fancybox=True,numpoints=1, loc='upper left',prop={'size':Fsize-3},ncol=1)#,
                     #bbox_to_anchor=anchorpos)  # add the legend
    leg.get_frame().set_alpha(0.7)
    #--------------------------

    plt.savefig(plotname)
    plt.clf()
    plt.close('all')




# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =