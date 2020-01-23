# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                       Utilities for stacking data arrays (spectra or images)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import glob
import os
import datetime
import kbsutilities as kbs
import stacking
import numpy as np
import sys
import MiGs
import pdb
import scipy
from scipy.interpolate import interp1d
import astropy.io.fits as afits
import matplotlib.pyplot as plt
import collections

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def stack_1D(wavelengths, fluxes, variances, stacktype='mean', wavemin=4500, wavemax=9500, deltawave=10.0,
             z_systemic=0.0, Nsigmaclip=None, outfile=None, verbose=True):
    """
    Stacking of 1D data arrays.

    --- INPUT ---
    wavelengths     List of wavelength grids for spectra to stack
    fluxes          List of fluxes of spectra to stack
    variances       List of variance spectra for the fluxes to stack
    stacktype       Type of stacing to perform. Choices are:
                        'mean'      Simple mean stack of spectra. Uncertainty returned as var_pix = Sum(var_i)/N
                                    where N are the number of pixels stacked
                        'median'    Median valie of stacked pixels
    wavemin         Minimum wavelength of output grid to interpolate spectra to
    wavemax         Maximum wavelength of output grid to interpolate spectra to
    deltawave       Wavelength resolution of output grid to interpolate spectra to
    z_systemic      Systemic redshifts if spectra to stack are to be moved to rest-frame
    Nsigmaclip      To clip the pixels to stack at N sigma provide the Nsigma value here. Default is no clipping.
    outfile         To save stack to a fits file (using the same format as TDOSE and FELIS)
                    provide path and name of output file
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---
    import stacking

    wave_out, flux_out, variance_out, Nspecstack = stacking.stack_1D()

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Defining output wavelength grid of stack')
    wave_out = np.arange(wavemin,wavemax+deltawave,deltawave)

    Nspec    = len(fluxes)
    if (Nspec != len(fluxes)) or (Nspec != len(variances)):
        sys.exit('Mis-match in number of wavelength ('+str(Nspec)+'), flux ('+str(len(fluxes))+
                 ') and variance ('+str(len(variances))+') vectors provided')

    fluxarr = np.zeros([len(wave_out),Nspec])
    fluxarr[:] = np.nan
    vararr  = np.zeros([len(wave_out),Nspec])
    vararr[:]  = np.nan

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Filling data arrays with spectra interpolated to output wavelength grid')
    for ww, wave in enumerate(wavelengths):
        newx            = wave_out
        oldx            = wave / (1.0 + z_systemic[ww])

        # interpolate fluxes (moving them to rest-frame if needed):
        oldy            = fluxes[ww] * (1.0 + z_systemic[ww])
        newy            = interp1d(oldx, oldy, kind='linear', fill_value=np.nan, bounds_error=False)(newx)
        fluxarr[:,ww]   = newy

        # interpolate variances (moving them to rest-frame if needed):
        oldy            = variances[ww] * (1.0 + z_systemic[ww])
        newy            = interp1d(oldx, oldy, kind='linear', fill_value=np.nan, bounds_error=False)(newx)
        vararr[:,ww]    = newy

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if Nsigmaclip is not None:
        if verbose:
            print(' - Clipping spectra to stack at '+str(Nsigmaclip)+' sigma (std) at each wavelength around '+stacktype.lower())
            print('   (setting flux values to NaN so they are handled by array masking)')
        for ss, spec in enumerate(wave_out):
            colvals_sort   = np.sort(fluxarr[ss,np.isfinite(fluxarr[ss,:])])
            Nval           = len(colvals_sort)
            if Nval > 0: # only clipping if there are values different from 0
                percentile_half = scipy.stats.norm(0, 1).cdf(-1*Nsigmaclip)
                Nclip           = int(Nval*percentile_half)
                if Nclip > 0:
                    value_low   = colvals_sort[Nclip-1]
                    value_high  = colvals_sort[-Nclip]
                    badent      = np.where( (fluxarr[ss,:] <= value_low) | (fluxarr[ss,:] >= value_high))[0]
                    fluxarr[ss,badent] = np.nan

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Performing stacking using stacktype = "'+str(stacktype)+'"')
    if stacktype.lower() == 'mean':
        mask_NaN       = np.ma.masked_invalid(fluxarr).mask
        mask_HighErr   = (np.sqrt(vararr) < 0.0)
        combined_mask  = (mask_NaN | mask_HighErr)

        fluxarr_ma      = np.ma.masked_array(fluxarr,mask=combined_mask,fill_value = np.nan)
        flux_out_ma     = np.mean(fluxarr_ma,axis=1)
        flux_out        = flux_out_ma.filled(fill_value=np.nan)

        fluxarr_ones    = fluxarr * 0.0 + 1.0
        Nspecstack_ma   = np.ma.masked_array(fluxarr_ones,mask=combined_mask,fill_value = np.nan)
        Nspecstack_ma   = np.sum(Nspecstack_ma,axis=1)
        Nspecstack      = Nspecstack_ma.filled(fill_value=0.0)

        vararr_ma       = np.ma.masked_array(vararr,mask=combined_mask,fill_value = np.nan)
        variance_out_ma = np.sum(vararr_ma,axis=1)
        # standard error on the mean, i.e., std scales as 1/sqrt(N) so
        # var_stack = Sum(var_i)/N    =>    std = sqrt(Sum(std_i)**2/N) = Sum(std_i) / sqrt(N)
        variance_out    = variance_out_ma.filled(fill_value=np.nan) / Nspecstack
    elif stacktype.lower() == 'median':
        fluxarr_ma      = np.ma.masked_array(fluxarr,mask=~np.isfinite(fluxarr),fill_value = np.nan)
        flux_out_ma     = np.median(fluxarr_ma,axis=1)
        flux_out        = flux_out_ma.filled(fill_value = np.nan)

        fluxarr_ones    = fluxarr * 0.0 + 1.0
        Nspecstack_ma   = np.ma.masked_array(fluxarr_ones,mask=~np.isfinite(fluxarr_ones),fill_value = np.nan)
        Nspecstack_ma   = np.sum(Nspecstack_ma,axis=1)
        Nspecstack      = Nspecstack_ma.filled(fill_value = np.nan)

        vararr_ma       = np.ma.masked_array(vararr,mask=~np.isfinite(fluxarr),fill_value = np.nan)
        variance_out_ma = np.median(vararr_ma,axis=1)
        variance_out    = variance_out_ma.filled(fill_value = np.nan)

        # alternatively, loop over columns and take confidence intervals, to indicate scatter of errors
        # Or use the error from the median value itself.

    else:
        if verbose: print('WARNING - stack_1D() did not recognize stacktype = "'+str(stacktype)+'"; returning None')
        flux_out        = None
        variance_out    = None
        Nspecstack      = None
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if outfile is not None:
        if verbose: print(' - Saving stack to '+outfile)
        stacking.save_stack_1D(outfile,wave_out,flux_out,np.sqrt(variance_out),Nspecstack,
                               headerinfo=None,overwrite=True,verbose=verbose)

    if verbose: print(' - Returning wavelengths, fluxes, variances of stack and '
                      'number of spectra contributing to each stacked pixel')
    return wave_out, flux_out, variance_out, Nspecstack
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def stack_1D_wrapper(spectra, redshifts, outfile, verbose=True,
                     stacktype='median', wavemin=1100, wavemax=3000, deltawave=0.1):
    """
    Wrapper around stack_1D() to load spectra on the TDOSE/FELIS format (and corresponding redshifts) and stack them

    --- INPUT ---
    spectra         List of spectra to stack in rest-frame
    redshifts       Redshifts used to move spectra to rest-frame
    outfile         File to save stacked spectrum to
    verbose         Toggle verbosity

    stacktype, wavemin, wavemax, deltawave are passed directly to stac_1D()

    --- EXAMPLE OF USE ---
    import stacking
    import glob

    spectra = glob.glob('/Volumes/DATABCKUP1/TDOSEextractions/180824_TDOSEextraction_LAEs60fields/modelimg/tdose_spectra/tdose_spectrum_candels-cdfs-16*.fits')
    stacktype = 'median'
    infofile  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo_UVemitters_3timesUDFcats.fits'
    objinfo   = afits.open(infofile)[1].data
    z_Lya     = objinfo['redshift'][129:140]

    wave_out, flux_out, variance_out, Nspecstack = stacking.stack_1D_wrapper(spectra,z_Lya, '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_spectrum_candels-cdfs-16_LAEs_stack_'+stacktype+'.fits',stacktype=stacktype)

    """
    if verbose: print(' - Will stack the '+str(len(spectra))+' provided spectra')

    wavelengths  = []
    fluxes       = []
    variances    = []

    for spectrum in spectra:
        data        = afits.open(spectrum)[1].data
        wavelengths.append(data['wave'])
        fluxes.append(data['flux'])
        variances.append(data['fluxerror']**2.0)

    if verbose: print(' - Start stacking of spectra in rest-frame')
    wave_out, flux_out, variance_out, Nspecstack = \
        stacking.stack_1D(wavelengths, fluxes, variances, z_systemic=redshifts,
                          stacktype=stacktype, wavemin=wavemin, wavemax=wavemax,
                          deltawave=deltawave, outfile=outfile, verbose=verbose)

    return wave_out, flux_out, variance_out, Nspecstack
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview(spectra, labels, wavecols, fluxcols, fluxerrcols, outname,
                        zoomwindows=[[1530,1570,'CIV'],[1620,1660,'HeII'],[1900,1920,'CIII'],[2780,2830,'MgII']],
                        redshift=0.0, plotSN=False, skyspectra=None, wavecols_sky=None, fluxcols_sky=None,
                        yrangefull=None, xrangefull=[0,10000], speccols=None, verbose=True):
    """

    Plotting overview with zoom-ins of 1D spectrum.
    Based on MUSEWidePlots.plot_1DspecOverview*()

    Can be used to plot multiple spectra for comparison purposes.

    --- INPUT ---
    spectra       List of spectra to include in plot
    labels        Labels to use for spectra
    wavecols      Column names of entries in "spectra" files containing wavelengths     [in angstrom]
    fluxcols      Column names of entries in "spectra" files containing fluxes          [1e-20 cgs]
    fluxerrcols   Column names of entries in "spectra" files containing flux errors     [1e-20 cgs]
    outname       Path and name of figure to generate
    zoomwindows   Regions of spectrum to show in sub-plot windows
    redshift      Redshift to position emission line markers. Default is rest-frame.
    plotSN        If true a S/N version (flux/fluxerr instead of flux) of the figure will be generated
    yrangefull    To fix the y-range for the overview panel including all spectra set it with this keyword
    verbose       Toggle verbosity

    --- EXAMPLE OF USE ---
    import stacking
    specfile = './MUSEWideStack_test181029.fits'
    zoomwindows=[[1180,1280,'Lya'],[1530,1570,'CIV'],[1620,1660,'HeII'],[1900,1920,'CIII'],[2780,2830,'MgII']]
    stacking.plot_1DspecOverview([specfile],['teststack'],['wave'],['flux'],['fluxerror'],specfile.replace('.fits','_overview.pdf'),xrangefull=[1100,2200],zoomwindows=zoomwindows)

    """
    if plotSN:
        specfigure = outname.replace('.pdf','_SN.pdf')
    else:
        specfigure = outname
    if verbose: print(' - 1D overview figure will be saved to:\n   '+specfigure)
    if verbose: print(' - Will plot emission line markers at redshift '+str("%.6f" % redshift))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Nspec = len(spectra)
    if verbose: print(' - Loading the '+str(Nspec)+' spectra provided for plotting')
    datadic = collections.OrderedDict()

    for ss, specname in enumerate(spectra):
        spec_dat      = afits.open(specname)[1].data
        spec_wave     = spec_dat[wavecols[ss]]
        spec_flux     = spec_dat[fluxcols[ss]]
        spec_ferr     = spec_dat[fluxerrcols[ss]]
        try:
            Nspecstack  =  spec_dat['nspecstack']
        except:
            Nspecstack  =  None

        spec_S2N      = spec_flux/spec_ferr
        spec_filllow  = spec_flux-spec_ferr
        spec_fillhigh = spec_flux+spec_ferr
        spec_wavecov  = [np.min(spec_wave),np.max(spec_wave)]

        datadic[specname] = {'spec_wave':spec_wave,
                             'spec_flux':spec_flux,
                             'spec_ferr':spec_ferr,
                             'spec_S2N':spec_S2N,
                             'spec_filllow':spec_filllow,
                             'spec_fillhigh':spec_fillhigh,
                             'spec_wavecov':spec_wavecov,
                             'Nspecstack':Nspecstack}
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading line lists to display in plotting windows')
    llistdic_em    = MiGs.linelistdic(listversion='full')
    llistdic_abs   = MiGs.linelistdic(listversion='full_absorption')
    llistdic_fs    = MiGs.linelistdic(listversion='full_finestructure')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Nwindows     = len(zoomwindows)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Plotting figure ')
    figuresize_x = 10
    figuresize_y = 7
    fig          = plt.figure(figsize=(figuresize_x,figuresize_y))
    Fsize        = 10
    LW           = 2
    plt.rc('text', usetex=True)                         # enabling LaTex rendering of text
    plt.rc('font', family='serif',size=Fsize)           # setting text font
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    left   = 0.08   # the left side of the subplots of the figure
    right  = 0.92   # the right side of the subplots of the figure
    bottom = 0.10   # the bottom of the subplots of the figure
    top    = 0.97   # the top of the subplots of the figure
    wspace = 0.60   # the amount of width reserved for blank space between subplots
    hspace = 0.40   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    if speccols is None:
        speccols          = ['black','green','red','magenta','cyan','orange','purple','yellow']

    xlabel            = '$\lambda$ / [\AA]'
    ylabel            = '$f_\lambda / [10^{-20}$erg/s/cm$^2$/\\AA]'
    if plotSN:
        ylabel        = 'S/N'
    color_marker_em   = 'blue'
    color_marker_abs  = 'gray'
    color_marker_fs   = 'skyblue'

    wavescale         = 1.0
    nrows             = 1+np.ceil(Nwindows/3.0)
    ncols             = 3

    xrangedic    = {}
    yrangedic    = {}
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for ww, wininfo in enumerate(zoomwindows):
        windowmin   = wininfo[0]
        windowmax   = wininfo[1]
        windowname  = wininfo[2]

        good_em     = {}
        for linekey in llistdic_em.keys():
            if (llistdic_em[linekey][1] > windowmin) & (llistdic_em[linekey][1] < windowmax):
                good_em[linekey] = llistdic_em[linekey]

        good_abs     = {}
        for linekey in llistdic_abs.keys():
            if (llistdic_abs[linekey][1] > windowmin) & (llistdic_abs[linekey][1] < windowmax):
                good_abs[linekey] = llistdic_abs[linekey]

        good_fs     = {}
        for linekey in llistdic_fs.keys():
            if (llistdic_fs[linekey][1] > windowmin) & (llistdic_fs[linekey][1] < windowmax):
                good_fs[linekey] = llistdic_fs[linekey]

        plotindex = ww+1
        xrangedic[windowname], yrangedic[windowname] = \
            stacking.plot_1DspecOverview_subplot(nrows,ncols,plotindex,datadic,spectra,speccols,plotSN,
                                                 [good_em,good_abs,good_fs],Fsize,
                                                 [color_marker_em,color_marker_abs,color_marker_fs],
                                                 LW,xlabel,ylabel,windowmin,windowmax,wavescale,redshift)

        plt.plot([windowmin,windowmax],[0,0],alpha=0.5,color='black',linestyle='--',linewidth=LW/2.)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ax1 = plt.subplot(nrows, ncols, (ncols*nrows-2,ncols*nrows)) # Full specs
    yrange = stacking.plot_1DspecOverview_plotspecs(datadic,spectra,speccols,xrangefull,plotSN=plotSN,
                                                    labels=labels,tightyrange=False)
    ax1.plot(xrangefull,[0,0],alpha=0.5,color='black',linestyle='--',linewidth=LW/2.)

    plot_1DspecOverview_plotlines(llistdic_em,xrangefull[1]-xrangefull[0],None,#Fsize-2.0,
                                  color_marker_em,xrangefull,yrange,yrange[1]*0.9,redshift,LW/2.0)

    plot_1DspecOverview_plotlines(llistdic_abs,xrangefull[1]-xrangefull[0],None,#Fsize-2.0,
                                  color_marker_abs,xrangefull,yrange,yrange[1]*0.9,redshift,LW/2.0)

    plot_1DspecOverview_plotlines(llistdic_fs,xrangefull[1]-xrangefull[0],None,#Fsize-2.0,
                                  color_marker_fs,xrangefull,yrange,yrange[1]*0.9,redshift,LW/2.0)


    if not plotSN:
        if yrangefull is None:
            yrangefull  = yrange
    else:
        yrangefull  = yrange
    Dyrangefull = yrangefull[1]-yrangefull[0]

    # --- "ZOOM BOXES" ---
    boxzorder = 150

    for ww, wininfo in enumerate(zoomwindows):
        windowname = wininfo[2]
        stacking.plot_1DspecOverview_genbox(windowname,xrangedic[windowname],yrangedic[windowname],
                                            LW,boxzorder,Dyrangefull,'black',Fsize)

    #--------- LEGEND ---------
    anchorpos = (0.5, 1.2)
    if len(labels) > 6:
        ncols = 6
    else:
        ncols = len(labels)
    leg = plt.legend(fancybox=True,numpoints=1, loc='upper center',prop={'size':Fsize},ncol=ncols)#,
                     #bbox_to_anchor=anchorpos)  # add the legend
    #leg.get_frame().set_alpha(0.7)
    #--------------------------

    for ss, specname in enumerate(spectra):
        if datadic[specname]['Nspecstack'] is not None:
            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
            ax2.plot(datadic[specname]['spec_wave'],
                     datadic[specname]['Nspecstack'],alpha=0.5, color=speccols[ss],zorder=1, drawstyle='steps' )
            ax2.tick_params(axis='y',labelcolor=speccols[ss])
            ax2.set_ylabel('Spectra per pixel',alpha=0.5,color=speccols[ss])
            ax2.set_xlabel(xlabel)
            ax2.set_xlim(xrangefull)


    ax1.set_xlim(xrangefull)
    ax1.set_ylim(yrangefull)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.savefig(specfigure) # dpi = dot per inch for rasterized points
    fig.clf()
    plt.close('all')
    if verbose: print(' - Saved figure to '+specfigure)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview_subplot(nrows,ncols,plotindex,datadic,spectra,speccols,plotSN,llistdics,Fsize,linecolors,LW,
                                xlabel,ylabel,windowmin,windowmax,wavescale,redshift):
    """
    Function generating sub plots for plot_1DspecOverview()
    """
    ax1 = plt.subplot(nrows,ncols,plotindex)
    xrange       = np.asarray([windowmin,windowmax])*(redshift+1)/wavescale
    windowwidth  = windowmax-windowmin

    yrange = stacking.plot_1DspecOverview_plotspecs(datadic,spectra,speccols,xrange,plotSN=plotSN,labels=None)

    for ll, llist in enumerate(llistdics):
        stacking.plot_1DspecOverview_plotlines(llist,windowwidth,Fsize-2,linecolors[ll],
                                               xrange,yrange,yrange[1]*(0.95-0.1*ll),redshift,LW,wavetype='vac')

    ax1.set_ylim(yrange)
    ax1.set_ylabel(ylabel)

    for ss, specname in enumerate(spectra):
        if datadic[specname]['Nspecstack'] is not None:
            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

            waveent = (datadic[specname]['spec_wave'] > xrange[0]) & (datadic[specname]['spec_wave'] < xrange[1])
            ax2.plot(datadic[specname]['spec_wave'][waveent],
                     datadic[specname]['Nspecstack'][waveent], alpha=0.5, color=speccols[ss],zorder=1, drawstyle='steps' )
            ax2.tick_params(axis='y',labelcolor=speccols[ss])
            ax2.set_ylabel('Spectra per pixel', alpha=0.5,color=speccols[ss])
            ax2.set_xlabel(xlabel)
            ax2.set_xlim(xrange)
    else:
        ax1.set_xlabel(xlabel)
        ax1.set_xlim(xrange)

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
def plot_1DspecOverview_plotspecs(datadic,spectra,colors,xrange,plotSN=False,labels=None,tightyrange=False):
    """

    """
    if labels is None:
        labels = [None]*len(spectra)
    yrangecomb = [0,10]

    for ss, specname in enumerate(spectra):

        waveent = (datadic[specname]['spec_wave'] > xrange[0]) & (datadic[specname]['spec_wave'] < xrange[1])

        if plotSN:
            plt.plot(datadic[specname]['spec_wave'][waveent],
                     datadic[specname]['spec_S2N'][waveent],
                     '-',alpha=0.8,color=colors[ss],label=labels[ss],zorder=100)

            try:
                specvals = datadic[specname]['spec_S2N'][waveent]
                fluxmin  = np.min(np.asarray( [0,  np.min(specvals[np.isfinite(specvals)])] ))
                fluxmax  = np.max(np.asarray( [5, np.max(specvals[np.isfinite(specvals)])] ))
                if tightyrange:
                    yrange  = [fluxmin,fluxmax*1.1]
                else:
                    dy      = fluxmax-fluxmin
                    yrange  = [fluxmin-np.abs(dy*0.2),fluxmax+np.abs(dy*0.5)]
            except:
                yrange = [0,10]

        else:
            try:
                specvals = datadic[specname]['spec_flux'][waveent]
                fluxmin  = np.min(np.asarray( [0,  np.min(specvals[np.isfinite(specvals)])] ))
                fluxmax  = np.max(np.asarray( [5, np.max(specvals[np.isfinite(specvals)])] ))
                if tightyrange:
                    yrange  = [fluxmin,fluxmax*1.1]
                else:
                    dy      = fluxmax-fluxmin
                    yrange  = [fluxmin-np.abs(dy*0.2),fluxmax+np.abs(dy*0.5)]
            except:
                yrange = [0,100]

            xvalues = datadic[specname]['spec_wave'][waveent]
            yvalues = datadic[specname]['spec_flux'][waveent]
            plt.plot(xvalues,yvalues,'-',alpha=0.8,color=colors[ss],label=labels[ss],zorder=100)

            ylow  = datadic[specname]['spec_filllow'][waveent]
            yhigh = datadic[specname]['spec_fillhigh'][waveent]
            yhigh[yhigh > yrange[1]] = yrange[1]
            ylow[ylow  < yrange[0]] = yrange[0]
            plt.fill_between(xvalues,ylow,yhigh,alpha=0.20,color=colors[ss],zorder=100)

        if yrange[0] < yrangecomb[0]:
            yrangecomb[0] = yrange[0]

        if yrange[1] > yrangecomb[1]:
            yrangecomb[1] = yrange[1]

    for ss, specname in enumerate(spectra):
        if plotSN:
            waveent = (datadic[specname]['spec_wave'] > xrange[0]) & (datadic[specname]['spec_wave'] < xrange[1])
            meanerr = np.mean(datadic[specname]['spec_ferr'])
        else:
            meanerr = 1.0

    # if yrangecomb == [0,100]: pdb.set_trace()
    return yrangecomb
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1DspecOverview_plotlines(llistdic,windowwidth,Fsize,col_linemarker,xrange,yrange,textypos,redshift,LW,wavetype='vac'):
    """

    """
    for ll in llistdic.keys():
        linedat      = llistdic[ll]
        linename     = linedat[0]
        if wavetype == 'vac':
            linewave = linedat[1]
        else:
            linewave = kbs.convert_wavelength(linedat[1],version='vac2air')
            pdb.set_trace()
        horalign     = linedat[2]
        lineposition = linewave*(redshift+1.0)

        if (lineposition > xrange[0]) & (lineposition < xrange[1]):
            plt.plot(np.zeros(2)+lineposition,yrange,color=col_linemarker,alpha=0.7,linestyle='-',linewidth=LW,zorder=20)
            if horalign == 'right':
                xpos = lineposition-0.03*windowwidth
            elif horalign == 'left':
                xpos = lineposition+0.03*windowwidth
            else:
                xpos = lineposition

            ypos = textypos
            if Fsize is not None:
                plt.text(xpos,ypos,linename,color=col_linemarker,size=Fsize,
                         rotation='horizontal',horizontalalignment=horalign,verticalalignment='top',zorder=20)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def save_stack_1D(outfile,wave,flux,fluxerr,nspecstack,headerinfo=None,overwrite=False,verbose=True):
    """
    Function saving spectrum to fits file

    --- INPUT ---
    outfile         Name of output file to generate
    wave            Wavelength in units of Angstrom
    flux            Flux to store to fits file
    fluxerr         Error on flux (sqrt(variance))
    nspecstack      Number of spectra contributing to each pixel in stacked spectrum
    headerinfo      To add info to the header provide it to this keyword as a dictionary on the format:
                       headerinfo = {'KEYNAME1':[VALUE1,INFOCOMMENT1], 'KEYNAME2':[VALUE2,INFOCOMMENT2]}
    overwrite       Overwrite existing file?
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---

    """
    S2N = flux/fluxerr
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Saving wavelenght and flux values to \n   '+outfile)
    mainHDU = afits.PrimaryHDU()       # primary HDU
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    c1 = afits.Column(name='wave',      format='D', unit='ANGSTROMS', array=wave)
    c2 = afits.Column(name='flux',      format='D', unit='', array=flux)
    c3 = afits.Column(name='fluxerror', format='D', unit='', array=fluxerr)
    c4 = afits.Column(name='s2n',       format='D', unit='', array=S2N)
    c5 = afits.Column(name='nspecstack',format='D', unit='', array=nspecstack)

    coldefs = afits.ColDefs([c1,c2,c3,c4,c5])
    tbHDU   = afits.BinTableHDU.from_columns(coldefs) # creating default header

    # writing hdrkeys:'---KEY--',                             '----------------MAX LENGTH COMMENT-------------'
    tbHDU.header.append(('EXTNAME ','SPEC1D'                     ,'cube containing source'),end=True)
    if headerinfo is not None:
        for key in headerinfo.keys():
            tbHDU.header.append((key,headerinfo[key][0],headerinfo[key][1]),end=True)

    hdulist = afits.HDUList([mainHDU,tbHDU])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    hdulist.writeto(outfile, overwrite=overwrite)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =