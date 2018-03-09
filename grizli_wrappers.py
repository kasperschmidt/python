# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Scripts, functions and wrappers for handling grizly reductions and simulations
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import glob
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
import sys
import grizli
import grizli.utils
from grizli.multifit import GroupFLT, MultiBeam, get_redshift_fit_defaults
import pdb
from astropy.io import ascii
from importlib import reload
import JADESutilities as ju

# plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec
mpl.rcParams['figure.figsize'] = (10.0, 6.0)
mpl.rcParams['font.size'] = 14
mpl.rcParams['savefig.dpi'] = 72

#importlib.reload(module)

# for grism sims
import grizli.fake_image
from grizli.utils import detect_with_photutils
from collections import OrderedDict
from astropy.io import fits
import astropy.wcs
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import drizzlepac
import pysynphot
import grizli_wrappers as gw

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def a2744_glassdatareduction(workdir='/Users/kschmidt/work/JWST/grizly_A2744/Prep',clobber=True):
    """
    Convenience wrapper to build and generate all the files needed for the TDOSE run

    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw
    gw.a2744_glassdatareduction(workdir='/Users/kschmidt/work/JWST/grizly_A2744/Prep180118',clobber=True)


    """
    print(' - Reducing A2744 GLASS data using ')
    print(' - Grizli version '+grizli.__version__)
    try:
        os.chdir(workdir)
        print(' - Moved to directory '+os.getcwd())
    except:
        sys.exit('Working directory '+workdir+' does not exists')

    files = glob.glob('../RAW/*flt.fits')
    info = grizli.utils.get_flt_info(files)
    ascii.write(info,'./rawFLTfileinfo.txt')
    visits, filters = grizli.utils.parse_flt_files(info=info, use_visit=False, uniquename=True)

    for i in range(len(visits)):
        print(dict(visits[i]))

    gotoextraction = True  # <------------------------------------------------------------------- Keyword for skipping
    if not gotoextraction:
        print(' ----------------------- Master Pre-processing ----------------------- ')
        from grizli.prep import process_direct_grism_visit
        print(' - defining visit pairs to process (NEEDS TO BE UPDATED TO INCLUE ALL OBS)')
        visitpairs = [[visits[0],visits[4]],[visits[2],visits[6]],[visits[7],visits[11]], [visits[9],visits[13]]]

        runmaster = False # <--------------------------------------------------------------------- Keyword for skipping
        if runmaster:
            for vv, visit in enumerate(visits):
                print(' ----- Visit No. '+str(vv+1)+' ----- ')
                visitfiles = visits[vv]['files']
                print(' Files: '+str(visitfiles))
                for vf in visitfiles:
                    infoent = np.where(info['FILE'] == vf)[0]
                    print('  '+vf+' infofilter = '+str(info['FILTER'][infoent][0]))

            for vp in visitpairs:
                status = process_direct_grism_visit(direct=vp[0], grism=vp[1],
                                                radec='../hst_a2744_60mas_merged_radec_F140Wlt24.txt',
                                                align_mag_limits=[14,23])
        else:
            print('Skipping master pre-processing')

        shiftlogs = glob.glob('*shifts.log')
        print(' - For info on shifts, check '+''.join(shiftlogs))
        # !ls *shifts.log
        # print('')
        # !cat *shifts.log

        print(' ----------------------- Grouping FLTs ----------------------- ')
        all_grism_files = []
        grismvisits     = [vp[1] for vp in visitpairs]
        for i in range(len(grismvisits)):
            all_grism_files.extend(grismvisits[i]['files'])

        print(' - Grism files (all_grism_files) to group: '+str(all_grism_files))
        refimgpath = '/Users/kschmidt/work/images_MAST/images_fullfov/'
        grp = GroupFLT(grism_files=all_grism_files, direct_files=[],
                      ref_file=refimgpath+'refimage_hlsp_frontier_hst_wfc3-60mas_abell2744_f140w_v1.0_drz.fits',
                      seg_file=refimgpath+'refimage_hlsp_frontier_hst_wfc3-60mas_abell2744_f140w_v1.0_drz_seg.fits',
                      catalog='/Users/kschmidt/work/catalogs/GLASScatalogs/GLASScatalog_A2744_150515.cat',
                      cpu_count=8)


        print(' ----------------------- Generate Continuum Model ----------------------- ')
        genmodels = False    # <------------------------------------------------------------------- Keyword for skipping
        if genmodels:
            grp.compute_full_model(mag_limit=25)

            print(' - Plotting: Show FLT residuals')
            fig = plt.figure(figsize=[12,6])
            ax = fig.add_subplot(121)
            ax.imshow(grp.FLTs[0].grism['SCI'] - grp.FLTs[0].model, vmin=-0.02, vmax=0.2, cmap='cubehelix_r',
                      interpolation='Nearest', origin='lower')
            ax.set_title('G102, %s' %(grp.FLTs[0].grism.parent_file))

            ax = fig.add_subplot(122)
            ax.imshow(grp.FLTs[4].grism['SCI'] - grp.FLTs[4].model, vmin=-0.02, vmax=0.2, cmap='cubehelix_r',
                      interpolation='Nearest', origin='lower')
            ax.set_title('G141, %s' %(grp.FLTs[4].grism.parent_file))

            for ax in fig.axes:
                #ax.set_xlim(500,700); ax.set_ylim(500,700)
                ax.set_xlim(100,1500); ax.set_ylim(100,1500)
            fig.savefig('./FLTresiduals_continuum_model.pdf')

        else:
            print('Skipping continuum model')


        print(' ----------------------- Generate Polynomial Model ----------------------- ')
        if genmodels:
            grp.refine_list(poly_order=2, mag_limits=[16, 24], verbose=False)


            print(' - Plotting: Show FLT residuals')
            fig = plt.figure(figsize=[12,6])
            ax = fig.add_subplot(121)
            ax.imshow(grp.FLTs[0].grism['SCI'] - grp.FLTs[0].model, vmin=-0.02, vmax=0.2, cmap='cubehelix_r',
                      interpolation='Nearest', origin='lower')
            ax.set_title('G102, %s' %(grp.FLTs[0].grism.parent_file))

            ax = fig.add_subplot(122)
            ax.imshow(grp.FLTs[4].grism['SCI'] - grp.FLTs[4].model, vmin=-0.02, vmax=0.2, cmap='cubehelix_r',
                      interpolation='Nearest', origin='lower')
            ax.set_title('G141, %s' %(grp.FLTs[4].grism.parent_file))

            for ax in fig.axes:
                #ax.set_xlim(500,700); ax.set_ylim(500,700)
                ax.set_xlim(100,1500); ax.set_ylim(100,1500)
            fig.savefig('./FLTresiduals_polynomial_model.pdf')
        else:
            print('Skipping Polynomial model')


        print(' ----------------------- Save Models ----------------------- ')
        if genmodels:
            grp.save_full_data()
        else:
            print(' did not generate models, so did not save models to disk')

    else:
        print('\n   NB\n - Going directly to spectral extraction...')

    print(' ----------------------- Prepare Fitting Spectra ----------------------- ')
    all_grism_files = ['ica501u3q_flt.fits', 'ica501uaq_flt.fits', 'ica501uhq_flt.fits', 'ica501uoq_flt.fits',
                       'ica501tbq_flt.fits', 'ica501tiq_flt.fits', 'ica501tpq_flt.fits', 'ica501twq_flt.fits',
                       'ica503fwq_flt.fits', 'ica503g3q_flt.fits', 'ica503gaq_flt.fits', 'ica503ghq_flt.fits',
                       'ica503ecq_flt.fits', 'ica503ejq_flt.fits', 'ica503eqq_flt.fits', 'ica503f5q_flt.fits']

    print(' - Grism files (all_grism_files) to group: '+str(all_grism_files))
    refimgpath = '/Users/kschmidt/work/images_MAST/images_fullfov/'
    grp = GroupFLT(grism_files=all_grism_files, direct_files=[],
                  ref_file=refimgpath+'refimage_hlsp_frontier_hst_wfc3-60mas_abell2744_f140w_v1.0_drz.fits',
                  seg_file=refimgpath+'refimage_hlsp_frontier_hst_wfc3-60mas_abell2744_f140w_v1.0_drz_seg.fits',
                  catalog='/Users/kschmidt/work/catalogs/GLASScatalogs/GLASScatalog_A2744_150515.cat',
                  cpu_count=8)

    print(' ----------------------- Setting up templates for fits ----------------------- ')
    # First is set with combined emission line complexes for the redshift fit
    # (don't allow infinite freedom) of the line ratios / fluxes
    templ0 = grizli.utils.load_templates(fwhm=1200, line_complexes=True, stars=False,
                                         full_line_list=None,  continuum_list=None,
                                         fsps_templates=True)

    # Second set has individual line templates for fitting the line fluxes
    templ1 = grizli.utils.load_templates(fwhm=1200, line_complexes=False, stars=False,
                                         full_line_list=None, continuum_list=None,
                                         fsps_templates=True)

    # Show the template names / dictionary keys
    fmt = '{0:<36s} {1:<36s}'
    print(fmt.format('templ0', 'templ1'))
    print(fmt.format('------', '------'))

    for i in range(len(templ1)):
        if i > len(templ0)-1:
            print(fmt.format('', list(templ1.keys())[i]))
        else:
            print(fmt.format(list(templ0.keys())[i], list(templ1.keys())[i]))

    # Parameters for drizzled line maps
    pline = {'kernel': 'point', 'pixfrac': 0.2, 'pixscale': 0.1, 'size': 8, 'wcs': None}


    print(' ----------------------- Pull out individual objects ----------------------- ')
    # grp `GroupFLT` object created as defined in WFC3IR_Reduction from the WFC3 ERS grism data
    target = 'glass_a2744'

    #      ELs      Cont      Cont
    ids = [161,     316,      694]

    for id in ids:
        # Pull out the 2D cutouts
        beams = grp.get_beams(id, size=80)
        mb = grizli.multifit.MultiBeam(beams, fcontam=0.5, group_name=target, psf=False)

        # Save a FITS file with the 2D cutouts (beams) from the individual exposures
        mb.write_master_fits()

        # Fit polynomial model for initial continuum subtraction
        wave = np.linspace(2000,2.5e4,100)
        poly_templates = grizli.utils.polynomial_templates(wave, order=7)
        pfit = mb.template_at_z(z=0, templates=poly_templates, fit_background=True,
                                fitter='lstsq', get_uncertainties=2)

        # Drizzle grisms / PAs
        hdu, fig = mb.drizzle_grisms_and_PAs(fcontam=0.2, flambda=False, kernel='point',
                                             size=32, zfit=pfit)

        # Save drizzle figure FITS file
        fig.savefig('{0}_{1:05d}.stack.png'.format(target, id))
        hdu.writeto('{0}_{1:05d}.stack.fits'.format(target, id), clobber=True)


        print(' ----------------------- Run wrapper on object '+str(id)+'----------------------- ')
        # High-level wrapper script for doing everything (redshift fits, line fluxes, drizzled line
        # maps).  More explanation of the details of individual steps TBD.
        #
        # Needs to be able to find {target}_{id:05d}.beams.fits and {target}_{id:05d}.stack.fits
        # generated above
        out = grizli.fitting.run_all(id, t0=templ0, t1=templ1, fwhm=1200,
                                     zr=[0.1, 1.7], dz=[0.004, 0.0005],
                                     fitter='nnls', group_name=target, prior=None, fcontam=0.,
                                     pline=pline, mask_sn_limit=7, fit_beams=True, fit_stacks=False,
                                     root=target+'_', fit_trace_shift=False, verbose=False,
                                     phot=None, scale_photometry=False, show_beams=True)


        print(" - For analyzing fit output see\n"
              "   https://github.com/gbrammer/grizli/blob/master/examples/NewSpectrumFits.ipynb")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def GOODSSERS_examplereduction():
    """

    Trying to follow the example at
        https://github.com/gbrammer/grizli/blob/master/examples/WFC3IR_Reduction.ipynb

    and (for the spectral extraction and analysis):
        https://github.com/gbrammer/grizli/blob/master/examples/NewSpectrumFits.ipynb

    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw
    gw.GOODSSERS_examplereduction()

    """
    print('Reducing GOODS-S ERS data using demo from Grizli GtiHub page')
    print(grizli.__version__)
    print(os.getcwd())

    files = glob.glob('../RAW/*flt.fits')
    info = grizli.utils.get_flt_info(files)
    visits, filters = grizli.utils.parse_flt_files(info=info, uniquename=True)

    for vv, visit in enumerate(visits):
        print(' ----- Visit No. '+str(vv+1)+' ----- ')
        visitfiles = visits[vv]['files']
        print(' Files: '+str(visitfiles))
        for vf in visitfiles:
            infoent = np.where(info['FILE'] == vf)[0]
            print('  '+vf+' infofilter = '+str(info['FILTER'][infoent][0]))

    for i in range(4):
        print(dict(visits[i]))

    gotofitting = True
    if not gotofitting:
        print(' ----------------------- Master Pre-processing ----------------------- ')
        from grizli.prep import process_direct_grism_visit

        runmaster = False
        if runmaster:
            for i in range(2):
                status = process_direct_grism_visit(direct=visits[i], grism=visits[i+2], # specific to the particular order for
                                                radec='../Catalog/goodss_radec.dat',     # the demo
                                                align_mag_limits=[14,23])
        else:
            print('Skipping master pre-processing')

        shiftlogs = glob.glob('*shifts.log')
        print(' - For info on shifts, check '+''.join(shiftlogs))
        # !ls *shifts.log
        # print('')
        # !cat *shifts.log


        print(' ----------------------- Grouping FLTs ----------------------- ')
        from grizli.multifit import GroupFLT, MultiBeam, get_redshift_fit_defaults


        all_grism_files = []
        for i in range(len(visits)):
            if '-g1' in visits[i]['product']:
                all_grism_files.extend(visits[i]['files'])

        grp = GroupFLT(grism_files=all_grism_files, direct_files=[],
                      ref_file='../Catalog/ERS_goodss_3dhst.v4.0.F160W_orig_sci.fits',
                      seg_file='../Catalog/ERS_GOODS-S_IR.seg.fits',
                      catalog='../Catalog/ERS_GOODS-S_IR.cat',
                      cpu_count=8)


        print(' ----------------------- Generate Continuum Model ----------------------- ')
        genmodels = True
        if genmodels:
            grp.compute_full_model(mag_limit=25)

            print(' - Plotting: Show FLT residuals')
            fig = plt.figure(figsize=[12,6])
            ax = fig.add_subplot(121)
            ax.imshow(grp.FLTs[0].grism['SCI'] - grp.FLTs[0].model, vmin=-0.02, vmax=0.2, cmap='cubehelix_r',
                      interpolation='Nearest', origin='lower')
            ax.set_title('G102, %s' %(grp.FLTs[0].grism.parent_file))

            ax = fig.add_subplot(122)
            ax.imshow(grp.FLTs[4].grism['SCI'] - grp.FLTs[4].model, vmin=-0.02, vmax=0.2, cmap='cubehelix_r',
                      interpolation='Nearest', origin='lower')
            ax.set_title('G141, %s' %(grp.FLTs[4].grism.parent_file))

            for ax in fig.axes:
                #ax.set_xlim(500,700); ax.set_ylim(500,700)
                ax.set_xlim(100,1500); ax.set_ylim(100,1500)
            fig.savefig('./FLTresiduals_continuum_model.pdf')

        else:
            print('Skipping continuum model')


        print(' ----------------------- Generate Polynomial Model ----------------------- ')
        if genmodels:
            grp.refine_list(poly_order=2, mag_limits=[16, 24], verbose=False)


            print(' - Plotting: Show FLT residuals')
            fig = plt.figure(figsize=[12,6])
            ax = fig.add_subplot(121)
            ax.imshow(grp.FLTs[0].grism['SCI'] - grp.FLTs[0].model, vmin=-0.02, vmax=0.2, cmap='cubehelix_r',
                      interpolation='Nearest', origin='lower')
            ax.set_title('G102, %s' %(grp.FLTs[0].grism.parent_file))

            ax = fig.add_subplot(122)
            ax.imshow(grp.FLTs[4].grism['SCI'] - grp.FLTs[4].model, vmin=-0.02, vmax=0.2, cmap='cubehelix_r',
                      interpolation='Nearest', origin='lower')
            ax.set_title('G141, %s' %(grp.FLTs[4].grism.parent_file))

            for ax in fig.axes:
                #ax.set_xlim(500,700); ax.set_ylim(500,700)
                ax.set_xlim(100,1500); ax.set_ylim(100,1500)
            fig.savefig('./FLTresiduals_polynomial_model.pdf')
        else:
            print('Skipping Polynomial model')


        print(' ----------------------- Save Models ----------------------- ')
        if genmodels:
            grp.save_full_data()
        else:
            print(' did not generate models, so did not save models to disk')

        grp = GroupFLT(grism_files=all_grism_files, direct_files=[],
                      ref_file='../Catalog/ERS_goodss_3dhst.v4.0.F160W_orig_sci.fits',
                      seg_file='../Catalog/ERS_GOODS-S_IR.seg.fits',
                      catalog='../Catalog/ERS_GOODS-S_IR.cat',
                      cpu_count=8)
    else:
        print('\n   NB\n - Going directly to the fitting part...')

    print(' ----------------------- Fitting Spectra ----------------------- ')
    all_grism_files = ['ib6o21qmq_flt.fits', 'ib6o21qoq_flt.fits', 'ib6o21r6q_flt.fits',
                       'ib6o21r8q_flt.fits', 'ib6o23rsq_flt.fits', 'ib6o23ruq_flt.fits',
                       'ib6o23ryq_flt.fits','ib6o23s0q_flt.fits']

    grp = grizli.multifit.GroupFLT(grism_files=all_grism_files, direct_files=[],
                  ref_file='../Catalog/ERS_goodss_3dhst.v4.0.F160W_orig_sci.fits',
                  seg_file='../Catalog/ERS_GOODS-S_IR.seg.fits',
                  catalog='../Catalog/ERS_GOODS-S_IR.cat',
                  cpu_count=8)

    print(' ----------------------- Fitting templates ----------------------- ')
    # First is set with combined emission line complexes for the redshift fit
    # (don't allow infinite freedom) of the line ratios / fluxes
    templ0 = grizli.utils.load_templates(fwhm=1200, line_complexes=True, stars=False,
                                         full_line_list=None,  continuum_list=None,
                                         fsps_templates=True)

    # Second set has individual line templates for fitting the line fluxes
    templ1 = grizli.utils.load_templates(fwhm=1200, line_complexes=False, stars=False,
                                         full_line_list=None, continuum_list=None,
                                         fsps_templates=True)

    # Show the template names / dictionary keys
    fmt = '{0:<36s} {1:<36s}'
    print(fmt.format('templ0', 'templ1'))
    print(fmt.format('------', '------'))

    for i in range(len(templ1)):
        if i > len(templ0)-1:
            print(fmt.format('', list(templ1.keys())[i]))
        else:
            print(fmt.format(list(templ0.keys())[i], list(templ1.keys())[i]))

    # Parameters for drizzled line maps
    pline = {'kernel': 'point', 'pixfrac': 0.2, 'pixscale': 0.1, 'size': 8, 'wcs': None}

    print(' ----------------------- Pull out individual objects ----------------------- ')
    # grp `GroupFLT` object created as defined in WFC3IR_Reduction from the WFC3 ERS grism data
    target = 'ers-grism'

    # Line-emitter
    id=40776

    # Weak continuum features
    #id=41147

    # Strong continuum
    #id=43114

    # Pull out the 2D cutouts
    beams = grp.get_beams(id, size=80)
    mb = grizli.multifit.MultiBeam(beams, fcontam=0.5, group_name=target, psf=False)

    # Save a FITS file with the 2D cutouts (beams) from the individual exposures
    mb.write_master_fits()

    # Fit polynomial model for initial continuum subtraction
    wave = np.linspace(2000,2.5e4,100)
    poly_templates = grizli.utils.polynomial_templates(wave, order=7)
    pfit = mb.template_at_z(z=0, templates=poly_templates, fit_background=True,
                            fitter='lstsq', get_uncertainties=2)

    # Drizzle grisms / PAs
    hdu, fig = mb.drizzle_grisms_and_PAs(fcontam=0.2, flambda=False, kernel='point',
                                         size=32, zfit=pfit)

    # Save drizzle figure FITS file
    fig.savefig('{0}_{1:05d}.stack.png'.format(target, id))
    hdu.writeto('{0}_{1:05d}.stack.fits'.format(target, id), clobber=True)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def NIRISSsim_UDF():
    """

    Trying to follow the example at
        https://github.com/gbrammer/grizli/blob/master/examples/NIRISS-simulation.ipynb

    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw
    gw.NIRISSsim_UDF()

    """
    cwd = os.getcwd()
    print('\n grizli version: %s' %(grizli.__version__))
    print('Simulating NIRISS exposure on UDF')

    print(' - The 3D-HST catalog http://www.stsci.edu/~brammer/Grizli/Demos/udf_3dhst_cat.fits was manually download to cwd = '+cwd)
    imagepath = '/Users/kschmidt/work/images_MAST/'

    generatesimulation = True  # <------------------------------------------------------------ Keyword for skipping
    if generatesimulation:
        print(' - Generating simulations')
        matchCatAndMakeSeg = True  # <------------------------------------------------------------ Keyword for skipping
        if matchCatAndMakeSeg:
            print(' - Make an object catalog / segmentation image with photutils')
            print('   Loading XDF F140 images from '+imagepath)
            sci = fits.open(imagepath+'hlsp_xdf_hst_wfc3ir-60mas_hudf_f140w_v1_sci.fits')
            wht = fits.open(imagepath+'hlsp_xdf_hst_wfc3ir-60mas_hudf_f140w_v1_wht.fits')
            rms = 1/np.sqrt(wht[0].data)
            rms[~np.isfinite(rms)] = 1000
            dq = rms > 900
            wcs = astropy.wcs.WCS(sci[0].header)

            print(' - Run the detection with the grizli / photutils wrapper')
            cat, seg = detect_with_photutils(sci[0].data, err=rms, dq=dq, seg=None,
                                             detect_thresh=1.4, npixels=6, grow_seg=3,
                                             gauss_fwhm=2.0, gsize=3, wcs=wcs, save_detection=False,
                                             root='udf_f140w_photutils', background=None, gain=None,
                                             AB_zeropoint=26.452,
                                             rename_columns={'ycentroid': 'y_flt', 'xcentroid': 'x_flt',
                                                             'dec_icrs_centroid': 'dec',
                                                             'ra_icrs_centroid': 'ra'},
                                             clobber=True, verbose=True)

            # code expects this columns later....
            cat['NUMBER'] = cat['id']

            print(' - Get matches from 3D-HST catalog')
            ref_3dhst = fits.open('udf_3dhst_cat.fits')
            ref_cat = Table.read(ref_3dhst[1])

            gs = SkyCoord(ra=ref_cat['ra']*u.degree, dec=ref_cat['dec']*u.degree)
            cat_rd = SkyCoord(ra=cat['ra'], dec=cat['dec'])

            gs_idx, d2d, d3d = cat_rd.match_to_catalog_sky(gs)
            has_gs_match = np.where(d2d < 2*u.arcsec)[0]

            # Use 3D-HST mags because quick photutils catalog has some issues,
            # perhaps related to background subtraction
            gs_mag = 25-2.5*np.log10(ref_cat['F204'])
            cat['MAG_AUTO'] = gs_mag[gs_idx]
            cat.write('udf_f140w_photutils.cat', format='ascii.commented_header')

            fits.writeto('udf_f140w_photutils_seg.fits', data=np.cast[int](seg),
                         header=sci[0].header, clobber=True)


        print(' - Setup fake (noise) images, centered in the UDF/XDF')
        ra, dec = 53.1592277508136, -27.782056346146
        pa_aper = 128.589

        np.random.seed(1)
        # Crude exposure parameters.  Rough read noise and backgrounds coded in
        # grizli.fake_image to make some reasonable noise estimate
        EXPTIME = 1.e4 # 10 ks ~ 4 HST orbits
        NEXP = 10      # divided between 10 exposures

        # JWST NIRISS, three filters & two orients
        for filt in ['F115W', 'F150W', 'F200W']:
            for theta in [0,90]:
                h, wcs = grizli.fake_image.niriss_header(filter=filt, ra=ra, dec=dec,
                                                         pa_aper=pa_aper+theta)
                print('Filter: {filter}, Background: {bg} e/s/pix, RN: {RN} e/exp'.format(filter=filt,
                                                                bg=h['BACKGR'], RN=h['READN']))
                output = 'niriss_{filt}_{theta:02d}_flt.fits'.format(filt=filt, theta=theta)
                grizli.fake_image.make_fake_image(h, output=output, exptime=EXPTIME, nexp=NEXP)


        print(' - Load GroupFLT for simulation, NB: input files are just noise')
        sim = grizli.multifit.GroupFLT(grism_files=glob.glob('niriss_*flt.fits'), direct_files=[],
                                       ref_file=imagepath+'hlsp_xdf_hst_wfc3ir-60mas_hudf_f140w_v1_sci.fits',
                                       seg_file='udf_f140w_photutils_seg.fits',
                                       catalog='udf_f140w_photutils.cat',
                                       cpu_count=0, pad=200)


        print('   First pass, flat continuum ')
        sim.compute_full_model(mag_limit=27)

        print('   Compute model grism spectra for 3D-HST matches based on full photo-z templates')
        detection_bp = pysynphot.ObsBandpass('wfc3,ir,f140w')
        for ix in has_gs_match:
            templ = ref_3dhst['WAVE'].data*(1+ref_cat['zbest'][gs_idx[ix]])
            tempf = ref_3dhst['FLAMBDA'].data[gs_idx[ix],:]
            # Needs to be normalized to unity in the detection band
            spec = pysynphot.ArraySpectrum(wave=templ, flux=tempf, waveunits='angstroms', fluxunits='flam')
            spec = spec.renorm(1., 'flam', detection_bp)

            id = cat['id'][ix]
            #print(id)
            sim.compute_single_model(id, mag=cat['MAG_AUTO'][ix], size=-1, store=False,
                                     spectrum_1d=[spec.wave, spec.flux], get_beams=None,
                                     in_place=True)

        print(' - Plot blotted reference image; "grism" exposures are still empty, just noise')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate([0,2,4,1,3,5]):
            ax = fig.add_subplot(2,3,ix+1)
            ax.imshow(sim.FLTs[i].direct['REF'], vmin=-0.01, vmax=0.05, cmap='viridis',
                      origin='lower')
            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.grid(color='w', alpha=0.8)
            ax.text(100,100,sim.FLTs[i].grism_file, color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./blotted_reference_image.pdf')

        print(' - Plot blotted segmentation image; "grism" exposures are still empty, just noise')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate([0,2,4,1,3,5]):
            ax = fig.add_subplot(2,3,ix+1)
            ax.imshow(sim.FLTs[i].seg, vmin=-0.01, vmax=3000, cmap='viridis', origin='lower')
            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.grid(color='w', alpha=0.8)
            ax.text(100,100,sim.FLTs[i].grism_file, color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./blotted_segmentation_image.pdf')

        print(' - Plot "grism" exposures that are still empty, just noise')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate([0,2,4,1,3,5]):
            ax = fig.add_subplot(2,3,ix+1)
            ax.imshow(sim.FLTs[i].grism['SCI'], vmin=-0.01, vmax=0.05, cmap='viridis', origin='lower')
            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.grid(color='w', alpha=0.8)
            ax.text(100,100,sim.FLTs[i].grism_file, color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./blotted_grismnoise_image.pdf')

        print(' - Plot model stored in FLTs[i].model attribute')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate([0,2,4,1,3,5]):
            ax = fig.add_subplot(2,3,ix+1)

            # show as if it were the rotated grism
            if (i % 2) > 0:
                ax.imshow(np.rot90(sim.FLTs[i].model,-1), vmin=-0.01, vmax=0.05, cmap='viridis',
                          origin='lower')
            else:
                ax.imshow(sim.FLTs[i].model, vmin=-0.01, vmax=0.05, cmap='viridis', origin='lower')

            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.grid(color='w', alpha=0.8)
            ax.text(100,100,sim.FLTs[i].grism_file, color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./flt_model_attribute.pdf')

        print(' - Update SCI extension of the fake FLT images with the models just computed')
        for flt in sim.FLTs:
            print('Update', flt.grism_file)
            orig_flt = fits.open(flt.grism_file, mode='update')
            orig_flt['SCI'].data += flt.model[flt.pad:-flt.pad, flt.pad:-flt.pad]
            orig_flt.flush()

        # to make sure that only sensible flux values are plotted us the following:
        # fig = plt.figure(figsize=[5,5])
        # ax = fig.add_subplot(1,1,1)
        #
        # #modimg  = sim.FLTs[i].direct['REF']
        # modimg  = sim.FLTs[i].model
        # goodent = (np.abs(modimg) < 1000.0)
        # imgarr  = modimg*0.0
        # imgarr[goodent]  = modimg[goodent]
        #
        # ax.imshow(imgarr, vmin=-0.01, vmax=0.05, cmap='viridis', origin='lower')
        # plt.show()

    else:
        print(' - Going directly to analysis of simulated data (assuming they exist)')

    print(' - Analyze simulated data ')
    print('   Now reloading simulations for fitting')
    grp = grizli.multifit.GroupFLT(grism_files=glob.glob('niriss_*flt.fits'), direct_files=[],
                                   ref_file=imagepath+'hlsp_xdf_hst_wfc3ir-60mas_hudf_f140w_v1_sci.fits',
                                   seg_file='udf_f140w_photutils_seg.fits',
                                   catalog='udf_f140w_photutils.cat',
                                   cpu_count=0, pad=200)



    print(' - First pass contamination model, flat continuum')
    grp.compute_full_model(mag_limit=26)

    print(' - Refine the (polynomial) continuum model for brighter objects')
    grp.refine_list(poly_order=2, mag_limits=[16, 24], verbose=False)

    # print(' - saving "grp" data')
    # grp.save_full_data()

    print(' - Fit parameters')
    pzfit, pspec2, pline = grizli.multifit.get_redshift_fit_defaults()

    # Redshift fit
    pzfit ['zr'] = [0.5, 2.4]
    pzfit['dz'] = [0.01, 0.001]

    # Drizzled line maps
    pline = {'kernel': 'square', 'pixfrac': 0.8, 'pixscale': 0.06, 'size': 10}

    # Full rectified 2D spectrum
    pspec2 = {'NY': 20, 'dlam': 50, 'spatial_scale': 1}

    ## emission line object
    id = 898 # ID in the detection catalog / segm image

    print(' - Extract spectrum cutouts from individual FLTs of object '+str(id))
    beams = grp.get_beams(id, size=40)

    print('   Put them in a MultiBeam object')
    mb = grizli.multifit.MultiBeam(beams, fcontam=1, group_name='niriss-udf')

    print('   Run the redshift fit and generate the emission line map')
    out = mb.run_full_diagnostics(pzfit=pzfit, pspec2=pspec2, pline=pline,
                                  GroupFLT=grp, prior=None, verbose=False)

    fit, fig, fig2, hdu2, hdu_line = out
    cmap = 'viridis_r'

    print(' - "Beams" are extracted for spectra of a given order.  Have attributes for contam, model etc.')
    fig = plt.figure(figsize=[9,9*1.2/3])
    for ix, i in enumerate([0,2,4,1,3,5]):
        ax = fig.add_subplot(2,3,ix+1)
        beam = mb.beams[i]
        ax.imshow(beam.grism['SCI'], vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower',
                  aspect='auto')
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                transform=ax.transAxes, size=10, ha='left', va='bottom')

    fig.axes[0].set_ylabel('Observed\nspectrum')
    fig.tight_layout(pad=0.1)
    plt.savefig('./extracted2Dregion.pdf')

    print(' - Each beam carries with it a static contamination model extracted from the full field')
    fig = plt.figure(figsize=[9,9*1.2/3])
    for ix, i in enumerate([0,2,4,1,3,5]):
        ax = fig.add_subplot(2,3,ix+1)
        beam = mb.beams[i]
        ax.imshow(beam.contam, vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                transform=ax.transAxes, size=10, ha='left', va='bottom')

    fig.axes[0].set_ylabel('Contamination model')
    fig.tight_layout(pad=0.1)
    plt.savefig('./contaminationModel.pdf')


    print(' - Under the hood, the fitting is done by specifying a single 1D template, which')
    print('   is used to generate model 2D spectra for each beam')
    fig = plt.figure(figsize=[9,9*1.2/3])
    for ix, i in enumerate([0,2,4,1,3,5]):
        ax = fig.add_subplot(2,3,ix+1)
        beam = mb.beams[i]
        ax.imshow(beam.model, vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                transform=ax.transAxes, size=10, ha='left', va='bottom')

    fig.axes[0].set_ylabel('Spectrum Model')
    fig.tight_layout(pad=0.1)
    plt.savefig('./observedSpectrum_model.pdf')

    print(' - Goodness of fit is computed by comparing the models in the full 2D pixel space')
    fig = plt.figure(figsize=[9,9*1.2/3])
    for ix, i in enumerate([0,2,4,1,3,5]):
        ax = fig.add_subplot(2,3,ix+1)
        beam = mb.beams[i]
        ax.imshow(beam.grism['SCI'] - beam.contam - beam.model, vmin=-0.01, vmax=0.05, cmap=cmap,
                  origin='lower', aspect='auto')
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                transform=ax.transAxes, size=10, ha='left', va='bottom')

    fig.axes[0].set_ylabel('Full residuals')
    fig.tight_layout(pad=0.1)
    plt.savefig('./fullresiduals.pdf')

    print(' - Trivial demo model, dropout in the middle of the F115W filter')
    xspec = np.arange(0.8,2.4,0.02)*1.e4
    yspec = (xspec/1.4e4)**(-2) # Beta=-2
    yspec[xspec < 1.2e4] = 0.
    plt.plot(xspec, yspec)
    mb.compute_model(spectrum_1d=[xspec, yspec])

    fig = plt.figure(figsize=[9,9*1.2/3])
    for ix, i in enumerate([0,2,4,1,3,5]):
        ax = fig.add_subplot(2,3,ix+1)
        beam = mb.beams[i]
        ax.imshow(beam.beam.model, vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                transform=ax.transAxes, size=10, ha='left', va='bottom')

    fig.axes[0].set_ylabel('Observed\nspectrum')
    fig.tight_layout(pad=0.1)
    plt.savefig('./simpeldropoutdemo.pdf')


    print(' - Emission line map')
    #line = fits.open('niriss-udf_zfit_00898.line.fits')
    line = fits.open('niriss-udf_00898.line.fits')
    print(line[0].header['HASLINES'])
    line.info()

    cmap = 'cubehelix_r'
    fig = plt.figure(figsize=[12,3])

    ax = fig.add_subplot(141)
    ax.imshow(line['DSCI'].data, vmin=-0.01, vmax=0.02, cmap=cmap, origin='lower')
    ax.text(5,5,'F140W direct image', ha='left', va='bottom')

    ax = fig.add_subplot(142)
    ax.imshow(line['LINE', 'Ha'].data, vmin=-0.01, vmax=0.02, cmap=cmap, origin='lower')
    ax.text(5,5,r'H$\alpha$', ha='left', va='bottom')

    ax = fig.add_subplot(143)
    ax.imshow(line['LINEWHT', 'Ha'].data, vmin=-0.01, vmax=20000, cmap='gray', origin='lower')
    ax.text(5,5,r'H$\alpha$, weight', ha='left', va='bottom', color='w')

    try:
        ax = fig.add_subplot(144)
        ax.imshow(line['LINE', 'OIII'].data, vmin=-0.03, vmax=0.06, cmap=cmap, origin='lower')
        ax.text(5,5,r'[OIII]$\lambda$4959,5007', ha='left', va='bottom')
    except:
        ax = fig.add_subplot(144)
        ax.imshow(line['LINE', 'Hd'].data, vmin=-0.03, vmax=0.06, cmap=cmap, origin='lower')
        ax.text(5,5,r'H$\delta$', ha='left', va='bottom')

    for ax in fig.axes:
        ax.set_xticklabels([]); ax.set_yticklabels([])

    fig.tight_layout(pad=0.1)

    plt.savefig('./emissionlinemap_Ha.pdf')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print('\n - Ran all commands successfully! ')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def NIRISSsim_A2744(generatesimulation=True):
    """

    NIRISS simulations of A2744 (based on NIRISSsim_UDF() above


    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw
    gw.NIRISSsim_A2744(generatesimulation=True)

    """
    cwd = os.getcwd()
    print('\n grizli version: %s' %(grizli.__version__))
    print('Simulating NIRISS exposure on A2744')

    print(' - Handle the GLASS object catalog and segmentation image ')
    print('   Define HST filter setup matching GLASS catalog')
    HFFimgfilter = 'f140w'

    imagepath  = '/Users/kschmidt/work/images_MAST/A2744/'
    ref_hffimg = imagepath+'hlsp_frontier_hst_wfc3-60mas_abell2744_'+HFFimgfilter+'_v1.0_drz.fits'

    if generatesimulation:
        print(' - Generating simulations')

        print('   Loading HFF '+HFFimgfilter+' images of A2744 from '+imagepath)

        # sci_hffimg = fits.open(ref_hffimg)
        # wht_hffimg = fits.open(ref_hffimg.replace('drz.fits','wht.fits'))
        #
        # rms = 1/np.sqrt(wht_hffimg[0].data)
        # rms[~np.isfinite(rms)] = 1000
        # dq = rms > 900
        # wcs = astropy.wcs.WCS(sci_hffimg[0].header)
        #
        # print(' - Run the detection with the grizli / photutils wrapper')
        # cat, seg = detect_with_photutils(sci_hffimg[0].data, err=rms, dq=dq, seg=None,
        #                                  detect_thresh=1.4, npixels=6, grow_seg=3,
        #                                  gauss_fwhm=2.0, gsize=3, wcs=wcs, save_detection=False,
        #                                  root='a2744_'+HFFimgfilter+'_photutils', background=None, gain=None,
        #                                  AB_zeropoint=AB_zeropoint,
        #                                  rename_columns={'ycentroid': 'y_flt', 'xcentroid': 'x_flt',
        #                                                  'dec_icrs_centroid': 'dec',
        #                                                  'ra_icrs_centroid': 'ra'},
        #                                  clobber=True, verbose=True)
        #
        # # code expects this columns later....
        # cat['NUMBER'] = cat['id']

        print(' - Loading GLASS A2744 source catalog and corresponding segmentation map')
        cat = np.genfromtxt('hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v001_glassmaster.txt',dtype=None,names=True)

        print(' - Get matches to ASTRODEEP catalog')
        ADpath  = '/Users/kschmidt/work/GLASS/LAEsearchFullGLASS/catalogs/ASTRODEEP/fullrelease/'
        AD_cat  = np.genfromtxt(ADpath+'A2744cl_26012016/A2744cl_A.cat',dtype=None,names=True)
        AD_cat  = AD_cat[AD_cat['MAG_JH140'] < 99.0]

        AD_radec  = SkyCoord(ra=AD_cat['RA']*u.degree, dec=AD_cat['DEC']*u.degree)
        cat_radec = SkyCoord(ra=cat['X_WORLD']*u.degree, dec=cat['Y_WORLD']*u.degree)

        matchtol = 2.0
        print('   Getting sources within the match toleracnce of 2 arc seconds')
        AD_idx, d2d, d3d = cat_radec.match_to_catalog_sky(AD_radec)
        has_AD_match = np.where(d2d < matchtol*u.arcsec)[0]
        Nmatches     = len(has_AD_match)
        print('   Found '+str(Nmatches)+' objects in the GLASS catalog with matches to the ASTRODEEP photometry')

        print('   Writing catalog and segmentation map to files')
        # Use ASTRODEEP mags because GLASS mags are sub-optimal
        AD_mag = AD_cat['MAG_JH140']
        cat['MAG_AUTO'] = AD_mag[AD_idx]
        np.savetxt('a2744_f140w_glassmodified.cat',cat,header=' '.join(cat.dtype.names))

        # cat.write('a2744_f140w_photutils.cat', format='ascii.commented_header')
        # fits.writeto('a2744_f140w_photutils_seg.fits', data=np.cast[int](seg),
        #              header=sci_hffimg[0].header, clobber=True)

        print(' - Setup fake (noise) images, centered in the UDF/XDF')
        ra, dec = 3.588197688, -30.39784202
        pa_aper = 135 # Using GLASS PA_V3

        np.random.seed(1)
        # Crude exposure parameters.  Rough read noise and backgrounds coded in
        # grizli.fake_image to make some reasonable noise estimate
        EXPTIME = 5218.07 # seconds exposure in GLASS-ERS grism setup per PA
        NEXP    = 6       # total Intergration (6 total dithers per setup, i.e. filter+C/R grism)

        # JWST NIRISS, three filters & two orients
        for filt in ['F115W', 'F150W', 'F200W']:
            for theta in [0,90]:
                h, wcs = grizli.fake_image.niriss_header(filter=filt, ra=ra, dec=dec,
                                                         pa_aper=pa_aper+theta)
                print('Filter: {filter}, Background: {bg} e/s/pix, RN: {RN} e/exp'.format(filter=filt,
                                                                bg=h['BACKGR'], RN=h['READN']))
                output = 'niriss_{filt}_{theta:02d}_flt.fits'.format(filt=filt, theta=theta)
                grizli.fake_image.make_fake_image(h, output=output, exptime=EXPTIME, nexp=NEXP)


        print(' - Load GroupFLT for simulation, NB: input files are just noise')
        sim = grizli.multifit.GroupFLT(grism_files=glob.glob('niriss_*flt.fits'), direct_files=[],
                                       ref_file=ref_hffimg, ref_ext=0,
                                       seg_file='hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v001_align-drz-seg.fits',
                                       catalog='a2744_f140w_glassmodified.cat',
                                       cpu_count=0, # <0 dont parallelize; =0 use all available; >0 CPUs to use
                                       pad=200)

        print(' - First pass, flat continuum ')
        sim.compute_full_model(mag_limit=27)

        print(' - Compute model grism spectra for ASTRODEEP matches based on full EAZY photo-z templates')
        detection_bp = pysynphot.ObsBandpass('wfc3,ir,'+HFFimgfilter)
        for ix in has_AD_match:
            AD_id = AD_cat['ID'][AD_idx[ix]]
            id    = cat['NUMBER'][ix]

            infostr = '   Computing model for id '+str(id)+' / ASTRODEEP '+str(AD_id)+' ('+\
                      str("%6.f" % (ix+1))+' / '+str("%6.f" % Nmatches)+')          '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

            temp_sed_dir = '/Users/kschmidt/work/GLASS/LAEsearchFullGLASS/EAZYruns/180105/' \
                           'A2744_180105/eazy_output_A2744_180105/'


            EAZYtempsed = glob.glob(temp_sed_dir+'eazy_output_*_'+str(AD_id)+'.temp_sed')
            if len(EAZYtempsed) == 1:
                temp_sed_dat = np.genfromtxt(EAZYtempsed[0],dtype=None,names=True,skip_header=0)
            else:
                print(' Did not just find 1 EAZY template SED for ASTRODEEP object '+str(AD_id))
                pdb.set_trace()

            temp_lambda = temp_sed_dat['lambda_zprior']
            temp_flux   = temp_sed_dat['tempflux_zprior']

            # Needs to be normalized to unity in the detection band
            spec = pysynphot.ArraySpectrum(wave=temp_lambda, flux=temp_flux, waveunits='angstroms', fluxunits='flam')
            spec = spec.renorm(1., 'flam', detection_bp)


            #print(id)
            sim.compute_single_model(id, mag=cat['MAG_AUTO'][ix], size=-1, store=False,
                                     spectrum_1d=[spec.wave, spec.flux], get_beams=None,
                                     in_place=True)

        print(' - Plot blotted reference image; "grism" exposures are still empty, just noise')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate([0,2,4,1,3,5]):
            ax = fig.add_subplot(2,3,ix+1)
            ax.imshow(sim.FLTs[i].direct['REF'], vmin=-0.01, vmax=0.05, cmap='viridis',
                      origin='lower')
            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.grid(color='w', alpha=0.8)
            ax.text(100,100,sim.FLTs[i].grism_file, color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./blotted_reference_image.pdf')

        print(' - Plot blotted segmentation image; "grism" exposures are still empty, just noise')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate([0,2,4,1,3,5]):
            ax = fig.add_subplot(2,3,ix+1)
            ax.imshow(sim.FLTs[i].seg, vmin=-0.01, vmax=3000, cmap='viridis', origin='lower')
            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.grid(color='w', alpha=0.8)
            ax.text(100,100,sim.FLTs[i].grism_file, color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./blotted_segmentation_image.pdf')

        print(' - Plot "grism" exposures that are still empty, just noise')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate([0,2,4,1,3,5]):
            ax = fig.add_subplot(2,3,ix+1)
            ax.imshow(sim.FLTs[i].grism['SCI'], vmin=-0.01, vmax=0.05, cmap='viridis', origin='lower')
            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.grid(color='w', alpha=0.8)
            ax.text(100,100,sim.FLTs[i].grism_file, color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./blotted_grismnoise_image.pdf')

        print(' - Plot model stored in FLTs[i].model attribute')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate([0,2,4,1,3,5]):
            ax = fig.add_subplot(2,3,ix+1)

            # show as if it were the rotated grism
            if (i % 2) > 0:
                ax.imshow(np.rot90(sim.FLTs[i].model,-1), vmin=-0.01, vmax=0.05, cmap='viridis',
                          origin='lower')
            else:
                ax.imshow(sim.FLTs[i].model, vmin=-0.01, vmax=0.05, cmap='viridis', origin='lower')

            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.grid(color='w', alpha=0.8)
            ax.text(100,100,sim.FLTs[i].grism_file, color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./flt_model_attribute.pdf')

        print(' - Update SCI extension of the fake FLT images with the models just computed')
        for flt in sim.FLTs:
            print('Update', flt.grism_file)
            orig_flt = fits.open(flt.grism_file, mode='update')
            orig_flt['SCI'].data += flt.model[flt.pad:-flt.pad, flt.pad:-flt.pad]
            orig_flt.flush()

        # to make sure that only sensible flux values are plotted us the following:
        # fig = plt.figure(figsize=[5,5])
        # ax = fig.add_subplot(1,1,1)
        #
        # #modimg  = sim.FLTs[i].direct['REF']
        # modimg  = sim.FLTs[i].model
        # goodent = (np.abs(modimg) < 1000.0)
        # imgarr  = modimg*0.0
        # imgarr[goodent]  = modimg[goodent]
        #
        # ax.imshow(imgarr, vmin=-0.01, vmax=0.05, cmap='viridis', origin='lower')
        # plt.show()

    else:
        print(' - Going directly to analysis of simulated data (assuming they exist)')

    print('   Now reloading simulations for fitting')
    sim = grizli.multifit.GroupFLT(grism_files=glob.glob('niriss_*flt.fits'), direct_files=[],
                                   ref_file=ref_hffimg, ref_ext=0,
                                   seg_file='hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v001_align-drz-seg.fits',
                                   catalog='a2744_f140w_glassmodified.cat',
                                   cpu_count=0, # <0 dont parallelize; =0 use all available; >0 CPUs to use
                                   pad=200)

    grp = sim

    print('\n - Analyze simulated data ')
    print(' - First pass contamination model, flat continuum')
    grp.compute_full_model(mag_limit=26)

    print(' - Refine the (polynomial) continuum model for brighter objects')
    grp.refine_list(poly_order=2, mag_limits=[16, 24], verbose=False)

    print(' - saving refined simulation output')
    grp.save_full_data()

    print(' - Fit parameters')
    pzfit, pspec2, pline = grizli.multifit.get_redshift_fit_defaults()

    # Redshift fit
    pzfit ['zr'] = [0.5, 2.4]
    pzfit['dz'] = [0.01, 0.001]

    # Drizzled line maps
    pline = {'kernel': 'square', 'pixfrac': 0.8, 'pixscale': 0.06, 'size': 10}

    # Full rectified 2D spectrum
    pspec2 = {'NY': 20, 'dlam': 50, 'spatial_scale': 1}

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## emission line object
    extractid = 793 # GLASS Line emitter at 1.34000 with zQ = 4.00000; ra,dec = 3.6061,-30.3920
    basename  = 'niriss-a2744_'+str("%.5d" % extractid)+'_'

    print(' - Extract spectrum cutouts from individual FLTs of object '+str(extractid))
    beams = grp.get_beams(extractid, size=40)

    print('   Put them in a MultiBeam object')
    mb = grizli.multifit.MultiBeam(beams, fcontam=1, group_name='niriss-a2744')

    mb.write_master_fits(get_hdu=False) # Saving fits file with all beams for "extractid"

    print('   Run the redshift fit and generate the emission line map')
    out = mb.run_full_diagnostics(pzfit=pzfit, pspec2=pspec2, pline=pline,
                                  GroupFLT=grp, prior=None, verbose=False)

    fit, fig, fig2, hdu2, hdu_line = out
    cmap = 'viridis_r'

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' - "Beams" are extracted for spectra of a given order.  Have attributes for contam, model etc.')
    fig = plt.figure(figsize=[9,9*1.2/3])
    for ix, i in enumerate([0,2,4,1,3,5]):
        ax = fig.add_subplot(2,3,ix+1)
        beam = mb.beams[i]
        ax.imshow(beam.grism['SCI'], vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                transform=ax.transAxes, size=10, ha='left', va='bottom')

    fig.axes[0].set_ylabel('Extraction')
    fig.tight_layout(pad=0.1)
    plt.savefig('./'+basename+'extracted2Dregion.pdf')
    plt.clf()
    plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' - Each beam carries with it a static contamination model extracted from the full field')
    fig = plt.figure(figsize=[9,9*1.2/3])
    for ix, i in enumerate([0,2,4,1,3,5]):
        ax = fig.add_subplot(2,3,ix+1)
        beam = mb.beams[i]
        ax.imshow(beam.contam, vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                transform=ax.transAxes, size=10, ha='left', va='bottom')

    fig.axes[0].set_ylabel('Contamination')
    fig.tight_layout(pad=0.1)
    plt.savefig('./'+basename+'contaminationModel.pdf')
    plt.clf()
    plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' - Under the hood, the fitting is done by specifying a single 1D template, which')
    print('   is used to generate model 2D spectra for each beam')
    fig = plt.figure(figsize=[9,9*1.2/3])
    for ix, i in enumerate([0,2,4,1,3,5]):
        ax = fig.add_subplot(2,3,ix+1)
        beam = mb.beams[i]
        ax.imshow(beam.model, vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                transform=ax.transAxes, size=10, ha='left', va='bottom')

    fig.axes[0].set_ylabel('Spec. Model')
    fig.tight_layout(pad=0.1)
    plt.savefig('./'+basename+'observedSpectrum_model.pdf')
    plt.clf()
    plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' - Goodness of fit is computed by comparing the models in the full 2D pixel space')
    fig = plt.figure(figsize=[9,9*1.2/3])
    for ix, i in enumerate([0,2,4,1,3,5]):
        ax = fig.add_subplot(2,3,ix+1)
        beam = mb.beams[i]
        ax.imshow(beam.grism['SCI'] - beam.contam - beam.model, vmin=-0.01, vmax=0.05, cmap=cmap,
                  origin='lower', aspect='auto')
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                transform=ax.transAxes, size=10, ha='left', va='bottom')

    fig.axes[0].set_ylabel('Residuals')
    fig.tight_layout(pad=0.1)
    plt.savefig('./'+basename+'fullresiduals.pdf')
    plt.clf()
    plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' - Emission line map')
    #line = fits.open('niriss-a2744_zfit_00898.line.fits')
    line = fits.open('niriss-a2744_'+str("%.5d" % extractid)+'.line.fits')
    print(line[0].header['HASLINES'])
    line.info()

    cmap = 'cubehelix_r'
    fig = plt.figure(figsize=[12,3])

    ax = fig.add_subplot(141)
    ax.imshow(line['DSCI'].data, vmin=-0.01, vmax=0.02, cmap=cmap, origin='lower')
    ax.text(5,5,'F140W direct image', ha='left', va='bottom')

    ax = fig.add_subplot(142)
    ax.imshow(line['LINE', 'Ha'].data, vmin=-0.01, vmax=0.02, cmap=cmap, origin='lower')
    ax.text(5,5,r'H$\alpha$', ha='left', va='bottom')

    ax = fig.add_subplot(143)
    ax.imshow(line['LINEWHT', 'Ha'].data, vmin=-0.01, vmax=20000, cmap='gray', origin='lower')
    ax.text(5,5,r'H$\alpha$, weight', ha='left', va='bottom', color='w')

    try:
        ax = fig.add_subplot(144)
        ax.imshow(line['LINE', 'OIII'].data, vmin=-0.03, vmax=0.06, cmap=cmap, origin='lower')
        ax.text(5,5,r'[OIII]$\lambda$4959,5007', ha='left', va='bottom')
    except:
        pass

    try:
        ax = fig.add_subplot(144)
        ax.imshow(line['LINE', 'Hd'].data, vmin=-0.03, vmax=0.06, cmap=cmap, origin='lower')
        ax.text(5,5,r'H$\delta$', ha='left', va='bottom')
    except:
        pass

    for ax in fig.axes:
        ax.set_xticklabels([]); ax.set_yticklabels([])

    fig.tight_layout(pad=0.1)

    plt.savefig('./'+basename+'emissionlinemap_Ha.pdf')
    plt.clf()
    plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print('\n - Ran all commands successfully! ')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def NIRCAMsim_A2744(generatesimulation=True, runfulldiagnostics=True, mockspec_type='manual_lines',
                    singlefilterrun=False, quickrun=False, plotsingleobj=False,
                    extractids=[65,476,726,233,754,755,793,848]):
    """

    NIRCAM simulations of A2744 (based on grizli_wrappers.NIRISSsim_A2744()

    --- INPUT ---
    generatesimulation    Generate the actual simulations (to save time for analysis and plotting,
                          only run this once to generate files etc.)
    runfulldiagnostics    Run the full redshift fits for objects to extract? Otherwise, only spectra extracted.
    mockspec_type         What mock spectra to use in the simulation. Choose between:
                             eazy_templates          The SED templates from the SED fits
                             manual_lines            Manually add lines to flat continuum and assign JADES spectra to
                                                     a few selected objects.
                             jades_templatematches   Assign the best-fitting JADES templates (in terms of mF140W and
                                                     redshift) to each of the objects reading matches from catalog
    singlefilterrun       Only run simulation for F356W to speed up things
    quickrun              If True corners are cut to make the modeling proceed faster (for testing)
    plotsingleobj         Plot the beams of the extracted objects?
    extractids            List of GLASS ids to extract.

    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw
    gw.NIRCAMsim_A2744(generatesimulation=True, runfulldiagnostics=False, mockspec_type='manual_lines')

    objects2extract = [476]

    gw.NIRCAMsim_A2744(generatesimulation=False, mockspec_type='jades_templatematches',singlefilterrun=False, quickrun=False, runfulldiagnostics=False,plotsingleobj=False,extractids=objects2extract)

    """
    Ncpu      = 1        # <0 dont parallelize; =0 use all available; >0 CPUs to use
    padval    = 2000     #

    # determine polynomial coeffficients for flat continuum model
    polymodelcoeffs = [0.1, 0.05] # coeffs=[1.2, -0.5]))
    polyxrange      = [2.0, 6.0]

    cwd = os.getcwd()+'/'
    if 'NIRCAM' not in cwd:
        sys.exit('"NIRCAM" is not part of workign directory name. Assumes this means the working location is incorrect')
    else:
        print('')
        print('\n --- Simulating NIRCAM grism obs of A2744 in: ')
        print(' --- '+cwd)
        print(' --- Started wrapper at:                       '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print('\n - Using grizli version: %s' %(grizli.__version__))

    if singlefilterrun:
        print(' \n\n          NB - running a "singlefilterrun" \n\n')
        filterloops  = ['F277W'] # ['F356W']
        orientations = [0,90]

    if quickrun:
        print(' \n\n          NB - using the "quickrun" setup \n\n')
        matchtol  = 0.01      # arcsec
        mag_limit  = 20
        mag_limits = [16, 19]
    else:
        matchtol  = 2.0      # arcsec
        mag_limit  = 26.0     # Magnitude limit for traces to simulate and account for
        mag_limits = [10, 26] # Magnitude range for bright objects to refine model for

    print(' - Handle the GLASS object catalog and segmentation image ')
    print('   Define HST filter setup matching GLASS catalog')
    HFFimgfilter    = 'f140w'
    imagepath       = '/Users/kschmidt/work/images_MAST/A2744/'
    ref_hffimg      = imagepath+'hlsp_frontier_hst_wfc3-60mas_abell2744_'+HFFimgfilter+'_v1.0_drz.fits'
    photcatout      = cwd+'a2744_f140w_glasscat_with_ASTRODEEP_f140w_mag.cat'
    segmentationmap = cwd+'hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v001_align-drz-seg.fits'

    if generatesimulation:
        ra, dec         = 3.588197688, -30.39784202 # Cluster center

        ra277R, dec277R = 3.574423698, -30.40970039 # SW ~1 arcmin from (ra, dec)

        #ra277C, dec277C = 3.601881571, -30.40965939 # Offset for 277 dispersion in column direction; ~1 arcmin from (ra, dec)
        #ra277C, dec277C = 3.574515015, -30.38604164 # Offset for 277 dispersion in column direction; ~1 arcmin from (ra, dec)
        ra277C, dec277C = 3.588197688, -30.39784202 # Cluster center
        #ra277C, dec277C = 3.602111525, -30.38638374 # Offset for 277 dispersion in row direction (180); ~1 arcmin from (ra, dec)
        #ra277C, dec277C = 3.594818398, -30.39189161 # NE ~30 arcsec from (ra, dec)
        #ra277C, dec277C = 3.595044930, -30.40375511 # SE ~30 arcsec from (ra, dec)
        #ra277C, dec277C = 3.580995897, -30.3917931  # NW ~30 arcsec from (ra, dec)

        pa_aper = 135      # Using GLASS PA_V3
        EXPTIME = 6012.591 # seconds exposure (GLASS NIRISS ERS is 5218.07)
        NEXP    = 16       # total Intergrations

        print(' - Generating simulations')
        print(' - Will loda HFF '+HFFimgfilter+' images of A2744 from \n   '+imagepath)
        print(' - Loading GLASS A2744 source catalog and corresponding segmentation map')
        cat = np.genfromtxt('hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v001_glassmaster.txt',dtype=None,names=True)

        print(' - Get matches to ASTRODEEP catalog')
        ADpath  = '/Users/kschmidt/work/GLASS/LAEsearchFullGLASS/catalogs/ASTRODEEP/fullrelease/'
        AD_cat  = np.genfromtxt(ADpath+'A2744cl_26012016/A2744cl_A.cat',dtype=None,names=True)

        print(' - Ignore objects with MAG_JH140 >= 99.0')
        AD_cat  = AD_cat[AD_cat['MAG_JH140'] < 99.0]

        AD_radec  = SkyCoord(ra=AD_cat['RA']*u.degree, dec=AD_cat['DEC']*u.degree)
        cat_radec = SkyCoord(ra=cat['X_WORLD']*u.degree, dec=cat['Y_WORLD']*u.degree)

        print(' - Getting sources within the match toleracnce of '+str(matchtol)+' arc seconds')
        AD_idx, d2d, d3d = cat_radec.match_to_catalog_sky(AD_radec)

        print(' - Writing modified catalog to file')
        print('   Using ASTRODEEP magnitudes because GLASS magnitudes are sub-optimal')
        AD_mag = AD_cat['MAG_JH140']
        cat['MAG_AUTO'] = AD_mag[AD_idx]
        np.savetxt(photcatout,cat,header=' '.join(cat.dtype.names))
        print('   Stored modified GLASS photometric catalog to\n   '+photcatout)

        has_AD_match = np.where(d2d < matchtol*u.arcsec)[0]
        # if quickrun:
        #     AD_idx, d2d, d3d = AD_idx[has_AD_match][10:40], d2d[has_AD_match][10:40], d3d[has_AD_match][10:40]
        #     has_AD_match     = has_AD_match[10:40]

        Nmatches     = len(has_AD_match)
        print('   Found '+str(Nmatches)+' objects in the GLASS catalog with matches to the ASTRODEEP photometry')

        print(' - Setup fake (noise) images, centered in the UDF/XDF')
        np.random.seed(1)
        print('   Read noise and background in images scaled to NEXP = '+str(NEXP)+' and EXPTIME = '+str(EXPTIME)+' seconds')
        # JWST NIRCAM, three filters & two orients
        if singlefilterrun:
            filterloops  = filterloops
            orientations = orientations
        else:
            filterloops  = ['F277W', 'F356W', 'F444W']
            orientations = [0,90]

        for filt in filterloops:
            for theta in orientations:
                if theta == 0:
                    grism_img = 'GRISMR'
                elif theta == 90:
                    grism_img = 'GRISMC'

                if (filt == 'F277W') & (theta == 0): # offset f277w row disperion
                    pointra  = ra277R
                    pointdec = dec277R
                elif (filt == 'F277W') & (theta == 90): # offset f277w column disperion
                    pointra  = ra277C
                    pointdec = dec277C
                else:
                    pointra  = ra
                    pointdec = dec

                h, wcs = grizli.fake_image.nircam_header(filter=filt, ra=pointra, dec=pointdec,
                                                         pa_aper=pa_aper+theta,grism=grism_img)
                print('Filter: {filter}, Background: {bg} e/s/pix, RN: {RN} e/exp'.format(filter=filt,
                                                                bg=h['BACKGR'], RN=h['READN']))
                output = 'nircam_{filt}_{theta:02d}_flt.fits'.format(filt=filt, theta=theta)
                grizli.fake_image.make_fake_image(h, output=output, exptime=EXPTIME, nexp=NEXP)


        print(' - Load GroupFLT for simulation, NB: input files still just noise')
        sim = grizli.multifit.GroupFLT(grism_files=glob.glob('nircam_*flt.fits'), direct_files=[],
                                       ref_file=ref_hffimg, ref_ext=0,
                                       seg_file=segmentationmap,
                                       catalog=photcatout,
                                       cpu_count=Ncpu, # <0 dont parallelize; =0 use all available; >0 CPUs to use
                                       pad=padval,
                                       polyx=polyxrange)

        print(' - Compute full model (for simulation); started at:                    '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) # First pass, flat continuum

        sim.compute_full_model(mag_limit=mag_limit, coeffs=polymodelcoeffs)

        print(' - Run compute_single_model with templates')
        temp_sed_dir = '/Users/kschmidt/work/GLASS/LAEsearchFullGLASS/EAZYruns/180105/' \
                       'A2744_180105/eazy_output_A2744_180105/'
        lamcol       = 'lambda_zprior'
        fluxcol      = 'tempflux_zprior'
        print('   Using wavelength and flux columns '+lamcol+' and '+fluxcol)

        detection_bp = pysynphot.ObsBandpass('wfc3,ir,'+HFFimgfilter)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        print('\n - Buildinfg and/or assembling models for individual objects in simulation ')
        if mockspec_type.lower() == 'eazy_templates':
            gw.compute_single_model_EAZY(cat,Nmatches,temp_sed_dir,sim,AD_cat,AD_idx,has_AD_match,lamcol,fluxcol,detection_bp)
        elif mockspec_type == 'manual_lines':
            gw.compute_single_model_MANUAL(sim,detection_bp,has_AD_match,AD_cat,AD_idx,cat,Nmatches)
        elif mockspec_type == 'jades_templatematches':
            if quickrun:
                quickrunIDs = AD_cat['ID'][AD_idx][has_AD_match]
            else:
                quickrunIDs = False
            gw.compute_single_model_JADES(sim,cat,detection_bp,quickrunIDs=quickrunIDs)
        print('   Done buildinfg and/or assembling models for individual objects in simulation \n')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        subplot_indices = [0,2,4,1,3,5]
        if singlefilterrun:
            subplot_indices = [0]
            if len(orientations) == 2:
                subplot_indices = [0,1]
        print(' - Plot blotted reference image; "grism" exposures are still empty, just noise')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate(subplot_indices):
            ax = fig.add_subplot(2,3,ix+1)
            ax.imshow(sim.FLTs[i].direct['REF'], vmin=-0.01, vmax=0.05, cmap='viridis',
                      origin='lower')
            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.grid(color='w', alpha=0.8)
            ax.text(100,100,sim.FLTs[i].grism_file, color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./blotted_reference_image.pdf')

        print(' - Plot blotted segmentation image; "grism" exposures are still empty, just noise')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate(subplot_indices):
            ax = fig.add_subplot(2,3,ix+1)
            ax.imshow(sim.FLTs[i].seg, vmin=-0.01, vmax=3000, cmap='viridis', origin='lower')
            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.grid(color='w', alpha=0.8)
            ax.text(100,100,sim.FLTs[i].grism_file, color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./blotted_segmentation_image.pdf')

        print(' - Plot "grism" exposures that are still empty, just noise')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate(subplot_indices):
            ax = fig.add_subplot(2,3,ix+1)
            ax.imshow(sim.FLTs[i].grism['SCI'], vmin=-0.01, vmax=0.05, cmap='viridis', origin='lower')
            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.grid(color='w', alpha=0.8)
            ax.text(100,100,sim.FLTs[i].grism_file, color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./blotted_grismnoise_image.pdf')

        print(' - Plot model stored in FLTs[i].model attribute')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate(subplot_indices):
            ax = fig.add_subplot(2,3,ix+1)

            # show as if it were the rotated grism
            if (i % 2) > 0:
                ax.imshow(np.rot90(sim.FLTs[i].model,-1), vmin=-0.01, vmax=0.05, cmap='viridis',
                          origin='lower')
            else:
                ax.imshow(sim.FLTs[i].model, vmin=-0.01, vmax=0.05, cmap='viridis', origin='lower')

            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.grid(color='w', alpha=0.8)
            ax.text(100,100,sim.FLTs[i].grism_file, color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./flt_model_attribute.pdf')

        print(' - Update SCI extension of the fake FLT images with the models just computed')
        for flt in sim.FLTs:
            print('Update', flt.grism_file)
            orig_flt = fits.open(flt.grism_file, mode='update')
            orig_flt['SCI'].data += flt.model[flt.pad:-flt.pad, flt.pad:-flt.pad]
            orig_flt.flush()

    else:
        print('\n - NB: Going directly to analysis of simulated data (assuming they exist)')

    print('   Reloading simulations to update the SCI extension')
    grp = grizli.multifit.GroupFLT(grism_files=glob.glob('nircam_*flt.fits'), direct_files=[],
                                   ref_file=ref_hffimg, ref_ext=0,
                                   seg_file=segmentationmap,
                                   catalog=photcatout,
                                   cpu_count=Ncpu, # <0 dont parallelize; =0 use all available; >0 CPUs to use
                                   pad=padval,
                                   polyx=polyxrange) # range the polynomial model is calculated over in microns

    print('\n - Computing model (for contam model); started at:                    '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print('   First pass to get contamination model, flat continuum')
    grp.compute_full_model(mag_limit=mag_limit, store=False, cpu_count=Ncpu, coeffs=polymodelcoeffs)
    # print('   Refine the (polynomial) continuum model for brighter objects')
    # grp.refine_list(poly_order=len(polymodelcoeffs), mag_limits=mag_limits, verbose=False)

    # print(' - saving refined contamination model')
    grp.save_full_data()

    print('\n - Analyze simulated data; started at:             '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print(' - Get default fit parameters')
    pzfit, pspec2, pline = grizli.multifit.get_redshift_fit_defaults()

    # Redshift fit
    pzfit ['zr'] = [0.5, 2.4]
    pzfit['dz'] = [0.01, 0.001]

    # Drizzled line maps
    pline = {'kernel': 'square', 'pixfrac': 0.8, 'pixscale': 0.063, 'size': 10}

    # Full rectified 2D spectrum
    pspec2 = {'NY': 20, 'dlam': 50, 'spatial_scale': 1}

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## emission line object
    id848 = 848 # Bright central galaxy:  ID_00848 (RA,DEC) = (3.5877,-30.3964) GLASS redshift = 0.316 | redshift quality = 2.0

    # ------ Vulcani Ha map (GLASS paper VII+VIII) ------
    id065 = 65  # Pa-a in F277 ID_00065 (RA,DEC) = (3.5770,-30.3795) GLASS redshift = 0.496 | redshift quality = 4.0
    id476 = 476 # BV Ha map; Pa-a in F277 ID_00476 (RA,DEC) = (3.6044,-30.3850) GLASS redshift = 0.300 | redshift quality = 4.0
    id726 = 726 # BV Ha map; Pa-a in F277 ID_00726 (RA,DEC) = (3.5809,-30.3908) GLASS redshift = 0.292 | redshift quality = 4.0

    # A274400065      3.57700058 -30.37947795  0.496   20.3983
    # A274400380      3.59327307 -30.38437748  0.2965  19.1761
    # A274400476      3.60440483 -30.38495383  0.3     21.152
    # A274400726      3.58094727 -30.39080141  0.2915  19.6931
    # A274400983      3.57986591 -30.39523146  0.45    24.0103
    # A274401110      3.58550641 -30.39715509  0.31    21.072
    # A274401222      3.59032801 -30.40038939  0.498   19.3537
    # A274401528      3.57572054 -30.40537791  0.6535  22.1135
    # A274401669      3.59866354 -30.40493392  0.5     24.5403

    id233 = 233 # MUSE LAE                ID_00233 (RA,DEC) = (3.5778434,-30.3812148) MUSE ID = 14518; MUSE redshift = 3.3756
    id754 = 754 # MUSE LAE 'arc'          ID_00754 (RA,DEC) = (3.5887928,-30.3938037) MUSE IDs = 8826 & 8789; MUSE redshift = 3.9803

    id755 = 755 # Line emitter:           ID_00755 (RA,DEC) = (3.5874,-30.3933) GLASS redshift = 0.600 | redshift quality = 3.0
    id793 = 793 # Line emitter:           ID_00793 (RA,DEC) = (3.6061,-30.3920) GLASS redshift = 1.340 | redshift quality = 4.0

    for ee, extractid in enumerate(extractids):
        print('\n - Extract object '+str("%.5d" % extractid)+' (object '+str(ee+1)+'/'+str(len(extractids))+')     ')
        group_name = 'nircam-a2744'
        base_name  = group_name+'_'+str("%.5d" % extractid)

        print(' - Get spectrum cutouts from individual FLTs ')
        beams = grp.get_beams(extractid, size=30) # , center_rd=radecs[ee]

        if len(beams) == 0:
            print('   WARNING: No beams for object '+str(extractid)+' so moving on ')
            continue

        print('   Put them in a MultiBeam object')
        mb = grizli.multifit.MultiBeam(beams, fcontam=1, group_name=group_name)

        mb.write_master_fits(get_hdu=False, strip=True) # Saving fits file with all beams for "extractid"
        beamfile  = './'+base_name+'.beams.fits'
        gw.gen_sci_nocontam_beams_file(beamfile,overwrite=True)

        print(' - Run the redshift fit and generate the emission line map?')
        if runfulldiagnostics:

            print('   -> Yes')
            out = mb.run_full_diagnostics(pzfit=pzfit, pspec2=pspec2, pline=pline,
                                          GroupFLT=grp, prior=None, verbose=True)

            fit, fig, fig2, hdu2, hdu_line = out
        else:
            print('   -> No')

        if plotsingleobj:
            cmap = 'viridis_r'
            try:
                subplot_indices = [0,2,4,1,3,5]
                mb.beams[5]
                print('   Checking for beams... six beams exist!')
            except:
                subplot_indices = [0,2,2,1,3,3]
                print('   Checking for beams... did not find six beams')

            if singlefilterrun:
                subplot_indices = [0,0,0,0,0,0]

            Ngrism          = len(subplot_indices)/2
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            print(' - "Beams" are extracted for spectra of a given order.  Have attributes for contam, model etc.')
            fig = plt.figure(figsize=[9,9*1.2/Ngrism])
            for ix, i in enumerate(subplot_indices):
                ax = fig.add_subplot(2,Ngrism,ix+1)
                beam = mb.beams[i]
                ax.imshow(beam.grism['SCI'], vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')
                ax.set_xticklabels([]); ax.set_yticklabels([])
                ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                        transform=ax.transAxes, size=10, ha='left', va='bottom')

            fig.axes[0].set_ylabel('Extraction')
            fig.tight_layout(pad=0.1)
            plt.savefig('./objectsextracted_plots/'+base_name+'_extracted2Dregion.pdf')
            plt.clf()
            plt.close('all')

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            print(' - Each beam carries with it a static contamination model extracted from the full field')
            fig = plt.figure(figsize=[9,9*1.2/Ngrism])
            for ix, i in enumerate(subplot_indices):
                ax = fig.add_subplot(2,Ngrism,ix+1)
                beam = mb.beams[i]
                ax.imshow(beam.contam, vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')
                ax.set_xticklabels([]); ax.set_yticklabels([])
                ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                        transform=ax.transAxes, size=10, ha='left', va='bottom')

            fig.axes[0].set_ylabel('Contamination')
            fig.tight_layout(pad=0.1)
            plt.savefig('./objectsextracted_plots/'+base_name+'_contaminationModel.pdf')
            plt.clf()
            plt.close('all')

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            print(' - "Beams" are extracted for spectra of a given order.  Have attributes for contam, model etc.')
            fig = plt.figure(figsize=[9,9*1.2/Ngrism])
            for ix, i in enumerate(subplot_indices):
                ax = fig.add_subplot(2,Ngrism,ix+1)
                beam = mb.beams[i]
                ax.imshow(beam.grism['SCI']-beam.contam, vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')
                ax.set_xticklabels([]); ax.set_yticklabels([])
                ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                        transform=ax.transAxes, size=10, ha='left', va='bottom')

            fig.axes[0].set_ylabel('Extraction')
            fig.tight_layout(pad=0.1)
            plt.savefig('./objectsextracted_plots/'+base_name+'_extracted2Dregion_contamremove.pdf')
            plt.clf()
            plt.close('all')

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            print(' - Under the hood, the fitting is done by specifying a single 1D template, which')
            print('   is used to generate model 2D spectra for each beam')
            fig = plt.figure(figsize=[9,9*1.2/Ngrism])
            for ix, i in enumerate(subplot_indices):
                ax = fig.add_subplot(2,Ngrism,ix+1)
                beam = mb.beams[i]
                ax.imshow(beam.model, vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')
                ax.set_xticklabels([]); ax.set_yticklabels([])
                ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                        transform=ax.transAxes, size=10, ha='left', va='bottom')

            fig.axes[0].set_ylabel('Spec. Model')
            fig.tight_layout(pad=0.1)
            plt.savefig('./objectsextracted_plots/'+base_name+'_observedSpectrum_model.pdf')
            plt.clf()
            plt.close('all')

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            print(' - Goodness of fit is computed by comparing the models in the full 2D pixel space')
            fig = plt.figure(figsize=[9,9*1.2/Ngrism])
            for ix, i in enumerate(subplot_indices):
                ax = fig.add_subplot(2,Ngrism,ix+1)
                beam = mb.beams[i]
                ax.imshow(beam.grism['SCI'] - beam.contam - beam.model, vmin=-0.01, vmax=0.05, cmap=cmap,
                          origin='lower', aspect='auto')
                ax.set_xticklabels([]); ax.set_yticklabels([])
                ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
                        transform=ax.transAxes, size=10, ha='left', va='bottom')

            fig.axes[0].set_ylabel('Residuals')
            fig.tight_layout(pad=0.1)
            plt.savefig('./objectsextracted_plots/'+base_name+'_fullresiduals.pdf')
            plt.clf()
            plt.close('all')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        print(' - Plot emission line map?')
        if runfulldiagnostics:
            print('   -> Yes')
            line = fits.open(base_name+'.line.fits')
            print(line[0].header['HASLINES'])
            line.info()

            cmap = 'cubehelix_r'
            fig = plt.figure(figsize=[12,3])

            ax = fig.add_subplot(141)
            ax.imshow(line['DSCI'].data, vmin=-0.01, vmax=0.02, cmap=cmap, origin='lower')
            ax.text(5,5,'F140W direct image', ha='left', va='bottom')

            try:
                ax = fig.add_subplot(142)
                ax.imshow(line['LINE', 'Ha'].data, vmin=-0.01, vmax=0.02, cmap=cmap, origin='lower')
                ax.text(5,5,r'H$\alpha$', ha='left', va='bottom')

                ax = fig.add_subplot(143)
                ax.imshow(line['LINEWHT', 'Ha'].data, vmin=-0.01, vmax=20000, cmap='gray', origin='lower')
                ax.text(5,5,r'H$\alpha$, weight', ha='left', va='bottom', color='w')
            except:
                pass

            try:
                ax = fig.add_subplot(144)
                ax.imshow(line['LINE', 'OIII'].data, vmin=-0.03, vmax=0.06, cmap=cmap, origin='lower')
                ax.text(5,5,r'[OIII]$\lambda$4959,5007', ha='left', va='bottom')
            except:
                pass

            try:
                ax = fig.add_subplot(144)
                ax.imshow(line['LINE', 'Hd'].data, vmin=-0.03, vmax=0.06, cmap=cmap, origin='lower')
                ax.text(5,5,r'H$\delta$', ha='left', va='bottom')
            except:
                pass

            for ax in fig.axes:
                ax.set_xticklabels([]); ax.set_yticklabels([])

            fig.tight_layout(pad=0.1)

            plt.savefig('./objectsextracted_plots/'+base_name+'emissionlinemap_Ha.pdf')
            plt.clf()
            plt.close('all')
        else:
            print('   -> No')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        print('\n v v v v v v v v v v DS9 command to display beams  v v v v v v v v v v ')
        beamfile_nocontam = './'+base_name+'.beams_nocontam.fits'
        beamhdu           = fits.open(beamfile)
        Nsetup            = beamhdu[0].header['COUNT']

        ds9cmd = '\n ds9 '
        for setup in np.arange(Nsetup):
            ds9cmd = ds9cmd + \
                     beamfile+'['+str(1+setup*7)+'] '+ beamfile+'['+str(2+setup*7)+'] '+ \
                     beamfile+'['+str(6+setup*7)+'] '+ beamfile+'['+str(3+setup*7)+'] '+ \
                     beamfile+'['+str(7+setup*7)+'] '+ beamfile_nocontam+'['+str(3+setup*7)+'] '
        ds9cmd = ds9cmd + ' -lock frame image -tile grid layout 6 '+str(Nsetup)+' & '
        print(ds9cmd)

        ds9cmd = '\n ds9 '
        for setup in np.arange(Nsetup):
            ds9cmd = ds9cmd + beamfile_nocontam+'['+str(3+setup*7)+'] '
        ds9cmd = ds9cmd + ' -lock frame image -tile grid layout 1 '+str(Nsetup)+' & '
        print(ds9cmd)

        print(' ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ')
    print('\n - Ran all commands successfully at:               '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def compute_single_model_EAZY(cat,Nmatches,temp_sed_dir,sim,AD_cat,AD_idx,has_AD_match,lamcol,fluxcol,detection_bp):
    """

    """
    print(' - Compute model spectra for ASTRODEEP matches based on full EAZY photo-z templates')
    print('   In template SED directory '+temp_sed_dir)

    for ii, ix in enumerate(has_AD_match):
        AD_id = AD_cat['ID'][AD_idx[ix]]
        id    = cat['NUMBER'][ix]

        infostr = '   Computing model for id '+str(id)+' / ASTRODEEP '+str(AD_id)+' ('+\
                  str("%6.f" % (ii+1))+' / '+str("%6.f" % Nmatches)+')          '
        sys.stdout.write("%s\r" % infostr)
        sys.stdout.flush()

        EAZYtempsed = glob.glob(temp_sed_dir+'eazy_output_*_'+str(AD_id)+'.temp_sed')
        if len(EAZYtempsed) == 1:
            temp_sed_dat = np.genfromtxt(EAZYtempsed[0],dtype=None,names=True,skip_header=0)
        else:
            print(' Did not just find 1 EAZY template SED for ASTRODEEP object '+str(AD_id))
            pdb.set_trace()

        temp_lambda = temp_sed_dat[lamcol]
        temp_flux   = temp_sed_dat[fluxcol]

        # Needs to be normalized to unity in the detection band
        spec = pysynphot.ArraySpectrum(wave=temp_lambda, flux=temp_flux, waveunits='angstroms', fluxunits='flam')
        spec = spec.renorm(1., 'flam', detection_bp)

        #print(id)
        sim.compute_single_model(id, mag=cat['MAG_AUTO'][ix], size=-1, store=False,
                                 spectrum_1d=[spec.wave, spec.flux], get_beams=None,
                                 in_place=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def compute_single_model_MANUAL(sim,detection_bp,has_AD_match,AD_cat,AD_idx,cat,Nmatches):
    """

    """
    print(' - Compute models using idealized template for all objects')
    temp_lambda   = np.arange(1000,10000,0.5)

    nocontinuum   = False
    if nocontinuum:
        print(' \n\n          NB - Removing continuum for all (but the JADES) sources \n\n')
        temp_flux     = temp_lambda*0.0 + 1e-20
        normval       = 1e-20
    else:
        temp_flux     = temp_lambda*0.0 + 0.25
        normval       = 1.0

    waveplateaus  = [[4860,4864],[4958,4960],[5006,5008]] # FWHM ~ 4A, 2A, 2A
    plateaulevels = [6,3,9]
    for ww, wp in enumerate(waveplateaus):
        waveent            = np.where( (temp_lambda > wp[0]) & (temp_lambda < wp[1]))[0]
        temp_flux[waveent] = plateaulevels[ww]

    redshift       = 6.0
    temp_lambda_z  = temp_lambda * (1.0+redshift)

    # Needs to be normalized to unity in the detection band
    spec = pysynphot.ArraySpectrum(wave=temp_lambda_z, flux=temp_flux, waveunits='angstroms', fluxunits='flam')
    spec = spec.renorm(normval, 'flam', detection_bp)

    insertsinglemodels = True

    if insertsinglemodels:
        for ii, ix in enumerate(has_AD_match):
            AD_id = AD_cat['ID'][AD_idx[ix]]
            id    = cat['NUMBER'][ix]

            #if id > 1:
            if id in [233,754]:
                Jobj = 189021 # @ z = 3.3424 with HST F140W AB mag = 24.06
                Jobj = 232519 # @ z = 4.4651 with HST F140W AB mag = 23.85

                print('\n - Loading JADES spectrum for id '+str(id))
                JADESinfo, temp_lambda, temp_flux = ju.get_JADESspecAndInfo(Jobj,observedframe=True,verbose=True)
                print('\n')

                # Needs to be normalized to unity in the detection band
                spec = pysynphot.ArraySpectrum(wave=temp_lambda, flux=temp_flux, waveunits='angstroms', fluxunits='flam')
                spec = spec.renorm(1., 'flam', detection_bp)

            infostr = '   Set idealized model for id '+str(id)+' / ASTRODEEP '+str(AD_id)+' ('+\
                      str("%6.f" % (ii+1))+' / '+str("%6.f" % Nmatches)+')          '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

            #print(id)
            objmag = cat['MAG_AUTO'][ix]
            sim.compute_single_model(id, mag=objmag, size=-1, store=False,
                                     spectrum_1d=[spec.wave, spec.flux], get_beams=None,
                                     in_place=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def compute_single_model_JADES(sim,cat,detection_bp,quickrunIDs=False):
    """

    """
    JADESmatches = '/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/A2744_JADESmatches_180305.txt'
    print(' - Assigning JADES templates as models for individual objects using \n'+JADESmatches)

    dat = np.genfromtxt(JADESmatches,names=True,dtype=None,skip_header=10)
    # for ii, id in enumerate(dat['JADESid']):
    #     if (id > 0) & (dat['JADESz'][ii] < 0.203):
    #         print(str(id)+' '+str(dat['JADESz'][ii])+' '+str(dat['JADESf140wmag'][ii]))

    for ii, id in enumerate(cat['NUMBER']):
        if quickrunIDs is not False:
            if id not in quickrunIDs: continue
        catent = np.where(dat['id_GLASS'].astype(int) == id)[0]

        if len(catent) == 0:
            sys.exit(' Match to JADES catalog not performed for id_GLASS = '+str(id))

        Jobj = dat['JADESid'][ii]

        if Jobj == -99:
            Jobj   = 274533
            modstr = 'JADES template'+str("%9s" % Jobj)+' (Manual: z=6.255 & m140=25.45)'
        elif Jobj == -9999:
            Jobj   = 4
            modstr = 'JADES template'+str("%9s" % Jobj)+' (Manual: z=0.200 & m140=27.26)'
        else:
            modstr = 'JADES template'+str("%9s" % Jobj)+' (ID assigned in z/JADES match)'

        #print(' - Loading JADES spectrum for id_GLASS = '+str(id))
        JADESinfo, temp_lambda, temp_flux = ju.get_JADESspecAndInfo(Jobj,observedframe=True,verbose=False)

        # Needs to be normalized to unity in the detection band
        spec = pysynphot.ArraySpectrum(wave=temp_lambda, flux=temp_flux, waveunits='angstroms', fluxunits='flam')
        spec = spec.renorm(1., 'flam', detection_bp)

        infostr = '   Use '+modstr+' for id_GLASS '+str(id)+\
                  ' ('+str("%6.f" % (ii+1))+' / '+str("%6.f" % len(cat['NUMBER']))+')             '
        sys.stdout.write("%s\r" % infostr)
        sys.stdout.flush()

        objmag = dat['f140wmag_GLASS'][ii]
        sim.compute_single_model(id, mag=objmag, size=-1, store=False,
                                 spectrum_1d=[spec.wave, spec.flux], get_beams=None,
                                 in_place=True)
    print('\n   ... done')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_sci_nocontam_beams_file(beamfile,overwrite=False):
    """

    Generate a fits file containing SCI-CONTAM instead of just SCI

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw
    gw.gen_sci_nocontam_beams_file('/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/nircam-a2744_00793.beams.fits')

    """
    outname = beamfile.replace('.fits','_nocontam.fits')
    beamhdu = fits.open(beamfile)
    Nsetup  = beamhdu[0].header['COUNT']

    for ext in np.arange(Nsetup):
        beamhdu[3+ext*7].data              = beamhdu[3+ext*7].data - beamhdu[6+ext*7].data
        beamhdu[3+ext*7].header['EXTNAME'] = 'SCINOCONTAM'

    beamhdu.writeto(outname,overwrite=overwrite)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def determine_JADESmatchForA2744obj(outfile, matchtol=0.1, overwrite=True, verbose=True):
    """
    Function generating file with pairings of the (GLASS) objects in the A2744 FoV with the mock spectra from
    JADES.

    MUSE redshift of objects are used when available, then GLASS redshifts and lastly ASTRODEEP photo-zs

    --- INPUT ---
    outfile         Name of output file to store list of matches to
    overwrite       Overwrite existing file?
    verbose         Toggle verbosity


    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw
    gw.determine_JADESmatchForA2744obj('A2744_JADESmatches.txt',overwrite=True,verbose=True)

    """
    if verbose: print(' - Loading GLASS catalog to get object IDs and coordinates')
    glasscat   = '/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/a2744_f140w_glasscat_with_ASTRODEEP_f140w_mag.cat'
    glassdat   = np.genfromtxt(glasscat,names=True,dtype=None)
    id_GLASS   = glassdat['NUMBER']
    mag_GLASS  = glassdat['MAG_AUTO']
    ra_GLASS   = glassdat['X_WORLD']
    dec_GLASS  = glassdat['Y_WORLD']

    if verbose: print(' - Loading JADES data ')
    JADESdir     = '/Users/kschmidt/work/catalogs/JADES_GTO/'
    jadesinfo    = fits.open(JADESdir+'JADES_SF_mock_r1_v1.0.fits')[1].data

    if verbose: print(' - Loading redshift catalogs ')
    catdir = '/Users/kschmidt/work/catalogs/'

    ADcat      = catdir + 'ASTRODEEP/fullrelease/A2744cl_26012016/A2744cl_ZPHOT.cat'
    ADdat      = np.genfromtxt(ADcat,names=True,dtype=None)

    ADcat_info = catdir + 'ASTRODEEP/fullrelease/A2744cl_26012016/A2744cl_A.cat'
    ADdat_info = np.genfromtxt(ADcat_info,names=True,dtype=None)

    glasszcat  = catdir + 'GLASSzcat/hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v002_redshiftcatalog.txt'
    glasszdat  = np.genfromtxt(glasszcat,names=True,dtype=None)
    glasszdat  = glasszdat[np.where(glasszdat['redshift'] > 0)]

    MUSEcat    = catdir + 'Mahler18_A2744_redshifts_cat_final.fits'
    MUSEdat    = fits.open(MUSEcat)[1].data


    if verbose: print(' - Loading GALFIT size estiamtesb ')
    sizedir = '/Users/kschmidt/work/observing/proposals_preparation/180406_JWSTcycle1_A2744/galfitFromTaka/'
    size_info = np.genfromtxt(sizedir + 'f160w_01.cat',names=True,dtype=None,skip_header=27)
    size_dat  = np.genfromtxt(sizedir + 'galfit_crs_01.cat',names=True,dtype=None)

    size_id_info  = size_info['NUMBER']
    size_ra       = size_info['X_WORLD']
    size_dec      = size_info['Y_WORLD']

    size_id  = size_dat['1ID']
    size_r50 = size_dat['27r50_sext']
    size_r90 = size_dat['28r90_sext']

    if verbose: print(' - Setting up output file '+outfile)
    if os.path.isfile(outfile) & (overwrite == False):
        sys.exit(outfile+' already exists and "overwrite"=False ')
    else:
        nowstr = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        fout   = open(outfile,'w')
        fout.write('# Best match JADES ids to GLASS objects detected in A2744. \n')
        fout.write('# Generated with grizli_wrappers.determine_JADESmatchForA2744obj() on '+nowstr+' \n')
        fout.write('#     \n')
        fout.write('#    -99      set if no match to id_GLASS in any of redshift catalogs searched \n')
        fout.write('#    -999     set if ASTRODEEP match has id > 100000 in whihc case there is not photo-z estimate \n')
        fout.write('#    -9999    set if match in redshift catalogs but zrange (z+/-0.1) is below JADES redshift, i.e., z_max < 0.2 \n')
        fout.write('#    -99999   set if no GALFIT results for matched (GALFIT) object \n')
        fout.write('#     \n')
        fout.write('# Catalog can be loaded with:   dat = np.genfromtxt("'+outfile+'",names=True,dtype=None,skip_header=10) \n')
        fout.write('#     \n')
        fout.write('# id_GLASS ra_GLASS dec_GLASS f140wmag_GLASS    '
                   'cat_match id_match ra_match dec_match r_match redshift    '
                   'JADESid  JADESz JADESf140wmag '
                   'size_id_info   size_ra_match  size_dec_match  size_r_match  size_id_galfit  size_r50obj_arcsec  size_r90_arcsec\n')

    if verbose: print(' - Loop over objects to get match to JADES mock catalog')
    for oo, objid in enumerate(id_GLASS):
        if (objid < 750) or (objid > 760): continue
        if verbose:
            infostr = '   matching and getting JADES info for id_GLASS = '+str(objid)+' ('+\
                      str("%6.f" % (oo+1))+' / '+str("%6.f" % len(id_GLASS))+')          '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        f140wmag  = mag_GLASS[oo]
        obj_radec = SkyCoord(ra=ra_GLASS[oo]*u.degree, dec=dec_GLASS[oo]*u.degree)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # if verbose: print(' - Crossmatch object to redshift catalogs ')
        AD_radec    = SkyCoord(ra=ADdat_info['RA']*u.degree, dec=ADdat_info['DEC']*u.degree)
        MUSE_radec  = SkyCoord(ra=MUSEdat['RA']*u.degree, dec=MUSEdat['DEC']*u.degree)
        glass_radec = SkyCoord(ra=glasszdat['RA']*u.degree, dec=glasszdat['DEC']*u.degree)
        size_radec = SkyCoord(ra=size_ra*u.degree, dec=size_dec*u.degree)

        # if verbose: print(' - Getting sources within '+str(matchtol)+' arc seconds of redshift catalogs ')
        AD_idx, AD_d2d, AD_d3d          = obj_radec.match_to_catalog_sky(AD_radec)
        MUSE_idx, MUSE_d2d, MUSE_d3d    = obj_radec.match_to_catalog_sky(MUSE_radec)
        glass_idx, glass_d2d, glass_d3d = obj_radec.match_to_catalog_sky(glass_radec)
        size_idx, size_d2d, size_d3d    = obj_radec.match_to_catalog_sky(size_radec)

        AD_idx = np.atleast_1d(AD_idx)[np.where(AD_d2d < matchtol*u.arcsec)[0]]
        AD_d2d = AD_d2d[np.where(AD_d2d < matchtol*u.arcsec)[0]]

        MUSE_idx = np.atleast_1d(MUSE_idx)[np.where(MUSE_d2d < matchtol*u.arcsec)[0]]
        MUSE_d2d = MUSE_d2d[np.where(MUSE_d2d < matchtol*u.arcsec)[0]]

        glass_idx = np.atleast_1d(glass_idx)[np.where(glass_d2d < matchtol*u.arcsec)[0]]
        glass_d2d = glass_d2d[np.where(glass_d2d < matchtol*u.arcsec)[0]]
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # if verbose: print(' - Extract redshift for object (MUSE over GLASS over ASTRODEEP')
        if len(MUSE_idx) > 0:
            bestent    = np.where(MUSE_d2d.arcsec == np.min(MUSE_d2d.arcsec))
            matchent   = MUSE_idx[bestent]
            cat_match  = 'MUSE'
            id_match   = MUSEdat['# ID'][matchent]
            ra_match   = MUSEdat['RA'][matchent]
            dec_match  = MUSEdat['DEC'][matchent]
            r_match    = MUSE_d2d.arcsec[bestent]
            redshift   = MUSEdat['Z'][matchent][0]
        elif len(glass_idx) > 0:
            bestent    = np.where(glass_d2d.arcsec == np.min(glass_d2d.arcsec))
            matchent   = glass_idx[bestent]
            cat_match  = 'GLASS'
            id_match   = glasszdat['ID'][matchent]
            ra_match   = glasszdat['RA'][matchent]
            dec_match  = glasszdat['DEC'][matchent]
            r_match    = glass_d2d.arcsec[bestent]
            redshift   = glasszdat['redshift'][matchent][0]
        elif len(AD_idx) > 0:
            bestent    = np.where(AD_d2d.arcsec == np.min(AD_d2d.arcsec))
            matchent   = AD_idx[bestent]
            cat_match  = 'ASTRODEEP'
            if ADdat_info['ID'][matchent] < 100000:
                if ADdat_info['ID'][matchent] != ADdat['ID'][matchent]:
                    sys.exit(' The IDs of the ASTRODEEP catalogs loaded do not match ')
                else:
                    id_match   = ADdat_info['ID'][matchent]
                redshift   = ADdat['ZBEST'][matchent][0]
            else:
                redshift   = -999
            ra_match   = ADdat_info['RA'][matchent]
            dec_match  = ADdat_info['DEC'][matchent]
            r_match    = AD_d2d.arcsec[bestent]
        else:
            # if verbose: print('   No match to id_GLASS = '+str(objid)+' in any of the redshift catalogs searched             ')
            cat_match  = 'NONE'
            id_match   = [-99]
            ra_match   = [-99]
            dec_match  = [-99]
            r_match    = [-99]
            redshift   = -99

        if redshift > 0:
            # if verbose: print(' - Select best JADES match (closest match to objects F140W magnitude within z_obj +/- 0.25)')
            zrange    = [redshift-0.1,redshift+0.1]
            JADESinfo = ju.get_JADESobjects(redshift=zrange,mag_f140w=[f140wmag,-99],jadesinfo=jadesinfo,verbose=False)
            if len(JADESinfo) == 0:
                JADESid   = [-9999]
                JADESz    = [-9999]
                JADESmag  = [-9999]
            else:
                JADESid   = JADESinfo['ID']
                JADESz    = JADESinfo['redshift']
                JADESmag  = -2.5*np.log10(JADESinfo['HST_F140W_fnu']/1e9)+8.90
        else:
            JADESid   = [-99]
            JADESz    = [-99]
            JADESmag  = [-99]


        # print(' - Add size of nearest object (to GLASS coordinates) from Takahiros estimates ')
        size_bestent    = np.where(size_d2d.arcsec == np.min(size_d2d.arcsec))
        size_matchent   = np.atleast_1d(size_idx)[size_bestent]
        size_ent        = np.where(size_id == size_id_info[size_matchent])[0]
        size_id1        = size_id_info[size_matchent]
        size_ra_match   = size_dec[size_matchent][0]
        size_dec_match  = size_ra[size_matchent][0]
        size_r_match    = size_d2d.arcsec[size_bestent][0]

        if len(size_ent) != 0:
            size_id2        = size_id[size_ent][0]
            size_r50obj     = size_r50[size_ent][0]
            size_r90obj     = size_r90[size_ent][0]
        else:
            size_id2        = -99999
            size_r50obj     = -99999
            size_r90obj     = -99999

        if size_r50obj < 0:
            size_r50obj = -99999
            size_r90obj = -99999
        else:
            size_r50obj = size_r50obj * 0.06
            size_r90obj = size_r90obj * 0.06

        # print(' - Store object information and JADES id to output file')
        outstring = str(objid)+'  '+str(ra_GLASS[oo])+'  '+str(dec_GLASS[oo])+'  '+str("%.4f" % f140wmag)+'  '+\
                    str("%10s" % cat_match)+'  '+str(id_match[0])+'  '+str(ra_match[0])+' '+\
                    str(dec_match[0])+'  '+str(r_match[0])+'  '+str("%.4f" % redshift)+'  '+\
                    str(JADESid[0])+'  '+str(JADESz[0])+'  '+str(JADESmag[0])+'  '+\
                    str(size_id1[0])+'  '+str(size_ra_match)+'  '+str(size_dec_match)+'  '+\
                    str(size_r_match)+'  '+str(size_id2)+'  '+str(size_r50obj)+'  '+str(size_r90obj)
        fout.write(outstring+'\n')
    if verbose: print('\n   ... done')

    fout.close()
    if verbose: print(' - Stored JADES ID matches to '+outfile)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =