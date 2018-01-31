# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Scripts, functions and wrappers for handling grizly reductions and simulations
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import grizli
import grizli.utils
from grizli.multifit import GroupFLT, MultiBeam, get_redshift_fit_defaults
import pdb
from astropy.io import ascii
from importlib import reload

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
    grp = grizli.multifit.GroupFLT(grism_files=glob.glob('niriss_*flt.fits'), direct_files=[],
                                   ref_file=ref_hffimg, ref_ext=0,
                                   seg_file='hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v001_align-drz-seg.fits',
                                   catalog='a2744_f140w_glassmodified.cat',
                                   cpu_count=0, # <0 dont parallelize; =0 use all available; >0 CPUs to use
                                   pad=200)

    print('\n - Analyze simulated data ')
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
    extractid = 793 # GLASS Line emitter at 1.34000 with zQ = 4.00000
    basename  = 'niriss-a2744_'+str("%.5d" % extractid)+'_'

    print(' - Extract spectrum cutouts from individual FLTs of object '+str(extractid))
    beams = grp.get_beams(extractid, size=40)

    print('   Put them in a MultiBeam object')
    mb = grizli.multifit.MultiBeam(beams, fcontam=1, group_name='niriss-a2744')

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
    plt.savefig('./'+basename+'extracted2Dregion.pdf')

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
    plt.savefig('./'+basename+'contaminationModel.pdf')


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
    plt.savefig('./'+basename+'observedSpectrum_model.pdf')

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
    plt.savefig('./'+basename+'fullresiduals.pdf')

    # print(' - Trivial demo model, dropout in the middle of the F115W filter')
    # xspec = np.arange(0.8,2.4,0.02)*1.e4
    # yspec = (xspec/1.4e4)**(-2) # Beta=-2
    # yspec[xspec < 1.2e4] = 0.
    # plt.plot(xspec, yspec)
    # mb.compute_model(spectrum_1d=[xspec, yspec])
    #
    # fig = plt.figure(figsize=[9,9*1.2/3])
    # for ix, i in enumerate([0,2,4,1,3,5]):
    #     ax = fig.add_subplot(2,3,ix+1)
    #     beam = mb.beams[i]
    #     ax.imshow(beam.beam.model, vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')
    #     ax.set_xticklabels([]); ax.set_yticklabels([])
    #     ax.text(0.1,0.1,beam.grism.parent_file, color='k', backgroundcolor='w',
    #             transform=ax.transAxes, size=10, ha='left', va='bottom')
    #
    # fig.axes[0].set_ylabel('Observed\nspectrum')
    # fig.tight_layout(pad=0.1)
    # plt.savefig('./simpeldropoutdemo.pdf')


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
        ax = fig.add_subplot(144)
        ax.imshow(line['LINE', 'Hd'].data, vmin=-0.03, vmax=0.06, cmap=cmap, origin='lower')
        ax.text(5,5,r'H$\delta$', ha='left', va='bottom')

    for ax in fig.axes:
        ax.set_xticklabels([]); ax.set_yticklabels([])

    fig.tight_layout(pad=0.1)

    plt.savefig('./'+basename+'emissionlinemap_Ha.pdf')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print('\n - Ran all commands successfully! ')


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def NIRCAMsim_UDF():
    """

    Simulation on UDF for NIRCAM grisms

    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw
    gw.NIRCAMsim_UDF()

    """
    cwd = os.getcwd()
    print('\n grizli version: %s' %(grizli.__version__))
    print('Simulating NIRCAM grism exposure on UDF')

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

        # JWST NIRCAM, three filters & two orients
        for filt in ['F277W', 'F356W', 'F444W']:
            for theta in [0,90]:
                h, wcs = grizli.fake_image.nircam_header(filter=filt, ra=ra, dec=dec,
                                                         pa_aper=pa_aper+theta)
                print('Filter: {filter}, Background: {bg} e/s/pix, RN: {RN} e/exp'.format(filter=filt,
                                                                bg=h['BACKGR'], RN=h['READN']))
                output = 'nircam_{filt}_{theta:02d}_flt.fits'.format(filt=filt, theta=theta)
                grizli.fake_image.make_fake_image(h, output=output, exptime=EXPTIME, nexp=NEXP)


        print(' - Load GroupFLT for simulation, NB: input files are just noise')
        sim = grizli.multifit.GroupFLT(grism_files=glob.glob('nircam_*flt.fits'), direct_files=[],
                                       ref_file=imagepath+'hlsp_xdf_hst_wfc3ir-60mas_hudf_f140w_v1_sci.fits',
                                       seg_file='udf_f140w_photutils_seg.fits',
                                       catalog='udf_f140w_photutils.cat',
                                       cpu_count=0, pad=200)

        pdb.set_trace() # here 180126

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