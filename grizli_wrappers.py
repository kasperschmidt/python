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
from astropy import wcs
from importlib import reload
import JADESutilities as ju

# plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import matplotlib.gridspec
mpl.rcParams['figure.figsize'] = (10.0, 6.0)
mpl.rcParams['font.size'] = 14
mpl.rcParams['savefig.dpi'] = 144

#importlib.reload(module)

# for grism sims
import grizli.fake_image
from grizli.utils import detect_with_photutils
from collections import OrderedDict
import astropy.io.fits as afits
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
            sci = afits.open(imagepath+'hlsp_xdf_hst_wfc3ir-60mas_hudf_f140w_v1_sci.fits')
            wht = afits.open(imagepath+'hlsp_xdf_hst_wfc3ir-60mas_hudf_f140w_v1_wht.fits')
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
            ref_3dhst = afits.open('udf_3dhst_cat.fits')
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

            afits.writeto('udf_f140w_photutils_seg.fits', data=np.cast[int](seg),
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
            orig_flt = afits.open(flt.grism_file, mode='update')
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
    #line = afits.open('niriss-udf_zfit_00898.line.fits')
    line = afits.open('niriss-udf_00898.line.fits')
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

        # sci_hffimg = afits.open(ref_hffimg)
        # wht_hffimg = afits.open(ref_hffimg.replace('drz.fits','wht.fits'))
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
        # afits.writeto('a2744_f140w_photutils_seg.fits', data=np.cast[int](seg),
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
            orig_flt = afits.open(flt.grism_file, mode='update')
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
    #line = afits.open('niriss-a2744_zfit_00898.line.fits')
    line = afits.open('niriss-a2744_'+str("%.5d" % extractid)+'.line.fits')
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
def NIRCAMsim_A2744(generatesimulation=True, runfulldiagnostics=True, zrangefit=[0.2, 8.0] , mockspec_type='manual_lines',
                    onlyreplaceextractIDs=False, useJADESz=True, fixJADEStemplate=None,
                    singlefilterrun=False, quickrun=False, plotsingleobj=False, workingdir=None,
                    extractids=[65,476,726,233,754,755,793,848],usesamecenter=True):
    """

    NIRCAM simulations of A2744 (based on grizli_wrappers.NIRISSsim_A2744()

    --- INPUT ---
    generatesimulation    Generate the actual simulations (to save time for analysis and plotting,
                          only run this once to generate files etc.)
    runfulldiagnostics    Run the full redshift fits for objects to extract? Otherwise, only spectra extracted.
    zrangefit             Redshift range to fit objects in, if runfulldiagnostics=True
    mockspec_type         What mock spectra to use in the simulation. Choose between:
                             eazy_templates          The SED templates from the SED fits
                             manual_lines            Manually add lines to flat continuum and assign JADES spectra to
                                                     a few selected objects.
                             jades_templatematches   Assign the best-fitting JADES templates (in terms of mF140W and
                                                     redshift) to each of the objects reading matches from catalog
    onlyreplaceextractIDs For the JADES template replacement, to only replace the ids in "extractids" set this
                          keyword to True.
    useJADESz             For the JADES template replacement, to replace the JADES template redshift with the object
                          GLASS/MUSE/ASTRODEEP redshift estimate, set to False. Otherwise JADES redshift is used.
    fixJADEStemplate      For the JADES template replacement, to fix the JADES template (ignoring matches)
                          provide JADES template id here
    singlefilterrun       Only run simulation for F356W to speed up things
    quickrun              If True corners are cut to make the modeling proceed faster (for testing)
    plotsingleobj         Plot the beams of the extracted objects?
    workingdir            Directory to work in. If not set the current working directory is used.
    extractids            List of GLASS ids to extract.
    usesamecenter         To use the same center for the different filters set this to True
                          Otherwise different centers will be used to align "full coverage regions"

    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw
    gw.NIRCAMsim_A2744(generatesimulation=True, runfulldiagnostics=False, mockspec_type='manual_lines')

    objects2extract = [476]

    gw.NIRCAMsim_A2744(generatesimulation=False, mockspec_type='jades_templatematches',singlefilterrun=False, quickrun=False, runfulldiagnostics=False,plotsingleobj=False,extractids=objects2extract)

    objids = [65,380,476,726,983,1110,1222,1528,1669]
    gw.NIRCAMsim_A2744(generatesimulation=True, singlefilterrun=False, quickrun=False, runfulldiagnostics=False,plotsingleobj=False, extractids=objids, mockspec_type='jades_templatematches', onlyreplaceextractIDs=True, useJADESz=False, fixJADEStemplate=695)

    """
    Ncpu      = 1        # <0 dont parallelize; =0 use all available; >0 CPUs to use
    padval    = 2000     #

    # determine polynomial coeffficients for flat continuum model
    polymodelcoeffs = [0.1, 0.05] # coeffs=[1.2, -0.5]))
    #polymodelcoeffs = [0.2, -0.025]
    #polymodelcoeffs = [0.05, 0.00]
    polyxrange      = [1.0, 7.0]

    if workingdir is None:
        cwd = os.getcwd()+'/'
    else:
        print('\n - Setting cwd = '+workingdir)
        if os.path.isdir(workingdir):
            answer = input('   -> that directory already exists, do you want to move there and continue?   ')
            if (answer.lower() == 'y') or (answer.lower() == 'yes'):
                os.chdir(workingdir)
                cwd = workingdir
            else:
                sys.exit('   You did not reply with (y)es so exiting')
        else:
            answer = input('   -> that directory does not exist. Do you want to create it, go there and continue?   ')
            if (answer.lower() == 'y') or (answer.lower() == 'yes'):
                os.mkdir(workingdir)
                os.chdir(workingdir)
                cwd = workingdir
            else:
                sys.exit('   You did not reply with (y)es so exiting')

    if 'nircam' not in cwd.lower():
        sys.exit('"NIRCAM" is not part of workign directory name. Assumes this means the working location is incorrect')
    else:
        print('')
        print('\n --- Simulating NIRCAM grism obs of A2744 in: ')
        print(' --- '+cwd)
        print(' --- Started wrapper at:                       '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print('\n - Using grizli version: %s' %(grizli.__version__))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if singlefilterrun:
        print(' \n\n          NB - Running a "singlefilterrun" ')

        if type(singlefilterrun) == str:
            singlefilt   = singlefilterrun


        else:
            singlefilt  = 'F356W'
        filterloops  = [singlefilt]
        print('               Using filter "'+singlefilt+'" \n\n')
        orientations = [0,90]
        pmc = {}
        pmc['F277W'] = [0.2, -0.025] # coeffs=[1.2, -0.5]))
        pmc['F356W'] = [0.2, -0.025]
        pmc['F444W'] = [0.05, 0.0]
        polymodelcoeffs = pmc[filterloops[0]]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    print(' - Polynomial coefficients to use for model: '+str(polymodelcoeffs))
    print(' - Polynomial range to use for model:        '+str(polyxrange))

    if quickrun:
        print(' \n\n          NB - using the "quickrun" setup \n\n')
        matchtol  = 0.01      # arcsec
        mag_limit  = 18
        mag_limits = [10, mag_limit]
    else:
        matchtol  = 2.0      # arcsec
        mag_limit  = 26.0     # Magnitude limit for traces to simulate and account for
        mag_limits = [10, 26] # Magnitude range for bright objects to refine model for

    print(' - Handle the GLASS object catalog and segmentation image ')
    print('   Define HST filter setup matching GLASS catalog')
    HFFimgfilter    = 'f140w'
    imagepath       = '/Users/kschmidt/work/images_MAST/A2744/'
    ref_hffimg      = imagepath+'hlsp_frontier_hst_wfc3-60mas_abell2744_'+HFFimgfilter+'_v1.0_drz.fits'
    photcatout      = cwd+'../referenceImageAndCat/a2744_f140w_glasscat_with_ASTRODEEP_f140w_mag.cat'
    segmentationmap = cwd+'../referenceImageAndCat/hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v001_align-drz-seg.fits'

    if generatesimulation:
        ra, dec         = 3.588197688, -30.39784202 # Cluster center
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Coordinates providing center of NIRCam module A putting cluster center
        # in "full coverage in both GrismR and GrismC full coverage" area
        radecdic          = {}
        radecdic['F277W'] = 3.5990425, -30.405018
        radecdic['F356W'] = 3.595451,  -30.401799
        radecdic['F444W'] = 3.5822446, -30.390036
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # pa_aper = 135      # Using GLASS PA_V3
        pa_aper = 0.0
        EXPTIME = 6012.591 # seconds exposure (GLASS NIRISS ERS is 5218.07)
        NEXP    = 16       # total Intergrations

        print(' - Generating simulations')
        print(' - Will load HFF '+HFFimgfilter+' images of A2744 from \n   '+imagepath)
        print(' - Loading GLASS A2744 source catalog and corresponding segmentation map')
        cat = np.genfromtxt(cwd+'../referenceImageAndCat/hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v001_glassmaster.txt',dtype=None,names=True)

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

        print(' - Looping over filters and orientations to generate fake images')
        for filt in filterloops:
            for theta in orientations:
                if theta == 0:
                    grism_img = 'GRISMR'
                elif theta == 90:
                    grism_img = 'GRISMC'

                if usesamecenter:
                    pointra  = ra
                    pointdec = dec
                elif (filt in ['F277W','F356W','F444W']) & (theta == 0): # offset f277w row disperion
                    # pointra  = ra277R
                    # pointdec = dec277R
                    pointra, pointdec  = radecdic[filt]
                elif (filt in ['F277W','F356W','F444W']) & (theta == 90): # offset f277w column disperion
                    # pointra  = ra277C
                    # pointdec = dec277C
                    pointra, pointdec  = radecdic[filt]
                else:
                    sys.exit('Nor coordinate set setup found for filt='+str(filt)+' and theta='+str(theta)+' combination')

                h, wcs = grizli.fake_image.nircam_header(filter=filt, ra=pointra, dec=pointdec,
                                                         pa_aper=pa_aper+theta,grism=grism_img)
                print('   Filter: {filter}, Orientation: {th}, Background: {bg} e/s/pix, RN: {RN}'
                      ' e/exp'.format(filter=filt, th=theta, bg=h['BACKGR'], RN=h['READN']))
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

        print(' - Compute full model (for simulation) for '+filt+'; started at:                    '+
              datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) # First pass, flat continuum

        sim.compute_full_model(mag_limit=mag_limit, coeffs=polymodelcoeffs)

        print(' - Run compute_single_model with templates')
        temp_sed_dir = '/Users/kschmidt/work/GLASS/LAEsearchFullGLASS/EAZYruns/180105/' \
                       'A2744_180105/eazy_output_A2744_180105/'
        lamcol       = 'lambda_zprior'
        fluxcol      = 'tempflux_zprior'
        print('   Using wavelength and flux columns '+lamcol+' and '+fluxcol)

        detection_bp = pysynphot.ObsBandpass('wfc3,ir,'+HFFimgfilter)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        print('\n - Building and/or assembling models for individual objects in simulation ')
        if mockspec_type.lower() == 'eazy_templates':
            gw.compute_single_model_EAZY(cat,Nmatches,temp_sed_dir,sim,AD_cat,AD_idx,has_AD_match,lamcol,fluxcol,detection_bp)
        elif mockspec_type == 'manual_lines':
            gw.compute_single_model_MANUAL(sim,detection_bp,has_AD_match,AD_cat,AD_idx,cat,Nmatches)
        elif mockspec_type == 'jades_templatematches':
            if quickrun:
                quickrunIDs = AD_cat['ID'][AD_idx][has_AD_match]
            else:
                quickrunIDs = False

            if onlyreplaceextractIDs:
                quickrunIDs = extractids

            gw.compute_single_model_JADES(sim,cat,detection_bp,selectIDs=quickrunIDs,
                                          useJADESz=useJADESz,fixJADEStemplate=fixJADEStemplate,
                                          JADESmatches='/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/A2744_JADESmatches_201021.txt')
        print('   Done buildinfg and/or assembling models for individual objects in simulation \n')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        subplot_indices = [0,2,4,1,3,5]
        if singlefilterrun:
            if singlefilt == 'F277W': subplot_indices = [0,-99,-99,1,-99,-99]
            if singlefilt == 'F356W': subplot_indices = [-99,0,-99,-99,1,-99]
            if singlefilt == 'F444W': subplot_indices = [-99,-99,0,-99,-99,1]
        #---------------------------------------------------------------------------
        Fsize = 8.0
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        #---------------------------------------------------------------------------
        print(' - Plot blotted reference image; "grism" exposures are still empty, just noise')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate(subplot_indices):
            if i != -99:
                ax = fig.add_subplot(2,3,ix+1)
                ax.imshow(sim.FLTs[i].direct['REF'], vmin=-0.01, vmax=0.05, cmap='viridis',
                          origin='lower')
                ax.set_xticklabels([]); ax.set_yticklabels([])
                ax.grid(color='w', alpha=0.8)
                ax.text(100,100,sim.FLTs[i].grism_file.replace('_','\_'), color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./blotted_reference_image.pdf')
        #---------------------------------------------------------------------------
        print(' - Plot blotted segmentation image; "grism" exposures are still empty, just noise')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate(subplot_indices):
            if i != -99:
                ax = fig.add_subplot(2,3,ix+1)
                ax.imshow(sim.FLTs[i].seg, vmin=-0.01, vmax=3000, cmap='viridis', origin='lower')
                ax.set_xticklabels([]); ax.set_yticklabels([])
                ax.grid(color='w', alpha=0.8)
                ax.text(100,100,sim.FLTs[i].grism_file.replace('_','\_'), color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./blotted_segmentation_image.pdf')
        #---------------------------------------------------------------------------
        print(' - Plot "grism" exposures that are still empty, just noise')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate(subplot_indices):
            if i != -99:
                ax = fig.add_subplot(2,3,ix+1)
                ax.imshow(sim.FLTs[i].grism['SCI'], vmin=-0.01, vmax=0.05, cmap='viridis', origin='lower')
                ax.set_xticklabels([]); ax.set_yticklabels([])
                ax.grid(color='w', alpha=0.8)
                ax.text(100,100,sim.FLTs[i].grism_file.replace('_','\_'), color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./blotted_grismnoise_image.pdf')
        #---------------------------------------------------------------------------
        print(' - Plot model stored in FLTs[i].model attribute')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate(subplot_indices):
            if i != -99:
                ax = fig.add_subplot(2,3,ix+1)
                # show as if it were the rotated grism
                if (i % 2) > 0:
                    ax.imshow(np.rot90(sim.FLTs[i].model,-1), vmin=-0.01, vmax=0.05, cmap='viridis',
                              origin='lower')
                else:
                    ax.imshow(sim.FLTs[i].model, vmin=-0.01, vmax=0.05, cmap='viridis', origin='lower')

                ax.set_xticklabels([]); ax.set_yticklabels([])
                ax.grid(color='w', alpha=0.8)
                ax.text(100,100,sim.FLTs[i].grism_file.replace('_','\_'), color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./flt_model_attribute.pdf')
        #---------------------------------------------------------------------------
        print(' - Update SCI extension of the fake FLT images with the models just computed')
        for flt in sim.FLTs:
            print('Update', flt.grism_file)
            orig_flt = afits.open(flt.grism_file, mode='update')
            orig_flt['SCI'].data += flt.model[flt.pad:-flt.pad, flt.pad:-flt.pad]
            orig_flt.flush()

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
        #print('   Refine the (polynomial) continuum model for brighter objects')
        # grp.refine_list(poly_order=len(polymodelcoeffs), mag_limits=mag_limits, verbose=False, wave = np.linspace(0.2,6.0e4,100))
        #grp.refine_list(poly_order=3, mag_limits=mag_limits, verbose=False, wave = np.linspace(0.2,6.0e4,100))

        # print(' - saving refined contamination model')
        grp.save_full_data()

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        GrismFLTimgs = glob.glob('nircam*.GrismFLT.fits')
        for gfi in GrismFLTimgs:
            outfile  = gfi.replace('GrismFLT.fits','GrismFLT_modeldiff.fits')
            gw.create_diffimg(gfi,'GSCI',gfi,'MODEL',outfile,overwrite=True,header=True)

        #---------------------------------------------------------------------------
        print(' - Plot GrismFLT_modeldiff.fits images')
        fig = plt.figure(figsize=[9,9*2./3])
        for ix, i in enumerate(subplot_indices):
            if i != -99:
                datfile = GrismFLTimgs[i].replace('GrismFLT.fits','GrismFLT_modeldiff.fits')
                ax = fig.add_subplot(2,3,ix+1)
                diffdata = afits.open(datfile)[0].data

                # show as if it were the rotated grism
                if (i % 2) > 0:
                    ax.imshow(np.rot90(diffdata,-1), vmin=-0.01, vmax=0.05, cmap='viridis', origin='lower')
                else:
                    ax.imshow(diffdata, vmin=-0.01, vmax=0.05, cmap='viridis', origin='lower')
                ax.set_xticklabels([]); ax.set_yticklabels([])
                ax.grid(color='w', alpha=0.8)
                ax.text(100,100,datfile.split('/')[-1].replace('_','\_'), color='w', size=10, ha='left', va='bottom')
        fig.tight_layout(pad=0.1)
        plt.savefig('./GrismFLT_modeldiff_files.pdf')
        #---------------------------------------------------------------------------
        if not singlefilterrun:
            splitstr = 'nircam_'
            basename = GrismFLTimgs[0].split(splitstr)[0]+splitstr
            gw.plot_grismFoV(basename,edgecut=1900,verbose=True)
        #---------------------------------------------------------------------------
    else:
        print('\n - NB: Going directly to analysis of simulated data (assuming they exist)')

        print('   Reloading simulations (and initializing grp class for analysis)')
        grp = grizli.multifit.GroupFLT(grism_files=glob.glob('nircam_*flt.fits'), direct_files=[],
                                       ref_file=ref_hffimg, ref_ext=0,
                                       seg_file=segmentationmap,
                                       catalog=photcatout,
                                       cpu_count=Ncpu, # <0 dont parallelize; =0 use all available; >0 CPUs to use
                                       pad=padval,
                                       polyx=polyxrange) # range the polynomial model is calculated over in microns

    if extractids is not None:
        print('\n - Analyze simulated data; started at:             '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ## emission line object
        # Bright central galaxy:  ID_00848 (RA,DEC) = (3.5877,-30.3964) GLASS redshift = 0.316 | redshift quality = 2.0
        # Pa-a in F277 ID_00065 (RA,DEC) = (3.5770,-30.3795) GLASS redshift = 0.496 | redshift quality = 4.0
        # BV Ha map; Pa-a in F277 ID_00476 (RA,DEC) = (3.6044,-30.3850) GLASS redshift = 0.300 | redshift quality = 4.0
        # BV Ha map; Pa-a in F277 ID_00726 (RA,DEC) = (3.5809,-30.3908) GLASS redshift = 0.292 | redshift quality = 4.0
        # MUSE LAE                ID_00233 (RA,DEC) = (3.5778434,-30.3812148) MUSE ID = 14518; MUSE redshift = 3.3756
        # MUSE LAE 'arc'          ID_00754 (RA,DEC) = (3.5887928,-30.3938037) MUSE IDs = 8826 & 8789; MUSE redshift = 3.9803
        # Line emitter:           ID_00755 (RA,DEC) = (3.5874,-30.3933) GLASS redshift = 0.600 | redshift quality = 3.0
        # Line emitter:           ID_00793 (RA,DEC) = (3.6061,-30.3920) GLASS redshift = 1.340 | redshift quality = 4.0
        # YD4:                    ID_00292 (RA,DEC) = (3.603853194,-30.382264279) confirmed z=8.32 galaxy (Schmidt+16 ID 2193)
        # XinZmapCand:            ID_00783  with zASTRODEEP = 6.2093

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

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            print(' - Run the redshift fit and generate the emission line map?')
            if runfulldiagnostics:
                print('   -> Yes')
                fit, fig, fig2, hdu2, hdu_line = gw.run_zfit(beamfile,zrangefit=zrangefit,grp=grp,verbose=True)
            else:
                print('   -> No')
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            print(' - Plot beams for individual objects?')
            if plotsingleobj:
                print('   -> Yes')
                gw.plot_beams(beamfile)
            else:
                print('   -> No')
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            print(' - Plot emission line map?')
            if runfulldiagnostics:
                print('   -> Yes')
                gw.plot_ELmaps(base_name+'.line.fits', map_vmin = -0.03, map_vmax = 0.08)
            else:
                print('   -> No')
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            print('\n v v v v v v v v v v DS9 command to display beams  v v v v v v v v v v ')
            beamfile_nocontam = './'+base_name+'.beams_nocontam.fits'
            beamhdu           = afits.open(beamfile)
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
def run_zfit(beamfile,zrangefit=[0.1,10.0],dzfit=[0.01, 0.001],grp=None,plotmaps=True,verbose=True):
    """
    Run redshift fit estimate on *.beams.fits files

    --- INPUT ---
    mb         a multi beam object

    --- EXAMPLE OF USE ---
    beamfile = '/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/180314_coeff0p2and-0p025/nircam-a2744_00726.beams.fits'
    fit, fig, fig2, hdu2, hdu_line = gw.run_zfit(beamfile,zrangefit=[0.1,0.5])

    """
    if verbose: print(' - Loading beams in \n   '+beamfile)
    base_name = beamfile.split('_0')[0]
    mb        = grizli.multifit.MultiBeam(beamfile, fcontam=1, group_name=base_name)

    if verbose: print(' - Get default fit parameters\n')
    pzfit, pspec2, pline = grizli.multifit.get_redshift_fit_defaults()

    # Redshift fit params
    pzfit ['zr'] = zrangefit
    pzfit['dz'] = dzfit

    # Drizzled line maps params
    pline = {'kernel': 'square', 'pixfrac': 0.8, 'pixscale': 0.063, 'size': 10}

    # Full rectified 2D spectrum params
    pspec2 = {'NY': 20, 'dlam': 50, 'spatial_scale': 1}

    out = mb.run_full_diagnostics(pzfit=pzfit, pspec2=pspec2, pline=pline,
                                  GroupFLT=grp, prior=None, verbose=True)

    fit, fig, fig2, hdu2, hdu_line = out

    if plotmaps:
        gw.plot_ELmaps(beamfile.replace('.beams.','.line.'), map_vmin = -0.03, map_vmax = 0.08)

    return out
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_grismFoV(basename,edgecut=0.0,markbeams=[None],cmap='viridis',verbose=True):
    """
    Function plotting simulation full-FoV output with direct image of GRISMR and GRISMC renderings

    --- INPUT ---
    basename      The path and main part of the output files. Will then look for the simulated GrismFLT files
                  corresponding to the 3 NIRCam filters.
    edgecut       Remove the edges of the images (correcting for the padval in the simulations)
    markbeams     Provide list of beam files to mark location in images

    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw, glob
    workdir   = '/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/201002_basicsimulation/'
    basename  = workdir+'nircam_'
    beamfiles = glob.glob(workdir+'*.beams.fits')
    gw.plot_grismFoV(basename,markbeams=beamfiles,edgecut=1900)

    """
    outfile = basename+'grismFoV.pdf'
    if verbose: print(' - Generating full-FoV figure of NIRCam grisms. Storing in\n   '+outfile)

    grismFLTs_F277W = glob.glob(basename+'F277W*.GrismFLT.fits')
    grismFLTs_F356W = glob.glob(basename+'F356W*.GrismFLT.fits')
    grismFLTs_F444W = glob.glob(basename+'F444W*.GrismFLT.fits')

    fig   = plt.figure(figsize=[12,4])
    Fsize = 10.0
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()


    datalist = [afits.open(grismFLTs_F277W[0])['DREF'].data[0+edgecut:-1*edgecut,0+edgecut:-1*edgecut]*5e19,
                afits.open(grismFLTs_F277W[0])['GSCI'].data[0+edgecut:-1*edgecut,0+edgecut:-1*edgecut],
                afits.open(grismFLTs_F356W[0])['DREF'].data[0+edgecut:-1*edgecut,0+edgecut:-1*edgecut]*5e19,
                afits.open(grismFLTs_F356W[0])['GSCI'].data[0+edgecut:-1*edgecut,0+edgecut:-1*edgecut],
                afits.open(grismFLTs_F444W[0])['DREF'].data[0+edgecut:-1*edgecut,0+edgecut:-1*edgecut]*5e19,
                afits.open(grismFLTs_F444W[0])['GSCI'].data[0+edgecut:-1*edgecut,0+edgecut:-1*edgecut],
                afits.open(grismFLTs_F277W[1])['DREF'].data[0+edgecut:-1*edgecut,0+edgecut:-1*edgecut]*5e19,
                afits.open(grismFLTs_F277W[1])['GSCI'].data[0+edgecut:-1*edgecut,0+edgecut:-1*edgecut],
                afits.open(grismFLTs_F356W[1])['DREF'].data[0+edgecut:-1*edgecut,0+edgecut:-1*edgecut]*5e19,
                afits.open(grismFLTs_F356W[1])['GSCI'].data[0+edgecut:-1*edgecut,0+edgecut:-1*edgecut],
                afits.open(grismFLTs_F444W[1])['DREF'].data[0+edgecut:-1*edgecut,0+edgecut:-1*edgecut]*5e19,
                afits.open(grismFLTs_F444W[1])['GSCI'].data[0+edgecut:-1*edgecut,0+edgecut:-1*edgecut]]

    namelist = [grismFLTs_F277W[0]+'DREF',
                grismFLTs_F277W[0]+'GSCI',
                grismFLTs_F356W[0]+'DREF',
                grismFLTs_F356W[0]+'GSCI',
                grismFLTs_F444W[0]+'DREF',
                grismFLTs_F444W[0]+'GSCI',
                grismFLTs_F277W[1]+'DREF',
                grismFLTs_F277W[1]+'GSCI',
                grismFLTs_F356W[1]+'DREF',
                grismFLTs_F356W[1]+'GSCI',
                grismFLTs_F444W[1]+'DREF',
                grismFLTs_F444W[1]+'GSCI']


    beamSCIheaders  = []
    coordlist       = []
    for beamfile in markbeams:
        if beamfile is not None:
            bhdu            = afits.open(beamfile)
            for hdu in bhdu[1:]:
                if (hdu.header['EXTNAME'] == 'SCI') or (hdu.header['EXTNAME'] == 'REF'):
                    beamSCIheaders.append(hdu.header)
                    coordlist.append([afits.open(beamfile)[0].header['RA'],afits.open(beamfile)[0].header['DEC']])

    for ii, datashow in enumerate(datalist):
        ax = fig.add_subplot(2,6,ii+1)
        # goodpix = datashow[np.isfinite(datashow) & (datashow!=0.0)]
        # pixsort = np.sort(goodpix)
        ax.imshow(datashow, vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')

        for bb, beamhdr in enumerate(beamSCIheaders):
            if (beamhdr['PARENT'].split('_flt.fits')[0] in namelist[ii]):
                if namelist[ii].endswith('REF') & (beamhdr['EXTNAME'] == 'REF'):
                    xpix, ypix = beamhdr['ORIGINX']-edgecut, beamhdr['ORIGINY']-edgecut
                    width      = beamhdr['NAXIS1']
                    height     = beamhdr['NAXIS2']

                    radius     = 5
                    objpatch   = patch.Circle((xpix+width/2.,ypix+height/2.),radius,color='red',fill=None)
                    ax.add_patch(objpatch)

                    objpatch   = patch.Rectangle((xpix,ypix), width, height, color='red', linewidth=1.0, linestyle='-', fill=None)
                    ax.add_patch(objpatch)

                    skipWCSpatch = True
                    if not skipWCSpatch:
                        # Base patch location on WCS and coordinates instead of header info
                        wcs_GrismFLT  = wcs.WCS(afits.open(namelist[ii][:-4])[namelist[ii][-4:]].header)
                        coords        = SkyCoord(coordlist[bb][0],coordlist[bb][1], unit="deg")
                        pixcoord      = wcs.utils.skycoord_to_pixel(coords,wcs_GrismFLT,origin=1)
                        xpix, ypix    = pixcoord[0]-edgecut, pixcoord[1]-edgecut
                        radius        = 30
                        objpatch      = patch.Circle((xpix,ypix),radius,color='white',fill=None)
                        ax.add_patch(objpatch)

                if namelist[ii].endswith('SCI') & (beamhdr['EXTNAME'] == 'SCI'):
                    # xpix, ypix = beamhdr['CRPIX1']-edgecut, beamhdr['CRPIX2']-edgecut
                    xpix, ypix = beamhdr['ORIGINX']-edgecut, beamhdr['ORIGINY']-edgecut
                    width      = beamhdr['NAXIS1']
                    height     = beamhdr['NAXIS2']
                    objpatch   = patch.Rectangle((xpix,ypix), width, height, color='red', linewidth=1.0, linestyle='-', fill=None)
                    ax.add_patch(objpatch)

        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xticklabels([]); ax.set_yticklabels([])

    fig.axes[0].set_title('F277W REF IMG')
    fig.axes[2].set_title('F356W REF IMG')
    fig.axes[4].set_title('F444W REF IMG')

    fig.axes[1].set_title('F277W GRISM')
    fig.axes[3].set_title('F356W GRISM')
    fig.axes[5].set_title('F444W GRISM')

    fig.axes[0].set_ylabel('GRISM R')
    fig.axes[6].set_ylabel('GRISM C')

    fig.tight_layout(pad=0.1)
    plt.savefig(outfile)
    plt.clf()
    plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_beams(beamfile,cmap='viridis_r',verbose=True):
    """
    Plot beams in *.beams.fits files

    --- INPUT ---
    mb         a multi beam object

    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw
    beamfile = '/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/200930_basicsimulation/nircam-a2744_00380.beams.fits'
    gw.plot_beams(beamfile)

    """
    if verbose: print(' - Loading beams in \n   '+beamfile)
    mb     = grizli.multifit.MultiBeam(beamfile) # providing file name invokes load_master_fits(beamfile)
    # mb.grisms
    # mb.PA
    Nbeams = len(mb.beams)
    if verbose: print(' - Found '+str(Nbeams)+' beams in beamfile')

    bhdu   = afits.open(beamfile)

    beamfiles   = []
    model_hdus  = []
    contam_hdus = []
    for hdu in bhdu[1:]:
        if hdu.header['EXTNAME'] == 'SCI':
            beamfiles.append(hdu.header['GPARENT'])
        if hdu.header['EXTNAME'] == 'MODEL':
            model_hdus.append(hdu)
        if hdu.header['EXTNAME'] == 'CONTAM':
            contam_hdus.append(hdu)

    subplot_indices = [-99,-99,-99,-99,-99,-99]
    for bb, bfile in enumerate(beamfiles):
        if bfile == 'nircam_F277W_00_flt.fits':
            subplot_indices[0] = bb
        elif bfile == 'nircam_F277W_90_flt.fits':
            subplot_indices[3] = bb
        elif bfile == 'nircam_F356W_00_flt.fits':
            subplot_indices[1] = bb
        elif bfile == 'nircam_F356W_90_flt.fits':
            subplot_indices[4] = bb
        elif bfile == 'nircam_F444W_00_flt.fits':
            subplot_indices[2] = bb
        elif bfile == 'nircam_F444W_90_flt.fits':
            subplot_indices[5] = bb
        else:
            print('\n\n WARNING - the beam parent file '+str(bfile)+
                  ' has not plot index assoiciated with it in gw.plot_beams(); hence ignored...\n\n')

    if Nbeams == 6:
        subplot_indices_ideal = [0,2,4,1,3,5]
        if subplot_indices != subplot_indices_ideal:
            print('\n\n WARNING - the 6 beams were not arranged as intended \n\n')
    # else:
    #     subplot_indices = [0,0,0,0,0,0]
    #     for ii, beam in enumerate(mb.beams):
    #         subplot_indices[ii] = ii

    beamtext     = [beam.grism.parent_file for beam in mb.beams]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - "Beams" are extracted for spectra of a given order.  Have attributes for contam, model etc.')
    outputfile   = beamfile.replace('.beams.fits','_spectrum2D_observed.pdf')
    beamdatalist = [beam.grism['SCI'] for beam in mb.beams]
    gw.plot_beams_add_figure(outputfile,beamdatalist,beamtext,subplot_indices,'sci',cmap=cmap)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Each beam carries with it a static contamination model extracted from the full field')
    outputfile   = beamfile.replace('.beams.fits','_contamination_model.pdf')
    # beamdatalist = [beam.contam for beam in mb.beams]
    contamdatalist = [contam.data for contam in contam_hdus]
    gw.plot_beams_add_figure(outputfile,contamdatalist,beamtext,subplot_indices,'contam',cmap=cmap)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Plotting difference between contam in multi beam and the beam fits file')
    outputfile   = beamfile.replace('.beams.fits','_contamination_model_multibeamVSfits.pdf')
    contamdiff     = [contam.data - mb.beams[cc].contam for cc,contam in enumerate(contam_hdus)]
    gw.plot_beams_add_figure(outputfile,contamdiff,beamtext,subplot_indices,'contam',cmap=cmap)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - "Beams" are extracted for spectra of a given order.  Have attributes for contam, model etc.')
    outputfile   = beamfile.replace('.beams.fits','_spectrum2D_observed_contamremove.pdf')
    beamdatalist = [beam.grism['SCI']-beam.contam for beam in mb.beams]
    gw.plot_beams_add_figure(outputfile,beamdatalist,beamtext,subplot_indices,'sci-contam',cmap=cmap)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Under the hood, the fitting is done by specifying a single 1D template, which')
    if verbose: print('   is used to generate model 2D spectra for each beam')
    outputfile    = beamfile.replace('.beams.fits','_spectrum2D_model.pdf')
    # beamdatalist  = [beam.model for beam in mb.beams]
    modeldatalist = [model.data for model in model_hdus]
    gw.plot_beams_add_figure(outputfile,modeldatalist,beamtext,subplot_indices,'model',cmap=cmap)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Plotting difference between model in multi beam and the beam fits file')
    outputfile    = beamfile.replace('.beams.fits','_spectrum2D_model_multibeamVSfits.pdf')
    modeldiff     = [mb.beams[cc].model - model.data for cc,model in enumerate(model_hdus)]
    gw.plot_beams_add_figure(outputfile,modeldiff,beamtext,subplot_indices,'model',cmap=cmap)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Goodness of fit is computed by comparing the models in the full 2D pixel space')
    outputfile   = beamfile.replace('.beams.fits','_fullresiduals.pdf')
    beamdatalist = [beam.grism['SCI'] - beam.contam - modeldatalist[bb] for bb, beam in enumerate(mb.beams)]
    gw.plot_beams_add_figure(outputfile,beamdatalist,beamtext,subplot_indices,'sci-contam-model',cmap=cmap)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_beams_add_figure(outputfile,beamdatalist,beamtext,panelindices,ylabel,cmap='viridis_r'):
    """
    Wrapper to add figure of individual beam content

    """
    Npanelcols = 3
    fig   = plt.figure(figsize=[9,9/Npanelcols])
    Fsize = 8.0
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    for ix, beamindex in enumerate(panelindices):
        ax = fig.add_subplot(2,Npanelcols,ix+1)
        if beamindex == -99:
            beamdata = beamdatalist[0]*0.0
            ax.imshow(beamdata, vmin=-0.01, vmax=0.05, cmap='Greys_r', origin='lower', aspect='auto')
            btext = 'No beam coverage'
        else:
            beamdata = beamdatalist[beamindex]
            ax.imshow(beamdata, vmin=-0.01, vmax=0.05, cmap=cmap, origin='lower', aspect='auto')
            btext = beamtext[beamindex].replace('_','\_')

        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(0.05,0.1,btext, color='k', backgroundcolor='w', transform=ax.transAxes, size=Fsize, ha='left', va='bottom')

    fig.axes[0].set_ylabel(ylabel)
    fig.axes[3].set_ylabel(ylabel)
    fig.tight_layout(pad=0.1)
    plt.savefig(outputfile)
    plt.clf()
    plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_ELmaps(linefile, map_vmin=-0.03, map_vmax=0.06, wht_vmin=-0.01, wht_vmax=20000, colormap='rainbow'):
    """
    Plot emission line maps in *.line.fits files

    --- EXAMPLE OF USE ---
    linefile = '/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/nircam-a2744_00783.line.fits'
    linefile = '/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/nircam-a2744_00726.line.fits'
    gw.plot_ELmaps(linefile, map_vmin = -0.03, map_vmax = 0.08)

    """
    line = afits.open(linefile)
    maps = line[0].header['HASLINES'].split()
    line.info()

    latexnames = {}
    latexnames['PaA']       = [r'Pa$\alpha$',              r'Weight', r'Continuum', r'Contam']
    latexnames['PaB']       = [r'Pa$\beta$',               r'Weight', r'Continuum', r'Contam']
    latexnames['PaG']       = [r'Pa$\gamma$',              r'Weight', r'Continuum', r'Contam']
    latexnames['BrA']       = [r'Br$\alpha$',              r'Weight', r'Continuum', r'Contam']
    latexnames['BrB']       = [r'Br$\beta$',               r'Weight', r'Continuum', r'Contam']
    latexnames['BrG']       = [r'Br$\gamma$',              r'Weight', r'Continuum', r'Contam']
    latexnames['FeII']      = [r'FeII$\lambda$16440',      r'Weight', r'Continuum', r'Contam']
    latexnames['HeI-1083']  = [r'HeI$\lambda$1083',        r'Weight', r'Continuum', r'Contam']
    latexnames['SIII']      = [r'SIII',                    r'Weight', r'Continuum', r'Contam']
    latexnames['SII']       = [r'SII',                     r'Weight', r'Continuum', r'Contam']
    latexnames['Ha']        = [r'H$\alpha$',               r'Weight', r'Continuum', r'Contam']
    latexnames['OI-6302']   = [r'OI$\lambda$6302',         r'Weight', r'Continuum', r'Contam']
    latexnames['OIII']      = [r'[OIII]$\lambda$5007',     r'Weight', r'Continuum', r'Contam']
    latexnames['Hb']        = [r'H$\beta$',                r'Weight', r'Continuum', r'Contam']
    latexnames['OIII-4363'] = [r'[OIII]$\lambda$4363',     r'Weight', r'Continuum', r'Contam']
    latexnames['Hg']        = [r'H$\gamma$',               r'Weight', r'Continuum', r'Contam']
    latexnames['Hd']        = [r'H$\delta$',               r'Weight', r'Continuum', r'Contam']
    latexnames['NeIII']     = [r'NeIII',                   r'Weight', r'Continuum', r'Contam']
    latexnames['OII']       = [r'OII$\lambda$3726',        r'Weight', r'Continuum', r'Contam']
    latexnames['NeVI']      = [r'NeVI',                    r'Weight', r'Continuum', r'Contam']
    latexnames['NeV']       = [r'NeV',                     r'Weight', r'Continuum', r'Contam']
    latexnames['MgII']      = [r'MgII',                    r'Weight', r'Continuum', r'Contam']
    latexnames['CIV-1549']  = [r'CIV$\lambda$1549',        r'Weight', r'Continuum', r'Contam']
    latexnames['CIII-1908'] = [r'CIV$\lambda$1908',        r'Weight', r'Continuum', r'Contam']
    latexnames['OIII-1663'] = [r'OIII$\lambda$1663',       r'Weight', r'Continuum', r'Contam']
    latexnames['HeII-1640'] = [r'HeII$\lambda$1640',       r'Weight', r'Continuum', r'Contam']
    latexnames['NIII-1750'] = [r'NIII$\lambda$1750',       r'Weight', r'Continuum', r'Contam']
    latexnames['NIV-1487']  = [r'NIV$\lambda$1487',        r'Weight', r'Continuum', r'Contam']
    latexnames['NV-1240']   = [r'NV$\lambda$1240',         r'Weight', r'Continuum', r'Contam']
    latexnames['Lya']       = [r'Ly$\alpha$',              r'Weight', r'Continuum', r'Contam']

    Nmaps    = len(maps)
    Nrows    = Nmaps+1
    Ncols    = 4
    FS       = 8

    # cmap = 'cubehelix_r'
    # cmap = 'rainbow'
    fig = plt.figure(figsize=[4,Nmaps+1])

    ax = fig.add_subplot(Nrows, Ncols, 1)
    ax.imshow(line['DSCI'].data, vmin=map_vmin, vmax=map_vmax, cmap=colormap, origin='lower')
    ax.text(5,5,'Direct '+line['DSCI'].header['FILTER'], ha='left', va='bottom', fontsize=FS)

    ax = fig.add_subplot(Nrows, Ncols, 2)
    ax.imshow(line['DWHT'].data, vmin=wht_vmin, vmax=wht_vmax, cmap='gray', origin='lower')
    ax.text(5,5,r'Direct weight', ha='left', va='bottom', color='w', fontsize=FS)

    for mm, map in enumerate(maps):
        ax = fig.add_subplot(Nrows, Ncols, 1+4*(mm+1))
        ax.imshow(line['LINE', map].data, vmin=map_vmin, vmax=map_vmax, cmap=colormap, origin='lower')
        ax.text(5,5,latexnames[map][0], ha='left', va='bottom', fontsize=FS)

        ax = fig.add_subplot(Nrows, Ncols, 2+4*(mm+1))
        ax.imshow(line['LINEWHT', map].data, vmin=wht_vmin, vmax=wht_vmax, cmap='gray', origin='lower')
        ax.text(5,5,latexnames[map][1], ha='left', va='bottom', color='w', fontsize=FS)

        ax = fig.add_subplot(Nrows, Ncols, 3+4*(mm+1))
        ax.imshow(line['CONTINUUM', map].data, vmin=map_vmin, vmax=map_vmax, cmap=colormap, origin='lower')
        ax.text(5,5,latexnames[map][2], ha='left', va='bottom', fontsize=FS)

        ax = fig.add_subplot(Nrows, Ncols, 4+4*(mm+1))
        ax.imshow(line['CONTAM', map].data, vmin=map_vmin, vmax=map_vmax, cmap=colormap, origin='lower')
        ax.text(5,5,latexnames[map][3], ha='left', va='bottom', fontsize=FS)

    for ax in fig.axes:
        ax.set_xticklabels([]); ax.set_yticklabels([])

    fig.tight_layout(pad=0.1)

    plt.savefig(linefile.replace('.fits','_emissionlinemaps.pdf'))
    plt.clf()
    plt.close('all')

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

    # adding a Gaussian emission line to spectrum
    # EL = pysynphot.GaussianSource(flux=line_flux,center=line_center,fwhm=fwhm_observed,)
    # spec += EL

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
def compute_single_model_JADES(sim,cat,detection_bp,useJADESz=True,fixJADEStemplate=None,selectIDs=False,
                               JADESmatches='/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/A2744_JADESmatches_201021.txt'):
    """

    """
    print(' - Assigning JADES templates as models for individual objects using \n   '+JADESmatches)

    dat = np.genfromtxt(JADESmatches,names=True,dtype=None,skip_header=10)
    # for ii, id in enumerate(dat['JADESid']):
    #     if (id > 0) & (dat['JADESz'][ii] < 0.203):
    #         print(str(id)+' '+str(dat['JADESz'][ii])+' '+str(dat['JADESf140wmag'][ii]))

    for ii, id in enumerate(cat['NUMBER']):
        if selectIDs is not False:
            if id not in selectIDs: continue
        catent = np.where(dat['id_GLASS'].astype(int) == id)[0]

        if len(catent) == 0:
            sys.exit(' Match to JADES catalog not performed for id_GLASS = '+str(id))

        Jobj = dat['JADESid'][ii]

        if fixJADEStemplate is not None:
            GLASSmapsIDs = [65,380,476,726,983,1110,1222,1528,1669]
            print('\n   Fixing JADES template to:')
            Jobj = fixJADEStemplate
            print('       idJADES  = '+str(Jobj))

        if Jobj == -99:
            Jobj   = 274533
            modstr = 'JADES template'+str("%9s" % Jobj)+' (Manual: z=6.255 & m140=25.45)'
        elif Jobj == -9999:
            Jobj   = 4
            modstr = 'JADES template'+str("%9s" % Jobj)+' (Manual: z=0.200 & m140=27.26)'
        else:
            modstr = 'JADES template'+str("%9s" % Jobj)+' (ID assigned in z/JADES match)'

        #print(' - Loading JADES spectrum for id_GLASS = '+str(id))
        if useJADESz:
            JADESinfo, temp_lambda, temp_flux = ju.get_JADESspecAndInfo(Jobj,observedframe=True,verbose=False)
        else:
            # fixing the JADES object to use
            fixredshift  = dat['redshift'][catent]
            print('\n   Not using JADES redshift (using z-obj) fixing spectrum for GLASS Ha-map obj to:')
            print('       idGLASS  = '+str(id))
            print('       z-obj    = '+str(fixredshift))
            JADESinfo, temp_lambda, temp_flux = ju.get_JADESspecAndInfo(Jobj,observedframe=False,verbose=False)

            temp_lambda  = temp_lambda  * (1 + fixredshift)
            temp_flux    = temp_flux    / (1 + fixredshift)

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
    beamhdu = afits.open(beamfile)
    Nsetup  = beamhdu[0].header['COUNT']

    for ext in np.arange(Nsetup):
        beamhdu[3+ext*7].data              = beamhdu[3+ext*7].data - beamhdu[6+ext*7].data
        beamhdu[3+ext*7].header['EXTNAME'] = 'SCINOCONTAM'

    beamhdu.writeto(outname,overwrite=overwrite)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def determine_JADESmatchForA2744obj(outfile, matchtol=0.1, overwrite=True, verbose=True, testing=False):
    """
    Function generating file with pairings of the (GLASS) objects in the A2744 FoV with the mock spectra from
    JADES.

    MUSE redshift of objects are used when available, then GLASS redshifts and lastly ASTRODEEP photo-zs

    --- INPUT ---
    outfile         Name of output file to store list of matches to
    overwrite       Overwrite existing file?
    verbose         Toggle verbosity
    testing         Set this to true if testing the code (more info printed and only limited number of objects matched)


    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw
    outfile = '/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/A2744_JADESmatches.txt'
    gw.determine_JADESmatchForA2744obj(outfile,overwrite=True,verbose=True,testing=True)

    """
    if verbose: print(' - Loading GLASS catalog to get object IDs and coordinates')
    glasscat   = '/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/referenceImageAndCat/a2744_f140w_glasscat_with_ASTRODEEP_f140w_mag.cat'
    glassdat   = np.genfromtxt(glasscat,names=True,dtype=None)
    id_GLASS   = glassdat['NUMBER']
    mag_GLASS  = glassdat['MAG_AUTO']
    ra_GLASS   = glassdat['X_WORLD']
    dec_GLASS  = glassdat['Y_WORLD']

    if verbose: print(' - Loading JADES data ')
    JADESdir     = '/Users/kschmidt/work/catalogs/JADES_GTO/'
    jadesinfo    = afits.open(JADESdir+'JADES_SF_mock_r1_v1.0.fits')[1].data

    if verbose: print(' - Loading redshift catalogs ')
    catdir = '/Users/kschmidt/work/catalogs/'

    ADcat      = catdir + 'ASTRODEEP/fullrelease/A2744cl_26012016/A2744cl_ZPHOT.cat'
    ADdat      = np.genfromtxt(ADcat,names=True,dtype=None)

    ADcat_info = catdir + 'ASTRODEEP/fullrelease/A2744cl_26012016/A2744cl_A.cat'
    ADdat_info = np.genfromtxt(ADcat_info,names=True,dtype=None)

    glasszcat  = catdir + 'GLASSzcat/hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v002_redshiftcatalog.txt'
    glasszdat  = np.genfromtxt(glasszcat,names=True,dtype=None)
    glasszdat  = glasszdat[np.where(glasszdat['redshift'] > 0)]

    MUSEcat    = catdir + 'richard20/A2744_DR_v1.0.fits'
    MUSEdat    = afits.open(MUSEcat)[1].data

    if verbose: print(' - Loading GALFIT size estiamtes ')
    sizedir = '/Users/kschmidt/work/observing/proposals_preparation/180406_JWSTcycle1_A2744/galfitFromTaka/'
    size_info = np.genfromtxt(sizedir + 'f160w_01.cat',names=True,dtype=None,skip_header=27)
    size_dat  = np.genfromtxt(sizedir + 'galfit_crs_01.cat',names=True,dtype=None)

    size_id_info  = size_info['NUMBER']
    size_ra       = size_info['X_WORLD']
    size_dec      = size_info['Y_WORLD']

    size_id  = size_dat['1ID']
    size_r50 = size_dat['27r50_sext']
    size_r90 = size_dat['28r90_sext']

    # make sure to only consider IDs in info list that are also in GALFIT results catalog
    for ii, sii in enumerate(size_id_info):
        if sii not in size_id:
            size_id_info[ii]  = -99
            size_ra[ii]       = 0.0
            size_dec[ii]      = 0.0

    if verbose: print(' - Setting up output file '+outfile)
    if os.path.isfile(outfile) & (overwrite == False):
        sys.exit(outfile+' already exists and "overwrite"=False ')
    else:
        nowstr = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        fout   = open(outfile,'w')
        fout.write('# Best match JADES ids to GLASS objects detected in A2744. \n')
        fout.write('# Generated with grizli_wrappers.determine_JADESmatchForA2744obj() on '+nowstr+' \n')
        fout.write('#     \n')
        fout.write('#    -99      set if no match to id_GLASS in any of redshift catalogs searched within '+str(matchtol)+' arcsec \n')
        fout.write('#    -999     set if ASTRODEEP match has id > 100000 in which case there is no photo-z estimate \n')
        fout.write('#    -9999    set if no good match in JADES catalog given z, AB mag, and potentially M* and SFR from ASTRODEEP fits\n')
        fout.write('#    -99999   set if no GALFIT results for matched (GALFIT) object \n')
        fout.write('#     \n')
        fout.write('# Catalog can be loaded with:   dat = np.genfromtxt("'+outfile+'",names=True,dtype=None,skip_header=10) \n')
        fout.write('#     \n')
        fout.write('# id_GLASS ra_GLASS dec_GLASS f140wmag_GLASS    '
                   'cat_match id_match ra_match dec_match r_match redshift    '
                   'AD_id AD_ra AD_dec AD_r_match AD_zbest AD_Mstar AD_mu AD_SFR    '
                   'JADESid  JADESz JADESf140wmag JADESMstar JADESsfr     '
                   'size_id_info   size_ra_match  size_dec_match  size_r_match  size_id_galfit  size_r50obj_arcsec  size_r90_arcsec\n')

    if verbose: print(' - Loop over objects to get match to JADES mock catalog')
    for oo, objid in enumerate(id_GLASS):
        if testing:
            verbose_loop = True
            verbose_ju   = True
            if (objid < 450) or (objid > 470): continue
        else:
            verbose_loop = False
            verbose_ju   = False
        if verbose:
            infostr = '   matching and getting JADES info for id_GLASS = '+str(objid)+' ('+\
                      str("%6.f" % (oo+1))+' / '+str("%6.f" % len(id_GLASS))+')          '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        f140wmag  = np.abs(mag_GLASS[oo])
        obj_radec = SkyCoord(ra=ra_GLASS[oo]*u.degree, dec=dec_GLASS[oo]*u.degree)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose_loop: print(' - Crossmatch object to redshift catalogs ')
        AD_radec    = SkyCoord(ra=ADdat_info['RA']*u.degree, dec=ADdat_info['DEC']*u.degree)
        MUSE_radec  = SkyCoord(ra=MUSEdat['RA']*u.degree, dec=MUSEdat['DEC']*u.degree)
        glass_radec = SkyCoord(ra=glasszdat['RA']*u.degree, dec=glasszdat['DEC']*u.degree)
        size_radec = SkyCoord(ra=size_ra*u.degree, dec=size_dec*u.degree)

        if verbose_loop: print(' - Getting sources within '+str(matchtol)+' arc seconds of redshift catalogs ')
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
        if verbose_loop: print(' - Extract redshift for object (MUSE over GLASS over ASTRODEEP)')
        if len(MUSE_idx) > 0:
            bestent    = np.where(MUSE_d2d.arcsec == np.min(MUSE_d2d.arcsec))
            matchent   = MUSE_idx[bestent]
            cat_match  = 'MUSE'
            id_match   = MUSEdat['iden'][matchent]
            ra_match   = MUSEdat['RA'][matchent]
            dec_match  = MUSEdat['DEC'][matchent]
            r_match    = MUSE_d2d.arcsec[bestent]
            redshift   = MUSEdat['z'][matchent][0]
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
            if verbose_loop: print('   No match to id_GLASS = '+str(objid)+' in any of the redshift catalogs searched             ')
            cat_match  = 'NONE'
            id_match   = [-99]
            ra_match   = [-99]
            dec_match  = [-99]
            r_match    = [-99]
            redshift   = -99

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose_loop: print(' - Assign mSFR to object to be used for selection if good AD match available.')
        if len(AD_idx) > 0:
            bestent    = np.where(AD_d2d.arcsec == np.min(AD_d2d.arcsec))
            matchent   = AD_idx[bestent]

            AD_id      = ADdat_info['ID'][matchent]
            AD_ra      = ADdat_info['RA'][matchent]
            AD_dec     = ADdat_info['DEC'][matchent]
            AD_r_match = AD_d2d.arcsec[bestent]

            if ADdat_info['ID'][matchent] < 100000:
                if ADdat_info['ID'][matchent] != ADdat['ID'][matchent]:
                    sys.exit(' The IDs of the ASTRODEEP catalogs loaded do not match ')
                else:
                    AD_zbest     = ADdat['ZBEST'][matchent][0]
                    AD_mu        = ADdat['MAGNIF'][matchent][0]

                    AD_Mstar     = np.log10(ADdat['MSTAR'][matchent][0]*1e9)
                    AD_Mstar_min = np.log10(ADdat['MSTAR_MIN'][matchent][0]*1e9)
                    AD_Mstar_max = np.log10(ADdat['MSTAR_MAX'][matchent][0]*1e9)

                    AD_SFR       = np.log10(ADdat['SFR'][matchent][0])
                    AD_SFR_min   = np.log10(ADdat['SFR_MIN'][matchent][0])
                    AD_SFR_max   = np.log10(ADdat['SFR_MAX'][matchent][0])
            else:
                AD_zbest   = -99
                AD_mu      = -99
                AD_Mstar   = -99
                AD_SFR     = -99
        else:
                AD_id      = [-99]
                AD_ra      = [-99]
                AD_dec     = [-99]
                AD_r_match = [-99]
                AD_zbest   = -99
                AD_mu      = -99
                AD_Mstar   = -99
                AD_SFR     = -99

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose_loop: print(' - Select best JADES match based on z, M* and F140W mag AB')
        if redshift > 0:
            zdiff     = 0.1
            magdiff   = 0.5

            if (AD_SFR > 0) & (AD_Mstar > 0):
                if verbose_loop: print('   M* and SFR available so including those in selection before matching in redshift')
                zrange     = [redshift,-99] # [redshift-zdiff,redshift+zdiff]
                mrange     = [f140wmag-magdiff,f140wmag+magdiff]
                Mstarrange = [AD_Mstar_min,AD_Mstar_max]
                sfrrange   = [AD_SFR_min,AD_SFR_max]
            # elif (AD_SFR > 0):
            #     if verbose_loop: print('   SFR (not M*) available so including that in selection before matching in redshift')
            #     zrange     = [redshift,-99]
            #     mrange     = [f140wmag-magdiff,f140wmag+magdiff]
            #     Mstarrange = None
            #     sfrrange   = [AD_SFR_min,AD_SFR_max]
            # elif AD_Mstar > 0:
            #     if verbose_loop: print('   M* (not SFR) available so including that in selection before matching in redshift')
            #     zrange     = [redshift,-99] # [redshift-zdiff,redshift+zdiff]
            #     mrange     = [f140wmag-magdiff,f140wmag+magdiff]
            #     Mstarrange = [AD_Mstar_min,AD_Mstar_max]
            #     sfrrange   = None
            else:
                if verbose_loop: print('   Neither M* not SFR so selecting best match to mF140W given z+/-'+str(zdiff))
                JADESzMin = 0.2
                if redshift < JADESzMin:
                    zrange     = [0.0,JADESzMin+zdiff]
                else:
                    zrange     = [redshift-zdiff,redshift+zdiff]
                mrange     = [f140wmag,-99]
                Mstarrange = None
                sfrrange   = None

            JADESinfo  = ju.get_JADESobjects(redshift=zrange,mag_f140w=mrange,mStar=Mstarrange,SFR=sfrrange,
                                            jadesinfo=jadesinfo,verbose=verbose_ju)

            if len(JADESinfo) == 0:
                JADESid    = [-9999]
                JADESz     = [-9999]
                JADESmag   = [-9999]
                JADESmStar = [-9999]
                JADESsfr   = [-9999]
            else:
                JADESid    = JADESinfo['ID']
                JADESz     = JADESinfo['redshift']
                JADESmag   = -2.5*np.log10(JADESinfo['HST_F140W_fnu']/1e9)+8.90
                JADESmStar = JADESinfo['mStar']
                JADESsfr   = JADESinfo['SFR_100']

        else:
            zrange     = [0,100]
            mrange     = [f140wmag,-99]
            Mstarrange = None
            sfrrange   = None

            JADESinfo  = ju.get_JADESobjects(redshift=zrange,mag_f140w=mrange,mStar=Mstarrange,SFR=sfrrange,
                                            jadesinfo=jadesinfo,verbose=verbose_ju)

            JADESid    = JADESinfo['ID']
            JADESz     = JADESinfo['redshift']
            JADESmag   = -2.5*np.log10(JADESinfo['HST_F140W_fnu']/1e9)+8.90
            JADESmStar = JADESinfo['mStar']
            JADESsfr   = JADESinfo['SFR_100']
            cat_match  = 'JADES'
            id_match   = JADESid
            ra_match   = [-99]
            dec_match  = [-99]
            r_match    = [-99]
            redshift   = JADESz

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose_loop: print(' - Add size of nearest object (to GLASS coordinates) from Takahiros estimates ')
        size_bestent    = np.where(size_d2d.arcsec == np.min(size_d2d.arcsec))
        size_matchent   = np.atleast_1d(size_idx)[size_bestent]
        size_ent        = np.where(size_id == size_id_info[size_matchent])[0]
        size_id1        = size_id[size_ent]
        size_ra_match   = size_ra[size_ent][0]
        size_dec_match  = size_dec[size_ent][0]
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

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose_loop: print(' - Store object information and JADES id to output file')
        outstring = str(objid)+'  '+str(ra_GLASS[oo])+'  '+str(dec_GLASS[oo])+'  '+str("%.4f" % f140wmag)+'       '+\
                    str("%10s" % cat_match)+'  '+str(id_match[0])+'  '+str(ra_match[0])+' '+\
                    str(dec_match[0])+'  '+str(r_match[0])+'  '+str("%.4f" % redshift)+'       '+\
                    str(AD_id[0])+'  '+str(AD_ra[0])+'  '+str(AD_dec[0])+'  '+str(AD_r_match[0])+'  '+\
                    str(AD_zbest)+'  '+str(AD_Mstar)+'  '+str(AD_mu)+'  '+str(AD_SFR)+'       '+\
                    str(JADESid[0])+'  '+str(JADESz[0])+'  '+str(JADESmag[0])+'  '+str(JADESmStar[0])+'  '+str(JADESsfr[0])+'       '+\
                    str(size_id1[0])+'  '+str(size_ra_match)+'  '+str(size_dec_match)+'  '+\
                    str(size_r_match)+'  '+str(size_id2)+'  '+str(size_r50obj)+'  '+str(size_r90obj)
        fout.write(outstring+'\n')

        # print(str(redshift)+'\n')
        # if redshift == -99.: pdb.set_trace()
        # if id_match == 1718: pdb.set_trace()
    if verbose: print('\n   ... done')

    fout.close()
    if verbose: print(' - Stored JADES ID matches of all '+str(len(id_GLASS))+' objects to '+outfile)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_matchAD2GLASS(matchtol=0.5,GLASSids=None,ASTRODEEPids=None,verbose=True):
    """

    --- INPUT ---
    matchtol         Tolerance for mataches to return
    GLASSids         List of GLASS IDs; matches only returned if GLASS id is in list.
    ASTRODEEPids     List of ASTRODEEP IDs; matches only returned if ASTRODEEP id is in list.
    verbose          Toggle verbosity

    --- EXAMPLE OF USE ---

    XWids = [12, 62, 73, 103, 145, 189, 203, 394, 437, 446, 466, 561, 585, 742, 834, 855, 952, 967, 1081, 1238, 1456, 1718, 1722, 2007, 2054, 2066, 2255, 2425, 2434, 2438, 2441, 2447, 2452, 2465, 2471, 2473, 2484, 2538, 2547, 2557, 2560, 2561, 2595, 20004, 20005, 20009, 20108, 20143, 20214, 20257, 20295, 20327, 20409, 20471, 20577, 20656, 20674, 20728, 20759, 20782, 20785, 20808, 20953, 21023, 21225, 21226, 21323, 21430, 21534, 21565, 21582, 21616, 21636, 21700, 21797, 21801, 21812, 21843, 21912, 22011, 22027, 22035, 22055, 22057, 22064, 22066, 22107, 22125, 22154, 22185, 22240, 22256, 22289, 22357, 22542, 23034, 23095, 23105, 23153, 23235, 23293, 23294, 23298, 23310, 23312, 23320, 23325, 23327, 23343, 23349, 23424, 23448, 23515, 23521]

    matches = gw.get_matchAD2GLASS(matchtol=0.5,ASTRODEEPids=XWids,verbose=False)

    """
    if verbose: print(' - Get matches to ASTRODEEP catalog')
    cat = np.genfromtxt('hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v001_glassmaster.txt',dtype=None,names=True)

    ADpath  = '/Users/kschmidt/work/GLASS/LAEsearchFullGLASS/catalogs/ASTRODEEP/fullrelease/'
    AD_cat  = np.genfromtxt(ADpath+'A2744cl_26012016/A2744cl_A.cat',dtype=None,names=True)

    if verbose: print(' - Ignore ASTRODEEP objects with MAG_JH140 >= 99.0')
    AD_cat  = AD_cat[AD_cat['MAG_JH140'] < 99.0]

    AD_radec  = SkyCoord(ra=AD_cat['RA']*u.degree, dec=AD_cat['DEC']*u.degree)
    cat_radec = SkyCoord(ra=cat['X_WORLD']*u.degree, dec=cat['Y_WORLD']*u.degree)

    if verbose: print(' - Getting sources within the match toleracnce of '+str(matchtol)+' arc seconds')
    AD_idx, d2d, d3d = cat_radec.match_to_catalog_sky(AD_radec)

    if verbose: print(' - Printing ID matches:')
    if verbose: print('   id_GLASS   id_ASTRODEEP   r_match_arcsec ')
    outarr = np.array([])
    for ii, id_GLASS in enumerate(cat['NUMBER']):
        if d2d[ii] < matchtol*u.arcsec:
            id_AD = AD_cat['ID'][AD_idx[ii]]

            if ASTRODEEPids is not None:
                if id_AD not in ASTRODEEPids: continue

            if GLASSids is not None:
                if id_GLASS not in GLASSids: continue

            if verbose: print(str("%10s" % id_GLASS)+' '+str("%10s" % id_AD)+
                              '        '+str("%.6f" % d2d[ii].arcsec))
            if len(outarr) == 0:
                outarr = np.array([int(id_GLASS),int(id_AD),float(d2d[ii].arcsec)])
            else:
                outarr = np.vstack([ outarr,np.array([int(id_GLASS),int(id_AD),float(d2d[ii].arcsec)]) ])

    return outarr

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def print_JADESoutputFromFile(GLASSids,verbose=True,
                              JADESmatches='/Users/kschmidt/work/JWST/grizly_A2744/Sim_A2744_NIRCAM/A2744_JADESmatches_201021.txt'):
    """
    GLASSids = [3.0,25.0,34.0,35.0,98.0,104.0,118.0,142.0,463.0,477.0,783.0,1014.0,1084.0,1517.0,1535.0,1627.0,1722.0,1744.0,1745.0,1787.0,1977.0,1987.0,2095.0,2113.0]

    gw.print_JADESoutputFromFile(GLASSids)

    """
    dat = np.genfromtxt(JADESmatches,names=True,dtype=None,skip_header=10)

    JADESids = []
    objents  = []
    for id in GLASSids:
        objent = np.where(dat['id_GLASS'] == float(id))[0]
        if len(objent) == 1:
            if verbose: print( str(dat[objent][0])[1:-1].replace(',',' ') )
            JADESids.append(dat['JADESid'][objent][0])
            objents.append(objent[0])
        else:
            if verbose: print(' WARNING: Found '++' matches to input GLASS id '+str(id)+' so skipping ')
            JADESids.append(np.nan)

    return JADESids, dat[np.asarray(objents)]
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def create_diffimg(fits1,ext1,fits2,ext2,outfile,overwrite=True,header=False,verbose=True):
    """

    --- EXAMPLE OF USE ---
    fitsfile = './nircam_F356W_00.01.GrismFLT.fits'
    outfile  = fitsfile.replace('GrismFLT.fits','GrismFLT_modeldiff.fits')
    gw.create_diffimg(fitsfile,'GSCI',fitsfile,'MODEL',outfile,overwrite=True,header=True)

    """
    img1 = afits.open(fits1)[ext1].data
    img2 = afits.open(fits2)[ext2].data

    if header:
        hdr = afits.open(fits1)[ext1].header
    else:
        hdr = None

    diffimage = img1-img2

    hdu = afits.PrimaryHDU(diffimage, header=hdr)
    hdu.writeto(outfile, overwrite=overwrite)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def check_coeffs(polycoeffs=[0.05, 0.01],polyxrange=[1.0, 7.0],verbose=True):
    """
    Function to generate and plot model polynomials

    --- EXAMPLE OF USE ---
    gw.check_coeffs(polycoeffs=[0.05, 0.00])

    """
    xspec_m1 = np.arange(polyxrange[0], polyxrange[1], 0.05)-1
    yspec_co = [xspec_m1**o*polycoeffs[o] for o in range(len(polycoeffs))]
    xspec    = (xspec_m1+1)*1.e4
    yspec    = np.sum(yspec_co, axis=0)

    plotname = './check_coeffs_plot.pdf'
    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.2, right=0.97, bottom=0.10, top=0.9)
    Fsize    = 10
    lthick   = 2
    marksize = 4
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    plt.title('polycoeffs: '+str(polycoeffs)+'    polyxrange: '+str(polyxrange))

    plt.plot(xspec,yspec)

    plt.xlabel(' Wavelength [A] ')
    plt.ylabel(' ')
    print(' - Saved plot to '+plotname)

    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def tracedisplacement(DISPL,DISPX,DISPY,waverange,verbose=True):
    """
    Tranlating Pirzkal & Ryan (2016) deispersion polynomial parameters to the grizli format.

    DISPL,DISPX, and DISPY dispersion polynomial parameters can all be taken from
    https://github.com/npirzkal/GRISM_NIRCAM/tree/master/V2

    --- INPUT ---
    DISPL0,DISPL1       The wavelength (lambda) DISPL_#_0 and DISPL_#_1 (# indicates spectral order) parameters
                        from *.conf files at https://github.com/npirzkal/GRISM_NIRCAM/tree/master/V2

    DISPX1,DISPX0       The x-direction DISPX_#_0 and DISPX_#_1 (# indicates spectral order) parameters
                        from *.conf files at https://github.com/npirzkal/GRISM_NIRCAM/tree/master/V2

    DISPY1,DISPY0       The y-direction DISPY_#_0 and DISPY_#_1 (# indicates spectral order) parameters
                        from *.conf files at https://github.com/npirzkal/GRISM_NIRCAM/tree/master/V2

    waverange           Wavelength range of dispersion filter to consider.
                        Can be obtained from senistivity curves at https://github.com/npirzkal/GRISM_NIRCAM/tree/master/V2

    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw
    gw.tracedisplacement([24000.,26000.],[-1530.8764939967652,2589.641434263754],[0,0],[22794.3334,32954.50519]) # OLD
    gw.tracedisplacement([23699.9735,8700.0146],[-1530.8764939967652,866.5338650000007],[19.97158346116109,1.5926025904834888],[23700.0,32400.0]) #v2

    """
    DLDP_1 = DISPL[1] / DISPX[1]
    DLDP_0 = DISPL[0] - DISPX[0] * DLDP_1

    pix_limit_min = (waverange[0] - DLDP_0) / DLDP_1
    pix_limit_max = (waverange[1] - DLDP_0) / DLDP_1

    if not np.array_equal(np.asarray(DISPY).astype(float), np.zeros(2)):
        if verbose: print(' WARNING: Ignoring provided dispersions in the y-direction')

    pixlimits = [int(pix_limit_min), int(pix_limit_max)]
    DLDP      = [DLDP_0, DLDP_1]

    return pixlimits, DLDP
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_grizli_dispersion_conf(outputdir,axeconfdir,filterlist=['F277W','F356W','F444W'],
                               grisms=['GRISMR','GRISMC'],modules=['A','B'],verbose=True):
    """
    Generating Tranlating Pirzkal & Ryan (2016) deispersion polynomial parameters to the grizli format.

    --- EXAMPLE OF USE ---
    import grizli_wrappers as gw

    outputdir   = '/Users/kschmidt/work/grizli/CONF/NIRCamConfPirzkal_grizliversion/'
    axeconfdir  = '/Users/kschmidt/work/grizli/CONF/NIRCamConfPirzkal/V2/'
    gw.gen_grizli_dispersion_conf(outputdir,axeconfdir,filterlist=['F277W','F356W','F444W'])

    """
    grizliorders = {'+1':'A','+2':'B'}


    if filterlist == 'all':
        filters = ['F430M', 'F460M', 'F250M', 'F335M', 'F322W2', 'F277W', 'F356W', 'F444W', 'F300M', 'F480M', 'F410M', 'F360M']
    else:
        filters = filterlist

    for filter in filters:
        for module in modules:
            for grism in grisms:
                confinput = axeconfdir+'NIRCAM_'+filter+'_mod'+module+'_'+grism[-1]+'.conf'
                if verbose: print( ' - Loading configuration in '+confinput)
                axeconf = grizli.grismconf.aXeConf(confinput)

                outconf   = outputdir+'NIRCam.'+module+'.'+grism+'.'+filter+'.conf'
                if verbose: print( '   Storing grizli-friendly configuration setup to \n   '+outconf)
                fout = open(outconf, 'w')

                for order in ['+1', '+2']:
                    outsens                = outputdir+'NIRCam.'+module+'.'+grism+'.'+filter+'.'+grizliorders[order]+'.sens.fits'
                    if verbose: print( '   Storing filter sensitivity curve to \n   '+outsens)
                    filtercurve            = afits.open(axeconfdir+axeconf.conf['SENSITIVITY_'+order])[1].data
                    outtab                 = grizli.utils.GTable()
                    outtab['WAVELENGTH']   = filtercurve['WAVELENGTH']*1.e4
                    outtab['SENSITIVITY']  = filtercurve['SENSITIVITY']
                    outtab['ERROR']        = filtercurve['ERROR']
                    outtab.write(outsens, format='fits', overwrite=True)
                    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    if verbose: print( '   Estimating dispersion parameters with gw.tracedisplacement()')
                    DISPL     = [axeconf.conf['DISPL_'+order+'_0']*1.e4,axeconf.conf['DISPL_'+order+'_1']*1.e4]
                    DISPX     = [axeconf.conf['DISPX_'+order+'_0'],axeconf.conf['DISPX_'+order+'_1']]
                    DISPY     = [axeconf.conf['DISPY_'+order+'_0'],axeconf.conf['DISPY_'+order+'_1']]

                    selectent = np.where(outtab['SENSITIVITY'] > 0.001*np.max(outtab['SENSITIVITY']))
                    lamrange  = [np.min(outtab['WAVELENGTH'][selectent]),
                                 np.max(outtab['WAVELENGTH'][selectent])]

                    pixlimits, DLDP = gw.tracedisplacement(DISPL,DISPX,DISPY,lamrange,verbose=True)

                    fout.write("""
#------------------------------
# Beam/order {o}={o_axe}

MMAG_EXTRACT_{o} 40

BEAM{o} {left} {right}

XOFF_{o} 0.0

YOFF_{o} 0.0

DYDX_ORDER_{o} 0

DYDX_{o}_0 0.

DISP_ORDER_{o} 1

DLDP_{o}_0 {DLDP_0:.2f}

DLDP_{o}_1 {DLDP_1:.2f}

SENSITIVITY_{o} {grismsens}

""".format(o=grizliorders[order], module=module, DLDP_0=DLDP[0], DLDP_1=DLDP[1], filt=filter, o_axe=order,
           left=pixlimits[0], right=pixlimits[1], grismsens=outsens.split('/Users/kschmidt/work/grizli/CONF/')[-1]))

                fout.close()
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =