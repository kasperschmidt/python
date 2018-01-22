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
#importlib.reload(module)
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


        print(" - For analyzing fit output see\n""
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
    import grizli
    import grizli.utils

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