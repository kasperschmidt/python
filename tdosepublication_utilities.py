# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import numpy as np
import pdb
import astropy.io.fits as afits
import os
import sys
import glob
import aplpy
import shutil
import MiGs
import tdose_utilities as tu
import MUSEWideUtilities as mu
import tdosepublication_utilities as tsu
import tdose_model_FoV as tmf
import kbsutilities as kbs
import matplotlib.pyplot as plt
import tdose_extract_spectra as tes
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_8685():
    """

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    tsu.plot_8685()

    """
    parentdir = '/Users/kschmidt/work/publications/TDOSE/TDOSEexampleruns/'
    specs = [parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_aperture_0000008685-0000008685.fits',
             parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_gauss_0000008685-0000008685.fits',
             parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_gauss_0008685000-0000008685.fits',
             parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_modelimg_0008685000-0000008685.fits',
             parentdir+'MWDR1_guo8685_mw102033149/aper_spectrum_candels-cdfs-02_102033149.fits',
             parentdir+'MWDR1_guo8685_mw102033149/emission_spectrum_candels-cdfs-02_102033149.fits',
             parentdir+'MWDR1_guo8685_mw102033149/tdose_spectrum_candels-cdfs-02_08685.fits']


    filelist    = [specs[2]]
    labels      = ['6-component Gauss model']

    compspec    = [specs[0],specs[1]]#,specs[3]]
    comp_labels = ['Aperture', '1-component Gauss model']#,'6-component Sersic model']
    comp_colors = ['blue','red']#,'green']

    plotname = parentdir+'/tdose_1Dspectra_comparison_8685.pdf'
    xrange   = [6050,6290]
    yrange   = [-100,650]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=False,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=True)


    plotname = parentdir+'/tdose_1Dspectra_comparison_8685_S2N.pdf'
    yrange   = [0,13]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=True,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=True)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    filelist    = [specs[0],specs[1],specs[2],specs[3],specs[6]]
    labels      = ['Aperture', '1-component Gauss model','6-component Gauss model','6-component Sersic model', 'MWDR1 TDOSE']
    colors      = ['blue','red','black','green','orange']

    compspec    = [specs[4],specs[5]]
    comp_labels = ['MWDR1 aperture', 'MWDR1 PSF weighted']
    comp_colors = ['magenta','cyan']

    comp_wavecol= 'WAVE_AIR'
    comp_fluxcol= 'FLUX'
    comp_errcol = 'FLUXERR'

    xrange   = [6160,6220]
    yrange   = [0,15]

    plotname = parentdir+'/tdose_1Dspectra_comparison_8685_S2N_AllwMWDR1specs.pdf'
    tes.plot_1Dspecs(filelist,plotname=plotname,colors=colors,labels=labels,plotSNcurve=True,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol=comp_wavecol,comp_fluxcol=comp_fluxcol,comp_errcol=comp_errcol,
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=False)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_9262():
    """

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    tsu.plot_9262()

    """
    parentdir = '/Users/kschmidt/work/publications/TDOSE/TDOSEexampleruns/'
    specs     = [parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_aperture_0000009262-0000009262.fits',
                 parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_gauss_0000009262-0000009262.fits',
                 parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_gauss_0000009262-0000009262.fits',
                 parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_modelimg_0000009262-0000009262.fits',
                 parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_gauss_0000009262-0000009262_bothcompmodelimg.fits',
                 parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_modelimg_0000009262-0000009262_bothcompmodelimg.fits',
                 parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_modelimg_0000009262-0000009262_oppcompmodelimg.fits']

    # LAE is oppcompmodelimg spectrum, i.e., specs[6] whihc corresponds to a source [0] extraction according to spectrum source cube header. Comparing that with the region file for the model reveles the LAE.

    filelist    = [specs[1]]
    labels      = ['Single Gauss model']

    compspec    = [specs[6],specs[3],specs[5]]#,specs[3]]
    comp_labels = ['Sersic source 1', 'Sersic source 2', 'Sersic both sources']
    comp_colors = ['red', 'blue', 'green']

    plotname = parentdir+'/tdose_1Dspectra_comparison_9262manually.pdf'
    xrange   = [5030,5080]
    yrange   = [-200,500]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=False,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=True)

    plotname = plotname.replace('.pdf','_S2N.pdf')
    yrange   = [-1,8]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=True,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_9640_MWDR1():
    """

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    tsu.plot_9640_MWDR1()

    """

    parentdir = '/Users/kschmidt/work/publications/TDOSE/MUSEWideDR1data/Guo9640/'

    specs     = [parentdir+'aper_spectrum_candels-cdfs-02_102009072.fits',
                 parentdir+'emission_spectrum_candels-cdfs-02_102009072.fits',
                 parentdir+'tdose_spectrum_candels-cdfs-02_09640.fits']

    filelist    = [specs[2]]
    labels      = ['TDOSE']

    compspec    = [specs[0],specs[1]]
    comp_labels = ['Aperture', 'PSF weighted']
    comp_colors = ['blue', 'red']

    plotname = parentdir+'/tdose_1Dspectra_comparison_9240_MWDR1_manually.pdf'
    xrange   = [6350,6740]
    yrange   = [-100,2000]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=False,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='WAVE_AIR',comp_fluxcol='FLUX',comp_errcol='FLUXERR',
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=True)

    plotname = plotname.replace('.pdf','_S2N.pdf')
    yrange   = [-2,35]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=True,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='WAVE_AIR',comp_fluxcol='FLUX',comp_errcol='FLUXERR',
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=True)


    plotname = parentdir+'/tdose_1Dspectra_comparison_fullrange_9240_MWDR1_manually.pdf'
    xrange   = [4900,9200]
    yrange   = [-100,2400]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=False,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='WAVE_AIR',comp_fluxcol='FLUX',comp_errcol='FLUXERR',
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=True)

    plotname = plotname.replace('.pdf','_S2N.pdf')
    yrange   = [-2,35]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=True,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='WAVE_AIR',comp_fluxcol='FLUX',comp_errcol='FLUXERR',
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=True)



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_UDF03_QSO(pubversion=False):
    """

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    tsu.plot_UDF03_QSO(pubversion=False)

    """

    parentdir = '/Users/kschmidt/work/MUSE/UDFquasar/tdose_extraction/tdose_modelbased_extractions/'

    specs     = [parentdir+'190201_DLAextraction_component4/tdose_spectrum_modelimg_0000000005-0000000005_DLAextraction_component4.fits',
                 parentdir+'190201_QSOextraction_component1and5/tdose_spectrum_modelimg_0000000005-0000000005_190201_QSOextraction_component1and5.fits',
                 parentdir+'190201_QSOextraction_component1only/tdose_spectrum_modelimg_0000000005-0000000005_190201_QSOextraction_component1only.fits']

    filelist    = [specs[0]]
    labels      = ['DLA?']

    compspec    = [specs[1]]#,specs[2]]
    comp_labels = ['QSO 2 components']#, 'QSO 1 component']
    comp_colors = ['blue']#, 'red']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MiGlines      = MiGs.linelistdic(listversion='full') # loading line list for plots
    line_waverest = np.asarray([MiGlines[key][1] for key in MiGlines.keys()])
    line_names    = np.asarray([MiGlines[key][0] for key in MiGlines.keys()])
    objredshift   = 3.187
    linelist      = [(1+objredshift)*line_waverest,line_names]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = parentdir+'/tdose_1Dspectra_UDF3_QSOandDLA_zoom1.pdf'
    xrange   = [6350,7000]
    yrange   = [-100,300]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=False,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',showlinelist=linelist,
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=pubversion)

    plotname = plotname.replace('.pdf','_S2N.pdf')
    yrange   = [-5,30]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=True,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',showlinelist=linelist,
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=pubversion)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = parentdir+'/tdose_1Dspectra_UDF3_QSOandDLA_zoom2.pdf'
    xrange   = [5050,5225]
    yrange   = [-100,150]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=False,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',showlinelist=linelist,
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=pubversion)

    plotname = plotname.replace('.pdf','_S2N.pdf')
    yrange   = [-5,10]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=True,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',showlinelist=linelist,
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=pubversion)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = parentdir+'/tdose_1Dspectra_UDF3_QSOandDLA_zoom3.pdf'
    xrange   = [7800,8050]
    yrange   = [-100,150]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=False,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',showlinelist=linelist,
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=pubversion)

    plotname = plotname.replace('.pdf','_S2N.pdf')
    yrange   = [-5,10]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=True,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',showlinelist=linelist,
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=pubversion)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = parentdir+'/tdose_1Dspectra_UDF3_QSOandDLA_fullrange.pdf'
    xrange   = [4900,9200]
    yrange   = [-100,400]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=False,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',showlinelist=linelist,
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=pubversion)

    plotname = plotname.replace('.pdf','_S2N.pdf')
    yrange   = [-5,30]

    tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=True,
                     comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                     comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',showlinelist=linelist,
                     xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=pubversion)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_specs():
    """

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    tsu.plot_specs()

    """
    parentdir  = '/Users/kschmidt/work/publications/TDOSE/TDOSEexampleruns/'

    specs9640  = [parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_aperture_0000009640-0000009640.fits',
                  parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_gauss_0000009640-0000009640.fits',
                  parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_gauss_0000009640-0000009640.fits',
                  parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_modelimg_0000009640-0000009640.fits']
    #ranges9640 = [[4800,9000],[-100,1200],[0,32]] # full range
    #ranges9640 = [[6400,6750],[-100,1200],[0,32]]  # [OIII]-doublet
    ranges9640 = [[4900,5030],[-100,1200],[-1,32]]  # [OII]-doublet
    #ranges9640 = [[8500,9000],[-100,1200],[0,32]]  # Ha @ ~8779
    pname9640  = parentdir+'/tdose_1Dspectra_comparison_9640.pdf'

    specs9262  = [parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_aperture_0000009262-0000009262.fits',
                  parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_gauss_0000009262-0000009262.fits',
                  parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_gauss_0000009262-0000009262.fits',
                  parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_modelimg_0000009262-0000009262.fits']
    ranges9262 = [[4950,5100],[-100,650],[-1,10]]
    pname9262  = parentdir+'/tdose_1Dspectra_comparison_9262.pdf'    #

    specs9093  = [parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_aperture_0000009093-0000009093.fits',
                  parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_gauss_0000009093-0000009093.fits',
                  parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_gauss_0000009093-0000009093.fits',
                  parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_modelimg_0000009093-0000009093.fits']
    ranges9093 = [[5200,5350],[-100,650],[-1,10]]
    pname9093  = parentdir+'/tdose_1Dspectra_comparison_9093.pdf'

    specs0102013086  = [parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_aperture_0102013086-0102013086.fits',
                        parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_gauss_0102013086-0102013086.fits',
                        parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_gauss_0102013086-0102013086.fits',
                        parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_modelimg_0102013086-0102013086.fits']
    ranges0102013086 = [[5050,5250],[-100,650],[-1,10]]
    pname0102013086  = parentdir+'/tdose_1Dspectra_comparison_0102013086.pdf'

    specs0102014087  = [parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_aperture_0102014087-0102014087.fits',
                        parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_gauss_0102014087-0102014087.fits',
                        parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_gauss_0102014087-0102014087.fits',
                        parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_modelimg_0102014087-0102014087.fits']
    ranges0102014087 = [[5050,5250],[-100,650],[-1,10]]
    pname0102014087  = parentdir+'/tdose_1Dspectra_comparison_0102014087.pdf'

    specs0102049176  = [parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_aperture_0102049176-0102049176.fits',
                        parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_gauss_0102049176-0102049176.fits',
                        parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_gauss_0102049176-0102049176.fits',
                        parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_modelimg_0102049176-0102049176.fits']
    ranges0102049176 = [[7350,7450],[-100,650],[-1,10]]
    pname0102049176  = parentdir+'/tdose_1Dspectra_comparison_0102049176.pdf'

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    specs9262bothcomp  = [parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_aperture_0000009262-0000009262.fits',
                          parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_gauss_0000009262-0000009262.fits',
                          parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_gauss_0000009262-0000009262_bothcompmodelimg.fits',
                          parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_modelimg_0000009262-0000009262_bothcompmodelimg.fits']
    ranges9262bothcomp = [[4950,5100],[-100,650],[-1,10]]
    pname9262bothcomp  = parentdir+'/tdose_1Dspectra_comparison_9262_bothcompmodelimg.pdf'


    specs9262oppcomp  = [parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_aperture_0000009262-0000009262.fits',
                          parentdir+'181008_gaussaperture/tdose_spectra/tdose_spectrum_gauss_0000009262-0000009262.fits',
                          parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_gauss_0000009262-0000009262_oppcompmodelimg.fits',
                          parentdir+'181008_modelimg/tdose_spectra/tdose_spectrum_modelimg_0000009262-0000009262_oppcompmodelimg.fits']
    ranges9262oppcomp = [[4950,5100],[-100,650],[-1,10]]
    pname9262oppcomp  = parentdir+'/tdose_1Dspectra_comparison_9262_oppcompmodelimg.pdf'

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    specvecs  = [specs9640,specs9262,specs9093,specs0102013086,specs0102014087,specs0102049176,
                 specs9262bothcomp,specs9262oppcomp]
    rangevecs = [ranges9640,ranges9262,ranges9093,ranges0102013086,ranges0102014087,ranges0102049176,
                 ranges9262bothcomp,ranges9262oppcomp]
    plotnames = [pname9640,pname9262,pname9093,pname0102013086,pname0102014087,pname0102049176,
                 pname9262bothcomp,pname9262oppcomp]
    for ss, specs in enumerate(specvecs):
        filelist    = [specs[3]]
        labels      = ['Modelimg']

        compspec    = [specs[0],specs[1],specs[2]]
        comp_labels = ['Aperture', 'Gauss', 'Gauss-modelimg']
        comp_colors = ['blue','red','green']

        plotname = plotnames[ss]
        xrange   = rangevecs[ss][0] # [6050,6290]
        yrange   = rangevecs[ss][1] # [-100,650]

        tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=False,
                         comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                         comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',
                         xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=True)


        plotname = plotname.replace('.pdf','_S2N.pdf')
        yrange   = rangevecs[ss][2] # [0,13]

        tes.plot_1Dspecs(filelist,plotname=plotname,colors=['black'],labels=labels,plotSNcurve=True,
                         comparisonspecs=compspec,comp_colors=comp_colors,comp_labels=comp_labels,
                         comp_wavecol='wave',comp_fluxcol='flux',comp_errcol='fluxerror',
                         xrange=xrange,yrange=yrange,showspecs=False,shownoise=True,verbose=True,pubversion=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def extract_spectrum(sourceIDs,layer_scale_arr,wavelengths,noise_cube=None,source_model_cube=None,
                     specname='tdose_extract_spectra_extractedspec.fits',obj_cube_hdr=None,data_cube=None,
                     clobber=False,verbose=True):
    return None

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def getOIIemitters(magcut=25,sepcut=0.3,verbose=True,
                   savefits='/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/OIIemitter_selection',
                   plotnamebase='/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/OIIemitter_selectionplot'):
    """

    Get sample of OII emitters to build GALFIT and TDOSE-Guass models of for comparison purposes.

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    MWsubcat, mags814 = tsu.getOIIemitters(magcut=25,sepcut=0.3)

    """
    if type(magcut) is not list:
        magcut = [0,magcut]

    DR1cat_el     = '/Users/kschmidt/work/catalogs/MUSE_GTO/MW_44fields_emline_table_v1.0.fits'
    DR1cat_main   = '/Users/kschmidt/work/catalogs/MUSE_GTO/MW_44fields_main_table_v1.0.fits'
    DR1dat_main   = afits.open(DR1cat_main)[1].data

    SkeltonCat    = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
    SkeltonDat    = afits.open(SkeltonCat)[1].data

    OIIent = np.where((DR1dat_main['LEAD_LINE'] == 'O2') & (DR1dat_main['CONFIDENCE'] == 3))[0]

    print(' - Found '+str(len(OIIent))+' object with OII being the leading line and a MW confidence of 3 ')

    ids_MW      = DR1dat_main['UNIQUE_ID'][OIIent]
    ids_Skelton = DR1dat_main['SKELTON_ID'][OIIent]
    sep_Skelton = DR1dat_main['SKELTON_SEP'][OIIent]

    goodent     = np.array([])
    mags814     = np.array([])
    for ii, id_MW in enumerate(ids_MW):
        if ids_Skelton[ii] > 0:
            ent_skelton = np.where(SkeltonDat['id'] == ids_Skelton[ii])[0]
            if len(ent_skelton) != 1:
                print('STOP! Number of matches in Skelton catalog to idMW='+str(id_MW)+' is not 1, but '+str(len(ent_skelton)))
                pdb.set_trace()
            else:
                f814 = SkeltonDat['f_f814wcand'][ent_skelton]
                m814 = 25.0-2.5*np.log10(f814)


                if (m814 > magcut[0]) & (m814 < magcut[1]) & (sep_Skelton[ii] < sepcut):
                    goodent = np.append(goodent,OIIent[ii])
                    mags814 = np.append(mags814,m814)

    print('   '+str(len(goodent))+' of these objects have a Skelton match within '+str(sepcut)+' arc seconds and '+str(magcut[0])+' < mag(f814w) < '+str(magcut[1]))
    goodsubcat = DR1dat_main[goodent.astype(int)]

    if savefits is not None:
        namestr = ('_'+str(magcut[0])+'LTm814LT'+str(magcut[1])+
                   '_SkelsepLT'+str(sepcut)+'_Nobj'+str(len(goodsubcat))).replace('.','p')+'.fits'
        hdu = afits.BinTableHDU.from_columns(goodsubcat._get_raw_data())
        hdu.writeto(savefits+namestr,overwrite=True)
        if verbose: print(' - Saved selection to '+savefits+namestr)

    settitle = True
    if settitle:
        titlestring = 'C=3 OII emitters w. '+str(magcut[0])+' $<$ mag(f814w) $<$ '+str(magcut[1])+\
                      ' and Skel\_sep$<$'+str(sepcut)+' (Nobj='+str(len(goodsubcat))+')'
        titlefsize  = 9

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = plotnamebase+'_selectionplot_zVSs2n.pdf'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=(5, 4))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.15, top=0.90)
    Fsize  = 12
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    if settitle: plt.title(titlestring,fontsize=titlefsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.plot(goodsubcat['Z'],goodsubcat['SN'],'ob',lw=lthick)
    plt.plot([np.min(goodsubcat['Z']),np.max(goodsubcat['Z'])],[5.0,5.0],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('Redshift', fontsize=Fsize)
    plt.ylabel('S/N of [OII] LSDCat detection', fontsize=Fsize)
    plt.yscale('log')
    if verbose: print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = plotnamebase+'_selectionplot_magVSs2n.pdf'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=(5, 4))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.15, top=0.90)
    Fsize  = 12
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    if settitle: plt.title(titlestring,fontsize=titlefsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.plot(mags814,goodsubcat['SN'],'ob',lw=lthick)
    plt.plot([np.min(mags814),np.max(mags814)],[5.0,5.0],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('magAB(814) from Skelton et al. (2014)', fontsize=Fsize)
    plt.ylabel('S/N of [OII] LSDCat detection', fontsize=Fsize)
    plt.yscale('log')
    if verbose: print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = plotnamebase+'_selectionplot_zVSsep.pdf'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=(5, 4))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.15, top=0.90)
    Fsize  = 12
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    if settitle: plt.title(titlestring,fontsize=titlefsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.plot(goodsubcat['Z'],goodsubcat['SKELTON_SEP'],'ob',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('Redshift', fontsize=Fsize)
    plt.ylabel('Separation to Skelton ID', fontsize=Fsize)
    if verbose: print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    return goodsubcat, mags814
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def extractGOODSS_cutoouts(dataarray,mwidcol='UNIQUE_ID',cutsize=5.0,racol='RA',deccol='DEC',clobber=True,
                           generatepng=False,outdir='/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/cutoutsTEMP/'):
    """
    Extract cutouts of a set of objects from an input data array
    Build to make cutouts of selection from tsu.getOIIemitters()

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    MWsubcat, mags814 = tsu.getOIIemitters(magcut=25,sepcut=0.3)
    tsu.extractGOODSS_cutoouts(MWsubcat,cutsize=10.0,generatepng=True)
    """
    DR1cat_main   = '/Users/kschmidt/work/catalogs/MUSE_GTO/MW_44fields_main_table_v1.0.fits'
    DR1dat_main   = afits.open(DR1cat_main)[1].data

    SkeltonCat    = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
    SkeltonDat    = afits.open(SkeltonCat)[1].data

    for ii, mwid in enumerate(dataarray[mwidcol]):
        # if ii > 2: continue
        fieldno    = str(mwid)[1:3]
        imgfile    = '/Users/kschmidt/work/images_MAST/MUSEWidePointings/acs_814w_candels-cdfs-'+fieldno+'_cut_v1.0.fits'
        ra         = dataarray[racol][ii]
        dec        = dataarray[deccol][ii]
        cutoutsize = [cutsize,cutsize]
        outname    = (outdir+'acs_814w_candels-cdfs-'+fieldno+'_cut_v1.0_MUSEWide'+\
                     str(int(mwid))+'_'+str(cutoutsize[0])+'X'+str(cutoutsize[1])+'arcseccut').replace('.','p')+'.fits'
        cutout     = tu.extract_subimage(imgfile,ra,dec,cutoutsize,outname=outname,clobber=clobber)

        if generatepng:
            fig = aplpy.FITSFigure(outname, figsize=(5, 5))
            fig.show_colorscale()
            fig.show_grayscale()
            fig.add_scalebar(1.0/3600.,color='white',label='1 arcsec')

            fig.add_label(ra, dec,mwid,color='red')
            fig.show_circles(ra, dec, 1.0/3600.*0.2,color='red')

            goodMW = np.where( (DR1dat_main['RA'] < ra + 1.0/3600. * cutsize) &
                               (DR1dat_main['RA'] > ra - 1.0/3600. * cutsize) &
                               (DR1dat_main['DEC'] < dec + 1.0/3600. * cutsize) &
                               (DR1dat_main['DEC'] > dec - 1.0/3600. * cutsize) )[0]
            for oo in goodMW:
                if DR1dat_main['UNIQUE_ID'][oo].astype(str) != str(mwid):
                    fig.add_label(DR1dat_main['RA'][oo], DR1dat_main['DEC'][oo],DR1dat_main['UNIQUE_ID'][oo].astype(str),color='white')
                    fig.show_circles(DR1dat_main['RA'][oo], DR1dat_main['DEC'][oo], 1.0/3600.*0.1,color='white')

            goodSK = np.where( (SkeltonDat['ra'] < ra + 1.0/3600. * cutsize) &
                               (SkeltonDat['ra'] > ra - 1.0/3600. * cutsize) &
                               (SkeltonDat['dec'] < dec + 1.0/3600. * cutsize) &
                               (SkeltonDat['dec'] > dec - 1.0/3600. * cutsize) )[0]
            for oo in goodSK:
                #fig.add_label(SkeltonDat['ra'][oo], SkeltonDat['dec'][oo],SkeltonDat['id'][oo].astype(str),color='green')
                fig.show_circles(SkeltonDat['ra'][oo], SkeltonDat['dec'][oo], 1.0/3600.*0.1,color='green')

            fig.tick_labels.set_font(size='large')
            fig.add_grid()
            fig.save(outname.replace('.fits','.png'))
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def OIIemitters_CopyAndRenameModels():
    """
    copy and rename galfit model imgblocks of OII emitters

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    tsu.OIIemitters_CopyAndRenameModels()

    """
    galfitresultsdir    = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/galfit_wrapper_results_final_181015/'
    models              = glob.glob(galfitresultsdir+'imgblocks/imgblock_acs_814w_*.fits')

    outdir              = '/Volumes/DATABCKUP1/TDOSEextractions/181016_MWDR1_OIIemitters/galfit_models/'
    for model in models:
        id       = model.split('814w_')[-1].split('.fit')[0]
        field    = id[1:3]
        newfile  = outdir+'model_acs_814w_candels-cdfs-CC_cut_v1.0_idIII_cutout4p0x4p0arcsec.fits'.replace('III',id).replace('CC',field)
        shutil.copy(model, newfile)
        tu.galfit_model_ds9region([newfile],clobber=True)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def OIIemitters_ConvertGALFITmodels2Cubes():
    """

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    tsu.OIIemitters_ConvertGALFITmodels2Cubes()

    """
    modeldir = '/Volumes/DATABCKUP1/TDOSEextractions/181016_MWDR1_OIIemitters/galfit_models/'
    models   = glob.glob(modeldir+'model_acs_814w_candels-cdfs-*_cut_v1.0_id*_cutout4p0x4p0arcsec.fits')
    Nmodels  = len(models)

    componentinfo = '/Users/kschmidt/work/TDOSE/181016_MWDR1_OIIemitters/OIIemitter_componentinfo_edited181016.txt'

    PSFmodel = afits.open('/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/PSF_models/acs_814w/imgblock_6475.fits')[2].data
    convkernel = np.zeros([133,133])
    convkernel[6:127,6:127] = PSFmodel


    print(' - Converting the '+str(Nmodels)+' models found in \n   '+modeldir)
    tu.galfit_convertmodel2cube(models,savecubesumimg=True,includewcs=True,sourcecat_compinfo=componentinfo,
                                convkernels=[convkernel]*Nmodels,normalizecomponents=False)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def OIIemitters_GenerateComponentInfo():
    """

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    tsu.OIIemitters_GenerateComponentInfo()

    """
    modeldir = '/Volumes/DATABCKUP1/TDOSEextractions/181016_MWDR1_OIIemitters/galfit_models/'
    models   = glob.glob(modeldir+'model_acs_814w_candels-cdfs-*_cut_v1.0_id*_cutout4p0x4p0arcsec.fits')

    compinfofile = '/Users/kschmidt/work/TDOSE/181016_MWDR1_OIIemitters/OIIemitter_componentinfo_RENAME.txt'
    fout = open(compinfofile,'w')
    fout.write("""# TDOSE source catalog components keys for 4x4 arcsec GALFIT models of MUSE-Wide DR OII emitters from
# /Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/galfit_wrapper_results_final_181015/OIIemitter_selection_23LTm814LT24_SkelsepLT0p3_Nobj153.fits
#
# --- TEMPLATE --- generated with tdosepublication_utilities.OIIemitters_GenerateComponentInfo() on %s
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# modelfilename        Path and name of model file
# id                   MUSE-Wide object ID
# componentinfo        Information on the model components given as ComponentNumber:InfoKey
#                      where the info keys are:  1 = object, 2 = contaminant and 3 = sky
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# modefilename  id  componentinfo
""" % tu.get_now_string())

    for mm, GFmodel in enumerate(models):
        objid       = GFmodel.split('/')[-1].split('_id')[-1].split('_cutout')[0]
        modelheader = afits.open(GFmodel)[2].header
        compstring  = ' '
        for key in modelheader.keys():
            if 'COMP_' in key:
                compNo = key.split('OMP_')[-1]
                if modelheader[key] == 'sky':
                    compstring = compstring + compNo + ':3  '
                else:
                    compstring = compstring + compNo + ':?  '

        outstring = GFmodel.split('/')[-1]+'  '+objid+'  '+compstring.ljust(50)
        fout.write(outstring+' \n')
    fout.close()
    print(' - Wrote component info to: '+compinfofile)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def OIIemitters_GenerateMWcoordinateSorcecats():
    """
    Generate source catalogs for the OII emitters based on the MW coordinates for generating cutouts

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    tsu.OIIemitters_GenerateMWcoordinateSorcecats()

    """
    outdir       = '/Volumes/DATABCKUP1/TDOSEextractions/181016_MWDR1_OIIemitters/tdose_sourcecats_MWcoord/'
    modeldir     = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/galfit_wrapper_results_final_181015/'
    fitscatalog  = modeldir+'OIIemitter_selection_23LTm814LT24_SkelsepLT0p3_Nobj153.fits'
    fitsdat      = afits.open(fitscatalog)[1].data

    for oo, objid in enumerate(fitsdat['UNIQUE_ID']):

        imgheader  = afits.open('/Users/kschmidt/work/images_MAST/MUSEWidePointings/acs_814w_'+
                                fitsdat['FIELD_NAME'][oo]+'_cut_v1.0.fits')[0].header
        sourcecatcenter = [fitsdat['RA'][oo],  fitsdat['DEC'][oo]]
        sourcecatradius = 45.0
        outname         = outdir+'tdose_sourcecat_MWcoordinates-'+str(int(objid))+'.txt'

        sourcecat    = tu.gen_sourcecat_from_FitsCat(fitscatalog,'UNIQUE_ID','RA','DEC',sourcecatcenter,
                                                     sourcecatradius ,imgheader,outname=outname,clobber=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def OIIemitters_plotcomparisons(showids=False, outliercut=1.25, figsize=(5, 5), fontsize=12, logaxes=True,
                                namebase='/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/OIIemitter_comparisons'):
    """
    Plotting comparison of extractions of OII emitters to evaluate performance.

    with the introduction of tsu.get_maxvalues(), tsu.extractMaxValues() and tsu.plot_maxvalues_OIIemitters()
    this routine is now obsolete...

    --- INPUT ---
    N/A

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    tsu.OIIemitters_plotcomparisons()

    """
    specdir        = '/Volumes/DATABCKUP1/TDOSEextractions/181016_MWDR1_OIIemitters/tdose_spectra/'
    files_modelimg = glob.glob(specdir+'tdose_spectrum_modelimg*.fits')
    files_migauss  = glob.glob(specdir+'tdose_spectrum_gauss*.fits')

    galfitdir      = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/galfit_wrapper_results_final_181015/'
    objcatalog     = afits.open(galfitdir+'OIIemitter_selection_23LTm814LT24_SkelsepLT0p3_Nobj153.fits')[1].data

    MWDR1dir       = '/Users/kschmidt/work/TDOSE/OIIemitter_MWDR1spectra/'

    print(' - Looping over objects and pulling out maximum OII flux and S/N from spectra')
    ids         = []
    SNsersic    = []
    SNgauss     = []
    Fsersic     = []
    Ferrsersic  = []
    Fgauss      = []
    Ferrgauss   = []

    SNlsdcat    = []

    SNmwdr1el   = []
    SNmwdr1td   = []
    Fmwdr1el    = []
    Fmwdr1td    = []
    Ferrmwdr1el = []
    Ferrmwdr1td = []

    Fsersic_filtered     = []
    Ferrsersic_filtered  = []
    SNsersic_filtered    = []

    for ff, file_mi in enumerate(files_modelimg):
        if '116004061' in file_mi:  # Ignore duplicated object
            continue
        objid       = file_mi.split('modelimg_0')[-1][:9]
        objid_check = files_migauss[ff].split('gauss_0')[-1][:9]
        if objid != objid_check:
            sys.exit('mismatch in ids... objid='+str(objid)+' and objid_check='+str(objid_check))
        else:
            ids.append(objid)

        infostr = ' - Getting info for object '+str(objid)+' ( file '+\
                  str("%3.f" % (ff+1))+' / '+str("%3.f" % len(files_modelimg)+')')
        sys.stdout.write("%s\r" % infostr)
        sys.stdout.flush()

        catent          = np.where(objcatalog['UNIQUE_ID'].astype(str) == objid)[0]
        objz            = objcatalog['Z'][catent]
        waveOII         = (1+objz)*3726.0
        SNlsdcat.append(objcatalog['SN'][catent][0])

        # - - - - - - - - - MODELIMG EXTRACTIONS - - - - - - - - -
        dat_modelimg    = afits.open(file_mi)[1].data
        dat_gauss       = afits.open(files_migauss[ff])[1].data

        wavediff        = np.abs(dat_modelimg['wave']-waveOII)
        waveent         = np.where(wavediff == np.min(wavediff))[0]

        dwave           = 10
        subdat_S2N_ser  = dat_modelimg['s2n'][waveent[0]-dwave:waveent[0]+dwave]
        subdat_S2N_gau  = dat_gauss['s2n'][waveent[0]-dwave:waveent[0]+dwave]
        subdat_F_ser    = dat_modelimg['flux'][waveent[0]-dwave:waveent[0]+dwave]
        subdat_F_gau    = dat_gauss['flux'][waveent[0]-dwave:waveent[0]+dwave]
        fmaxent_ser     = np.where(subdat_F_ser == np.max(subdat_F_ser))[0]
        fmaxent_gau     = np.where(subdat_F_gau == np.max(subdat_F_gau))[0]
        subdat_Ferr_ser = dat_modelimg['fluxerror'][waveent[0]-dwave:waveent[0]+dwave][fmaxent_ser]
        subdat_Ferr_gau = dat_gauss['fluxerror'][waveent[0]-dwave:waveent[0]+dwave][fmaxent_gau]

        SNgauss.append(np.max(subdat_S2N_ser))
        SNgauss.append(np.max(subdat_S2N_gau))
        Fsersic.append(np.max(subdat_F_ser))
        Fgauss.append(np.max(subdat_F_gau))
        Ferrsersic.append(subdat_Ferr_ser[0])
        Ferrgauss.append(subdat_Ferr_gau[0])

        # v v v v LSDCatify spectrum and pull out S/N, flux and variance v v v v
        wave, flux_filtered, variance_filtered = tsu.LSDCatify_spectrum(file_mi,verbose=False)
        s2n_filtered         = flux_filtered / np.sqrt(variance_filtered)
        subdat_S2N_serfilt   = s2n_filtered[waveent[0]-dwave:waveent[0]+dwave]
        subdat_F_serfilt     = flux_filtered[waveent[0]-dwave:waveent[0]+dwave]
        fmaxent_serfilt      = np.where(subdat_F_serfilt == np.max(subdat_F_serfilt))[0]
        subdat_Ferr_serfilt  = np.sqrt(variance_filtered)[waveent[0]-dwave:waveent[0]+dwave][fmaxent_serfilt]

        SNsersic_filtered.append(np.max(subdat_S2N_serfilt))
        Fsersic_filtered.append(np.max(subdat_F_serfilt))
        Ferrsersic_filtered.append(subdat_Ferr_serfilt[0])
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^

        # - - - - - - - - - MWDR1 EXTRACTIONS - EMISSION LINES - - - - - - - - -
        file_MWDR1el  = MWDR1dir+'spectrum_'+objid+'.fits'
        spec_MWDR1el  = afits.open(file_MWDR1el)[1].data

        wavediff        = np.abs(spec_MWDR1el['WAVE_AIR']-waveOII)
        waveent         = np.where(wavediff == np.min(wavediff))[0]
        dwave           = 10

        SNvec           = spec_MWDR1el['FLUX']/spec_MWDR1el['FLUXERR']
        subdat_S2N      = SNvec[waveent[0]-dwave:waveent[0]+dwave]
        subdat_F        = spec_MWDR1el['FLUX'][waveent[0]-dwave:waveent[0]+dwave]
        fmaxent         = np.where(subdat_F == np.max(subdat_F))[0]
        Ferrmax         = spec_MWDR1el['FLUXERR'][waveent[0]-dwave:waveent[0]+dwave][fmaxent]

        SNmwdr1el.append(np.max(subdat_S2N))
        Fmwdr1el.append(np.max(subdat_F))
        Ferrmwdr1el.append(Ferrmax[0])

        # - - - - - - - - - MWDR1 EXTRACTIONS - TDOSE GUO - - - - - - - - -
        guoidstr        = str("%.10d" % objcatalog['GUO_ID'][catent])
        file_MWDR1td    = MWDR1dir+'tdose_spectrum_'+objcatalog['FIELD_NAME'][catent][0]+'_gauss_'+guoidstr+'.fits'
        if guoidstr == '0000005408':
            file_MWDR1td = file_MWDR1td.replace(objcatalog['FIELD_NAME'][catent][0],'candels-cdfs-05')
        if guoidstr == '0000006034':
            file_MWDR1td = file_MWDR1td.replace(objcatalog['FIELD_NAME'][catent][0],'candels-cdfs-09')
        if guoidstr == '0000006561':
            file_MWDR1td = file_MWDR1td.replace(objcatalog['FIELD_NAME'][catent][0],'candels-cdfs-14')
        if guoidstr == '0000005499':
            file_MWDR1td = file_MWDR1td.replace(objcatalog['FIELD_NAME'][catent][0],'candels-cdfs-16')
        if guoidstr == '0000013882':
            file_MWDR1td = file_MWDR1td.replace(objcatalog['FIELD_NAME'][catent][0],'candels-cdfs-46')
        if guoidstr == '0000013882':
            file_MWDR1td = file_MWDR1td.replace(objcatalog['FIELD_NAME'][catent][0],'candels-cdfs-45')

        spec_MWDR1td    = afits.open(file_MWDR1td)[1].data

        wavediff        = np.abs(spec_MWDR1td['wave']-waveOII)
        waveent         = np.where(wavediff == np.min(wavediff))[0]
        dwave           = 10

        subdat_S2N      = spec_MWDR1td['s2n'][waveent[0]-dwave:waveent[0]+dwave]
        subdat_F        = spec_MWDR1td['flux'][waveent[0]-dwave:waveent[0]+dwave]
        fmaxent         = np.where(subdat_F == np.max(subdat_F))[0]
        Ferrmax         = spec_MWDR1td['fluxerror'][waveent[0]-dwave:waveent[0]+dwave][fmaxent]

        SNmwdr1td.append(np.max(subdat_S2N))
        Fmwdr1td.append(np.max(subdat_F))
        Ferrmwdr1td.append(Ferrmax[0])

    print('\n   ... done')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    odrfit_F = kbs.fit_function_to_data_with_errors_on_both_axes(Fsersic,Fmwdr1td,Ferrsersic,Ferrmwdr1td,
                                                               fitfunction='linear',plotresults=namebase+'_ODRfit2data.pdf')
    print(' - A bit of stats on the input data for the ODR fit:')
    print('   -> np.median(np.asarray(Fsersic)/np.asarray(Fmwdr1td))  = '+
          str(np.median(np.asarray(Fsersic)/np.asarray(Fmwdr1td))))
    print('   -> np.std(np.asarray(Fsersic)/np.asarray(Fmwdr1td))  = '+
          str(np.std(np.asarray(Fsersic)/np.asarray(Fmwdr1td))))

    # pdb.set_trace() #... 181107 why does the ODR not work? Funky shape of input arrays?
    odrfit_SN = kbs.fit_function_to_data_with_errors_on_both_axes(SNsersic_filtered,SNlsdcat,
                                                                  np.asarray(SNlsdcat)*0.0+1.0,np.asarray(SNlsdcat)*0.0+1.0,
                                                                  fitfunction='linear',
                                                                  plotresults=namebase+'_ODRfit2data_SN.pdf')
    print(' - A bit of stats on the input data for the ODR fit:')
    print('   -> np.median(np.asarray(SNsersic_filtered)/np.asarray(SNlsdcat))  = '+
          str(np.median(np.asarray(SNsersic_filtered)/np.asarray(SNlsdcat))))
    print('   -> np.std(np.asarray(SNsersic_filtered)/np.asarray(SNlsdcat))  = '+
          str(np.std(np.asarray(SNsersic_filtered)/np.asarray(SNlsdcat))))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2N.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.plot(SNsersic,SNgauss,'ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if (SNsersic[ii]/SNgauss[ii] > outliercut):
                plt.text(SNsersic[ii],SNgauss[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            if (SNgauss[ii]/SNsersic[ii] > outliercut):
                plt.text(SNsersic[ii],SNgauss[ii],id,fontsize=Fsize-5,horizontalalignment='right')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [2,100]
    else:
        axisrange = [0,70]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot([5,5],axisrange,'--k',lw=lthick)
    plt.plot(axisrange,[5,5],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('max(S/N) TDOSE multicomponent Gauss model', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    # leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1,
    #                  bbox_to_anchor=(0.01, 0.99))  # add the legend
    # leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2Nlsdcat.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    scalefactor = 1.5
    plt.plot(SNsersic,np.asarray(SNlsdcat)/scalefactor,'.r',label='S/N LSDcat scaled by '+str(scalefactor))

    plt.plot(SNsersic,SNlsdcat,'ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if (SNsersic[ii]/(SNlsdcat[ii]/scalefactor) > outliercut):
                plt.text(SNsersic[ii],SNlsdcat[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            if ((SNlsdcat[ii]/scalefactor)/SNsersic[ii] > outliercut):
                plt.text(SNsersic[ii],SNlsdcat[ii],id,fontsize=Fsize-5,horizontalalignment='right')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [2,100]
    else:
        axisrange = [0,70]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot([5,5],axisrange,'--k',lw=lthick)
    plt.plot(axisrange,[5,5],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('S/N MW DR1 LSDcat measurements', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1,
                     bbox_to_anchor=(0.01, 0.99))  # add the legend
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2Nlsdcat_filtered.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # scalefactor = 1.5
    # plt.plot(SNsersic,np.asarray(SNlsdcat)/scalefactor,'.r',label='S/N LSDcat scaled by '+str(scalefactor))
    # plt.plot(SNsersic_filtered,SNlsdcat,'.g',label='LSDcatified TDOSE S/N')

    plt.plot(SNsersic_filtered,SNlsdcat,'ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if (SNsersic[ii]/(SNlsdcat[ii]/scalefactor) > outliercut):
                plt.text(SNsersic[ii],SNlsdcat[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            if ((SNlsdcat[ii]/scalefactor)/SNsersic[ii] > outliercut):
                plt.text(SNsersic[ii],SNlsdcat[ii],id,fontsize=Fsize-5,horizontalalignment='right')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # overplot ODR fit
    popt = odrfit_SN.beta
    perr = odrfit_SN.sd_beta
    nstd = 3. # draw 5-sigma intervals
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    x_fit = np.linspace(min(SNsersic_filtered), max(SNsersic_filtered), 100)
    fit = kbs.odr_fitfunction_linear(popt, x_fit)
    fit_up = kbs.odr_fitfunction_linear(popt_up, x_fit)
    fit_dw= kbs.odr_fitfunction_linear(popt_dw, x_fit)

    plt.fill_between(x_fit, fit_up, fit_dw, alpha=0.25, label='3-sigma interval',color='red')
    plt.plot(x_fit, fit, 'r', lw=2, label='Best fit curve')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [2,100]
    else:
        axisrange = [0,70]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot([5,5],axisrange,'--k',lw=lthick)
    plt.plot(axisrange,[5,5],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT multicomponent Sersic model (filtered)', fontsize=Fsize)
    plt.ylabel('S/N MW DR1 LSDcat measurements', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    # leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1,
    #                  bbox_to_anchor=(0.01, 0.99))  # add the legend
    # leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2Nmwdr1el.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    scalefactor = 1.0
    if scalefactor != 1.0:
        plt.plot(SNsersic,np.asarray(SNmwdr1el)/scalefactor,'.r',label='S/N EL spectrum scaled by '+str(scalefactor))

    plt.plot(SNsersic,SNmwdr1el,'ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if (SNsersic[ii]/SNmwdr1el[ii] > outliercut):
                plt.text(SNsersic[ii],SNmwdr1el[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            if (SNmwdr1el[ii]/SNsersic[ii] > outliercut):
                plt.text(SNsersic[ii],SNmwdr1el[ii],id,fontsize=Fsize-5,horizontalalignment='right')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [2,100]
    else:
        axisrange = [0,50]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot([5,5],axisrange,'--k',lw=lthick)
    plt.plot(axisrange,[5,5],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('max(S/N) MW DR1 EL spectrum', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    if scalefactor != 1.0:
        leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1,
                         bbox_to_anchor=(0.01, 0.99))  # add the legend
        leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2Nmwdr1td.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    scalefactor = 1.0
    if scalefactor != 1.0:
        plt.plot(SNsersic,np.asarray(SNmwdr1td)/scalefactor,'.r',label='S/N TDOSE spectrum scaled by '+str(scalefactor))

    plt.plot(SNsersic,SNmwdr1td,'ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if (SNsersic[ii]/SNmwdr1td[ii] > outliercut):
                plt.text(SNsersic[ii],SNmwdr1td[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            if (SNmwdr1td[ii]/SNsersic[ii] > outliercut):
                plt.text(SNsersic[ii],SNmwdr1td[ii],id,fontsize=Fsize-5,horizontalalignment='right')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [2,100]
    else:
        axisrange = [0,50]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot([5,5],axisrange,'--k',lw=lthick)
    plt.plot(axisrange,[5,5],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('max(S/N) MW DR1 TDOSE spectrum', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    if scalefactor != 1.0:
        leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1,
                         bbox_to_anchor=(0.01, 0.99))  # add the legend
        leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_flux.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.errorbar(Fsersic,Fgauss, yerr=Ferrgauss, xerr=Ferrsersic,fmt='ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if (Fsersic[ii]/Fgauss[ii] > outliercut):
                plt.text(Fsersic[ii],Fgauss[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            if (Fgauss[ii]/Fsersic[ii] > outliercut):
                plt.text(Fsersic[ii],Fgauss[ii],id,fontsize=Fsize-5,horizontalalignment='right')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [50,6000]
    else:
        axisrange = [0,5000]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(F/[cgs]) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('max(F/[cgs]) TDOSE multicomponent Gauss model', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    # leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1,
    #                  bbox_to_anchor=(0.01, 0.99))  # add the legend
    # leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_flux_mwdr1el.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    scalefactor = 2.0
    if scalefactor != 1.0:
        plt.plot(Fsersic,np.asarray(Fmwdr1el)/scalefactor,'.r',label='S/N MW DR1 EL spectrum scaled by '+str(scalefactor))

    plt.errorbar(Fsersic,Fmwdr1el, yerr=Ferrgauss, xerr=Ferrsersic,fmt='ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if (Fsersic[ii]/(Fmwdr1el[ii]/scalefactor) > outliercut):
                plt.text(Fsersic[ii],Fmwdr1el[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            if ((Fmwdr1el[ii]/scalefactor)/Fsersic[ii] > outliercut):
                plt.text(Fsersic[ii],Fmwdr1el[ii],id,fontsize=Fsize-5,horizontalalignment='right')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [100,4000]
    else:
        axisrange = [0,3000]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(F/[cgs]) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('max(F/[cgs]) MW DR1 EL spectrum', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    if scalefactor != 1.0:
        leg = plt.legend(fancybox=True, loc='lower right',prop={'size':Fsize},ncol=1,numpoints=1)
        leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_flux_mwdr1td.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.errorbar(Fsersic,Fmwdr1td, yerr=Ferrgauss, xerr=Ferrsersic,fmt='ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if (Fsersic[ii]/Fmwdr1td[ii] > outliercut):
                plt.text(Fsersic[ii],Fmwdr1td[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            if (Fmwdr1td[ii]/Fsersic[ii] > outliercut):
                plt.text(Fsersic[ii],Fmwdr1td[ii],id,fontsize=Fsize-5,horizontalalignment='right')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # overplot ODR fit
    popt = odrfit_F.beta
    perr = odrfit_F.sd_beta
    nstd = 3. # draw 5-sigma intervals
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    x_fit = np.linspace(min(Fsersic), max(Fsersic), 100)
    fit = kbs.odr_fitfunction_linear(popt, x_fit)
    fit_up = kbs.odr_fitfunction_linear(popt_up, x_fit)
    fit_dw= kbs.odr_fitfunction_linear(popt_dw, x_fit)

    plt.fill_between(x_fit, fit_up, fit_dw, alpha=0.25, label='3-sigma interval',color='red')
    plt.plot(x_fit, fit, 'r', lw=2, label='Best fit curve')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [100,4000]
    else:
        axisrange = [0,1500]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(F/[cgs]) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('max(F/[cgs]) MW DR1 TDOSE spectrum', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def OIIemitters_WhiteLightImages(outputdir,overwrite=True,skipImageCreation=False,printfilenames=False,verbose=True):
    """
    Function to generate white light images for the OII emitters

    --- INPUT ---
    outputdir       Directory to contain white light images

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    outputdir = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/DATE_WhiteLightImages/'
    tsu.OIIemitters_WhiteLightImages(outputdir)

    """
    modeldir        = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/galfit_wrapper_results_final_181015/'
    fitscatalog     = modeldir+'OIIemitter_selection_23LTm814LT24_SkelsepLT0p3_Nobj153.fits'
    datacubedir     = '/Volumes/DATABCKUP1/MUSE-Wide/DATACUBES/'
    tdosemodeldir   = '/Volumes/DATABCKUP1/TDOSEextractions/181016_MWDR1_OIIemitters/tdose_models/'
    tdosecutoutdir  = '/Volumes/DATABCKUP1/TDOSEextractions/181016_MWDR1_OIIemitters/tdose_cutouts/'
    tdosespecdir    = '/Volumes/DATABCKUP1/TDOSEextractions/181016_MWDR1_OIIemitters/tdose_spectra/'
    galfitmodeldir  = '/Volumes/DATABCKUP1/TDOSEextractions/181016_MWDR1_OIIemitters/galfit_models/'

    objdat   = afits.open(fitscatalog)[1].data

    for oo, id in enumerate(objdat['UNIQUE_ID']):
        if verbose: print(' >>> Generating images for object '+str(oo+1)+' / '+str(len(objdat['UNIQUE_ID']))+' <<<')
        pointing        = objdat['FIELD_NAME'][oo]
        ra              = [objdat['RA'][oo]]
        Dra             = 4.0 # arcsec
        dec             = [objdat['DEC'][oo]]
        Ddec            = Dra
        objid           = str(objdat['UNIQUE_ID'][oo])

        datacubestr     = datacubedir+'DATACUBE_'+pointing+'_v1.0_dcbgc_effnoised.fits'
        datacube        = glob.glob(datacubestr)
        if len(datacube) != 1:
            print('\n\n WARNING - the datacube '+datacubestr+' was not found... continuing\n\n')
            continue

        sourcemodelstr  = tdosespecdir+'tdose_spectrum_modelimg*'+objid+'*.fits'
        sourcemodel     = glob.glob(sourcemodelstr)
        if len(sourcemodel) != 1:
            print('\n\n WARNING - the spectrum with source model '+sourcemodelstr+' was not found... continuing\n\n')
            continue

        if verbose: print(' ---- Extracting white light image and cube for '+objid+' ---- ')
        wcenter   = [[7050]]
        dwave     = [[2250]]
        names     = [objid+'DATACUBE_whitelight',objid+'SOURCEMODEL_whitelight']
        for cc, dcube in enumerate([datacube[0],sourcemodel[0]]):
            if not skipImageCreation:
                mu.create_narrowband_subcube(dcube,ra,dec,Dra,Ddec,wcenter,dwave,outputdir,
                                             cube_ext=['DATA_DCBGC'],names=[names[cc]],clobber=overwrite)
        imgwhiteDC  = glob.glob(outputdir+'*'+names[0]+'*narrowbandimage*fits')
        imgwhiteSM  = glob.glob(outputdir+'*'+names[1]+'*narrowbandimage*fits')

        if verbose: print(' ---- Extracting OII image and cube for '+objid+' ---- ')
        redshift  = float(str(objdat['Z'][oo]))
        linewave  = 3727.5
        wcenter   = [[linewave*(redshift+1.0)]]
        dwave     = [[500.0/299792.0 * linewave * (redshift+1.0)]]# narrowband width is 2x500=1000 km/s rest-frame
        names     = [objid+'DATACUBE_OIIdoubletWidth1000kmsRest',objid+'SOURCEMODEL_OIIdoubletWidth1000kmsRest']
        for cc, dcube in enumerate([datacube[0],sourcemodel[0]]):
            if not skipImageCreation:
                mu.create_narrowband_subcube(dcube,ra,dec,Dra,Ddec,wcenter,dwave,outputdir,
                                             cube_ext=['DATA_DCBGC'],names=[names[cc]],clobber=overwrite)
        imgOIIDC    = glob.glob(outputdir+'*'+names[0]+'*narrowbandimage*fits')
        imgOIISM    = glob.glob(outputdir+'*'+names[1]+'*narrowbandimage*fits')

        if verbose: print(' ---- Globbing for other images of object ---- ')
        modelstr       = tdosemodeldir+'acs_814w_'+pointing+'*'+objid+'_cutout*tdose_modelimage_gauss.fits'
        acsmodel       = glob.glob(modelstr)

        if len(acsmodel) != 1:
            print('\n\n WARNING - the ACS model image '+modelstr+' was not found... continuing\n\n')
            acsmodel          = [' ']
            modelGaussCubeWCS = [' ']
        else:
            cubeWCSstr        = acsmodel[0].replace('_gauss.fits','_cubeWCS_gauss.fits')
            modelGaussCubeWCS = glob.glob(cubeWCSstr)

        galfitmodelstr = galfitmodeldir+'model_acs_814w_'+pointing+'*'+objid+'_cutout4p0x4p0arcsec_cubesum.fits'
        galfitmodel    = glob.glob(galfitmodelstr)

        if len(galfitmodel) != 1:
            print('\n\n WARNING - the ACS GALFIT model cubesum image '+galfitmodelstr+' was not found... continuing\n\n')
            galfitmodel      = [' ']
            galfitModCubeWCS = [' ']
        else:
            cubeWCSstr   = acsmodel[0].replace('_gauss.fits','_cubeWCS_modelimg.fits')
            galfitModCubeWCS = glob.glob(cubeWCSstr)

        refimage  = glob.glob(tdosecutoutdir+'acs_814w_'+pointing+'_cut_*'+objid+'*.fits')

        ds9cmd = 'ds9 -scale mode minmax '+\
                 imgwhiteDC[0]+' '+imgwhiteSM[0]+' '+imgOIIDC[0]+' '+imgOIISM[0]+' '+modelGaussCubeWCS[0]+' '+acsmodel[0]+' '+refimage[0]+' '+galfitModCubeWCS[0]+' '+galfitmodel[0]+\
                 ' -lock frame wcs -tile grid layout 9 1 -geometry 1400x500 -zoom to fit &'
        if verbose: print(ds9cmd)

        if printfilenames:
            print(imgwhiteDC[0])
            print(imgwhiteSM[0])
            print(imgOIIDC[0])
            print(imgOIISM[0])
            print(modelGaussCubeWCS[0])
            print(acsmodel[0])
            print(refimage[0])
            print(galfitModCubeWCS[0])
            print(galfitmodel[0])

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def OIIemitters_WhiteLightImages_estimatesize(outputdir,datestr='XXXXXX',overwrite=True,verbose=True):
    """
    Estimating size of white light image (and others incl. ref image) content for the OII emitters

    --- INPUT ---
    outputdir       Directory to contain white light images

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    outputdir  = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/190212_WhiteLightImages/sizeestimates_1902XX/'
    tsu.OIIemitters_WhiteLightImages_estimatesize(outputdir,datestr='1902XX')

    """
    imagefiles = outputdir+'../filenames.txt'
    filenames  = np.genfromtxt(imagefiles,skip_header=4,dtype="S190",comments='#')

    resultstxt = outputdir+'fitparameters'+datestr+'.txt'
    resultsout = open(resultstxt,'w')
    resultsout.write('# The fit parameters obtained with tdosepublication_utilities.OIIemitters_WhiteLightImages_estimatesize()\n')
    resultsout.write('# on '+datestr+'\n')
    resultsout.write('# \n')
    resultsout.write('# objid xpos ypos fluxscale xsigma ysigma angle '
                  'xpos_init ypos_init fluxscale_init xsigma_init ysigma_init angle_init imagefit \n')

    for fn in filenames:
        dataimg    = afits.open(fn)[0].data
        try:
            objid = int(fn.split('/')[-1].split('_')[1])
        except:
            objid = int(fn.split('/')[-1].split('_id')[-1].split('_')[0])

        txtcat    = outputdir+fn.split('/')[-1].replace('.fits','_sourcecat_automatic.txt')
        fout      = open(txtcat,'w')
        fout.write('# object xpos  ypos  fluxscale  \n')
        imgcenter = int(dataimg.shape[1]/2.)
        icstring  = str(imgcenter)
        fout.write(str(objid)+'  '+icstring+'  '+icstring+'  '+str(dataimg[imgcenter,imgcenter])+' \n')
        fout.close()
        fitscat = tu.ascii2fits(txtcat,asciinames=True,skip_header=0,outpath=outputdir,verbose=verbose)

        param_init, fit_output = \
            tmf.gen_fullmodel(dataimg,fitscat,xpos_col='xpos',ypos_col='ypos',
                              sigysigxangle=None, #np.array([[imgcenter/5.,imgcenter/10.,45.0]]),
                              fluxscale='fluxscale',
                              generateimage=outputdir+fn.split('/')[-1].replace('.fits','_model.fits'),
                              generateresidualimage=outputdir+fn.split('/')[-1].replace('.fits','_modelresidual.fits'),
                              max_centroid_shift=int(imgcenter*0.5),
                              verbose=verbose,datanoise=None,clobber=overwrite)

        outputstr = str(fit_output[0][1::6][0])+' '+str(fit_output[0][0::6][0])+' '+str(fit_output[0][2::6][0])+' '+\
                    str(fit_output[0][4::6][0])+' '+str(fit_output[0][3::6][0])+' '+str(fit_output[0][5::6][0])+'       '+\
                    str(param_init[1::6][0])+' '+str(param_init[0::6][0])+' '+str(param_init[2::6][0])+' '+\
                    str(param_init[4::6][0])+' '+str(param_init[3::6][0])+' '+str(param_init[5::6][0])
        resultsout.write(str(objid)+'      '+outputstr+'       '+fn+'\n')
    resultsout.close()
    fitsfmt     = ['K']+['D']*12+['A300']
    resultsfits = tu.ascii2fits(resultstxt,asciinames=True,skip_header=3,outpath=outputdir,fitsformat=fitsfmt,verbose=verbose)

    print('\n\n - Check the residuals with')
    print('ds9 -scale minmax '+outputdir+'*residual.fits  -lock frame image & ')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def OIIemitters_compare_modeparams(paramsummary,plotdir,verbose=True):
    """
    Diagnostic plots and calculations comparing OII emitters to continuum images

    --- INPUT ---
    paramsummary      Fits file generated by tsu.OIIemitters_WhiteLightImages_estimatesize() summarizing fits

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu

    maindir   = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/190212_WhiteLightImages/sizeestimates/'
    plotdir   = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/190212_WhiteLightImages/sizeestimates_plots/'
    paramfile = maindir+'fitparameters190213_153obj.fits'

    maindir   = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/190212_WhiteLightImages/sizeestimates_190213/'
    plotdir   = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/190212_WhiteLightImages/sizeestimates_190213_plots/'
    paramfile = maindir+'fitparameters190213.fits'

    tsu.OIIemitters_compare_modeparams(paramfile,plotdir)

    """
    paramdat   = afits.open(paramsummary)[1].data
    Nparamset  = len(paramdat['objid'])
    unique_ids = np.unique(paramdat['objid'])
    Nobj       = len(unique_ids)
    if verbose:
        print(' - Found '+str(Nobj)+' unique IDs in parameter summary file to evaluate and plot.')
        print(' - With '+str(Nparamset)+' entries in file, that makes for '+str(float(Nparamset)/float(Nobj))+' files/object')

    dat_whitelight = paramdat[0::5]
    dat_OIIband    = paramdat[1::5]
    dat_CubeWCSMod = paramdat[2::5]
    dat_ACSmodel   = paramdat[3::5]
    dat_ACSref     = paramdat[4::5]

    ACSpix2arcsec  = 0.03 # arcsec/pix
    MUSEpix2arcsec = 0.20 # arcsec/pix
    # ---------------------------------------------------------------------------------------------------------
    # ---------------- MINOR AXES ----------------
    WLminoraxes         = np.zeros(Nobj)
    OIIminoraxes        = np.zeros(Nobj)
    CubeWCSModminoraxes = np.zeros(Nobj)
    ACSMODminoraxes     = np.zeros(Nobj)
    REFminoraxes        = np.zeros(Nobj)

    for ii, id in enumerate(unique_ids):
        WLminoraxes[ii]         = np.min([dat_whitelight['xsigma'][ii],dat_whitelight['ysigma'][ii]]) * MUSEpix2arcsec
        OIIminoraxes[ii]        = np.min([dat_OIIband['xsigma'][ii],dat_OIIband['ysigma'][ii]])       * MUSEpix2arcsec
        CubeWCSModminoraxes[ii] = np.min([dat_CubeWCSMod['xsigma'][ii],dat_CubeWCSMod['ysigma'][ii]])   * MUSEpix2arcsec
        ACSMODminoraxes[ii]     = np.min([dat_ACSmodel['xsigma'][ii],dat_ACSmodel['ysigma'][ii]])     * ACSpix2arcsec
        REFminoraxes[ii]        = np.min([dat_ACSref['xsigma'][ii],dat_ACSref['ysigma'][ii]])         * ACSpix2arcsec

    diff_WLvsOII_minor           = WLminoraxes-OIIminoraxes
    med_WLvsOII_minor            = np.median(diff_WLvsOII_minor)
    std_WLvsOII_minor            = np.std(diff_WLvsOII_minor)
    diff_OIIvsCubeWCSMod_minor   = OIIminoraxes-CubeWCSModminoraxes
    med_OIIvsCubeWCSMod_minor    = np.median(diff_OIIvsCubeWCSMod_minor)
    std_OIIvsCubeWCSMod_minor    = np.std(diff_OIIvsCubeWCSMod_minor)
    diff_OIIvsACSMOD_minor       = OIIminoraxes-ACSMODminoraxes
    med_OIIvsACSMOD_minor        = np.median(diff_OIIvsACSMOD_minor)
    std_OIIvsACSMOD_minor        = np.std(diff_OIIvsACSMOD_minor)
    diff_OIIvsREF_minor          = OIIminoraxes-REFminoraxes
    med_OIIvsREF_minor           = np.median(diff_OIIvsREF_minor)
    std_OIIvsREF_minor           = np.std(diff_OIIvsREF_minor)
    diff_WLvsCubeWCSMod_minor    = WLminoraxes-CubeWCSModminoraxes
    med_WLvsCubeWCSMod_minor     = np.median(diff_WLvsCubeWCSMod_minor)
    std_WLvsCubeWCSMod_minor     = np.std(diff_WLvsCubeWCSMod_minor)
    diff_WLvsACSMOD_minor        = WLminoraxes-ACSMODminoraxes
    med_WLvsACSMOD_minor         = np.median(diff_WLvsACSMOD_minor)
    std_WLvsACSMOD_minor         = np.std(diff_WLvsACSMOD_minor)
    diff_WLvsREF_minor           = WLminoraxes-REFminoraxes
    med_WLvsREF_minor            = np.median(diff_WLvsREF_minor)
    std_WLvsREF_minor            = np.std(diff_WLvsREF_minor)

    # ---------------------------------------------------------------------------------------------------------
    # ---------------- MAJOR AXES ----------------
    WLmajoraxes         = np.zeros(Nobj)
    OIImajoraxes        = np.zeros(Nobj)
    CubeWCSModmajoraxes = np.zeros(Nobj)
    ACSMODmajoraxes     = np.zeros(Nobj)
    REFmajoraxes        = np.zeros(Nobj)

    for ii, id in enumerate(unique_ids):
        WLmajoraxes[ii]         = np.max([dat_whitelight['xsigma'][ii],dat_whitelight['ysigma'][ii]]) * MUSEpix2arcsec
        OIImajoraxes[ii]        = np.max([dat_OIIband['xsigma'][ii],dat_OIIband['ysigma'][ii]])       * MUSEpix2arcsec
        CubeWCSModmajoraxes[ii] = np.max([dat_CubeWCSMod['xsigma'][ii],dat_CubeWCSMod['ysigma'][ii]])   * MUSEpix2arcsec
        ACSMODmajoraxes[ii]     = np.max([dat_ACSmodel['xsigma'][ii],dat_ACSmodel['ysigma'][ii]])     * ACSpix2arcsec
        REFmajoraxes[ii]        = np.max([dat_ACSref['xsigma'][ii],dat_ACSref['ysigma'][ii]])         * ACSpix2arcsec

    diff_WLvsOII_major           = WLmajoraxes-OIImajoraxes
    med_WLvsOII_major            = np.median(diff_WLvsOII_major)
    std_WLvsOII_major            = np.std(diff_WLvsOII_major)
    diff_OIIvsCubeWCSMod_major   = OIImajoraxes-CubeWCSModmajoraxes
    med_OIIvsCubeWCSMod_major    = np.median(diff_OIIvsCubeWCSMod_major)
    std_OIIvsCubeWCSMod_major    = np.std(diff_OIIvsCubeWCSMod_major)
    diff_OIIvsACSMOD_major       = OIImajoraxes-ACSMODmajoraxes
    med_OIIvsACSMOD_major        = np.median(diff_OIIvsACSMOD_major)
    std_OIIvsACSMOD_major        = np.std(diff_OIIvsACSMOD_major)
    diff_OIIvsREF_major          = OIImajoraxes-REFmajoraxes
    med_OIIvsREF_major           = np.median(diff_OIIvsREF_major)
    std_OIIvsREF_major           = np.std(diff_OIIvsREF_major)
    diff_WLvsCubeWCSMod_major    = WLmajoraxes-CubeWCSModmajoraxes
    med_WLvsCubeWCSMod_major     = np.median(diff_WLvsCubeWCSMod_major)
    std_WLvsCubeWCSMod_major     = np.std(diff_WLvsCubeWCSMod_major)
    diff_WLvsACSMOD_major        = WLmajoraxes-ACSMODmajoraxes
    med_WLvsACSMOD_major         = np.median(diff_WLvsACSMOD_major)
    std_WLvsACSMOD_major         = np.std(diff_WLvsACSMOD_major)
    diff_WLvsREF_major           = WLmajoraxes-REFmajoraxes
    med_WLvsREF_major            = np.median(diff_WLvsREF_major)
    std_WLvsREF_major            = np.std(diff_WLvsREF_major)

    # ---------------------------------------------------------------------------------------------------------
    # ---------------- PRINTING ----------------
    if verbose:
        print('\n - dat_CubeWCSModel axis estimates: ')
        print('   median(whitelight - OIIband                 ) = '+
              str(med_WLvsOII_minor)+'arcsec +/- '+str(std_WLvsOII_minor))
        print(' ')
        print('   median(OIIband    - ACS model in CUBE WCS) = '+
              str(med_OIIvsCubeWCSMod_minor)+'arcsec +/- '+str(std_OIIvsCubeWCSMod_minor))
        print('   median(OIIband    - ACS model               ) = '+
              str(med_OIIvsACSMOD_minor)+'arcsec +/- '+str(std_OIIvsACSMOD_minor))
        print('   median(OIIband    - ACS reference image     ) = '+
              str(med_OIIvsREF_minor)+'arcsec +/- '+str(std_OIIvsREF_minor))
        print(' ')
        print('   median(White Light Img - ACS model in CUBE WCS) = '+
              str(med_WLvsCubeWCSMod_minor)+'arcsec +/- '+str(std_WLvsCubeWCSMod_minor))
        print('   median(White Light Img - ACS model               ) = '+
              str(med_WLvsACSMOD_minor)+'arcsec +/- '+str(std_WLvsACSMOD_minor))
        print('   median(White Light Img - ACS reference image     ) = '+
              str(med_WLvsREF_minor)+'arcsec +/- '+str(std_WLvsREF_minor))

        print('\n - Major axis estimates: ')
        print('   median(whitelight - OIIband                 ) = '+
              str(med_WLvsOII_major)+'arcsec +/- '+str(std_WLvsOII_major))
        print(' ')
        print('   median(OIIband    - ACS model in CUBE WCS) = '+
              str(med_OIIvsCubeWCSMod_major)+'arcsec +/- '+str(std_OIIvsCubeWCSMod_major))
        print('   median(OIIband    - ACS model               ) = '+
              str(med_OIIvsACSMOD_major)+'arcsec +/- '+str(std_OIIvsACSMOD_major))
        print('   median(OIIband    - ACS reference image     ) = '+
              str(med_OIIvsREF_major)+'arcsec +/- '+str(std_OIIvsREF_major))
        print(' ')
        print('   median(White Light Img - ACS model in CUBE WCS) = '+
              str(med_WLvsCubeWCSMod_major)+'arcsec +/- '+str(std_WLvsCubeWCSMod_major))
        print('   median(White Light Img - ACS model               ) = '+
              str(med_WLvsACSMOD_major)+'arcsec +/- '+str(std_WLvsACSMOD_major))
        print('   median(White Light Img - ACS reference image     ) = '+
              str(med_WLvsREF_major)+'arcsec +/- '+str(std_WLvsREF_major))

    # ---------------------------------------------------------------------------------------------------------
    # ---------------- PLOTTING ----------------
    histvals = ['WLminoraxes','OIIminoraxes','CubeWCSModminoraxes','ACSMODminoraxes','REFminoraxes',
                'WLmajoraxes','OIImajoraxes','CubeWCSModmajoraxes','ACSMODmajoraxes','REFmajoraxes']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = plotdir+'histograms.pdf'
    # if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=(7,4))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.05, right=0.98, bottom=0.15, top=0.98)
    Fsize  = 10
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for hh, histlabel in enumerate(histvals):
        plt.hist(eval(histvals[hh]),bins=np.arange(0,2.2,0.06),label=histlabel.replace('_','\_'),histtype='step')

    plt.plot([MUSEpix2arcsec,MUSEpix2arcsec],[0,60],'--k',lw=lthick,label='MUSE pixel size')
    plt.plot([MUSEpix2arcsec*1./2.355,MUSEpix2arcsec*1./2.355],[0,60],
             linestyle='--',color='gray',lw=lthick,label='MUSE min. $\sigma$ allowed')

    plt.plot([ACSpix2arcsec,ACSpix2arcsec],[0,60],':k',lw=lthick,label='HST/ACS pixel size')
    plt.plot([ACSpix2arcsec*1./2.355,ACSpix2arcsec*1./2.355],[0,60],
             linestyle=':',color='gray',lw=lthick,label='HST min. $\sigma$ allowed')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel("Gaussian $\sigma$ [$''$]", fontsize=Fsize)
    leg = plt.legend(fancybox=True, loc='upper right',prop={'size':Fsize},ncol=1,numpoints=1,
                     bbox_to_anchor=(0.99, 0.99))  # add the legend
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xvals  = ['WLminoraxes','OIIminoraxes','OIIminoraxes','OIIminoraxes',
              'WLminoraxes','WLminoraxes','WLminoraxes','WLmajoraxes',
              'OIImajoraxes','OIImajoraxes','OIImajoraxes','WLmajoraxes',
              'WLmajoraxes','WLmajoraxes']
    yvals  = ['diff_WLvsOII_minor','diff_OIIvsCubeWCSMod_minor','diff_OIIvsACSMOD_minor','diff_OIIvsREF_minor',
              'diff_WLvsCubeWCSMod_minor','diff_WLvsACSMOD_minor','diff_WLvsREF_minor','diff_WLvsOII_major',
              'diff_OIIvsCubeWCSMod_major','diff_OIIvsACSMOD_major','diff_OIIvsREF_major','diff_WLvsCubeWCSMod_major',
              'diff_WLvsACSMOD_major','diff_WLvsREF_major']
    for pp in xrange(len(xvals)):
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plotname = plotdir+'axesVSdiff_'+yvals[pp]+'.pdf'
        # if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        fig = plt.figure(figsize=(7,4))
        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.10, right=0.98, bottom=0.15, top=0.98)
        Fsize  = 10
        lthick = 1
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.plot(eval(xvals[pp]),eval(yvals[pp]),'.r')#,label='S/N LSDcat scaled by '+str(scalefactor))

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # if logaxes:
        #     axisrange = [2,100]
        # else:
        #     axisrange = [0,70]
        plt.plot([np.min(eval(xvals[pp])),np.max(eval(xvals[pp]))],[0,0],'--k',lw=lthick)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.xlabel(xvals[pp].replace('_','\_')+" [$''$]", fontsize=Fsize)
        plt.ylabel(yvals[pp].replace('_','\_')+" [$''$]", fontsize=Fsize)
        # plt.ylim(axisrange)
        # plt.xlim(axisrange)
        # if logaxes:
        if 'major' in xvals[pp]:
            plt.xscale('log')
            plt.yscale('log')
        # leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1,
        #                  bbox_to_anchor=(0.01, 0.99))  # add the legend
        # leg.get_frame().set_alpha(0.7)
        print(' - Saving plot to '+plotname)
        plt.savefig(plotname)
        fig.clf()
        plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xvals  = ['WLminoraxes','OIIminoraxes','OIIminoraxes','OIIminoraxes',
              'WLminoraxes','WLminoraxes','WLminoraxes','WLmajoraxes',
              'OIImajoraxes','OIImajoraxes','OIImajoraxes','WLmajoraxes',
              'WLmajoraxes','WLmajoraxes',
              'CubeWCSModminoraxes','CubeWCSModminoraxes','ACSMODminoraxes',
              'CubeWCSModmajoraxes','CubeWCSModmajoraxes','ACSMODmajoraxes']
    yvals  = ['OIIminoraxes','CubeWCSModminoraxes','ACSMODminoraxes','REFminoraxes',
              'CubeWCSModminoraxes','ACSMODminoraxes','REFminoraxes','OIImajoraxes',
              'CubeWCSModmajoraxes','ACSMODmajoraxes','REFmajoraxes','CubeWCSModmajoraxes',
              'ACSMODmajoraxes','REFmajoraxes',
              'REFminoraxes','ACSMODminoraxes','REFminoraxes',
              'REFmajoraxes','ACSMODmajoraxes','REFmajoraxes']

    for pp in xrange(len(xvals)):
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plotname = plotdir+'axes_'+xvals[pp]+'_vs_'+yvals[pp]+'.pdf'
        # if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        fig = plt.figure(figsize=(7,4))
        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.10, right=0.98, bottom=0.15, top=0.98)
        Fsize  = 10
        lthick = 1
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.plot(eval(xvals[pp]),eval(yvals[pp]),'.r')#,label='S/N LSDcat scaled by '+str(scalefactor))

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # if logaxes:
        #     axisrange = [2,100]
        # else:
        #     axisrange = [0,70]
        plt.plot([np.min(eval(xvals[pp])),np.max(eval(xvals[pp]))],[np.min(eval(xvals[pp])),np.max(eval(xvals[pp]))],
                 '--k',lw=lthick)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.xlabel(xvals[pp].replace('_','\_')+" [$''$]", fontsize=Fsize)
        plt.ylabel(yvals[pp].replace('_','\_')+" [$''$]", fontsize=Fsize)
        # plt.ylim(axisrange)
        # plt.xlim(axisrange)
        # if logaxes:
        if 'major' in xvals[pp]:
            plt.xscale('log')
            plt.yscale('log')
        # leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1,
        #                  bbox_to_anchor=(0.01, 0.99))  # add the legend
        # leg.get_frame().set_alpha(0.7)
        print(' - Saving plot to '+plotname)
        plt.savefig(plotname)
        fig.clf()
        plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def MW102009072Guo9640_WhiteLightImage(outputdir,overwrite=True,verbose=True):
    """
    Function to generate white light image for the Star-Galaxy blend of MW=102009072/Guo=9640

    --- INPUT ---
    outputdir       Directory to contain white light images

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    outputdir = '/Users/kschmidt/work/publications/TDOSE/TDOSEexampleruns/'
    tsu.MW102009072Guo9640_WhiteLightImage(outputdir)

    """
    cubedir   = '/Volumes/DATABCKUP1/MUSE-Wide/DATACUBES/'
    pointing  = 'candels-cdfs-02'
    # from /Users/kschmidt/work/publications/TDOSE/TDOSEexampleruns/tdose_sourcecats/catalog_photometry_candels-cdfs-02_tdose_sourcecat.txt
    ra        = [53.0624529]
    Dra       = 4.0 # arcsec
    dec       = [-27.8225088]
    Ddec      = Dra
    objid     = str(102009072)
    dcube     = cubedir+'DATACUBE_'+pointing+'_v1.0_dcbgc_effnoised.fits'
    datacube  = glob.glob(dcube)

    print(' ---- Extracting white light image and cube for '+objid+' ---- ')
    wcenter   = [[7050]]
    dwave     = [[2250]]
    name      = [objid+'_whitelight']
    mu.create_narrowband_subcube(datacube[0],ra,dec,Dra,Ddec,wcenter,dwave,outputdir,names=name,clobber=overwrite)
    imgwhite  = glob.glob(outputdir+'*'+name[0]+'*narrowbandimage*fits')

    print(' ---- Extracting OIII image and cube for '+objid+' ---- ')
    redshift  = 0.3377
    linewave  = 5007.0
    wcenter   = [[linewave*(redshift+1.0)]]
    dwave     = [[500.0/299792.0 * linewave * (redshift+1.0)]]# narrowband width is 2x500=1000 km/s rest-frame
    name      = [objid+'_OIII5007Width1000kmsRest']
    mu.create_narrowband_subcube(datacube[0],ra,dec,Dra,Ddec,wcenter,dwave,outputdir,names=name,clobber=overwrite)
    imgOIII   = glob.glob(outputdir+'*'+name[0]+'*narrowbandimage*fits')

    refimgcut = '/Users/kschmidt/work/publications/TDOSE/TDOSEexampleruns/181008_modelimg/tdose_cutouts/acs_814w_candels-cdfs-02_cut_v1.0_id9640_cutout4p0x4p0arcsec.fits'

    ds9cmd = 'ds9 -scale mode minmax '+\
             imgwhite[0]+' '+imgOIII[0]+' '+refimgcut+\
             ' -lock frame wcs -tile grid layout 4 1 &'
    print(ds9cmd)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def print_TDOSEsetup_PSFparameter(fields,printresults=True):
    """
    Function loading main PSF catalog and printing the 3 paramaters needed for the TDOSE setup files

    --- INPUT ---
    fields       list of fields to print PSF parameters for

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    fields   = ['cdfs-02','cdfs-03','cdfs-04','cdfs-04','cdfs-04','cdfs-42']
    psfparam = tsu.print_TDOSEsetup_PSFparameter(fields)

    """
    PSFcat = '/Users/kschmidt/work/MUSE/MUSEWide_PSFs/MUSEWide_PSF_Catalog_180723.fits'
    PSFdat = afits.open(PSFcat)[1].data

    psfparam = np.zeros([len(fields)],
                        dtype=[('field', 'S20'),('psf_FWHMp0', 'f4'), ('psf_FWHMp1', 'f4'), ('psf_FWHMp2', 'f4')])
    for ff, field in enumerate(fields):
        psfent = np.where(PSFdat['field'] == field)[0]

        if len(psfent) != 1:
            print(' WARNING: Found '+str(len(psfent))+' entries in PSF catalog for field = '+str(field))
        else:
            psfparam['field'][ff]      = field
            psfparam['psf_FWHMp0'][ff] = PSFdat['p0_sel_g'][psfent]
            psfparam['psf_FWHMp1'][ff] = PSFdat['p1_sel_g'][psfent]
            psfparam['psf_FWHMp2'][ff] = 7050.00

    if printresults:
        print('#            field              psf_FWHMp0            psf_FWHMp1              psf_FWHMp2 ')
        for pp in xrange(len(psfparam)):
            pstr = str("%20s" % psfparam['field'][pp])+'  '+\
                   str("%20s" % psfparam['psf_FWHMp0'][pp])+'  '+\
                   str("%20s" % psfparam['psf_FWHMp1'][pp])+'  '+\
                   str("%20s" % psfparam['psf_FWHMp2'][pp])
            print(pstr)

    return psfparam

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def OIIemitters_minosDownloadCMDs(outdir='/Users/kschmidt/work/TDOSE/OIIemitter_MWDR1spectra/'):
    """
    Function generating text file with download commands for minos.aip.de to grab MW DR1 spectra.

    --- INPUT ---


    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    tsu.OIIemitters_minosDownloadCMDs()

    """
    modeldir     = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/galfit_wrapper_results_final_181015/'
    fitscatalog  = modeldir+'OIIemitter_selection_23LTm814LT24_SkelsepLT0p3_Nobj153.fits'

    objdat = afits.open(fitscatalog)[1].data

    outcmdfile = outdir+'minos_downloadcommands.txt'
    fout = open(outcmdfile,'w')
    fout.write('# Command for grabbing spectra of OII emitters in:\n')
    fout.write('# '+fitscatalog+'\n')
    fout.write('# from minos.aip.de. These are the spectra in the MW DR1 at https://musewide.aip.de/dr1\n')

    for ii, mwid in enumerate(objdat['UNIQUE_ID']):
        fout.write('# >>>>>>>>>> SPECTRA FOR MUSE-Wide '+str(mwid)+' <<<<<<<<<<\n')

        ELcmd = 'scp kasper@minos.aip.de:/store/vo/musewide/musewide/emission_spectra_dr1/*'+\
                str(mwid)+'*  '+outdir

        TDcmd = 'scp kasper@minos.aip.de:/store/vo/musewide/musewide/guo_photometric_catalog_spectra/*'+\
                str("%.10d" % objdat['GUO_ID'][ii])+'*  '+outdir

        if str(mwid) == '116004061': # Ignore download of duplicate object
            ELcmd = '# '+ELcmd
            TDcmd = '# '+TDcmd

        fout.write(ELcmd+'\n')
        fout.write(TDcmd+'\n')

    fout.close()
    print(' - Wrote commands to:\n   '+outcmdfile)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def LSDCatify_spectrum(specfile,plot_comparison=False,xrange=None,yranges=None,verbose=True):
    """
    Apply matrix filtering to spectrum, to make the S/N values comparable to the LSDCat detections

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu

    specdir  = '/Volumes/DATABCKUP1/TDOSEextractions/181016_MWDR1_OIIemitters/tdose_spectra/'
    specfile = specdir+'tdose_spectrum_modelimg_0101014030-0101014030.fits'
    specfile = specdir+'tdose_spectrum_modelimg_0103024094-0103024094.fits'
    wave, flux_filtered, variance_filtered = tsu.LSDCatify_spectrum(specfile,plot_comparison='/Users/kschmidt/work/MUSE/notes/181119/LSDCatify_comparison.pdf',xrange=[6150,6350],yranges=[[-220,320],[-3,13]])

    """
    spechdu = afits.open(specfile)
    specdat = spechdu[1].data
    # if verbose: print(' - Convert spectrum to observed frame (if already observed frame set "redshift=0.0" ')

    if 'wave' in spechdu[1].columns.names:
        wavevec = specdat['wave']
    else:
        wavevec = specdat['wave_air']

    if 'flux' in spechdu[1].columns.names:
        fluxvec = specdat['flux']
    else:
        fluxvec = specdat['FLUX']

    if 'fluxerror' in spechdu[1].columns.names:
        fluxerrvec = specdat['fluxerror']
    else:
        fluxerrvec = specdat['FLUXERR']

    wave      = wavevec #* (1.0 + redshift)
    flux      = fluxvec
    variance  = fluxerrvec ** 2.0

    if verbose: print(' - Defining the filter matrix for the LDSCat filtering ')
    velocity       = 250.00 # FWHM[km/s] cf. Urrutia et al. MWDR1 section 4.1
    lambda_start   = wave[0]
    cdelt          = np.median(np.diff(wave))
    lMax           = len(wave)
    filter_matrix  = tsu.create_filter_matrix_vel(velocity,lambda_start=lambda_start,cdelt=cdelt,lMax=lMax)

    if verbose: print(' - Inserting flux and variance into dummy cube to filter spectrum and variance ')
    fluxcube               = np.zeros([lMax,1])
    fluxcube[:,0]          = flux
    variancecube           = np.zeros([lMax,1])
    variancecube[:,0]      = variance
    fluxcube_filtered      = tsu.filter_spectrum(filter_matrix,fluxcube)
    variancecube_filtered  = tsu.filter_spectrum(filter_matrix**2,variancecube)
    flux_filtered          = fluxcube_filtered[:,0]
    variance_filtered      = variancecube_filtered[:,0]

    if plot_comparison:
        figsize  = (5, 5)
        fontsize = 12
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plotname = plot_comparison
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        fig = plt.figure(figsize=figsize)
        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
        Fsize  = fontsize
        lthick = 1
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        plt.title('Gaussian filter with FWHM = '+str(velocity)+'km/s applied',fontsize=Fsize)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.plot(wave,flux,'-k',alpha=1.0,lw=lthick,label='Intrinsic')
        plt.fill_between(wave,flux-np.sqrt(variance),flux+np.sqrt(variance), alpha=0.4, label='+/-sqrt(variance)',
                         color='black')

        #plt.plot(wave,flux_filtered[:-10],'-g',alpha=1.0,lw=lthick,label='Filtered[:-10]')
        plt.plot(wave,flux_filtered[8:-8],'-r',alpha=1.0,lw=lthick,label='Filtered[8:-8]')
        plt.fill_between(wave,flux_filtered[8:-8]-np.sqrt(variance_filtered[8:-8]),
                         flux_filtered[8:-8]+np.sqrt(variance_filtered[8:-8]),
                         alpha=0.4, label='+/-sqrt(variance)', color='red')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.xlabel('wavelength', fontsize=Fsize)
        plt.ylabel('flux', fontsize=Fsize)
        leg = plt.legend(fancybox=True, loc='lower center',prop={'size':Fsize},ncol=2,numpoints=1)
        leg.get_frame().set_alpha(0.7)
        if xrange is not None:
            plt.xlim(xrange)
        if yranges is not None:
            plt.ylim(yranges[0])

        print(' - Saving plot to '+plotname)
        plt.savefig(plotname)
        fig.clf()
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plotname = plot_comparison.replace('.pdf','_S2N.pdf')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        fig = plt.figure(figsize=figsize)
        fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
        Fsize  = fontsize
        lthick = 1
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        plt.title('Gaussian filter with FWHM = '+str(velocity)+'km/s applied',fontsize=Fsize)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.plot(wave,flux/np.sqrt(variance),'-k',alpha=1.0,lw=lthick,label='Intrinsic')
        s2n_filtered = flux_filtered/np.sqrt(variance_filtered)
        #plt.plot(wave,s2n_filtered[:-10],'-g',alpha=1.0,lw=lthick,label='Filtered[:-10]')
        plt.plot(wave,s2n_filtered[8:-8],'-r',alpha=1.0,lw=lthick,label='Filtered[8:-8]')
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.xlabel('wavelength', fontsize=Fsize)
        plt.ylabel('flux', fontsize=Fsize)
        leg = plt.legend(fancybox=True, loc='lower center',prop={'size':Fsize},ncol=2,numpoints=1)
        leg.get_frame().set_alpha(0.7)
        if xrange is not None:
            plt.xlim(xrange)
        if yranges is not None:
            plt.ylim(yranges[1])

        print(' - Saving plot to '+plotname)
        plt.savefig(plotname)
        fig.clf()
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    return wave, flux_filtered[8:-8], variance_filtered[8:-8]

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def rename_MWDR1specs(objcatalog,speclist,outdir,idswap='GUO_ID',verbose=True):
    """
    Function to replace Guo ID with the MUSE-Wide id in spectra extracted with TDOSE based on Guo IDs (MW DR1)

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu
    import glob

    galfitdir      = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/galfit_wrapper_results_final_181015/'
    objcatalog     = galfitdir+'OIIemitter_selection_23LTm814LT24_SkelsepLT0p3_Nobj153.fits'

    speclist       = glob.glob('/Users/kschmidt/work/TDOSE/OIIemitter_MWDR1spectra/tdose_spectrum_*.fits')

    outdir         = '/Users/kschmidt/work/TDOSE/OIIemitter_MWDR1spectra/renamed_TDOSEGuoSpec/'

    tsu.rename_MWDR1specs(objcatalog,speclist,outdir,idswap='GUO_ID')

    """
    objdat   = afits.open(objcatalog)[1].data
    MWids    = objdat['UNIQUE_ID']
    swapids  = objdat[idswap]

    if verbose: print(' - Copying spectra to '+outdir+'\n   replaceing the '+idswap+
                      ' ids with the MUSE-Wide ids in the catalog\n'+objcatalog)
    for spec in speclist:
        specid  = int(spec.split('_')[-1].split('.fit')[0])
        swapent = np.where(swapids.astype('int') == specid)[0]
        MWid    = MWids[swapent]

        if len(swapent) >1:
            print ('\n\nWARNING multiple matches to swapid '+str(specid)+
                   ' using first MUSE-Wide ID for renaming: '+str(MWids[swapent]))
            MWid    = MWid[0]
        newname = outdir+spec.split('/')[-1].replace(str("%.10d" % specid),str("%.10d" % MWid))
        #if verbose: print(' - copying '+spec+' to '+newname)
        shutil.copyfile(spec,newname)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_maxvalues(tablename,specdirstrings,objcatalog,IDandRedshiftCol=['UNIQUE_ID','Z'],catalog=None,
                  colnames='SN',waveline=3726.0,dwave=10,namerefs='default',overwrite=False,verbose=True):
    """
    Function to pull out the max fluxes and S/N for a given range around an emission feature (line)
    for comparison of different extraction methods.

    --- INPUT ---
    tablename           Path and name of fits binary table to contain the results
    specdirstrings      List of strings to be used for collecting the spectra.
                        III will be replaced with IDs from objcatalog.
    objcatalog          Catalog containing IDs and redshifts of objects to include in table
    IDandRedshiftCol    Name of columns containing the ID and redshifts of objects in "objcatalog"
    catalog             Catalog with further information (e.g. LSDCat information) to add to output fits table
    colnames            Names of columns in "catalog" to include in output fits table
    waveline            Rest-frame wavelength of emission feature (emission line) to extract value from.
    dwave               Rest-frame wavelength size of region to extract values from (waveline +/- dwave).
    namerefs            Reference names to use for column name of the spectra found when searching the
                        "specdirstrings".
    overwrite           Overwrite existing output fits table
    verbose             Toggle verbosity

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu

    specdir        = '/Volumes/DATABCKUP1/TDOSEextractions/181016_MWDR1_OIIemitters/tdose_spectra/'
    MWDR1dir       = '/Users/kschmidt/work/TDOSE/OIIemitter_MWDR1spectra/'
    specdirstrings = [specdir+'tdose_spectrum_modelimg_0III-0III.fits',specdir+'tdose_spectrum_gauss_0III-0III.fits',MWDR1dir+'spectrum_III.fits',MWDR1dir+'renamed_TDOSEGuoSpec/tdose_spectrum_*III.fits']

    galfitdir      = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/galfit_wrapper_results_final_181015/'
    objcatalog     = galfitdir+'OIIemitter_selection_23LTm814LT24_SkelsepLT0p3_Nobj153.fits'

    tablename      = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/test181116_maxvalues.fits'
    tsu.get_maxvalues(tablename,specdirstrings,objcatalog,IDandRedshiftCol=['UNIQUE_ID','Z'],catalog=objcatalog, colnames=['UNIQUE_ID','SN'],waveline=3726.0,dwave=10,namerefs=['modelimg','gauss','mwdr1el','mwdr1td'])

    """
    Nstrings = len(specdirstrings)
    objdat   = afits.open(objcatalog)[1].data
    obj_ids  = objdat[IDandRedshiftCol[0]]
    obj_zs   = objdat[IDandRedshiftCol[1]]
    Nobj     = len(obj_ids)
    if verbose: print(' - Found '+str(Nobj)+' objects in the provided catalog to pull out spectral values for ')

    if catalog is None:
        Nxtracol = 0
    else:
        Nxtracol    = len(colnames)
        catalogdat  = afits.open(catalog)[1].data

    # outarray will have columns:
    # ID, z, [f,ferr,fs2n,s2n]*Nstrings, [f_filtered,ferr_filtered,fs2n_filt,s2n_filtered]*Nstrings, colnames
    Ncolperspec = 8
    outarray    = np.zeros([Nobj,2 + Nstrings*8 + Nxtracol])

    if verbose: print(' - Looping over objects to fill output array with extracted and (LSDCat) filtered values')
    for oo, objid in enumerate(obj_ids):
        #if oo < 605: continue
        infostr = ' - Getting info for object '+str(objid)+' ( file '+\
                  str("%3.f" % (oo+1))+' / '+str("%3.f" % len(obj_ids)+')')
        sys.stdout.write("%s\r" % infostr)
        sys.stdout.flush()

        waveobs  = (1+obj_zs[oo])*waveline
        outarray[oo,0] = objid
        outarray[oo,1] = obj_zs[oo]

        for ss, sstring in enumerate(specdirstrings):
            spec  = glob.glob(sstring.replace('III',str(objid)))
            Nspec = len(spec)

            if Nspec == 1:
                specfile = spec[0]
                spechdu  = afits.open(specfile)
                specdat  = spechdu[1].data
                if 's2n' in spechdu[1].columns.names:
                    s2nvec = specdat['s2n']
                else:
                    s2nvec = specdat['flux']/specdat['fluxerr']

                if 'wave' in spechdu[1].columns.names:
                    wavevec = specdat['wave']
                else:
                    wavevec = specdat['wave_air']

                if 'flux' in spechdu[1].columns.names:
                    fluxvec = specdat['flux']
                else:
                    fluxvec = specdat['FLUX']

                if 'fluxerror' in spechdu[1].columns.names:
                    fluxerrvec = specdat['fluxerror']
                else:
                    fluxerrvec = specdat['FLUXERR']

                if len(fluxvec[np.isfinite(fluxvec)]) == 0:
                    maxF, maxFerr, maxFs2n, maxs2n, maxF_filt, maxFerr_filt, maxFs2n_filt, maxs2n_filt = [np.nan]*8
                    if verbose: print ('\n   WARNING: No finite flux spectrum pixels for \n   '+specfile+
                                       '\n   Returning NaNs for the "max values" and "filtered max values"')
                elif (waveobs < wavevec[0]) or (waveobs > wavevec[-1]):
                    maxF, maxFerr, maxFs2n, maxs2n, maxF_filt, maxFerr_filt, maxFs2n_filt, maxs2n_filt = [np.nan]*8
                    if verbose: print ('\n   WARNING: Line wavelength ('+str(waveobs)+'A @ z='+
                                       str(obj_zs[oo])+') for object '+str(objid)+' outside spectral range '
                                       '\n   Returning NaNs for the "max values" and "filtered max values"')
                else:
                    maxF, maxFerr, maxFs2n, maxs2n = \
                        tsu.extractMaxValues(wavevec,fluxvec,fluxerrvec,s2nvec,waveobs,dwave)

                    wave_filt, flux_filt, variance_filt = tsu.LSDCatify_spectrum(specfile,verbose=False)
                    s2n_filt       = flux_filt / np.sqrt(variance_filt)

                    maxF_filt, maxFerr_filt, maxFs2n_filt, maxs2n_filt = \
                        tsu.extractMaxValues(wave_filt,flux_filt,np.sqrt(variance_filt),s2n_filt,waveobs,dwave)

                outarray[oo,2+ss*Ncolperspec] = maxF
                outarray[oo,3+ss*Ncolperspec] = maxFerr
                outarray[oo,4+ss*Ncolperspec] = maxFs2n
                outarray[oo,5+ss*Ncolperspec] = maxs2n

                outarray[oo,6+ss*Ncolperspec] = maxF_filt
                outarray[oo,7+ss*Ncolperspec] = maxFerr_filt
                outarray[oo,8+ss*Ncolperspec] = maxFs2n_filt
                outarray[oo,9+ss*Ncolperspec] = maxs2n_filt
            else:
                if verbose: print('\n   WARNING: Found '+str(Nspec)+' spectra globbing for \n   '+
                                  str(sstring.replace('III',str(objid)))+
                                  '\n   for object '+str(objid)+' so storing NaNs ')
                outarray[oo,2:2+8*Nstrings] = [None]*8*Nstrings
            if colnames is not None:
                objcatent = np.where(catalogdat[IDandRedshiftCol[0]] == objid)[0]
                if len(objcatent) == 1:
                    outarray[oo,2+8*Nstrings:] = [catalogdat[cc][objcatent][0] for cc in colnames]
                else:
                    outarray[oo,2+8*Nstrings:] = [None]*len(colnames)
    if verbose: print('\n   ... done')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Saving output array to fits file')
    mainHDU = afits.PrimaryHDU()       # primary HDU

    c0 = afits.Column(name='id',                    format='D', unit='' , array=outarray[:,0])
    c1 = afits.Column(name='redshift',              format='D',   unit='' , array=outarray[:,1])

    collist = [c0,c1]

    for ss, sstring in enumerate(specdirstrings):
        c2 = afits.Column(name='max_flux_'+namerefs[ss],          format='D',   unit='' , array=outarray[:,2+ss*8])
        c3 = afits.Column(name='max_fluxerror_'+namerefs[ss],     format='D',   unit='' , array=outarray[:,3+ss*8])
        c4 = afits.Column(name='max_flux_s2n_'+namerefs[ss],      format='D',   unit='' , array=outarray[:,4+ss*8])
        c5 = afits.Column(name='max_s2n_'+namerefs[ss],           format='D',   unit='' , array=outarray[:,5+ss*8])
        c6 = afits.Column(name='max_flux_filt_'+namerefs[ss],     format='D',   unit='' , array=outarray[:,6+ss*8])
        c7 = afits.Column(name='max_fluxerror_filt_'+namerefs[ss],format='D',   unit='' , array=outarray[:,7+ss*8])
        c8 = afits.Column(name='max_flux_s2n_filt_'+namerefs[ss], format='D',   unit='' , array=outarray[:,8+ss*8])
        c9 = afits.Column(name='max_s2n_filt_'+namerefs[ss],      format='D',   unit='' , array=outarray[:,9+ss*8])

        collist = collist+[c2,c3,c4,c5,c6,c7,c8,c9]
    if colnames is not None:
        xtra_collist = []
        for cc, cname in enumerate(colnames):
            cX = afits.Column(name=cname, format='D', unit='' , array=outarray[:,2+Nstrings*8+cc])
            xtra_collist.append(cX)

        collist = collist + xtra_collist

    tbHDU   = afits.BinTableHDU.from_columns(collist)#, header=head)
    hdulist = afits.HDUList([mainHDU,tbHDU])
    hdulist.writeto(tablename, overwrite=overwrite)
    if verbose: print('   Saved to '+tablename)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def extractMaxValues(wave,flux,fluxerror,s2n,wavecenter,dwave):
    """
    Function to return the extracted max-values
    """
    wavediff    = np.abs(wave-wavecenter)
    waveent     = np.where(wavediff == np.min(wavediff))[0]

    if len(waveent[np.isfinite(waveent)]) == 0:
        maxF, maxFerr, maxFs2n, maxs2n, maxF_filt, maxFerr_filt, maxFs2n_filt, maxs2n_filt = [np.nan]*8
        if verbose: print ('\n   WARNING: No finite "waveent" entries so returning NaNs for the max values')
        maxF, maxFerr, maxFs2n, maxs2n = [np.nan]*4
    else:

        try:
            entmin = waveent[0]-dwave
            if entmin < 0: entmin = 0
            entmax = waveent[0]+dwave

            subdat_S2N  = s2n[entmin:entmax]
            subdat_F    = flux[entmin:entmax]
            fmaxent     = np.where(subdat_F == np.max(subdat_F))[0]
            subdat_Ferr = fluxerror[entmin:entmax][fmaxent]

            maxF        = np.max(subdat_F)
            maxFerr     = subdat_Ferr[0]
            maxFs2n     = np.max(subdat_F)/subdat_Ferr[0]
            maxs2n      = np.max(subdat_S2N)
        except:
            pdb.set_trace()
    return maxF, maxFerr, maxFs2n, maxs2n
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_maxvalues_OIIemitters(maxvalfitstable, outliercut=1.25, figsize=(5, 5), fontsize=12, logaxes=True, showids=False,
                   nstd = 3.0, # The number of sigma shown as region around for ODR fitted lines
                   namebase='/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/maxval_comparison_fitscat'):
    """
    Plotting comparison of extractions of max values.
    Based on tsu.OIIemitters_plotcomparisons() making this function obsolete.

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu

    fitstab  = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/OIIemitter_spectra_maxvalues_181116.fits'
    namebase = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/maxval_comparison_fitscat'
    tsu.plot_maxvalues_OIIemitters(fitstab,namebase=namebase,logaxes=False)

    """
    maxvaldat   = afits.open(maxvalfitstable)[1].data
    goodent     = np.isfinite(maxvaldat['max_s2n_modelimg'])
    maxvaldat   = maxvaldat[goodent]

    ids         = maxvaldat['id']
    selectedids = maxvaldat['id'] #[10233149., 10209072., 10212085.]
    # for id in ids:
    #     if id in selectedids:
    #         print(id)
    redshift    = maxvaldat['redshift']

    SNsersic    = maxvaldat['max_s2n_modelimg']
    Fsersic     = maxvaldat['max_flux_modelimg']
    Ferrsersic  = maxvaldat['max_fluxerror_modelimg']

    SNsersic_filtered    = maxvaldat['max_s2n_filt_modelimg']
    Fsersic_filtered     = maxvaldat['max_flux_filt_modelimg']
    Ferrsersic_filtered  = maxvaldat['max_fluxerror_filt_modelimg']

    SNgauss    = maxvaldat['max_s2n_gauss']
    Fgauss     = maxvaldat['max_flux_gauss']
    Ferrgauss  = maxvaldat['max_fluxerror_gauss']

    SNgauss_filtered    = maxvaldat['max_s2n_filt_gauss']
    Fgauss_filtered     = maxvaldat['max_flux_filt_gauss']
    Ferrgauss_filtered  = maxvaldat['max_fluxerror_filt_gauss']

    SNmwdr1td    = maxvaldat['max_s2n_mwdr1td']
    Fmwdr1td     = maxvaldat['max_flux_mwdr1td']
    Ferrmwdr1td  = maxvaldat['max_fluxerror_mwdr1td']

    SNmwdr1td_filtered    = maxvaldat['max_s2n_filt_mwdr1td']
    Fmwdr1td_filtered     = maxvaldat['max_flux_filt_mwdr1td']
    Ferrmwdr1td_filtered  = maxvaldat['max_fluxerror_filt_mwdr1td']

    SNmwdr1el    = maxvaldat['max_s2n_mwdr1el']
    Fmwdr1el     = maxvaldat['max_flux_mwdr1el']
    Ferrmwdr1el  = maxvaldat['max_fluxerror_mwdr1el']

    SNmwdr1el_filtered    = maxvaldat['max_s2n_filt_mwdr1el']
    Fmwdr1el_filtered     = maxvaldat['max_flux_filt_mwdr1el']
    Ferrmwdr1el_filtered  = maxvaldat['max_fluxerror_filt_mwdr1el']


    IDlsdcat    = maxvaldat['UNIQUE_ID']
    SNlsdcat    = maxvaldat['SN']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_ODRfit2data_fluxtd.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    odrfit_Fsersic_Fmwdr1td = kbs.fit_function_to_data_with_errors_on_both_axes(Fsersic,Fmwdr1td,Ferrsersic,Ferrmwdr1td,
                                                                              fitfunction='linear',
                                                                              plotresults=plotname)
    print(' - A bit of stats on the input data for the ODR fit:')
    print('   -> np.median(np.asarray(Fsersic)/np.asarray(Fmwdr1td))  = '+
          str(np.median(np.asarray(Fsersic)/np.asarray(Fmwdr1td))))
    print('   -> np.std(np.asarray(Fsersic)/np.asarray(Fmwdr1td))  = '+
          str(np.std(np.asarray(Fsersic)/np.asarray(Fmwdr1td))))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_ODRfit2data_fluxel.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    odrfit_Fsersic_Fmwdr1el = kbs.fit_function_to_data_with_errors_on_both_axes(Fsersic,Fmwdr1el,Ferrsersic,Ferrmwdr1el,
                                                                              fitfunction='linear',
                                                                              plotresults=plotname)
    print(' - A bit of stats on the input data for the ODR fit:')
    print('   -> np.median(np.asarray(Fsersic)/np.asarray(Fmwdr1el))  = '+
          str(np.median(np.asarray(Fsersic)/np.asarray(Fmwdr1el))))
    print('   -> np.std(np.asarray(Fsersic)/np.asarray(Fmwdr1el))  = '+
          str(np.std(np.asarray(Fsersic)/np.asarray(Fmwdr1el))))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_ODRfit2data_SN.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    odrfit_SNsersic_SNlsdcat = kbs.fit_function_to_data_with_errors_on_both_axes(SNsersic_filtered,SNlsdcat,
                                                                                 np.asarray(SNlsdcat)*0.0+1.0,
                                                                                 np.asarray(SNlsdcat)*0.0+1.0,
                                                                                 fitfunction='linear',
                                                                                 plotresults=namebase+'_ODRfit2data_SN.pdf')

    print(' - A bit of stats on the input data for the ODR fit:')
    print('   -> np.median(np.asarray(SNsersic_filtered)/np.asarray(SNlsdcat))  = '+
          str(np.median(np.asarray(SNsersic_filtered)/np.asarray(SNlsdcat))))
    print('   -> np.std(np.asarray(SNsersic_filtered)/np.asarray(SNlsdcat))  = '+
          str(np.std(np.asarray(SNsersic_filtered)/np.asarray(SNlsdcat))))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_ODRfit2data_SNtd.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    odrfit_SNsersic_SNmwdr1td = kbs.fit_function_to_data_with_errors_on_both_axes(SNsersic,SNmwdr1td,
                                                                                 np.asarray(SNmwdr1td)*0.0+1.0,
                                                                                 np.asarray(SNmwdr1td)*0.0+1.0,
                                                                                 fitfunction='linear',
                                                                                 plotresults=plotname)
    print(' - A bit of stats on the input data for the ODR fit:')
    print('   -> np.median(np.asarray(SNsersic)/np.asarray(SNmwdr1td))  = '+
          str(np.median(np.asarray(SNsersic)/np.asarray(SNmwdr1td))))
    print('   -> np.std(np.asarray(SNsersic)/np.asarray(SNmwdr1td))  = '+
          str(np.std(np.asarray(SNsersic)/np.asarray(SNmwdr1td))))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_ODRfit2data_SNel.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    odrfit_SNsersic_SNmwdr1el = kbs.fit_function_to_data_with_errors_on_both_axes(SNsersic,SNmwdr1el,
                                                                                 np.asarray(SNmwdr1el)*0.0+1.0,
                                                                                 np.asarray(SNmwdr1el)*0.0+1.0,
                                                                                 fitfunction='linear',
                                                                                 plotresults=plotname)
    print(' - A bit of stats on the input data for the ODR fit:')
    print('   -> np.median(np.asarray(SNsersic)/np.asarray(SNmwdr1el))  = '+
          str(np.median(np.asarray(SNsersic)/np.asarray(SNmwdr1el))))
    print('   -> np.std(np.asarray(SNsersic)/np.asarray(SNmwdr1el))  = '+
          str(np.std(np.asarray(SNsersic)/np.asarray(SNmwdr1el))))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2N.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.plot(SNsersic,SNgauss,'ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if id in selectedids:
                if (SNsersic[ii]/SNgauss[ii] > outliercut):
                    plt.text(SNsersic[ii],SNgauss[ii],id,fontsize=Fsize-5,horizontalalignment='left')
                if (SNgauss[ii]/SNsersic[ii] > outliercut):
                    plt.text(SNsersic[ii],SNgauss[ii],id,fontsize=Fsize-5,horizontalalignment='right')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [2,100]
    else:
        axisrange = [0,70]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot([5,5],axisrange,'--k',lw=lthick)
    plt.plot(axisrange,[5,5],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('max(S/N) TDOSE multicomponent Gauss model', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    # leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1,
    #                  bbox_to_anchor=(0.01, 0.99))  # add the legend
    # leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2Nlsdcat.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    scalefactor = 1.5
    plt.plot(SNsersic,np.asarray(SNlsdcat)/scalefactor,'.r',label='S/N LSDcat scaled by '+str(scalefactor))

    plt.plot(SNsersic,SNlsdcat,'ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if id in selectedids:
                if (SNsersic[ii]/(SNlsdcat[ii]/scalefactor) > outliercut):
                    plt.text(SNsersic[ii],SNlsdcat[ii],id,fontsize=Fsize-5,horizontalalignment='left')
                if ((SNlsdcat[ii]/scalefactor)/SNsersic[ii] > outliercut):
                    plt.text(SNsersic[ii],SNlsdcat[ii],id,fontsize=Fsize-5,horizontalalignment='right')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [2,100]
    else:
        axisrange = [0,70]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot([5,5],axisrange,'--k',lw=lthick)
    plt.plot(axisrange,[5,5],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('S/N MW DR1 LSDcat measurements', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1,
                     bbox_to_anchor=(0.01, 0.99))  # add the legend
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2Nlsdcat_filtered.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    scalefactor = 1.0
    # plt.plot(SNsersic,np.asarray(SNlsdcat)/scalefactor,'.r',label='S/N LSDcat scaled by '+str(scalefactor))
    # plt.plot(SNsersic_filtered,SNlsdcat,'.g',label='LSDcatified TDOSE S/N')

    plt.plot(SNsersic_filtered,SNlsdcat,'ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if id in selectedids:
                if (SNlsdcat[ii]/SNsersic_filtered[ii] > 1.0): # showing IDs of objects where LSDCat is highest
                    print('   ID, S/N_LSDCat, S/N_LSDCat / S/N_TDOSEfiltered '+str(id)+', '+str(SNlsdcat[ii])+
                          ', '+str(SNlsdcat[ii]/SNsersic_filtered[ii]))
                    plt.text(SNsersic_filtered[ii],SNlsdcat[ii],id,fontsize=Fsize-5,horizontalalignment='right')

            # if (SNsersic_filtered[ii]/(SNlsdcat[ii]/scalefactor) > outliercut):
            #     plt.text(SNsersic_filtered[ii],SNlsdcat[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            # if ((SNlsdcat[ii]/scalefactor)/SNsersic_filtered[ii] > outliercut):
            #     plt.text(SNsersic_filtered[ii],SNlsdcat[ii],id,fontsize=Fsize-5,horizontalalignment='right')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # overplot ODR fit
    popt = odrfit_SNsersic_SNlsdcat.beta
    perr = odrfit_SNsersic_SNlsdcat.sd_beta
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    x_fit = np.linspace(min(SNsersic_filtered), max(SNsersic_filtered), 100)
    fit = kbs.odr_fitfunction_linear(popt, x_fit)
    fit_up = kbs.odr_fitfunction_linear(popt_up, x_fit)
    fit_dw= kbs.odr_fitfunction_linear(popt_dw, x_fit)

    plt.fill_between(x_fit, fit_up, fit_dw, alpha=0.25, label='3-sigma interval',color='red')
    plt.plot(x_fit, fit, 'r', lw=2, label='Best fit curve')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [2,100]
    else:
        axisrange = [0,70]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot([5,5],axisrange,'--k',lw=lthick)
    plt.plot(axisrange,[5,5],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT multicomponent Sersic model (filtered)', fontsize=Fsize)
    plt.ylabel('S/N MW DR1 LSDcat measurements', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_flux_mwdr1td.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.errorbar(Fsersic,Fmwdr1td, yerr=Ferrgauss, xerr=Ferrsersic,fmt='ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if id in selectedids:
                if (Fsersic[ii]/Fmwdr1td[ii] > outliercut):
                    plt.text(Fsersic[ii],Fmwdr1td[ii],id,fontsize=Fsize-5,horizontalalignment='left')
                if (Fmwdr1td[ii]/Fsersic[ii] > outliercut):
                    plt.text(Fsersic[ii],Fmwdr1td[ii],id,fontsize=Fsize-5,horizontalalignment='right')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # overplot ODR fit
    popt = odrfit_Fsersic_Fmwdr1td.beta
    perr = odrfit_Fsersic_Fmwdr1td.sd_beta
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    x_fit = np.linspace(min(Fsersic), max(Fsersic), 100)
    fit = kbs.odr_fitfunction_linear(popt, x_fit)
    fit_up = kbs.odr_fitfunction_linear(popt_up, x_fit)
    fit_dw= kbs.odr_fitfunction_linear(popt_dw, x_fit)

    plt.fill_between(x_fit, fit_up, fit_dw, alpha=0.25, label='3-sigma interval',color='red')
    plt.plot(x_fit, fit, 'r', lw=2, label='Best fit curve')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [100,4000]
    else:
        axisrange = [0,1500]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(F/[cgs]) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('max(F/[cgs]) MW DR1 TDOSE spectrum', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_flux_mwdr1el.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    scalefactor = 1/0.5156
    if scalefactor != 1.0:
        plt.plot(Fsersic,Fmwdr1el/scalefactor,'.r',label='Scaled PSF weighted flux') #/ '+str(scalefactor)

    plt.errorbar(Fsersic,Fmwdr1el, yerr=Ferrgauss, xerr=Ferrsersic,fmt='ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if id in selectedids:
                if (Fsersic[ii]/Fmwdr1el[ii] > outliercut):
                    plt.text(Fsersic[ii],Fmwdr1el[ii],id,fontsize=Fsize-5,horizontalalignment='left')
                if (Fmwdr1el[ii]/Fsersic[ii] > outliercut):
                    plt.text(Fsersic[ii],Fmwdr1el[ii],id,fontsize=Fsize-5,horizontalalignment='right')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # overplot ODR fit
    popt = odrfit_Fsersic_Fmwdr1el.beta
    perr = odrfit_Fsersic_Fmwdr1el.sd_beta
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    x_fit = np.linspace(min(Fsersic), max(Fsersic), 100)
    fit = kbs.odr_fitfunction_linear(popt, x_fit)
    fit_up = kbs.odr_fitfunction_linear(popt_up, x_fit)
    fit_dw= kbs.odr_fitfunction_linear(popt_dw, x_fit)

    plt.fill_between(x_fit, fit_up, fit_dw, alpha=0.25, label='3-sigma interval',color='red')
    plt.plot(x_fit, fit, 'r', lw=2, label='Best fit curve')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        xaxisrange = [100,2000]
        yaxisrange = [100,6000]
    else:
        xaxisrange = [0,1500]
        yaxisrange = [0,1500]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(F/[cgs]) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('max(F/[cgs]) MW DR1 PSF weighted spectrum', fontsize=Fsize)
    plt.ylim(yaxisrange)
    plt.xlim(xaxisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2N_mwdr1td.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.plot(SNsersic,SNmwdr1td,'ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if id in selectedids:
                if (SNsersic[ii]/SNmwdr1td[ii] > outliercut):
                    plt.text(SNsersic[ii],SNmwdr1td[ii],id,fontsize=Fsize-5,horizontalalignment='left')
                if (SNmwdr1td[ii]/SNsersic[ii] > outliercut):
                    plt.text(SNsersic[ii],SNmwdr1td[ii],id,fontsize=Fsize-5,horizontalalignment='right')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # overplot ODR fit
    popt = odrfit_SNsersic_SNmwdr1td.beta
    perr = odrfit_SNsersic_SNmwdr1td.sd_beta
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    x_fit = np.linspace(min(SNsersic), max(SNsersic), 100)
    fit = kbs.odr_fitfunction_linear(popt, x_fit)
    fit_up = kbs.odr_fitfunction_linear(popt_up, x_fit)
    fit_dw= kbs.odr_fitfunction_linear(popt_dw, x_fit)

    plt.fill_between(x_fit, fit_up, fit_dw, alpha=0.25, label='3-sigma interval',color='red')
    plt.plot(x_fit, fit, 'r', lw=2, label='Best fit curve')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [3,60]
    else:
        axisrange = [0,70]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('max(S/N) MW DR1 TDOSE spectrum', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2N_mwdr1el.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.plot(SNsersic,SNmwdr1el,'ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if id in selectedids:
                if (SNsersic[ii]/SNmwdr1el[ii] > outliercut):
                    plt.text(SNsersic[ii],SNmwdr1el[ii],id,fontsize=Fsize-5,horizontalalignment='left')
                if (SNmwdr1el[ii]/SNsersic[ii] > outliercut):
                    plt.text(SNsersic[ii],SNmwdr1el[ii],id,fontsize=Fsize-5,horizontalalignment='right')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # overplot ODR fit
    popt = odrfit_SNsersic_SNmwdr1el.beta
    perr = odrfit_SNsersic_SNmwdr1el.sd_beta
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    x_fit = np.linspace(min(SNsersic), max(SNsersic), 100)
    fit = kbs.odr_fitfunction_linear(popt, x_fit)
    fit_up = kbs.odr_fitfunction_linear(popt_up, x_fit)
    fit_dw= kbs.odr_fitfunction_linear(popt_dw, x_fit)

    plt.fill_between(x_fit, fit_up, fit_dw, alpha=0.25, label='3-sigma interval',color='red')
    plt.plot(x_fit, fit, 'r', lw=2, label='Best fit curve')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [3,60]
    else:
        axisrange = [0,70]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT multicomponent Sersic model', fontsize=Fsize)
    plt.ylabel('max(S/N) MW DR1 PSF weighted spectrum', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_maxvalues_LAEs(maxvalfitstable, outliercut=1.25, figsize=(5, 5), fontsize=12, logaxes=True, showids=False,
                        nstd = 3.0, # The number of sigma shown as region around for ODR fitted lines
                        namebase='/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/maxval_comparison_fitscat'):
    """
    Plotting comparison of extractions of max values.
    Based on tsu.OIIemitters_plotcomparisons() making this function obsolete.

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import tdosepublication_utilities as tsu

    fitstab  = '/Users/kschmidt/work/publications/TDOSE/fig/LAEcomparison/LAE_spectra_maxvalues_181120.fits'
    namebase = '/Users/kschmidt/work/publications/TDOSE/fig/LAEcomparison/LAE_maxval_comparison_fitscat'
    tsu.plot_maxvalues_LAEs(fitstab,namebase=namebase,logaxes=False)

    """
    maxvaldat   = afits.open(maxvalfitstable)[1].data
    goodent     = np.isfinite(maxvaldat['max_s2n_modelimg'])
    maxvaldat   = maxvaldat[goodent]

    ids         = maxvaldat['id']
    redshift    = maxvaldat['redshift']

    SNmodelimg    = maxvaldat['max_s2n_modelimg']
    Fmodelimg     = maxvaldat['max_flux_modelimg']
    Ferrmodelimg  = maxvaldat['max_fluxerror_modelimg']

    SNmodelimg_filtered    = maxvaldat['max_s2n_filt_modelimg']
    Fmodelimg_filtered     = maxvaldat['max_flux_filt_modelimg']
    Ferrmodelimg_filtered  = maxvaldat['max_fluxerror_filt_modelimg']

    # SNgauss    = maxvaldat['max_s2n_gauss']
    # Fgauss     = maxvaldat['max_flux_gauss']
    # Ferrgauss  = maxvaldat['max_fluxerror_gauss']
    #
    # SNgauss_filtered    = maxvaldat['max_s2n_filt_gauss']
    # Fgauss_filtered     = maxvaldat['max_flux_filt_gauss']
    # Ferrgauss_filtered  = maxvaldat['max_fluxerror_filt_gauss']

    SNaper      = maxvaldat['max_s2n_aperture0p5']
    Faper       = maxvaldat['max_flux_aperture0p5']
    Ferraper    = maxvaldat['max_fluxerror_aperture0p5']

    SNaper_filtered    = maxvaldat['max_s2n_filt_aperture0p5']
    Faper_filtered     = maxvaldat['max_flux_filt_aperture0p5']
    Ferraper_filtered  = maxvaldat['max_fluxerror_filt_aperture0p5']

    fwhm        = maxvaldat['fwhm_A']
    fwhm        = maxvaldat['fwhm_A_std']
    peaksep     = maxvaldat['peak_sep_kms']
    peakseperr  = maxvaldat['peak_sep_kms_std']
    EW0         = maxvaldat['EW_0']
    EW0err      = maxvaldat['EW_0_err']

    SNleadline  = maxvaldat['leadlineS2N']
    leadline    = []
    LAEinfodat  = afits.open('/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo.fits')[1].data
    for objid in ids:
        objent = np.where(LAEinfodat['id'] == objid)[0]
        leadline.append(LAEinfodat['leadline'][objent][0])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_ODRfit2data_flux_modelaper.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    odrfit_Fmodelimg_Faper = kbs.fit_function_to_data_with_errors_on_both_axes(Fmodelimg,Faper,Ferrmodelimg,Ferraper,
                                                                              fitfunction='linear',
                                                                              plotresults=plotname)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_ODRfit2data_SN_modelaper.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    odrfit_SNmodelimg_SNaper = kbs.fit_function_to_data_with_errors_on_both_axes(SNmodelimg,SNaper,
                                                                                 np.asarray(SNmodelimg)*0.0+1.0,
                                                                                 np.asarray(SNaper)*0.0+1.0,
                                                                                 fitfunction='linear',
                                                                                 plotresults=plotname)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_ODRfit2data_SN_modelfiltVSleadline.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    odrfit_SNmodelimgfilt_SNleadline = kbs.fit_function_to_data_with_errors_on_both_axes(SNmodelimg_filtered,SNleadline,
                                                                                         np.asarray(SNmodelimg)*0.0+1.0,
                                                                                         np.asarray(SNaper)*0.0+1.0,
                                                                                         fitfunction='linear',
                                                                                         plotresults=plotname)

    print(' - A bit of stats on the input data for the ODR fit:')
    print('   -> np.median(np.asarray(SNleadline)/np.asarray(SNmodelimg_filtered))  = '+
          str(np.median(np.asarray(SNleadline)/np.asarray(SNmodelimg_filtered))))
    print('   -> np.mean(np.asarray(SNleadline)/np.asarray(SNmodelimg_filtered))  = '+
          str(np.mean(np.asarray(SNleadline)/np.asarray(SNmodelimg_filtered))))
    print('   -> np.std(np.asarray(SNleadline)/np.asarray(SNmodelimg_filtered))  = '+
          str(np.std(np.asarray(SNleadline)/np.asarray(SNmodelimg_filtered))))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # plotname = namebase+'_ODRfit2data_Flux_modelfiltVSleadline.pdf'
    # if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # odrfit_Fmodelimgfilt_Fleadline = kbs.fit_function_to_data_with_errors_on_both_axes(Fmodelimg_filtered,Fleadline,
    #                                                                                      np.asarray(SNmodelimg)*0.0+1.0,
    #                                                                                      np.asarray(SNaper)*0.0+1.0,
    #                                                                                      fitfunction='linear',
    #                                                                                      plotresults=plotname)
    #
    # print(' - A bit of stats on the input data for the ODR fit:')
    # print('   -> np.median(np.asarray(Fleadline)/np.asarray(Fmodelimg_filtered))  = '+
    #       str(np.median(np.asarray(Fleadline)/np.asarray(Fmodelimg_filtered))))
    # print('   -> np.mean(np.asarray(Fleadline)/np.asarray(Fmodelimg_filtered))  = '+
    #       str(np.mean(np.asarray(Fleadline)/np.asarray(Fmodelimg_filtered))))
    # print('   -> np.std(np.asarray(Fleadline)/np.asarray(Fmodelimg_filtered))  = '+
    #       str(np.std(np.asarray(Fleadline)/np.asarray(Fmodelimg_filtered))))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_ODRfit2data_SN_modelVSleadline.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    odrfit_SNmodelimg_SNleadline = kbs.fit_function_to_data_with_errors_on_both_axes(SNmodelimg,SNleadline,
                                                                                     np.asarray(SNmodelimg)*0.0+1.0,
                                                                                     np.asarray(SNaper)*0.0+1.0,
                                                                                     fitfunction='linear',
                                                                                     plotresults=plotname)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2N_modelaper.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.plot(SNmodelimg,SNaper,'ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if (SNmodelimg[ii]/SNaper[ii] > outliercut):
                plt.text(SNmodelimg[ii],SNaper[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            if (SNaper[ii]/SNmodelimg[ii] > outliercut):
                plt.text(SNmodelimg[ii],SNaper[ii],id,fontsize=Fsize-5,horizontalalignment='right')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # overplot ODR fit
    popt = odrfit_SNmodelimg_SNaper.beta
    perr = odrfit_SNmodelimg_SNaper.sd_beta
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    x_fit = np.linspace(min(SNmodelimg), max(SNmodelimg), 100)
    fit = kbs.odr_fitfunction_linear(popt, x_fit)
    fit_up = kbs.odr_fitfunction_linear(popt_up, x_fit)
    fit_dw= kbs.odr_fitfunction_linear(popt_dw, x_fit)

    plt.fill_between(x_fit, fit_up, fit_dw, alpha=0.25, label='3-sigma interval',color='red')
    plt.plot(x_fit, fit, 'r', lw=2, label='Best fit curve')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [2,100]
    else:
        axisrange = [0,70]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot([5,5],axisrange,'--k',lw=lthick)
    plt.plot(axisrange,[5,5],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT Sersic model', fontsize=Fsize)
    plt.ylabel('max(S/N) TDOSE 0.5'' aperture extraction', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2N_modelVSleadline.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.plot(SNmodelimg,SNleadline,'ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if (SNmodelimg[ii]/SNleadline[ii] > outliercut):
                plt.text(SNmodelimg[ii],SNleadline[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            if (SNleadline[ii]/SNmodelimg[ii] > outliercut):
                plt.text(SNmodelimg[ii],SNleadline[ii],id,fontsize=Fsize-5,horizontalalignment='right')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # overplot ODR fit
    popt = odrfit_SNmodelimg_SNleadline.beta
    perr = odrfit_SNmodelimg_SNleadline.sd_beta
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    x_fit = np.linspace(min(SNmodelimg), max(SNmodelimg), 100)
    fit = kbs.odr_fitfunction_linear(popt, x_fit)
    fit_up = kbs.odr_fitfunction_linear(popt_up, x_fit)
    fit_dw= kbs.odr_fitfunction_linear(popt_dw, x_fit)

    plt.fill_between(x_fit, fit_up, fit_dw, alpha=0.25, label='3-sigma interval',color='red')
    plt.plot(x_fit, fit, 'r', lw=2, label='Best fit curve')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [1,100]
    else:
        axisrange = [0,70]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot([5,5],axisrange,'--k',lw=lthick)
    plt.plot(axisrange,[5,5],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT Sersic model (intrinsic)', fontsize=Fsize)
    plt.ylabel('S/N MW DR1 LSDCat measurement', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_S2N_modelfiltVSleadline.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.plot(SNmodelimg_filtered,SNleadline,'ok',alpha=0.5    ,markersize=3.0)
    print("---> Number of LAE points plotted"+str(len(SNmodelimg)))
    if showids:
        for ii, id in enumerate(ids):
            if (SNmodelimg_filtered[ii]/SNleadline[ii] > outliercut):
                plt.text(SNmodelimg_filtered[ii],SNleadline[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            if (SNleadline[ii]/SNmodelimg_filtered[ii] > outliercut):
                plt.text(SNmodelimg_filtered[ii],SNleadline[ii],id,fontsize=Fsize-5,horizontalalignment='right')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # overplot ODR fit
    popt = odrfit_SNmodelimgfilt_SNleadline.beta
    perr = odrfit_SNmodelimgfilt_SNleadline.sd_beta
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    x_fit = np.linspace(min(SNmodelimg_filtered), max(SNmodelimg_filtered), 100)
    fit = kbs.odr_fitfunction_linear(popt, x_fit)
    fit_up = kbs.odr_fitfunction_linear(popt_up, x_fit)
    fit_dw= kbs.odr_fitfunction_linear(popt_dw, x_fit)

    plt.fill_between(x_fit, fit_up, fit_dw, alpha=0.25, label='3-sigma interval',color='red')
    plt.plot(x_fit, fit, 'r', lw=2, label='Best fit curve')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        xaxisrange = [1,100]
        yaxisrange = [4.,100]
    else:
        xaxisrange = [0,70]
        yaxisrange = xaxisrange
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    plt.plot([5,5],axisrange,'--k',lw=lthick)
    plt.plot(axisrange,[5,5],'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(S/N) GALFIT Sersic model (filtered)', fontsize=Fsize)
    plt.ylabel('max(S/N) MW DR1 LSDCat measurement', fontsize=Fsize)
    plt.ylim(yaxisrange)
    plt.xlim(xaxisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = namebase+'_flux_modelaper.pdf'
    if logaxes: plotname = plotname.replace('.pdf','_log.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.1, top=0.95)
    Fsize  = fontsize
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('TDOSE 1D spectra'),fontsize=Fsize)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # scalefactor = 2.0
    # if scalefactor != 1.0:
    #     plt.plot(Fsersic,Fmwdr1el/scalefactor,'.r',label='S/N MW DR1 EL spectrum scaled by '+str(scalefactor))

    plt.errorbar(Fmodelimg,Faper, yerr=Ferraper, xerr=Ferrmodelimg,fmt='ok',alpha=0.5)
    if showids:
        for ii, id in enumerate(ids):
            if (Fmodelimg[ii]/Faper[ii] > outliercut):
                plt.text(Fmodelimg[ii],Faper[ii],id,fontsize=Fsize-5,horizontalalignment='left')
            if (Faper[ii]/Fmodelimg[ii] > outliercut):
                plt.text(Fmodelimg[ii],Faper[ii],id,fontsize=Fsize-5,horizontalalignment='right')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # overplot ODR fit
    popt = odrfit_Fmodelimg_Faper.beta
    perr = odrfit_Fmodelimg_Faper.sd_beta
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    x_fit = np.linspace(min(Fmodelimg), max(Fmodelimg), 100)
    fit = kbs.odr_fitfunction_linear(popt, x_fit)
    fit_up = kbs.odr_fitfunction_linear(popt_up, x_fit)
    fit_dw= kbs.odr_fitfunction_linear(popt_dw, x_fit)

    plt.fill_between(x_fit, fit_up, fit_dw, alpha=0.25, label='3-sigma interval',color='red')
    plt.plot(x_fit, fit, 'r', lw=2, label='Best fit curve')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if logaxes:
        axisrange = [100,4000]
    else:
        axisrange = [0,1500]
    plt.plot(axisrange,axisrange,'--k',lw=lthick)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlabel('max(F/[cgs]) GALFIT Sersic model', fontsize=Fsize)
    plt.ylabel('max(F/[cgs]) TDOSE 0.5'' aperture extraction', fontsize=Fsize)
    plt.ylim(axisrange)
    plt.xlim(axisrange)
    if logaxes:
        plt.xscale('log')
        plt.yscale('log')
    leg = plt.legend(fancybox=True, loc='upper left',prop={'size':Fsize},ncol=1,numpoints=1)
    leg.get_frame().set_alpha(0.7)
    print(' - Saving plot to '+plotname)
    plt.savefig(plotname)
    fig.clf()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #////////////////////////////////////////////////////////////////////////////////////////////////


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v
#
# A few functions to smooth spectrum as done by LSDCat when filtering cubes to determine S/N of emission lines
# Taken from wavelength_smooth_lib.py (author C. Herenz) available at:
# https://bitbucket.org/Knusper2000/lsdcat/src/master/lib/wavelength_smooth_lib.py
# For details on create_filter_matrix_vel see Appendix A.2 in:
# https://ttt.astro.su.se/~ehere/pdf/Herenz_MASTER_THESIS.pdf
import math as m
import pylab as p
from scipy import signal
from scipy.sparse import csr_matrix
import getopt
import sys
#from . import line_em_funcs as lef # my own library with convenience functions
import multiprocessing
import warnings
warnings.filterwarnings("ignore",category=FutureWarning)

def create_filter_matrix_vel(velocity,lambda_start=4800,cdelt=1.3,lMax=3463):
    """
    filter_matrix = create_filter_matrix_vel(velocity,
                                             lambda_start=4800,
                                             cdelt=1.3,lMax=3463)
    ---

    Creates the filter matrix which will then be multiplied to
    a vector containing the spectrum.

    <velocity> = linewidth in velocity space (i.e. constant rest-frame)
                 of emission line (FWHM, in km/s) for which the filter
                 will be optimized
    <lambda_start> =  wavelength corresponding to 1st spectral pixel (in Angstrom)
                      (default: 4800 Angstrom)
    <cdelt> = wavelength interval corresponding to each pixel (in Angstrom)
              (stored in CD3_3 in QSim Datacube-Header), default 1.3 Angstrom
    <lMax> = No. of spectral elements (stored in NAXIS3 in Qsim Datacube-Header)
    """
    # same as create_filter_matrix but this time uses velocity as input -
    # instead of \delta_\lambda and \lambda_0
    # to understand the source see create_filter_matrix function
    speedoflight = 299792.458 # speed of light in km/s (with
                              # sufficient accuracy...)
    C = 4 #truncation length
    sqrt2pi_reci = 1./m.sqrt(2*m.pi)
    fwhm_factor = 2*m.sqrt(2*m.log(2))
    velocity=float(velocity); lambda_start=float(lambda_start);
    cdelt = float(cdelt); lMax = float(lMax)

    # this is the the first line which differs from  create_filter_matrix
    #  lambda / delta_lambda -> velocity / speedoflight
    tMax = abs(
        ((velocity / speedoflight) * ((lambda_start / cdelt)+lMax))\
            /fwhm_factor
              )

    M = C*m.sqrt(tMax)+1  # should be big enough...
    M = int(m.ceil(2*M))
    lambda_start = lambda_start - cdelt*(M/2)
    if M % 2 == 0:
        M = M + 1
    h_l = p.zeros((int(lMax)+M-1,M))
    for i in range(int(lMax)+M-1):
        t = abs(
            ((velocity / speedoflight) * ((lambda_start / cdelt) + i))\
                /fwhm_factor
               )
        h_l[i,:] = (sqrt2pi_reci/t)*signal.gaussian(M,t)
    filter_matrix = p.zeros((int(lMax)+M-1,int(lMax)))
    for i in range(int(lMax)+M-1):
        if i < M:
            filter_matrix[i,0:i+1] = h_l[i,M-i-1:M]
        elif i >= M and i < int(lMax):
            filter_matrix[i,i+1-M:i+1] = h_l[i,:]
        elif i >= int(lMax):
            filter_matrix[i,i+1-M:int(lMax)] = h_l[i,0:M-(i+1-int(lMax))]
    return filter_matrix

# the loop which does the actual filtering on the array of spectra
def filter_spectrum(filter_matrix,data):
    """
    filtered_data = filter_spectrum(filter_matrix,data):
    the loop which does the actual filtering on the array of spectra
    """
    length = data.shape[1]
    filtered_data = p.zeros((filter_matrix.shape[0],length))
    filter_matrix_sparse = csr_matrix(filter_matrix)
    for i in range(length):
        # filtering implemented as matrix multiplication see my
        # Masters-Thesis for more description
        # https://ttt.astro.su.se/~ehere/pdf/Herenz_MASTER_THESIS.pdf
        # - Appendix A.2
        filtered_data[:,i] = filter_matrix_sparse.dot(data[:,i])

    return filtered_data

#
# ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =