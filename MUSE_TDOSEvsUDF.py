# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Scripts, functions and routines to for the comparison between TDOSE and "Lyon" extractions on UDF 10
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import os
import commands
import sys
import glob
import MiGs
import pyfits
import datetime
import numpy as np
import collections
import matplotlib
import shutil
import fits2ascii as f2a
from astropy.coordinates import SkyCoord
import MUSEWideUtilities as mu
import MUSEWidePlots as mwp
import kbsutilities as kbs
import tdose_utilities as tu
from astropy import wcs
import MUSE_TDOSEvsUDF as mtu
import ciiiEmitterCandidates as cec
from astropy.io.votable import parse_single_table
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def buildANDgenerate(clobber=True):
    """
    Convenience wrapper to build and generate all the files needed for the TDOSE run

    --- EXAMPLE OF USE ---
    import MUSE_TDOSEvsUDF as mtu
    mtu.buildANDgenerate()

    """
    sourcecatdir = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_sourcecats/'
    mtu.gen_sourcecat(sourcecatdir,LAEinfofile,modelcoord=True)

    SETUPinfofile = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/tdose_setupfiles/MUSEWide_infofile_arche_PSFupdate_LAEs.txt'
    mtu.gen_TDOSEsetupfiles(SETUPinfofile,clobber=clobber)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def run_TDOSEextraction():
    """
    Command (to copy-paste into ipython) to run TDOSE on setup files generated with mtu.gen_TDOSEsetupfiles()

    --- EXAMPLE OF USE ---
    copy-paste into ipython

    """
    # mkdir tdose_models tdose_cutouts tdose_spectra
    # ipython
    import tdose, glob
    import numpy as np
    Nsessions = 2

    setupfiles = ['/Users/kschmidt/work/TDOSE/UDFextractionComparison/tdose_setupfiles/tdose_setup_udf_udf10_06698.txt']
    #setupfiles = glob.glob('/Users/kschmidt/work/TDOSE/UDFextractionComparison/tdose_setupfiles/tdose_setup_udf_udf10_*.txt')

    bundles, paralleldic = tdose.perform_extractions_in_parallel(setupfiles,Nsessions=Nsessions,clobber=True,performcutout=True,store1Dspectra=True,plot1Dspectra=True,generateFullFoVmodel=True,generateOverviewPlots=True,skipextractedobjects=False,logterminaloutput=True,verbosePE=True,verbosefull=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_UDFobj_sourcecats(sourcefiledir='/Users/kschmidt/work/TDOSE/UDFextractionComparison/udf10_c042_e031_withz_iter4/',
                          outputdir='/Users/kschmidt/work/TDOSE/UDFextractionComparison/tdose_sourcecats_individual/',
                          IMGEXT='IMA_HST_F814W_DATA',verbose=True):
    """
    Generating source catalogs for TDOSE based on UDF source files

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import MUSE_TDOSEvsUDF as mtu
    mtu.gen_UDFobj_sourcecats(verbose=True)

    """
    sourcefiles = glob.glob(sourcefiledir+'udf*.fits')
    Nsources    = len(sourcefiles)
    if verbose: print ' - Found '+str(Nsources)+' source files to generate TDOSE source catalogs for'

    for ss, sourcefile in enumerate(sourcefiles):
        source_mainhdr = pyfits.open(sourcefile)[0].header
        source_imghdr  = pyfits.open(sourcefile)[IMGEXT].header

        id       = source_mainhdr['ID']
        parentid = id
        ra       = source_mainhdr['RA']
        dec      = source_mainhdr['DEC']
        ximg     = int(source_imghdr['NAXIS1']/2.)+1
        yimg     = int(source_imghdr['NAXIS2']/2.)+1

        basename    = sourcefile.split('/')[-1].split('.fits')[0]
        sourcecat   = outputdir+'tdose_sourcecat_'+basename+'.txt'
        fout = open(sourcecat,'w')
        fout.write('# TDOSE Source catalog for '+sourcefile.split('/')[-1]+' generated with MUSE_TDOSEvsUDF.gen_UDFobj_sourcecats() \n')
        fout.write('# \n')
        fout.write('# parent_id id ra dec x_image y_image fluxscale \n')

        objstr = str(parentid) + ' '
        objstr = objstr + str(id) + ' '
        objstr = objstr + str(ra) + ' '
        objstr = objstr + str(dec) + ' '
        objstr = objstr + str(ximg) + ' '
        objstr = objstr + str(yimg) + ' '
        objstr = objstr + ' 1.0000 ' + ' \n'

        fout.write(objstr)
        fout.close()

        sourcecat_fits = f2a.ascii2fits(sourcecat,asciinames=True,skip_header=2,fitsformat='D',verbose=verbose)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def match_CUBEtoRafelski2015(cat_raf2015='/Users/kschmidt/work/catalogs/uvudf_rafelski_2015.fits',pixorigin=1,
                              fitscube='/Users/kschmidt/work/TDOSE/UDFextractionComparison/DATACUBE-UDF10-v0.32.fits',verbose=True):
    """
    match Rafelski et al. 2015 UDF catalog to a MUSE cube FoV

    --- INPUT ---


    --- EXAMPLE OF USE ---
    import MUSE_TDOSEvsUDF as mtu
    mtu.match_CUBEtoRafelski2015()

    """
    cubehdr   = pyfits.open(fitscube)['DATA'].header
    cubewcs   = wcs.WCS(tu.strip_header(cubehdr.copy()))

    raf15      = pyfits.open(cat_raf2015)[1].data
    ids_raf15  = raf15['ID']
    ras_raf15  = raf15['RA']
    decs_raf15 = raf15['DEC']

    goodids    = []
    goodras    = []
    gooddecs   = []

    pixcoords  = wcs.utils.skycoord_to_pixel(SkyCoord(ras_raf15, decs_raf15, unit="deg"),cubewcs,origin=pixorigin)

    sourcecat   = fitscube.replace('.fits','_Rafelski2015match.txt')
    fout = open(sourcecat,'w')
    fout.write('# Object from the Rafelski et al. (2015) catalog within the FoV of '+fitscube+'\n')
    fout.write('# Generated with MUSE_TDOSEvsUDF.match_CUBEtoRafelski2015() \n')
    fout.write('# objectnumber id ra dec x_cube y_cube fluxscale \n')

    for oo, id_raf15 in enumerate(ids_raf15):
        xpix      = pixcoords[0][oo]
        ypix      = pixcoords[1][oo]

        if (xpix >= 1) & (xpix <= cubehdr['NAXIS1']) & (ypix >= 1) & (ypix <= cubehdr['NAXIS2']):
            conclusionstring = '   goodobj: (x,y) = ('+str("%.5d" % xpix)+','+str("%.5d" % ypix)+')   '
            objstr = str(oo+1) + ' '
            objstr = objstr + str(id_raf15) + ' '
            objstr = objstr + str(ras_raf15[oo]) + ' '
            objstr = objstr + str(decs_raf15[oo]) + ' '
            objstr = objstr + str(xpix) + ' '
            objstr = objstr + str(ypix) + ' '
            objstr = objstr + ' 1.0000 ' + ' \n'
            fout.write(objstr)
        else:
            conclusionstring = '   badobj:  (x,y) = ('+str("%.5d" % xpix)+','+str("%.5d" % ypix)+')   '

        if verbose:
            infostr = '   Checked '+str("%.6d" % id_raf15)+' ('+str("%.6d" % oo)+' / '+str("%.6d" % len(ids_raf15))+')  '+conclusionstring
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

    fout.close()
    if verbose: print '\n - Found '+str(len(goodids))+' objects in the Rafelski et al. (2015) ' \
                                                      'catalog wihtin the cube FoV \n   (cube = '+fitscube+')'
    if verbose: print '   - Stored the info in '+sourcecat

    sourcecat_fits = f2a.ascii2fits(sourcecat,asciinames=True,skip_header=2,fitsformat='D',verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_TDOSEsetupfiles(infofile,namebase='MUSEWide_tdose_setup_LAEs',clobber=False,
                        outputdir='/Users/kschmidt/work/TDOSE/UDFextractionComparison/tdose_setupfiles/',verbose=True):
    """
    Generate TDOSE setupfiles for the LAE extractions

    --- INPUT ---


    --- EXAMPLE OF USE ---
    import MUSE_TDOSEvsUDF as mtu
    infofile = '/Users/kschmidt/work/TDOSE/UDFextractionComparison/tdose_setupfiles/UDF10_TDOSEsetupbuild_infofile.txt'
    mtu.gen_TDOSEsetupfiles(infofile)

    """
    tu.duplicate_setup_template(outputdir,infofile,namebase=namebase,clobber=clobber,loopcols='all',infofmt="S250",infohdr=2)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def load_sourcefile_spectra(sourcefiles,spectypes=['SPE_MUSE_TOT_SKYSUB','SPE_MUSE_PSF_SKYSUB',
                                                   'SPE_MUSE_WHITE_SKYSUB','SPE_MUSE_APER_0.4_SKYSUB'],
                            verbose=True):
    """
    Loading source file an pulling out spectra

    --- INPUT ---
    import MUSE_TDOSEvsUDF as mtu

    sourcefiles = ['/Volumes/DATABCKUP1/UDFvsMUSEWide/udf10_c042_e031_withz_iter6/udf_udf10_00006.fits','/Volumes/DATABCKUP1/UDFvsMUSEWide/udf10_c042_e031_withz_iter6/udf_udf10_00533.fits']
    redshifts, specdic = mtu.load_sourcefile_spectra(sourcefiles)

    """
    Nfiles = len(sourcefiles)
    if verbose: print(' - Loading the spectra from the '+str(Nfiles)+' source files provided and returning dictionary')

    specdic_all = collections.OrderedDict()
    redshifts   = collections.OrderedDict()

    for sourcefile in sourcefiles:
        sfhdu   = pyfits.open(sourcefile)
        specdic = collections.OrderedDict()
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for spectype in spectypes:
            flux                 = sfhdu[spectype+'_DATA'].data
            error                = sfhdu[spectype+'_STAT'].data
            wave                 = (sfhdu[spectype+'_DATA'].header['CRVAL1'] -
                                    sfhdu[spectype+'_DATA'].header['CDELT1'] * (sfhdu[spectype+'_DATA'].header['CRPIX1'] - 1)) +\
                                   np.arange(len(flux)) * sfhdu[spectype+'_DATA'].header['CDELT1']
            specarr              = np.zeros(len(flux),dtype=[('lambda',np.float32),('flux',np.float32),('fluxerror',np.float32)])
            specarr['lambda']    = wave
            specarr['flux']      = flux
            specarr['fluxerror'] = np.sqrt(error)

            specdic[spectype]    = specarr
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        specdic_all[sourcefile] = specdic
        redshifts[sourcefile]   = sfhdu['Z'].data
    return redshifts, specdic_all

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_UDFspec_fits(sourcefiles,outputdir,clobber=False,verbose=True):
    """

    Load the 1D spectra from a set of source files, and save the indivdual spectra as
    stand alone fits files (makes comparison plotting with MUSEWidePlots.plot_1DspecOverview() straight forward)

    --- EXAMPLE OF USE ---
    import MUSE_TDOSEvsUDF as mtu

    sourcefiles = ['/Volumes/DATABCKUP1/UDFvsMUSEWide/udf10_c042_e031_withz_iter6/udf_udf10_00004.fits','/Volumes/DATABCKUP1/UDFvsMUSEWide/udf10_c042_e031_withz_iter6/udf_udf10_00006.fits','/Volumes/DATABCKUP1/UDFvsMUSEWide/udf10_c042_e031_withz_iter6/udf_udf10_00533.fits']
    outputdir  = '/Volumes/DATABCKUP1/UDFvsMUSEWide/udf10_c042_e031_withz_iter6_1Dspecs/'

    mtu.gen_UDFspec_fits(sourcefiles,outputdir,clobber=False)

    """
    if not os.path.isdir(outputdir):
        sys.exit('Output directory specified ('+outputdir+') does not appear to exist')

    redshifts, specdic = mtu.load_sourcefile_spectra(sourcefiles,verbose=verbose)

    if verbose: print(' - Storing spectra extracted from source files to individual fits files')

    for sourcefile in sourcefiles:
        for spec in specdic[sourcefile].keys():
            if spec.startswith('SPE_'):
                specarr  = specdic[sourcefile][spec]
                fitsname = outputdir+sourcefile.split('/')[-1].replace('.fits','_'+spec+'.fits')
                if verbose: print('   Generating '+fitsname)
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                c1  = pyfits.Column(name='lambda',    format='D', unit='ang',                 array=specarr['lambda'])
                c2  = pyfits.Column(name='flux',      format='D', unit='1e-20 erg/ang/cm2/s', array=specarr['flux'])
                c3  = pyfits.Column(name='fluxerror', format='D', unit='1e-20 erg/ang/cm2/s', array=specarr['fluxerror'])

                coldefs = pyfits.ColDefs([c1,c2,c3])
                tb      = pyfits.new_table(coldefs) # creating default header
                head    = tb.header
                tbHDU   = pyfits.new_table(coldefs, header=head)
                tbHDU.writeto(fitsname, clobber=clobber)
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def load_UDFmastertable(votablename):
    """

    Loading the master ID list VOTable
    NB that there is now also a fits version of this table

    --- EXAMPLE OF USE ---
    import MUSE_TDOSEvsUDF as mtu
    votablename = '/Volumes/DATABCKUP1/UDFvsMUSEWide/master_idlist_20170214.vot'
    table = mtu.load_UDFmastertable(votablename)

    """
    table = parse_single_table(votablename)

    return table.array
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1Doverviewcomparison_crossobject_wrapper(plotdir,specdir_tdose,specdir_udf,idset='nonLAEs',
                                                  specs2show=['TDOSE gauss','PSF SKYSUB'], # ,'TDOSE aperture'
                                                  yrangefull=[-200,1000], xrangefull=[4500,9500],clobber=False):
    """

    """

    #### getting id redshifts
    if idset == 'LAEs':
        iddat   = np.genfromtxt('/Users/kschmidt/work/MUSE/UDFvsMUSEWide/blends_UDF10.txt',names=True,comments='#')
    elif idset == 'nonLAEs':
        iddat   = np.genfromtxt('/Users/kschmidt/work/MUSE/UDFvsMUSEWide/blends_UDF10_pt2.txt',names=True,comments='#')
    else:
        sys.exit('Unknown idset in plot_1Doverviewcomparison_crossobject_wrapper()')

    MUSEids = iddat['id_MUSE']
    zdat    = pyfits.open('/Users/kschmidt/work/MUSE/UDFvsMUSEWide/master_idlist_20170214.fits')[1].data
    for id in MUSEids:
        objent = np.where(zdat['ID'] == id)[0]
        print str(int(id))+'       '+str(zdat['RAF_ID'][objent[0]])+'       # z = '+str(zdat['UDF10_Z'][objent[0]])

    idsetdic = collections.OrderedDict()

    if idset == 'LAEs':
        idsetdic['01'] = [[50,208,183], 5.7810, 1, None]
        idsetdic['02'] = [[50,208,183], 3.3178, 1, None]
        idsetdic['03'] = [[50,208,183], 3.1863, 1, None]
        idsetdic['04'] = [[559,218,44], 4.5120, 1, None]
        idsetdic['05'] = [[559,218,44], 3.0483, 1, None]
        idsetdic['06'] = [[559,218,44], 1.6103, 1, None]
        idsetdic['07'] = [[587,720,6336], 5.9765, 2, [1,0,0]]
        idsetdic['08'] = [[587,720,6336], 3.71472, 2, [1,0,0]]
        idsetdic['09'] = [[6300,277], 4.5221, 1, None]
        idsetdic['10'] = [[6300,277], 3.7148, 1, None]
        idsetdic['11'] = [[6312,6769,6680], 6.3097, 1, None]
        idsetdic['12'] = [[6312,6769,6680], 4.5046, 1, None]
        idsetdic['13'] = [[6311], 5.4200, 1, None]
        idsetdic['14'] = [[109,748,6327], 3.0863, 1, None]
        idsetdic['15'] = [[109,748,6327], 5.1351, 1, None]
        idsetdic['16'] = [[287], 3.7242, 1, None]
        idsetdic['17'] = [[319], 6.5658, 1, None]
        idsetdic['18'] = [[304,173,46], 5.9405, 1, None]
        idsetdic['19'] = [[304,173,46], 1.4135, 1, None]
        idsetdic['20'] = [[459,823], 6.1852, 1, None]
        idsetdic['21'] = [[6291], 3.2994, 1, None]
    elif idset == 'nonLAEs':
        idsetdic['01'] = [[4,11,12,14], 0.7650, 1, None]
        idsetdic['02'] = [[4,11,12,14], 1.0376, 1, None]
        idsetdic['03'] = [[4,11,12,14], 0.9967, 1, None]
        idsetdic['04'] = [[21,79],  0.8939, 1, None]
        idsetdic['05'] = [[21,79],  0.3631, 1, None]
        idsetdic['06'] = [[26,40],  1.1348, 1, None]
        idsetdic['07'] = [[26,40],  0.9798, 1, None]
        idsetdic['08'] = [[32,121], 1.3066, 1, None]
        idsetdic['09'] = [[32,121], 1.3060, 1, None]
        idsetdic['10'] = [[38,627], 0.7655, 1, None]
        idsetdic['11'] = [[38,627], 5.1352, 1, None]
        idsetdic['12'] = [[46,92],  1.4135, 1, None]
        idsetdic['13'] = [[46,92],  1.4137, 1, None]
        idsetdic['14'] = [[65,605], 1.3069, 1, None]
        idsetdic['15'] = [[65,605], 3.4390, 1, None]
        idsetdic['16'] = [[69,72],  0.5439, 1, None]
        idsetdic['17'] = [[69,72],  1.0969, 1, None]
        idsetdic['18'] = [[93,156], 0.6181, 1, None]
        idsetdic['19'] = [[93,156], 0.8588, 1, None]
        idsetdic['20'] = [[96,114], 0.6220, 2, [0,1]]
        idsetdic['21'] = [[96,114], 0.6777, 2, [0,1]]

        idsetdic['22'] = [[2,6368], 0.4193, 1, None]
        idsetdic['23'] = [[2,6368], 0.2754, 1, None]
        idsetdic['24'] = [[33,6324], 1.4153, 1, None]
        idsetdic['25'] = [[33,6324], 5.1321, 1, None]
        idsetdic['26'] = [[48,6322], 0.4681, 1, None]
        idsetdic['27'] = [[48,6322], 3.6056, 1, None]


    for key in idsetdic.keys():
        idset        = idsetdic[key][0]
        plotredshift = idsetdic[key][1]
        Nrafids      = idsetdic[key][2]


        if Nrafids == 2:
            namestring   = 'ids'+'-'.join(np.asarray(idset).astype(str))+'_rafid2'
            rafent       = idsetdic[key][3]
            mtu.plot_1Doverviewcomparison_crossobject(idset, plotredshift, plotdir, specdir_tdose=specdir_tdose, rafident=rafent,
                                                      namestring=namestring, specdir_udf=specdir_udf, clobber=clobber,
                                                      yrangefull=yrangefull, xrangefull=xrangefull, specs2show=specs2show)

        namestring   = 'ids'+'-'.join(np.asarray(idset).astype(str))
        mtu.plot_1Doverviewcomparison_crossobject(idset, plotredshift, plotdir, specdir_tdose=specdir_tdose,
                                                  rafident=0,
                                                  namestring=namestring, specdir_udf=specdir_udf, clobber=clobber,
                                                  yrangefull=yrangefull, xrangefull=xrangefull, specs2show=specs2show)



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1Doverviewcomparison_crossobject(ids,plotredshift,outputdir,namestring='_RENAME_',rafident=0,
                                          specs2show=['TDOSE gauss','TDOSE aperture','PSF SKYSUB'],
                                          yrangefull=[-1000,2000],xrangefull=[4500,18000],
                                          specdir_udf='/Volumes/DATABCKUP1/UDFvsMUSEWide/udf10_c042_e031_withz_iter6_1Dspecs/',
                                          specdir_tdose='/Volumes/DATABCKUP1/UDFvsMUSEWide/tdose_spectra/',
                                          clobber=False,verbose=True):
    """
    Plot 1D overview plots of UDF, MUSE-Wide and 3D-HST spectra

    --- INPUT ---
    ids      MUSE ids to plot

    --- EXAMPLE OF USE ---
    import MUSE_TDOSEvsUDF as mtu

    mtu.plot_1Doverviewcomparison([4,6,533],1.2,'/Volumes/DATABCKUP1/UDFvsMUSEWide/overviewplots/',clobber=False)

    """
    try:
        float(rafident)
        rafident = [rafident]*len(ids)
    except:
        pass

    votablename = '/Volumes/DATABCKUP1/UDFvsMUSEWide/master_idlist_20170214.vot'
    table = mtu.load_UDFmastertable(votablename)

    colormaps    = ['Reds','Greens','Blues','Oranges','Purples','Greys']
    speccolmaps  = collections.OrderedDict()
    for cmapuse in colormaps:
        cmap = matplotlib.cm.get_cmap(cmapuse)
        norm = matplotlib.colors.Normalize(vmin=0, vmax=len(specs2show)*2.5)
        speccolmaps[cmapuse]    = [ cmap(norm(objnumner+len(specs2show))) for objnumner in np.arange(len(specs2show))]

    # speccols['red']    = ['darkred','red','tomato','pink','lightred']
    # speccols['green']  = ['darkgreen','green','seagreen','lime','lightgreen']
    # speccols['blue']   = ['darkblue','blue','steelblue','cyan','lightblue']
    # speccols['orange'] = ['darkorange','orange','gold','yellow','khaki']

    Nids = len(ids)
    if verbose: print(' - Will attempt to generate plots for the spectra found for the '+str(Nids)+' IDs in list')
    spectra, labels, wavecols, fluxcols, ferrcols, skyspecs, wavecols_sky, fluxcols_sky = [],[],[],[],[],[],[],[]
    speccols_all = []
    for ii, id in enumerate(ids):
        speccols = speccolmaps[speccolmaps.keys()[ii]]
        objent_udf = np.where(table['ID'] == int(id))[0]
        Nent_udf   = len(objent_udf)
        if Nent_udf == 0:
            print(' WARNING No match to id='+str(id)+' found in master_idlist_20170214.vot')
        elif Nent_udf > 1:
            print(' WARNING More than one match to id='+str(id)+' found in master_idlist_20170214.vot; using the first match')
            objent_udf = objent_udf[0]

        rafID            = table['RAF_ID'][objent_udf]
        obj_UDFfile      = table['UDF10_FILENAME'][objent_udf]

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #                             UDF-10
        spec_udf         = glob.glob(specdir_udf+obj_UDFfile.data[0].replace('.fits','*.fits'))

        spec_udf_show    = []
        for su in spec_udf:
            lab_udf          = ' '.join(su.split('udf_udf10_')[-1].split('.fit')[0][6:].split('_')[2:])
            if lab_udf in specs2show:
                spec_udf_show.append(su)
        Nspecudf         = len(spec_udf_show)

        wave_udf         = ['lambda']*Nspecudf
        flux_udf         = ['flux']*Nspecudf
        ferr_udf         = ['fluxerror']*Nspecudf
        sky_udf          = ['/Users/kschmidt/work/MUSE/skyspectra/SKY_SPECTRUM_candels-cdfs-35_av.fits']+[None]*(Nspecudf-1)
        sky_wc_udf       = ['lambda']+[None]*(Nspecudf-1)
        sky_fc_udf       = ['data']+[None]*(Nspecudf-1)

        lab_udf          = [' '.join(ss.split('udf_udf10_')[-1].split('.fit')[0][6:].split('_')[2:])+' id'+str(int(id))
                            for ss in spec_udf_show]
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #                             TDOSE
        spec_tdose       = glob.glob(specdir_tdose+'tdose_spectrum_*'+
                                     str("%.10d" % int(rafID.data[0].split(',')[rafident[ii]]))+'*.fits')

        spec_tdose_show  = []
        for st in spec_tdose:
            lab_tdose          = 'TDOSE '+st.split('tdose_spectrum_UDF-10_')[-1].split('_')[0]
            if lab_tdose in specs2show:
                spec_tdose_show.append(st)
        Nspectdose       = len(spec_tdose_show)

        wave_tdose       = ['wave']*Nspectdose
        flux_tdose       = ['flux']*Nspectdose
        ferr_tdose       = ['fluxerror']*Nspectdose
        sky_tdose        = [None]*Nspectdose
        sky_wc_tdose     = [None]*Nspectdose
        sky_fc_tdose     = [None]*Nspectdose

        lab_tdose        = ['TDOSE '+ss.split('tdose_spectrum_UDF-10_')[-1].split('_')[0]+' id'+str(int(id))
                            for ss in spec_tdose_show]


        spectra      = spectra + spec_udf_show + spec_tdose_show
        speccols_all = speccols_all + speccols
        if len(spectra) == 0:
            if verbose: print(' WARNING No spectra found for object id='+str(id)+' with source file '+obj_UDFfile.data[0])
        else:
            labels       = labels + lab_udf  + lab_tdose

            wavecols     = wavecols + wave_udf + wave_tdose
            fluxcols     = fluxcols + flux_udf + flux_tdose
            ferrcols     = ferrcols + ferr_udf + ferr_tdose

            skyspecs     = skyspecs + sky_udf  + sky_tdose
            wavecols_sky = wavecols_sky + sky_wc_udf + sky_wc_tdose
            fluxcols_sky = fluxcols_sky + sky_fc_udf + sky_fc_tdose

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #                            PLOTTING
    voffset      = 0.0
    outputfig    = outputdir+'UDF10_1Doverviewplots_crossobj_'+namestring+'.pdf'
    for plotSN in [True,False]:
        serachstr   = outputfig.replace('plots.pdf','plots*.pdf')
        filesstored = glob.glob(serachstr)
        if (len(filesstored) > 0) & (clobber==False):
            if verbose: print(' WARNING Skipping, as clobber=False and plot exists for:\n   '+serachstr)
        else:
            mwp.plot_1DspecOverview(spectra, labels, wavecols, fluxcols, ferrcols, plotredshift, voffset=voffset,
                                    skyspectra=skyspecs,wavecols_sky=wavecols_sky, fluxcols_sky=fluxcols_sky,
                                    outputfigure=outputfig,yrangefull=yrangefull, xrangefull=xrangefull,
                                    plotSN=plotSN,speccols=speccols_all,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_1Doverviewcomparison(ids,outputdir,yrangefull=[-1000,2000],xrangefull=[4500,18000],rafident=0,
                              specdir_udf='/Volumes/DATABCKUP1/UDFvsMUSEWide/udf10_c042_e031_withz_iter6_1Dspecs/',
                              specdir_tdose='/Volumes/DATABCKUP1/UDFvsMUSEWide/tdose_spectra/',
                              clobber=False,verbose=True):
    """
    Plot 1D overview plots of UDF, MUSE-Wide and 3D-HST spectra

    --- EXAMPLE OF USE ---
    import MUSE_TDOSEvsUDF as mtu

    mtu.plot_1Doverviewcomparison([4,6,533],'/Volumes/DATABCKUP1/UDFvsMUSEWide/overviewplots/',clobber=False)

    """
    votablename = '/Volumes/DATABCKUP1/UDFvsMUSEWide/master_idlist_20170214.vot'
    table = mtu.load_UDFmastertable(votablename)

    if ids == 'all':
        ids = table['ID']

    Nids = len(ids)
    if verbose: print(' - Will attempt to generate plots for the spectra found for the '+str(Nids)+' IDs in list')
    for id in ids:
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #                             UDF-10
        objent_udf = np.where(table['ID'] == int(id))[0]
        Nent_udf   = len(objent_udf)
        if Nent_udf == 0:
            print(' WARNING No match to id='+str(id)+' found in master_idlist_20170214.vot')
        elif Nent_udf > 1:
            print(' WARNING More than one match to id='+str(id)+' found in master_idlist_20170214.vot; using the first match')
            objent_udf = objent_udf[0]

        objRA            = table['DEC'][objent_udf]
        objDEC           = table['RA'][objent_udf]
        rafID            = table['RAF_ID'][objent_udf]
        obj_UDFfile      = table['UDF10_FILENAME'][objent_udf]
        obj_Mosaicfil    = table['MOSAIC_FILENAME'][objent_udf]

        sourcefiledir    = '/Volumes/DATABCKUP1/UDFvsMUSEWide/udf10_c042_e031_withz_iter6/'
        sourcefile       = sourcefiledir+obj_UDFfile.data[0]

        zs, specdic      = mtu.load_sourcefile_spectra([sourcefile])
        obj_z_ent_muse   = np.where(zs[sourcefile]['Z_DESC'] == 'MUSE')[0]
        if len(obj_z_ent_muse) == 1:
            obj_z = zs[sourcefile]['Z'][obj_z_ent_muse][0]
        else:
            obj_z_ent_eazy   = np.where(zs[sourcefile]['Z_DESC'] == 'EAZY')[0]
            obj_z            = zs[sourcefile]['Z'][obj_z_ent_eazy][0]

        spec_udf         = glob.glob(specdir_udf+obj_UDFfile.data[0].replace('.fits','*.fits'))

        wave_udf         = ['lambda']*len(spec_udf)
        flux_udf         = ['flux']*len(spec_udf)
        ferr_udf         = ['fluxerror']*len(spec_udf)
        sky_udf          = ['/Users/kschmidt/work/MUSE/skyspectra/SKY_SPECTRUM_candels-cdfs-35_av.fits']+[None]*(len(spec_udf)-1)
        sky_wc_udf       = ['lambda']+[None]*len(spec_udf)
        sky_fc_udf       = ['data']+[None]*len(spec_udf)

        lab_udf          = [' '.join(ss.split('udf_udf10_')[-1].split('.fit')[0][6:].split('_')[2:]) for ss in spec_udf]
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #                             TDOSE
        spec_tdose       = glob.glob(specdir_tdose+'tdose_spectrum_*'+
                                     str("%.10d" % int(rafID.data[0].split(',')[rafident]))+'*.fits')
        wave_tdose       = ['wave']*len(spec_tdose)
        flux_tdose       = ['flux']*len(spec_tdose)
        ferr_tdose       = ['fluxerror']*len(spec_tdose)
        sky_tdose        = [None]*len(spec_tdose)
        sky_wc_tdose     = [None]*len(spec_tdose)
        sky_fc_tdose     = [None]*len(spec_tdose)

        lab_tdose        = ['TDOSE '+ss.split('tdose_spectrum_UDF-10_')[-1].split('_')[0] for ss in spec_tdose]
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #                             3D-HST
        spec_3dhst     = []
        if len(spec_3dhst) > 0:
            sky_3dhst      = ['/Users/kschmidt/work/MUSE/skytable.fits']+[None]*(len(spec_3dhst)-1)
        else:
            sky_3dhst      = []
        sky_wc_3dhst   = ['lam']+[None]*len(spec_udf)
        sky_fc_3dhst   = ['flux']+[None]*len(spec_udf)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #                            PLOTTING
        spectra      = spec_udf + spec_tdose + spec_3dhst
        if len(spectra) == 0:
            if verbose: print(' WARNING No spectra found for object id='+str(id)+' with source file '+obj_UDFfile.data[0])
        else:
            labels       = lab_udf  + lab_tdose  #+ lab_3dhst

            wavecols     = wave_udf + wave_tdose #+ wave_3dhst
            fluxcols     = flux_udf + flux_tdose #+ flux_3dhst
            ferrcols     = ferr_udf + ferr_tdose #+ ferr_3dhst

            skyspecs     = sky_udf  + sky_tdose  + sky_3dhst
            wavecols_sky = sky_wc_udf + sky_wc_tdose + sky_wc_3dhst
            fluxcols_sky = sky_fc_udf + sky_fc_tdose + sky_fc_3dhst

            redshift     = obj_z
            voffset      = 0.0

            outputfig    = outputdir+obj_UDFfile.data[0].replace('.fits','_1Doverviewplots.pdf')
            for plotSN in [True,False]:
                serachstr   = outputfig.replace('plots.pdf','plots*.pdf')
                filesstored = glob.glob(serachstr)
                if (len(filesstored) > 0) & (clobber==False):
                    if verbose: print(' WARNING Skipping, as clobber=False and plot exists for:\n   '+serachstr)
                else:
                    mwp.plot_1DspecOverview(spectra, labels, wavecols, fluxcols, ferrcols, redshift, voffset=voffset,
                                            skyspectra=skyspecs,wavecols_sky=wavecols_sky, fluxcols_sky=fluxcols_sky,
                                            outputfigure=outputfig,yrangefull=yrangefull, xrangefull=xrangefull,
                                            plotSN=plotSN,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_sourcecat_FromFitsCats(outputdir,centralcoords=[[53.160554,-27.778928]],catnames=['UDF-10'],addMUSEcat=True,
                               fluxcol='FLUX_F775W',fluxfactor=1.0,returnSeparations=False,clobber=False,verbose=True):
    """
    Generate source catalogs from UDF 10 fits catalog and Rafelski catalog

    --- INPUT ---
    outputdir            Output directory to contain object source catalogs
    centralcoords        list of coordinates [ra,deg] in degrees to generate catalogs around
    catnames             Names of the catalogs generated around centralcoords
    addMUSEcat           Adding the MUSE catalog (including ORION-only detections)
    fluxcol              Flux column to use as for flux scaling of sources
    fluxfactor           Value to use as flux factor in source catalog.
                         If returnSeparations=True this is overwritten.
    returnSeparations    Return the separations in the fluxscale column of the source catalog
    clobber              Overwrite existing files
    Verbose              Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSE_TDOSEvsUDF as mtu

    outdir = '/Users/kschmidt/work/MUSE/UDFvsMUSEWide/tdose_sourcecats/'
    udf10coord = [53.160554,-27.778928]

    mtu.gen_sourcecat_FromFitsCats(outdir,centralcoords=[udf10coord],catnames=['UDF-10'],addMUSEcat=False,clobber=False)

    """
    scat_radius  = 70.0 # using a radius of 70 arc sec to ensure all relevant objects are included in source catalog
    cat_raf      = '/Users/kschmidt/work/catalogs/uvudf_rafelski_2015.fits'

    cat_udf      = '/Users/kschmidt/work/MUSE/UDFvsMUSEWide/udf10_c042_e031_withz_iter6.fits'
    dat_udf      = pyfits.open(cat_udf)[1].data

    refimg       = '/Users/kschmidt/work/images_MAST/hlsp_xdf_hst_acswfc-30mas_hudf_f814w_v1_sci.fits'
    imgheader    = pyfits.open(refimg)[0].header

    for cc, ccoord in enumerate(centralcoords):
        ra_cen, dec_cen = ccoord
        outname      = outputdir+'tdose_sourcecat_from_fitscat_'+catnames[cc]+'.txt'
        sourcelist = []

        if addMUSEcat:
            for ss, udfid in enumerate(dat_udf['ID']):
                sourcelist.append([int(udfid+1e9), int(udfid+1e9), dat_udf['RA'][ss], dat_udf['DEC'][ss], 1])

        if type(fluxfactor) == str:
            fluxfactor = dat_udf

        if returnSeparations:
            fluxfactor = 'separation'

        sourcecat    = tu.gen_sourcecat_from_FitsCat(cat_raf,'ID','RA','DEC',[ra_cen,dec_cen],scat_radius,imgheader,
                                                     fluxfactor=fluxfactor,fluxcol=fluxcol,outname=outname,newsources=sourcelist,
                                                     clobber=clobber,verbose=verbose)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =