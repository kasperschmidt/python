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
import shutil
import fits2ascii as f2a
from astropy.coordinates import SkyCoord
import MUSEWideUtilities as mu
import kbsutilities as kbs
import tdose_utilities as tu
from astropy import wcs
import MUSE_TDOSEvsUDF as mtu
import ciiiEmitterCandidates as cec
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
    specdic = mtu.load_sourcefile_spectra(sourcefiles)

    """
    Nfiles = len(sourcefiles)
    if verbose: print(' - Loading the spectra from the '+str(Nfiles)+' source files provided and returning dictionary')

    specdic_all = collections.OrderedDict()

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
            specarr['fluxerror'] = error

            specdic[spectype]    = specarr
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        specdic_all[sourcefile] = specdic
    return specdic_all

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_UDFspec_fits(sourcefiles,outputdir,clobber=False,verbose=True):
    """

    Load the 1D spectra from a set of source files, and save the indivdual spectra as
    stand alone fits files (makes comparison plotting with MUSEWidePlots.plot_1DspecOverview() straight forward)

    --- EXAMPLE OF USE ---
    import MUSE_TDOSEvsUDF as mtu

    sourcefiles = ['/Volumes/DATABCKUP1/UDFvsMUSEWide/udf10_c042_e031_withz_iter6/udf_udf10_00006.fits','/Volumes/DATABCKUP1/UDFvsMUSEWide/udf10_c042_e031_withz_iter6/udf_udf10_00533.fits']
    outputdir  = '/Volumes/DATABCKUP1/UDFvsMUSEWide/udf10_c042_e031_withz_iter6_1Dspecs/'

    mtu.gen_UDFspec_fits(sourcefiles,outputdir,clobber=False)

    """
    if not os.path.isdir(outputdir):
        sys.exit('Output directory specified ('+outputdir+') does not appear to exist')

    specdic = mtu.load_sourcefile_spectra(sourcefiles,verbose=verbose)

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