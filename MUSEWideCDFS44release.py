# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                     Utilities for preparing MUSE-Wide CDFS first 44 fields release
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pyfits
import numpy as np
import MUSEWideUtilities as mwu
import MUSEWideCDFS44release as mw44
import tdose_utilities as tu
import commands
import glob
import tdose
import pdb
import os
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def run_arche_notes():
    """
    The commands described in the MUSE notes 170721 section 22.3.7

    """
    #First I created the setup catalogs (on my laptop) with:

    outputdir = '/Users/kschmidt/work/TDOSE/muse_tdose_setups/'
    infofile  = outputdir+'musewide_infofile_arche_PSFupdate.txt'
    tu.duplicate_setup_template(outputdir,infofile,namebase='MUSEWide_tdose_setup_arche',loopcols='all',clobber=True)

    #Here the PSFupdate comes from running
    mwu.insert_PSFfits_into_TDOSEinfofile(duplicatefirstline=True)

    # convert the SExtractor catalogs into fits catalogs to use them for initial guesses of the Gaussian modeling:
    catalogs = glob.glob('/Volumes/DATABCKUP2/MUSE-Wide/catalogs_photometry/catalog_photometry_candels-cdfs-*.cat')
    tu.SExtractorCat2fits(catalogs,stringcols=[1],header=73,verbose=True)

    # generate the source catalogs for the CDFS fields (as already specified in the setup files) with:
    # catfile = '/store/data/musewide/candels-cdfs-01/catalog_photometry_candels_cdfs_01.cat'

    catfiles = glob.glob('/store/data/musewide/candels*/catalog_photometry_candels*.cat')
    for catfile in catfiles:
        field = catfile.split('/')[-2]
        output = '/store/data/musewide/TDOSE/tdose_sourcecats/'+catfile.split('/')[-1].split('.')[0]+'_tdose_sourcecat.txt'
        output = output.replace('candels-cdfs-','candels_cdfs_')
        imgfile = '/store/data/musewide/'+field+'/acs_814w_'+field+'_cut_v2.0.fits'

        if os.path.isfile(imgfile):
            print ' - Creating source catalog for '+field+' using: \n   '+catfile+' \n   '+imgfile
            imgheader = pyfits.open(imgfile)[0].header
            sourcecat     = tu.gen_sourcecat_from_SExtractorfile(catfile,output,imgheader=imgheader,idcol=0,racol=2,deccol=3,fluxcol=22,fluxfactor=100.,generateDS9reg=True,verbose=False,clobber=True)

        else:
            print ' - WARNING Imgfile not found ('+imgfile+'); skipping generating catalog for '+field

    # And from there simply run TDOSE with:
    setupfile = '/store/data/musewide/TDOSE/tdose_setupfiles/MUSEWide_tdose_setup_arche_candels-cdfs-01.txt'

    tdose.perform_extraction(setupfile=setupfile, verbose=True, verbosefull=True, clobber=True, performcutout=True, logterminaloutput=True)

    tdose.gen_fullFoV_from_cutouts(setupfile)

    tu.gen_overview_plot('all',setupfile,skipobj=True)

    #Or in parallel for multiple setup files with:
    Nsessions = 4

    setupfiles = ['/store/data/musewide/TDOSE/tdose_setupfiles/MUSEWide_tdose_setup_arche_candels-cdfs-03.txt', '/store/data/musewide/TDOSE/tdose_setupfiles/MUSEWide_tdose_setup_arche_candels-cdfs-04.txt', '/store/data/musewide/TDOSE/tdose_setupfiles/MUSEWide_tdose_setup_arche_candels-cdfs-05.txt', '/store/data/musewide/TDOSE/tdose_setupfiles/MUSEWide_tdose_setup_arche_candels-cdfs-06.txt']

    bundles, paralleldic = tdose.perform_extractions_in_parallel(setupfiles,Nsessions=Nsessions,clobber=True,performcutout=True,store1Dspectra=True,plot1Dspectra=True,generateFullFoVmodel=True,generateOverviewPlots=True,skipextractedobjects=True,logterminaloutput=True)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_repeat_files(password,scpdir='/Users/kschmidt/work/MUSE/MUSEWide_CFDS44release/repeat_overviewplots_170714/',
                     searchstr='tdose_spectra/*FFFFF*IIIII*source_overview.pdf',fieldnameunderscore=False,verbose=True):
    """
    scp files for the ids in the list of object to repeat.
    Default is to copy over the overview plots from the 170714 extraction to laptop

    --- INPUT ---
    password                password for login to arche
    scpdir                  directory to copy found files to
    searchstr               string to use for locating files.
                                FFFFF is replaced with the fieldname
                                IIIII is replaced with the object id
    fieldnameunderscore     To replace '-' with '_' in field name set to True
    verbose                 toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWideCDFS44release as mw44

    # - - - - Overviewplots - - - -
    scpdir    = '/Users/kschmidt/work/MUSE/MUSEWide_CFDS44release/repeat_overviewplots_170714/'
    searchstr = 'tdose_spectra/*FFFFF*IIIII*source_overview.pdf'
    mw44.get_repeat_files(password,scpdir=scpdir,searchstr=searchstr,fieldnameunderscore=False)

    # - - - - Sourcecats - - - -
    scpdir    = '/Users/kschmidt/work/MUSE/MUSEWide_CFDS44release/repeat_sourcecats_170714/'
    searchstr = 'tdose_sourcecats/catalog_photometry_*FFFFF*_tdose_sourcecat_id*IIIII*_cutout5p0x5p0arcsec.fits'
    mw44.get_repeat_files(password,scpdir=scpdir,searchstr=searchstr,fieldnameunderscore=True)

    """
    archedir = '/store/data/musewide/TDOSE/170714_cdfs44fields_TDOSEextraction/'
    if verbose: print(' - Locating files in\n   '+archedir+'\n   on arche.aip.de and copying them to\n   '+scpdir)
    iddic    = mw44.get_repeast_ids()

    for fieldname in iddic.keys():
        idlist = iddic[fieldname]
        for objid in idlist:
            filestring = searchstr.replace('IIIII',str(objid))
            if fieldnameunderscore:
                filestring = filestring.replace('FFFFF',fieldname.replace('-','_'))
            else:
                filestring = filestring.replace('FFFFF',fieldname)
            scpout     = commands.getoutput('sshpass -p '+password+' scp -r kasper@arche.aip.de:'+archedir+filestring+'  '+scpdir)
            if verbose & (scpout != ''): print(scpout)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_repeast_ids(idlist='/Users/kschmidt/work/MUSE/MUSEWide_CFDS44release/repeat_tdose_kasper.txt'):
    """

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import MUSEWideCDFS44release as mw44
    iddic = mw44.get_repeast_ids()

    """
    iddic     = {}
    fieldname = 'None'

    for idline in open(idlist,'r').readlines():
        if idline.startswith('#'):
            pass
        elif idline == ('\n'):
            pass
        elif idline.startswith('c'):
            fieldname = idline.split('\n')[0]
            iddic[fieldname] = []
        else:
            try:
                iddic[fieldname].append(int(idline.split(' ')[0].split('\n')[0]))
            except:
                pdb.set_trace()

    return iddic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def prepare_repeat_setup():
    """
    Generate the TDOSE setups based on an infofile

    --- INPUT ---

    --- EXAMPLE OF USE ---
    mw44.prepare_repeat_setup()

    """
    infofile  = '/Users/kschmidt/work/MUSE/MUSEWide_CFDS44release/musewide_infofile_DATABCKUP1.txt'
    outputdir = '/Volumes/DATABCKUP1/TDOSEextractions/tdose_setupfiles/'
    tu.duplicate_setup_template(outputdir,infofile,namebase='MUSEWide_CDFS44repeat_tdose_setup_DATABCKUP1_candels',
                                loopcols='all',clobber=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def prepare_repeat_photocat():
    """
    convert the SExtractor catalogs into fits catalogs to use them for initial guesses of the Gaussian modeling:

    --- INPUT ---

    --- EXAMPLE OF USE ---
    mw44.prepare_repeat_photocat()

    """
    catalogs = glob.glob('/Volumes/DATABCKUP1/MUSE-Wide/catalogs_photometry/catalog_photometry_candels-cdfs-*.cat')
    tu.SExtractorCat2fits(catalogs,stringcols=[1],header=73,verbose=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_repeat_sourcecats():
    """
    Generate the source catalogs

    --- INPUT ---

    --- EXAMPLE OF USE ---
    mw44.gen_repeat_sourcecats()

    """

    catfiles = glob.glob('/Volumes/DATABCKUP1/MUSE-Wide/catalogs_photometry/catalog_photometry_candels*.cat')
    for catfile in catfiles:
        field = catfile.split('/')[-1].split('_')[-1].split('.ca')[0]
        output = '/Volumes/DATABCKUP1/TDOSEextractions/tdose_sourcecats/'+catfile.split('/')[-1].split('.')[0]+'_tdose_sourcecat.txt'
        output = output.replace('candels-cdfs-','candels_cdfs_')

        if (field == 'candels-cdfs-60') or (field == 'candels-cdfs-61'):
            imgfile = '/Volumes/DATABCKUP1/MUSE-Wide/hst_cutouts/acs_814w_'+field+'_cut_v1.0.fits'
        else:
            imgfile = '/Volumes/DATABCKUP1/MUSE-Wide/hst_cutouts/acs_814w_'+field+'_cut_v2.0.fits'

        if os.path.isfile(imgfile):
            print ' - Creating source catalog for '+field+' using: \n   '+catfile+' \n   '+imgfile
            imgheader  = pyfits.open(imgfile)[0].header
            sourcecat  = tu.gen_sourcecat_from_SExtractorfile(catfile,output,imgheader=imgheader,idcol=0,racol=2,deccol=3,
                                                              fluxcol=22,fluxfactor=100.,generateDS9reg=True,verbose=False,
                                                              clobber=True)
        else:
            print ' - WARNING Imgfile not found ('+imgfile+'); skipping generating catalog for '+field

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_repeat_IDlists(infofile='/Users/kschmidt/work/MUSE/MUSEWide_CFDS44release/musewide_infofile_DATABCKUP1.txt',
                       verbose=True):
    """
    Print IDs to be inserted into the setup info file used to generate the setup files with prepare_repeat_setup

    --- INPUT ---
    infofile      used to get the list of fields to return ID lists for

    --- EXAMPLE OF USE ---
    mw44.get_repeat_IDlists()

    """
    iddic = mw44.get_repeast_ids()

    for infoline in open(infofile,'r').readlines():
        if infoline.startswith('#'):
            pass
        else:
            fieldname  = infoline.split()[0]
            try:
                idlist = iddic[fieldname]
            except:
                idlist = 'IDlistNotFound'
            idliststr = fieldname.ljust(20)+str(idlist).replace(', ',',')
            if verbose: print (idliststr.ljust(100)+'END')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def run_repeat_TDOSE(skipextractedobjects=True,gencutputsandscats=True,verbose=True):
    """

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import MUSEWideCDFS44release as mw44
    mw44.run_repeat_TDOSE(skipextractedobjects=True,gencutputsandscats=False)


    """
    setupdir   = '/Volumes/DATABCKUP1/TDOSEextractions/tdose_setupfiles/'
    setupfiles = glob.glob(setupdir+'MUSEWide_CDFS44repeat_tdose_setup_DATABCKUP1_candels_candels-cdfs-*.txt')

    files2run  = setupfiles[-2:]

    for sfile in files2run:
        tdose.perform_extraction(setupfile=sfile, modelrefimage=True,refimagemodel2cubewcs=True,definePSF=True,
                                 modeldatacube=True,createsourcecube=True,store1Dspectra=True,plot1Dspectra=True,
                                 plotS2Nspectra=True,save_init_model_output=False,clobber=True,verbose=True,verbosefull=True,
                                 logterminaloutput=False,skipextractedobjects=skipextractedobjects,skipspecificobjects=None,
                                 performcutout=gencutputsandscats,generatesourcecat=gencutputsandscats)

    for sfile in files2run:
        tu.gen_overview_plot('all',sfile,skipobj=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =