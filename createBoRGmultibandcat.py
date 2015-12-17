#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# createBoRGmultibandcat.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Turning the individual SExtractor catalogs of the BoRG fields into
# multiband catalogs on the format required by the BPZ scripts
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# parentdir        : The path to the top level directory where all the catalogs recide.
#                    The individual SExtractor catalogs should be in subdirectories named
#                    after the fields. These excpected field names are hardcoded below.
#                    End string with '/'
# outputdir        : Directory to put the resulting *multiband.cat in (assumes that it 
#                    already excists) - end string with '/'
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --nofieldprefix  : Keyword to be set when final_sources*.cat files have now field name as
#                    prefix. This is the case when the catalogs have just been created with
#                    the process_borg_field.sh script
# --F105W          : set --F105W to add columns of F105W data as new columns to the end of output
# --verbose        : set --verbose to get info/messages printed to the screen
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# 
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x createBoRGmultibandcat.py       (only required once)
# bash> createBoRGmultibandcat.py '/Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p1/' '/Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p1/multiband/' --verbose 

# bash> createBoRGmultibandcat.py /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/ /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2_forcephot/multiband/ --F105W --nofieldprefix --verbose 

#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-11-06  started by K. B. Schmidt (UCSB)
# 2013-01-06  added --nofieldprefix and fixed division with 0. K. B. Schmidt (UCSB)
# 2013-04-17  added --F105W keyword K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import kbsutilities as kbs #
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import commands            # get output from spawned command line processes
import getopt              # used to extract/obtain the optional input
import numpy as np         # enable opening with genfromtxt
import matplotlib as plt
import pdb                 # for debugging
#----------------------------
#   FUNCTIONS
#----------------------------
#
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("parentdir", type=str, help="Path of parent directory for data (end with '/')")
parser.add_argument("outputdir", type=str, help="Excisting output directory for new catalogs (end with '/')")
# ---- optional arguments ----
parser.add_argument("--nofieldprefix", action="store_true", help="Set this keyword if the final_sources* files have no field prefix in maes.")
parser.add_argument("--F105W", action="store_true", help="Add F105W data to output")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# list of expected subdirectories to search through
specialcase = 1  # use to distinguish between the standard run and other runs
if specialcase == 0: 
    fields  = ['borg_0110-0224','borg_0214+1255','borg_0228-4102','borg_0240-1857','borg_0427+2538','borg_0436-5259','borg_0439-5317','borg_0440-5244','borg_0540-6409','borg_0553-6405','borg_0624-6432','borg_0624-6440','borg_0637-7518','borg_0751+2917','borg_0756+3043','borg_0808+3946','borg_0819+4911','borg_0820+2332','borg_0835+2456','borg_0846+7654','borg_0906+0255','borg_0909+0002','borg_0914+2822','borg_0922+4505','borg_0926+4000','borg_0926+4426','borg_1010+3001','borg_1014-0423','borg_1031+3804','borg_1031+5052','borg_1033+5051','borg_1051+3359','borg_1103-2330','borg_1111+5545','borg_1119+4026','borg_1131+3114','borg_1152+5441','borg_1153+0056','borg_1209+4543','borg_1230+0750','borg_1242+5716','borg_1245+3356','borg_1301+0000','borg_1337+0028','borg_1341+4123','borg_1408+5503','borg_1437+5043','borg_1510+1115','borg_1524+0954','borg_1555+1108','borg_1632+3733','borg_1632+3737','borg_1815-3244','borg_2057-4412','borg_2132+1004','borg_2155-4411','borg_2203+1851','borg_2345+0054','borg_2351-4332'] # version 2.0

    Nfields = len(fields)

    str606 = '_final_sources_cleaned_F606.cat'
    str600 = '_final_sources_cleaned_F600.cat'
    str098 = '_final_sources_cleaned_F098.cat'
    str125 = '_final_sources_cleaned_F125.cat'
    str160 = '_final_sources_cleaned_F160.cat'
else:
#    fields  = ['borg_1033+5051']
    fields  = ['borg_1437+5043']
#    fields  = ['borg_0456-2203','borg_0951+3304','borg_1059+0519','borg_1118-1858','borg_1358+4334','borg_1459+7146','borg_1510+1115','borg_2132-1202','borg_2313-2243'] # Cycle 19
#    fields  = ['borg_0952+5304','borg_1358+4326','borg_1416+1638','borg_1429-0331'] # ramining Cycle 19 cats

    Nfields = len(fields)
    str606 = '_final_sources_F606.cat'
    str600 = '_final_sources_F600.cat'
    str098 = '_final_sources_F098.cat'
    str125 = '_final_sources_F125.cat'
    str160 = '_final_sources_F160.cat'
    if args.F105W: str105 = '_final_sources_F105.cat'

#-------------------------------------------------------------------------------------------------------------
# list of expected subdirectories to search through
if args.verbose: print ':: '+sys.argv[0]+' :: Will loop over and search in '+str(Nfields)+' subdirectories in '+args.parentdir

for ii in range(Nfields):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # reading the V band catalog
    catfileVband = args.parentdir+fields[ii]+'/'+fields[ii]+str606           # catalogfile to look for
    if args.nofieldprefix: catfileVband = args.parentdir+fields[ii]+'/'+str606[1:] # catalogfile to look for
    if os.path.isfile(catfileVband):                                         # if catalog exists then read it
        catVband     = np.genfromtxt(catfileVband, dtype=None, comments='#')
        NobjVband    = len(catVband['f0'])
        colheadVband = 'f606w'
    else:
        catfileVband = args.parentdir+fields[ii]+'/'+fields[ii]+str600        # catalogfile to look for
        if args.nofieldprefix: catfileVband = args.parentdir+fields[ii]+'/'+str600[1:] # catalogfile to look for
        if os.path.isfile(catfileVband):                                      # if catalog exists then read it
            catVband     = np.genfromtxt(catfileVband, dtype=None, comments='#')
            NobjVband    = len(catVband['f0'])
            colheadVband = 'f600lp'
        else:                                                                           # otherwise abort
            if args.verbose: print ':: '+sys.argv[0]+' :: ERROR Cannot find the file '+catfileVband+' (wanna use --nofieldprefix?) --> ABORTING' 
            sys.exit()
    print ':: '+sys.argv[0]+' :: Read: '+catfileVband
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # reading the Y band catalog
    catfileYband = args.parentdir+fields[ii]+'/'+fields[ii]+str098           # catalogfile to look for
    if args.nofieldprefix: catfileYband = args.parentdir+fields[ii]+'/'+str098[1:] # catalogfile to look for
    if os.path.isfile(catfileYband):                                         # if catalog exists then read it
        catYband     = np.genfromtxt(catfileYband, dtype=None, comments='#')
        NobjYband    = len(catYband['f0'])
        colheadYband = 'f098m'
    else:
        if args.verbose: print ':: '+sys.argv[0]+' :: ERROR Cannot find the file '+catfileYband+' (wanna use --nofieldprefix?) --> ABORTING' 
        sys.exit()
    print ':: '+sys.argv[0]+' :: Read: '+catfileYband
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # reading the J band catalog
    catfileJband = args.parentdir+fields[ii]+'/'+fields[ii]+str125           # catalogfile to look for
    if args.nofieldprefix: catfileJband = args.parentdir+fields[ii]+'/'+str125[1:] # catalogfile to look for
    if os.path.isfile(catfileJband):                                         # if catalog exists then read it
        catJband     = np.genfromtxt(catfileJband, dtype=None, comments='#')
        NobjJband    = len(catJband['f0'])
        colheadJband = 'f125w'
    else:
        if args.verbose: print ':: '+sys.argv[0]+' :: ERROR Cannot find the file '+catfileJband+' (wanna use --nofieldprefix?) --> ABORTING' 
        sys.exit()
    print ':: '+sys.argv[0]+' :: Read: '+catfileJband
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # reading the H band catalog
    catfileHband = args.parentdir+fields[ii]+'/'+fields[ii]+str160           # catalogfile to look for
    if args.nofieldprefix: catfileHband = args.parentdir+fields[ii]+'/'+str160[1:] # catalogfile to look for
    if os.path.isfile(catfileHband):                                         # if catalog exists then read it
        catHband     = np.genfromtxt(catfileHband, dtype=None, comments='#')
        NobjHband    = len(catHband['f0'])
        colheadHband = 'f160w'
    else:
        if args.verbose: print ':: '+sys.argv[0]+' :: ERROR Cannot find the file '+catfileHband+' (wanna use --nofieldprefix?) --> ABORTING' 
        sys.exit()
    print ':: '+sys.argv[0]+' :: Read: '+catfileHband
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    NobjF105W = NobjVband # default Number of objects for F105W filter (ensures no break if --F105W not given)
    if args.F105W:
        # reading the F105W catalog
        catfileF105W = args.parentdir+fields[ii]+'/'+fields[ii]+str105           # catalogfile to look for
        if args.nofieldprefix: catfileF105W = args.parentdir+fields[ii]+'/'+str105[1:] # catalogfile to look for
        if os.path.isfile(catfileF105W):                                         # if catalog exists then read it
            catF105W     = np.genfromtxt(catfileF105W, dtype=None, comments='#')
            NobjF105W    = len(catF105W['f0'])
            colheadF105W = 'f105w'
        else:
            if args.verbose: print ':: '+sys.argv[0]+' :: ERROR Cannot find the file '+catfileF105W+' (wanna use --nofieldprefix?) --> ABORTING' 
            sys.exit()
        print ':: '+sys.argv[0]+' :: Read: '+catfileF105W
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Checking that there is the same number of objects in all catalogs
    if (NobjVband == NobjYband) and (NobjVband == NobjJband) and (NobjVband == NobjHband) and (NobjVband == NobjF105W):
        if args.verbose: print ':: '+sys.argv[0]+' :: Found '+str(NobjVband)+' objects in all catalogs'
    else:
        if args.verbose: print ':: '+sys.argv[0]+' :: ERROR The number of objects in the different bands are different --> ABORTING' 
        sys.exit()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Defining the various columns of the output catalog 
    id                  = catJband['f0']
    cid                 = []
    for jj in range(NobjJband):
        cid.append(fields[ii]+'_'+str(id[jj]).zfill(4))
    x                   = catJband['f31']
    y                   = catJband['f32']
    ra                  = catJband['f33']
    dec                 = catJband['f34']
    kron_radius         = catJband['f29']
    a_image             = catJband['f35']
    b_image             = catJband['f36']
    theta_image         = catJband['f37']
    fwhm_image          = catJband['f38']
    fwhm_world          = catJband['f39']
    stellarity          = catJband['f41']
    seflags             = catJband['f40']
    # ------------------------------- V band ------------------------------- 
    f606w_bkgrd         = catVband['f30']
    f606w_isoflux       = catVband['f1']
    f606w_isofluxerr    = catVband['f2']
    f606w_isosnr        = kbs.divideWith0(f606w_isoflux,f606w_isofluxerr,-99)
    f606w_autoflux      = catVband['f25']
    f606w_autofluxerr   = catVband['f26']
    f606w_autosnr       = kbs.divideWith0(f606w_autoflux,f606w_autofluxerr,-99)
    f606w_isomag        = catVband['f3']
    f606w_isomagerr     = catVband['f4']
    f606w_automag       = catVband['f27']
    f606w_automagerr    = catVband['f28']
    # ------------------------------- Y band ------------------------------- 
    f098m_bkgrd         = catYband['f30']
    f098m_isoflux       = catYband['f1']
    f098m_isofluxerr    = catYband['f2']
    f098m_isosnr        = kbs.divideWith0(f098m_isoflux,f098m_isofluxerr,-99)
    f098m_autoflux      = catYband['f25']
    f098m_autofluxerr   = catYband['f26']
    f098m_autosnr       = kbs.divideWith0(f098m_autoflux,f098m_autofluxerr,-99)
    f098m_isomag        = catYband['f3']
    f098m_isomagerr     = catYband['f4']
    f098m_automag       = catYband['f27']
    f098m_automagerr    = catYband['f28']
    # ------------------------------- J band ------------------------------- 
    f125w_bkgrd         = catJband['f30']
    f125w_isoflux       = catJband['f1']
    f125w_isofluxerr    = catJband['f2']
    f125w_isosnr        = kbs.divideWith0(f125w_isoflux,f125w_isofluxerr,-99)
    f125w_autoflux      = catJband['f25']
    f125w_autofluxerr   = catJband['f26']
    f125w_autosnr       = kbs.divideWith0(f125w_autoflux,f125w_autofluxerr,-99)
    f125w_isomag        = catJband['f3']
    f125w_isomagerr     = catJband['f4']
    f125w_automag       = catJband['f27']
    f125w_automagerr    = catJband['f28']
    # ------------------------------- H band ------------------------------- 
    f160w_bkgrd         = catHband['f30']
    f160w_isoflux       = catHband['f1']
    f160w_isofluxerr    = catHband['f2']
    f160w_isosnr        = kbs.divideWith0(f160w_isoflux,f160w_isofluxerr,-99)
    f160w_autoflux      = catHband['f25']
    f160w_autofluxerr   = catHband['f26']
    f160w_autosnr       = kbs.divideWith0(f160w_autoflux,f160w_autofluxerr,-99)
    f160w_isomag        = catHband['f3']
    f160w_isomagerr     = catHband['f4']
    f160w_automag       = catHband['f27']
    f160w_automagerr    = catHband['f28']
    # ------------------------------- F105W ------------------------------- 
    if args.F105W:
        f105w_bkgrd         = catF105W['f30']
        f105w_isoflux       = catF105W['f1']
        f105w_isofluxerr    = catF105W['f2']
        f105w_isosnr        = kbs.divideWith0(f105w_isoflux,f105w_isofluxerr,-99)
        f105w_autoflux      = catF105W['f25']
        f105w_autofluxerr   = catF105W['f26']
        f105w_autosnr       = kbs.divideWith0(f105w_autoflux,f105w_autofluxerr,-99)
        f105w_isomag        = catF105W['f3']
        f105w_isomagerr     = catF105W['f4']
        f105w_automag       = catF105W['f27']
        f105w_automagerr    = catF105W['f28']
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Writing the new catalog
    outname    = args.outputdir+fields[ii]+'_multiband.cat'                  # name of output catalog
    outputcat  = open(outname,"w")

    outputcat.write('# Merging the following four catalogs together. \n')
    outputcat.write('# '+catfileVband+' \n')
    outputcat.write('# '+catfileYband+' \n')
    outputcat.write('# '+catfileJband+' \n')
    outputcat.write('# '+catfileHband+' \n')
    if args.F105W: outputcat.write('# '+catfileF105W+' \n')
    outputcat.write('# \n')
    outputcat.write('# Columns 3-14 are taken from catfileJband \n')
    outputcat.write('#   1 id                   \n')
    outputcat.write('#   2 cid                  \n')
    outputcat.write('#   3 x                    \n')
    outputcat.write('#   4 y                    \n')
    outputcat.write('#   5 ra                   \n')
    outputcat.write('#   6 dec                  \n')
    outputcat.write('#   7 kron_radius          \n')
    outputcat.write('#   8 a_image              \n')
    outputcat.write('#   9 b_image              \n')
    outputcat.write('#  10 theta_image          \n')
    outputcat.write('#  11 fwhm_image           \n')
    outputcat.write('#  12 fwhm_world           \n')
    outputcat.write('#  13 stellarity           \n')
    outputcat.write('#  14 seflags              \n')
    outputcat.write('#  15 '+colheadVband+'_bkgrd         \n')
    outputcat.write('#  16 '+colheadVband+'_isoflux       \n')
    outputcat.write('#  17 '+colheadVband+'_isofluxerr    \n')
    outputcat.write('#  18 '+colheadVband+'_isosnr        \n')
    outputcat.write('#  19 '+colheadVband+'_autoflux      \n')
    outputcat.write('#  20 '+colheadVband+'_autofluxerr   \n')
    outputcat.write('#  21 '+colheadVband+'_autosnr       \n')
    outputcat.write('#  22 '+colheadVband+'_isomag        \n')
    outputcat.write('#  23 '+colheadVband+'_isomagerr     \n')
    outputcat.write('#  24 '+colheadVband+'_automag       \n')
    outputcat.write('#  25 '+colheadVband+'_automagerr    \n')
    outputcat.write('#  26 f098m_bkgrd          \n')
    outputcat.write('#  27 f098m_isoflux        \n')
    outputcat.write('#  28 f098m_isofluxerr     \n')
    outputcat.write('#  29 f098m_isosnr         \n')
    outputcat.write('#  30 f098m_autoflux       \n')
    outputcat.write('#  31 f098m_autofluxerr    \n')
    outputcat.write('#  32 f098m_autosnr        \n')
    outputcat.write('#  33 f098m_isomag         \n')
    outputcat.write('#  34 f098m_isomagerr      \n')
    outputcat.write('#  35 f098m_automag        \n')
    outputcat.write('#  36 f098m_automagerr     \n')
    outputcat.write('#  37 f125w_bkgrd          \n')
    outputcat.write('#  38 f125w_isoflux        \n')
    outputcat.write('#  39 f125w_isofluxerr     \n')
    outputcat.write('#  40 f125w_isosnr         \n')
    outputcat.write('#  41 f125w_autoflux       \n')
    outputcat.write('#  42 f125w_autofluxerr    \n')
    outputcat.write('#  43 f125w_autosnr        \n')
    outputcat.write('#  44 f125w_isomag         \n')
    outputcat.write('#  45 f125w_isomagerr      \n')
    outputcat.write('#  46 f125w_automag        \n')
    outputcat.write('#  47 f125w_automagerr     \n')
    outputcat.write('#  48 f160w_bkgrd          \n')
    outputcat.write('#  49 f160w_isoflux        \n')
    outputcat.write('#  50 f160w_isofluxerr     \n')
    outputcat.write('#  51 f160w_isosnr         \n')
    outputcat.write('#  52 f160w_autoflux       \n')
    outputcat.write('#  53 f160w_autofluxerr    \n')
    outputcat.write('#  54 f160w_autosnr        \n')
    outputcat.write('#  55 f160w_isomag         \n')
    outputcat.write('#  56 f160w_isomagerr      \n')
    outputcat.write('#  57 f160w_automag        \n')
    outputcat.write('#  58 f160w_automagerr \n')
    if args.F105W:
        outputcat.write('#  59 f105w_bkgrd          \n')
        outputcat.write('#  60 f105w_isoflux        \n')
        outputcat.write('#  61 f105w_isofluxerr     \n')
        outputcat.write('#  62 f105w_isosnr         \n')
        outputcat.write('#  63 f105w_autoflux       \n')
        outputcat.write('#  64 f105w_autofluxerr    \n')
        outputcat.write('#  65 f105w_autosnr        \n')
        outputcat.write('#  66 f105w_isomag         \n')
        outputcat.write('#  67 f105w_isomagerr      \n')
        outputcat.write('#  68 f105w_automag        \n')
        outputcat.write('#  69 f105w_automagerr \n')
    outputcat.write('# \n')

    columnstring = '# id   cid   x   y   ra   dec   kron_radius   a_image   b_image   theta_image   fwhm_image   fwhm_world   stellarity   seflags   '+colheadVband+'_bkgrd   '+colheadVband+'_isoflux   '+colheadVband+'_isofluxerr   '+colheadVband+'_isosnr   '+colheadVband+'_autoflux   '+colheadVband+'_autofluxerr   '+colheadVband+'_autosnr   '+colheadVband+'_isomag   '+colheadVband+'_isomagerr   '+colheadVband+'_automag   '+colheadVband+'_automagerr   f098m_bkgrd   f098m_isoflux   f098m_isofluxerr   f098m_isosnr   f098m_autoflux   f098m_autofluxerr   f098m_autosnr   f098m_isomag   f098m_isomagerr   f098m_automag   f098m_automagerr   f125w_bkgrd   f125w_isoflux   f125w_isofluxerr   f125w_isosnr   f125w_autoflux   f125w_autofluxerr   f125w_autosnr   f125w_isomag   f125w_isomagerr   f125w_automag   f125w_automagerr   f160w_bkgrd   f160w_isoflux   f160w_isofluxerr   f160w_isosnr   f160w_autoflux   f160w_autofluxerr   f160w_autosnr   f160w_isomag   f160w_isomagerr   f160w_automag   f160w_automagerr \n'

    if args.F105W: columnstring.replace(' \n','   f105w_bkgrd   f105w_isoflux   f105w_isofluxerr   f105w_isosnr   f105w_autoflux   f105w_autofluxerr   f105w_autosnr   f105w_isomag   f105w_isomagerr   f105w_automag   f105w_automagerr \n')

    outputcat.write(columnstring)

    for kk in range(NobjJband):
        objectstring = str(id[kk])+'  '+cid[kk]+'  '+str(x[kk])+'  '+str(y[kk])+'  '+str(ra[kk])+'  '+str(dec[kk])+'  '+str(kron_radius[kk])+'  '+str(a_image[kk])+'  '+str(b_image[kk])+'  '+str(theta_image[kk])+'  '+str(fwhm_image[kk])+'  '+str(fwhm_world[kk])+'  '+str(stellarity[kk])+'  '+str(seflags[kk])+'  '+str(f606w_bkgrd[kk])+'  '+str(f606w_isoflux[kk])+'  '+str(f606w_isofluxerr[kk])+'  '+str(f606w_isosnr[kk])+'  '+str(f606w_autoflux[kk])+'  '+str(f606w_autofluxerr[kk])+'  '+str(f606w_autosnr[kk])+'  '+str(f606w_isomag[kk])+'  '+str(f606w_isomagerr[kk])+'  '+str(f606w_automag[kk])+'  '+str(f606w_automagerr[kk])+'  '+str(f098m_bkgrd[kk])+'  '+str(f098m_isoflux[kk])+'  '+str(f098m_isofluxerr[kk])+'  '+str(f098m_isosnr[kk])+'  '+str(f098m_autoflux[kk])+'  '+str(f098m_autofluxerr[kk])+'  '+str(f098m_autosnr[kk])+'  '+str(f098m_isomag[kk])+'  '+str(f098m_isomagerr[kk])+'  '+str(f098m_automag[kk])+'  '+str(f098m_automagerr[kk])+'  '+str(f125w_bkgrd[kk])+'  '+str(f125w_isoflux[kk])+'  '+str(f125w_isofluxerr[kk])+'  '+str(f125w_isosnr[kk])+'  '+str(f125w_autoflux[kk])+'  '+str(f125w_autofluxerr[kk])+'  '+str(f125w_autosnr[kk])+'  '+str(f125w_isomag[kk])+'  '+str(f125w_isomagerr[kk])+'  '+str(f125w_automag[kk])+'  '+str(f125w_automagerr[kk])+'  '+str(f160w_bkgrd[kk])+'  '+str(f160w_isoflux[kk])+'  '+str(f160w_isofluxerr[kk])+'  '+str(f160w_isosnr[kk])+'  '+str(f160w_autoflux[kk])+'  '+str(f160w_autofluxerr[kk])+'  '+str(f160w_autosnr[kk])+'  '+str(f160w_isomag[kk])+'  '+str(f160w_isomagerr[kk])+'  '+str(f160w_automag[kk])+'  '+str(f160w_automagerr[kk])+'  \n'

        if args.F105W: 
            repstring = '  '+str(f105w_bkgrd[kk])+'  '+str(f105w_isoflux[kk])+'  '+str(f105w_isofluxerr[kk])+'  '+str(f105w_isosnr[kk])+'  '+str(f105w_autoflux[kk])+'  '+str(f105w_autofluxerr[kk])+'  '+str(f105w_autosnr[kk])+'  '+str(f105w_isomag[kk])+'  '+str(f105w_isomagerr[kk])+'  '+str(f105w_automag[kk])+'  '+str(f105w_automagerr[kk])+'  \n'
            objectstring = objectstring.replace(' \n',repstring)

        outputcat.write(objectstring)

    outputcat.close()

#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
