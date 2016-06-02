# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import os
import commands
import sys
import glob
import pyfits
import numpy as np
import ciiiEmitterCandidates as cec
import CDFS_MUSEvs3DHST as cm3
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_candidates(verbose=True):
    """

    Get the sample(s) of candidates to investigate

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    cec.get_candidates()

    """
    catMUSE     = '/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1_MUSEsubcat.fits'
    cat3DHST    = '/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1_3DHSTinfo.fits'
    catMUSEorig = '/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits'


    if verbose: print ' ------------ OBJECT WITH CIII] FALLING IN 3D-HST G141 WAVELENGTH RANGE ------------ '
    cm3.get_candidates(zrange=[4.7,7.9],matchtol=0.5,catMUSE=catMUSE,cat3DHST=cat3DHST)

    if verbose: print ' ------------ OBJECTS WITH CIII] FALLING IN THE MUSE WAVELENGTH RANGE ------------ '
    dat = pyfits.open(catMUSEorig)[1].data
    zrange = [2.8,3.9]
    goodent = np.where((dat['REDSHIFT'] > zrange[0]) & (dat['REDSHIFT'] < zrange[1]))[0]

    if verbose: print ' - Found '+str(len(goodent))+' candidate objects in total'
    print dat['redshift'][goodent]
    for id in dat['UNIQUE_ID'][goodent]: print id

    if verbose: print ' ------------ OBJECTS WITH CIII] FALLING IN THE MUSE WAVELENGTH RANGE (3D-HST MATCH) ------------ '
    cm3.get_candidates(zrange=zrange,matchtol=0.5,catMUSE=catMUSE,cat3DHST=cat3DHST)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def copy_spectra(outputdir='./',verbose=True,
                 specdir_MUSE='/Users/kschmidt/work/MUSE/emission_spectra_0.2/',
                 specdir_3DHST='/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1/'):
    """

    Get the sample(s) of candidates to investigate

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    cec.copy_spectra(outputdir='./spectra_CIIIcandidates/')

    """
    if os.path.isdir(outputdir):
        if verbose: print ' - Copying spectra to '+outputdir
    else:
        if verbose: print ' - The outputdir does not exists; doing nothing. \n   outputdir='+outputdir
        return

    cand3DHST           = '/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/candidateCIIIemitters.txt'
    candMUSE            = '/Users/kschmidt/work/MUSE/notes/160612/MUSE_CIII_potentialemitters.txt'
    candMUSE_3DHSTmatch = '/Users/kschmidt/work/MUSE/notes/160612/MUSE_CIII_potentialemitters_matchto3DHST.txt'

    dat_3DHST            = np.genfromtxt(cand3DHST          ,dtype=None,comments='#',names=True)
    dat_MUSE             = np.genfromtxt(candMUSE           ,dtype=None,comments='#',names=True)
    dat_MUSE_3DHSTmatch  = np.genfromtxt(candMUSE_3DHSTmatch,dtype=None,comments='#',names=True)


    if verbose: print ' - Copying MUSE spectra over'
    ids_MUSEall = np.append(dat_3DHST['ID_MUSE'],dat_MUSE['ID_MUSE'])
    for objid in ids_MUSEall:
        cpcmd = 'cp '+specdir_MUSE+'spectrum*'+str(objid)+'*.fits  '+outputdir
        cpout = commands.getoutput(cpcmd)
        if not cpout == '':
            print cpout

    if verbose: print ' - Reformatting and copying 3D-HST spectra over'
    ids_3DHSTall = np.append(dat_3DHST['ID_3DHST'],dat_MUSE_3DHSTmatch['ID_3DHST'])
    ids_MUSEall  = np.append(dat_3DHST['ID_MUSE'],dat_MUSE_3DHSTmatch['ID_MUSE'])
    for oo, objid_3dhst in enumerate(ids_3DHSTall):
        objid_muse = ids_MUSEall[oo]
        if verbose:
            idno    = oo+1
            infostr = '   Reformat 3D-HST 1D spectrum for '+str(objid_3dhst)+' / '+str(objid_muse)+\
                      '   ('+str(idno)+'/'+str(len(ids_3DHSTall))+')'
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        specs       = glob.glob(specdir_3DHST+'*'+str("%.5d" % objid_3dhst)+'*1D.fits')
        for spec in specs: # looping over the 1D spectra found (from multiple PAs)
            if len(spec) == 0:
                sys.exit('Something is not right; Did not find any spectrum for '+str(objid_3dhst)+
                         ' ('+str(objid_muse)+')\n'+spec)

            reformatname = 'spectrum_'+str(objid_muse)+'_'+spec.split('/')[-1].split('.')[0]+'_MiG1Dreformat.fits'

            cec.prep_3DHST_1Dspec4MiG1D(spec,outname=outputdir+reformatname,verbose=False)

    if verbose: print '\n - Done...'
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def prep_3DHST_1Dspec4MiG1D(specname,outname='DEFAULT',verbose=True):
    """

    Convert a 3DHST *1D.fits spectrum to

    --- INPUT ---
    sepcname         1D grism spectrum to reformat for MiG1D inspection.
    outname          Name (incl. paht) of file to save reformatted spectrum to.

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    newspec = cec.prep_3DHST_1Dspec4MiG1D('../candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1/goodss-31-G141_24503.1D.fits')

    """
    if outname == 'DEFAULT':
        outfile = './'+specname.split('/')[-1].split('.')[0]+'_MiG1Dreformat.fits'
    else:
        outfile = outname

    if verbose: print ' - Copying original spectrum ('+specname+')\n   to file containing output ('+outfile+')'
    cpcmd = 'cp '+specname+' '+outfile
    cpout = commands.getoutput(cpcmd)
    if not cpout == '':
        print cpout

    if verbose: print ' - Loading data and reformating:'
    hdulist = pyfits.open(outfile, mode='update')

    if verbose: print '   o Change header keywords'
    hdulist[1].header['TTYPE1']   = 'WAVE_AIR'

    hdulist[1].header['TTYPE2']   = 'FLUX'
    hdulist[1].header['TUNIT2']   = '1e-20 erg/s/cm2/A'

    hdulist[1].header['TTYPE3']   = 'FLUXERR'
    hdulist[1].header['TUNIT3']   = '1e-20 erg/s/cm2/A'

    hdulist[1].header['TTYPE4']   = 'CONTAM'
    hdulist[1].header['TUNIT4']   = '1e-20 erg/s/cm2/A'

    hdulist[1].header['TTYPE5']   = 'TRACE'
    hdulist[1].header['TUNIT5']   = '1e-20 erg/s/cm2/A'

    hdulist[1].header['TTYPE6']   = 'ETRACE'
    hdulist[1].header['TUNIT6']   = '1e-20 erg/s/cm2/A'

    hdulist[1].header['TTYPE7']   = 'SENSITIVITY'
    hdulist[1].header['TUNIT7']   = '1e-20 erg/s/cm2/A  /  E/S'

    if verbose: print '   o Adjust flux levels to be in cgs/[1e-20] units'
    hdulist[1].data['FLUX']       = hdulist[1].data['FLUX']/hdulist[1].data['SENSITIVITY']    *1.0e3
    hdulist[1].data['FLUXERR']    = hdulist[1].data['FLUXERR']/hdulist[1].data['SENSITIVITY'] *1.0e3
    hdulist[1].data['TRACE']      = hdulist[1].data['TRACE']/hdulist[1].data['SENSITIVITY']   *1.0e3
    hdulist[1].data['ETRACE']     = hdulist[1].data['ETRACE']/hdulist[1].data['SENSITIVITY']  *1.0e3

    if verbose: print '   o Correct sensititivity conversion'
    hdulist[1].data['SENSITIVITY']= 1.0e3 / hdulist[1].data['SENSITIVITY']

    if verbose: print ' - Storing correctiosn to output file'
    hdulist.flush()
    return outfile
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def inspect_3DHST_2Dspec(matchcatalog='/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/candidateCIIIemitters.txt',
                         smooth=True,verbose=True):
    """

    Inspect the 2D grism spectra

    --- INPUT ---
    matchcatalog    A catalog with MUSE objects matched to 3D-HST; expects the format printed
                    to the screen by cm3.get_candidates()
    smooth          Smooth the 2D spectra in DS9 with the default Gaussian kernel?
    verbose         Toggle verbosity

    """
    specdir = '/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1/'
    cm3.inspect_G141specs(matchcatalog,specdir=specdir,smooth=smooth,
                          ds9circlename='CIII]',verbose=verbose,oneobj=False)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =