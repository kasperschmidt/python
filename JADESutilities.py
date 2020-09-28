# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#            Utilities for handling the JADES mock catalogs by Williams et al. (2018)
#                             http://fenrir.as.arizona.edu/jwstmock/
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import JADESutilities as ju
import glob
import os
import datetime
import numpy as np
import sys
import pdb
from astropy.io import fits
import matplotlib.pyplot as plt

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_JADESspecAndInfo(JADESid,observedframe=True,showspec=False,verbose=True):
    """
    Assemble JADES spectrum and info

    --- INPUT ---
    JADESid        JADES mock object ID (Williams et al. 2018) to return info and spectrum for
    observedframe  Return the spectrum in observed frame (shifting lambda and scaling Flambda)?
    showspec       To show the JADES spectrum (in puplot window), set to true.
    verbose        Toggle verbosity

    --- EXAMPLE OF USE ---
    import JADESutilities as ju

    # 274533 z=6.255485 F140Wmag = 25.4527

    JADESinfo, spec_lam, spec_flux = ju.get_JADESspecAndInfo(189021,verbose=True)

    """
    JADESdir     = '/Users/kschmidt/work/catalogs/JADES_GTO/'
    JADEScat     = fits.open(JADESdir+'JADES_SF_mock_r1_v1.0.fits')[1].data

    objent       = np.where(JADEScat['ID'] == JADESid)[0]
    JADESinfo    = JADEScat[objent]

    if (JADESid > 0) & (JADESid < 50002):
        specfile = 'JADES_SF_mock_r1_v1.0_spec_5A_ID_1_50001.fits'
    elif (JADESid > 50001) & (JADESid < 100003):
        specfile = 'JADES_SF_mock_r1_v1.0_spec_5A_ID_50002_100002.fits'
    elif (JADESid > 100002) & (JADESid < 150004):
        specfile = 'JADES_SF_mock_r1_v1.0_spec_5A_ID_100003_150003.fits'
    elif (JADESid > 150003) & (JADESid < 200005):
        specfile = 'JADES_SF_mock_r1_v1.0_spec_5A_ID_150004_200004.fits'
    elif (JADESid > 200004) & (JADESid < 250013):
        specfile = 'JADES_SF_mock_r1_v1.0_spec_5A_ID_200005_250012.fits'
    elif (JADESid > 250012) & (JADESid < 300018):
        specfile = 'JADES_SF_mock_r1_v1.0_spec_5A_ID_250013_300017.fits'
    elif (JADESid > 300018) & (JADESid < 302515):
        specfile = 'JADES_SF_mock_r1_v1.0_spec_5A_ID_300018_302514.fits'
    else:
        sys.exit(' Provided ID ('+str(JADESid)+') is outside JADES id range ')

    JADEShdu     = fits.open(JADESdir+specfile)
    objent       = np.where(JADEShdu['OBJECT PROPERTIES'].data['ID'] == JADESid)
    spec_lam     = JADEShdu['FULL SED WL'].data
    spec_flux    = JADEShdu['FULL SED'].data[objent][0]

    if observedframe:
        spec_lam     = spec_lam  * (1 + JADESinfo['redshift'])
        spec_flux    = spec_flux / (1 + JADESinfo['redshift'])

    HSTF140Wmag  = -2.5 * np.log10(JADESinfo['HST_F140W_fnu']/1e9) + 8.90

    if verbose:
        print(' - Returning info and spec for JADESid = '+str(JADESid)+
              ' @ z = '+str("%.4f" % JADESinfo['redshift'])+
              ' with HST F140W AB mag = '+str("%.2f" % HSTF140Wmag))

    if showspec:
        import pylab as plt
        plt.plot(spec_lam,spec_flux)
        plt.show()

    return JADESinfo, spec_lam, spec_flux

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_JADESobjects(redshift=[3.2,3.6],mag_f140w=[23.5,24.5],MUV=None,mStar=None,jadesinfo=None,verbose=True):
    """
    Function to return JADES mock objects given a set of parameter ranges.

    If one (and only one) paramter has -99 as the upper bound, all other selections will be performed first,
    and then the best match to that single parameter, among the remaining objects is returned. E.g. the
    selections redshift=[3.4,3.45],mag_f140w=[23.5,-99] will return the object with the F140W magnitude
    closest to 23.5 among all objects between redshift 3.4 and 3.45.


    --- INPUT ---
    redshift         Range of redshifts to return objects for. If None all values are accepted
                     If redshift[1]=-99 best match to redshift[0] is returned
    mag_f140w        Range of HST F140W AB mags to return objects for. If None all values are accepted
                     If mag_f140w[1]=-99 best match to mag_f140w[0] is returned
    MUV              Range of MUV to return objects for. If None all values are accepted
                     If MUV[1]=-99 best match to MUV[0] is returned
    mStar            Range of mStar to return objects for. If None all values are accepted
                     If mStar[1]=-99 best match to mStar[0] is returned
    jadesinfo        If JADES catalog should not be loaded in-function, provide the catalog in here
    verbose          Toggle verbosity

    --- EXAMPLE OF USE ---
    import JADESutilities as ju
    JADESinfo = ju.get_JADESobjects(redshift=[3.2,3.6],mag_f140w=[23.5,24.5],MUV=None,mStar=None,jadesinfo=None,verbose=True)

    JADESinfo = ju.get_JADESobjects(redshift=[3.4,-99],mag_f140w=[23.5,24.5],MUV=None,mStar=None,jadesinfo=None,verbose=True)

    JADESinfo = ju.get_JADESobjects(redshift=[3.4,3.45],mag_f140w=[23.5,-99],MUV=None,mStar=None,jadesinfo=None,verbose=True)

    """
    if jadesinfo is None:
        JADESdir     = '/Users/kschmidt/work/catalogs/JADES_GTO/'
        jadesinfo    = fits.open(JADESdir+'JADES_SF_mock_r1_v1.0.fits')[1].data

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if mag_f140w is not None:
        if verbose: print(' - Converting F140W magnitudes to nJu fluxes')
        HST_F140W_fnu    = [0,0]
        if mag_f140w[1] == -99:
            HST_F140W_fnu[1] = -99
            HST_F140W_fnu[0] = 10**( (mag_f140w[0]-8.90) / -2.5 ) * 1e9
        else:
            HST_F140W_fnu[1] = 10**( (mag_f140w[0]-8.90) / -2.5 ) * 1e9
            HST_F140W_fnu[0] = 10**( (mag_f140w[1]-8.90) / -2.5 ) * 1e9
    else:
        HST_F140W_fnu = mag_f140w
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    inputs  = {'redshift':redshift, 'HST_F140W_fnu':HST_F140W_fnu, 'MUV':MUV, 'mStar':mStar}

    if verbose: print(' - Performing selection for objects with provided ranges ')
    goodindices    = np.arange(len(jadesinfo))
    bestmatchinput = None
    for inputkey in inputs.keys():
        if inputs[inputkey] is not None:
            if inputs[inputkey][1] != -99:
                goodi = ju.get_subcat(jadesinfo,inputkey,inputs[inputkey])
                goodindices = np.intersect1d(goodindices,goodi)
            else:
                bestmatchinput = inputkey, inputs[inputkey]

    if (bestmatchinput is not None) & (len(goodindices) > 0):
        if verbose: print(' - Finding best match to "'+bestmatchinput[0]+'" value among the '+str(len(goodindices))+' remaining objects ')
        goodi = ju.get_subcat(jadesinfo[goodindices],bestmatchinput[0],bestmatchinput[1])
        outputinfo  = jadesinfo[goodindices.astype(int)][goodi]
    else:
        outputinfo  = jadesinfo[goodindices.astype(int)]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if len(outputinfo) == 0:
        if verbose:
            print(' - WARNING: No objects were found satisfying the selections:')
            print(inputs)
    else:
        if verbose:
            print(' - Returning the indices for the '+str(len(outputinfo))+' JADES mock objects satisfying the selections:')
            print(inputs)

    return outputinfo
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_subcat(jadesinfo,col,range):
    """
    Return info for a certain set of objects give a seleciton on a single catalog column

    if range[1] == -99 the index for the best match is returned

    """
    if range[1] == -99:
        dval = np.abs(jadesinfo[col] - range[0])
        goodindices = np.where(dval == np.min(dval))
    else:
        goodindices = np.where( (jadesinfo[col] > range[0]) & (jadesinfo[col] < range[1]))

    return goodindices
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =