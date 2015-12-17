#!/usr/bin/env python
# 2012.02.24  L. Bradley
# 2013.02.05  Using IDL version of dust_eval to get E(B-V) value. K. B. Schmidt (UCSB)
# 2013.04.17  Adding F105W filter info. K. B. Schmidt (UCSB)

import sys
import os
import glob
import pyfits
import numpy as np
from optparse import OptionParser
import dust_getval as dust
from astropysics.coords import ICRSCoordinates, GalacticCoordinates
import pdb # for debugging

def get_fsig(filtname, uvis='borg'):
    depth_fn = 'sigma_empty.dat'
    if os.path.isfile(depth_fn):
        fm = open(depth_fn)
        sigs = fm.readlines()
        fm.close()
        filtname = filtname.upper()
    else:
        return None
    sig = {}
    if uvis == 'borg':
        sig['F606W'] = sigs[0]
    if uvis == 'hippies':
        sig['F600LP'] = sigs[0]
    sig['F098M'] = sigs[1]
    sig['F125W'] = sigs[2]
    sig['F160W'] = sigs[3]

    if filtname not in sig:
        return None
    return float(sig[filtname])

def field_info(fits_fn, uvis='borg'):
    if '.fits' not in fits_fn:
        print 'Skipping %s....not a fits file' % fits_fn
        return
    field = fits_fn[0:14]
    av_ebv = {}
    av_ebv['F300X'] = 6.78362003559
    av_ebv['F475X'] = 3.79441819047
    av_ebv['F475W'] = 3.82839055809
    av_ebv['F606W'] = 3.01882984135
    av_ebv['F600LP'] = 2.24159324026
    av_ebv['F098M'] = 1.29502816006
    av_ebv['F105W'] = 1.18148250758
    av_ebv['F125W'] = 0.893036743585
    av_ebv['F160W'] = 0.633710427959

    hstzpt = {}
    hstzpt['F300X'] = 24.9651978604
    hstzpt['F475X'] = 26.1590189073
    hstzpt['F475W'] = 25.6883822362
    hstzpt['F606W'] = 26.0814041237
    hstzpt['F600LP'] = 25.8814691531
    hstzpt['F098M'] = 25.68102452
    hstzpt['F105W'] = 26.27
    hstzpt['F125W'] = 26.2473632068
    hstzpt['F160W'] = 25.9558992372
    fobj = pyfits.open(fits_fn)
    hdr = fobj[0].header
    filtname = hdr['FILTER']
    ra = hdr['CRVAL1']
    dec = hdr['CRVAL2']
    hst_zpt = hstzpt[filtname]
    gcoords = ICRSCoordinates(ra, dec).convert(GalacticCoordinates)
    gall = gcoords.l.degrees
    galb = gcoords.b.degrees
    ebv = dust.get_dust(gall, galb, interp=1)[0]
    ebv = ebv[0]
    # --- KBS Using idl to get ebv value ---
    import commands
    cmdstring = 'idl -e '+'"print,dust_getval('+str(gall)+','+str(galb)+",/interp,ipath='~/idl711mac/itt/idl71/lib/schlegelmapsIDL/SFD_4096/')xxx".replace('xxx','"')
    ebv = float(commands.getoutput(cmdstring).split('\n')[-1])  # overwriting ebv value
    #import pdb                 # KBS for debugging
    #pdb.set_trace() # KBS
    # --------------------------------------
    if filtname in av_ebv:
        av = av_ebv[filtname] * ebv
    else:
        print 'filter %s not found'
        sys.exit()
    hstzpt_ext = hst_zpt - av
    exptime = hdr['EXPTIME']
    area_fn = 'search_area.dat'
    if os.path.isfile(area_fn):
        fa = open(area_fn)
        npix = fa.readlines()[0]
        scale = 0.08    # arcsec/pixel
        area = np.int(npix) * scale**2 / 3600.0      # sq arcmin
        fa.close()
    else:
        area = -99

    fsig = get_fsig(filtname, uvis=uvis)
    if fsig:
        mlim = hstzpt_ext - (2.5 * np.log10(fsig * 5.0))     # 5 sigma lim
    else:
        mlim = -99.0
    return field, filtname, exptime, area, mlim, ra, dec, gall, galb, ebv, hst_zpt, av, hstzpt_ext


def main():
    usage = '%prog [options] <filelist>'
    version = '%prog 0.1'
    parser = OptionParser(usage=usage, version=version)
    (options, args) = parser.parse_args()

    fn_out = 'zzzzzz_borg_field_info.txt'
    fo = open(fn_out, mode='a')
    str = """# 1 Field
# 2 Filter
# 3 Exptime (s)
# 4 F125W Search Area (square arcmin)
# 5 5-sigma Limiting Magnitude (including Galactic extinction)
# 6 RA
# 7 Dec
# 8 Galactic l
# 9 Galactic b
# 10 E(B-V) from Schlegel et al. (1998)
# 11 Nominal HST zero-point
# 12 A_V
# 13 Corrected HST zero-point
"""
    fo.write(str)
    dirs = glob.glob('borg_*')
    for dir in dirs:
        if os.path.isdir(dir):
            print 'Processing %s' % dir
            os.chdir(dir)
            fits_fns = glob.glob('*_drz.fits')
            for fits_fn in fits_fns:
                if 'f600lp' in fits_fn:
                    uvis = 'hippies'
                if 'f606w' in fits_fn:
                    uvis = 'borg'
            for fits_fn in fits_fns:
                (field, filtname, exptime, area, mlim, ra, dec, gall, galb, ebv, hst_zpt, av, hstzpt_ext) = field_info(fits_fn, uvis=uvis)
                line = '%s %s %s %s %s %s %s %s %s %s %s %s %s' % (field, filtname, exptime, area, mlim, ra, dec, gall, galb, ebv, hst_zpt, av, hstzpt_ext)
                fo.write('%s\n' % line)
            fo.write('\n')
            os.chdir('../')
    print 'Writing %s' % fn_out
    fo.close()

if __name__ == '__main__':
    main()

