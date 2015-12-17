#!/usr/bin/env python

import sys
import pyfits
from optparse import OptionParser
import cosmics


def cr_clean(fits_fn):
    det = pyfits.getval(fits_fn, 'DETECTOR')

    if det == 'IR':
        print('CR-cleaning IR FLT file: %s' % fits_fn)

        readnoise = 20.0     # for WFC3/IR
        #sigclip = 5.0
        sigclip = 5.0
        sigfrac = 0.3     # small means crs are grown more  0.3 is good!
        #objlim = 5.0
        objlim = 12.0
        #maxiter = 4
        maxiter = 3
        exptime = pyfits.getval(fits_fn, 'EXPTIME')
        gain = exptime

        (image, header) = cosmics.fromfits(fits_fn)
        c = cosmics.cosmicsimage(image, gain=gain, readnoise=readnoise, sigclip=sigclip, sigfrac=sigfrac, objlim=objlim)
        c.run(maxiter=maxiter)

        clean_fn = fits_fn.replace('.fits', '_crclean.fits')
        mask_fn = fits_fn.replace('.fits', '_crmask.fits')

        fo = pyfits.open(fits_fn)
        fo[1].data = c.cleanarray
        fo.writeto(clean_fn)
        fo.close()

        #cosmics.tofits(clean_fn, c.cleanarray, header)
        cosmics.tofits(mask_fn, c.mask, header)

    if det == 'UVIS':
        print('CR-cleaning UVIS FLT file: %s' % fits_fn)
        readnoise = 3.105     # for WFC3/UVIS

        # good
        #sigclip = 5.0
        #sigfrac = 0.1     # small means crs are grown more  0.3 is good!
        #objlim = 5.0

        # some residuals
        #sigclip = 4.5
        #sigfrac = 0.1     # small means crs are grown more  0.3 is good!
        #objlim = 4.0

        sigclip = 4.5
        sigfrac = 0.1     # small means crs are grown more  0.3 is good!
        objlim = 2.0

        maxiter = 4
        exptime = pyfits.getval(fits_fn, 'EXPTIME')
        gain = 1.0
        ff = pyfits.open(fits_fn)
        header = ff[0].header
        image = ff[1].data

        #(image, header) = cosmics.fromfits(fits_fn)
        c = cosmics.cosmicsimage(image, gain=gain, readnoise=readnoise, sigclip=sigclip, sigfrac=sigfrac, objlim=objlim)
        c.run(maxiter=maxiter)

        image2 = ff[4].data
        d = cosmics.cosmicsimage(image2, gain=gain, readnoise=readnoise, sigclip=sigclip, sigfrac=sigfrac, objlim=objlim)
        d.run(maxiter=maxiter)

        clean_fn = fits_fn.replace('.fits', '_crclean.fits')
        mask_fn = fits_fn.replace('.fits', '_crmask.fits')

        fo = pyfits.open(fits_fn)
        fo[1].data = c.cleanarray
        fo[4].data = d.cleanarray
        fo.writeto(clean_fn)
        fo.close()

        #cosmics.tofits(clean_fn, c.cleanarray, header)
        cosmics.tofits(mask_fn, c.mask, header)


def main():
    usage = '%prog [options] <filelist>'
    version = '%prog 0.1'
    parser = OptionParser(usage=usage, version=version)
    parser.add_option('-e', '--exten', type='int', default=0, help='FITS extension number')
    parser.add_option('-o', action='store_true', default=False, help='Save fits header to a file')
    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
        sys.exit()
    for fits_fn in args:
        cr_clean(fits_fn)

if __name__ == '__main__':
    main()


