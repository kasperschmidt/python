#!/usr/bin/env python

import sys
import os
import os.path as op
import glob
import numpy as np
import drizzlepac as dp
import pyfits

def classify_flts(flt_fns):
    fndict = {}
    for fn in flt_fns:
        filtname = pyfits.getheader(fn)['FILTER'].lower()
        fndict.setdefault(filtname, []).append(fn)
    return fndict

def deg2sex2(deg, hours=False, verbose=False):
    sgn = 1
    if deg < 0:
        sgn = -1
        deg *= sgn
    if hours:
        dd = deg / 15.0
        fmt_d = '%02d'
        d_or_h = 'h'
        dms_hms = '(hrs, min, sec)'
    else:
        dd = deg
        fmt_d = '%+03d'
        d_or_h = 'd'
        dms_hms = '(deg, min, sec)'
    fmt_m = '%02d'
    fmt_s = '%10.8f'
    d = int(dd)
    min = (dd - d) * 60.0
    min = round(min)
    dms1 = ('%s%s' % (fmt_d, fmt_m)) % (sgn * d, min)
    return dms1

def borgfieldname(fn):
    ff = pyfits.open(fn)
    ra = ff[0].header['RA_TARG']        # exactly the same as CRVAL1
    dec = ff[0].header['DEC_TARG']      # exactly the same as CRVAL2
    ra_sex = deg2sex2(ra, hours=True)
    dec_sex = deg2sex2(dec)
    fieldname = '%s_%s%s' % ('borg', ra_sex, dec_sex)
    return fieldname

def align_flt(filtname, flt_fns, conv_width=None, threshold=10.0, minobj=10, searchrad=2.0, final=False, wcsname='fltaln'):
    if not conv_width:
        if filtname in ['f600lp', 'f606w']:
            conv_width = 3.5   # UVIS, ACS
        else:
            conv_width = 2.5   # IR
    peakmax = 50000      # to exclude saturated stars
    if final:
        see2dplot = False
        residplot = 'None'
        updatehdr = True
    else:
        see2dplot = True
        residplot = 'both'
        updatehdr = False
    fns_str = ','.join(flt_fns)
    dp.tweakreg.TweakReg(fns_str, updatehdr=updatehdr, conv_width=conv_width, peakmax=peakmax, minobj=minobj, threshold=threshold, searchrad=searchrad, see2dplot=see2dplot, residplot=residplot, wcsname=wcsname)

def align_drz(drz_fns, refimage, conv_width=3.5, threshold=5.0, minobj=10, searchrad=2.0, final=True):
    peakmax = 50000    # exclude saturated stars
    if final:
        see2dplot = False
        residplot = 'None'
        updatehdr = True
    else:
        see2dplot = True
        residplot = 'both'
        updatehdr = False
    fns_str = ','.join(drz_fns)
    dp.tweakreg.TweakReg(fns_str, refimage=refimage, updatehdr=updatehdr, conv_width=conv_width, peakmax=peakmax, minobj=minobj, threshold=threshold, searchrad=searchrad, see2dplot=see2dplot, residplot=residplot)

def drizzle_flt(filtname, flt_fns, outname, refimage=None):
    """ drizzle flt files """
    final_wcs = True       # required to set final_scale, final_rot
    final_scale = 0.08
    final_pixfrac = 0.75
    final_rot = 0.0
    if filtname in ['f600lp', 'f606w']:
        #bits = '64,32'   # allow warm pixels and dark CTE tails
        bits = 64        # allow warm pixels
    else:
        bits = 512     # allow IR blobs
    wht_type = 'IVM'
    num_cores = 1
    dp.astrodrizzle.AstroDrizzle(flt_fns, output=outname, num_cores=num_cores, driz_sep_bits=bits, driz_cr_corr=True, final_bits=bits, final_wht_type=wht_type, final_wcs=final_wcs, final_scale=final_scale, final_pixfrac=final_pixfrac, final_rot=final_rot, final_refimage=refimage, final_fillval=0)

