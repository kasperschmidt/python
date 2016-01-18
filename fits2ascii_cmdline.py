#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# fits2ascii.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Reading fits file and turning it into an ascii file
#----------------------------
#   COMMENTS
#----------------------------
# 
#----------------------------
#   INPUTS:
#----------------------------
# fitsfile         : Fits file to read and turn into ascii
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --columns        : columns to write to ascii file.
#                    DEFAULT is to write all columns to ascii
# --ascii2fits     : The 'inverted' method, i.e. converts ascii file to fits file.
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# asciifile        : Ascii file outputted. Same name (and loaction) as fitsfile but with
#                    .fit* replaced with .ascii
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> fits2ascii.py borg_1510+1115_130426_Y_borg_1510+1115_1705_eps_1Dspec_CALIBRATED.fits --verbose
#
# For multiple files do:
# for i in path/*.fits ; do fits2ascii.py $i --verbose ; done
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-07-08  started by K. B. Schmidt (UCSB)
# 2013-12-09  added ascii2fits keyword. K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import numpy as np         # enable opening with genfromtxt
import pdb                 # for debugging with pdb.set_trace()
import pyfits
import asciitable
import re
import datetime            # print currecnt date and time with datetime.datetime.now()
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("fitsfile", type=str, help="Fitsfile to turn into ascii")
# ---- optional arguments ----
parser.add_argument("--columns", type=str, nargs='+', help="Columns to write to ascii instead of all")
parser.add_argument("--ascii2fits", action="store_true", help="The 'inverted' method, i.e. converts ascii to fits file.")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--show", action="store_true", help="Showing plots on screen")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
if not args.ascii2fits:
    #-------------------------------------------------------------------------------------------------------------
    #                                               FITS 2 ASCII
    #-------------------------------------------------------------------------------------------------------------
    # reading fits file
    datfits  = pyfits.open(args.fitsfile)
    fitstab  = datfits[1].data
    fitscol  = fitstab.columns
    #-------------------------------------------------------------------------------------------------------------
    # the list of keys (columns) to write to ascii
    if not args.columns:
        keys = [fitscol[ii].name for ii in range(len(fitscol[:]))]
    else:
        keys = args.columns
    #-------------------------------------------------------------------------------------------------------------
    # Initialize and fill dictionary with data
    asciidata = {}
    for kk in keys:
        asciidata[kk] = []
        asciidata[kk][:] = fitstab[kk][:]
    #-------------------------------------------------------------------------------------------------------------
    # write dictionary to ascii file
    head = re.compile('\.fi', re.IGNORECASE).split(args.fitsfile)[0]  # making sure case is ignored
    asciiname = head+'.ascii'
    #commentstr = '# File created on ',str(datetime.datetime.now()),' with fits2ascii.py from ',args.fitsfile
    asciitable.write(asciidata, asciiname, Writer=asciitable.CommentedHeader, names=keys)
    #-------------------------------------------------------------------------------------------------------------
    if args.verbose:
        print 'Wrote the columns: ',keys
        print 'From fits file   : ',args.fitsfile
        print 'to the ascii file: ',asciiname
else:
    #-------------------------------------------------------------------------------------------------------------
    #                                             ASCII 2 FITS
    #-------------------------------------------------------------------------------------------------------------
    # reading ascii file
    hdr     = open(args.fitsfile).readlines()[0]
    hdrkeys = hdr.split()[1:]
    data    = np.genfromtxt(args.fitsfile,comments='#')
    #-------------------------------------------------------------------------------------------------------------
    # the list of keys (columns) to write to fits file
    if not args.columns:
        keys = hdr.split()[1:]
    else:
        keys = args.columns
    #-------------------------------------------------------------------------------------------------------------
    # Initialize and fill dictionary with data
    datadic = {}
    for kk in keys:
        datadic[kk] = []
        colent         = np.where(np.asarray(hdrkeys) == kk)[0]
        datadic[kk][:] = data[:,colent]
    #-------------------------------------------------------------------------------------------------------------
    # writing to fits table
    tail = args.fitsfile.split('.')[-1]# reomove extension
    outputfile = args.fitsfile.replace('.'+tail,'.fits')

    columndefs = []
    for key in keys:
        columndefs.append(pyfits.Column(name=key  , format='F', array=datadic[key]))

    cols     = pyfits.ColDefs(columndefs)
    tbhdu    = pyfits.new_table(cols)          # creating table header
    hdu      = pyfits.PrimaryHDU()             # creating primary (minimal) header
    thdulist = pyfits.HDUList([hdu, tbhdu])    # combine primary and table header to hdulist
    thdulist.writeto(outputfile,clobber=True)  # write fits file (clobber=True overwrites excisting file)
    #-------------------------------------------------------------------------------------------------------------
    if args.verbose:
        print 'Wrote the columns: ',keys
        print 'From fits file   : ',args.fitsfile
        print 'to the ascii file: ',outputfile


#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------

