# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Converting fits files to ascii and ascii files to fits from within python
# For command line compatability see fits2ascii_cmdline.py
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
#import asciitable
import sys
import pyfits
import numpy as np
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def fits2ascii(fitsfile,outpath=None,columns=['all'],verbose=True):
    """
    Convert fits file into ascii (see ascii2fits for the reverse command)

    --- INPUT

    --- EXAMPLE OF USE ---


    NB! not tested as of 160118... download asciitable
    """
    #-------------------------------------------------------------------------------------------------------------
    if verbose: print ' - Loading fits file '+fitsfile
    datfits  = pyfits.open(fitsfile)
    fitstab  = datfits[1].data
    fitscol  = fitstab.columns

    pdb.set_trace()
    #-------------------------------------------------------------------------------------------------------------
    if verbose: print ' - Will write the following columns to the ascii file:'
    if 'all' in columns:
        keys = fitscol.names
        if verbose: print '   all ("all" was found in list of columns)'
    else:
        keys = columns
        if verbose: print '   '+','.join(keys)
    #-------------------------------------------------------------------------------------------------------------
    if verbose: print ' - Initializing and fillling dictionary with data'
    asciidata = {}
    for kk in keys:
        asciidata[kk] = []
        asciidata[kk][:] = fitstab[kk][:]
    #-------------------------------------------------------------------------------------------------------------
    if verbose: print ' - Write dictionary to ascii file:'
    head = fitsfile.split('.fit')[0]
    asciiname = head+'.ascii'
    if outpath != None:
        asciibase = outpath+asciiname.split('/')[-1]
    asciitable.write(asciidata, asciiname, Writer=asciitable.CommentedHeader, names=keys)
    #-------------------------------------------------------------------------------------------------------------
    if verbose: print ' - Wrote data to: ',asciiname

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def ascii2fits(asciifile,asciinames=True,skip_header=0,outpath=None,verbose=True):
    """
    Convert ascii file into fits (see fits2ascii for the reverse command)

    --- INPUT ---
    asciifile        Ascii file to convert
    asciinames       Do the ascii file contain the column names in the header?
    skip_header      The number of header lines to skip when reading the ascii file.
    outpath          Alternative destination for the resulting fits file.

    --- EXAMPLE OF USE ---


    NB! not tested as of 160118... download asciitable
    """
    #-------------------------------------------------------------------------------------------------------------
    if verbose: print ' - Reading ascii file ',asciifile
    data    = np.genfromtxt(asciifile,names=asciinames,skip_header=skip_header,comments='#')
    keys    = data.dtype.names
    #-------------------------------------------------------------------------------------------------------------
    if verbose: print ' - Initialize and fill dictionary with data'
    datadic = {}
    for kk in keys:
        datadic[kk] = []
        datadic[kk] = data[kk]
    if verbose: print ' - found the columns '+','.join(keys)
    #-------------------------------------------------------------------------------------------------------------
    # writing to fits table
    tail = asciifile.split('.')[-1]# reomove extension
    outputfile = asciifile.replace('.'+tail,'.fits')
    if outpath != None:
        outputfile = outpath+outputfile.split('/')[-1]

    columndefs = []
    for key in keys:
        columndefs.append(pyfits.Column(name=key  , format='f', array=datadic[key]))

    cols     = pyfits.ColDefs(columndefs)
    tbhdu    = pyfits.new_table(cols)          # creating table header
    hdu      = pyfits.PrimaryHDU()             # creating primary (minimal) header
    thdulist = pyfits.HDUList([hdu, tbhdu])    # combine primary and table header to hdulist
    thdulist.writeto(outputfile,clobber=True)  # write fits file (clobber=True overwrites excisting file)
    #-------------------------------------------------------------------------------------------------------------
    if verbose: print ' - Wrote the data to: ',outputfile
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
