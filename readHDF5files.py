#!/opt/local/bin/python2.7
#+
#----------------------------
#   NAME
#----------------------------
# readHDF5files.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Reading HDF5 file (*.hdf5) and writing the content to a python dictionary (*.npz) 
# if 'writedict' keyword is set.
#
# Created to read BoRG selection and completeness functions from Pascal Oesch
#----------------------------
#   COMMENTS
#----------------------------
# Note that the python version #!/opt/local/bin/python2.7 (installed via macports)
# is used instead of the default #!/usr/bin/env python2.7 (containing pip installations)
# as this is the architecture I managed to buld HDF5 and h5py on (KBS: 130628)
#
# To execute on mutliple files simply do something along the lines:
# bash> for i in ~/path/to/files*/*.hdf5 ; do readHDF5files.py $i --writedict --verbose ; done
#----------------------------
#   INPUTS:
#----------------------------
# hdf5file         : HDF5 file containing data to read
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --writedict      : write dictionary with HDF5 file content.
#                    The dictionary will be put in the same directory as the HDF5 file
#                    and will have the extention .npz
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
#
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> readHDF5files.py ~/work/BoRG/borgdata/pascalsim/DELIVERY_v1/borg_1437+5043/Szm_borg_1437+5043_SN5.hdf5 --writedict --verbose
# bash> readHDF5files.py ~/work/BoRG/borgdata/pascalsim/DELIVERY_v1/borg_1437+5043/C_borg_1437+5043_SN5.hdf5 --writedict --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-06-28  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import getopt              # used to extract/obtain the optional input
import numpy as np
import h5py
import pdb
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("hdf5file", type=str, help="HDF5 file to read")
# ---- optional arguments ----
parser.add_argument("--writedict", action="store_true", help="Setting this keyword will write output to dictionary")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print 'Reading',args.hdf5file
f     = h5py.File(args.hdf5file, 'r')
keys  = f.keys()
if args.verbose: print 'Found the keys: ',keys
Nkeys = len(keys)
#-------------------------------------------------------------------------------------------------------------
if args.writedict:
    # Building dictionary
    dict = {}
    for ii in range(Nkeys):
        dict[str(keys[ii])] = np.asarray(f[keys[ii]])

    outfile = args.hdf5file.replace('.hdf5','.npz')
    np.savez(outfile,**dict) # save array as binary file
    if args.verbose: print '\nSaved output dictionary to: \n',outfile
#-------------------------------------------------------------------------------------------------------------
f.close() # closing HDF5 file
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

