#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# modifyBoRGcats.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Procedure for modifying a list of borg catalogs.
#----------------------------
#   COMMENTS
#----------------------------
# After running BPZ on the mofidied catalogs modifyBoRGcats_postBPZ.py can be used to
# modify catalogs according to the BPZ results.
#----------------------------
#   INPUTS:
#----------------------------
# catlist          : Path AND name of file containing names of catalogs to modify. These catalogs
#                    should be *_multiband* catalogs containing V, Y, J, and H information.
# outputdir        : Path of directory to create (i.e. cannot already excist) where the
#                    modified catlogs are put.
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --stellaritycut  : The maximum stellarity value alloved
# --sextractorflag : The maximum flag value allowed, i.e., flag <= flagvalue
#                    see section 9.1 here https://www.astromatic.net/pubsvn/software/sextractor/trunk/doc/sextractor.pdf
#                    for description of flags.
# --sn606          : minimum S/N on flux defined as autoflux/autoflux_err for band F606W or F600LP (V) band (Negative value returns SN < sn606)
# --sn098          : minimum S/N on flux defined as autoflux/autoflux_err for band F098W (Y) band (Negative value returns SN < sn098) 
# --sn125          : minimum S/N on flux defined as autoflux/autoflux_err for band F125W (J) band (Negative value returns SN < sn125) 
# --sn160          : minimum S/N on flux defined as autoflux/autoflux_err for band F160W (H) band (Negative value returns SN < sn160) 
# --mag606         : Magnitude range allowed for automag. Expects the two arguments min(mag) max(mag),
#                    for the F606W or F600LP (V) band, i.e., bright limit followed by faint limit
# --mag098         : Magnitude range allowed for automag. Expects the two arguments min(mag) max(mag),
#                    for the F098W (Y) band, i.e., bright limit followed by faint limit
# --mag125         : Magnitude range allowed for automag. Expects the two arguments min(mag) max(mag),
#                    for the F125W (J) band, i.e., bright limit followed by faint limit
# --mag160         : Magnitude range allowed for automag. Expects the two arguments min(mag) max(mag),
#                    for the F160W (H) band, i.e., bright limit followed by faint limit
# --stars          : set this flag to invert the stellarity cut, i.e., search for stars instead of 
#                    galaxies such that objects are required to have stellarity >= stellaritycut instead
#                    of the default stellarity <= stellaritycut
#                    This can be used to create stellar catalogs of the field used for observation planning
# --YmJ            : Cut on Y-J color (i.e. selecting objects to be > YmJ)
# --JmH            : Cut on J-H color (i.e. selecting objects to be < JmH)
# --covcheck       : checking that the object does not fall in a UVIS gap (and therfore has no V detection)
#                    The keyword expects the path and name to a file containg two columns: 1) the rms map
#                    and 2) the segmentation maps.
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
#
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x modifyBoRGcats.py       (only required once)
# bash> modifyBoRGcats.py '/Users/kasperborelloschmidt/work/BoRG/borg_cats/BoRGcatalognames.txt' '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121024' --verbose --stellaritycut 0.8 --sn125 10 --sn160 5 --mag125 1 26 --sextractorflag 3
#
# ---- TESTNG COVCHECK ----
# bash> modifyBoRGcats.py '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/catalognames.txt' '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcatsTESTcovcheck' --verbose --stellaritycut 0.8 --sextractorflag 3 --YmJ 1.75 --JmH 0.5 --sn606 -1.5 --sn160 2.5 --covcheck '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/rmsANDsegmMAPS.txt'
#
# ---- OBSERVATIONS DECEMBER 2012 ----
# bash> modifyBoRGcats.py '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/catalognames.txt' '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/modifiedcatalogs121121' --verbose --stellaritycut 0.8 --sextractorflag 3 --YmJ 1.75 --JmH 0.5 --sn606 -1.5 --sn160 2.5

# bash> modifyBoRGcats.py '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/catalognames.txt' '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/modifiedcatalogs131301_brightobj' --verbose --stellaritycut 0.8 --sextractorflag 3 --mag125 10 23

# bash> modifyBoRGcats.py '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/catalognames.txt' '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/modifiedcatalogs121121stars' --verbose --stellaritycut 0.8 --stars --sextractorflag 3 --sn125 5
#
# ---- CATALOG VERSION 2.1 ----
# bash> modifyBoRGcats.py '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121106/BoRGcatalognames.txt' '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121106/' --verbose --stellaritycut 0.8 --sn125 10 --sn160 5 --sextractorflag 3
#
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-10-24  started by K. B. Schmidt (UCSB)
# 2012-11-21  added --covcheck. K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import subprocess # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
import numpy as np
import pyfits
#----------------------------
#   FUNCTIONS
#----------------------------
def pathAname(str):                         # splitting string with path and name in to
    strsplit = str.split('/')               # splitting string
    name     = strsplit[-1]                 # saving filename (last entry of catsplit)
    slash    = '/'
    path=slash.join(strsplit[0:-1])         # putting path back together
    return [path,name]
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("catlist", type=str, help="Path and name of file containing list of catalog names")
parser.add_argument("outputdir", type=str, help="Path of directory to create and put new catalogs in (cannot exist already)")
# ---- optional arguments ----
parser.add_argument("--stellaritycut", type=float,  help="Maximum value on the object stellarity allowed")
parser.add_argument("--sextractorflag", type=int,  help="Maximum allowed value of the sextractor flags")
parser.add_argument("--flag", type=float, help="allowed SExtractor flag values")
parser.add_argument("--sn606", type=float, help="S/N limit in the V (F606W or F600LP) band")
parser.add_argument("--sn098", type=float, help="S/N limit in the Y (F098W) band")
parser.add_argument("--sn125", type=float, help="S/N limit in the J (F125W) band")
parser.add_argument("--sn160", type=float, help="S/N limit in the H (F160W) band")
parser.add_argument("--mag606", type=float, nargs=2, help="Limits on the F606W (V) magnitude: min(mag) max(mag)")
parser.add_argument("--mag098", type=float, nargs=2, help="Limits on the F098W (Y) magnitude: min(mag) max(mag)")
parser.add_argument("--mag125", type=float, nargs=2, help="Limits on the F125W (J) magnitude: min(mag) max(mag)")
parser.add_argument("--mag160", type=float, nargs=2, help="Limits on the F160W (H) magnitude: min(mag) max(mag)")
parser.add_argument("--YmJ", type=float, help="Lower limit on Y-J color.")
parser.add_argument("--JmH", type=float, help="Upper limit on J-H color.")
parser.add_argument("--covcheck", type=str, help="Path and name of file with rms and segmentation maps for V band")
parser.add_argument("--stars", action="store_true", help="Swap stellaritycut to be stellarity >= stellaritycut")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# reading input file
cats    = np.genfromtxt(args.catlist, dtype=str, comments='#')
Nfields = len(cats)
#-------------------------------------------------------------------------------------------------------------
# preparing output directory
if os.path.exists(args.outputdir): # checking if directory exists
    print ':: '+sys.argv[0]+' :: Selected output directory already exist. Please select another --> ABORTING'
    sys.exit()
else:
    import datetime
    os.mkdir(args.outputdir)

#-------------------------------------------------------------------------------------------------------------
# looping over and manipulating catalogs
for ii in range(Nfields): # looping over the fields found in the input file
    if args.verbose:
        print ':: '+sys.argv[0]+' :: Modifying the catalog:'
        print cats[ii]
    dat      = np.genfromtxt(cats[ii], dtype=None, comments='#')            # reading catalog
    Nobjinit = len(dat['f0'])                                               # The intitial number of objects in catalog
    Nobjmodified = Nobjinit                                                 # initializing Nobjmodified
    PNcats  = pathAname(cats[ii])                                           # splitting catalog in path and name
    if args.stars:
        catname = args.outputdir+'/'+PNcats[1].replace('.cat','_modifiedSTARS.cat')
    else:
        catname = args.outputdir+'/'+PNcats[1].replace('.cat','_modified.cat')

    outcat  = open(catname,"w")
    outcat.write('# Modified version of '+cats[ii]+' \n')
    outcat.write('# For column definitions see header of that file. \n')
    outcat.write('# Created with '+sys.argv[0]+' on '+str(datetime.datetime.now())+' \n')
    outcat.write('# \n')
    outcat.write('# The performed modifications are the following: \n')
    outcat.write('# \n')
    # ---- checking for keywords, manipulating catalog and writing header ----
    # ..................................................................................................................
    # ..................................................................................................................
    if args.stellaritycut:
        stellarity = dat['f12']
        if args.stars:
            goodent    = np.where(stellarity >= args.stellaritycut)
            dat        = dat[goodent]
            Nobjcut    = len(dat['f0'])
            outcat.write('# stellarity     >=   '+str(args.stellaritycut)+' \n')
        else:
            goodent    = np.where(stellarity <= args.stellaritycut)
            dat        = dat[goodent]
            Nobjcut    = len(dat['f0'])
            outcat.write('# stellarity     <=   '+str(args.stellaritycut)+' \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# stellarity     <=   n/a \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')
    # .........................................................
    if args.sextractorflag:
        sextractorflag = dat['f13']
        goodent    = np.where(sextractorflag <= args.sextractorflag)
        dat        = dat[goodent]
        Nobjcut    = len(dat['f0'])

        outcat.write('# SEXtractorflag <=   '+str(args.sextractorflag)+' \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# SEXtractorflag <=   n/a \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')
    # .........................................................
    if args.sn606:
        if args.sn606 > 0:
            #sn606 = dat['f17'] # isomag
            sn606 = dat['f21'] # automag
            goodent    = np.where(sn606 >= args.sn606)
            dat        = dat[goodent]
            Nobjcut    = len(dat['f0'])
            outcat.write('# S/N V-band (automag)     >=   '+str(args.sn606)+' \n')
        else:
            #sn606 = dat['f17'] # isomag
            sn606 = dat['f20'] # automag
            goodent    = np.where(sn606 <= (-1)*args.sn606)
            dat        = dat[goodent]
            Nobjcut    = len(dat['f0'])
            outcat.write('# S/N V-band (automag)     <=   '+str((-1)*args.sn606)+' \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# S/N V-band (automag)     >=   n/a \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')
    # .........................................................
    if args.sn098:
        if args.sn098 > 0:
            #sn098 = dat['f28'] # isomag
            sn098 = dat['f31'] # automag
            goodent    = np.where(sn098 >= args.sn098)
            dat        = dat[goodent]
            Nobjcut    = len(dat['f0'])
            outcat.write('# S/N Y-band (automag)     >=   '+str(args.sn098)+' \n')
        else:
            #sn098 = dat['f28'] # isomag
            sn098 = dat['f31'] # automag
            goodent    = np.where(sn098 <= (-1)*args.sn098)
            dat        = dat[goodent]
            Nobjcut    = len(dat['f0'])
            outcat.write('# S/N Y-band (automag)     <=   '+str((-1)*args.sn098)+' \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# S/N Y-band (automag)      >=   n/a \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')
    # .........................................................
    if args.sn125:
        if args.sn125 > 0:
            #sn125 = dat['f39'] # isomag
            sn125 = dat['f42'] # automag
            goodent    = np.where(sn125 >= args.sn125)
            dat        = dat[goodent]
            Nobjcut    = len(dat['f0'])
            outcat.write('# S/N J-band (automag)     >=   '+str(args.sn125)+' \n')
        else:
            #sn125 = dat['f39'] # isomag
            sn125 = dat['f42'] # automag
            goodent    = np.where(sn125 <= (-1)*args.sn125)
            dat        = dat[goodent]
            Nobjcut    = len(dat['f0'])
            outcat.write('# S/N J-band (automag)     <=   '+str((-1)*args.sn125)+' \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# S/N J-band (automag)      >=   n/a \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')
    # .........................................................
    if args.sn160:
        if args.sn160 > 0:
            #sn160 = dat['f50'] # isomag
            sn160 = dat['f53'] # automag
            goodent    = np.where(sn160 >= args.sn160)
            dat        = dat[goodent]
            Nobjcut    = len(dat['f0'])
            outcat.write('# S/N H-band (automag)     >=   '+str(args.sn160)+' \n')
        else:
            #sn160 = dat['f50'] # isomag
            sn160 = dat['f53'] # automag
            goodent    = np.where(sn160 <= (-1)*args.sn160)
            dat        = dat[goodent]
            Nobjcut    = len(dat['f0'])
            outcat.write('# S/N H-band (automag)     <=   '+str((-1)*args.sn160)+' \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# S/N H-band (automag)     >=   n/a \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')
    # .........................................................
    if args.mag606:
        mag606 = dat['f23']
        goodent    = np.where((mag606 >= args.mag606[0]) & (mag606 <= args.mag606[1]))
        dat        = dat[goodent]
        Nobjcut    = len(dat['f0'])
        outcat.write('# Automag F606W   =   '+str(args.mag606)+' \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# Automag F606W   =   n/a \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')
    # .........................................................
    if args.mag098:
        mag098 = dat['f34']
        goodent    = np.where((mag098 >= args.mag098[0]) & (mag098 <= args.mag098[1]))
        dat        = dat[goodent]
        Nobjcut    = len(dat['f0'])
        outcat.write('# Automag F098W   =   '+str(args.mag098)+' \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# Automag F098W   =   n/a \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')
    # .........................................................
    if args.mag125:
        mag125 = dat['f45']
        goodent    = np.where((mag125 >= args.mag125[0]) & (mag125 <= args.mag125[1]))
        dat        = dat[goodent]
        Nobjcut    = len(dat['f0'])
        outcat.write('# Automag F125W   =   '+str(args.mag125)+' \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# Automag F125W   =   n/a \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')
    # .........................................................
    if args.mag160:
        mag160 = dat['f56']
        goodent    = np.where((mag160 >= args.mag160[0]) & (mag160 <= args.mag160[1]))
        dat        = dat[goodent]
        Nobjcut    = len(dat['f0'])
        outcat.write('# Automag F160W   =   '+str(args.mag160)+' \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# Automag F160W   =   n/a \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')
    # .........................................................
    if args.YmJ:
        #YJcol      = dat['f32']-dat['f43'] # isomag
        YJcol      = dat['f34']-dat['f45'] # automag
        goodent    = np.where(YJcol > args.YmJ)
        dat        = dat[goodent]
        Nobjcut    = len(dat['f0'])
        outcat.write('# Y-J (automag)   >   '+str(args.YmJ)+' \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# Y-J (automag)   =   n/a \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')
    # .........................................................
    if args.JmH:
        #JHcol      = dat['f43']-dat['f54'] # isomag
        JHcol      = dat['f45']-dat['f56'] # automag
        goodent    = np.where(JHcol < args.JmH)
        dat        = dat[goodent]
        Nobjcut    = len(dat['f0'])
        outcat.write('# J-H (automag)   <   '+str(args.JmH)+' \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# J-H (automag)   =   n/a \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')



    # .........................................................
    if args.covcheck:
        RMSandSEGM    = np.genfromtxt(args.covcheck, dtype=None, comments='#')
        Nmaps         = len(RMSandSEGM)
        if Nmaps != Nfields: sys.exit('Number of RMS and SEGMENTATION maps is not the same as Nfields --> ABORTING')

        if args.verbose:
            print ':: '+sys.argv[0]+' :: Reading the RMS and SEGMENTATION maps:'
            print RMSandSEGM[ii,0]
            print RMSandSEGM[ii,1]
        rmsmap     = pyfits.open(RMSandSEGM[ii,0])
        rmsmap     = rmsmap[0].data
        segmap     = pyfits.open(RMSandSEGM[ii,1])
        segmap     = segmap[0].data
        if rmsmap.size != segmap.size: sys.exit('Npixels of the RMS map and the SEGMENTATION map do not agree --> ABORTING')

        ID         = dat['f0']              # IDs for remaining objects
        covvec     = np.array([0]*len(ID))  # nparray for keeping track of objects without coverage

        for kk in range(Nobjmodified): # looping over objects
            IDobj    = ID[kk]
            pixelent = np.where(segmap == IDobj) # pixels for object in segmentation map
            rmsvals  = rmsmap[pixelent]          # the rms values of pixels belogning to object
            if np.isfinite(np.sum(rmsvals)):   # checking for any pixels = inf (would make the sum inf)
                covvec[kk] = ID[kk]
            else:
                covvec[kk] = -99
 
        #print ID[np.where(covvec == -99)]  # printing the IDs of objects falling in UVIS gap

        dat        = dat[np.where(covvec != -99)]
        Nobjcut    = len(dat['f0'])
        outcat.write('# Checked for objects in UVIS gap \n')
        outcat.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat.write('# Did not check for objects in UVIS gap \n')
        outcat.write('# (n/a / n/a passed cut) \n')
        outcat.write('# \n')
    # ..................................................................................................................
    # ..................................................................................................................
    outcat.write('# Hence, the number of objects was reduced from '+str(Nobjinit)+' to '+str(Nobjmodified)+' \n')
    if args.verbose: print '      Number of objects in modified catalog:',Nobjmodified
    outcat.write('# \n')
    for jj in range(len(dat['f0'])): # looping over remaining lines of catalog
        objdata = str(dat[jj])
        objdata = objdata.replace(',','')
        objdata = objdata.replace("'","")
        objdata = objdata.replace('(','')
        objdata = objdata.replace(')','')
        outcat.write(objdata+' \n')
    outcat.close()
    if args.verbose: print ' '
#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

