#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# modifyBoRGcats_postBPZ.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Procedure similar to modifyBoRGcats.py but here modifying BoRG catalogs after the
# BPZ photometric redshifts have been obtained.
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# dirlist          : Path AND name of file containing names of directories containing BPZ output.
#                    If the results have been created with runBPZmultiplecats.py the directories
#                    contain names including a time stamp yymmddhhmmss, e.g., 121019120545_
#                    NOTE: same input as web4BPZgifs.py and plotBPZhistcounts.py use.
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --zb             : The allowed range in redshift, i.e., minz maxz
# --zberr          : The maximum allowed error on the redshift. The same value is applied to zerr+ and zerr-
# --odds           : The minimum redshift odds allowed
# --chisq2         : the largest reduced chi^2 values accepted 
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# new catalogs     : Catalogs after applying specified cuts. Named as the found catalogs (*BPZinput_bpz.cat 
#                    and *BPZinput.cat) but added the extension 'postBPZmodified'
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x modifyBoRGcats_postBPZ.py       (only required once)
# bash> modifyBoRGcats_postBPZ.py '/Users/kasperborelloschmidt/work/BoRG/modifiedBoRGcats121025BPZ/resultdirectorylist121025run.txt' --verbose --zb 0.4 3.5 --zberr 3.0 --odds 0.5 --chisq2 1.0
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-10-29  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import subprocess # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
import numpy as np
import commands
import datetime
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
parser.add_argument("dirlist", type=str, help="Path and name of file containing list of directories")
# ---- optional arguments ----
parser.add_argument("--zb",type=float,nargs=2, help="Range of photo z")
parser.add_argument("--zberr",type=float,help="Largest allowed error on phot z (same cut on zerr+ and zerr-)")
parser.add_argument("--odds",type=float,help="The smallest odds value accepted")
parser.add_argument("--chisq2",type=float,help="The largest chisq2 (reduce chi squared) value accepted")
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
dirALL  = np.genfromtxt(args.dirlist, dtype=str, comments='#')
Ndirs    = len(dirALL)
#-------------------------------------------------------------------------------------------------------------
# Finding BPZ otput file and reading it

for ii in range(Ndirs):                                             # loop over directories
    dir       = dirALL[ii]

    # --- looking for BPZ input catalog ---
    borgpath  = dir+'/*_BPZinput.cat'                               # search string
    borgcat   = commands.getoutput('ls '+borgpath)
    if os.path.isfile(borgcat):                                     # if file excists read it
        catBORG     = np.genfromtxt(borgcat, dtype=None, comments='#')
        PNborgcat = pathAname(borgcat)
        catname     = PNborgcat[-1].split('_multiband')             # the BoRG field name
    else:
        print ':: '+sys.argv[0]+' :: ERROR :: The '+borgcat+' was not found in '+dir[0]+'  --> ABPRTING'
        sys.exit()

    # --- looking for BPZ output catalog ---
    bpzpath   = dir+'/*BPZinput_bpz.cat'                                    # search string
    bpzoutput = commands.getoutput('ls '+bpzpath)
    if os.path.isfile(bpzoutput):                                   # if file excists read it
        catBPZ      = np.genfromtxt(bpzoutput, dtype=None, comments='#')
        PNbpzoutput = pathAname(bpzoutput)
    else:
        print ':: '+sys.argv[0]+' :: ERROR :: The '+bpzoutput+' was not found in '+dir[0]+'  --> ABPRTING'
        sys.exit()


    if len(catBPZ['f0']) != len(catBORG['f0']): 
        print ':: '+sys.argv[0]+' :: ERROR :: The two read catalogs:'
        print ':: '+sys.argv[0]+'             '+borgcat
        print ':: '+sys.argv[0]+'             '+bpzoutput
        print ':: '+sys.argv[0]+'             do not have the same number of objects --> ABORTING'
        sys.exit()
    #-------------------------------------------------------------------------------------------------------------
    Nobjinit = len(catBPZ['f0'])                                    # The intitial number of objects in catalog
    Nobjmodified = Nobjinit                                         # initializing Nobjmodified

    catname = borgcat.replace('.cat','_postBPZmodified.cat')
    outcat1 = open(catname,"w")

    outcat1.write('# Modified version of '+borgcat+' \n')
    outcat1.write('# For column definitions see header of that file. \n')
    outcat1.write('# Created with '+sys.argv[0]+' on '+str(datetime.datetime.now())+' \n')
    outcat1.write('# \n')
    outcat1.write('# The performed modifications are the following: \n')
    outcat1.write('# \n')

    catname = bpzoutput.replace('_bpz.cat','_postBPZmodified_bpz.cat')
    outcat2 = open(catname,"w")

    outcat2.write('# Modified version of '+bpzoutput+' \n')
    outcat2.write('# For column definitions see header of that file. \n')
    outcat2.write('# Created with '+sys.argv[0]+' on '+str(datetime.datetime.now())+' \n')
    outcat2.write('# \n')
    outcat2.write('# The performed modifications are the following: \n')
    outcat2.write('# \n')
    # ---- checking for keywords, manipulating catalog and writing header ----
    # ..................................................................................................................
    # ..................................................................................................................
    if args.odds:
        odds       = catBPZ['f5']
        goodent    = np.where(odds >= args.odds)
        catBPZ     = catBPZ[goodent]
        Nobjcut    = len(catBPZ['f0'])

        outcat1.write('# zb odds        >=   '+str(args.odds)+' \n')
        outcat1.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat1.write('# \n')

        outcat2.write('# zb odds        >=   '+str(args.odds)+' \n')
        outcat2.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat2.write('# \n')

        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat1.write('# zb odds        >=   n/a \n')
        outcat1.write('# (n/a / n/a passed cut) \n')
        outcat1.write('# \n')

        outcat2.write('# zb odds        >=   n/a \n')
        outcat2.write('# (n/a / n/a passed cut) \n')
        outcat2.write('# \n')
    # .........................................................
    if args.zb:
        zb = catBPZ['f1']
        goodent    = np.where((zb >= args.zb[0]) & (zb <= args.zb[1]))
        catBPZ        = catBPZ[goodent]
        Nobjcut    = len(catBPZ['f0'])

        outcat1.write('# zb (photo z)    =   '+str(args.zb)+' \n')
        outcat1.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat1.write('# \n')

        outcat2.write('# zb (photo z)    =   '+str(args.zb)+' \n')
        outcat2.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat2.write('# \n')

        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat1.write('# zb (photo z)    =   n/a \n')
        outcat1.write('# (n/a / n/a passed cut) \n')
        outcat1.write('# \n')

        outcat2.write('# zb (photo z)    =   n/a \n')
        outcat2.write('# (n/a / n/a passed cut) \n')
        outcat2.write('# \n')
    # .........................................................
    if args.zberr:
        zberrmin   = catBPZ['f2']
        zberrmax   = catBPZ['f3']
        goodent    = np.where((zberrmin <= args.zberr) & (zberrmax <= args.zberr))
        catBPZ     = catBPZ[goodent]
        Nobjcut    = len(catBPZ['f0'])

        outcat1.write('# zberr          <=   '+str(args.zberr)+' \n')
        outcat1.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat1.write('# \n')

        outcat2.write('# zberr          <=   '+str(args.zberr)+' \n')
        outcat2.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat2.write('# \n')

        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat1.write('# zberr          <=   n/a \n')
        outcat1.write('# (n/a / n/a passed cut) \n')
        outcat1.write('# \n')

        outcat2.write('# zberr          <=   n/a \n')
        outcat2.write('# (n/a / n/a passed cut) \n')
        outcat2.write('# \n')
    # .........................................................
    if args.chisq2:
        chisq2     = catBPZ['f10']
        goodent    = np.where(chisq2 <= args.chisq2)
        catBPZ     = catBPZ[goodent]
        Nobjcut    = len(catBPZ['f0'])

        outcat1.write('# zb chisq2      >=   '+str(args.chisq2)+' \n')
        outcat1.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat1.write('# \n')

        outcat2.write('# zb chisq2      >=   '+str(args.chisq2)+' \n')
        outcat2.write('# ('+str(Nobjcut)+' / '+str(Nobjmodified)+' passed cut) \n')
        outcat2.write('# \n')
        Nobjmodified = Nobjcut    # update Nobjmodified
    else:
        outcat1.write('# zb chisq2      >=   n/a \n')
        outcat1.write('# (n/a / n/a passed cut) \n')
        outcat1.write('# \n')

        outcat2.write('# zb chisq2      >=   n/a \n')
        outcat2.write('# (n/a / n/a passed cut) \n')
        outcat2.write('# \n')
    # ..................................................................................................................
    # ..................................................................................................................
    outcat1.write('# Hence, the number of objects was reduced from '+str(Nobjinit)+' to '+str(Nobjmodified)+' \n')
    outcat1.write('# \n')

    outcat2.write('# Hence, the number of objects was reduced from '+str(Nobjinit)+' to '+str(Nobjmodified)+' \n')
    outcat2.write('# \n')

    for jj in range(len(catBPZ['f0'])): # looping over remaining lines of catalog
        IDbpz   = catBPZ['f0']
        goodent = np.where(catBORG['f0'] == IDbpz[jj]) # matching on ID
        if len(goodent) == 1:
            objdata1 = str(catBORG[goodent])
            objdata1 = objdata1.replace(',','')
            objdata1 = objdata1.replace("'","")
            objdata1 = objdata1.replace('(','')
            objdata1 = objdata1.replace(')','')
            objdata1 = objdata1.replace('[','')
            objdata1 = objdata1.replace(']','')
            outcat1.write(objdata1+' \n')
        else:
           print ':: '+sys.argv[0]+' :: ERROR :: Object with ID '+str(IDbpz)+' missing in '+borgcat+'  --> ABPRTING'
           sys.exit()

        objdata2 = str(catBPZ[jj])
        objdata2 = objdata2.replace(',','')
        objdata2 = objdata2.replace("'","")
        objdata2 = objdata2.replace('(','')
        objdata2 = objdata2.replace(')','')
        outcat2.write(objdata2+' \n')

    outcat1.close()
    outcat2.close()
#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

