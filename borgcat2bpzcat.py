#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# borgcat2bpzcat.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Routine turning (multiple) standard BoRG catalogs into (multiple) 
# calags readable by the BPZ code (*.cat). For each *.cat file a *.columns
# file also needed by BPZ is created as well.
#----------------------------
#   COMMENTS
#----------------------------
#
# NB! The content of the *.columns file is hardcoaded and the same 
# for all catalogs - i.e. it is tailored at BoRG catalogs only. To change
# the content edit the code below.
#
#----------------------------
#   INPUTS:
#----------------------------
# borgcatlist     : string containing name and path to a text file containing
#                   names and paths of the borg catalogs to transform
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# -nondetection   : if this keyword is set the magnitudes and errors of non-detections (S/N < 1) 
#                   are set to 99 and the limiting magnitude respectively as suggested on
#                   http://acs.pha.jhu.edu/%7Etxitxo/bpzdoc.html
# -vband*         : Use this keyword to indicate which V-filter the fields have. Set either:
#                       -vbandF606W   (DEFAULT)
#                       -vbandF600LP
# -mag*           : Use this keyword to indicate what magnitudes to use. Set either:
#                       -magauto      (DEFAULT)
#                       -magiso
# -F105W          : Set this keyword if the input catalogs contain F105W data
#                   Assumes it's added as the last columns (58-69) to the catalgs
#                   and includes the nescessary lines in output
# -fieldinfo      : provide thepath an name to the borg infofile for fields. Defaut is 
#                   /Users/kasperborelloschmidt/work/BoRG/BoRG_field_info_all130313.txt 
# -EPS            : set -EPS to create an .eps version of plots
# -VERBOSE        : set -VERBOSE to get info/messages printed to the screen
#----------------------------
#   OUTPUTS:
#----------------------------
# outputfile      : New Fits file containing SDSS objects, same format
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x borgcat2bpzcat.py       (only required once)
# bash> ./borgcat2bpzcat.py '/Users/kasperborelloschmidt/work/BoRG/borg_catsBPZ/BoRGcatalogs.txt' -VERBOSE

# bash> ./borgcat2bpzcat.py '/Users/kasperborelloschmidt/work/BoRG/observationDec2012/BPZrun121217/catalognames_modified.txt' -VERBOSE
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-10-05  started by K. B. Schmidt (UCSB)
# 2013-01-25  added nondetection keyword. K. B. Schmidt (UCSB)
# 2013-03-13  added -vband* and -mag* keywords. K. B. Schmidt (UCSB)
# 2013-04-17  added F105W and -fieldinfo keywords. K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import pyfits
import numpy as np
import sys    # enabling arguments to code
import pdb                 # for debugging
#----------------------------
#   FUNCTIONS
#----------------------------
def test4float(str):
    try:
        float(str)
        return 1
    except ValueError:
        return 0
def pathAname(str):
    strsplit = str.split('/')                  # splitting string
    name     = strsplit[-1]                    # saving filename (last entry of catsplit)
    slash    = '/'
    path=slash.join(strsplit[0:-1])+'/'        # putting path back together
    return [path,name]
def mlim5sigTO1sig(mlim5sig):                  # converting 5 sigma limiting magnitude to 1 sigma limiting magnitude
    mlim1sig = -2.5*np.log10(1./5.)+mlim5sig   # calculating 1 sigma limit from 5 sigma limit
    return mlim1sig
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# optional input
VB = 0
if '-VERBOSE' in sys.argv: VB = 1
PS = 0
if '-EPS' in sys.argv: PS = 1
ND = 0
if '-nondetection' in sys.argv: ND = 1

vband = 'F606W' # default V band set to F606W
if '-vbandF600LP' in sys.argv: vband = 'F600LP'
if '-vbandF606LP' in sys.argv: vband = 'F606W'

maguse = 'auto' # default magnitudes to use
if '-magauto' in sys.argv: maguse = 'auto'
if '-magiso'  in sys.argv: maguse = 'iso'
#-------------------------------------------------------------------------------------------------------------
if VB == 1:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
Narg = 1
if len(sys.argv) < Narg: # checking that at lest Narg arguments are provided
    errorSTR = "ERROR: Must provide at least "+str(Narg-1)+" arguments to "+sys.argv[0]+" (besides optional arguments) "
    sys.exit(errorSTR)   # stop the program and print an error message
#-------------------------------------------------------------------------------------------------------------
# loading info about fields
infofile = '/Users/kasperborelloschmidt/work/BoRG/BoRG_field_info_all130313.txt'
if '-fieldinfo' in sys.argv: 
    ent = np.where(np.asarray(sys.argv) == '-fieldinfo')
    infofile = sys.argv[ent[0]+1]

if VB == 1: print ':: '+sys.argv[0]+' :: Using the field info from ',infofile
fldinfo  = np.genfromtxt(infofile,dtype=None, comments='#')
#-------------------------------------------------------------------------------------------------------------
# reading input file
buffer    = open(sys.argv[1]).readlines()   # read file into buffer
pAnINPUT  = pathAname(sys.argv[1])          # splitting input file in name and path
Nlines    = len(buffer)                     # total number of lines in file
catlist   = []                              # defining list to put catalog names in
Ncatalogs = 0                               # resetting counter
for line in buffer:
    if Ncatalogs>=Nlines : break            # if all lines read, break out of loop
    if line[0]=='#': continue               # if comment continue to next line
    Ncatalogs=Ncatalogs+1
    catlist.append(line.rstrip())           # appending catalog to list; rstrip removes the \n for new lines
if VB == 1: print ':: '+sys.argv[0]+' :: Found',Ncatalogs,'catalogs in input file'
#-------------------------------------------------------------------------------------------------------------
for cat in catlist:                         # looping over catalogs an transforming them into BPZ cat files
    pAn = pathAname(cat)                    # getting path and string
    catnamesplit = pAn[1].split('.') 
    if len(catnamesplit)>2:
        errorSTR = "ERROR: "+sys.argv[0]+" The filename "+catname+" contains multiple dots (.)"
        sys.exit(errorSTR)                  # stop the program and print an error message
    bpzcat  = pAnINPUT[0]+catnamesplit[0]+'_BPZinput.cat'      # putting to gether name of BPZ catalog
    bpzcol  = pAnINPUT[0]+catnamesplit[0]+'_BPZinput.columns'  # putting to gether name of BPZ columns file
    if VB == 1:
        print ':: '+sys.argv[0]+' :: Will write output to ',bpzcat
        print '                      and                  ',bpzcol

    # --- Editing individual lines and writing them to new catalog ---
    bufferIN = open(cat).readlines()                                 # read file into buffer
    fileOUT  = open(bpzcat,"w")                                      # open new file (overwrites or creates)
    Nobj     = 0                                                     # resetting counter
    for line in bufferIN:                                            # looping over lines
        if line[0]=='#': 
            fileOUT.write(line)                                      # writing header lines to BPZ catalog
            continue
        Nobj = Nobj+1
        if Nobj==1: fileOUT.write('# \n')                            # making sure the header ends with an empty line starting with #
        linelist = line.split()                                      # turning line into list of values
        newline  = ' '
        for ii in range(len(linelist)):                              # looping over columns and getting rid of strings
            col = linelist[ii]
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            SNnondetection = 1.0                                     # the limit between detection and non-detection
            if ND == 1:                                              # check if non-detection. Marking as m=99 and dm=mlim 
                field    = linelist[1][0:14]
                if ii in [17,20]: # Vband
                    mliment  = np.where((fldinfo['f0'] == field) & (fldinfo['f1'] == vband)) 
                    mlim5sig = fldinfo['f4'][mliment]                # hstzpt_ext-2.5*np.log10(fsig*5.0) from borg_fldinfo.py
                    mlim1sig = mlim5sigTO1sig(mlim5sig)              # getting 1 sigma limiting magnitude

                if ii in [28,31]: # Yband
                    yband = 'F098M'
                    mliment  = np.where((fldinfo['f0'] == field) & (fldinfo['f1'] == yband)) 
                    mlim5sig = fldinfo['f4'][mliment]                # hstzpt_ext-2.5*np.log10(fsig*5.0) from borg_fldinfo.py
                    mlim1sig = mlim5sigTO1sig(mlim5sig)              # getting 1 sigma limiting magnitude

                if ii in [39,42]: # Jband
                    jband = 'F125W'
                    mliment  = np.where((fldinfo['f0'] == field) & (fldinfo['f1'] == jband)) 
                    mlim5sig = fldinfo['f4'][mliment]                # hstzpt_ext-2.5*np.log10(fsig*5.0) from borg_fldinfo.py
                    mlim1sig = mlim5sigTO1sig(mlim5sig)              # getting 1 sigma limiting magnitude

                if ii in [50,53]: # Hband
                    hband = 'F160W'
                    mliment  = np.where((fldinfo['f0'] == field) & (fldinfo['f1'] == hband)) 
                    mlim5sig = fldinfo['f4'][mliment]                # hstzpt_ext-2.5*np.log10(fsig*5.0) from borg_fldinfo.py
                    mlim1sig = mlim5sigTO1sig(mlim5sig)              # getting 1 sigma limiting magnitude

                if '-F105W' in sys.argv:
                    if ii in [61,64]: # F105W
                        F105band = 'F105W'
                        mliment  = np.where((fldinfo['f0'] == field) & (fldinfo['f1'] == F105band)) 
                        mlim5sig = fldinfo['f4'][mliment]
                        mlim1sig = mlim5sigTO1sig(mlim5sig)

                isosncol = [17,28,39,50]                             # list of S/N_iso columns (col number-1)
                if '-F105W' in sys.argv: isosncol = [17,28,39,50,61]
                if ii in isosncol and float(col) < SNnondetection:   # check if S/N is high enough
                    linelist[ii+4] = str(99)                         # setting iso magnitude of non-detection to 99
                    linelist[ii+5] = str(mlim1sig[0])                # setting iso magnitude error of non-detection to limiting mag

                autosncol = [20,31,42,53]                            # list of S/N_auto columns (col number-1)
                if '-F105W' in sys.argv: autosncol = [20,31,42,53,64]
                if ii in autosncol and float(col) < SNnondetection:  # check if S/N is high enough
                    linelist[ii+3] = str(99)                         # setting auto magnitude of non-detection to 99
                    linelist[ii+4] = str(mlim1sig[0])                # setting auto magnitude error of non-detection to limiting mag
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            if test4float(col)==1: newline=newline+'  '+col          # if column is float append
            if test4float(col)==0: newline=newline+'  '+str(9999.99) # writing strings as float 9999.99

        #if Nobj == 460: pdb.set_trace()

        newline=newline+'\n'                # advancing to next line
        fileOUT.write(newline)
    fileOUT.close()
    
    # --- Writing the corresponding *.columns file ---
    fileCOL  = open(bpzcol,"w")
    fileCOL.write("# Filter               columns  AB/Vega  zp_error  zp_offset   \n")

    if maguse == 'auto':
        if vband == 'F600LP': fileCOL.write("f600lp_UVIS1     24, 25   AB        0.01      0.0   \n")  # automag
        if vband == 'F606W':  fileCOL.write("HST_ACS_WFC_F606W     24, 25   AB        0.01      0.0   \n") # automag
        fileCOL.write("f098m_IR     35, 36   AB        0.01      0.0   \n") # automag
        fileCOL.write("HST_WFC3_IR_F125W     46, 47   AB        0.01      0.0   \n") # automag
        fileCOL.write("HST_WFC3_IR_F160W     57, 58   AB        0.01      0.0   \n") # automag
        if '-F105W' in sys.argv: fileCOL.write("HST_WFC3_IR_F105W     68, 69   AB        0.01      0.0   \n") # automag
        fileCOL.write("M_0                   46   \n") # automag


    if maguse == 'iso':
        if vband == 'F600LP': fileCOL.write("f600lp_UVIS1     22, 23   AB        0.01      0.0   \n")  # isomag
        if vband == 'F606W':  fileCOL.write("HST_ACS_WFC_F606W     22, 23   AB        0.01      0.0   \n") # isomag
        fileCOL.write("f098m_IR     33, 34   AB        0.01      0.0   \n") # isomag
        fileCOL.write("HST_WFC3_IR_F125W     44, 45   AB        0.01      0.0   \n") # isomag
        fileCOL.write("HST_WFC3_IR_F160W     55, 56   AB        0.01      0.0   \n") # isomag
        if '-F105W' in sys.argv: fileCOL.write("HST_WFC3_IR_F105W     66, 67   AB        0.01      0.0   \n") # automag
        fileCOL.write("M_0                   44   \n") # isomag

    fileCOL.write("#Z_S                 14   \n")
    fileCOL.write("ID                    1   \n")
    fileCOL.close
#-------------------------------------------------------------------------------------------------------------
if VB == 1:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
