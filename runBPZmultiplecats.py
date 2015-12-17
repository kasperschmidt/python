#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# runBPZmultiplecats.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Routine for running the BPZ code on multiple input catalogs. The
# catalogs are expected to be on the form required by BPZ, i.e., a
# *.cat file and a *.column file for each catalog
#----------------------------
#   COMMENTS
#----------------------------
# In the case of BoRG catalogs the files can be prepared using 
# ~/work/BoRG/borgcat2bpzcat.py
#
# code has to be run in directory where *.cat and *.columns files are located.
#----------------------------
#   INPUTS:
#----------------------------
# catlist          : Path AND name of file containing names of catalogs to run BPZ for.
#                    The file should contain the names used for the .cat and .column files.
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --noplots        : set this keyword to skip the webpage and plot creation (to save time)
# --keywords       : --keywords 'string' to specify the BPZ keywords to run the code with.
#                    If none given a default setup will be used.
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program at specified place
# --help           : printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# Directoryfile    : Output file containing the name of the directories the output has been
#                    put into. Directories will be moved to directory of catlist
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x runBPZmultiplecats.py       (only required once)
# bash> runBPZmultiplecats.py '/Users/kasperborelloschmidt/work/BoRG/BPZmultitest121015/BoRGcatalognamesTEST.txt' --verbose --keywords '-INTERP 5 -SOMETHING more'

# Running simulated high-z galaxies from Michele 130223
# bash> runBPZmultiplecats.py /Users/kasperborelloschmidt/work/BoRG/simobjects_redshift9/catalognames.txt  --verbose --noplots
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-10-15  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import subprocess # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
import gc
import pdb        # for debugging
#----------------------------
#   FUNCTIONS
#----------------------------
def replaceAll(file,searchExp,replaceExp):
    # search and replace in file
    import fileinput
    import sys
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)
def pathAname(str):
    # splitting string with path and name in to
    strsplit = str.split('/')               # splitting string
    name     = strsplit[-1]                 # saving filename (last entry of catsplit)
    slash    = '/'
    path=slash.join(strsplit[0:-1])+'/'     # putting path back together
    return [path,name]
def DandTstr():
    # Creating a string with date and time on the format yymmddhhmmss
    import datetime
    now    = datetime.datetime.now()
    now    = str(now)
    now0   = now.split(' ')
    date   = now0[0].split('-')
    date1  = date[0].split('20')
    time   = now0[1].split(':')
    HHMM   = ''.join(time[0:2])
    SS     = time[2].split('.')
    HHMMSS = HHMM+str(SS[0])
    DandT = str(date1[1])+''.join(date[1:3])+HHMMSS
    return str(DandT)
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("catlist", type=str, help="Provide path and name of file containing list of catalog names")
# ---- optional arguments ----
parser.add_argument("-np", "--noplots", action="store_true", help="Skip website and plot creation")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("-k", "--keywords", type=str, help="Provide list of keywords to run BPZ with in string")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
# reading input file
buffer    = open(sys.argv[1]).readlines()   # read file into buffer
pAnINPUT  = pathAname(sys.argv[1])          # splitting input file in name and path
Nlines    = len(buffer)                     # total number of lines in file
catnames  = []                              # defining list to put catalog names in
Ncatalogs = 0                               # resetting counter
for line in buffer:
    if Ncatalogs>=Nlines : break            # if all lines read, break out of loop
    if line[0]=='#': continue               # if comment continue to next line
    Ncatalogs=Ncatalogs+1
    catnames.append(line.rstrip())          # appending catalog to list; rstrip removes the \n for new lines
if args.verbose: print ':: '+sys.argv[0]+' :: Found',Ncatalogs,'catalogs in input file'
#-------------------------------------------------------------------------------------------------------------
# Running BPZ on catalogs
dirnames  = []                              # defining list to put directory names in
for cat in catnames:                         # looping over catalogs
    # checking if BPZ keywords are given on the command line - if not set to default
    if args.keywords: 
        keywords = args.keywords
    else:
        keywords = ' -INTERP 2'

    cmdstr = 'python $BPZPATH/bpz.py '+cat+'.cat '+keywords
    os.system(cmdstr)  # running BPZ on catalog with keywords/run-flags
        # KBS ------- Here one should be able to adjust the priors --------
    os.system('python $BPZPATH/bpzfinalize.py '+cat) # finalizing BPZ run and creating output

    # Making sure that no values of exact zero excist in output from BPZ (makes webpage.py choke)
    searchAreplaceFILE = cat+'.flux_comparison'
    searchExp          = '0.000e+00'
    replaceExp         = '1.000e-16'
    replaceAll(searchAreplaceFILE,searchExp,replaceExp)

    # --- CREATE PLOTS IF REQUESTED ---
    if args.noplots:
        if args.verbose: print ':: '+sys.argv[0]+' :: Skipping website and plot creation'
    else:
        if args.verbose: print ':: '+sys.argv[0]+' :: Creating website, plots and gifs for objects'
        # creating SED plots and webpage for BPZ run
        os.system('webpageBPZ.py '+cat+' --verbose')                    

        # creating probability plots for webpage
        os.system('plotBPZprobs.py '+pAnINPUT[0]+cat+'.probs ./html/probplots/ --verbose')  

        # creating ra,dec,zspec plots and gif
        os.system('plotBPZresults.py '+pAnINPUT[0]+cat+'_bpz.cat --verbose --png --Nzbin 50 --savecounts')  
    # ---------------------------------

    # Copying output to a seperate directory
    cmdstr  = 'saveBPZoutput.py '+cat+' -VERBOSE'
    proc    = subprocess.Popen(cmdstr, stdout=subprocess.PIPE, shell=True)
    output  = proc.stdout.read()
    outlist = output.split('\n')  # getting the individual lines of the output

    for line in outlist:
        if line.startswith(" Files moved to the directory: "):
            linesplit = line.split(' ')
            outdir = linesplit[-1]
            if os.path.exists(pAnINPUT[0]+outdir): # checking that the directory was created in the right place
                if args.verbose: print ':: '+sys.argv[0]+' :: Ouput directory was put in the right place'
            else:                                  # if not move it
                os.system('mv '+outdir+' '+pAnINPUT[0]+outdir)
                outdir = pAnINPUT[0]+outdir        # and change outdir value
            break # exiting loop
    dirnames.append(outdir)
#-------------------------------------------------------------------------------------------------------------
# Writing directory names to output file
datetime      = DandTstr()
pAnCODE       = pathAname(sys.argv[1])          # splitting input file in name and path
codename      = (pAnCODE[1]).split('.')
outname       = pAnINPUT[0]+codename[0]+'_'+datetime+'OUTPUT.txt'
directoryfile = open(outname,"w")             # open output file (overwrites if already exists - creates if doesn't exist)
directoryfile.write('# Directories created running '+sys.argv[0]+' on '+datetime+' \n')
directoryfile.write('# \n')
for dir in dirnames:                   # looping over directories
    directoryfile.write(dir+'\n')
directoryfile.close()

if args.verbose: print ':: '+sys.argv[0]+' :: Wrote list of directories created to '+outname
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
