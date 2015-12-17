#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# web4BPZgifs.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Creating a webpage for the gifs created with plotBPZresults.py
#----------------------------
#   COMMENTS
#----------------------------
# Note: the gifs put
#----------------------------
#   INPUTS:
#----------------------------
# dirlist          : Path AND name of file containing names of directories containing BPZ output.
#                    If the results have been created with runBPZmultiplecats.py the directories
#                    contain names including a time stamp yymmddhhmmss, e.g., 121019120545_
#                    NOTE: same input as plotBPZhistcounts.py uses
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
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
# bash> chmod +x web4BPZgifs.py       (only required once)
# bash> web4BPZgifs.py '/Users/kasperborelloschmidt/work/BoRG/BPZmultitest121015/resultdirectorylist.txt' --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-10-17  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import subprocess # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
import commands   # enable easy capture of commandline output
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
parser.add_argument("dirlist", type=str, help="List of directories to search for plots/*deltaz.gif and plots/*zHIST.gif")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
pnDIRLIST = pathAname(args.dirlist)
buffer    = open(args.dirlist).readlines()  # read file into buffer
Nlines    = len(buffer)                     # total number of lines in file
gifnames  = []
histnames = []
for ii in range(Nlines):
    dirsplit = buffer[ii].split('\n')
    dir      = dirsplit[0]

    gifpath  = dir+'/plots/*deltaz.gif'
    giffile  = commands.getoutput('ls '+gifpath)
    if os.path.isfile(giffile):
        gifnames.append(giffile)
    else:
        print ':: '+sys.argv[0]+' :: The '+giffile+' was not found in '+dir[0]

    histpath  = dir+'/plots/*zHIST.gif'
    histfile = commands.getoutput('ls '+histpath)
    if os.path.isfile(histfile):
        histnames.append(histfile)
    else:
        print ':: '+sys.argv[0]+' :: The '+histfile+' was not found in '+dir[0]

Ngif      = len(gifnames)
Nhist     = len(histnames)
if Ngif != Nhist:
    ERRORMSG = "Didn't find the same number of *deltaz.gif ("+str(Ngif)+') and *zHIST.gif ('+str(Nhist)+') --> ABORTING'
    sys.exit(ERRORMSG)    
if Ngif == 0:
    ERRORMSG = "No files found --> ABORTING"
    sys.exit(ERRORMSG)    

#-------------------------------------------------------------------------------------------------------------
# creating webpage. 

outfile = pnDIRLIST[0]+'/gifsANDhistBPZ_'+pnDIRLIST[1].replace('.txt','.html')
fout    = open(outfile, 'w')

fout.write('<HTML>')
fout.write('<HEAD>')
fout.write('<TITLE>GIF and HIST</TITLE>')
fout.write('</HEAD>')
fout.write('<BODY>')
fout.write('<style type="text/css">')
fout.write('<!--')
fout.write('body,td,th {color: #000000;font-family: Arial, Helvetica, sans-serif;font-size: medium;}')
fout.write('body {background-color: #999999;}')
fout.write('-->')
fout.write('</style>')
fout.write('<HEAD><BODY>')

fout.write('\n')
fout.write('<h2>histograms found in directories in %s</h2>\n\n\n' % args.dirlist)

for ii in range(Ngif):
    dir      = buffer[ii]
    fout.write('Directory: %s' % dir)

    fout.write('<br>\n')

    fout.write('<a href="%s">\n' % gifnames[ii])
    fout.write('  <img src="%s"   border=0 height=400 width=500>' % gifnames[ii])
    fout.write('<a>\n')

    fout.write('<a href="%s">\n' % histnames[ii])
    fout.write('  <img src="%s"   border=0 height=400 width=500>' % histnames[ii])
    fout.write('<a>\n')

    fout.write('<br>\n')
    fout.write('<br>\n\n')
    
fout.write('</BODY>')
fout.write('</HTML>')

fout.close()
#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

