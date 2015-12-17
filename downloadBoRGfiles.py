#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# downloadBoRGfiles.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Script for downloading all files of a certaing kind from the BoRG data repository
#----------------------------
#   COMMENTS
#----------------------------
# As asterisk are not working with wget on the BoRG server only full filenames are supported.
# We circumvented this by simply downloading the vull version 2.1 catalog to 
# /Users/kasperborelloschmidt/work/BoRG/borgdata/version_2.1
# with 
# nohup scp -r "trenti@wolf359.colorado.edu:/data/borg_web_archive/full_dataset/version_2.1/" .
# ^z
# bg

#
# JUST USE:
#
# bash> wget --http-user=BoRG_CoI --http-password=Loc7o91n3! --no-check-certificate -r --no-host-directories https://wolf359.colorado.edu/borg_web_archive/full_dataset/version_2.1/borg_0110-0224/
#
#----------------------------
#   INPUTS:
#----------------------------
# filestring       : Part of string to search for in the BoRG directories. For instance 'final_sources*_F125.cat'
#                    will download all the source catalogs from SExtractor - both cleaned and non-cleaned.
#                    NB! has to be the full filename... asterisks are not supported.
# outputdir        : Directories to download files to (assumes that it already exists) - end with '/'
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program at specified place
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# Directoryfile    : Output file containing the name of the directories the output has been
#                    put into. Directories will be moved to directory of catlist
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x downloadBoRGfiles.py       (only required once)
# bash> downloadBoRGfiles.py '*final_sources*_F125.cat' '/Users/kasperborelloschmidt/work/BoRG/borg_cats_v2p1/' --verbose
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-11-06  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse   # argument managing
import sys        # enabling arguments to code
import os         # enabling command line runs and executing other Python scripts with os.system('string')
import commands   # get output from spawned command line processes
import getopt     # used to extract/obtain the optional input
#----------------------------
#   FUNCTIONS
#----------------------------
#
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("filesting", type=str, help="String contained in file to download")
parser.add_argument("outputdir", type=str, help="Directory to donwload files to")
# ---- optional arguments ----
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program at specified place")
parser.add_argument("-k", "--keywords", type=str, help="Provide list of keywords to run BPZ with in string")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
borgurl = 'https://wolf359.colorado.edu/borg_web_archive/full_dataset/version_2.1/'
usrpsw  = ' --http-user=BoRG_CoI --http-password=Loc7o91n3! --no-check-certificate  '

# list of subdirectories to search through
fields  = ['borg_0110-0224/','borg_0214+1255/','borg_0228-4102/','borg_0240-1857/','borg_0427+2538/','borg_0436-5259/','borg_0439-5317/','borg_0440-5244/','borg_0540-6409/','borg_0553-6405/','borg_0624-6432/','borg_0624-6440/','borg_0637-7518/','borg_0751+2917/','borg_0756+3043/','borg_0808+3946/','borg_0819+4911/','borg_0820+2332/','borg_0835+2456/','borg_0846+7654/','borg_0906+0255/','borg_0909+0002/','borg_0914+2822/','borg_0922+4505/','borg_0926+4000/','borg_0926+4426/','borg_1010+3001/','borg_1014-0423/','borg_1031+3804/','borg_1031+5052/','borg_1033+5051/','borg_1051+3359/','borg_1103-2330/','borg_1111+5545/','borg_1119+4026/','borg_1131+3114/','borg_1152+5441/','borg_1153+0056/','borg_1209+4543/','borg_1230+0750/','borg_1242+5716/','borg_1245+3356/','borg_1301+0000/','borg_1337+0028/','borg_1341+4123/','borg_1408+5503/','borg_1437+5043/','borg_1510+1115/','borg_1524+0954/','borg_1555+1108/','borg_1632+3733/','borg_1632+3737/','borg_1815-3244/','borg_2057-4412/','borg_2132+1004/','borg_2155-4411/','borg_2203+1851/','borg_2345+0054/','borg_2351-4332/']
Nfields = len(fields)
if args.verbose: print ':: '+sys.argv[0]+' :: Will loop over and search in '+str(Nfields)+' subdirectories at '+borgurl

presentdir = commands.getoutput('pwd')
os.system('cd '+args.outputdir) # changing directory to donwload directory
 
for ii in range(Nfields):  # looping over the borg fields in the data repository
    cmdstring = 'wget '+usrpsw+borgurl+fields[ii]+args.filesting  # string to parse
    os.system(cmdstring) # downloading files

os.system('cd '+presentdir) # changing directory back after download

#-------------------------------------------------------------------------------------------------------------
if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------

