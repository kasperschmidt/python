#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# saveBPZoutput.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Script copying the output from a BPZ run to a seperate folder with the 
# name containing the date and time and a user specified keyword for easy
# recognition. 
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# keyword        : String containing the prefix to the .cat and .columns BPZ input files 
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# -VERBOSE       : set -VERBOSE to get info/messages printed to the screen
#----------------------------
#   OUTPUTS:
#----------------------------
#
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> saveBPZoutput.py 'borg_0110-0224_multiband_BPZinput' -VERBOSE
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2012-10-12  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import sys    # enabling arguments to code
import os     # enabling command line runs and executing other Python scripts with os.system('string')
#----------------------------
#   FUNCTIONS
#----------------------------
#
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
    SS     = str(round(float(time[2]))).split('.')
    HHMMSS = HHMM+str(SS[0])
    DandT = str(date1[1])+''.join(date[1:3])+HHMMSS
    return str(DandT)
#----------------------------
#-
#-------------------------------------------------------------------------------------------------------------
# optional input
VB = 0
if '-VERBOSE' in sys.argv: VB = 1
#-------------------------------------------------------------------------------------------------------------
if VB == 1:
    print ' '
    print ':: '+sys.argv[0]+' :: -- START OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------
Narg = 2
if len(sys.argv) < Narg: # checking that at least Narg arguments are provided
    errorSTR = "Must provide at least "+str(Narg-1)+" arguments to "+sys.argv[0]+" (besides optional arguments) "
    sys.exit(errorSTR)   # stop the program and print an error message
#-------------------------------------------------------------------------------------------------------------
keyword = sys.argv[1]
DT = DandTstr()
dirname = DT+'_BPZrun_'+keyword

# running bash commands to create dir and move dir/files
os.system('mkdir '+dirname)
if os.path.isdir('html'):                      os.system('mv html '+dirname)
if os.path.isdir('plots'):                     os.system('mv plots '+dirname)
if os.path.isfile(keyword+'.bpz.bak'):         os.system('mv '+keyword+'.bpz.bak '+dirname)
if os.path.isfile(keyword+'.probs'):           os.system('mv '+keyword+'.probs '+dirname)
if os.path.isfile(keyword+'.bpz'):             os.system('mv '+keyword+'.bpz '+dirname)
if os.path.isfile(keyword+'_bpz.cat'):         os.system('mv '+keyword+'_bpz.cat '+dirname)
if os.path.isfile(keyword+'.flux_comparison'): os.system('mv '+keyword+'.flux_comparison '+dirname)
if os.path.isfile(keyword+'_bpz_DzCOUNTS.txt'):    os.system('mv '+keyword+'_bpz_DzCOUNTS.txt '+dirname)
if os.path.isfile(keyword+'.cat'):             os.system('cp '+keyword+'.cat '+dirname)           # cp photo catalog
if os.path.isfile(keyword+'.columns'):         os.system('cp '+keyword+'.columns '+dirname)       # cp photo catalog


if VB == 1:
    print ' Files moved to the directory: '+dirname
#-------------------------------------------------------------------------------------------------------------
if VB == 1:
    print ' '
    print ':: '+sys.argv[0]+' :: -- END OF PROGRAM -- '
    print ' '
#-------------------------------------------------------------------------------------------------------------


