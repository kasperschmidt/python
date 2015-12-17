#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# webpageBPZ.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Creating plots of BPZ results and putting them together in a webpage
#----------------------------
#   COMMENTS
#----------------------------
# In the case of BoRG catalogs the files can be prepared using 
# borgcat2bpzcat.py
#----------------------------
#   INPUTS:
#----------------------------
# catalog          : The output catalog from bpzfinalize.py, i.e., *_bpz.cat
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
# bash> chmod +x webpageBPZ.py       (only required once)
# bash> webpageBPZ.py borg_0110-0224_multiband_BPZinputTESTfull --verbose
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
sys.path.append('/Users/kasperborelloschmidt/work/BoRG/BPZ/bpz-1.99.3/plots/')  # adding path instead of to PYTHONPATH
import sedplotAB
import probplot
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
parser.add_argument("catalog", type=str, help="The output catalog from bpzfinalize.py, i.e., *_bpz.cat")
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
buffer    = open(args.catalog+'_bpz.cat').readlines()  # read file into buffer
Nlines    = len(buffer)                     # total number of lines in file
objid     = []                              # empty list for good indexes
zb        = []                              # empty list for good indexes
zbmin     = []                              # empty list for good indexes
zbmax     = []                              # empty list for good indexes
tb        = []                              # empty list for good indexes
odds      = []                              # empty list for good indexes
chisq2    = []                              # empty list for good indexes
for ll in range(Nlines):
#for ll in range(18):
    line  = buffer[ll]
    if line[0] != '#':
        linesplit = line.split()
        # devning lists of variables
        objid.append(linesplit[0])
        zb.append(linesplit[1])
        zbmin.append(linesplit[2])
        zbmax.append(linesplit[3])
        tb.append(linesplit[4])
        odds.append(linesplit[5])
        chisq2.append(linesplit[10])
Nobj      = len(objid)
pAnINPUT  = pathAname(args.catalog)          # splitting input file in name and path
#-------------------------------------------------------------------------------------------------------------
outdir    = 'html/'
if os.path.exists(pAnINPUT[0]+outdir): # checking if the directory excists
    if args.verbose: print ':: '+sys.argv[0]+' :: html directory already exists'
else:                                  # if not move it
    if args.verbose: print ':: '+sys.argv[0]+' :: creating html/ directory'
    os.mkdir(outdir)

SEDoutdir = outdir+'/sedplots'
if os.path.exists(pAnINPUT[0]+SEDoutdir): # checking if the directory excists
    if args.verbose: print ':: '+sys.argv[0]+' :: html/sedplots directory already exists'
else:                                  # if not move it
    if args.verbose: print ':: '+sys.argv[0]+' :: creating html/sedplot directory'
    os.mkdir(SEDoutdir)

PROBoutdir= outdir+'/probplots'
if os.path.exists(pAnINPUT[0]+PROBoutdir): # checking if the directory excists
    if args.verbose: print ':: '+sys.argv[0]+' :: html/probplots directory already exists'
else:                                  # if not move it
    if args.verbose: print ':: '+sys.argv[0]+' :: creating html/probplots directory'
    os.mkdir(PROBoutdir)
#-------------------------------------------------------------------------------------------------------------
b=sedplotAB.bpzPlots(pAnINPUT[1], objid, probs=None) # no need to loop over sedplotAB.py
b.flux_comparison_plots(show_plots=0, save_plots='png', colors={}, nomargins=0, outdir=SEDoutdir)

# note - comment out this and use plotBPZprobs.py to create the probability plots instead
#for jj in range(Nobj): # 
    #os.system('python $BPZPATH/plots/probplot.py '+pAnINPUT[1]+' '+objid[jj]+' "'+PROBoutdir+'" -SAVE') # SLOOOOW
    #probplot.probplot(pAnINPUT[1], int(objid[jj]), nomargins=0, outdir=PROBoutdir) # QUICKER; but chockes after 475 objects - memory leak when looping...

#-------------------------------------------------------------------------------------------------------------
# creating the actual webpage. 

outfile = outdir+'index.html'
fout    = open(outfile, 'w')

fout.write('<HTML>')
fout.write('<HEAD>')
fout.write('<TITLE>PLOTS</TITLE>')
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
fout.write('<h2>BPZ results for %s.cat</h2>\n\n' % pAnINPUT[1])

for i in range(Nobj):
    id = objid[i]
    fout.write('Object #%s' % str(id))
    fout.write(' &nbsp; BPZ = %.2f' % float(zb[i]))
    fout.write(' [%.2f--%.2f]' % (float(zbmin[i]), float(zbmax[i])))
    fout.write(' &nbsp; type = %.2f' % float(tb[i]))
    fout.write(' &nbsp; ODDS = %.2f' % float(odds[i]))
    fout.write(' &nbsp; chisq2 = %.2f' % float(chisq2[i]))
    #fout.write(' &nbsp; spec-z = %.2f' % float(zspec[i]))

    fout.write('<br>\n')

    fout.write('<a href="sedplots/%s_sed_%d.png">\n' % (pAnINPUT[1], int(id)))
    fout.write('  <img src="sedplots/%s_sed_%d.png"   border=0 height=300 width=400>' % (pAnINPUT[1], int(id)))
    fout.write('<a>\n')

    fout.write('<a href="probplots/probplot%d.png">\n ' % (int(id)))
    fout.write(' <img src="probplots/probplot%d.png" border=0 height=300 width=400>' % (int(id)))
    fout.write('</a>\n ')
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

