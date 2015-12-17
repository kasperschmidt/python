#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# extractNoiseApertures.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Creating and extracting a sample of noise (empty) apertures for a given fits file.
# Empty regions are determined via the sextractor segmentation map.
# Searches for other fits images (by globbing fitsimage...*_drz.fits and creates
# numpy arrays for these as well.
#
# Extracts randomly positioned apertures as opposed to extractNoiseApertures_ongrid.py
# which extracts the maximum number of apertures on a grid.
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   INPUTS:
#----------------------------
# fitsimage        : Image to extract apertures from
# rmsimage         : RMS image of fitsimage
# sexsegm          : Sextractor segmentation map of fitsimage defining 'empty' regions
# Naper            : The number of apertures to extract
# Raper            : Radius of apertures to extract given in arcsec
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --outputdir      : If set all the outputs wil be written to de given directory instead of './'
# --verbose        : set -verbose to get info/messages printed to the screen
# --show           : showning plots on screen for manipulation and saving
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# aperarray        : Saving numpy data array containing the total count/flux in aperture
#                    for the Naper apertures to file
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> extractNoiseApertures.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_0456-2203/borg_0456-2203_f125w_wfc3ir_drz.fits /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_0456-2203/borg_0456-2203_f125w_wfc3ir_rms.fits  /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_0456-2203/final_sources_F125_segm.fits 10 0.4 --verbose --outputdir testdir2

# bash> extractNoiseApertures.py /Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p0/borg_1033+5051/borg_1033+5051_f125w_wfc3ir_drz.fits /Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p0/borg_1033+5051/borg_1033+5051_f125w_wfc3ir_rms.fits /Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p0/segm/borg_1033+5051_final_sources_F125_segm.fits  10 0.4 --verbose --outputdir testdir2

# bash> extractNoiseApertures.py /Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p0/borg_1245+3356/borg_1245+3356_f125w_wfc3ir_drz.fits /Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p0/borg_1245+3356/borg_1245+3356_f125w_wfc3ir_rms.fits /Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p0/segm/borg_1245+3356_final_sources_F125_segm.fits  10 0.4 --verbose --outputdir testdir

# bash> extractNoiseApertures.py /Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p0/borg_0240-1857/borg_0240-1857_f125w_wfc3ir_drz.fits /Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p0/borg_0240-1857/borg_0240-1857_f125w_wfc3ir_rms.fits /Users/kasperborelloschmidt/work/BoRG/borgdata/version_2p0/segm/borg_0240-1857_final_sources_F125_segm.fits  10 0.4 --verbose --outputdir testdir

#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-07-18  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import numpy as np         # enable opening with genfromtxt
import pdb                 # for debugging with pdb.set_trace()
import pyfits
import pywcs
import glob
import commands
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("fitsimage", type=str, help="Image to extract apertures from")
parser.add_argument("rmsimage", type=str, help="RMS map of fitsimage")
parser.add_argument("sexsegm", type=str, help="Sextractor segmentation map for fitsimage")
parser.add_argument("Naper", type=int, help="Number of apertures to extract")
parser.add_argument("Raper", type=float, help="Radius of apertures in arcsec")
# ---- optional arguments ----
parser.add_argument("--outputdir", type=str, help="Outputdirectory to use instead of './'")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
# set output directory
outdir = './'
if args.outputdir:
    if os.path.isdir(args.outputdir):
        outdir = args.outputdir
    else:
        out = commands.getoutput('mkdir '+args.outputdir)
        if args.verbose: print ' - Created output directory "'+args.outputdir+'" as it did not excist.'
        outdir = args.outputdir
if outdir[-1] != '/': outdir = outdir+'/'
#-------------------------------------------------------------------------------------------------------------
# Get string-base for naming output.
namebase = args.fitsimage.split('/')[-1].split('.fi')[0]
#-------------------------------------------------------------------------------------------------------------
# reading fits files
imgHDU = pyfits.open(args.fitsimage)
img    = imgHDU[0].data
imghdr = imgHDU[0].header  
simg   = img.shape
if args.verbose: print ' - shape of fits image is       ',simg

rmsHDU = pyfits.open(args.rmsimage)
rms    = rmsHDU[0].data
srms   = rms.shape
if args.verbose: print ' - shape of rms image is        ',srms

segHDU = pyfits.open(args.sexsegm)
seg    = segHDU[0].data
sseg   = seg.shape
if args.verbose: print ' - shape of segmentaion map is  ',sseg
#-------------------------------------------------------------------------------------------------------------
# looking for other images
otherimg = glob.glob('/'.join( args.fitsimage.split('/')[0:-1] )+'/*drz.fits')
del otherimg[np.where(np.asarray(otherimg) == args.fitsimage)[0]] # removing the args.fitsimage itself
Nother = len(otherimg)
otherarr = np.zeros((Nother,simg[0],simg[1]))
if Nother > 0:
    if args.verbose: print ' - Found ',Nother,' other images... loading them'
    for oo in range(Nother):
        imgotherHDU = pyfits.open(otherimg[oo])
        imgother    = imgotherHDU[0].data
        simgother   = imgother.shape
        otherarr[oo,:,:] = imgother # filling array
        if args.verbose: print '   shape of fits image is       ',simgother         
#-------------------------------------------------------------------------------------------------------------
# get pixelsize
wcs          = pywcs.WCS(imghdr)
xcoord       = round(simg[1]/2.)
ycoord       = round(simg[0]/2.)
pixcrd       = np.array([[xcoord,ycoord],[xcoord+1,ycoord+1]], np.float_) 
radecs       = wcs.wcs_pix2sky(pixcrd,1)
arcsecperpix = [np.abs(np.diff(radecs[:,1])*3600),np.abs(np.diff(radecs[:,1])*3600)] # arcsec per pix for x (ra) and y (dec)
#-------------------------------------------------------------------------------------------------------------
# aperture size in pixels
Raperpix = [round(args.Raper/arcsecperpix[0]),round(args.Raper/arcsecperpix[1])]
Rpix     = np.max(Raperpix) # using the largest dimension as radius for apertures
if Raperpix[0] != Raperpix[1]:
    print 'WARNING - the pixel radius in x and y is not the same - has the image been squeezed?'
    print '          Dimensions of apertures was changed from ',Raperpix,' pixels to ',[Rpix,Rpix]
    Raperpix = [Rpix,Rpix]
#-------------------------------------------------------------------------------------------------------------
# creating mask of good noise space
if args.verbose: print ' - Create mask of empty regions'
empty             = np.ones(simg)
empty[rms>=1]     = 0
empty[seg>=1]     = 0
#-------------------------------------------------------------------------------------------------------------
# grow empty mask
Ngrow = 5
emptygrow = np.ones(simg)

if args.verbose: print ' - Growing mask by '+str(Ngrow)+' pixels in x and y'
for yy in np.arange(-Ngrow,Ngrow,1):
    for xx in np.arange(-Ngrow,Ngrow,1):
            emptygrow = emptygrow * np.roll(np.roll(empty,yy,axis=0),xx,axis=1)
Ngoodpix = len(np.where(emptygrow == 1)[0])
Nbadpix  = len(np.where(emptygrow == 1)[0])

hdu = pyfits.PrimaryHDU(emptygrow)
emptymaskname = outdir+namebase+'_emptymaskgrow.fits'
hdu.writeto(emptymaskname,clobber=True)
if args.verbose: print ' - Wrote empty pixel mask to '+emptymaskname
#-------------------------------------------------------------------------------------------------------------
# put down and find apertures
aperarr  = np.zeros((args.Naper,3+Nother)) # array to contain output
apermask = np.zeros(simg)
pixratio = np.float(len(np.where(emptygrow == 0)[0]))/len(np.where(emptygrow == 1)[0]) # ratio of bad to good pixels

Nrand = args.Naper * int(np.ceil(pixratio)) * 1000
xval  = np.random.rand(Nrand)*(simg[1]-int(2*Rpix))+int(Rpix) # the -2Rpix + Rpix makes sure edges are avoided
xent  = xval.astype(int)
yval  = np.random.rand(Nrand)*(simg[0]-int(2*Rpix))+int(Rpix) # the -2Rpix + Rpix makes sure edges are avoided
yent  = yval.astype(int)

counter = 0 # resetting counter
if args.verbose: print ' - Putting down test apertures; saving good ones '
for ii in xrange(Nrand):
    if counter == args.Naper: continue
    subarr = emptygrow[yent[ii]-Rpix:yent[ii]+Rpix,xent[ii]-Rpix:xent[ii]+Rpix]

    if (subarr == 1).all():
        y,x = np.ogrid[-yent[ii]:simg[0]-yent[ii],-xent[ii]:simg[1]-xent[ii]] #creating mesh vectos around position
        mask = x*x + y*y <= Rpix*Rpix  # creating circle mask
        
        sumpix = np.sum(img[mask])
        # ---- saving values to array ----
        if Nother == 0:
            aperarr[counter,:] = xent[ii],yent[ii],sumpix
        else:
            for jj in xrange(Nother):
                sumpix = np.append(sumpix,np.sum(otherarr[jj,mask]))
            aperarr[counter,:] = np.append(np.array([xent[ii],yent[ii]]),sumpix[:])
        # --------------------------------
        counter = counter + 1 # keeping track of save apertures
        apermask[mask] = 1 # drawing aperture in apermask image
        
    if (ii/float(args.Naper) == round(ii/float(args.Naper))) and (ii > 0) and args.verbose: print '     aperture ',ii,' tested (found',counter,'so far)'
hdu = pyfits.PrimaryHDU(apermask)
apermaskname = outdir+namebase+'_aperturemask.fits'
hdu.writeto(apermaskname,clobber=True)
if args.verbose: print ' - Wrote aperture mask to '+apermaskname
#-------------------------------------------------------------------------------------------------------------
# save aperture array
if counter == args.Naper: # only create npy file if the requested number of good apertures were found
    aperarrname = outdir+namebase+'_aperturearray.txt'
    hdrout = 'columns are:  \nxpix  ypix  apertureval_x_Nfiles \nwhere the files are: \nx y '+args.fitsimage+' '+' '.join(otherimg)
    np.savetxt(aperarrname,aperarr,fmt='%.5f',header=hdrout)
    if args.verbose: print ' - Saved aperture array to '+aperarrname
else:
    if args.verbose: 
        print ' - Did not find the requested ',args.Naper,' apertures among the ',Nrand,' created'
        print '       --> no npz file outputted'
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------

