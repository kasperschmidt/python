#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# extractNoiseApertures_ongrid.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Creating and extracting a sample of noise (empty) apertures for a given fits file.
# Empty regions are determined via the sextractor segmentation map.
#
# Extracts the maximum number of apertures on a grid as opposed to extractNoiseApertures.py
# which extracts randomly positioned apertures.
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
# bash> extractNoiseApertures_ongrid.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_0456-2203/borg_0456-2203_f125w_wfc3ir_drz.fits /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_0456-2203/borg_0456-2203_f125w_wfc3ir_rms.fits  /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_0456-2203/final_sources_F125_segm.fits 0.32 --verbose --outputdir testdir2

#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-08-30  started by K. B. Schmidt (UCSB) (based on extractNoiseApertures.py)
# 2013-09-13  added constraint on good pixels from all images. K. B. Schmidt
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
from time import localtime, strftime
#-------------------------------------------------------------------------------------------------------------
def createsuperrms(imgname,dataarr,rms,verbose=False):
    '''
    Script creating the combined rms of the V, Y, J and H band for the particular field

    imgname     the names of the images put in dataarr
    dataarr     a cube of the data of the individual rms maps
    '''
    superrms = rms.copy()
    if verbose: print ' - Added initial (F125W) rms map to super rms '        
    flag606  = 0 #flag to indicate if a f606w image is present
    
    for ii in xrange(len(imgname)): # looping over images
        if '_f606w' in imgname[ii]:
            superrms = rms + otherrmsarr[ii,:,:]
            if verbose: print ' - Added ',imgname[ii],' to super rms '
            flag606  = 1 # set flag so we know a f606w (V) image has been added
        if '_f160w' in imgname[ii]:
            superrms = rms + otherrmsarr[ii,:,:]
            if verbose: print ' - Added rms for ',imgname[ii],' to super rms '
        if '_f098m' in imgname[ii]:
            superrms = rms + otherrmsarr[ii,:,:]
            if verbose: print ' - Added rms for ',imgname[ii],' to super rms '

    if flag606 == 0: #In case no f606w image was found add the f600lp as V-band
        for ii in xrange(len(imgname)): # looping over images
            if '_f600lp' in imgname[ii]:
                superrms = rms + otherrmsarr[ii,:,:]
                if verbose: print ' - Added ',imgname[ii],' to super rms '        

    return superrms
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("fitsimage", type=str, help="Image to extract apertures from")
parser.add_argument("rmsimage", type=str, help="RMS map of fitsimage")
parser.add_argument("sexsegm", type=str, help="Sextractor segmentation map for fitsimage")
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
Nother      = len(otherimg)
otherarr    = np.zeros((Nother,simg[0],simg[1]))
otherrmsarr = np.zeros((Nother,simg[0],simg[1]))
if Nother > 0:
    if args.verbose: print ' - Found ',Nother,' other images... loading them and their rms maps'
    for oo in range(Nother):
        imgotherHDU      = pyfits.open(otherimg[oo])
        imgother         = imgotherHDU[0].data
        simgother        = imgother.shape
        otherarr[oo,:,:] = imgother # filling array
        if args.verbose: print '   shape of fits image is       ',simgother

        rmsotherHDU         = pyfits.open(otherimg[oo].replace('drz.fits','rms.fits'))
        rmsother            = rmsotherHDU[0].data
        srmsother           = rmsother.shape
        otherrmsarr[oo,:,:] = rmsother # filling array
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
Npixaper = Rpix*Rpix*np.pi
if Raperpix[0] != Raperpix[1]:
    print 'WARNING - the pixel radius in x and y is not the same - has the image been squeezed?'
    print '          Dimensions of apertures was changed from ',Raperpix,' pixels to ',[Rpix,Rpix]
    Raperpix = [Rpix,Rpix]
#-------------------------------------------------------------------------------------------------------------
# Constraining good pixels based on all images by creating super-rms map of V, Y, J and H rms maps
if Nother != 0:
    if args.verbose: print ' - Making mask to ensure coverage of aperture in all images'
    superrms = createsuperrms(otherimg,otherrmsarr,rms,verbose=True)
    #for ii in xrange(Nother):
    #    ent0 = np.where(otherarr[ii,:,:] == 0)  # where image is 0 (no coverage)
    #    img[ent0] = 0.0 

#hdu = pyfits.PrimaryHDU(superrms)
#testname = outdir+namebase+'_testimgsuperrms.fits'
#hdu.writeto(testname,clobber=True)
#if args.verbose: print ' - Wrote test image to '+testname
#-------------------------------------------------------------------------------------------------------------
# creating mask of good noise space
if args.verbose: print ' - Create mask of empty regions (1 = good, 0 = bad)'
empty             = np.ones(simg)
empty[seg>=1]     = 0
#empty[superrms>=np.max(superrms)] = 0 
empty[~np.isfinite(superrms)] = 0 # regions where no bad pixels exist and there is coverage in both V,Y,J and H
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
apermask = np.zeros(simg)
pixratio = np.float(len(np.where(emptygrow == 0)[0]))/len(np.where(emptygrow == 1)[0]) # ratio of bad to good pixels

edgeskip = 2 # number of apertures to skip at edges of img array

Nxaper = np.int(np.floor(simg[1]/(2*Rpix)))-edgeskip # number of apertures in x (subtracting 4 to clear edges of img array)
Nyaper = np.int(np.floor(simg[0]/(2*Rpix)))-edgeskip # number of apertures in y (subtracting 4 to clear edges of img array)

xent = np.arange(Nxaper)*2.0*Rpix+2*edgeskip*Rpix # x centers of apertures
yent = np.arange(Nyaper)*2.0*Rpix+2*edgeskip*Rpix # y centers of apertures

xent = xent.astype(int)
yent = yent.astype(int)

Nmaxaper = Nxaper*Nyaper
aperarr  = np.zeros((Nmaxaper,3+Nother))-9999 # array to fill

counter = 0 # reset counter
yogrid, xogrid = np.ogrid[0:simg[0],0:simg[1]] #creating mesh vectos around position

yogrid, xogrid = np.ogrid[-Rpix:Rpix+1,-Rpix:Rpix+1]
circlemask     = xogrid**2 + yogrid**2 <= Rpix*Rpix         # creating circle mask in subarray
circledraw     = np.zeros(circlemask.shape)
circledraw[circlemask] = 1.0 

for xx in xrange(Nxaper):
#for xx in [79]:
    if (xx/float(100) == round(xx/float(100))) and args.verbose: print '   Aperture column ',xx
    for yy in xrange(Nyaper):
    #for yy in np.arange(5)+100:
        empty_sub = emptygrow[yent[yy]-Rpix:yent[yy]+Rpix+1,xent[xx]-Rpix:xent[xx]+Rpix+1]
        img_sub   = img[yent[yy]-Rpix:yent[yy]+Rpix+1,xent[xx]-Rpix:xent[xx]+Rpix+1]
        
        Nbadpix   = len(np.where(empty_sub[circlemask] == 0)[0]) # counting pixels that are bad
        N0spix    = len(np.where(img_sub[circlemask] == 0)[0])   # counting pixels with values 0 (to check for aper outside img)
        if (Nbadpix == 0):          # making sure all pixels are good    # and < 10% are 0s : "and (N0spix < Npixaper*0.10)"
            sumpix  = np.sum(img_sub[circlemask])
            # ---- saving values to array ----
            if Nother == 0:
                aperarr[counter,:] = xent[xx],yent[yy],sumpix
            else:
                for jj in xrange(Nother): # if other files detected get aperture flux from those
                    oarr_sub = otherarr[jj,yent[yy]-Rpix:yent[yy]+Rpix+1,xent[xx]-Rpix:xent[xx]+Rpix+1]
                    sumoarr  = np.sum(oarr_sub[circlemask])
                    sumpix = np.append(sumpix,sumoarr)
                aperarr[counter,:] = np.append(np.array([xent[xx],yent[yy]]),sumpix[:])
            # --------------------------------
            counter = counter + 1 # keeping track of saved apertures
            apermask[yent[yy]-Rpix:yent[yy]+Rpix+1,xent[xx]-Rpix:xent[xx]+Rpix+1] = circledraw # drawing aperture in apermask image
            
if args.verbose: print ' - Found ',counter,' apertures in total '    
hdu = pyfits.PrimaryHDU(apermask)
apermaskname = outdir+namebase+'_aperturemask.fits'
hdu.writeto(apermaskname,clobber=True)
if args.verbose: print ' - Wrote aperture mask to '+apermaskname
#-------------------------------------------------------------------------------------------------------------
# save aperture array
aperarr     = aperarr[np.where(aperarr[:,0] != -9999)[0],:]
aperarrname = outdir+namebase+'_aperturearray.txt'

hdrout = 'columns are:  \nxpix  ypix  apertureval_x_Nfiles \nwhere the files are: \nx y '+args.fitsimage+' '+' '.join(otherimg)
np.savetxt(aperarrname,aperarr,fmt='%.5f',header=hdrout)
if args.verbose: print ' - Saved aperture array to '+aperarrname
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------

