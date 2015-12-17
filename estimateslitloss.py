#!/usr/bin/env python2.7
#+
#----------------------------
#   NAME
#----------------------------
# estimateslitloss.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Estimating slitloss of object given a fits image, a specified slit
# (extraction aperture/box) and an object poistion.
#----------------------------
#   COMMENTS
#----------------------------
# 
#----------------------------
#   INPUTS:
#----------------------------
# fitsimage        : Fits image containing object
# seeing           : Seeing FWHM of slit observations used to covolve the input fits image
# slit             : The slit/extraction box's dimensions in pixels (of input image). 
#                    Expects: width length angle. Angle corresponds to the DS9 region box angle in degrees.
#                    (The pixel scale of drizzled BoRG data is 0.0800 ''/pixel)
#                    (The spatial pixel scale of MOSFIRE is 0.1799 ''/pixel)
#                    A 0.7''x15(MOSFIREpix) slit has HST dimensions 0.7/0.08 x 15*0.1799/0.08 ~ 9 x 34
# pixelpos         : Position of object (where to put down slit) in pixels. Expects: x y (counting from 1)
#                    These can come from e.g. a sextractor catalog of the fits image.
# magobj           : the magnitude of the object and the zeropoint of the band it's calucalted in.
#                    The zeropoints used in BoRG are:  ZP098, ZP125, ZP160 = 25.68, 26.25, 25.96
#                    Used to get the actual slitloss, i.e. the fraction of flux lost due to the slit/aperture
#----------------------------
#   OPTIONAL INPUTS:
#----------------------------
# --relative       : Give the fraction of light _lost_ for objects (dropouts) which the Magnitude adjustment
#                    has to be estimated relative to.
# --verbose        : set -verbose to get info/messages printed to the screen
# --stop           : stoppping program before end for de-bugging
# --help           : Printing help menu
#----------------------------
#   OUTPUTS:
#----------------------------
# Fraction of estimated flux lost due to the slit
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# --- borg_1510+1115_0745 ---
# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_1510+1115/borg_1510+1115_f125w_wfc3ir_drz.fits 0.9 9 34 -3.000 1704 1210 19.29 26.25 --verbose  --relative 0.36

# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_1510+1115/borg_1510+1115_f098m_wfc3ir_drz.fits 0.9 9 34 -3.000 1704 1210 19.57 25.68 --verbose --relative 0.36

# --- 1510+1115_0354 Y >27.83  J 27.03 \pm 0.22 ---
# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_1510+1115/borg_1510+1115_f125w_wfc3ir_drz.fits 0.9 9 34 -3.000 722.287 701.621 27.03 26.25 --verbose

# --- 1510+1115_1218 Y >27.83  J 26.87 \pm 0.22 ---
# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_1510+1115/borg_1510+1115_f125w_wfc3ir_drz.fits 0.9 9 34 -3.000 916.363  2054.621 26.87 26.25 --verbose

# --- 1510+1115_1487 Y >27.83  J 27.60 \pm 0.24 ---
# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_1510+1115/borg_1510+1115_f125w_wfc3ir_drz.fits 0.9 9 34 -3.000 1398.763  1650.497 27.60 26.25 --verbose

# --- 1510+1115_1524 Y >27.83  J 26.63 \pm 0.15 ---
# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_1510+1115/borg_1510+1115_f125w_wfc3ir_drz.fits 0.9 9 34 -3.000 1117 1785 26.63 26.25 --verbose

# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_1510+1115/borg_1510+1115_f098m_wfc3ir_drz.fits 0.9 9 34 -3.000 1117 1785 27.83 25.68 --verbose

# --- 1510+1115_1705 Y >27.83  J 27.00 \pm 0.19 ---
# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_1510+1115/borg_1510+1115_f125w_wfc3ir_drz.fits 0.9 9 34 -3.000 1030.535  1586.088 27.00 26.25 --verbose


# --- borg_0951+3304_0204 --- bright no. 1
# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_0951+3304/borg_0951+3304_f098m_wfc3ir_drz.fits 0.9 9 34 -62.000 2201 1069 18.91 25.68 --verbose  --relative 0.36

# --- borg_0951+3304_0281 --- bright no. 2
# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_0951+3304/borg_0951+3304_f098m_wfc3ir_drz.fits 0.9 9 34 -62.000  153.8  1078.6 20.43 25.68 --verbose --relative 0.36



# --- borg_0951+3304_0277 ---
# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/cycle19_new_nov14_cat/borg_0951+3304/borg_0951+3304_f098m_wfc3ir_drz.fits 0.9 9 34 -62.000 1758 983 26.83 25.68 --verbose


# --- borg_1437+5043_0145_r3 ---
# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region3/borg_1437+5043/borg_1437+5043_f098m_wfc3ir_drz.fits 0.9 9 34 -34.500 2512 1633 21.26 25.68 --verbose  --relative 0.36

# --- borg_1437+5043_r2_0637_T12a ---
# bash> estimateslitloss.py /Users/kasperborelloschmidt/work/BoRG/borgdata/BoRG_1437+5043_followup/cat_byhandmasking130418/region2/borg_1437+5043/borg_1437+5043_f098m_wfc3ir_drz.fits 0.9 9 34 -34.500 1784 1741 28.05 25.68 --verbose


#----------------------------
#   BUGS
#----------------------------
# 130814 - the commands.getoutput for some reason doesn't work for spawning the command that convolves the fits image.
#          Hence, the need to run it manually if non-excisting. 
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-08-28  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import argparse            # argument managing
import sys                 # enabling arguments to code
import os                  # enabling command line runs and executing other Python scripts with os.system('string')
import commands            # get output from spawned command line processes
import getopt              # used to extract/obtain the optional input
import numpy as np         # enable opening with genfromtxt
import pdb                 # for debugging with pdb.set_trace()
import pyfits
#-------------------------------------------------------------------------------------------------------------
def Rotate2D(pts,cnt,ang=3.141593/4.):
    '''
    pts = {} Rotates points(nx2) about center cnt(2) by angle ang(1) in radian
    from http://gis.stackexchange.com/questions/23587/how-do-i-rotate-the-polygon-about-an-anchor-point-using-python-script
    '''
    import scipy
    ar     = scipy.array
    dot    = scipy.dot
    sin    = scipy.sin
    cos    = scipy.cos
    return dot(pts-cnt,ar([[cos(ang),sin(ang)],[-sin(ang),cos(ang)]]))+cnt
def pointsinpoly(corners,imgshape):
    '''
    Creating mask of polygon. Corners give the corners of the polygon in order (e.g. clockwise)
    from http://stackoverflow.com/questions/3654289/scipy-create-2d-polygon-mask
    '''
    import numpy as np
    from matplotlib.nxutils import points_inside_poly
    
    nx, ny = imgshape[1], imgshape[0]
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x,y)).T

    grid = points_inside_poly(points, corners)
    grid = grid.reshape((ny,nx))

    #print grid
    return grid
#-------------------------------------------------------------------------------------------------------------
# Managing arguments with argparse (see http://docs.python.org/howto/argparse.html)
parser = argparse.ArgumentParser()
# ---- required arguments ---- :
parser.add_argument("fitsimage", type=str, help="Fits image containing object")
parser.add_argument("seeing", type=float, help="Seeing FWHM of slit observations")
parser.add_argument("slit", type=float, nargs=3, help="The slit/extraction box dimension in pixels (width length)")
parser.add_argument("pixelpos", type=float, nargs=2, help="Position of object (center of slit) in fits image (x y)")
parser.add_argument("magobj", type=float, nargs=2, help="Magnitude and zeropoint (for band mag is in) of object")
# ---- optional arguments ----
parser.add_argument("--relative", type=float, help="The _lost_ flux the corrected mag should be relative to")
parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose comments")
parser.add_argument("--stop", action="store_true", help="Stopping program before end for debugging")
parser.add_argument("-k", "--keywords", type=str, help="Provide list of keywords to run BPZ with in string")

args = parser.parse_args()
#-------------------------------------------------------------------------------------------------------------
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- START OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------
# reading fits iamge
if args.verbose: print ' - Loading ',args.fitsimage
imgHDU = pyfits.open(args.fitsimage)
img    = imgHDU[0].data
imghdr = imgHDU[0].header  
simg   = img.shape
if args.verbose: print ' - shape of fits image is       ',simg
#-------------------------------------------------------------------------------------------------------------
# Convolve fits image with psf via convolveFITSimages.py
if args.verbose: print ' - Looking for colvoved image'
replacestr   = '_convolvedFWHM'+str(args.seeing).replace('.','p')+'arcsec.fits'
convolvename = args.fitsimage.replace('.fits',replacestr)  # creating name of convolve image

if os.path.isfile(convolvename):
    print ' - Found the concolved fits image : ',convolvename
else:
    pwd = commands.getoutput('pwd')
    dummyfile = pwd+'/dummyimagelist.txt' 
    f = open(dummyfile, 'w')
    f.write(args.fitsimage)
    f.close

    cmd = ' convolveFITSimages.py "'+dummyfile+'" --verbose --kernelsize '+str(args.seeing)+' '
    #cmdout = commands.getoutput(cmd) # Not working???
    #print cmdout
    #os.remove(dummyfile)
    
    print ' - Manually run the command:'
    print '   bash> '+cmd
    print '   and then rerun... --> ABORTING'
    pdb.set_trace()
#-------------------------------------------------------------------------------------------------------------
# Get pixels within slit
if args.verbose: print ' - Getting pixels falling within slit '
center = np.array(args.pixelpos)
width  = args.slit[0]
length = args.slit[1]
angrad = args.slit[2]*np.pi/180.0
points = np.array([[center[0]-width/2.,center[1]+length/2.],
                   [center[0]+width/2.,center[1]+length/2.],
                   [center[0]+width/2.,center[1]-length/2.],
                   [center[0]-width/2.,center[1]-length/2.]])

pointsnew = Rotate2D(points,center,ang=angrad) # rotating slit
mask      = pointsinpoly(pointsnew,simg)       # getting mask of pixels in slit
pixent    = np.where(mask == True)             # getting the corresponding pixel entries

writetest = 1
if writetest == 1:
    test = np.zeros(simg)
    test[pixent] = 1.0
    hdu = pyfits.PrimaryHDU(test)
    hdu.writeto('testofmask.fits',clobber=True)
#-------------------------------------------------------------------------------------------------------------
# Get flux in slit
if args.verbose: print ' - Compare flux in slit with magnitude '+str(args.magobj[0])+' given the zeropoint '+str(args.magobj[1])

flux    = np.sum(img[pixent])
mag     = -2.5*np.log10(flux) + args.magobj[1]
fluxobj = 10**( (args.magobj[0]-args.magobj[1]) / -2.5 )
floss   = (1.0-flux/fluxobj)

magcorr = args.magobj[0] - 2.5*np.log10(flux/fluxobj)
if args.relative:
    magcorrrel = args.magobj[0] - 2.5*np.log10(flux/fluxobj/(1-args.relative))

if args.verbose:
    print '\n - The flux of the input mag is         ',str("%.2f" % fluxobj)
    print ' - The flux in the slit was found to be ',str("%.2f" % flux)
    print '   corresponding to a magnitude of      ',str("%.2f" % mag)
    print '\n   Hence;'
    print ' - The flux ratio between the input and slit value is ',fluxobj/flux
    print ' - I.e. due to the slit a flux loss of                ',str("%.2f" % (floss*100.))+'%'
    print ' - corresponding to a corrected magnitude of          ',str("%.2f" % magcorr)
    if args.relative:
        print '\n - The relative slitloss assuming loosing '+str(args.relative)+'% of light from dropouts'
        print '   then corresponds to an adjusted magnitude of   ',str("%.2f" % magcorrrel)        
#-------------------------------------------------------------------------------------------------------------
if args.stop: pdb.set_trace()
#if args.stop: sys.exit('STOPPED PROGRAM AS REQUESTED')
if args.verbose: print '\n:: '+sys.argv[0]+' :: -- END OF PROGRAM -- \n'
#-------------------------------------------------------------------------------------------------------------

