# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                           Utilities for MUSE-Wide related stuff
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pyfits
import numpy as np
import MUSEWideUtilities as mwu
import sys
import scipy.ndimage
import tdose_utilities as tu
from astropy import wcs
from astropy.coordinates import SkyCoord
import glob
import pdb
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_pixelpos(ra,dec,pointingname,imgdir='/Users/kschmidt/work/images_MAST/MUSEWidePointings/',imgext=0,
                 radecunit="deg",pixorigin=1,verbose=True):
    """
    Get pixel positions in a muse cube for

    """
    searchstr = imgdir+'*'+pointingname+'*.fits'
    img = glob.glob(searchstr)
    if len(img) == 1:
        imghdr = pyfits.open(img[0])[imgext].header
    else:
        sys.exit(' Found '+str(len(img))+' images looking for '+searchstr)

    if verbose: print ' Getting pixel position in '+img[0]
    imgwcs    = wcs.WCS(tu.strip_header(imghdr.copy()))

    pixcoord  = wcs.utils.skycoord_to_pixel(SkyCoord(ra, dec, unit=radecunit),imgwcs,origin=pixorigin)
    xpix      = pixcoord[0]
    ypix      = pixcoord[1]

    return xpix[0], ypix[0]
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_pointingname(id):
    """
    Returning the pointing name for a MUSE-Wide id

    """
    idstr = str(id)
    pointingname = 'candels-'

    # get field
    if idstr[0] == '1':
        pointingname = pointingname+'cdfs-'
    elif idstr[0] == '2':
        pointingname = pointingname+'cosmos-'
    else:
        sys.exit(' mu.gen_pointingname(): Field corresponding to ID[0]='+str(idstr[0])+' not known')

    # get pointing number
    pointingname = pointingname+idstr[1:3]

    return pointingname
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def convert_wht2sigma(wht_file,wht_units='ivar',wht_ext=0,
                      smooth=False,smooth_type='gauss',smooth_param=4,
                      output_file='default',clobber=False,verbose=True):
    """
    Convert an wheight map to a sigma image.
    Add optional smoohting before saving with the smooth_* keywords.

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu
    imgdir     = '/Volumes/DATABCKUP2/TDOSEextractions/tdose_cutouts/'
    wht_file   = imgdir+'acs_814w_candels-cdfs-02_wht_cut_v1.0_id9262_cutout6p0x6p0arcsec.fits'
    outputfile = mwu.convert_wht2sigma(wht_file,verbose=True,clobber=True,smooth=True,smooth_param=3)
    outputfile = mwu.convert_wht2sigma(wht_file,verbose=True,clobber=True,smooth=False)


    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading inverse variance wheight map'
    wht_data            = pyfits.open(wht_file)[wht_ext].data
    wht_data_shape      = wht_data.shape

    if wht_units == 'ivar':
        if verbose: print ' - Weight image is in units "ivar", i.e., inverser variamce. Turning map into sigmas'
        sigma_img = np.sqrt(1.0/wht_data)
    elif (wht_units == 'sigma') or (wht_units == 'stddev'):
        if verbose: print ' - Weight image already in sigma units "sigma" or "stddev". No conversion applied'
        sigma_img = wht_data
    else:
        sys.exit(' ---> Invalid choice of wht_units="'+wht_units+'" ')

    outstr = 'sigma'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if smooth:
        outstr = outstr+'_smooth'
        if smooth_type.lower() == 'gauss':
            outputimg = scipy.ndimage.filters.gaussian_filter(sigma_img, smooth_param,cval=0.0)
            smoothkey = 'gaussSig'+str(smooth_param).replace('.','p')+'pix'
            outstr    = outstr+'_'+smoothkey

        else:
            sys.exit(' ---> Invalid choice of smooth_type="'+smooth_type+'" ')
    else:
        smoothkey = 'None'
        outputimg = sigma_img

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if output_file.lower() == 'default':
        outfile = wht_file.replace('.fits','_'+outstr+'.fits')
    else:
        outfile = output_file

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Updating fits header with conversion information '
    sigma_hdr  = pyfits.open(wht_file)[wht_ext].header
    removekeys = ['FILENAME']
    for key in removekeys:
        if key in sigma_hdr.keys():
            sigma_hdr.remove(key)

    # writing hdrkeys:'---KEY--',                  '----------------MAX LENGTH COMMENT-------------'
    sigma_hdr.append(('WHTTYPE ',wht_units        ,'Assumed weight map units when converting '),end=True)
    sigma_hdr.append(('SMOOTH  ',smoothkey        ,'Information about sigma map smoothing'),end=True)
    #sigma_hdr['OBJECT'] = sigma_hdr['OBJECT'].replace('EXP','VARCUBE')

    mwu.save_img(outputimg,sigma_hdr,outfile,clobber=clobber,verbose=verbose)

    return outfile
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def save_img(imagedata,header,filename,clobber=False,verbose=True):
    """
    Saving image to fits file

    """
    if verbose: print ' - Saving image to \n   '+filename
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if 'XTENSION' in header.keys():
        hduprim        = pyfits.PrimaryHDU()  # default HDU with default minimal header
        hducube        = pyfits.ImageHDU(imagedata,header=header)
        hdus           = [hduprim,hducube]
    else:
        hducube = pyfits.PrimaryHDU(imagedata,header=header)
        hdus           = [hducube]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hdulist = pyfits.HDUList(hdus)             # turn header into to hdulist
    hdulist.writeto(filename,clobber=clobber)  # write fits file (clobber=True overwrites excisting file)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def insert_PSFfits_into_TDOSEinfofile(TDOSEinfofile='/Users/kschmidt/work/TDOSE/muse_tdose_setups/MUSEWide_infofile_arche.txt',
                                      duplicatefirstline=False,verbose=True):
    """
    Generating new TDOSE infofile where the PSF information is updated

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu
    mwu.insert_PSFfits_into_TDOSEinfofile(duplicatefirstline=True)

    """
    infofile   = np.genfromtxt(TDOSEinfofile,skip_header=2,names=True,dtype=None)
    Nsetups    = len(infofile[0])
    fmt        = ['S150']*Nsetups
    infofile   = np.genfromtxt(TDOSEinfofile,skip_header=2,names=True,dtype=fmt)

    outfile    = TDOSEinfofile.replace('.txt','_PSFupdate.txt')
    if verbose: print ' - Will store the updated infofile in \n   '+outfile
    fout       = open(outfile,'w')

    PSFinfo    = pyfits.open('/Users/kschmidt/work/catalogs/MUSE_GTO/psf_all_Converted_cleaned.fits')[1].data

    fieldnames = [nn.split('candels-')[-1] for nn in infofile['setupname']]

    fout.write('# File updated with mwu.insert_PSFfits_into_TDOSEinfofile() \n#\n# '+
               '  '.join(list(infofile.dtype.names))+'\n' )

    if duplicatefirstline:
        fieldnames = PSFinfo['field']

    for ff, fieldname in enumerate(fieldnames):
        if 'hudf' in fieldname: continue
        fieldname = fieldname.split()[0]

        fieldent  = np.where(PSFinfo['field'] == fieldname.ljust(9) )[0]
        p0        = PSFinfo['p0_comb_g'][fieldent][0]
        p1        = PSFinfo['p1_comb_g'][fieldent][0]

        if (p0 == 0.0) & (p1 == 0.0):
            p0        = PSFinfo['p0_PM_g'][fieldent][0]
            p1        = PSFinfo['p1_PM_g'][fieldent][0]
            if (p0 == 0.0) & (p1 == 0.0):
                p0        = PSFinfo['p0_FFTS_g'][fieldent][0]
                p1        = PSFinfo['p1_FFTS_g'][fieldent][0]
                if (p0 == 0.0) & (p1 == 0.0):
                    p0        = PSFinfo['p0_FFT_g'][fieldent][0]
                    p1        = PSFinfo['p1_FFT_g'][fieldent][0]

        if ff == 0:
            infoent   = np.where(infofile['setupname'] == 'candels-'+fieldname)[0]
            infofile['psf_FWHMp0'][infoent] = str(p0)
            infofile['psf_FWHMp1'][infoent] = str(p1)
            stringout = '  '.join(list(infofile[infoent][0]))+'\n'

            if duplicatefirstline:
                basename   = infofile['setupname'][infoent]
                basestring = stringout
                basep0     = str(p0)
                basep1     = str(p1)
        else:
            if duplicatefirstline:
                stringout = basestring.replace(basename[0],'candels-'+fieldname)
                stringout = stringout.replace(basename[0].replace('-','_'),('candels-'+fieldname).replace('-','_'))
                stringout = stringout.replace(basep0,str(p0))
                stringout = stringout.replace(basep1,str(p1))
            else:
                infoent   = np.where(infofile['setupname'] == 'candels-'+fieldname)[0]
                infofile['psf_FWHMp0'][infoent] = str(p0)
                infofile['psf_FWHMp1'][infoent] = str(p1)
                stringout = '  '.join(list(infofile[infoent][0]))+'\n'

        fout.write(stringout)

    fout.close()


    #     if PSFinfo['p0_FFTS_g'][fieldent] != 0.0:
    #         p0 = PSFinfo['p0_FFTS_g'][fieldent]
    #     else:
    #         if PSFinfo['p0_PM_g'][fieldent] != 0.0:
    #             p0 = PSFinfo['p0_PM_g'][fieldent]
    #         else:
    #             if PSFinfo['p0_FFT_g'][fieldent] != 0.0:
    #                 p0 = PSFinfo['p0_FFT_g'][fieldent]
    #             else:
    #                 sys.exit(' No desired PSF fit found for object ')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_specstat(wave,flux,fluxerr):
    """
    Return statistics on a spectrum

    """
    statdic        = {}
    statdic['mea_f']          = np.mean(flux[np.isfinite(flux)])
    statdic['med_f']          = np.median(flux[np.isfinite(flux)])
    statdic['max_f']          = np.max(flux[np.isfinite(flux)])
    statdic['std_f']          = np.std(flux[np.isfinite(flux)])
    statdic['max_f_wave']     = wave[np.where(flux == statdic['max_f'])[0]][0]

    s2n                       = flux/fluxerr
    statdic['s2n']            = s2n
    statdic['mea_s2n']        = np.mean(s2n[np.isfinite(s2n)])
    statdic['med_s2n']        = np.median(s2n[np.isfinite(s2n)])
    statdic['max_s2n']        = np.max(s2n[np.isfinite(s2n)])
    statdic['std_s2n']        = np.std(s2n[np.isfinite(s2n)])
    statdic['max_s2n_wave']   = wave[np.where(s2n == statdic['max_s2n'])[0]][0]

    return statdic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def launch_bluebumpGUI():
    """
    Wrapper to set up and launch Josies GUI to serach for blue bumps in LAE spectra

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mu
    mu.launch_bluebumpGUI()

    """
    print 'ERROR: Cannot be run within ipython using import. Run by simply typing "python determine_bb.py',pdb.set_trace()
    print ' - Assuming positioned in /Users/kschmidt/work/MUSE/Josie_DoublePeakInspection/'
    print ' - Importing determine_bb.py (remember to edit file and catalog info '
    import determine_bb
    print ' - Launching inspection GUI; happy inspecting... '
    determine_bb.main()

    print ' - GUI existed after inspections'

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def create_narrowband_subcube(datacube,ras,decs,dras,ddecs,wavecenters,dwaves,outputdir,
                              cube_ext=['DATA_DCBGC','EFF_STAT'],
                              names=None,clobber=False,verbose=True):
    """
    Create narrow band images (sum over predefined set of wavelength slices) from a given data cube.
    Will generate subcube cutouts in the porcess using tu.extract_subcube()

    --- INPUT ---
    datacube     Data cube to extract narrow band images from
    ras          List of R.A.s to extract narrow band images around
    decs         List of Dec.s to extract narrow band images around
    dras         Size of narrow band images in R.A. direction (provided in arcsec)
    ddecs        Size of narrow band images in Dec. direction (provided in arcsec)
    wavecenters  List of lists of central wavelengths of narrow band images to generate for the N locations/objects
                 provided with ras and decs input. Use format:
                     [[obj1line1,obj1line2,obj1line3],[obj2line1,obj2line2],...]
    dwaves       The wavelength ranges to collabse into narrow band images (by summing) for the wavecenters provided.
                 Use format:
                     [[obj1dwave1,obj1dwave2,obj1dwave3],[obj2dwave1,obj2dwave2],...]
    outputdir    Directory to stroe outputs to
    cube_ext     list of extensions to extract from cube into sub-cube. The header of the first extension mentioned
                 in the list will be used as WCS and Header template for the narrow-band image generated
    names        Names or IDs of locations/objects indicated by the ras and decs input
    clobber      Overwrite existing files
    verbose      Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mu

    pointing = 'candels-cdfs-01'
    outdir   = '/Volumes/DATABCKUP2/TDOSEextractions/MW_LAEs_JKgalfitmodels/'
    ras      = [53.0612466535,53.0540222689]
    decs     = [-27.8047374553,-27.8037275129]
    wcenters = [[5830.0,7426.5,9152.5],[5869.5,7476.9,11000]]
    dwaves   = [[10,10,10],[10,10,10]]
    names    = ['101011026','101012027']
    datacube = glob.glob('/Volumes/DATABCKUP2/MUSE-Wide/datacubes_dcbgc_effnoised/DATACUBE_'+pointing+'_v1.0_dcbgc_effnoised.fits')

    mu.create_narrowband_subcube(datacube[0],ras,decs,5.0,5.0,wcenters,dwaves,outdir,names=names,clobber=True)

    """
    Nsubcubes = len(ras)
    if len(decs) != Nsubcubes:
        sys.exit(' The number of RAs does not match the number of Decs. Make sure they do ')
    if verbose: print ' - Found '+str(Nsubcubes)+' coordinate sets to generate sub cubes for '

    if (type(dras) == float) & (type(ddecs) == float):
        dras  = [dras] * Nsubcubes
        ddecs = [ddecs] * Nsubcubes

    cubehdr = pyfits.open(datacube)[cube_ext[0]].header
    wavevec = np.arange(cubehdr['NAXIS3'])*cubehdr['CD3_3']+cubehdr['CRVAL3']

    for ii in xrange(Nsubcubes):
        subcubestr  = (str(dras[ii])+'x'+str(ddecs[ii]) ).replace('.','p')
        subcubename = outputdir+'subcube_iiiiii_'+subcubestr+'arcsec.fits'
        if names is not None:
            subcubename = subcubename.replace('iiiiii',names[ii])
        else:
            subcubename = subcubename.replace('iiiiii',str("%.5d" % (ii+1)))

        tu.extract_subcube(datacube,ras[ii],decs[ii],[dras[ii],ddecs[ii]],subcubename,cubeext=cube_ext,
                           clobber=clobber,imgfiles=None,imgexts=None,imgnames=None,verbose=verbose)

        subcube = pyfits.open(subcubename)['DATA_DCBGC'].data

        for ww, cwave in enumerate(wavecenters[ii]):
            dwave = dwaves[ii][ww]
            narrowbandname  = subcubename.replace('.fits','narrowbandimage_cwave'+str(cwave).replace('.','p')+
                                                  'dwave'+str(dwave).replace('.','p')+'.fits')

            wavemin = cwave-dwave
            wavemax = cwave+dwave

            goodent = np.where((wavevec < wavemax) & (wavevec > wavemin))[0]

            if len(goodent) >= 2:
                if verbose: print ' - Saving model cube to \n   '+narrowbandname
                narrowbandimage = np.sum(subcube[goodent,:,:],axis=0)

                imghdr = tu.strip_header(pyfits.open(subcubename)[cube_ext[0]].header.copy())
                for key in imghdr.keys():
                    if '3' in key:
                        del imghdr[key]
                    if 'ZAP' in key:
                        del imghdr[key]

                hduimg   = pyfits.PrimaryHDU(narrowbandimage,header=imghdr)
                hdus     = [hduimg]
                hdulist  = pyfits.HDUList(hdus)                  # turn header into to hdulist
                hdulist.writeto(narrowbandname,clobber=clobber)  # write fits file (clobber=True overwrites excisting file)
            else:
                if verbose: print ' - WARNING: less than 2 slices in narrowband extraction from subcube trying to generate'
                if verbose: print '   '+narrowbandname

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =