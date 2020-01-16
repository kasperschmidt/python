# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                           Utilities for MUSE-Wide related stuff
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#import pyfits     # KBS180118 -> import astropy.io.fits as pyfits
import astropy.io.fits as afits
import numpy as np
import MUSEWideUtilities as mwu
import sys
import stacking
import scipy.ndimage
import tdose_utilities as tu
import kbsutilities as kbs
from astropy import wcs
from astropy.coordinates import SkyCoord
import datetime
import felis
import glob
import pdb
import collections
import os
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_pixelpos(ra,dec,pointingname,imgdir='/Users/kschmidt/work/images_MAST/MUSEWidePointings/',imgext=0,
                 radecunit="deg",ignorestr='_wht_',pixorigin=1,verbose=True):
    """
    Get pixel positions in a muse cube for

    --- INPUT ---
    ra                  Right ascension of coordinate set to convert.
    dec                 Declination of coordinate set to convert.
    pointingname        Name of pointing (or similar) to combine with imgdir when searching for images
    imgdir              directory to serach for images
    imgext              Image extension to use WCS from
    radecunit           Units of ra and dec
    ignorestr           Ignore images contianing the provided string when searching through imgdir
    pixorigin           Pixel orgigin when coinvertin to pixcels
    verbose             Toggle verbosity

    """
    searchstr = imgdir+'*'+pointingname+'*.fits'
    img_in = glob.glob(searchstr)

    img = []
    for imname in img_in:
        if ignorestr in imname.split('/')[-1]:
            continue
        else:
            img.append(imname)

    if len(img) == 1:
        imghdr = afits.open(img[0])[imgext].header
    else:
        sys.exit(' Found '+str(len(img))+' images looking for '+searchstr)

    if verbose: print(' Getting pixel position in '+img[0])
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
        pointingname = 'candels-cdfs-'
    elif idstr[0] == '2':
        pointingname = 'candels-cosmos-'
    elif idstr[0] == '3':
        pointingname = 'hudf09-1-'
    elif idstr[0] == '4':
        pointingname = 'hudf09-2-'
    elif idstr[0] == '5':
        pointingname = 'udf-'
    elif len(idstr) <= 8:
        pdb.set_trace('Length of ID string of non-UDF-mosaic objects should be longer than 8;'
                      ' stopping in MUSEWideUtilities.gen_pointingname()...')
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
    if verbose: print(' - Loading inverse variance wheight map')
    wht_data            = afits.open(wht_file)[wht_ext].data
    wht_data_shape      = wht_data.shape

    if wht_units == 'ivar':
        if verbose: print(' - Weight image is in units "ivar", i.e., inverser variamce. Turning map into sigmas')
        sigma_img = np.sqrt(1.0/wht_data)
    elif (wht_units == 'sigma') or (wht_units == 'stddev'):
        if verbose: print(' - Weight image already in sigma units "sigma" or "stddev". No conversion applied')
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
    if verbose: print(' - Updating fits header with conversion information ')
    sigma_hdr  = afits.open(wht_file)[wht_ext].header
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
    if verbose: print(' - Saving image to \n   '+filename)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if 'XTENSION' in header.keys():
        hduprim        = afits.PrimaryHDU()  # default HDU with default minimal header
        hducube        = afits.ImageHDU(imagedata,header=header)
        hdus           = [hduprim,hducube]
    else:
        hducube = afits.PrimaryHDU(imagedata,header=header)
        hdus           = [hducube]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hdulist = afits.HDUList(hdus)             # turn header into to hdulist
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
    if verbose: print(' - Will store the updated infofile in \n   '+outfile)
    fout       = open(outfile,'w')

    PSFinfo    = afits.open('/Users/kschmidt/work/catalogs/MUSE_GTO/psf_all_Converted_cleaned.fits')[1].data

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
    print('ERROR: Cannot be run within ipython using import. Run by simply typing "python determine_bb.py'),pdb.set_trace()
    print(' - Assuming positioned in /Users/kschmidt/work/MUSE/Josie_DoublePeakInspection/')
    print(' - Importing determine_bb.py (remember to edit file and catalog info ')
    import determine_bb
    print(' - Launching inspection GUI; happy inspecting... ')
    determine_bb.main()

    print(' - GUI existed after inspections')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def create_manualimgFromCube(datacube,ras,decs,dras,ddecs,wavecenters,dwaves,outputdir,
                             cube_ext=['DATA_DCBGC','EFF_STAT'],skip_cube_extraction=False,
                             invertselection=False,
                             names=None,overwrite=False,verbose=True):
    """
    Create a manually defined image by removing specified layers of a data cube.
    This is a compliment to mu.create_narrowband_subcube() which only handles a continuous set of wavelength layers.
    However, for cutting out line emission from data cubes, it can be useful to be able to remove multiple
    pre-defined wavelength ranges. This is what this routine is for.

    --- INPUT ---
    datacube                Data cube to extract narrow band images from
    ras                     List of R.A.s to extract narrow band images around
    decs                    List of Dec.s to extract narrow band images around
    dras                    Size of narrow band images in R.A. direction (provided in arcsec)
    ddecs                   Size of narrow band images in Dec. direction (provided in arcsec)
    wavecenters             List of lists of central wavelengths of narrow bands to remove from data cube.
                            Use format:
                                [[obj1center1,obj1center2,...,obj1centerN],[obj2center1,obj2center2,...,obj2centerN],...]
    dwaves                  The width of the band to remove from the data cube given in Angstrom.
                            Use format:
                                [[obj1dwave1,obj1dwave2,...,obj1dwaveN],[obj2dwave1,obj2dwave2,...,obj2dwaveN],...]
    outputdir               Directory to store outputs in
    cube_ext                list of extensions to extract from cube into sub-cube. The header of the first extension mentioned
                            in the list will be used as WCS and Header template for the narrow-band image generated
    names                   Names or IDs of locations/objects indicated by the ras and decs input
    overwrite               Overwrite existing files
    skip_cube_extraction    If sub-cubes have already been extracted, re-generating them can be skipped
                            setting this keywoprd to True
    invertselection         If True, images of an inverted selection (i.e., keeping instead of removing the defined
                            wavelength ranges) will be generated in addition to the original picture.
                            For instance, if emission lines are removed to define a continuum white light image, an
                            emission line image will be genrated if invertselection = True.
    verbose                 Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mu

    # remove Lya emission lines and look for continuum
    pointing = 'candels-cdfs-01'
    outdir   = '/Users/kschmidt/work/TDOSE/190215_Test_create_manualimgFromCube/'
    ras      = [53.0612466535,53.0540222689]
    decs     = [-27.8047374553,-27.8037275129]
    wcenters = [[5830.0,7426.5,9152.5],[5869.5,7476.9,11000]]
    dwaves   = [[10,10,10],[10,10,10]]
    names    = ['101011026','101012027']
    datacube = glob.glob('/Volumes/DATABCKUP1/MUSE-Wide/DATACUBES/DATACUBE_'+pointing+'_v1.0_dcbgc_effnoised.fits')

    mu.create_manualimgFromCube(datacube[0],ras,decs,5.0,5.0,wcenters,dwaves,outdir,names=names,overwrite=True)

    # Generate image of continuum (removing emission lines) image of OII and OIII lines
    pointing = 'candels-cdfs-03'
    outdir   = '/Users/kschmidt/work/TDOSE/190215_Test_create_manualimgFromCube/'
    ras      = [53.0697905498543,53.0697905498543]
    decs     = [-27.8390839474992,-27.8390839474992]
    zobj     = 0.560494468420294
    linepos  = np.array([3726.0,3729.0,4861.33,4959.0,10000.0])*(1.0+zobj)
    dwave    = 500.0/299792.0 * linepos  # 500km/s x 2 to give width of 1000km/s bands
    wcenters = [linepos,[5300,6777,8564]]
    dwaves   = [dwave,[500,947,736]]
    names    = ['103013064_cont','103013064_oxygenlines']
    datacube = glob.glob('/Volumes/DATABCKUP1/MUSE-Wide/DATACUBES/DATACUBE_'+pointing+'_v1.0_dcbgc_effnoised.fits')

    mu.create_manualimgFromCube(datacube[0],ras,decs,15.0,15.0,wcenters,dwaves,outdir,names=names,overwrite=True,invertselection=True,skip_cube_extraction=True)

    """
    Nsubcubes = len(ras)
    if len(decs) != Nsubcubes:
        sys.exit(' The number of RAs does not match the number of Decs. Make sure they do ')
    if verbose: print(' - Found '+str(Nsubcubes)+' coordinate sets to generate sub cubes for ')

    if (type(dras) == float) & (type(ddecs) == float):
        dras  = [dras] * Nsubcubes
        ddecs = [ddecs] * Nsubcubes

    cubehdr = afits.open(datacube)[cube_ext[0]].header
    wavevec = np.arange(cubehdr['NAXIS3'])*cubehdr['CD3_3']+cubehdr['CRVAL3']

    if verbose: print(' - Looping of (ra,dec) coordinate pairs in data cube to generate images for ')
    for ii, ra in enumerate(ras):
        subcubestr  = (str(dras[ii])+'x'+str(ddecs[ii]) ).replace('.','p')
        subcubename = outputdir+'subcube_iiiiii_'+subcubestr+'arcsec.fits'
        if names is not None:
            subcubename = subcubename.replace('iiiiii',names[ii])
        else:
            subcubename = subcubename.replace('iiiiii',str("%.5d" % (ii+1)))

        if os.path.isfile(subcubename) & (overwrite == False):
            if verbose: print(' - '+subcubename+' exists and overwrwite=False so skipping...')
        else:
            if not skip_cube_extraction:
                tu.extract_subcube(datacube,ras[ii],decs[ii],[dras[ii],ddecs[ii]],subcubename,cubeext=cube_ext,
                                   clobber=overwrite,imgfiles=None,imgexts=None,imgnames=None,verbose=verbose)

            for cc, cext in enumerate(cube_ext):
                subcube = afits.open(subcubename)[cext].data
                narrowbandname  = subcubename.replace('.fits','_narrowbandimage_'+cext+'_manualimgFromCube.fits')
                if verbose: print(' - Found '+str(len(wavecenters[ii]))+' regions to remove from cube before collapsing ')
                if verbose: print(' - Defining layers to extract based on ranges to '
                                  'exclude provide via "wavecenters" and "dwaves" keywords.')
                badwaveent_all = np.array([])
                for ww, cwave in enumerate(wavecenters[ii]):
                    dwave      = dwaves[ii][ww]
                    wavemin    = cwave-dwave
                    wavemax    = cwave+dwave
                    badwaveent     = np.where((wavevec < wavemax) & (wavevec > wavemin))[0]
                    if len(badwaveent) == 0:
                        if verbose: print('\n - WARNING: wavelength range '+str(cwave)+'+/-'+str(dwave)+'A is outside cube\n')
                    badwaveent_all = np.append(badwaveent_all,badwaveent)

                badwaveent_all = np.unique(np.sort(badwaveent_all))
                entlist        = np.arange(len(wavevec))
                goodent        = np.delete(entlist,badwaveent_all.astype(int))

                if len(goodent) >= 2:
                    if verbose: print(' - Saving narrowband image to \n   '+narrowbandname)
                    narrowbandimage = np.sum(subcube[goodent,:,:],axis=0)

                if len(goodent) >= 2:
                    cubehdr = tu.strip_header(afits.open(subcubename)[cube_ext[0]].header.copy())

                    if invertselection:
                        mwu.collapsecube(narrowbandname,subcube,cubehdr,layers=badwaveent_all.astype(int),
                                         overwrite=overwrite,verbose=verbose)
                    else:
                        mwu.collapsecube(narrowbandname,subcube,cubehdr,layers=goodent.astype(int),
                                         overwrite=overwrite,verbose=verbose)
                else:
                    if verbose: print(' - WARNING: less than 2 slices in narrowband extraction from subcube '+cext+' in')
                    if verbose: print('   '+narrowbandname)
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
    if verbose: print(' - Found '+str(Nsubcubes)+' coordinate sets to generate sub cubes for ')

    if (type(dras) == float) & (type(ddecs) == float):
        dras  = [dras] * Nsubcubes
        ddecs = [ddecs] * Nsubcubes

    cubehdr = afits.open(datacube)[cube_ext[0]].header
    wavevec = np.arange(cubehdr['NAXIS3'])*cubehdr['CD3_3']+cubehdr['CRVAL3']

    for ii in xrange(Nsubcubes):
        subcubestr  = (str(dras[ii])+'x'+str(ddecs[ii]) ).replace('.','p')
        subcubename = outputdir+'subcube_iiiiii_'+subcubestr+'arcsec.fits'
        if names is not None:
            subcubename = subcubename.replace('iiiiii',names[ii])
        else:
            subcubename = subcubename.replace('iiiiii',str("%.5d" % (ii+1)))

        if os.path.isfile(subcubename) & (clobber == False):
            if verbose: print(' - '+subcubename+' exists and clobber=False so skipping...')
        else:
            tu.extract_subcube(datacube,ras[ii],decs[ii],[dras[ii],ddecs[ii]],subcubename,cubeext=cube_ext,
                               clobber=clobber,imgfiles=None,imgexts=None,imgnames=None,verbose=verbose)

            subcube = afits.open(subcubename)[cube_ext[0]].data

            for ww, cwave in enumerate(wavecenters[ii]):
                dwave = dwaves[ii][ww]
                narrowbandname  = subcubename.replace('.fits','_narrowbandimage_cwave'+str(cwave).replace('.','p')+
                                                      'dwave'+str(dwave).replace('.','p')+'.fits')

                wavemin = cwave-dwave
                wavemax = cwave+dwave

                goodent = np.where((wavevec < wavemax) & (wavevec > wavemin))[0]

                if len(goodent) >= 2:
                    cubehdr = tu.strip_header(afits.open(subcubename)[cube_ext[0]].header.copy())
                    mwu.collapsecube(narrowbandname,subcube,cubehdr,layers=goodent,overwrite=clobber,verbose=verbose)
                else:
                    if verbose: print(' - WARNING: less than 2 slices in narrowband extraction from subcube trying to generate')
                    if verbose: print('   '+narrowbandname)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def collapsecube(outname,dataarray,fitsheader,layers='all',normalize=False,overwrite=False,verbose=True):
    """
    Function to callapse specified layers of data cube and store them as image

    --- INPUT ---
    outname    Name of image fits file to generate
    dataarray       The data cube to collapse. Will be collapsed along the 0th axis
    cubeheader      Fits header of data cube (in whihc case all third dimension coomponents are removed) or
                    image to use when storing the output image to a fits file.
    layers          The layers of the cube to collaps. If layers='all' the full cube is collapsed.
    normalize       Normalize output by the number of layers collapsed
    overwrite       Overwrite existing image?
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu
    mwu.collapsecube(narrowbandname,subcube,cubehdr,layers=goodent,verbose=verbose)
    """
    if verbose: print(' - Saving narrowband image to \n   '+outname)

    if layers == 'all':
        layers = np.arange(dataarray.shape[0])

    dataarray[np.where(~np.isfinite(dataarray))] = 0.0
    narrowbandimage = np.sum(dataarray[layers,:,:],axis=0)

    for key in fitsheader.keys():
        if '3' in key:
            try:
                del fitsheader[key]
            except:
                pass
        if 'ZAP' in key:
            try:
                del fitsheader[key]
            except:
                pass

    hduimg   = afits.PrimaryHDU(narrowbandimage,header=fitsheader)
    hdus     = [hduimg]
    hdulist  = afits.HDUList(hdus)                  # turn header into to hdulist
    hdulist.writeto(outname,overwrite=overwrite)  # write fits file (clobber=True overwrites excisting file)

    if normalize:
        imghdu = afits.open(outname, mode='update')
        imghdu[0].data = imghdu[0].data / len(layers)
        imghdu.flush()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_PSFasciiselectiontemplate(outname='./PSFselectionAsciiTemplate_RENAME_.txt',campaign='E40',clobber=False,verbose=True):
    """

    --- INPUT ---
    outname    Name of ascii template to generate
    campaign   Defines what inspection campaign to add fields names to ascii template for
    clobber    Overwrite existing file?
    verbose    Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu
    campaigns = ['E24','E36','E40']
    for campaign in campaigns:
        output    = '/Users/kschmidt/work/MUSE/MUSEWide_PSFs/PSFselectionAsciiTemplate_'+campaign+'_TEMPLATE.txt'
        mwu.gen_PSFasciiselectiontemplate(outname=output,campaign=campaign,clobber=True)
    """

    if os.path.isfile(outname) & (clobber == False):
        sys.exit(' The output '+outname+' exists and clobber=False ')
    else:
        if verbose: print(' - Getting list of fields in campaign = '+campaign)
        if campaign == 'E24':
            fields = mwu.E24fields()
        elif campaign == 'E36':
            fields = mwu.E36fields()
        elif campaign == 'E40':
            fields = mwu.E40fields()
        else:
            sys.exit(' Fields still need to be set up for campaign = '+campaign)

        fout = open(outname,'w')
        fout.write('# Template for selection of PSF parameter generated on '+
                   datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+
                   ' with MUSEWideUtilities.gen_PSFasciiselectiontemplate()\n')
        fout.write('# field                PSF_FullConv         PSF_PampelMuse        PSF_GalaxyFitting  # \n')

        if verbose: print(' - Writing template to '+outname)
        for field in fields:
            fout.write(str("%15s" %field)+'              0                     0                       0         # \n')

        fout.close()
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def check_PSFasciiselectiontemplate(selectionascii):
    """
    Check that 1 (and exactly 1) option is chosen for each of the fields in a PSF selection file

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu
    data = mwu.check_PSFasciiselectiontemplate('/Users/kschmidt/work/MUSE/MUSEWide_PSFs/PSFselectionAsciiTemplate_E40_TEST.txt')

    """
    data = np.genfromtxt(selectionascii,dtype=None,names=True,skip_header=1)
    cols = data.dtype.names

    for ff, field in enumerate(data['field']):
        selectionvec = np.asarray([data[col][ff] for col in cols[1:]])
        Nsel = int(np.sum(selectionvec))
        if Nsel == 0:
            print(' - WARNING: Field '+field+' has NO selected PSF')
        elif Nsel > 1:
            print(' - WARNING: Field '+field+' has '+str(Nsel)+' selected PSFs')
        else:
            continue

    return data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_masterPSFcat(fitscatalogs,outcatname,psfselections=None,clobber=False,verbose=True):
    """

    --- INPUT ---
    fitscatalogs        Fits catalogs with PSF parameters from campaigns (E24, E36 and E40 for instance) to combine
    psfselections       PSF selection ascii files corresponding to the fits catalogs. If None values are set to -99
                        for the selection columns in output catalog
    outcatname          The name of the output catalog to generate
    clobber             Overwrite existing output?
    verbose             Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu
    import glob
    fitscatalogs  = ['/Users/kschmidt/work/catalogs/MUSE_GTO/psf_all_Converted_cleaned.fits','/Users/kschmidt/work/MUSE/MUSEWide_PSFs/psf_E40_final.fits']
    outcatname    = '/Users/kschmidt/work/MUSE/MUSEWide_PSFs/PSFmastercatalogTEST.fits'
    psfselections = glob.glob('/Users/kschmidt/work/MUSE/MUSEWide_PSFs/PSFselectionAsciiTemplate_E*_180112.txt')
    mwu.gen_masterPSFcat(fitscatalogs,outcatname,psfselections=psfselections,clobber=False)

    """
    import pyfits
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Gettting list of MUSE-Wide fields to look for')
    MWfields = mwu.E24fields() + mwu.E36fields() + mwu.E40fields()

    columns = ['field', 'AG',
               'p0_FFT_m', 'p1_FFT_m', 'beta_FFT', 'p0_FFT_g', 'p1_FFT_g',
               'p0_FFTS_m', 'p1_FFTS_m', 'beta_FFTS', 'p0_FFTS_g', 'p1_FFTS_g',
               'p0_comb_m', 'p1_comb_m', 'beta_comb', 'p0_comb_g', 'p1_comb_g',
               'p0_PM_m', 'p1_PM_m', 'beta_PM', 'p0_PM_g', 'p1_PM_g',
               'p0_S_g', 'p1_S_g',
               'psf_sel_m','p0_sel_m','p1_sel_m','beta_sel',
               'psf_sel_g','p0_sel_g','p1_sel_g']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    psfseldatadic = {}
    if psfselections is not None:
        if verbose: print(' - Loading data from PSF selection files')
        for psfselfile in psfselections:
            seldata = mwu.check_PSFasciiselectiontemplate(psfselfile)
            psfseldatadic[psfselfile] = seldata

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Extract information from fitscatalogs: ')
    for cc, cat in enumerate(fitscatalogs):
        if verbose: print('   -> '+cat)
        data     = pyfits.open(cat)[1].data
        datacols = data.dtype.names

        for ff, fieldname in enumerate(data['field']):
            fieldname = fieldname.split()[0] #removing trailing spaces

            if ('psf_all_Converted_cleaned.fits' in cat) & \
                    (fieldname in ['cdfs-47','cdfs-48','cdfs-49','cdfs-51','cdfs-52']):
                continue

            if fieldname in MWfields:
                fieldentry = [0]*len(columns)
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                for oo, col in enumerate(columns):
                    val = -99
                    if col in datacols:
                        val = data[col][ff]
                    fieldentry[oo] = val

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Filling entries for Gauss PSF selection
                if psfselections is not None:
                    for psfselfile in psfselections:
                        data_psfsel = psfseldatadic[psfselfile]
                        fieldent    = np.where(data_psfsel['field'] == fieldname)[0]

                        if len(fieldent) == 1:
                            for selcol in ['PSF_FullConv', 'PSF_PampelMuse', 'PSF_GalaxyFitting', 'PSF_Sextractor']:
                                try:
                                    selval = data_psfsel[selcol][fieldent]
                                except:
                                    selval = 0

                                if selval == 1:
                                    if selcol == 'PSF_FullConv':
                                        selstring = 'FFT'
                                    elif selcol == 'PSF_PampelMuse':
                                        selstring = 'PM'
                                    elif selcol == 'PSF_GalaxyFitting':
                                        selstring = 'comb'
                                    elif selcol == 'PSF_Sextractor':
                                        selstring = 'S'

                                    fieldentry[-3] = selstring
                                    fieldentry[-2] = data['p0_'+selstring+'_g'][ff]# 'p0_sel_g'
                                    fieldentry[-1] = data['p1_'+selstring+'_g'][ff]# 'p1_sel_g'
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Filling entries for Moffat PSF selection
                if ('psf_all_Converted_cleaned.fits' in cat) & ('cosmos' not in fieldentry[0]):
                    if fieldentry[17] != 0: # Checking if PampelMUSE
                        fieldentry[-7] = 'PM'
                        fieldentry[-6] = data['p0_PM_m'][ff]# 'p0_sel_m'
                        fieldentry[-5] = data['p1_PM_m'][ff]# 'p1_sel_m'
                        fieldentry[-4] = data['beta_PM'][ff]# 'beta_sel'
                    else:
                        fieldentry[-7] = 'FFT'
                        fieldentry[-6] = data['p0_FFT_m'][ff]# 'p0_sel_m'
                        fieldentry[-5] = data['p1_FFT_m'][ff]# 'p1_sel_m'
                        fieldentry[-4] = data['beta_FFT'][ff]# 'beta_sel'

                elif 'psf_E40_final.fits' in cat:
                    #if fieldentry[17] != 0: # Checking if PampelMUSE
                    fieldentry[-7] = 'TBD'
                    fieldentry[-6] = -99
                    fieldentry[-5] = -99
                    fieldentry[-4] = -99
                else:
                    #print(' - No PampelMUSE fits to check for in catalog '+cat)
                    fieldentry[-7] = 'TBD'
                    fieldentry[-6] = -99
                    fieldentry[-5] = -99
                    fieldentry[-4] = -99

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                if fieldname == 'cdfs-28':
                    fieldentry[-3] = 'MeanSlopeAndAGfix'
                    fieldentry[-2] = 0.963# 'p0_sel_g'
                    fieldentry[-1] = -4.4394e-05# 'p1_sel_g'
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                if fieldname == 'cdfs-37':
                    fieldentry[-7] = 'FFTS'# 'psf_sel_m'
                    fieldentry[-6] = data['p0_FFTS_m'][ff]# 'p0_sel_m'
                    fieldentry[-5] = data['p1_FFTS_m'][ff]# 'p1_sel_m'
                    fieldentry[-4] = data['beta_FFTS'][ff]# 'beta_sel'

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                for vv, val in enumerate(fieldentry):
                    if val == 0.0:
                        fieldentry[vv] = -99

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                try:
                    masterarray = np.vstack([masterarray, np.asarray(fieldentry)])
                except:
                    masterarray = np.asarray(fieldentry)

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            else:
                print(' WARNING The field '+fieldname+' from '+cat+' is not in the MUSE-Wide field list! ')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Check for duplicate entries in catalog')
    arrayfields = [field.split()[0] for field in masterarray[:,0]]
    ufields     = np.unique(arrayfields)
    if len(arrayfields) != len(ufields):
        for uf in ufields:
            if len(np.where(np.asarray(arrayfields) == uf)[0]) > 1:
                if verbose: print('   WARNING Duplicate entry for '+uf)
    else:
        if verbose: print('   Good - only one entry per field.')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Writing assembled array to fitsfile')
    columndic = collections.OrderedDict()

    for cc, col in enumerate(columns):
        if col in ['field','psf_sel_g','psf_sel_m']:
            columndic[col]    =  pyfits.Column(name=col, format='A20', unit='', array=masterarray[:,cc].astype(str))
        elif 'beta' in col:
            columndic[col]    =  pyfits.Column(name=col+'_m', format='D', unit='', array=masterarray[:,cc].astype(float))
        else:
            columndic[col]    =  pyfits.Column(name=col, format='D', unit='', array=masterarray[:,cc].astype(float))

    collist   = [columndic[key] for key in columndic.keys()]
    coldefs   = pyfits.ColDefs(collist)
    th        = pyfits.new_table(coldefs) # creating default header

    # writing hdrkeys:'---KEY--',                             '----------------MAX LENGTH COMMENT-------------'
    #th.header.append(('MAG     ' , spec2D[0].header['MAG']   ,'MAG_AUTO from interlaced catalog'),end=True)

    head    = th.header
    tbHDU   = pyfits.new_table(coldefs, header=head)

    tbHDU.writeto(outcatname, clobber=clobber)
    if verbose: print('   Fits table stored in \n   '+outcatname)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def E24fields():
    """
    Returning list of fields in E40
    """
    fields = ['cdfs-01', 'cdfs-02', 'cdfs-03', 'cdfs-04', 'cdfs-05', 'cdfs-06',
               'cdfs-07', 'cdfs-08', 'cdfs-09', 'cdfs-10', 'cdfs-11', 'cdfs-12',
               'cdfs-13', 'cdfs-14', 'cdfs-15', 'cdfs-16', 'cdfs-17', 'cdfs-18',
               'cdfs-19', 'cdfs-20', 'cdfs-21', 'cdfs-22', 'cdfs-23', 'cdfs-24']
    return fields
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def E36fields():
    """
    Returning list of fields in E40
    """
    fields = ['cdfs-25', 'cdfs-26', 'cdfs-28', 'cdfs-29', 'cdfs-30', 'cdfs-31', 'cdfs-32',
              'cdfs-33', 'cdfs-34', 'cdfs-35', 'cdfs-36', 'cdfs-37', 'cdfs-39', 'cdfs-40',
              'cdfs-41', 'cdfs-42', 'cdfs-43', 'cdfs-44', 'cdfs-45', 'cdfs-46',
              'cosmos-01', 'cosmos-02', 'cosmos-03', 'cosmos-04', 'cosmos-05', 'cosmos-06',
              'cosmos-07', 'cosmos-08', 'cosmos-09', 'cosmos-10', 'cosmos-11', 'cosmos-12',
              'cosmos-13', 'cosmos-14', 'cosmos-15', 'cosmos-16']

    return fields
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def E40fields():
    """
    Returning list of fields in E40
    """
    fields = ['udf-01', 'udf-02', 'udf-03', 'udf-04', 'udf-05',
              'udf-06', 'udf-07', 'udf-08', 'udf-09',
              'cdfs-47', 'cdfs-48', 'cdfs-49', 'cdfs-50', 'cdfs-51', 'cdfs-52',
              'cdfs-53', 'cdfs-54', 'cdfs-55', 'cdfs-56', 'cdfs-57', 'cdfs-58',
              'cdfs-59', 'cdfs-60', 'cdfs-61', 'cdfs-62',
              'cosmos-17', 'cosmos-18', 'cosmos-19', 'cosmos-20', 'cosmos-21', 'cosmos-58', 'cosmos-59',
              'hudf09-1-01', 'hudf09-1-02', 'hudf09-1-03', 'hudf09-1-04',
              'hudf09-2-01', 'hudf09-2-02', 'hudf09-2-03', 'hudf09-2-04']
    return fields
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def match_templates2JosiesStack(spec,outfile, templatedir='/Users/kschmidt/work/MUSE/uvEmissionlineSearch/felis_templates/',
                                verbose=True):
    """
    Wrapper around felis.match_templates2specs() to template to the potential lines in Josie's stack of LyC candidates

    --- EXAMPLE OF USE ---

    spec    = '/Users/kschmidt/work/MUSE/KeruttLyCstack/comb_median_all_zcomb_with_cands.fits'
    outfile = spec.replace('.fits','_tdoseformat.fits')

    import MUSEWideUtilities as mwu
    mwu.match_templates2JosiesStack(spec,outfile)


    """

    JKspec  = afits.open(spec)[1].data
    wave    = JKspec['Wavelength']
    flux    = JKspec['Spectrum']
    fluxerr = JKspec['Error']

    print(' - Saving spectrum in the FELIS/TDOSE format: \n   '+outfile)
    felis.save_spectrum(outfile,wave,flux,fluxerr,headerinfo=None,overwrite=True,verbose=verbose)
    datestr = '180907'

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Grabbing templates from '+templatedir)
    waverests   = {'CIII':1908.0, 'CIV':1549.0, 'HEII':1640.420}


    for line in ['CIII','CIV','HEII']:
        temps    = glob.glob(templatedir+'uves_felis_template_'+line+'doublet_sig_*_flux'+line+'1_1p0_flux*.fits')

        plotdir  = '/Users/kschmidt/work/MUSE/KeruttLyCstack/'
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print(' - Crosscorrelating templates to spectra using FELIS')
        picklefile = '/Users/kschmidt/work/MUSE/KeruttLyCstack/CCresults'+datestr+'_'+line+'_RENAME_.pkl'

        ccdic      = felis.match_templates2specs(temps,[outfile],[0.0],picklefile,wavewindow=[15],plotdir=plotdir,
                                                 wavecen_restframe=[waverests[line]],vshift=[0],min_template_level=1e-4,
                                                 plot_allCCresults=True,subtract_spec_median=True)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def build_galfit_wrapper_setups(outdir,build='DR1objWOskeltonID'):
    """
    Building setup files that can be loaded by /Users/kschmidt/work/galfit_wrapper/galfit_wrapper_obj.py

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu
    outdir    = '/Users/kschmidt/work/MUSE/MWDR1galfitmodeling/galfit_wrapper_setups/'
    mwu.build_galfit_wrapper_setups(outdir,build='DR1objWOskeltonID')

    import MUSEWideUtilities as mwu
    outdir    = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/galfit_wrapper_setups/'
    mwu.build_galfit_wrapper_setups(outdir,build='OIIemitters')

    """
    if build == 'DR1objWOskeltonID':
        maindir   = '/Users/kschmidt/work/MUSE/MWDR1galfitmodeling/'
        imagedir  = '/Users/kschmidt/work/images_MAST/MUSEWidePointings/'
        sigmadir  = '/Users/kschmidt/work/MUSE/MWDR1galfitmodeling/hst_galfit_sigma_images/'
        pstampdir = 'PPPPPP'

        cat     = 'MW_44fields_emission_woskeltoncp_no_LAEs.fits'
        catpath = maindir+cat
        fitsdat = afits.open(catpath)[1].data

        areas  = [fn.split('-')[0] for fn in np.unique(fitsdat['field_name'])]
        fields = [fn.replace(fn.split('-')[0]+'-','') for fn in np.unique(fitsdat['field_name'])]
        #bands  = ['acs_775w','acs_814w']
        bands  = ['acs_435w', 'acs_606w', 'acs_775w', 'acs_814w', 'wfc3_105w', 'wfc3_125w', 'wfc3_160w'] # no PSF for 'acs_850lp'

        for ff, field in enumerate(fields):
            for band in bands:
                setupfile = outdir+'galfit_wrapper_setup_'+field+'_'+band+'.txt'

                if os.path.isfile(setupfile):
                    sys.exit(setupfile+' already exists')

                inputimage    = imagedir+band+'_'+areas[ff]+'-'+field+'_cut_v1.0.fits'
                sigmaimage    = sigmadir+band+'_'+areas[ff]+'-'+field+'_wht_cut_sigma_smooth.fits'
                outlogdir     = 'galfit_wrapper_output_logs/'
                outputdir     = 'galfit_wrapper_output_fits/'
                outputimg     = outputdir+'imgblock_'+band # .fits automatically appended
                outputcat     = 'galfit_wrapper_output_cat_'+areas[ff]+'-'+field+'_'+band+'_'+build+'.fits'
                resultsdir    = 'galfit_wrapper_results/'
                psfimage      = 'PSF_models/'+band+'/imgblock_6475.fits[2]'
                resultsfinal  = 'galfit_wrapper_results_final/' #+areas[ff]+'-'+band+'_parameters_with_sky/'
                inputparamdir = 'galfit_wrapper_input_params/'
                inputparams   = inputparamdir+'params_'+band

                if band == 'acs_814w':
                    shapeparam    = 'None'
                else:
                    shapeparam    = maindir+resultsfinal+outputcat.replace(band,'acs_814w')

                fout = open(setupfile,'w')
                fout.write("""# Setup file to load in galfit_wrapper_obj.py to overwrite hard-coded parameters
# name  parameter
use_for                    %s
field                      %s
band                       %s
catalogue                  %s
area                       %s
input_catalogue            %s
input_image                %s
sigma_image                %s
out_log_dir                %s
image_dir                  %s
output_image               %s
output_catalogue           %s
Guo_catalogue              /Users/kschmidt/work/galfit_wrapper/examples/CANDELS.GOODSS.F160W.v1_1.photom.cat
counterparts_catalogue     None
counterparts_imgs          None # XXX/work1/josie/galfit/Tanya_cutouts/postage/
counterparts_suff          None # XXX_cutout.jpg
galfit_path                /usr/local/bin/galfit
result_folder              %s
psf_image                  %s
psf_sampling               1
bad_pixel                  None
param_constr               None
save_stuff                 True
results_final              %s
file_comments              None
use_shape_params           %s
placeholder_bombed         /Users/kschmidt/work/galfit_wrapper/examples/empty.fits
input_params_dir           %s
input_params               %s
                """ % ('obj',field,band,cat,areas[ff],catpath,inputimage,sigmaimage,
                       outlogdir,outputdir,outputimg,outputcat,resultsdir,
                       psfimage,resultsfinal,shapeparam,inputparamdir,inputparams))
                fout.close()

    elif build == 'OIIemitters':
        maindir      = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/'
        imagedir     = '/Users/kschmidt/work/images_MAST/MUSEWidePointings/'
        sigmadir     = '/Users/kschmidt/work/MUSE/MWDR1galfitmodeling/hst_galfit_sigma_images/'
        pstampdir    = '/Users/kschmidt/work/TDOSE/OIIemitterGalfitModels/cutouts/'

        cat     = 'OIIemitter_selection_23LTm814LT24_SkelsepLT0p3_Nobj153.fits'
        catpath = maindir+cat
        fitsdat = afits.open(catpath)[1].data

        areas  = [fn.split('-')[0] for fn in np.unique(fitsdat['field_name'])]
        fields = [fn.replace(fn.split('-')[0]+'-','') for fn in np.unique(fitsdat['field_name'])]
        bands  = ['acs_814w']
        #bands  = ['acs_435w', 'acs_606w', 'acs_775w', 'acs_814w', 'wfc3_105w', 'wfc3_125w', 'wfc3_160w']

        for ff, field in enumerate(fields):
            for band in bands:
                setupfile = outdir+'galfit_wrapper_setup_'+field+'_'+band+'.txt'

                if os.path.isfile(setupfile):
                    sys.exit(setupfile+' already exists')

                inputimage    = imagedir+band+'_'+areas[ff]+'-'+field+'_cut_v1.0.fits'
                sigmaimage    = sigmadir+band+'_'+areas[ff]+'-'+field+'_wht_cut_sigma_smooth.fits'
                outlogdir     = 'galfit_wrapper_output_logs/'
                outputdir     = 'galfit_wrapper_output_fits/'
                outputimg     = outputdir+'imgblock_'+band # .fits automatically appended
                outputcat     = 'galfit_wrapper_output_cat_'+areas[ff]+'-'+field+'_'+band+'_'+build+'.fits'
                resultsdir    = 'galfit_wrapper_results/'
                psfimage      = 'PSF_models/'+band+'/imgblock_6475.fits[2]'
                resultsfinal  = 'galfit_wrapper_results_final/' #+areas[ff]+'-'+band+'_parameters_with_sky/'
                inputparamdir = 'galfit_wrapper_input_params/'
                inputparams   = inputparamdir+'params_'+band
                pstampsuffix  = band+'_'+areas[ff]+'-'+field+'_cut_v1p0_MUSEWideIII_4p0X4p0arcseccut.png'

                if band == 'acs_814w':
                    shapeparam    = 'None'
                else:
                    shapeparam    = maindir+resultsfinal+outputcat.replace(band,'acs_814w')

                fout = open(setupfile,'w')
                fout.write("""# Setup file to load in galfit_wrapper_obj.py to overwrite hard-coded parameters
# name  parameter
use_for                    %s
field                      %s
band                       %s
catalogue                  %s
area                       %s
#
out_log_dir                %s     # empty dir to fill
image_dir                  %s     # empty dir to fill
result_folder              %s     # empty dir to fill
results_final              %s     # empty dir to fill
input_params_dir           %s     # empty dir to fill
#
output_image               %s
input_params               %s
output_catalogue           %s
input_catalogue            %s
#
input_image                %s
sigma_image                %s
Guo_catalogue              /Users/kschmidt/work/galfit_wrapper/examples/CANDELS.GOODSS.F160W.v1_1.photom.cat
psf_image                  %s  # THE SINNER FOR THE ABORT TRAP BREAK IF FULL PATH INCLUDED
placeholder_bombed         /Users/kschmidt/work/galfit_wrapper/examples/empty.fits
#
counterparts_catalogue     %s
counterparts_imgs          %s
counterparts_suff          %s
galfit_path                /usr/local/bin/galfit
psf_sampling               1
bad_pixel                  None
param_constr               None
save_stuff                 True
file_comments              None
use_shape_params           %s
                """ % ('obj',field,band,cat,areas[ff],
                       outlogdir,outputdir,resultsdir,resultsfinal,inputparamdir,
                       outputimg,inputparams,outputcat,catpath,
                       inputimage,sigmaimage,psfimage,
                       catpath,pstampdir,pstampsuffix,shapeparam))
                fout.close()
    elif build == 'JKexample':
        sys.exit('build option "'+build+'" not enabled yet')
        # write 'psf' (for fitting a star with a moffat) or 'obj'
        use_for = 'obj'

        field = 'cdfs-16' # the PSF star 6475 is in field 16
        band  = 'acs_435w'

        catalogue = 'MW_1-24_v3.1_LAEs.fits'
        #catalogue = 'e36_emline_master_v1.0_LAEs.fits'
        area = 'candels'
        #KBS input_catalogue = '/work1/josie/MUSE_WIDE/catalogs/'+catalogue
        input_catalogue = '/Users/kschmidt/work/galfit_wrapper/examples/'+catalogue
        #input_catalogue = 'PSF_objs.fits'

        # directory with the HST data
        #KBS input_image_pref = '/work1/josie/MUSE_WIDE/HST_images/cut/'+band+'_'+area+'-'
        input_image_pref = '/Users/kschmidt/work/galfit_wrapper/examples/'+band+'_'+area+'-'
        input_image_suff = '_cut.fits'
        input_image = input_image_pref+field+input_image_suff

        # KBS ->  /Users/kschmidt/work/images_MAST/MUSEWidePointings

        # write a file name or None, but it is always better to have sigma images!
        sigma_image = '/Users/kschmidt/work/galfit_wrapper/examples/acs_435w_candels-cdfs-16_wht_cut_sigma_smooth.fits'
        # sigma_image_pref = '/work1/josie/MUSE_WIDE/HST_images/sigma/smooth/'\
        #                    +band+'_'+area+'-'
        # if 'candels' in area and 'cdfs' in field or 'cosmos' in field:
        #     sigma_image_suff = '_wht_cut_sigma_smooth.fits'
        # else:
        #     sigma_image_suff = '_cut_wht_sigma_smooth.fits'
        # sigma_image = sigma_image_pref+field+sigma_image_suff

        output_image = 'output_fits/imgblock' # .fits will be added automatically

        if area == 'candels':
            #KBS Guo_catalogue = '/home/josie/MUSE_WIDE/Guo_catalogue/CANDELS.GOODSS.F160W.v1_1.photom.cat'
            Guo_catalogue = '/Users/kschmidt/work/galfit_wrapper/examples/CANDELS.GOODSS.F160W.v1_1.photom.cat'

        # catalogue with counterparts (by Tanya), or write None
        #KBS: counterparts_catalogue ='/work1/josie/galfit/Tanya_cutouts/counterparts2.fits'
        counterparts_catalogue ='/Users/kschmidt/work/galfit_wrapper/examples/counterparts2.fits'

        # cutouts of HST showing the counterparts (by Tanya)
        #KBS counterparts_imgs = '/work1/josie/galfit/Tanya_cutouts/postage/'
        counterparts_imgs = '/Users/kschmidt/work/galfit_wrapper/examples/'
        counterparts_suff = '_cutout.jpg'

        output_catalogue = 'cat_ident_'+area+'-'+field+'_rid_galfit_'+band+'.fits'

        # path to where you installed galfit
        #KBS galfit_path = '/work1/josie/Programs/./galfit'
        galfit_path = '/usr/local/bin/galfit'
        #KBS galfit_path = '/Users/kschmidt/work/galfit/galfit'

        # directory where you want to store your ouput images
        image_dir = 'output_fits/'
        # directory where you want to store the output log files
        out_log_dir = 'fit_logs/'
        # directory where you want to store the final output images
        result_folder = image_dir+'results_'+band+"/"

        # PSF image, optional, write a file name or None
        #KBS psf_image    = 'output_fits/PSF_models/'+band+'/imgblock_6475.fits[2]'#'psfJ.fits'
        psf_image    = '/Users/kschmidt/work/galfit_wrapper/examples/imgblock_6475.fits[2]'
        if '336' in band:
            psf_image = 'output_fits/PSF_models/'+band+'/imgblock_9836.fits[2]'
        psf_sampling = '1' # if no file name was given, this is ignored

        bad_pixel    = None # Bad pixel mask (FITS image or ASCII coord list)
        param_constr = None # File with parameter constraints (ASCII file)

        save_stuff = False

        # folder for all final results
        results_final = 'results/acs_814w_parameters_with_sky/'+band+'/'

        # file with comments, write None if there is none yet
        #KBS file_comments = 'comments_objects.txt'
        file_comments = 'comments_objectsKBS.txt'

        # if you want to use parameters from 814, enter directory, otherwise write None
        #KBS use_shape_params = 'results/acs_814w_parameters/acs_814w/cat_ident_candels-'+field+'_rid_galfit_acs_814w.fits'
        use_shape_params = None #'/Users/kschmidt/work/galfit_wrapper/examples/cat_ident_candels-'+field+'_rid_galfit_acs_814w.fits'

        # in case galfit crashed you need something to display
        #KBS placeholder_bombed = "empty.fits"
        placeholder_bombed = "/Users/kschmidt/work/galfit_wrapper/examples/empty.fits"

        # folder with input parameters
        input_params_dir = 'input_params/'
        input_params = input_params_dir+'params' # file name for input parameters
    else:
        sys.exit('Invalid build option "'+build+'"')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def stack_MUSEWideSpecs(stacktype='mean',stackobjects='TDOSEgalfitextractionsE60',runallspecs=False,
                        outfile = './MUSEWideStack_RENAME.fits', verbose=True):
    """
    Wrapper around stacking.stack_1D() to stack the MUSE-Wide spectra

    --- INPUT ---
    stacktype       keyword tellong stacking.stack_1D() how to stack
    stackobjects    keyword deciding what spectra to load and stack
    runallspecs     whether to run on all spectra, or just a small subset of selection (good for testing/setting up)
    verbose         toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu
    outfile = '/Users/kschmidt/work/MUSE/spectral_stacking/TDOSEgalfitextractionsE60_meanstack_181030.fits'
    mwu.stack_MUSEWideSpecs(stacktype='mean',stackobjects='TDOSEgalfitextractionsE60',outfile=outfile,runallspecs=True)

    """

    if verbose: print(' - Grabbing spectra for "'+stackobjects+'" to be stacked ')
    if verbose: print('   and building lists of wavelengths, fluxes and variances')
    if stackobjects == 'TDOSEgalfitextractionsE60':
        infofile  = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/LAEinfo_100fields_UDFshallow.fits'
        LAEinfo   = afits.open(infofile)[1].data
        z_Lya     = LAEinfo['redshift']
        z_sys_V18 = LAEinfo['z_sys_V18']
        id_all    = LAEinfo['id']
        specdir   = '/Volumes/DATABCKUP1/TDOSEextractions/180824_TDOSEextraction_LAEs60fields/modelimg/tdose_spectra/'
        if runallspecs:
            spectra  = glob.glob(specdir+'tdose_spectrum_candels-*.fits')
        else:
            # spectra  = [specdir+'tdose_spectrum_candels-cdfs-04_modelimg_0104014050-0104014050.fits',
            #             specdir+'tdose_spectrum_candels-cdfs-15_modelimg_0115003085-0115003085.fits',
            #             specdir+'tdose_spectrum_candels-cdfs-06_modelimg_0106004019-0106004019.fits',
            #             specdir+'tdose_spectrum_candels-cdfs-25_modelimg_0125042115-0125042115.fits']

            spectra  = [specdir+'tdose_spectrum_candels-cdfs-15_modelimg_0115003085-0115003085.fits']*100

        wavelengths  = []
        fluxes       = []
        variances    = []
        z_systemic   = []
        if verbose: print(' - Found '+str(len(spectra))+' spectra to stack')
        for spectrum in spectra:
            data        = afits.open(spectrum)[1].data
            wavelengths.append(data['wave'])
            fluxes.append(data['flux'])
            variances.append(data['fluxerror']**2.0)

            objid       = spectrum.split('-0')[-1].split('.fit')[0]
            objinfoent  = np.where(id_all == int(objid))[0]
            z_sys_obj   = z_sys_V18[objinfoent]
            if z_sys_obj == 0.0:
                z_sys_obj   = z_Lya[objinfoent]
            z_systemic  = np.append(z_systemic , z_sys_obj[0])
    else:
        sys.exit('No setup available for stackobjects = "TDOSEgalfitextractionsE60"')

    if verbose: print(' - Stack spectra')
    wave_out, flux_out, variance_out, Nspecstack = \
        stacking.stack_1D(wavelengths, fluxes, variances, z_systemic=z_systemic,
                          stacktype=stacktype, wavemin=1100, wavemax=3000,
                          deltawave=0.1, outfile=outfile, verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def stack_MUSEWideSpecs_plot(plotmultiple=False,plotSN=False):
    """
    Plotting MUSEWide LAE stack created with stack_MUSEWideSpecs_plot()

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu
    mwu.stack_MUSEWideSpecs_plot()

    """
    stackfile_test = '/Users/kschmidt/work/MUSE/spectral_stacking/MUSEWideStack_test181030.fits'
    spec_AGN       = '/Volumes/DATABCKUP1/TDOSEextractions/180824_TDOSEextraction_LAEs60fields/modelimg/tdose_spectra/tdose_spectrum_candels-cdfs-15_modelimg_0115003085-0115003085.fits'
    stackfile      = '/Users/kschmidt/work/MUSE/spectral_stacking/TDOSEgalfitextractionsE60_meanstack_181030.fits'
    if plotmultiple:
        spectra = [stackfile,stackfile_test,spec_AGN]
        labels  = ['674 obj stack','4 obj stack test','115003085 (AGN)']
    else:
        spectra = [stackfile]
        labels  = ['MUSE-Wide 60 Fields Ly$\\alpha$ emitter stack']

    Nspec   = len(spectra)
    plotname = stackfile.replace('.fits','_overview.pdf')

    xrangefull = [np.min(afits.open(stackfile)[1].data['wave']),np.max(afits.open(stackfile)[1].data['wave'])]
    yrangefull = [-100,900]

    zoomwindows=[[1190,1250,'Lya'],
                 [1300,1360,'CII'],
                 [1375,1410,'OIV'],
                 [1540,1555,'CIV'],
                 [1520,1540,'SiII'],
                 [1590,1650,'HeII'],
                 [1650,1675,'OIII'],
                 [1880,1900,'SiIII'],
                 [1900,1915,'CIII']]
                 # [2750,2900,'MgII']

    stacking.plot_1DspecOverview(spectra,labels,['wave']*Nspec,['flux']*Nspec,['fluxerror']*Nspec,plotname,
                                 plotSN=plotSN,xrangefull=xrangefull,yrangefull=yrangefull,zoomwindows=zoomwindows)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def TDOSE_sourcecat_from_Rafelski(outputnamebase,refimage,minRaper=0.5,minCutwidth=4.5,overwrite=False,verbose=True):
    """
    Generate a TDOSE source catatalog (and catalog with intitial guesses for Gaussian modeling) from
    the Rafelski+2015 UVUDF catalog.

    --- INPUT ---
    outputnamebase    The name base (prepended the path) of the source catalog to generate
    refimage          Provide fits image with WCS information to use for the pixel positions in the source cata
    overwrite         Overwrite output files if they exists?
    verbose           Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu

    imgpath        = '/Users/kschmidt/work/MUSE/TDOSE_fullextractions/Rafelski-UDF-full-v1p0/refimages/'
    refimage       = imgpath+'acs_775w_udf-04_cut_rot.fits'
    outputnamebase = refimage.replace('acs_775w','tdose_sourcecat_Rafelski_in_acs_775w').replace('.fits','')

    mwu.TDOSE_sourcecat_from_Rafelski(outputnamebase,refimage)

    """
    if verbose: print(' - Generating TDOSE source catalog from Rafelski+2015 catalog\n'
                      '   Restricting source to FoV of reference image '+refimage+'\n'
                                                                                  '   Output will be saved to '+outputnamebase+'.txt/fits')
    rafelskidat      = afits.open('/Users/kschmidt/work/catalogs/rafelski/uvudf_rafelski_2015.fits')[1].data
    refimgdata       = afits.open(refimage)[0].data
    refimghdr        = afits.open(refimage)[0].header
    arcsecPerPix_Raf = 0.03

    ids_all     = rafelskidat['ID']
    ras_all     = rafelskidat['RA']
    decs_all    = rafelskidat['DEC']
    fluxf_all   = rafelskidat['FLUX_ISO_F775W']
    ellip_all   = rafelskidat['ELLIPTICITY']
    theta_all   = rafelskidat['THETA']
    area_all    = rafelskidat['AREAF']
    a_area_all  = np.sqrt(area_all/np.pi/(1-ellip_all))
    b_area_all  = (1.0-ellip_all)*a_area_all
    x_raf_all   = rafelskidat['X']
    y_raf_all   = rafelskidat['Y']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Creating cutoutsizes.txt and aperturesizes.txt for the TDOSE extractions along the way ')
    if verbose: print('   (storing both files together with source catalog output) ')

    refimgpath  = '/'.join(refimage.split('/')[:-1])+'/'
    cut_outtxt  = refimgpath+'cutoutsizes_Rafelski_4xA_AREA.txt'
    aper_outtxt = refimgpath+'apertureradii_Rafelski_1xA_AREA.txt'

    fout = open(cut_outtxt,'w')
    fout.write('# Cutout sizes of 4*A_AREA on each side (corresponding to twice the aperture diameter '
               'used for the aperture extractions) estimated with '
               'MUSEWideUtilities.TDOSE_sourcecat_from_Rafelski() on '+kbs.DandTstr2()+'\n')
    fout.write('# \n')
    fout.write('# id xsize ysize \n')

    cutoutwidth = 4.0 * a_area_all * arcsecPerPix_Raf
    cutoutwidth[cutoutwidth < minCutwidth]  = minCutwidth
    cutoutwidth[cutoutwidth > 10.0] = 10.0000

    for ii, id in enumerate(ids_all):
        fout.write(str(ids_all[ii])+' '+
                   str(cutoutwidth[ii])+' '+
                   str(cutoutwidth[ii])+'  \n')
    fout.close()

    fout = open(aper_outtxt,'w')
    fout.write('# Aperturesizes of A_AREA estimated with '
               'MUSEWideUtilities.TDOSE_sourcecat_from_Rafelski() on '+kbs.DandTstr2()+'\n')
    fout.write('# \n')
    fout.write('# id aperturesize \n')

    Raper_arcsec = a_area_all * arcsecPerPix_Raf
    Raper_arcsec[Raper_arcsec < minRaper]  = minRaper

    for ii, id in enumerate(ids_all):
        fout.write(str(ids_all[ii])+' '+str(Raper_arcsec[ii])+'  \n')
    fout.close()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Finding objects within reference image field-of-view')
    if verbose: print('   (Estimating pixel positions using wcs info from header)')
    striphdr   = tu.strip_header(refimghdr,verbose=verbose)
    wcs_in     = wcs.WCS(striphdr)
    arcsecPerPix_refimg = wcs.utils.proj_plane_pixel_scales(wcs_in)*60*60

    skycoord   = SkyCoord(ras_all, decs_all, frame='fk5', unit='deg')
    pixcoord   = wcs.utils.skycoord_to_pixel(skycoord,wcs_in,origin=1)
    xpos       = pixcoord[0]
    ypos       = pixcoord[1]
    goodent    = np.where((xpos < refimghdr['NAXIS1']) & (xpos > 0) &
                          (ypos < refimghdr['NAXIS2']) & (ypos > 0)   )[0]

    if verbose: print('   (Make sure no 0s exist in a 6x6 pixel region around position, i.e., ignoring 0-edges of ref images.)')
    Ngoodinit  = len(goodent)
    intxpos    = np.round(xpos).astype(int)
    intypos    = np.round(ypos).astype(int)
    for ent in goodent:
        if (refimgdata[np.max([intypos[ent]-5,0]):np.min([intypos[ent]+5,refimghdr['NAXIS2']-1]),
            np.max([intxpos[ent]-5,0]):np.min([intxpos[ent]+5,refimghdr['NAXIS1']-1])] == 0).any():
            goodent[np.where(goodent == ent)[0]] = -99
    Nedge     = len(np.where(goodent == -99)[0])
    goodent   = goodent[np.where(goodent != -99)[0]]
    if verbose: print('   (Ended up removing '+str(Nedge)+'/'+str(Ngoodinit)+
                      ' objects that fall in the edge region but are within the image FoV)')

    ids             = ids_all[goodent]
    ras             = ras_all[goodent]
    decs            = decs_all[goodent]
    x_image         = xpos[goodent]
    y_image         = ypos[goodent]
    flux_iso_f775w  = fluxf_all[goodent]
    a_area          = a_area_all[goodent] * arcsecPerPix_Raf / np.mean(arcsecPerPix_refimg)
    b_area          = b_area_all[goodent] * arcsecPerPix_Raf / np.mean(arcsecPerPix_refimg)
    theta           = theta_all[goodent]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    outtxt          = outputnamebase+'.txt'
    if (overwrite == False) & os.path.isfile(outtxt):
        sys.exit('Output ('+outtxt+') already exists and clobber=False')
    else:
        if verbose: print(' - Will save source catalog to '+outtxt+' (overwriting any existing files)')
        fout = open(outtxt,'w')
        fout.write('# TDOSE Source catalog generated with '
                   'MUSEWideUtilities.TDOSE_sourcecat_from_Rafelski() on '+kbs.DandTstr2()+'\n')
        fout.write('# \n')
        fout.write('# parent_id id ra dec x_image y_image flux_iso_f775w a_area b_area theta \n')
        for ii, id in enumerate(ids):
            fout.write(str(ids[ii])+' '+str(ids[ii])+' '+str(ras[ii])+' '+str(decs[ii])+' '+
                       str(x_image[ii])+' '+str(y_image[ii])+' '+str(flux_iso_f775w[ii])+' '+
                       str(a_area[ii])+' '+str(b_area[ii])+' '+str(theta[ii])+'  \n')

        fout.close()

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        regionfile = outtxt.replace('.txt','.reg')
        if verbose: print(' - Storing DS9 region file to '+regionfile)
        idsstr     = [str(id) for id in ids]
        tu.create_simpleDS9region(regionfile,ras,decs,color='green',circlesize=a_area*np.mean(arcsecPerPix_refimg),
                                  textlist=idsstr,clobber=overwrite)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        outnamefits = outtxt.replace('.txt','.fits')
        if verbose: print(' - Saving fits version of source catalog to '+outnamefits)
        fitsfmt       = ['D']*10
        sourcecatfits = tu.ascii2fits(outtxt,asciinames=True,skip_header=2,fitsformat=fitsfmt,verbose=verbose)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def TDOSE_sourcecat_from_Skelton(outputnamebase,refimage,minRaper=0.5,minCutwidth=4.5,overwrite=False,verbose=True):
    """
    Generate a TDOSE source catatalog (and catalog with intitial guesses for Gaussian modeling) from
    the Skelton+2014 CDFS and COSMOS catalogs.

    --- INPUT ---
    outputnamebase    The name base (prepended the path) of the source catalog to generate
    refimage          Provide fits image with WCS information to use for the pixel positions in the source cata
    overwrite         Overwrite output files if they exists?
    verbose           Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu

    imgpath        = '/Users/kschmidt/work/images_MAST/MUSEWidePointings/'
    refimage       = imgpath+'acs_814w_XXX'
    outputnamebase = refimage.replace('acs_814w','tdose_sourcecat_skelton_in_acs_775w').replace('.fits','')

    mwu.TDOSE_sourcecat_from_Skelton(outputnamebase,refimage)

    """
    if verbose: print(' - Generating TDOSE source catalog from Skelton+2014 catalog\n'
                      '   Restricting source to FoV of reference image '+refimage+'\n'
                      '   Output will be saved to '+outputnamebase+'.txt/fits')

    if '-cosmos-' in refimage:
        skeltoncat   = '/Users/kschmidt/work/catalogs/skelton/cosmos_3dhst.v4.1.cats/Catalog/cosmos_3dhst.v4.1.cat.FITS'
    else:
        skeltoncat   = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
    skeltondat       = afits.open(skeltoncat)[1].data
    refimgdata       = afits.open(refimage)[0].data
    refimghdr        = afits.open(refimage)[0].header
    arcsecPerPix_Ske = 0.06

    ids_all     = skeltondat['id']
    ras_all     = skeltondat['ra']
    decs_all    = skeltondat['dec']
    fluxf_all   = skeltondat['f_f160w']
    theta_all   = skeltondat['theta_j2000']
    a_image_all = skeltondat['a_image']
    b_image_all = skeltondat['b_image']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Creating cutoutsizes.txt and aperturesizes.txt for the TDOSE extractions along the way ')
    if verbose: print('   (storing both files together with source catalog output) ')

    outpath     = '/'.join(outputnamebase.split('/')[:-1])+'/'

    if '-cosmos-' in refimage:
        cut_outtxt  = outpath+'../'+'cutoutsizes_Skelton_cosmos_4xA_IMAGE.txt'
        aper_outtxt = outpath+'../'+'apertureradii_Skelton_cosmos_1xA_IMAGE.txt'
    else:
        cut_outtxt  = outpath+'../'+'cutoutsizes_Skelton_goodss_4xA_IMAGE.txt'
        aper_outtxt = outpath+'../'+'apertureradii_Skelton_goodss_1xA_IMAGE.txt'

    fout = open(cut_outtxt,'w')
    fout.write('# Cutout sizes of 4*A_IMAGE on each side (corresponding to twice the aperture diameter '
               'used for the aperture extractions) estimated with '
               'MUSEWideUtilities.TDOSE_sourcecat_from_Skelton() on '+kbs.DandTstr2()+'\n')
    fout.write('# \n')
    fout.write('# id xsize ysize \n')

    cutoutwidth = 4.0 * a_image_all * arcsecPerPix_Ske
    cutoutwidth[cutoutwidth < minCutwidth]  = minCutwidth
    cutoutwidth[cutoutwidth > 10.0] = 10.0000

    for ii, id in enumerate(ids_all):
        fout.write(str(ids_all[ii])+' '+
                   str(cutoutwidth[ii])+' '+
                   str(cutoutwidth[ii])+'  \n')
    fout.close()

    fout = open(aper_outtxt,'w')
    fout.write('# Aperturesizes of A_AREA estimated with '
               'MUSEWideUtilities.TDOSE_sourcecat_from_Skelton() on '+kbs.DandTstr2()+'\n')
    fout.write('# \n')
    fout.write('# id aperturesize \n')

    Raper_arcsec = a_image_all * arcsecPerPix_Ske
    Raper_arcsec[Raper_arcsec < minRaper]  = minRaper

    for ii, id in enumerate(ids_all):
        fout.write(str(ids_all[ii])+' '+str(Raper_arcsec[ii])+'  \n')
    fout.close()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Finding objects within reference image field-of-view')
    if verbose: print('   (Estimating pixel positions using wcs info from header)')
    striphdr   = tu.strip_header(refimghdr,verbose=verbose,delkeys=['COMMENT','HISTORY','','A_ORDER','B_ORDER'])
    wcs_in     = wcs.WCS(striphdr)
    arcsecPerPix_refimg = wcs.utils.proj_plane_pixel_scales(wcs_in)*60*60

    skycoord   = SkyCoord(ras_all, decs_all, frame='fk5', unit='deg')
    pixcoord   = wcs.utils.skycoord_to_pixel(skycoord,wcs_in,origin=1)
    xpos       = pixcoord[0]
    ypos       = pixcoord[1]
    goodent    = np.where((xpos < refimghdr['NAXIS1']) & (xpos > 0) &
                          (ypos < refimghdr['NAXIS2']) & (ypos > 0) &
                          (a_image_all > 0.0) )[0]

    if verbose: print('   (Make sure no 0s exist in a 6x6 pixel region around position, i.e., ignoring 0-edges of ref images.)')
    Ngoodinit  = len(goodent)
    intxpos    = np.round(xpos).astype(int)
    intypos    = np.round(ypos).astype(int)
    for ent in goodent:
        if (refimgdata[np.max([intypos[ent]-5,0]):np.min([intypos[ent]+5,refimghdr['NAXIS2']-1]),
            np.max([intxpos[ent]-5,0]):np.min([intxpos[ent]+5,refimghdr['NAXIS1']-1])] == 0).any():
            goodent[np.where(goodent == ent)[0]] = -99
    Nedge     = len(np.where(goodent == -99)[0])
    goodent   = goodent[np.where(goodent != -99)[0]]
    if verbose: print('   (Ended up removing '+str(Nedge)+'/'+str(Ngoodinit)+
                      ' objects that fall in the edge region but are within the image FoV)')

    ids             = ids_all[goodent]
    ras             = ras_all[goodent]
    decs            = decs_all[goodent]
    x_image         = xpos[goodent]
    y_image         = ypos[goodent]
    flux_f160w      = fluxf_all[goodent]
    a_image         = a_image_all[goodent] * arcsecPerPix_Ske / np.mean(arcsecPerPix_refimg)
    b_image         = b_image_all[goodent] * arcsecPerPix_Ske / np.mean(arcsecPerPix_refimg)
    theta           = theta_all[goodent]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    outtxt          = outputnamebase+'.txt'
    if (overwrite == False) & os.path.isfile(outtxt):
        sys.exit('Output ('+outtxt+') already exists and clobber=False')
    else:
        if verbose: print(' - Will save source catalog to '+outtxt+' (overwriting any existing files)')
        fout = open(outtxt,'w')
        fout.write('# TDOSE Source catalog generated with '
                   'MUSEWideUtilities.TDOSE_sourcecat_from_Skelton() on '+kbs.DandTstr2()+'\n')
        fout.write('# \n')
        fout.write('# parent_id id ra dec x_image y_image flux_f160w a_image b_image theta \n')
        for ii, id in enumerate(ids):
            fout.write(str(ids[ii])+' '+str(ids[ii])+' '+str(ras[ii])+' '+str(decs[ii])+' '+
                       str(x_image[ii])+' '+str(y_image[ii])+' '+str(flux_f160w[ii])+' '+
                       str(a_image[ii])+' '+str(b_image[ii])+' '+str(theta[ii])+'  \n')

        fout.close()

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        regionfile = outtxt.replace('.txt','.reg')
        if verbose: print(' - Storing DS9 region file to '+regionfile)
        idsstr     = [str(id) for id in ids]
        tu.create_simpleDS9region(regionfile,ras,decs,color='green',circlesize=a_image*np.mean(arcsecPerPix_refimg),
                                  textlist=idsstr,clobber=overwrite)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        outnamefits = outtxt.replace('.txt','.fits')
        if verbose: print(' - Saving fits version of source catalog to '+outnamefits)
        fitsfmt       = ['D']*10
        sourcecatfits = tu.ascii2fits(outtxt,asciinames=True,skip_header=2,fitsformat=fitsfmt,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def TDOSE_sourcecat_from_Laigle(outputnamebase,refimage,minRaper=0.5,minCutwidth=4.5,overwrite=False,verbose=True):
    """
    Generate a TDOSE source catatalog (and catalog with intitial guesses for Gaussian modeling) from
    the Laigle+2016 COSMOS catalog.

    --- INPUT ---
    outputnamebase    The name base (prepended the path) of the source catalog to generate
    refimage          Provide fits image with WCS information to use for the pixel positions in the source cata
    overwrite         Overwrite output files if they exists?
    verbose           Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu

    imgpath        = '/Users/kschmidt/work/images_MAST/MUSEWidePointings/'
    refimage       = imgpath+'acs_814w_candels-cosmos'
    outputnamebase = refimage.replace('acs_814w','tdose_sourcecat_laigle_in_acs_775w').replace('.fits','')

    mwu.TDOSE_sourcecat_from_Laigle(outputnamebase,refimage)

    """
    if verbose: print(' - Generating TDOSE source catalog from Laigle+2016 catalog\n'
                      '   Restricting source to FoV of reference image '+refimage+'\n'
                      '   Output will be saved to '+outputnamebase+'.txt/fits')


    laiglecat        = '/Users/kschmidt/work/catalogs/laigle/COSMOS2015_Laigle_v1.1_candelsregion.fits'
    laigledat        = afits.open(laiglecat)[1].data
    refimgdata       = afits.open(refimage)[0].data
    refimghdr        = afits.open(refimage)[0].header
    arcsecPerPix_Lai = 0.15 # pixel scale from Table 9 of Laigle+16 (4.16666*10-5 * 3600 arcsec/deg)

    ids_all     = laigledat['NUMBER']
    ras_all     = laigledat['ALPHA_J2000']
    decs_all    = laigledat['DELTA_J2000']
    fluxf_all   = laigledat['FLUX_814W']
    r_img_all   = laigledat['FLUX_RADIUS']   # radius enclosing 0.5 of the total flux (FLUX_AUTO) in pixels

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Creating cutoutsizes.txt and aperturesizes.txt for the TDOSE extractions along the way ')
    if verbose: print('   (storing both files together with source catalog output) ')

    outpath     = '/'.join(outputnamebase.split('/')[:-1])+'/'


    cut_outtxt  = outpath+'../'+'cutoutsizes_Laigle_cosmos_4xFLUX_RADIUS.txt'
    aper_outtxt = outpath+'../'+'apertureradii_Laigle_cosmos_1xFLUX_RADIUS.txt'

    fout = open(cut_outtxt,'w')
    fout.write('# Cutout sizes of 4*FLUX_RADIUS on each side taken from Laigle catalog with '
               'MUSEWideUtilities.TDOSE_sourcecat_from_Laigle() on '+kbs.DandTstr2()+'\n')
    fout.write('# \n')
    fout.write('# id xsize ysize \n')

    cutoutwidth = 4.0 * r_img_all * arcsecPerPix_Lai
    cutoutwidth[cutoutwidth < minCutwidth]  = minCutwidth
    cutoutwidth[cutoutwidth > 10.0] = 10.0000

    for ii, id in enumerate(ids_all):
        fout.write(str(ids_all[ii])+' '+
                   str(cutoutwidth[ii])+' '+
                   str(cutoutwidth[ii])+'  \n')
    fout.close()

    fout = open(aper_outtxt,'w')
    fout.write('# Aperturesizes of FLUX_RADIUS taken from Laigle catalog with '
               'MUSEWideUtilities.TDOSE_sourcecat_from_Laigle() on '+kbs.DandTstr2()+'\n')
    fout.write('# \n')
    fout.write('# id aperturesize \n')

    Raper_arcsec = r_img_all * arcsecPerPix_Lai
    Raper_arcsec[Raper_arcsec < minRaper]  = minRaper

    for ii, id in enumerate(ids_all):
        fout.write(str(ids_all[ii])+' '+str(Raper_arcsec[ii])+'  \n')
    fout.close()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Finding objects within reference image field-of-view')
    if verbose: print('   (Estimating pixel positions using wcs info from header)')
    striphdr   = tu.strip_header(refimghdr,verbose=verbose,delkeys=['COMMENT','HISTORY','','A_ORDER','B_ORDER'])
    wcs_in     = wcs.WCS(striphdr)
    arcsecPerPix_refimg = wcs.utils.proj_plane_pixel_scales(wcs_in)*60*60

    skycoord   = SkyCoord(ras_all, decs_all, frame='fk5', unit='deg')
    pixcoord   = wcs.utils.skycoord_to_pixel(skycoord,wcs_in,origin=1)
    xpos       = pixcoord[0]
    ypos       = pixcoord[1]
    goodent    = np.where((xpos < refimghdr['NAXIS1']) & (xpos > 0) &
                          (ypos < refimghdr['NAXIS2']) & (ypos > 0) &
                          (r_img_all > 0.0) )[0]

    if verbose: print('   (Make sure no 0s exist in a 6x6 pixel region around position, i.e., ignoring 0-edges of ref images.)')
    Ngoodinit  = len(goodent)
    intxpos    = np.round(xpos).astype(int)
    intypos    = np.round(ypos).astype(int)
    for ent in goodent:
        if (refimgdata[np.max([intypos[ent]-5,0]):np.min([intypos[ent]+5,refimghdr['NAXIS2']-1]),
            np.max([intxpos[ent]-5,0]):np.min([intxpos[ent]+5,refimghdr['NAXIS1']-1])] == 0).any():
            goodent[np.where(goodent == ent)[0]] = -99
    Nedge     = len(np.where(goodent == -99)[0])
    goodent   = goodent[np.where(goodent != -99)[0]]
    if verbose: print('   (Ended up removing '+str(Nedge)+'/'+str(Ngoodinit)+
                      ' objects that fall in the edge region but are within the image FoV)')

    ids             = ids_all[goodent]
    ras             = ras_all[goodent]
    decs            = decs_all[goodent]
    x_image         = xpos[goodent]
    y_image         = ypos[goodent]
    flux_f814w      = fluxf_all[goodent]
    r_img_all       = r_img_all[goodent] * arcsecPerPix_Lai / np.mean(arcsecPerPix_refimg)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    outtxt          = outputnamebase+'.txt'
    if (overwrite == False) & os.path.isfile(outtxt):
        sys.exit('Output ('+outtxt+') already exists and clobber=False')
    else:
        if verbose: print(' - Will save source catalog to '+outtxt+' (overwriting any existing files)')
        fout = open(outtxt,'w')
        fout.write('# TDOSE Source catalog generated with (a_image = b_image = FLUX_RADIUS) '
                   'MUSEWideUtilities.TDOSE_sourcecat_from_Laigle() on '+kbs.DandTstr2()+'\n')
        fout.write('# \n')
        fout.write('# parent_id id ra dec x_image y_image flux_f814w a_image b_image theta \n')
        for ii, id in enumerate(ids):
            fout.write(str(ids[ii])+' '+str(ids[ii])+' '+str(ras[ii])+' '+str(decs[ii])+' '+
                       str(x_image[ii])+' '+str(y_image[ii])+' '+str(flux_f814w[ii])+' '+
                       str(r_img_all[ii])+'  '+str(r_img_all[ii])+'  '+str(r_img_all[ii]*0.0)+' \n')

        fout.close()

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        regionfile = outtxt.replace('.txt','.reg')
        if verbose: print(' - Storing DS9 region file to '+regionfile)
        idsstr     = [str(id) for id in ids]
        tu.create_simpleDS9region(regionfile,ras,decs,color='green',circlesize=r_img_all*np.mean(arcsecPerPix_refimg),
                                  textlist=idsstr,clobber=overwrite)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        outnamefits = outtxt.replace('.txt','.fits')
        if verbose: print(' - Saving fits version of source catalog to '+outnamefits)
        fitsfmt       = ['D']*10
        sourcecatfits = tu.ascii2fits(outtxt,asciinames=True,skip_header=2,fitsformat=fitsfmt,verbose=verbose)



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def TDOSE_sourcecat_from_Guo(outputnamebase,refimage,minRaper=0.5,minCutwidth=4.5,overwrite=False,verbose=True):
    """
    Generate a TDOSE source catatalog (and catalog with intitial guesses for Gaussian modeling) from
    the Guo+2013 CDFS catalog.

    --- INPUT ---
    outputnamebase    The name base (prepended the path) of the source catalog to generate
    refimage          Provide fits image with WCS information to use for the pixel positions in the source cata
    overwrite         Overwrite output files if they exists?
    verbose           Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu

    imgpath        = '/Users/kschmidt/work/images_MAST/MUSEWidePointings/'
    refimage       = imgpath+'acs_814w_candels-cdfs'
    outputnamebase = refimage.replace('acs_814w','tdose_sourcecat_guo_in_acs_775w').replace('.fits','')

    mwu.TDOSE_sourcecat_from_Guo(outputnamebase,refimage)

    """
    if verbose: print(' - Generating TDOSE source catalog from Guo+2013 catalog\n'
                      '   Restricting source to FoV of reference image '+refimage+'\n'
                      '   Output will be saved to '+outputnamebase+'.txt/fits')


    guocat        = '/Users/kschmidt/work/catalogs/guo/hlsp_candels_hst_wfc3_goodss-tot-multiband_f160w_v1_cat.fits'
    guodat        = afits.open(guocat)[1].data
    refimgdata    = afits.open(refimage)[0].data
    refimghdr     = afits.open(refimage)[0].header
    arcsecPerPix_Guo = 0.06 # pixel scale from page 4 of Guo+13

    ids_all     = guodat['ID']
    ras_all     = guodat['RA']
    decs_all    = guodat['DEC']
    fluxf_all   = guodat['ACS_F814W_FLUX']
    theta_all   = guodat['theta_image']
    a_image_all = guodat['a_image']
    b_image_all = guodat['b_image']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Creating cutoutsizes.txt and aperturesizes.txt for the TDOSE extractions along the way ')
    if verbose: print('   (storing both files together with source catalog output) ')

    outpath     = '/'.join(outputnamebase.split('/')[:-1])+'/'

    cut_outtxt  = outpath+'../'+'cutoutsizes_Guo_goodss_4xA_IMAGE.txt'
    aper_outtxt = outpath+'../'+'apertureradii_Guo_goodss_1xA_IMAGE.txt'

    fout = open(cut_outtxt,'w')
    fout.write('# Cutout sizes of 4*A_IMAGE on each side taken from Guo+13 catalog with '
               'MUSEWideUtilities.TDOSE_sourcecat_from_Guo() on '+kbs.DandTstr2()+'\n')
    fout.write('# \n')
    fout.write('# id xsize ysize \n')

    cutoutwidth = 4.0 * a_image_all * arcsecPerPix_Guo
    cutoutwidth[cutoutwidth < minCutwidth]  = minCutwidth
    cutoutwidth[cutoutwidth > 10.0] = 10.0000

    for ii, id in enumerate(ids_all):
        fout.write(str(ids_all[ii])+' '+
                   str(cutoutwidth[ii])+' '+
                   str(cutoutwidth[ii])+'  \n')
    fout.close()

    fout = open(aper_outtxt,'w')
    fout.write('# Aperturesizes of A_IMAGE taken from Guo+13 catalog with  '
               'MUSEWideUtilities.TDOSE_sourcecat_from_Guo() on '+kbs.DandTstr2()+'\n')
    fout.write('# \n')
    fout.write('# id aperturesize \n')

    Raper_arcsec = a_image_all * arcsecPerPix_Guo
    Raper_arcsec[Raper_arcsec < minRaper]  = minRaper

    for ii, id in enumerate(ids_all):
        fout.write(str(ids_all[ii])+' '+str(Raper_arcsec[ii])+'  \n')
    fout.close()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Finding objects within reference image field-of-view')
    if verbose: print('   (Estimating pixel positions using wcs info from header)')
    striphdr   = tu.strip_header(refimghdr,verbose=verbose,delkeys=['COMMENT','HISTORY','','A_ORDER','B_ORDER'])
    wcs_in     = wcs.WCS(striphdr)
    arcsecPerPix_refimg = wcs.utils.proj_plane_pixel_scales(wcs_in)*60*60

    skycoord   = SkyCoord(ras_all, decs_all, frame='fk5', unit='deg')
    pixcoord   = wcs.utils.skycoord_to_pixel(skycoord,wcs_in,origin=1)
    xpos       = pixcoord[0]
    ypos       = pixcoord[1]
    goodent    = np.where((xpos < refimghdr['NAXIS1']) & (xpos > 0) &
                          (ypos < refimghdr['NAXIS2']) & (ypos > 0) &
                          (a_image_all > 0.0) )[0]

    if verbose: print('   (Make sure no 0s exist in a 6x6 pixel region around position, i.e., ignoring 0-edges of ref images.)')
    Ngoodinit  = len(goodent)
    intxpos    = np.round(xpos).astype(int)
    intypos    = np.round(ypos).astype(int)
    for ent in goodent:
        if (refimgdata[np.max([intypos[ent]-5,0]):np.min([intypos[ent]+5,refimghdr['NAXIS2']-1]),
            np.max([intxpos[ent]-5,0]):np.min([intxpos[ent]+5,refimghdr['NAXIS1']-1])] == 0).any():
            goodent[np.where(goodent == ent)[0]] = -99
    Nedge     = len(np.where(goodent == -99)[0])
    goodent   = goodent[np.where(goodent != -99)[0]]
    if verbose: print('   (Ended up removing '+str(Nedge)+'/'+str(Ngoodinit)+
                      ' objects that fall in the edge region but are within the image FoV)')

    ids             = ids_all[goodent]
    ras             = ras_all[goodent]
    decs            = decs_all[goodent]
    x_image         = xpos[goodent]
    y_image         = ypos[goodent]
    flux_f814w      = fluxf_all[goodent]
    a_image         = a_image_all[goodent] * arcsecPerPix_Guo / np.mean(arcsecPerPix_refimg)
    b_image         = b_image_all[goodent] * arcsecPerPix_Guo / np.mean(arcsecPerPix_refimg)
    theta           = theta_all[goodent]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    outtxt          = outputnamebase+'.txt'
    if (overwrite == False) & os.path.isfile(outtxt):
        sys.exit('Output ('+outtxt+') already exists and clobber=False')
    else:
        if verbose: print(' - Will save source catalog to '+outtxt+' (overwriting any existing files)')
        fout = open(outtxt,'w')
        fout.write('# TDOSE Source catalog generated with '
                   'MUSEWideUtilities.TDOSE_sourcecat_from_Guo() on '+kbs.DandTstr2()+'\n')
        fout.write('# \n')
        fout.write('# parent_id id ra dec x_image y_image flux_f814w a_image b_image theta \n')
        for ii, id in enumerate(ids):
            fout.write(str(ids[ii])+' '+str(ids[ii])+' '+str(ras[ii])+' '+str(decs[ii])+' '+
                       str(x_image[ii])+' '+str(y_image[ii])+' '+str(flux_f814w[ii])+' '+
                       str(a_image[ii])+' '+str(b_image[ii])+' '+str(theta[ii])+'  \n')

        fout.close()

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        regionfile = outtxt.replace('.txt','.reg')
        if verbose: print(' - Storing DS9 region file to '+regionfile)
        idsstr     = [str(id) for id in ids]
        tu.create_simpleDS9region(regionfile,ras,decs,color='green',circlesize=a_image*arcsecPerPix_refimg,
                                  textlist=idsstr,clobber=overwrite)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        outnamefits = outtxt.replace('.txt','.fits')
        if verbose: print(' - Saving fits version of source catalog to '+outnamefits)
        fitsfmt       = ['D']*10
        sourcecatfits = tu.ascii2fits(outtxt,asciinames=True,skip_header=2,fitsformat=fitsfmt,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def TDOSE_sourcecat_from_Whitaker(outputnamebase,refimage,minRaper=0.5,minCutwidth=4.5,overwrite=False,verbose=True):
    """
    Generate a TDOSE source catatalog (and catalog with intitial guesses for Gaussian modeling) from
    the Whitaker+2019 Goods-S catalog.

    --- INPUT ---
    outputnamebase    The name base (prepended the path) of the source catalog to generate
    refimage          Provide fits image with WCS information to use for the pixel positions in the source cata
    overwrite         Overwrite output files if they exists?
    verbose           Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu

    imgpath        = '/Users/kschmidt/work/images_MAST/MUSEWidePointings/'
    refimage       = imgpath+'acs_814w_XXX'
    outputnamebase = refimage.replace('acs_814w','tdose_sourcecat_Whitaker_acs_814w').replace('.fits','')

    mwu.TDOSE_sourcecat_from_Whitaker(outputnamebase,refimage)

    ===> but see /Users/kschmidt/work/MUSE/TDOSE_fullextractions/readme.txt

    """
    if verbose: print(' - Generating TDOSE source catalog from Whitaker+2019 catalog\n'
                      '   Restricting source to FoV of reference image '+refimage+'\n'
                      '   Output will be saved to '+outputnamebase+'.txt/fits')

    whitcat          = '/Users/kschmidt/work/catalogs/whitaker/XXX.fits'
    whitdat          = afits.open(whitcat)[1].data
    refimgdata       = afits.open(refimage)[0].data
    refimghdr        = afits.open(refimage)[0].header
    arcsecPerPix_Whi = 0.06 # 60 mas

    ids_all     = whitdat['id']
    ras_all     = whitdat['ra']
    decs_all    = whitdat['dec']
    fluxf_all   = whitdat['f_f160w']
    theta_all   = whitdat['theta_j2000']
    a_image_all = whitdat['a_image']
    b_image_all = whitdat['b_image']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Creating cutoutsizes.txt and aperturesizes.txt for the TDOSE extractions along the way ')
    if verbose: print('   (storing both files together with source catalog output) ')

    outpath     = '/'.join(outputnamebase.split('/')[:-1])+'/'
    cut_outtxt  = outpath+'../'+'cutoutsizes_Whitaker_goodss_4xA_IMAGE.txt'
    aper_outtxt = outpath+'../'+'apertureradii_Whitaker_goodss_1xA_IMAGE.txt'

    fout = open(cut_outtxt,'w')
    fout.write('# Cutout sizes of 4*A_IMAGE on each side (corresponding to twice the aperture diameter '
               'used for the aperture extractions) estimated with '
               'MUSEWideUtilities.TDOSE_sourcecat_from_Whitaker() on '+kbs.DandTstr2()+'\n')
    fout.write('# \n')
    fout.write('# id xsize ysize \n')

    cutoutwidth = 4.0 * a_image_all * arcsecPerPix_Whi
    cutoutwidth[cutoutwidth < minCutwidth]  = minCutwidth
    cutoutwidth[cutoutwidth > 10.0] = 10.0000

    for ii, id in enumerate(ids_all):
        fout.write(str(ids_all[ii])+' '+
                   str(cutoutwidth[ii])+' '+
                   str(cutoutwidth[ii])+'  \n')
    fout.close()

    fout = open(aper_outtxt,'w')
    fout.write('# Aperturesizes of A_AREA estimated with '
               'MUSEWideUtilities.TDOSE_sourcecat_from_Whitaker() on '+kbs.DandTstr2()+'\n')
    fout.write('# \n')
    fout.write('# id aperturesize \n')

    Raper_arcsec = a_image_all * arcsecPerPix_Whi
    Raper_arcsec[Raper_arcsec < minRaper]  = minRaper

    for ii, id in enumerate(ids_all):
        fout.write(str(ids_all[ii])+' '+str(Raper_arcsec[ii])+'  \n')
    fout.close()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Finding objects within reference image field-of-view')
    if verbose: print('   (Estimating pixel positions using wcs info from header)')
    striphdr   = tu.strip_header(refimghdr,verbose=verbose,delkeys=['COMMENT','HISTORY','','A_ORDER','B_ORDER'])
    wcs_in     = wcs.WCS(striphdr)
    arcsecPerPix_refimg = wcs.utils.proj_plane_pixel_scales(wcs_in)*60*60

    skycoord   = SkyCoord(ras_all, decs_all, frame='fk5', unit='deg')
    pixcoord   = wcs.utils.skycoord_to_pixel(skycoord,wcs_in,origin=1)
    xpos       = pixcoord[0]
    ypos       = pixcoord[1]
    goodent    = np.where((xpos < refimghdr['NAXIS1']) & (xpos > 0) &
                          (ypos < refimghdr['NAXIS2']) & (ypos > 0) &
                          (a_image_all > 0.0) )[0]

    if verbose: print('   (Make sure no 0s exist in a 6x6 pixel region around position, i.e., ignoring 0-edges of ref images.)')
    Ngoodinit  = len(goodent)
    intxpos    = np.round(xpos).astype(int)
    intypos    = np.round(ypos).astype(int)
    for ent in goodent:
        if (refimgdata[np.max([intypos[ent]-5,0]):np.min([intypos[ent]+5,refimghdr['NAXIS2']-1]),
            np.max([intxpos[ent]-5,0]):np.min([intxpos[ent]+5,refimghdr['NAXIS1']-1])] == 0).any():
            goodent[np.where(goodent == ent)[0]] = -99
    Nedge     = len(np.where(goodent == -99)[0])
    goodent   = goodent[np.where(goodent != -99)[0]]
    if verbose: print('   (Ended up removing '+str(Nedge)+'/'+str(Ngoodinit)+
                      ' objects that fall in the edge region but are within the image FoV)')

    ids             = ids_all[goodent]
    ras             = ras_all[goodent]
    decs            = decs_all[goodent]
    x_image         = xpos[goodent]
    y_image         = ypos[goodent]
    flux_f160w      = fluxf_all[goodent]
    a_image         = a_image_all[goodent] * arcsecPerPix_Whi / np.mean(arcsecPerPix_refimg)
    b_image         = b_image_all[goodent] * arcsecPerPix_Whi / np.mean(arcsecPerPix_refimg)
    theta           = theta_all[goodent]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    outtxt          = outputnamebase+'.txt'
    if (overwrite == False) & os.path.isfile(outtxt):
        sys.exit('Output ('+outtxt+') already exists and clobber=False')
    else:
        if verbose: print(' - Will save source catalog to '+outtxt+' (overwriting any existing files)')
        fout = open(outtxt,'w')
        fout.write('# TDOSE Source catalog generated with '
                   'MUSEWideUtilities.TDOSE_sourcecat_from_Whitaker() on '+kbs.DandTstr2()+'\n')
        fout.write('# \n')
        fout.write('# parent_id id ra dec x_image y_image flux_f160w a_image b_image theta \n')
        for ii, id in enumerate(ids):
            fout.write(str(ids[ii])+' '+str(ids[ii])+' '+str(ras[ii])+' '+str(decs[ii])+' '+
                       str(x_image[ii])+' '+str(y_image[ii])+' '+str(flux_f160w[ii])+' '+
                       str(a_image[ii])+' '+str(b_image[ii])+' '+str(theta[ii])+'  \n')

        fout.close()

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        regionfile = outtxt.replace('.txt','.reg')
        if verbose: print(' - Storing DS9 region file to '+regionfile)
        idsstr     = [str(id) for id in ids]
        tu.create_simpleDS9region(regionfile,ras,decs,color='green',circlesize=a_image*np.mean(arcsecPerPix_refimg),
                                  textlist=idsstr,clobber=overwrite)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        outnamefits = outtxt.replace('.txt','.fits')
        if verbose: print(' - Saving fits version of source catalog to '+outnamefits)
        fitsfmt       = ['D']*10
        sourcecatfits = tu.ascii2fits(outtxt,asciinames=True,skip_header=2,fitsformat=fitsfmt,verbose=verbose)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =