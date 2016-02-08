# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import os
import commands
import sys
import glob
import pyfits
import numpy as np
import manipulate_3DHSTdata as m3d
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def rename_3DHST2GLASS(outputdirectory,id,id_GLASSmock=999,
                       parentdatadir='/Volumes/DATABCKUP3/GOODSdata/3D-HST/GOODSS_WFC3_V4.1.5/',
                       clobber=False,testing=False,verbose=True):
    """
    Rename a set of 3D-HST files to the GLASS naming scheme.
    This will enable inspection with GiG and GiGz

    --- INPUT ---
    outputdirectory     Name of directory where the renamed files will be put
    id                  3D-HST ID of object to rename files for (including cluster name, e.g.,
    id_GLASSmock        Mock GLASS id to assign to object (GiG and GiGz expect id on the format %.5d)
                        Use IDs < 1000 to enable including the 3D-HST field number in the ID.
    parentdatadir       Top level directory of 3D-HST data
    clobber             Overwrite files in directory?
    verbose             Toggle verbosity

    --- EXAMPLE OF USE ---
    import manipulate_3DHSTdata as m3d
    filenames = m3d.rename_3DHST2GLASS('./data_test/',43386)

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Searching for files with ID = '+str(id)
    if testing: # searching for files takes a while, so only do it if "testing" not set
        if verbose: print '   NB! "testing=True" so not searching for files, returning hardcoded list'
        filenames   = m3d.get_filenames(id=id)
    else:
        filenames   = glob.glob(parentdatadir+'*/*/*/*'+str(id)+'*')

    Nfiles      = len(filenames)
    if verbose: print ' - Found '+str(Nfiles)+' files in the sub-directories of \n   '+parentdatadir
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Setting up output directory = '+outputdirectory
    if os.path.isdir(outputdirectory):
        if clobber == False:
            sys.exit(' Directory already exists and "clobber = False" so aborting ')
        else:
            if verbose: print '   Directory already exists and "clobber = True" so will replace files'
    else:
        os.makedirs(outputdirectory)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for fn in filenames:
        fnbase           = fn.split('_'+str(id))[0]
        fnext            = fn.split('_'+str(id))[1]
        field, field_no  = fn.split('/')[-4].split('-')
        IDnew            = str("%.3d" % id_GLASSmock)+str("%.2d" % int(field_no))

        GLASSnamebase = field.upper()+'3DHST-555-'+str("%.5d" % int(field_no))+'_'+IDnew+'-G141'

        if fnext.endswith('.png'):
            if fn.endswith(str(id)+'.1D.png'):
                cpcmd = ' cp '+fn+' '+outputdirectory+GLASSnamebase+'.1D.png'
            elif fn.endswith(str(id)+'.2D.png'):
                cpcmd = ' cp '+fn+' '+outputdirectory+GLASSnamebase+'.2D.png'
            else:
                cpcmd = ' cp '+fn+' '+outputdirectory+GLASSnamebase+fnext
        elif fnext.endswith('.1D.fits'):
            cpcmd = ' cp '+fn+' '+outputdirectory+GLASSnamebase+'.1D.fits'
        elif fnext.endswith('.2D.fits'):
            if '-big' in fnbase: GLASSnamebase = GLASSnamebase.replace('G141','G102')
            cpcmd = ' cp '+fn+' '+outputdirectory+GLASSnamebase+'.2D.fits'
        else:
            if verbose: print '   No copy-case set up for '+fn.split('/')[-1]+' so not copied'
            cpcmd = None

        if cpcmd != None:
            cpout = commands.getoutput(cpcmd)
            if cpout != '': print cpout

            if '.2D.fits' not in cpcmd:
                cpcmd = cpcmd.replace(IDnew+'-G141',IDnew+'-G102')
                cpout = commands.getoutput(cpcmd)
                if cpout != '': print cpout

    return filenames

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def caruanaLAEs(outputdir,verbose=True):
    """
    copy/rename files for the list of LAEs J. Caruana sent on 151221 so they can be inspected with GiG and GiGz

    --- INPUT ---
    outputdir  Directory to store renamed data in
    verbose    Toggle verbosity

    --- EXAMPLE OF USE ---
    import manipulate_3DHSTdata as m3d
    outdir = '/Users/kschmidt/work/MUSE/Caruana_LAE_sources/CaruanaLAEs_data/'
    m3d.caruanaLAEs(outdir,verbose=True)

    """
    ID_skelton, ID_GLASSmock, Lyawaves, redshifts = m3d.caruanaLAE_info()

    for oo, obj in enumerate(ID_skelton):
        if verbose: print '\n - Copying/renaming data of object '+str(obj)+'/"ID_GLASS"='+\
                          str(ID_GLASSmock[oo])+' (Lya at z = '+str("%.2f" % redshifts[oo])+')'
        obj_files = m3d.rename_3DHST2GLASS(outputdir,ID_skelton[oo],id_GLASSmock=ID_GLASSmock[oo],
                                           testing=False,clobber=True,verbose=verbose)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def caruanaLAE_info():
    """
    Retuning ingormation about the Caruane LAEs

    """
    ID_skelton    = np.array([25614,16492,19953,19906,16007,15419,17484,15601,19546,12210])
    ID_GLASSmock  = np.array([1    ,2    ,3    ,4    ,5    ,6    ,7    ,8    ,9    ,10   ])
    Lyawaves      = np.array([7111.10742188,6935.99023438,6684.72607422,6684.72607422,6545.99560547,
                              6323.49023438,5762.24121094,5192.36279297,5055.98535156,4980.97900391])
    redshifts     = Lyawaves/1216.-1.

    return  ID_skelton, ID_GLASSmock, Lyawaves, redshifts
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def generateCaruanaLAEregions(verbose=True):
    """

    --- EXAMPLE OF USE ---
    import manipulate_3DHSTdata as m3d
    m3d.generateCaruanaLAEregions()

    """
    outputpath    = './CaruanaLAEs_data/'
    ID_skelton, ID_GLASSmock, Lyawaves, redshifts = m3d.caruanaLAE_info()

    for ll, lw in enumerate(Lyawaves):
        outputfile    = str("%.3d" % ID_GLASSmock[ll])+'_emissionlines.reg'
        m3d.create_ELds9region(lw, outname=outputfile, outpath=outputpath)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def create_ELds9region(lyawave,outname='ds9region_emissionlines.reg',outpath='./',markersize=200,verbose=True):
    """
    Region files showing Lya/OII, OIII, CIII and CIV

    --- INPUT ---
    lyawave    Lya wavelength of object
    verbose    Toggle verbosity

    --- EXAMPLE OF USE ---
    import manipulate_3DHSTdata as m3d
    create_ELds9region(7111.10742188)

    """
    zobj      = lyawave/1216.-1.
    zobj_low  = lyawave/3726.-1.
    if verbose: print ' - Generating region file for object at z='+str("%.2f" % zobj)+\
                      ' (or z='+str("%.2f" % zobj_low)+'):\n   '+outpath+outname

    oiiiwave  = 5007.*(zobj_low+1.)
    civwave   = 1549.*(zobj+1.)
    ciiiwave  = 1909.*(zobj+1.)
    mgiiwaveH  = 2788.*(zobj+1.)
    mgiiwaveL  = 2788.*(zobj_low+1.)
    halphawave = 6563.*(zobj_low+1.)

    regtext = """
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
wcs;
circle(%s,0.0,%s) # color=green width=2 text={Lya @ z~%s (or [OII] @ z~%s)}
circle(%s,0.0,%s) # color=magenta width=2 text={[OIII] @ z~%s}
circle(%s,0.0,%s) # color=green width=2 text={CIV @ z~%s}
circle(%s,0.0,%s) # color=green width=2 text={CIII] @ z~%s}
circle(%s,0.0,%s) # color=green width=2 text={MgII @ z~%s}
circle(%s,0.0,%s) # color=magenta width=2 text={MgII @ z~%s}
circle(%s,0.0,%s) # color=magenta width=2 text={Halpha @ z~%s}
""" % (lyawave, markersize, str("%.2f" % zobj),        str("%.2f" % zobj_low),
       oiiiwave, markersize, str("%.2f" % zobj_low),
       civwave,  markersize, str("%.2f" % zobj),
       ciiiwave, markersize, str("%.2f" % zobj),
       mgiiwaveH, markersize, str("%.2f" % zobj),
       mgiiwaveL, markersize, str("%.2f" % zobj_low),
       halphawave, markersize, str("%.2f" % zobj_low)
       )

    fout = open(outpath+outname,'w')
    fout.write(regtext)
    fout.close()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_filenames(id=43386,verbose=True):
    """
    Returning list of files for a given id

    --- INPUT ---
    id          id of object to return files for
    verbose     toggle verbosity

    """
    parentpath = '/Volumes/DATABCKUP3/GOODSdata/3D-HST/GOODSS_WFC3_V4.1.5/goodss-01/'
    files43386 = [parentpath+'1D/ASCII/goodss-01-G141_43386.1D.ascii',
                  parentpath+'1D/FITS/goodss-01-G141_43386.1D.fits',
                  parentpath+'1D/PNG/goodss-01-G141_43386.1D.png',
                  parentpath+'2D/FITS/goodss-01-G141_43386.2D.fits',
                  parentpath+'2D/PNG/goodss-01-G141_43386.2D.png',
                  parentpath+'BIG/2D/goodss-01-G141-big_43386.2D.fits',
                  parentpath+'BIG/ZFIT/goodss-01-G141-big_43386.new_zfit.fits',
                  parentpath+'LINE/DAT/goodss-01-G141_43386.linefit.dat',
                  parentpath+'LINE/FITS/goodss-01-G141_43386.linefit.fits',
                  parentpath+'LINE/PNG/goodss-01-G141_43386.linefit.chain.png',
                  parentpath+'LINE/PNG/goodss-01-G141_43386.linefit.png',
                  parentpath+'ZFIT/DAT/goodss-01-G141_43386.new_zfit.dat',
                  parentpath+'ZFIT/FITS/goodss-01-G141_43386.new_zfit.fits',
                  parentpath+'ZFIT/PNG/goodss-01-G141_43386.new_zfit.2D.png',
                  parentpath+'ZFIT/PNG/goodss-01-G141_43386.new_zfit.png',
                  parentpath+'ZFIT/PZ/goodss-01-G141_43386.new_zfit.pz.fits',
                  parentpath+'ZFIT/TILT/goodss-01-G141_43386.new_zfit_tilt.dat',
                  parentpath+'ZFIT/TILT/goodss-01-G141_43386.new_zfit_tilt.png']

    objdirectory = {}
    objdirectory[str(id)] = files43386

    try:
        if verbose: print ' - Returning file names for ID = '+str(id)+' from hardcoded dictionary'
        fns = objdirectory[str(id)]
    except:
        if verbose: print ' - No file names found for ID = '+str(id)+' in hardcoded dictionary'
        fns = None


    return fns

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =