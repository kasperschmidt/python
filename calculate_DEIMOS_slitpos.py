#----------------------------
#   NAME
#----------------------------
# calculate_DEIMOS_slitpos.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Set of functions to calculate the expected pixel postion of objects in a DEIMOS mask
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2014-03-28  started by K. B. Schmidt (UCSB)
#----------------------------
#   MODULES
#----------------------------
import calculate_DEIMOS_slitpos as cds
import pyfits
import numpy as np
import os
import sys
import pdb
#-------------------------------------------------------------------------------------------------------------
def estimatepos_all(maskfits,pixscale=0.1185,verbose=True):
    """
    Load the info fits file for a mask and calculate the expected position
    of each object in the slits.

    -- INPUT --
    maskfits       Information fits file for DEIMOS mask.
    pixscale       The spatial pixel scale to use for conversion
                   Default is 0.1185 taken from http://www2.keck.hawaii.edu/inst/deimos/specs.html

    -- EXAMPLE OF USE --
    import calculate_DEIMOS_slitpos as cds
    slitpos = cds.estimatepos_all(maskfits='/Users/kasperborelloschmidt/work/observing/131030_DEIMOS/GOODS-S_for_Kasper/goodss_mask3.fits')

    cds.printextractcmds(slitpos,aperturewidth=12)

    """
    if not os.path.isfile(maskfits):
        sys.exit('Could not find '+maskfits+' --> ABORTING')

    if verbose: print ' - Loading data from fits tables in ',maskfits
    hdulist = pyfits.open(maskfits)

    objectcat  = hdulist[1].data
    DesiSlits  = hdulist[4].data
    SlitObjMap = hdulist[5].data
    Nspec      = len(DesiSlits['dSlitId'])
    if verbose: print ' - Found ',Nspec,' objects with spectra'

    slitpos = {}
    for ss in xrange(Nspec):
        sID   = SlitObjMap['dSlitId'][ss]                               # Slit ID
        sname = str("%.3d" % sID)                                       # Slit name
        oID   = SlitObjMap['ObjectId'][sID == SlitObjMap['dSlitId']][0] # Object ID
        oname = objectcat['OBJECT'][oID == objectcat['ObjectID']][0]      # Object name

        slitlen = DesiSlits['slitLen'][ss]
        topdist = SlitObjMap['TopDist'][ss]
        botdist = SlitObjMap['BotDist'][ss]

        tolerance = 0.005
        if topdist+botdist-slitlen > tolerance:
            print '   WARNING The length of the slit '+sname+' and the summed offsets differ by >',tolerance

        toppix = cds.calcpixpos(topdist,pixscale,verbose=verbose)
        botpix = cds.calcpixpos(botdist,pixscale,verbose=verbose)

        if verbose: print '   [sname,oname,slitlen,top",bot",toppix,botpix] = ['+sname+', '+oname+', '+\
                          str("%.3f" % slitlen)+', '+\
                          str("%.3f" % topdist)+', '+\
                          str("%.3f" % botdist)+', '+\
                          str("%.3f" % toppix)+', '+\
                          str("%.3f" % botpix)+']'

        slitpos[sname] = [sname,oname,slitlen,topdist,botdist,toppix,botpix]

    return slitpos
#-------------------------------------------------------------------------------------------------------------
def calcpixpos(position,pixscale,verbose=True):
    """
    Convertin position in arcseconds to pixels given a pixel scale
    """
    pixpos = position/pixscale

    return pixpos
#-------------------------------------------------------------------------------------------------------------
def printextractcmds(slitpos,aperturewidth=12):
    """
    Print the extraction commands for each objects in the slitpos dictionary
    """
    aw = str(int(aperturewidth))

    for key in np.sort(slitpos.keys()):
        values = slitpos[key]
        offset = str(int(np.round(values[5])))
        print """
print, '---------------------- SLIT %s (OBJECT %s) ----------------------'
spec = extract1d('reduction/slit.gsm3.%sB.fits',%s,%s)
mwrfits,spec,'spectra/slit.gsm3.%sB_specIDL.fits'
spawn, 'python ~/PYTHONCODE/mattspec.py spectra/slit.gsm3.%sB_specIDL.fits'

spec = extract1d('reduction/slit.gsm3.%sR.fits',%s,%s)
mwrfits,spec,'spectra/slit.gsm3.%sR_specIDL.fits'
spawn, 'python ~/PYTHONCODE/mattspec.py spectra/slit.gsm3.%sR_specIDL.fits'

spawn, 'python ~/PYTHONCODE/merge_sides.py %s'
spawn, 'python ~/PYTHONCODE/resp_correct.py spectra/slit.gsm3.%s_merged.fits'
        """ % (key,values[1],key,offset,aw,key,key,key,offset,aw,key,key,key,key)
#-------------------------------------------------------------------------------------------------------------
def getobjinfo(objects,verbose=True,field='GOODS-S'):
    """
    Extract the info for a given object from the corresponding catalog.

    Note: returns a numby void where the column names are used to access the data

    -- EXAMPLE OF USE --
    import calculate_DEIMOS_slitpos as cds
    objdata = cds.getobjinfo(['J033242.23-27553'],verbose=True,field='GOODS-S')

    """

    if field == 'GOODS-S':
        catalog = '/Users/kasperborelloschmidt/work/observing/131030_DEIMOS/idropcandidatesGOODSS_130821.txt'
        columns = ['Name','RA','Dec','zmagauto','zmagautoerr',
                   'zflux','zfluxerr','zmagaper','zmagapererr','sn_z850']
        dat     = np.genfromtxt(catalog,dtype='40a,f,f,f,f,f,f,f,f,f',names=columns)

        namecut = [name[:-3] for name in dat['Name']] # removing decimals of las part of name
        namevec = namecut
    else:
        sys.exit(' The field '+field+' is invalid --> ABORTING')

    Nobj   = len(objects)
    objent = np.atleast_1d(np.int_(np.zeros(Nobj)-99.0))
    for ii in xrange(Nobj):
        ent = np.where(objects[ii] == np.asarray(namevec))[0]
        if len(ent) == 0:
            if verbose: print ' WARNING Found no match in catalog for object "',objects[ii],'"'
        else:
            objent[ii] = ent

    if len(objent[objent != -99]) == 0:
        sys.exit(' Did not find any matches to the '+str(Nobj)+' objects provided --> ABORTING')
    else:
        objdat = dat[objent[objent != -99]]

    if verbose:
        print ' - The following data was returned from the catalog \n',catalog
        print '   ',columns
        print '   ',objdat

    return objdat
#-------------------------------------------------------------------------------------------------------------
def plotandinspect(specid,specdir='./',verbose=True,ds9=True):
    """
    Plot 1D spectrum and diusplay the 2D spectra in ds9

    -- EXAMPLE OF USE --
    << %cpaste full code into torrone iPython session >>
    specidsGSM1 =
    objdata = plotandinspect(specids,'/data1/homedirs/kschmidt/DEIMOS/gsm1_n12/',verbose=True)

    """
    import glob
    import commands
    import ale_specplot
    import pylab
    import sys
    import pdb
    import pyfits
    pylab.ion()

    Nspec = len(specid)
    for ii in xrange(Nspec):
        if verbose: print ' - Inspecting spectrum '+str("%.3d" % specid[ii])

        spec2DB = glob.glob(specdir+'reduction/*'+str("%.3d" % specid[ii])+'B_bgsub.fits')
        if len(spec2DB) != 1:
            sys.exit(' Found '+str(len(spec2DB))+' blue 2D spectra for spectrum '+str("%.3d" % specid[ii]))
        else:
            spec2DB = spec2DB[0]

        spec2DR = glob.glob(specdir+'reduction/*'+str("%.3d" % specid[ii])+'R_bgsub.fits')
        if len(spec2DR) != 1:
            sys.exit(' Found '+str(len(spec2DR))+' red 2D spectra for spectrum '+str("%.3d" % specid[ii]))
        else:
            spec2DR = spec2DR[0]

        spec1D  = glob.glob(specdir+'spectra/*'+str("%.3d" % specid[ii])+'_merged_fluxcalib.fits')
        if len(spec1D) != 1:
            sys.exit(' Found '+str(len(spec1D))+' 1D spectra for spectrum '+str("%.3d" % specid[ii]))
        else:
            spec1D = spec1D[0]

        sizeblue = pyfits.open(spec2DB)[0].shape
        sizered  = pyfits.open(spec2DR)[0].shape
        if sizeblue != sizered:
            if verbose: print '   WARNING Size of blue and red 2Dspec for ID '+\
                              str("%.3d" % specid[ii])+' differ: B=',sizeblue,' R=',sizered
        ds9cmd = 'ds9 -scale zscale -tile mode row '+spec2DB+' '+spec2DR+' '
        if ds9: ds9out = commands.getoutput(ds9cmd)

        # plot 1D spectrum
        spectrum = spec1D
        Spec = ale_specplot.SpecPlot(spectrum)

        plotcmd = ''
        while plotcmd != 'close':
            plotcmd = raw_input(" - What to do with the plot? (ready to move on then type 'close')? ")
            # - - - - - - - - - - - - - - - - - - - - - - - -
            if plotcmd.startswith('gauss'):
                if len(plotcmd) > 5:
                    val = plotcmd.split('(')[-1].split(')')[0]
                    if verbose: print '   Executing:   Spec.gauss('+val+')'
                    Spec.gauss(width=float(val))
                else:
                    if verbose: print '   Executing:   Spec.gauss()'
                    Spec.gauss()
            # - - - - - - - - - - - - - - - - - - - - - - - -
            if plotcmd.startswith('boxcar'):
                if len(plotcmd) > 6:
                    val = plotcmd.split('(')[-1].split(')')[0]
                    if verbose: print '   Executing:   Spec.boxcar('+val+')'
                    Spec.boxcar(width=float(val))
                else:
                    if verbose: print '   Executing:   Spec.boxcar()'
                    Spec.boxcar()
            # - - - - - - - - - - - - - - - - - - - - - - - -
            if plotcmd.startswith('id'):
                wave = plotcmd.split('(')[-1].split(',')[0]
                line = plotcmd.split(')')[0].split(',')[-1]
                if verbose: print '   Executing:   Spec.id('+wave+','+line+')'
                Spec.id(float(wave),float(line))
            # - - - - - - - - - - - - - - - - - - - - - - - -
            if plotcmd.startswith('redraw'):
                if verbose: print '   Executing:   Spec.redraw()'
                Spec.redraw()
            # - - - - - - - - - - - - - - - - - - - - - - - -
            if plotcmd.startswith('sky'):
                if verbose: print '   Executing:   Spec.sky()'
                Spec.sky()
            # - - - - - - - - - - - - - - - - - - - - - - - -
            if plotcmd.startswith('ds9'):
                if verbose: print '   Executing:   ds9 command, i.e. showing 2D spec'
                ds9out = commands.getoutput(ds9cmd)
            # - - - - - - - - - - - - - - - - - - - - - - - -
            if plotcmd == 'close':
                Spec.close()


#-------------------------------------------------------------------------------------------------------------
#                                                  END
#-------------------------------------------------------------------------------------------------------------
