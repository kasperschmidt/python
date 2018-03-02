# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                           Utilities for M3G studies
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import M3G_utilities as m3g
import astropy.io.fits as fits
import numpy as np
import tdose_utilities as tu
import tdose_modify_cube as tmc
from astropy import units as u
from astropy.coordinates import SkyCoord
import glob
import pdb
import os
import sys
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gen_tdosesourcecat(coordfile,outputfile,fluxscales=1.0,clobber=False,verbose=True):
    """
    Convert coord file to TDOSE-friendly source list

    --- INPUT ---
    coordfile     Input file containing information about sources.
    outputfile    TDOSE source catalog to generate
    fluxscales    Flux scales for sources; if float, same flux scale for all sources is used
    clobber       Overwrite existing output file?
    verbose       Toggle verbosity

    --- EXAMPLE OF USE ---
    import M3G_utilities as m3g

    coordfile  = '/Users/kschmidt/work/MUSE/M3G/data/PGC099522_sat_coord.txt'
    outputfile = '/Users/kschmidt/work/MUSE/M3G/tdose_sourcecats/tdose_sourcecats_PGC099522.txt'
    fluxscales = [1300.,300.,800.,300.,300.,130.,50.,5.]

    coordfile  = '/Users/kschmidt/work/MUSE/M3G/data/PGC065588_sat_coord.txt'
    outputfile = '/Users/kschmidt/work/MUSE/M3G/tdose_sourcecats/tdose_sourcecats_PGC065588.txt'
    fluxscales = [5000.,11000.,45000.,2800.,960.]

    outfits    = m3g.gen_tdosesourcecat(coordfile,outputfile,fluxscales=fluxscales,clobber=False)

    """
    dat = np.genfromtxt(coordfile,names=True,skip_header=0,dtype=None)

    if os.path.isfile(outputfile) & (clobber == False):
        sys.exit(' Output ('+outputfile+') exists and clobber=False')
    else:
        fout = open(outputfile,'w')
    fout.write('# id  x_image  y_image  ra  dec  fluxscale \n')

    ralist  = []
    declist = []
    for oo, object in enumerate(dat['num']):
        coord = SkyCoord(dat['RA'][oo]+' '+dat['DEC'][oo], unit=(u.hourangle, u.deg))

        if type(fluxscales) == float:
            fluxscale = fluxscales
        else:
            fluxscale = fluxscales[oo]

        objstr = str(dat['num'][oo])+'  '+\
                 str(dat['X'][oo])+'  '+\
                 str(dat['Y'][oo])+'  '+\
                 str(coord.ra.degree)+'  '+\
                 str(coord.dec.degree)+'  '+\
                 str(fluxscale)+'  \n'

        ralist.append(coord.ra.degree)
        declist.append(coord.dec.degree)

        fout.write(objstr)
    fout.close()
    if verbose: print(' - Wrote TDOSE-friendly source catalog to '+outputfile+' (and corresponding region file)')
    tu.create_simpleDS9region(outputfile.split('.')[0]+'.reg',ralist,declist,
                              color='red',circlesize=0.5,textlist=dat['num'].astype(str),clobber=clobber)
    if verbose: print('   and the corresponding DS9 region file ')

    outputfile_fits = tu.ascii2fits(outputfile,asciinames=True,skip_header=0,outpath=None,verbose=verbose)
    return outputfile_fits

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def convertSex2Fits(sexcat,verbose=True):
    """
    Convert a SExtractor catalog into fits file.

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import M3G_utilities as m3g
    m3g.convertSex2Fits('/Users/kschmidt/work/MUSE/M3G/sextractorruns/PCG065588_idedit.cat')

    """
    tu.SExtractorCat2fits([sexcat],stringcols=[],header=12,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def removesatellites(datacube,sourcemodelcube,outnameext='tdose_nosattelites',mainsources=[],dataext='DATA',
                     clobber=False,verbose=True):
    """
    Wrapper for removing the satellites from the original data cube

    --- INPUT ---
    datacube            Data cube to remove sources from
    sourcemodelcube     The source model cube to use for removing sources
    outnameext          Extension to append modified data cube
    mainsources         The source(s) index to keep in the data cube. Everything else (assumed to be contaminating
                        satellites) is removed. A link between source number (index+1) and ID is provided in
                        *tdose_modelimage_ds9_gauss.reg or by going through the source model cube
    dataext             Extension containing data to modify in datacube
    clobber             Overwrite existing files
    verbose             Toggle verbosity

    --- EXAMPLE OF USE ---
    import M3G_utilities as m3g

    datacube        = '/Users/kschmidt/work/MUSE/M3G/tdose_cutouts/PGC065588_poor_sat_20150713_id65588_cutout62p0x62p0arcsec.fits'
    sourcemodelcube = '/Users/kschmidt/work/MUSE/M3G/tdose_models/PGC065588_poor_sat_20150713_id65588_cutout62p0x62p0arcsec_tdose_source_modelcube_gauss.fits'
    sources         = [1]

    modified_cube   = m3g.removesatellites(datacube,sourcemodelcube,mainsources=sources,clobber=False)

    """
    tmc.remove_object(datacube,sourcemodelcube,objects=mainsources,remove=False,dataext=dataext,sourcemodelext=1,
                      savecube=outnameext,clobber=clobber,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =