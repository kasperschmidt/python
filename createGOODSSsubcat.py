#-------------------------------------------------------------------------------------------------------------
import numpy as np
import kbsutilities as kbs
import pdb
#-------------------------------------------------------------------------------------------------------------
def createsubcat(originalcat,outputcat,objlist,verbose=True,clobber=True):
    """
    Extract a list of objects from a GOODSS catalog

    -- EXAMPLE OF USE --
    incat   = 'hlsp_candels_hst_wfc3_goodss-tot-multiband_f160w_v1_cat.txt.reformat.original'
    outcat  = 'hlsp_candels_hst_wfc3_goodss-tot-multiband_f160w_v1_cat.txt.reformat'
    objlist = 'deflectors2.txt'
    import createGOODSSsubcat as cc
    cc.createsubcat(incat,outcat,objlist,verbose=True)

    """
    if verbose: print ' - Loading object list'
    obj   = np.genfromtxt(objlist,dtype='I',names='objid')
    objid = obj['objid']
    Nobj  = len(objid)
    if verbose: print '   Found ',Nobj,' objects to extract from catalog'

    if verbose: print ' - Loading goods catalog'
    dat = np.genfromtxt(originalcat,dtype=None,names=True)
    Ngoodsobj = len(dat['id'])
    if verbose: print '   Found ',Ngoodsobj,' objects in the catalog to extract data from'

    if verbose: print ' - Extracting information for objects '
    goodent = []
    for ii in xrange(Nobj):
        ent = np.where(dat['id'] == objid[ii])
        if len(ent[0]) == 1:
            goodent.append(ent[0])
    goodent = np.asarray(goodent)
    subdat  = dat[:][goodent]

    if verbose: print ' - Saving sub-catalog to ',outputcat
    kbs.void2file(subdat,outputcat,clobber=clobber)
#-------------------------------------------------------------------------------------------------------------
