# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Functions and scripts to analyze and deal with NEOGAL ionization models
# http://www.iap.fr/neogal/models.html
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import numpy as np
import glob
import astropy.io.fits as afits
import sys
import pdb
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import matplotlib
import NEOGALmodels as nm
import itertools
import literaturecollection_emissionlinestrengths as lce
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def load_model(Zgas,filepath='/Users/kschmidt/work/catalogs/NEOGALlines/nebular_emission/',verbose=True):
    """
    Loat the content of a model output

    --- INPUT ---
    Zgas       The interstellar gas metallicity of the model to load. If set to 'combined'
               the combined model output [generated with nm.combine_modeloutputs()] will be loaded.


    --- EXAMPLE OF USE ---
    import NEOGALmodels as nm
    modeldata = nm.load_model(0.0001)

    modeldata = nm.load_model('combined')

    modeldata  = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/nebular_emission/')
    modeldata2 = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/AGN_NLR_nebular_feltre16/')

    """
    if Zgas == 'combined':
        Zgasstring = Zgas
    else:
        Zgasstring = str(Zgas).split('.')[-1]
    if 'feltre' in filepath:
        filename   = filepath+'nlr_nebular_Z'+Zgasstring+'.txt'
    else:
        filename   = filepath+'nebular_emission_Z'+Zgasstring+'.txt'
    if verbose: print(' - Attempting to load data from:\n   '+filename)
    modeldata  = np.genfromtxt(filename,names=True,dtype=None)
    if verbose: print('   ... successful')
    return modeldata
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def linenames():
    """
    Return dictionary that can convert NEGOAL column names to flux catalog LINENAMEs

    --- EXAMPLE OF USE ---
    import NEOGALmodels as nm
    linenamesdic = nm.linenames()

    """
    linenamesdic = {}

    linenamesdic['ovi1']   = ['OVI $\\lambda$1032'                , 1031.9261,         'right'      , 'Morton1991tab2']
    linenamesdic['ovi2']   = ['OVI $\\lambda$1038'                , 1037.6167,         'left'       , 'Morton1991tab2']
    linenamesdic['lyb']    = ['Ly$\\beta$ $\\lambda$1025'         , 1025.7219,         'right'      , 'Morton1991tab5']
    linenamesdic['lya']    = ['Ly$\\alpha$ $\\lambda$1216'        , 1215.6737,         'right'      , 'Morton1991tab5']
    linenamesdic[ 'NV1240']    = ['NV $\\lambda$1239'                 , 1238.821 ,         'right'      , 'Morton1991tab5']
    linenamesdic['nv2']    = ['NV $\\lambda$1243'                 , 1242.804 ,         'left'       , 'Morton1991tab5']
    linenamesdic['cii']    = ['CII $\\lambda$1336'                , 1335.6627,         'right'      , 'Morton1991tab5']
    linenamesdic['Siiv1']  = ['SiIV $\\lambda$1394'               , 1393.755 ,         'right'      , 'Morton1991tab5']
    linenamesdic['oiv1']   = ['OIV $\\lambda$1397'                , 1397.232 ,         'right'      , 'Morton1991tab5']
    linenamesdic['oiv2']   = ['OIV $\\lambda$1400'                , 1399.780 ,         'left'       , 'Morton1991tab5']
    linenamesdic['Siiv2']  = ['SiIV $\\lambda$1403'               , 1402.770 ,         'left'       , 'Morton1991tab5']
    linenamesdic['CIV1548']   = ['CIV $\\lambda$1548'                , 1548.195 ,         'right'      , 'Morton1991tab5']
    linenamesdic['CIV1551']   = ['CIV $\\lambda$1551'                , 1550.770 ,         'left'       , 'Morton1991tab5']
    linenamesdic['HeII1640']   = ['HeII $\\lambda$1640'               , 1640.420 ,         'right'      , 'vandenberk+2001']
    linenamesdic['OIII1661'] = ['OIII] $\\lambda$1661'              , 1660.809 ,         'right'      , 'Morton1991tab2']
    linenamesdic['OIII1666'] = ['OIII] $\\lambda$1666'              , 1666.150 ,         'left'       , 'Morton1991tab2']
    linenamesdic['ciii1']  = ['[CIII] $\\lambda$1907'             , 1907.    ,         'right'      , 'stark+2015']
    linenamesdic['CIII1908']  = ['CIII] $\\lambda$1909'              , 1909.    ,         'left'       , 'stark+2015']
    linenamesdic['ciib']   = ['CII] $\\lambda$2326'               , 2326.113 ,         'right'      , 'Morton1991tab5']
    linenamesdic['mgii1']  = ['MgII] $\\lambda$2796'              , 2795.528 ,         'right'      , 'Morton1991tab5']
    linenamesdic['mgii2']  = ['MgII] $\\lambda$2803'              , 2802.705 ,         'left'       , 'Morton1991tab5']
    linenamesdic['OII3727']   = ['[OII] $\\lambda$3726'              , 3726.    ,         'right'      , 'Pradhan2006']
    linenamesdic['oii2']   = ['[OII] $\\lambda$3729'              , 3729.    ,         'left'       , 'Pradhan2006']

    return linenamesdic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def load_fluxes_LSDCatForcedRun(line1='CIV1548',line2='CIII1908',fluxcol='F_3KRON',
                                datadir='/Users/kschmidt/work/MUSE/ciii_candidates/fluxAndEWmeasurements/'
                                        'ForceFluxC3inMUSE_fullrun161031/',
                                redshiftcat='/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits',
                                convert_Flux2Lbol=False,verbose=True):
    """
    Loading the flux output from an LSDcat run and potentially turning it into bolometric luminoisity.
    Returning data array to be plotted on NEGOAL diagrams

    --- EXAMPLE OF USE ---
    import NEOGALmodels as nm
    lumarray = nm.load_fluxes_LSDCatForcedRun(line1='CIV1548',line2='CIII1908',fluxcol='F_3KRON',convert_Flux2Lbol=True)


    """
    if verbose: print(' - Loading redshift catalog: \n   '+redshiftcat)
    z_data       = afits.open(redshiftcat)[1].data
    linenamesdic = nm.linenames()
    if verbose: print(' - Grabbing files with flux measurements in data directory: \n   '+datadir)
    fluxfiles = glob.glob(datadir+'*_linelist_fluxes.fits')
    Nfiles    = len(fluxfiles)
    if Nfiles == 0:
        sys.exit("Didn't find any *_linelist_fluxes.fits files in datadir="+datadir)

    outputarray = np.ones([Nfiles,11])*-99

    for ff, ffile in enumerate(fluxfiles):
        f_data   = afits.open(ffile)[1].data
        objid    = ffile.split('/')[-1][:8]
        objent   = np.where(z_data['UNIQUE_ID'] == objid)[0]
        if len(objent) != 1:
            if verbose: print(' - WARNING Found '+str(len(objent))+' matches to '+str(objid)+' in redshift catalog')
            continue

        line1name = linenamesdic[line1][0]
        line2name = linenamesdic[line2][0]
        line1ent = np.where(f_data['LINENAME'] == line1name)[0]
        line2ent = np.where(f_data['LINENAME'] == line2name)[0]
        if (len(line1ent) != 1) or (len(line2ent) != 1):
            if verbose: print(' - WARNING No match in flux table for '+line1+' and '+line2+' for '+str(objid))
            continue

        redshift     = z_data['REDSHIFT'][objent]
        redshifterr  = z_data['REDSHIFT_ERR'][objent]
        lineflux1    = f_data[fluxcol][line1ent]
        linefluxerr1 = f_data[fluxcol+'_ERR'][line1ent]
        lineflux2    = f_data[fluxcol][line2ent]
        linefluxerr2 = f_data[fluxcol+'_ERR'][line2ent]

        if convert_Flux2Lbol:
            Lbol1, Lbolerr1 = nm.convert_Fline2Lbol(lineflux1,linefluxerr1,redshift,verbose=False)
            Lbol2, Lbolerr2 = nm.convert_Fline2Lbol(lineflux2,linefluxerr2,redshift,verbose=False)
        else:
            Lbol1, Lbolerr1 = -99, -99
            Lbol2, Lbolerr2 = -99, -99

        outputarray[ff,:] = int(objid), redshift, redshifterr, \
                            lineflux1, linefluxerr1, lineflux2, linefluxerr2, \
                            Lbol1, Lbolerr1, Lbol2, Lbolerr2,

    return outputarray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_LvsL(modeldata,line1='CIV1548',line2='CIII1908',plotname='./TESTPLOT.pdf',
              Zgas=False,logU=False,xid=0.3,nh=100,COratio=0.38,Mcutoff=100,
              logx=False,logy=False,logp1=False,logp2=False,fixxrange=False,fixyrange=False,
              showobs=None,noobserr=False,verbose=True):
    """
    Plotting the model grids (in luminoisity) for two lines against each other

    --- INPUT ---
    modeldata The model data to use for plotting.
    line1     line to plot on x-axis (column name in nebular emission model output file without brackets)
    line2     line to plot on y-axis (column name in nebular emission model output file without brackets)
    plotname  Name of plot to generate

    Zgas      The interstellar gas metallicity in units of solar;
              Choose between: False, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008,
                              0.01, 0.014, 0.017, 0.02, 0.03, 0.04
    logU      Logarithmic value of the ionization parameter;
              Choose between: False, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0
    xid       Dust-to-metal mass ratio
              Choose between: False, 0.1, 0.3, 0.5
    nh        Hydrogen gas density (per cubic cm)
              Choose between: False, 10, 100, 1000, 10000
    COratio   (C/O)/(C/O)sol, carbon-to-oxygen abundance ratio in units of the solar value, (C/O)sol=0.44
              Choose between: False,  0.1, 0.14, 0.2, 0.27, 0.38, 0.52, 0.72, 1.0, 1.4
    Mcutoff   upper mass cutoff of the Chabrier IMF
              Choose between: False, 100, 300

              The 6 model parameters are used to define what model grid to plot. Hence, two (and only two)
              of them should _not_ be specified, i.e., set to false. The grid of these two paramters will
              then be plotted while fixing the remaining 4 parameters to the provided values.
              See the default values for an example of a Zgas vs. logU grid with ficed xid, nh, COratio and
              Mcutoff.

    logx      Plot log x axis
    logy      Plot log y axis
    logp1     Plot colors of varying model parameter 1 in log
    logp2     Plot colors of varying model parameter 2 in log
    fixxrange To fix the x plotting range provide [ymin,ymax] with this keyword
    fixyrange To fix the y plotting range provide [ymin,ymax] with this keyword
    showobs   To overlay/show observed measurements provide an array with these data with shape (Nobjects,11)
              where the 11 columns contain:
                  objid             ID of object
                  redshift          Redshift of object           (can be a dummy value if not known - no used here)
                  redshifterr       Uncertainty on the redshift  (can be a dummy value if not known - no used here)
                  lineflux1         Flux of line1                (can be a dummy value if not known - no used here)
                  linefluxerr1      Uncertainty on flux of line1 (can be a dummy value if not known - no used here)
                  lineflux2         Flux of line2                (can be a dummy value if not known - no used here)
                  linefluxerr2      Uncertainty on flux of line2 (can be a dummy value if not known - no used here)
                  Lbol1             Bolometric luminoisity of line1 for the given flux and redshift
                  Lbolerr1          Uncertainty on bolometric luminoisity of line1
                  Lbol2             Bolometric luminoisity of line2 for the given flux and redshift
                  Lbolerr2          Uncertainty on bolometric luminoisity of line1
              An array on this format can be obtained from the LSDCat output from a forced flux run with
              nm.load_fluxes_LSDCatForcedRun()
    noobserr  To not show the error bars on the observations set noobserr=True
    verbose   Toggle verbosity

    --- EXAMPLE OF USE ---
    import NEOGALmodels as nm
    modeldata = nm.load_model('combined',verbose=True)
    nm.plot_LvsL(modeldata,line1='CIV1548',line2='CIII1908',Zgas=False,logU=False,xid=0.3,nh=100,COratio=0.38,Mcutoff=100,logx=True,logy=True,logp1=True)

    nm.plot_LvsL(modeldata,line1='CIV1548',line2='CIII1908',Zgas=False,logU=False,xid=0.1,nh=10,COratio=0.1,Mcutoff=100,logx=True,logy=True,logp1=True)


    line1     = 'CIV1548'
    line2     = 'CIII1908'
    modeldata = nm.load_model('combined',verbose=True)
    obsdata   = nm.load_fluxes_LSDCatForcedRun(line1=line1,line2=line2,fluxcol='F_3KRON',convert_Flux2Lbol=True)
    nm.plot_LvsL(modeldata,line1=line1,line2=line2,Zgas=False,logU=False,xid=0.3,nh=100,COratio=0.38,Mcutoff=100,logx=True,logy=True,logp1=True,showobs=obsdata,fixxrange=[1e3,1e9],fixyrange=[1e0,1e9],plotname='./TESTPLOT_wdata.pdf',noobserr=True)

    """
    NFalse    = 0
    freeparam = []
    inforstr  = ""
    # - - - - - - - - - - - - - - - - - - - - - - - -
    legenddic = {}
    legenddic['Zgas']     = r'Z$_\textrm{gas}$'
    legenddic['logUs']    = r'log$_\textrm{10}$(U)'
    legenddic['xid']      = r'$\xi_\textrm{d}$'
    legenddic['nh']       = r'n$_\textrm{H}$ / [cm$^3$]'
    legenddic['COCOsol']  = r'C/O / [C/O]$_\textrm{sun}$'
    legenddic['mup']      = r'M$_\textrm{cut IMF}$ / [M$_\textrm{sun}]$'
    # - - - - - - - - - - - - - - - - - - - - - - - -
    if not Zgas:
        Zgasrange = [0.0,1.0]
        NFalse    = NFalse + 1.0
        #inforstr  = inforstr+' Zgas:vary, '
        freeparam.append('Zgas')
    else:
        Zgasrange = [Zgas-1e-6,Zgas+1e-6]
        inforstr  = inforstr+' '+legenddic['Zgas']+'='+str(Zgas)+', '
    # - - - - - - - - - - - - - - - - - - - - - - - -
    if not logU:
        logUrange = [-5.0,0.0]
        NFalse    = NFalse + 1.0
        #inforstr  = inforstr+' logU:vary, '
        freeparam.append('logUs')
    else:
        logUrange = [logU-0.1,logU+0.1]
        inforstr  = inforstr+' '+legenddic['logUs']+'='+str(logU)+', '
    # - - - - - - - - - - - - - - - - - - - - - - - -
    if not xid:
        xidrange  = [0.0,0.6]
        NFalse    = NFalse + 1.0
        #inforstr  = inforstr+' xid:vary, '
        freeparam.append('xid')
    else:
        xidrange  = [xid-0.01,xid+0.01]
        inforstr  = inforstr+' '+legenddic['xid']+'='+str(xid)+', '
    # - - - - - - - - - - - - - - - - - - - - - - - -
    if not nh:
        nhrange   = [0.0,1.0e6]
        NFalse    = NFalse + 1.0
        #inforstr  = inforstr+' nH:vary, '
        freeparam.append('nh')
    else:
        nhrange   = [nh-1.0,nh+1.0]
        inforstr  = inforstr+' '+legenddic['nh']+'='+str(nh)+', '
    # - - - - - - - - - - - - - - - - - - - - - - - -
    if not COratio:
        COratiorange  = [0.0,2.0]
        NFalse        = NFalse + 1.0
        #inforstr      = inforstr+' C/O:vary, '
        freeparam.append('COCOsol')
    else:
        COratiorange  = [COratio-0.001,COratio+0.001]
        inforstr      = inforstr+' '+legenddic['COCOsol']+'='+str(COratio)+', '
    # - - - - - - - - - - - - - - - - - - - - - - - -
    if not Mcutoff:
        Mcutoffrange  = [0.0,400.0]
        NFalse        = NFalse + 1.0
        #inforstr      = inforstr+' Mcutoff:vary, '
        freeparam.append('mup')
    else:
        Mcutoffrange  = [Mcutoff-1.0,Mcutoff+1.0]
        inforstr      = inforstr+' '+legenddic['mup']+'='+str(Mcutoff)+', '
    # - - - - - - - - - - - - - - - - - - - - - - - -

    if NFalse != 2:
        sys.exit(' Two and only two of the model parameters (Zgas,logU,xid,nh,COratio,Mcutoff) '
                 'should be set to Flase to define the model grid; however it appears '+str(NFalse)+
                 ' parameters where not set')

    goodent   = np.where( (modeldata['Zgas']    > Zgasrange[0])    & (modeldata['Zgas']    < Zgasrange[1]) &
                          (modeldata['logUs']   > logUrange[0])    & (modeldata['logUs']   < logUrange[1]) &
                          (modeldata['xid']     > xidrange[0])     & (modeldata['xid']     < xidrange[1]) &
                          (modeldata['nh']      > nhrange[0])      & (modeldata['nh']      < nhrange[1]) &
                          (modeldata['COCOsol'] > COratiorange[0]) & (modeldata['COCOsol'] < COratiorange[1]) &
                          (modeldata['mup']     > Mcutoffrange[0]) & (modeldata['mup']     < Mcutoffrange[1]) )

    Ngoodent  = len(goodent[0])

    if Ngoodent > 1:
        if verbose: print(' - Getting data for '+str(Ngoodent)+' data points satisfying model selection ')
        param1 = modeldata[freeparam[0]][goodent]
        if logp1:
            param1 = np.log10(param1)

        param2 = modeldata[freeparam[1]][goodent]
        if logp2:
            param2 = np.log10(param2)

        lum1   = modeldata[line1][goodent]
        lum2   = modeldata[line2][goodent]
    else:
        if verbose: print(' WARNING: Less than 2 model grid points to plot; no output generated')
        return

    # - - - - - - - - - - - PLOTTING - - - - - - - - - - -
    if verbose: print(' - Setting up and generating plot')
    plotname = plotname
    fig = plt.figure(figsize=(9, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.99, bottom=0.10, top=0.95)
    Fsize    = 10
    lthick   = 1
    marksize = 3
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    plt.title(inforstr[:-2],fontsize=Fsize)

    margin  = 0.1
    dx      = np.abs(np.max(lum1)-np.min(lum1))
    dy      = np.abs(np.max(lum2)-np.min(lum2))


    if fixxrange:
        xrange      = fixxrange
    else:
        if logx:
            xrange  = [np.min(lum1)-np.min(lum1)/2.,np.max(lum1)+np.max(lum1)/2.]
        else:
            xrange  = [np.min(lum1)-dx*margin,np.max(lum1)+dx*margin]

    if fixyrange:
        yrange      = fixyrange
    else:
        if logy:
            yrange  = [np.min(lum2)-np.min(lum2)/2.,np.max(lum2)+np.max(lum2)/2.]
        else:
            yrange  = [np.min(lum2)-dy*margin,np.max(lum2)+dy*margin]

    # ------------ PARAM1 ------------
    cmap    = plt.cm.get_cmap('winter')
    cmin    = np.min(param1)
    cmax    = np.max(param1)
    colnorm = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, 30) #cmax-cmin)
    mm      = plt.cm.ScalarMappable(cmap=cmap)
    mm.set_array(cmaparr)
    cb1     = plt.colorbar(mm)#shrink=0.25

    pstr1   = legenddic[freeparam[0]]
    if logp1:
        pstr1 = r'log$_\textrm{10}$('+pstr1+')'

    cb1.set_label(pstr1+' (outer circle) - Fixed: black line')

    for p1 in np.unique(param1):
        p1col   = cmap(colnorm(p1))
        p1ent   = np.where(param1 == p1)

        plt.plot(lum1[p1ent],lum2[p1ent],'-',lw=lthick, color='k',zorder=1)

        plt.errorbar(lum1[p1ent],lum2[p1ent],xerr=None,yerr=None,
                     marker='o',lw=0, markersize=marksize*3,
                     markerfacecolor=p1col,ecolor=p1col,markeredgecolor = 'k',zorder=10)

    # ------------ PARAM2 ------------
    cmap    = plt.cm.get_cmap('spring')
    cmin    = np.min(param2)
    cmax    = np.max(param2)
    colnorm = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, 30) #cmax-cmin)
    mm      = plt.cm.ScalarMappable(cmap=cmap)
    mm.set_array(cmaparr)
    cb2     = plt.colorbar(mm)#shrink=0.25

    pstr2   = legenddic[freeparam[1]]
    if logp2:
        pstr2 = 'log10('+pstr2+')'

    cb2.set_label(pstr2+' (inner circle) - Fixed: gray line')

    for p2 in np.unique(param2):
        p2col   = cmap(colnorm(p2))
        p2ent   = np.where(param2 == p2)

        plt.plot(lum1[p2ent],lum2[p2ent],'-',lw=lthick, color='gray',zorder=2)

        plt.errorbar(lum1[p2ent],lum2[p2ent],xerr=None,yerr=None,
                     marker='o',lw=0, markersize=marksize*1.5,
                     markerfacecolor=p2col,ecolor=p2col,markeredgecolor = 'k',zorder=20)

    if showobs != None:
        for ii, objid in enumerate(showobs[:,0]):
            if (showobs[:,7][ii] > xrange[0]) & (showobs[:,7][ii] < xrange[1]) & \
                    (showobs[:,9][ii] > yrange[0]) & (showobs[:,9][ii] < yrange[1]):

                if noobserr:
                    obsxerr = None
                    obsyerr = None
                else:
                    obsxerr = showobs[:,8][ii]
                    obsyerr = showobs[:,10][ii]
                plt.errorbar(showobs[:,7][ii],showobs[:,9][ii],xerr=obsxerr,yerr=obsyerr,
                             marker='*',lw=lthick, markersize=marksize*2,
                             markerfacecolor='k',ecolor='k',markeredgecolor = 'k',zorder=30)

    plt.xlabel(r'L$_\textrm{'+line1+r'}$ / [3.826$\times$1e33 erg/s]/[Msun/yr]')
    plt.ylabel(r'L$_\textrm{'+line2+r'}$ / [3.826$\times$1e33 erg/s]/[Msun/yr]')

    plt.xlim(xrange)
    plt.ylim(yrange)

    if logx:
        plt.xscale('log')
    if logy:
        plt.yscale('log')

    #--------- LEGEND ---------
    # plt.errorbar(-1,-1,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='white', markersize=marksize*2,
    #              markerfacecolor='white',markeredgecolor = 'k',label='Ground-based spec')
    #
    # leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize},ncol=1,numpoints=1)
    #                  #bbox_to_anchor=(1.25, 1.03))  # add the legend
    # leg.get_frame().set_alpha(0.7)
    #--------------------------

    if verbose: print('   Saving plot to'+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_LvsL_multiple(line1='CIV1548',line2='CIII1908',line1range=[1e3,1e8],line2range=[1e0,1e8],
                       outputdir='./',verbose=True):
    """
    Plotting the model grids for all possible combinations of free and fixed parameters
    for two specific lines

    --- EXAMPLE OF USE ---
    import NEOGALmodels as nm
    nm.plot_LvsL_multiple(outputdir='NEOGALplots_clean1611XX/')

    """
    modeldata = nm.load_model('combined',verbose=verbose)

    if verbose: print(' - Putting together permutations of chosen setups for plotting')
    infodic          = {}
    infodic['Zgas']  = [False,0.0001,0.006,0.040], True
    infodic['logUs'] = [False,-1.0,-2.5,-4.0]    , False
    infodic['xid']   = [False,0.1,0.3,0.5]       , False
    infodic['nh']    = [False,10,100,1000,10000] , False
    infodic['CO']    = [False,0.1,0.38,1.4]      , False
    infodic['Mcut']  = [False,100,300]           , False

    variables = [infodic['Zgas'][0],infodic['logUs'][0],infodic['xid'][0],
                 infodic['nh'][0],infodic['CO'][0],infodic['Mcut'][0]]

    permutations            = list(itertools.product(*variables))
    permutations_with2false = [sublist for sublist in permutations if sublist.count(False) == 2.]
    Nplots                  = len(permutations_with2false)

    if verbose: print(' - With the restriction Nfalse=2 the setup will results in '+str(Nplots)+\
                      ' plots (if model data allows)')
    if verbose: print(' - These will be saved to the output directory: '+outputdir)
    for pp, perm in enumerate(permutations_with2false):
        Zval     = perm[0]
        Uval     = perm[1]
        Xival    = perm[2]
        Nhval    = perm[3]
        COval    = perm[4]
        Mval     = perm[5]

        plotname = outputdir+'NEOGALmodelgrid_Zgas'+str(Zval).replace('.','p')+\
                   '_logU'+str(Uval).replace('.','p')+\
                   '_xid'+str(Xival).replace('.','p')+\
                   '_nH'+str(Nhval).replace('.','p')+\
                   '_CO'+str(COval).replace('.','p')+\
                   '_Mcut'+str(Mval).replace('.','p')+'.pdf'

        plotname = plotname.replace('False','Free')

        if verbose:
            plotno    = pp+1
            infostr = ' - Generating plot '+str("%.4d" % plotno)+'/'+str("%.4d" % Nplots)+': '+plotname.split('/')[-1]+'    '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        if not Zval:
            logp1 = True
        else:
            logp1 = False

        nm.plot_LvsL(modeldata,line1=line1,line2=line2,logx=True,logy=True,logp1=logp1,logp2=False,verbose=False,
                     Zgas=Zval,logU=Uval,xid=Xival,nh=Nhval,COratio=COval,Mcutoff=Mval,
                     fixxrange=line1range,fixyrange=line2range,plotname=plotname)

    print('\n   ... done')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_lineratios(modeldata,modeldata2='None',line1='CIV1551',line2='CIII1908',line3='CIV1551',line4='HeII1640',
                    plotname='./TESTPLOT.pdf',Zgas=False,logU=False,xid=0.3,nh=100,COratio=0.38,Mcutoff=100,
                    logx=True,logy=True,logp1=False,logp2=False,fixxrange=False,fixyrange=False,
                    showobs=None,noobserr=False,verbose=True):
    """
    Plotting the model grids for lines ratios line1/line2 and line3/line4 against each other

    --- INPUT ---
    modeldata The model data to use for plotting.
    modeldata2 TO also plot the narrow-loine region models from Feltre et al. provide these here
    line1     line to plot on x-axis (column name in nebular emission model output file without brackets)
    line2     line to plot on y-axis (column name in nebular emission model output file without brackets)
    plotname  Name of plot to generate

    Zgas      The interstellar gas metallicity in units of solar;
              Choose between: False, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008,
                              0.01, 0.014, 0.017, 0.02, 0.03, 0.04 (0.05, 0.06, 0.07)
    logU      Logarithmic value of the ionization parameter;
              Choose between: False, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0
    xid       Dust-to-metal mass ratio
              Choose between: False, 0.1, 0.3, 0.5
    nh        Hydrogen gas density (per cubic cm)
              Choose between: False, 10, 100, 1000, 10000
    COratio   (C/O)/(C/O)sol, carbon-to-oxygen abundance ratio in units of the solar value, (C/O)sol=0.44
              Choose between: False,  0.1, 0.14, 0.2, 0.27, 0.38, 0.52, 0.72, 1.0, 1.4
    Mcutoff   upper mass cutoff of the Chabrier IMF
              Choose between: False, 100, 300

              The 6 model parameters are used to define what model grid to plot. Hence, two (and only two)
              of them should _not_ be specified, i.e., set to false. The grid of these two paramters will
              then be plotted while fixing the remaining 4 parameters to the provided values.
              See the default values for an example of a Zgas vs. logU grid with ficed xid, nh, COratio and
              Mcutoff.

    logx      Plot log x axis
    logy      Plot log y axis
    logp1     Plot colors of varying model parameter 1 in log
    logp2     Plot colors of varying model parameter 2 in log
    fixxrange To fix the x plotting range provide [ymin,ymax] with this keyword
    fixyrange To fix the y plotting range provide [ymin,ymax] with this keyword
    showobs   To overlay/show observed measurements provide an array with these data with shape (Nobjects,11)
              where the 11 columns contain:
                  objid             ID of object
                  redshift          Redshift of object           (can be a dummy value if not known - no used here)
                  redshifterr       Uncertainty on the redshift  (can be a dummy value if not known - no used here)
                  lineflux1         Flux of line1                (can be a dummy value if not known - no used here)
                  linefluxerr1      Uncertainty on flux of line1 (can be a dummy value if not known - no used here)
                  lineflux2         Flux of line2                (can be a dummy value if not known - no used here)
                  linefluxerr2      Uncertainty on flux of line2 (can be a dummy value if not known - no used here)
                  Lbol1             Bolometric luminoisity of line1 for the given flux and redshift
                  Lbolerr1          Uncertainty on bolometric luminoisity of line1
                  Lbol2             Bolometric luminoisity of line2 for the given flux and redshift
                  Lbolerr2          Uncertainty on bolometric luminoisity of line1
              An array on this format can be obtained from the LSDCat output from a forced flux run with
              nm.load_fluxes_LSDCatForcedRun()
    noobserr  To not show the error bars on the observations set noobserr=True
    verbose   Toggle verbosity

    --- EXAMPLE OF USE ---
    import NEOGALmodels as nm

    modeldata  = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/nebular_emission/')
    modeldata2 = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/AGN_NLR_nebular_feltre16/')

    line1='CIV1551'
    line2='HeII1640'
    line3='CIII1908'
    line4='HeII1640'

    nm.plot_lineratios(modeldata,line1=line1,line2=line2,line3=line3,line4=line4,Zgas=False,logU=False,xid=0.3,nh=100,COratio=0.38,Mcutoff=100,logx=True,logy=True,logp1=False,logp2=False,fixxrange=[1e-3,5e2],fixyrange=[1e-3,5e2],plotname='./TESTPLOT_CIVHeIIvsCIIIHeII.pdf',modeldata2=modeldata2)

    line1='CIV1551'
    line2='CIII1908'
    line3='CIII1908'
    line4='HeII1640'

    nm.plot_lineratios(modeldata,line1=line1,line2=line2,line3=line3,line4=line4,Zgas=False,logU=False,xid=0.3,nh=100,COratio=0.38,Mcutoff=100,logx=True,logy=True,logp1=False,logp2=False,fixxrange=[1e-3,5e2],fixyrange=[1e-3,5e2],plotname='./TESTPLOT_CIVCIIIvsCIIIHeII.pdf',modeldata2=modeldata2)


    line1='CIV1551'
    line2='HeII1640'
    line3='CIV1551'
    line4='CIII1908'

    nm.plot_lineratios(modeldata,line1=line1,line2=line2,line3=line3,line4=line4,Zgas=False,logU=False,xid=0.3,nh=100,COratio=0.38,Mcutoff=100,logx=True,logy=True,logp1=False,logp2=False,fixxrange=[1e-3,5e2],fixyrange=[1e-3,5e2],plotname='./TESTPLOT_CIVHeIIvsCIVCIII.pdf',modeldata2=modeldata2)


    """
    NFalse    = 0
    freeparam = []
    inforstr  = ""
    # - - - - - - - - - - - - - - - - - - - - - - - -
    legenddic = {}
    legenddic['Zgas']     = r'Z$_\textrm{gas}$'
    legenddic['logUs']    = r'log$_\textrm{10}$(U)'
    legenddic['xid']      = r'$\xi_\textrm{d}$'
    legenddic['nh']       = r'n$_\textrm{H}$ / [cm$^3$]'
    legenddic['COCOsol']  = r'C/O / [C/O]$_\textrm{sun}$'
    legenddic['mup']      = r'M$_\textrm{cut IMF}$ / [M$_\textrm{sun}]$'
    # - - - - - - - - - - - - - - - - - - - - - - - -
    if not Zgas:
        Zgasrange = [0.0,1.0]
        NFalse    = NFalse + 1.0
        #inforstr  = inforstr+' Zgas:vary, '
        freeparam.append('Zgas')
    else:
        Zgasrange = [Zgas-1e-6,Zgas+1e-6]
        inforstr  = inforstr+' '+legenddic['Zgas']+'='+str(Zgas)+', '
    # - - - - - - - - - - - - - - - - - - - - - - - -
    if not logU:
        logUrange = [-5.0,0.0]
        NFalse    = NFalse + 1.0
        #inforstr  = inforstr+' logU:vary, '
        freeparam.append('logUs')
    else:
        logUrange = [logU-0.1,logU+0.1]
        inforstr  = inforstr+' '+legenddic['logUs']+'='+str(logU)+', '
    # - - - - - - - - - - - - - - - - - - - - - - - -
    if not xid:
        xidrange  = [0.0,0.6]
        NFalse    = NFalse + 1.0
        #inforstr  = inforstr+' xid:vary, '
        freeparam.append('xid')
    else:
        xidrange  = [xid-0.01,xid+0.01]
        inforstr  = inforstr+' '+legenddic['xid']+'='+str(xid)+', '
    # - - - - - - - - - - - - - - - - - - - - - - - -
    if not nh:
        nhrange   = [0.0,1.0e6]
        NFalse    = NFalse + 1.0
        #inforstr  = inforstr+' nH:vary, '
        freeparam.append('nh')
    else:
        nhrange   = [nh-1.0,nh+1.0]
        inforstr  = inforstr+' '+legenddic['nh']+'='+str(nh)+', '
    # - - - - - - - - - - - - - - - - - - - - - - - -
    if not COratio:
        COratiorange  = [0.0,2.0]
        NFalse        = NFalse + 1.0
        #inforstr      = inforstr+' C/O:vary, '
        freeparam.append('COCOsol')
    else:
        COratiorange  = [COratio-0.001,COratio+0.001]
        inforstr      = inforstr+' '+legenddic['COCOsol']+'='+str(COratio)+', '
    # - - - - - - - - - - - - - - - - - - - - - - - -
    if not Mcutoff:
        Mcutoffrange  = [0.0,400.0]
        NFalse        = NFalse + 1.0
        #inforstr      = inforstr+' Mcutoff:vary, '
        freeparam.append('mup')
    else:
        Mcutoffrange  = [Mcutoff-1.0,Mcutoff+1.0]
        inforstr      = inforstr+' '+legenddic['mup']+'='+str(Mcutoff)+', '
    # - - - - - - - - - - - - - - - - - - - - - - - -

    if NFalse != 2:
        sys.exit(' Two and only two of the model parameters (Zgas,logU,xid,nh,COratio,Mcutoff) '
                 'should be set to Flase to define the model grid; however it appears '+str(NFalse)+
                 ' parameters where not set')

    # - - - - - - - - - - - - - - - - - - - - - - - -
    goodent   = np.where( (modeldata['Zgas']    > Zgasrange[0])    & (modeldata['Zgas']    < Zgasrange[1]) &
                          (modeldata['logUs']   > logUrange[0])    & (modeldata['logUs']   < logUrange[1]) &
                          (modeldata['xid']     > xidrange[0])     & (modeldata['xid']     < xidrange[1]) &
                          (modeldata['nh']      > nhrange[0])      & (modeldata['nh']      < nhrange[1]) &
                          (modeldata['COCOsol'] > COratiorange[0]) & (modeldata['COCOsol'] < COratiorange[1]) &
                          (modeldata['mup']     > Mcutoffrange[0]) & (modeldata['mup']     < Mcutoffrange[1]) )

    Ngoodent  = len(goodent[0])

    if Ngoodent > 1:
        if verbose: print(' - Getting data for '+str(Ngoodent)+' data points satisfying (SFR)model selection ')
        param1_1 = modeldata[freeparam[0]][goodent]
        if logp1:
            param1_1 = np.log10(param1_1)

        param1_2 = modeldata[freeparam[1]][goodent]
        if logp2:
            param1_2 = np.log10(param1_2)

        ratio1_1   = modeldata[line1][goodent]/modeldata[line2][goodent]
        ratio1_2   = modeldata[line3][goodent]/modeldata[line4][goodent]
    else:
        if verbose: print(' WARNING: Less than 2 (SFR)model grid points to plot; no output generated')
        return

    # - - - - - - - - - - - - - - - - - - - - - - - -
    if modeldata2 != 'None':
        goodent2   = np.where( (modeldata2['Zgas']    > Zgasrange[0])    & (modeldata2['Zgas']    < Zgasrange[1]) &
                              (modeldata2['logUs']   > logUrange[0])    & (modeldata2['logUs']   < logUrange[1]) &
                              (modeldata2['xid']     > xidrange[0])     & (modeldata2['xid']     < xidrange[1]) &
                              (modeldata2['nh']      > nhrange[0])      & (modeldata2['nh']      < nhrange[1])  )

        Ngoodent2  = len(goodent2[0])

        if Ngoodent > 1:
            if verbose: print(' - Getting data for '+str(Ngoodent2)+' data points satisfying (AGN)model selection ')
            param2_1 = modeldata2[freeparam[0]][goodent2]
            if logp1:
                param2_1 = np.log10(param2_1)

            param2_2 = modeldata2[freeparam[1]][goodent2]
            if logp2:
                param2_2 = np.log10(param2_2)

            l2s = ['x','x','x','x'] # line names to use for Feltre+16 file
            for ll, linestr in enumerate([line1,line2,line3,line4]):
                if '1908' in linestr:
                    l2 = linestr.replace('1908','1907')
                else:
                    l2 = linestr

                l2s[ll] = l2

            ratio2_1   = modeldata2[l2s[0]][goodent2]/modeldata2[l2s[1]][goodent2]
            ratio2_2   = modeldata2[l2s[2]][goodent2]/modeldata2[l2s[3]][goodent2]
        else:
            if verbose: print(' WARNING: Less than 2 (AGN)model grid points to plot; no output generated')
            return

    # - - - - - - - - - - - PLOTTING - - - - - - - - - - -
    if verbose: print(' - Setting up and generating plot')
    plotname = plotname
    fig = plt.figure(figsize=(9, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.99, bottom=0.10, top=0.95)
    Fsize    = 10
    lthick   = 1
    marksize = 3
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    plt.title(inforstr[:-2],fontsize=Fsize)

    margin  = 0.1
    dx      = np.abs(np.max(ratio1_1)-np.min(ratio1_1))
    dy      = np.abs(np.max(ratio1_2)-np.min(ratio1_2))


    if fixxrange:
        xrange      = fixxrange
    else:
        if logx:
            xrange  = [np.min(ratio1_1)-np.min(ratio1_1)/2.,np.max(ratio1_1)+np.max(ratio1_1)/2.]
        else:
            xrange  = [np.min(ratio1_1)-dx*margin,np.max(ratio1_1)+dx*margin]

    if fixyrange:
        yrange      = fixyrange
    else:
        if logy:
            yrange  = [np.min(ratio1_2)-np.min(ratio1_2)/2.,np.max(ratio1_2)+np.max(ratio1_2)/2.]
        else:
            yrange  = [np.min(ratio1_2)-dy*margin,np.max(ratio1_2)+dy*margin]

    # ------------ PARAM1 ------------
    cmap    = plt.cm.get_cmap('winter')
    cmin    = np.min(param1_1)
    cmax    = np.max(param1_1)
    colnorm = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, 30) #cmax-cmin)
    mm      = plt.cm.ScalarMappable(cmap=cmap)
    mm.set_array(cmaparr)
    cb1     = plt.colorbar(mm)#shrink=0.25

    pstr1   = legenddic[freeparam[0]]
    if logp1:
        pstr1 = r'log$_\textrm{10}$('+pstr1+')'

    cb1.set_label(pstr1+' (outer circle) - Fixed: black line')

    for p1 in np.unique(param1_1):
        p1col   = cmap(colnorm(p1))
        p1ent   = np.where(param1_1 == p1)

        plt.plot(ratio1_1[p1ent],ratio1_2[p1ent],'-',lw=lthick, color='k',zorder=1)

        plt.errorbar(ratio1_1[p1ent],ratio1_2[p1ent],xerr=None,yerr=None,
                     marker='o',lw=0, markersize=marksize*3,
                     markerfacecolor=p1col,ecolor=p1col,markeredgecolor = 'k',zorder=10)

        if modeldata2 is not 'None':
            p1ent   = np.where(param2_1 == p1)

            plt.plot(ratio2_1[p1ent],ratio2_2[p1ent],'-',lw=lthick, color='k',zorder=1)

            plt.errorbar(ratio2_1[p1ent],ratio2_2[p1ent],xerr=None,yerr=None,
                         marker='D',lw=0, markersize=marksize*3,
                         markerfacecolor=p1col,ecolor=p1col,markeredgecolor = 'k',zorder=10)


    # ------------ PARAM2 ------------
    cmap    = plt.cm.get_cmap('spring')
    cmin    = np.min(param1_2)
    cmax    = np.max(param1_2)
    colnorm = matplotlib.colors.Normalize(vmin=cmin,vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, 30) #cmax-cmin)
    mm      = plt.cm.ScalarMappable(cmap=cmap)
    mm.set_array(cmaparr)
    cb2     = plt.colorbar(mm)#shrink=0.25

    pstr2   = legenddic[freeparam[1]]
    if logp2:
        pstr2 = 'log10('+pstr2+')'

    cb2.set_label(pstr2+' (inner circle) - Fixed: gray line')

    for p2 in np.unique(param1_2):
        p2col   = cmap(colnorm(p2))
        p2ent   = np.where(param1_2 == p2)

        plt.plot(ratio1_1[p2ent],ratio1_2[p2ent],'-',lw=lthick, color='gray',zorder=2)

        plt.errorbar(ratio1_1[p2ent],ratio1_2[p2ent],xerr=None,yerr=None,
                     marker='o',lw=0, markersize=marksize*1.5,
                     markerfacecolor=p2col,ecolor=p2col,markeredgecolor = 'k',zorder=20)

        if modeldata2 is not 'None':
            p2ent   = np.where(param2_2 == p2)

            plt.plot(ratio2_1[p2ent],ratio2_2[p2ent],'-',lw=lthick, color='gray',zorder=2)

            plt.errorbar(ratio2_1[p2ent],ratio2_2[p2ent],xerr=None,yerr=None,
                         marker='D',lw=0, markersize=marksize*1.5,
                         markerfacecolor=p2col,ecolor=p2col,markeredgecolor = 'k',zorder=20)


    if showobs != None:
        for ii, objid in enumerate(showobs[:,0]):
            if (showobs[:,7][ii] > xrange[0]) & (showobs[:,7][ii] < xrange[1]) & \
                    (showobs[:,9][ii] > yrange[0]) & (showobs[:,9][ii] < yrange[1]):

                if noobserr:
                    obsxerr = None
                    obsyerr = None
                else:
                    obsxerr = showobs[:,8][ii]
                    obsyerr = showobs[:,10][ii]
                plt.errorbar(showobs[:,7][ii],showobs[:,9][ii],xerr=obsxerr,yerr=obsyerr,
                             marker='*',lw=lthick, markersize=marksize*2,
                             markerfacecolor='k',ecolor='k',markeredgecolor = 'k',zorder=30)

    plt.xlabel(line1+'/'+line2)
    plt.ylabel(line3+'/'+line4)

    plt.xlim(xrange)
    plt.ylim(yrange)

    if logx:
        plt.xscale('log')
    if logy:
        plt.yscale('log')

    #--------- LEGEND ---------
    # plt.errorbar(-1,-1,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='white', markersize=marksize*2,
    #              markerfacecolor='white',markeredgecolor = 'k',label='Ground-based spec')
    #
    # leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize},ncol=1,numpoints=1)
    #                  #bbox_to_anchor=(1.25, 1.03))  # add the legend
    # leg.get_frame().set_alpha(0.7)
    #--------------------------

    if verbose: print('   Saving plot to'+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def combine_modeloutputs(outputname='xxRENAMExx_Zcombined.txt',
                         data='sfr',
                         verbose=True):
    """
    Combine the model outputs to have a single 'master-model' with all variables included

    --- EXAMPLE OF USE ---
    import NEOGALmodels as nm
    nm.combine_modeloutputs(outputname='nebular_emission_Zcombined.txt',data='sfr')
    nm.combine_modeloutputs(outputname='nlr_nebular_Zcombined.txt',data='agn')

    """
    if data == 'sfr':
        filepath     = '/Users/kschmidt/work/catalogs/NEOGALlines/nebular_emission/'
        modelfilestr = filepath+'nebular_emission_Z0*.txt'
        splitstr     = 'emission_Z'
    elif data == 'agn':
        filepath     = '/Users/kschmidt/work/catalogs/NEOGALlines/AGN_NLR_nebular_feltre16/'
        modelfilestr = filepath+'nlr_nebular_Z0*.txt'
        splitstr     = 'nebular_Z'
    else:
        sys.exit('Inavlid value of data="'+data+'"')

    output     = filepath+outputname
    if verbose: print(' - Setting up output for:\n   '+output)
    modelfiles = glob.glob(modelfilestr)
    header     = open(modelfiles[0]).readline().rstrip()
    if data == 'sfr':
        header     = header.replace('##','#  Zgas  ')
    elif data == 'agn':
        header     = header.replace('#','#  Zgas  ')
        header     = header+'\n'

    fout = open(output, 'w')
    fout.write(header)
    if verbose: print(' - Writing the following files to ouput:')
    for mf in modelfiles:
        if verbose: print('   '+mf)
        Zgasstring   = mf.split('/')[-1].split(splitstr)[-1].split('.txt')[0]

        with open(mf, 'r') as f:
            linesall = f.readlines()

        for linestring in linesall:
            if linestring.startswith('#'):
                pass
            elif linestring == ' \n':
                fout.write(linestring)
            else:
                fout.write('0.'+Zgasstring+'    '+linestring)

    fout.close()
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def convert_Fline2Lbol(lineflux,linefluxerr,redshift,verbose=True):
    """
    Converting an observed integrated line flux [erg/s/cm2] to bolometric luminoisity [erg/s]

    --- EXAMPLE OF USE ---
    import NEOGALmodels as nm
    LbolLsun, LbolLsunerr = nm.convert_Fline2Lbol(2000,100,4.1)

    """
    cosmo       = FlatLambdaCDM(H0=70, Om0=0.3)
    if verbose: print(' - Estimating bolometric luminoisity for flat standard cosmology (H0=70, Om0=0.3, OL0=0.7)')
    DL          = cosmo.luminosity_distance(redshift).value
    # DLplus    = cosmo.luminosity_distance(redshift+redshifterr).value
    # DLminus   = cosmo.luminosity_distance(redshift-redshifterr).value
    Mpc2cm      = 3.086                          # 10**24 cm/Mpc
    Asphere     = 4*np.pi*(DL*Mpc2cm)**2         # 10**48 cm2
    Lbol        = lineflux*Asphere               # 10**28 erg/s ; assuming line fluxes are in 10**-20 erg/s/cm2
    Lbolerr     = linefluxerr*Asphere            # 10**28 erg/s ; assuming line fluxes are in 10**-20 erg/s/cm2
    LbolLsun    = Lbol/3.826*10**-5              # in units of Lbol_sun = 3.826*10**33 erg/s
    LbolLsunerr = Lbolerr/3.826*10**-5           # in units of Lbol_sun = 3.826*10**33 erg/s
    if verbose: print(' - Retunring luminoisity in units of Lbol_sun = 3.826e33 erg/s')
    if verbose: print(' - Result is: '+str(LbolLsun)+' +/- '+str(LbolLsunerr)+' [3.826e33 erg/s]')
    return LbolLsun, LbolLsunerr
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def convert_Lbol2Fline(Lbol,redshift,verbose=True):
    """
    Converting bolometric luminoisity [erg/s] to integrated line flux [erg/s/cm2]
    """
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def estimate_object_PDFs(fluxratiodictionarylist,generatePDFplots=False,AGNcol='blue',SFcol='red',verbose=True):
    """
    Function to estimate "PDFs" from the SF and AGN NEOGAL models given a set of flux ratio measurements.

    --- INPUT ---
    fluxratiodictionarylist   List of flux ratio dictionaris to get "PDFs" for.
                              Length of list corresponds to objects with measurements to get "PDFs" for.
    AGNcol                    Color of AGN model related data in plots
    SFcol                     Color of SF model related data in plots
    verbose                   Toggle verbosity


    --- EXAMPLE OF USE ---
    import NEOGALmodels as nm
    FRdic = [{'id':111111, 'OIII1663/HeII1640':[1e-3,1e-2],'CIII1908/CIV1550':[1.0,10.0]}, {'id':222222, 'OIII1663/HeII1640':[1e-1,1.0],'CIII1908/CIV1550':[0.1,10.0]}, {'id':333333, 'OIII1663/HeII1640':[1e2,1e3],'CIII1908/CIV1550':[1e-3,1e-2]}, {'id':444444, 'OIII1663/HeII1640':[1e-2,1e-1],'CIII1908/CIV1550':[5e-1,1e-0]}]
    parametercollection_SF, parametercollection_AGN, stat_SF, stat_AGN = nm.estimate_object_PDFs(FRdic)

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Loading NEOGAL models ')
    SF_models  = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/nebular_emission/')

    AGN_models = nm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/NEOGALlines/AGN_NLR_nebular_feltre16/')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Define all possible line ratios from the lines:\n    '
                      'NV1240, CIV1550, CIII1908, HeII1640, OIII1663, and SiIII1888')
    fluxratiodic = {}                   # [[SF range], [AGN range]]
    fluxratiodic['NV1240/CIV1550']     = [[-999,-999],[-999,-999]]
    fluxratiodic['NV1240/CIII1908']    = [[-999,-999],[-999,-999]]
    fluxratiodic['NV1240/HeII1640']    = [[-999,-999],[-999,-999]]
    fluxratiodic['NV1240/OIII1663']    = [[-999,-999],[-999,-999]]
    fluxratiodic['NV1240/SiIII1888']   = [[-999,-999],[-999,-999]]

    fluxratiodic['CIV1550/NV1240']     = [[-999,-999],[-999,-999]]
    fluxratiodic['CIV1550/CIII1908']   = [[-999,-999],[-999,-999]]
    fluxratiodic['CIV1550/HeII1640']   = [[-999,-999],[-999,-999]]
    fluxratiodic['CIV1550/OIII1663']   = [[-999,-999],[-999,-999]]
    fluxratiodic['CIV1550/SiIII1888']  = [[-999,-999],[-999,-999]]

    fluxratiodic['CIII1908/NV1240']     = [[-999,-999],[-999,-999]]
    fluxratiodic['CIII1908/CIV1550']    = [[-999,-999],[-999,-999]]
    fluxratiodic['CIII1908/HeII1640']   = [[-999,-999],[-999,-999]]
    fluxratiodic['CIII1908/OIII1663']   = [[-999,-999],[-999,-999]]
    fluxratiodic['CIII1908/SiIII1888']  = [[-999,-999],[-999,-999]]

    fluxratiodic['HeII1640/NV1240']     = [[-999,-999],[-999,-999]]
    fluxratiodic['HeII1640/CIV1550']    = [[-999,-999],[-999,-999]]
    fluxratiodic['HeII1640/CIII1908']   = [[-999,-999],[-999,-999]]
    fluxratiodic['HeII1640/OIII1663']   = [[-999,-999],[-999,-999]]
    fluxratiodic['HeII1640/SiIII1888']  = [[-999,-999],[-999,-999]]

    fluxratiodic['OIII1663/NV1240']     = [[-999,-999],[-999,-999]]
    fluxratiodic['OIII1663/CIV1550']    = [[-999,-999],[-999,-999]]
    fluxratiodic['OIII1663/CIII1908']   = [[-999,-999],[-999,-999]]
    fluxratiodic['OIII1663/HeII1640']   = [[-999,-999],[-999,-999]]
    fluxratiodic['OIII1663/SiIII1888']  = [[-999,-999],[-999,-999]]

    fluxratiodic['SiIII1888/NV1240']    = [[-999,-999],[-999,-999]]
    fluxratiodic['SiIII1888/CIV1550']   = [[-999,-999],[-999,-999]]
    fluxratiodic['SiIII1888/CIII1908']  = [[-999,-999],[-999,-999]]
    fluxratiodic['SiIII1888/HeII1640']  = [[-999,-999],[-999,-999]]
    fluxratiodic['SiIII1888/OIII1663']  = [[-999,-999],[-999,-999]]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Set up mode line flux vectors')
    fluxdic = {}         # [[SF flxu],                                   [AGN flux]]
    fluxdic['NV1240']    = [SF_models['NV1240'],                           AGN_models['NV1240']]
    fluxdic['CIV1550']   = [SF_models['CIV1548']+SF_models['CIV1551'],     AGN_models['CIV1548']+AGN_models['CIV1551']]
    fluxdic['CIII1908']  = [SF_models['CIII1908'],                         AGN_models['CIII1907']+AGN_models['CIII1910']]
    fluxdic['HeII1640']  = [SF_models['HeII1640'],                         AGN_models['HeII1640']]
    fluxdic['OIII1663']  = [SF_models['OIII1661']+SF_models['OIII1666'],  AGN_models['OIII1661']+AGN_models['OIII1666']]
    fluxdic['SiIII1888'] = [SF_models['SiIII1888'],  AGN_models['SiIII1888']]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Get ranges of model flux ratios')
    for FR in fluxratiodic.keys():
        numerator   = FR.split('/')[0]
        denominator = FR.split('/')[1]
        fluxratiodic[FR][0] = [np.min(fluxdic[numerator][0]/fluxdic[denominator][0]),
                               np.max(fluxdic[numerator][0]/fluxdic[denominator][0])]
        fluxratiodic[FR][1] = [np.min(fluxdic[numerator][1]/fluxdic[denominator][1]),
                               np.max(fluxdic[numerator][1]/fluxdic[denominator][1])]


    Nobj = len(fluxratiodictionarylist)
    if verbose: print(' - Get model selection given flux ratio ranges according to '+
                      str(Nobj)+" object's data provided ")

    parametercollection_SF  = [{'id':0, 'Zgas':[],'logUs':[],'xid':[],'nh':[],'COCOsol':[],'mup':[]}]*Nobj
    parametercollection_AGN = [{'id':0, 'Zgas':[],'logUs':[],'xid':[],'nh':[],'alpha':[]}]*Nobj

    for oo, FRdic_input in enumerate(fluxratiodictionarylist):
        fluxratiodic_obj = fluxratiodic    # resetting flux ratio dictionary for object
        for FR in FRdic_input.keys():
            if FR in fluxratiodic.keys():
                fluxratiodic_obj[FR] = [FRdic_input[FR],FRdic_input[FR]]
            elif FR == 'id':
                pass
            else:
                print(' WARNING nm.estimate_object_PDFs(): The flux ratio entry '+FR+' is not availble in the \n'
                      '                                    dictionary from the NEOGAL models. Define that flux \n'
                      '                                    ratio or correct input data.')

        goodent_SF   = np.arange(len(SF_models))
        goodent_AGN  = np.arange(len(AGN_models))
        for FR in fluxratiodic.keys():
            numerator   = FR.split('/')[0]
            denominator = FR.split('/')[1]

            goodent_FR_SF  = np.where( (fluxdic[numerator][0]/fluxdic[denominator][0] >= fluxratiodic_obj[FR][0][0]) &
                                       (fluxdic[numerator][0]/fluxdic[denominator][0] <= fluxratiodic_obj[FR][0][1]))[0]
            goodent_SF = np.intersect1d(goodent_SF,goodent_FR_SF)

            goodent_FR_AGN = np.where( (fluxdic[numerator][1]/fluxdic[denominator][1] >= fluxratiodic_obj[FR][1][0]) &
                                       (fluxdic[numerator][1]/fluxdic[denominator][1] <= fluxratiodic_obj[FR][1][1]))[0]
            goodent_AGN    = np.intersect1d(goodent_AGN,goodent_FR_AGN)


        parametercollection_SF[oo]  = {'id'     : FRdic_input['id'],
                                       'Zgas'   : SF_models['Zgas'][goodent_SF],
                                       'logUs'  : SF_models['logUs'][goodent_SF],
                                       'xid'    : SF_models['xid'][goodent_SF],
                                       'nh'     : SF_models['nh'][goodent_SF],
                                       'COCOsol': SF_models['COCOsol'][goodent_SF],
                                       'mup'    : SF_models['mup'][goodent_SF]}

        parametercollection_AGN[oo] = {'id'     : FRdic_input['id'],
                                       'Zgas'   : AGN_models['Zgas'][goodent_AGN],
                                       'logUs'  : AGN_models['logUs'][goodent_AGN],
                                       'xid'    : AGN_models['xid'][goodent_AGN],
                                       'nh'     : AGN_models['nh'][goodent_AGN],
                                       'alpha'  : AGN_models['alpha'][goodent_AGN]}

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Getting distribution ranges (percentiles) for parameter collections ')
    stat_SF   = []
    stat_AGN  = []

    for oo in np.arange(Nobj):
        stat_SF.append({'id':parametercollection_SF[oo]['id'], 'Zgas':[],'logUs':[],'xid':[],'nh':[],'COCOsol':[],'mup':[]})
        for key in stat_SF[oo].keys():
            if key == 'id': continue

            if len(parametercollection_SF[oo][key]) > 0:

                meanval_SF   = np.mean(parametercollection_SF[oo][key])
                std_SF       = np.std(parametercollection_SF[oo][key])
                medianval_SF = np.median(parametercollection_SF[oo][key])
                perc2p5_SF   = np.sort(parametercollection_SF[oo][key])[int(len(parametercollection_SF[oo][key])*0.025)]
                perc16_SF    = np.sort(parametercollection_SF[oo][key])[int(len(parametercollection_SF[oo][key])*0.16)]
                perc25_SF    = np.sort(parametercollection_SF[oo][key])[int(len(parametercollection_SF[oo][key])*0.25)]
                perc50_SF    = np.sort(parametercollection_SF[oo][key])[int(len(parametercollection_SF[oo][key])*0.50)]
                perc75_SF    = np.sort(parametercollection_SF[oo][key])[int(len(parametercollection_SF[oo][key])*0.75)]
                perc84_SF    = np.sort(parametercollection_SF[oo][key])[int(len(parametercollection_SF[oo][key])*0.84)]
                perc97p5_SF  = np.sort(parametercollection_SF[oo][key])[int(len(parametercollection_SF[oo][key])*0.975)]

                stat_SF[oo][key] = [meanval_SF,std_SF,medianval_SF,perc2p5_SF,perc16_SF,perc25_SF,
                                    perc50_SF,perc75_SF,perc84_SF,perc97p5_SF]
            else:
                stat_SF[oo][key] = [np.nan]*10

        stat_AGN.append({'id':parametercollection_AGN[oo]['id'], 'Zgas':[],'logUs':[],'xid':[],'nh':[],'alpha':[]})
        for key in stat_AGN[oo].keys():
            if key == 'id': continue

            if len(parametercollection_AGN[oo][key]) > 0:
                meanval_AGN   = np.mean(parametercollection_AGN[oo][key])
                std_AGN       = np.std(parametercollection_AGN[oo][key])
                medianval_AGN = np.median(parametercollection_AGN[oo][key])
                perc2p5_AGN   = np.sort(parametercollection_AGN[oo][key])[int(len(parametercollection_AGN[oo][key])*0.025)]
                perc16_AGN    = np.sort(parametercollection_AGN[oo][key])[int(len(parametercollection_AGN[oo][key])*0.16)]
                perc25_AGN    = np.sort(parametercollection_AGN[oo][key])[int(len(parametercollection_AGN[oo][key])*0.25)]
                perc50_AGN    = np.sort(parametercollection_AGN[oo][key])[int(len(parametercollection_AGN[oo][key])*0.50)]
                perc75_AGN    = np.sort(parametercollection_AGN[oo][key])[int(len(parametercollection_AGN[oo][key])*0.75)]
                perc84_AGN    = np.sort(parametercollection_AGN[oo][key])[int(len(parametercollection_AGN[oo][key])*0.84)]
                perc97p5_AGN  = np.sort(parametercollection_AGN[oo][key])[int(len(parametercollection_AGN[oo][key])*0.975)]

                stat_AGN[oo][key] = [meanval_AGN,std_AGN,medianval_AGN,perc2p5_AGN,perc16_AGN,perc25_AGN,
                                     perc50_AGN,perc75_AGN,perc84_AGN,perc97p5_AGN]
            else:
                stat_AGN[oo][key] = [np.nan]*10

    stat_idlist                =  [stat_AGN[oo]['id'] for oo in np.arange(len(stat_AGN))]
    parametercollection_idlist =  [parametercollection_AGN[oo]['id'] for oo in np.arange(len(parametercollection_AGN))]

    if stat_idlist != parametercollection_idlist:
        sys.exit(' NEOGALmodels.estimate_object_PDFs(): Wait a minute... the ID lists are not identical between \n'
                 '                                      the parameter collection ('+str(parametercollection_idlist)+') and \n'
                 '                                      the stats ('+str(stat_idlist)+')')

    plotname = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/NEOGALpdffigures/NEOGALobjectStats.pdf'
    nm.plot_stat(plotname,stat_SF,stat_AGN,SFcol=SFcol,AGNcol=AGNcol,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if generatePDFplots:
        if verbose: print(' - Plotting the extracted model parameter collections')
        plotname = '/Users/kschmidt/work/MUSE/uvEmissionlineSearch/NEOGALpdffigures/NEOGALobjectPDF.pdf'
        nm.plot_modelparametercollections(plotname, parametercollection_SF, parametercollection_AGN,
                                          stat_SF, stat_AGN, AGNcol=AGNcol,SFcol=SFcol,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return parametercollection_SF, parametercollection_AGN, stat_SF, stat_AGN

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_stat(plotname, stat_SF, stat_AGN, SFcol='red', AGNcol='blue', verbose=True):
    """

    """
    if verbose: print(' - Plot stats of parameter selection')

    xkey   = 'Zgas'
    xlabel = 'Z'
    ykey   = 'logUs'
    ylabel = 'log(U)'

    xrange = [0.00001,0.2]
    yrange = [-5.2,-0.4]

    plotname_obj = plotname.replace('.pdf','_meanANDstd_ZgasVSlogUs.pdf')
    nm.plot_stat_diagram(plotname_obj,stat_SF,stat_AGN,xkey,ykey,xlabel,ylabel,SFcol,AGNcol,xrange,yrange,
                         xlog=True,ylog=False,datsetup='meanstd',verbose=verbose)


    plotname_obj = plotname.replace('.pdf','_median68percent_ZgasVSlogUs.pdf')
    nm.plot_stat_diagram(plotname_obj,stat_SF,stat_AGN,xkey,ykey,xlabel,ylabel,SFcol,AGNcol,xrange,yrange,
                         xlog=True,ylog=False,datsetup='med68',verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_stat_diagram(plotname,stat_SF,stat_AGN,xkey,ykey,xlabel,ylabel,SFcol,AGNcol,xrange,yrange,
                      xlog=False,ylog=False,datsetup='med68',verbose=True):

    if datsetup == 'med68':
        cen_ent     = 2
        errlow_ent  = 4
        errhigh_ent = 8
    elif datsetup == 'medqart':
        cen_ent     = 2
        errlow_ent  = 5
        errhigh_ent = 7
    elif datsetup == 'med95':
        cen_ent     = 2
        errlow_ent  = 3
        errhigh_ent = 9
    elif datsetup == 'meanstd':
        cen_ent     = 0
        errlow_ent  = 1
        errhigh_ent = 1

    xvalues_SF  = []
    xerrlow_SF  = []
    xerrhigh_SF = []

    yvalues_SF  = []
    yerrlow_SF  = []
    yerrhigh_SF = []

    xvalues_AGN  = []
    xerrlow_AGN  = []
    xerrhigh_AGN = []

    yvalues_AGN  = []
    yerrlow_AGN  = []
    yerrhigh_AGN = []

    for oo in np.arange(len(stat_SF)):
        xvalues_SF.append(stat_SF[oo][xkey][cen_ent])
        xerrlow_SF.append(stat_SF[oo][xkey][errlow_ent])
        xerrhigh_SF.append(stat_SF[oo][xkey][errhigh_ent])

        yvalues_SF.append(stat_SF[oo][ykey][cen_ent])
        yerrlow_SF.append(stat_SF[oo][ykey][errlow_ent])
        yerrhigh_SF.append(stat_SF[oo][ykey][errhigh_ent])

        xvalues_AGN.append(stat_AGN[oo][xkey][cen_ent])
        xerrlow_AGN.append(stat_AGN[oo][xkey][errlow_ent])
        xerrhigh_AGN.append(stat_AGN[oo][xkey][errhigh_ent])

        yvalues_AGN.append(stat_AGN[oo][ykey][cen_ent])
        yerrlow_AGN.append(stat_AGN[oo][ykey][errlow_ent])
        yerrhigh_AGN.append(stat_AGN[oo][ykey][errhigh_ent])

    xvalues_SF  = np.asarray(xvalues_SF )
    xvalues_AGN = np.asarray(xvalues_AGN)
    yvalues_SF  = np.asarray(yvalues_SF )
    yvalues_AGN = np.asarray(yvalues_AGN)

    fig = plt.figure(1, figsize=(6, 6))
    # definitions for the axes
    left, width = 0.15, 0.60
    bottom, height = 0.15, 0.60
    bottom_h = left_h = left + width + 0.01

    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=left, right=left+width, bottom=bottom, top=bottom+height)
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    Fsize    = 14
    lthick   = 2
    marksize = 6
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    plt.xlim(xrange)
    xminsys, xmaxsys = plt.xlim()

    plt.ylim(yrange)
    yminsys, ymaxsys = plt.ylim()

    for oo in np.arange(len(stat_SF)): # loop necessary for coloring and markers
        objid = stat_SF[oo]['id']
        mfc   = True
        if (objid < 6e8): # CDFS and COSMOS
            markersym   = 'o'
        elif (objid < 7e8) & (objid > 6e8): # UDF
            markersym   = 'D'
        elif (objid < 9e8) & (objid > 7e8): # UDF10
            markersym   = 'X'
        elif (objid > 1e9): # Literature objects
            markersym   = lce.get_reference_fromID(objid,verbose=False)[4]
            mfc         = False
        else:
            print(' WARNING - stopped as could not assing a marker symbol to the id '+str(objid))
            pdb.set_trace()

        ms          = marksize
        limsizefrac = 0.05


        # change color of limits
        markerzorder = 20

        if mfc:
            markerfacecolor = SFcol
        else:
            markerfacecolor = 'None'

        if datsetup == 'meanstd':
            plt.errorbar(xvalues_SF[oo],yvalues_SF[oo],xerr=xerrlow_SF[oo],yerr=yerrlow_SF[oo],capthick=0.5,
                         marker=markersym,lw=lthick/2., markersize=ms, alpha=0.5,
                         markerfacecolor=SFcol,ecolor=SFcol,
                         markeredgecolor=SFcol,zorder=markerzorder)
        else:
            plt.errorbar(xvalues_SF[oo],yvalues_SF[oo],
                         xerr=[[xvalues_SF[oo]-xerrlow_SF[oo],xerrhigh_SF[oo]-xvalues_SF[oo]]],
                         yerr=[[yvalues_SF[oo]-yerrlow_SF[oo],yerrhigh_SF[oo]-yvalues_SF[oo]]],capthick=0.5,
                         marker=markersym,lw=lthick/2., markersize=ms, alpha=0.5,
                         markerfacecolor=SFcol,ecolor=SFcol,
                         markeredgecolor=SFcol,zorder=markerzorder)


        if mfc:
            markerfacecolor = AGNcol
        else:
            markerfacecolor = 'None'

        if datsetup == 'meanstd':
            plt.errorbar(xvalues_AGN[oo],yvalues_AGN[oo],xerr=xerrlow_AGN[oo],yerr=yerrlow_AGN[oo],capthick=0.5,
                         marker=markersym,lw=lthick/2., markersize=ms, alpha=0.5,
                         markerfacecolor=AGNcol,ecolor=AGNcol,
                         markeredgecolor=AGNcol,zorder=markerzorder)
        else:
            plt.errorbar(xvalues_AGN[oo],yvalues_AGN[oo],
                         xerr=[[xvalues_AGN[oo]-xerrlow_AGN[oo],xerrhigh_AGN[oo]-xvalues_AGN[oo]]],
                         yerr=[[yvalues_AGN[oo]-yerrlow_AGN[oo],yerrhigh_AGN[oo]-yvalues_AGN[oo]]],capthick=0.5,
                         marker=markersym,lw=lthick/2., markersize=ms, alpha=0.5,
                         markerfacecolor=AGNcol,ecolor=AGNcol,
                         markeredgecolor=AGNcol,zorder=markerzorder)


    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if ylog:
        plt.yscale('log')
    if xlog:
        plt.xscale('log')

    # --------- HISTOGRAMS ---------
    Nbins   = 50
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    axHistx.xaxis.set_major_formatter(NullFormatter())
    axHisty.yaxis.set_major_formatter(NullFormatter())

    binwidth_x = np.diff([xminsys,xmaxsys])/Nbins
    bindefs    = np.arange(xminsys, xmaxsys+binwidth_x, binwidth_x)
    if xlog:
        bindefs = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))
        axHistx.set_xscale('log')

    axHistx.hist(xvalues_SF[np.isfinite(xvalues_SF)],  bins=bindefs,histtype='step',color=SFcol)
    axHistx.hist(xvalues_AGN[np.isfinite(xvalues_AGN)], bins=bindefs,histtype='step',color=AGNcol)
    axHistx.set_xticks([])
    axHistx.set_xlim([xminsys,xmaxsys])

    binwidth_y = np.diff([yminsys,ymaxsys])/Nbins
    bindefs    = np.arange(yminsys, ymaxsys+binwidth_y, binwidth_y)
    if ylog:
        bindefs = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))
        axHisty.set_yscale('log')

    axHisty.hist(yvalues_SF[np.isfinite(yvalues_SF)],  bins=bindefs,histtype='step',color=SFcol, orientation='horizontal')
    axHisty.hist(yvalues_AGN[np.isfinite(yvalues_AGN)], bins=bindefs,histtype='step',color=AGNcol, orientation='horizontal')
    axHisty.set_yticks([])
    axHisty.set_ylim([yminsys,ymaxsys])

    # if colorcode:
    #     cb      = plt.colorbar(m,extend=cextend,orientation='vertical',
    #                            pad=0.01,aspect=10,shrink=0.35,anchor=(-15.0,1.58),use_gridspec=False)
    #     cb.set_label(clabel)

    #--------- LEGEND ---------
    # plt.errorbar(-5000,-5000,xerr=None,yerr=1,marker='o',lw=0, markersize=marksize,alpha=1.0,
    #              markerfacecolor='k',ecolor='k',markeredgecolor='black',zorder=1,label='MUSE-Wide LAE')
    # plt.errorbar(-5000,-5000,xerr=None,yerr=None,marker='*',lw=0, markersize=marksize*2,alpha=1.0,
    #              markerfacecolor='None',ecolor='None',markeredgecolor='black',zorder=1,label='AGN')
    # plt.errorbar(-5000,-5000,xerr=None,yerr=None,marker='D',lw=0, markersize=marksize,alpha=1.0,
    #              markerfacecolor='None',ecolor='None',markeredgecolor='black',zorder=1,label='AGN candidate')
    #
    # leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize/1.0},ncol=5,numpoints=1,
    #                  bbox_to_anchor=(0.5, 1.1),)  # add the legend
    # leg.get_frame().set_alpha(0.7)
    #--------------------------

    if verbose: print('   Saving plot to '+plotname)
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_modelparametercollections(plotname, parametercollection_SF, parametercollection_AGN,
                                   stat_SF, stat_AGN, AGNcol='blue',SFcol='red',verbose=True):
    """

    --- INPUT ---
    plotname
    parametercollection_SF
    parametercollection_AGN
    stat_SF
    stat_AGN
    AGNcol
    SFcol
    verbose

    --- EXAMPLE OF USE ---

    """

    Nobj = len(parametercollection_SF)
    if verbose: print(' - Will generate plots of NEOGAL "PDFs" for all '+str(Nobj)+' objects in parameter collections')
    for oo in np.arange(len(parametercollection_SF)):
        objid        = parametercollection_SF[oo]['id']
        plotname_obj = plotname.replace('.pdf','_id'+str(objid)+'.pdf')
        if verbose: print(' - Generating the figure '+plotname_obj)
        figuresize_x = 6
        figuresize_y = 5
        fig          = plt.figure(figsize=(figuresize_x,figuresize_y))
        Fsize        = 9
        LW           = 2
        plt.rc('text', usetex=True)                         # enabling LaTex rendering of text
        plt.rc('font', family='serif',size=Fsize)           # setting text font
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()

        left   = 0.10   # the left side of the subplots of the figure
        right  = 0.95   # the right side of the subplots of the figure
        bottom = 0.10   # the bottom of the subplots of the figure
        top    = 0.93   # the top of the subplots of the figure
        wspace = 1.50   # the amount of width reserved for blank space between subplots
        hspace = 0.50   # the amount of height reserved for white space between subplots
        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        Nrows, Ncols      = 3, 6
        ylabel            = 'Number of NEOGAL SF ('+str(SFcol)+') and AGN ('+str(AGNcol)+') models'

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Nmodels_SF  = float(len(parametercollection_SF[oo]['Zgas']))
        Nmodels_AGN = float(len(parametercollection_AGN[oo]['Zgas']))

        titlestr = 'Models satisfying ID='+str(objid)+' cuts: SF='+str(Nmodels_SF)+'; AGN='+str(Nmodels_AGN)
        if (Nmodels_AGN > 0) & (Nmodels_SF > 0):
            Nmodels_ratio = Nmodels_SF/Nmodels_AGN
            titlestr_addition = '; SF/AGN='+str("%.4f" % Nmodels_ratio)
        else:
            titlestr_addition = ' '
            Nmodels_ratio = np.nan
        fig.suptitle(titlestr+titlestr_addition)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Zgas
        plt.subplot(Nrows, Ncols, (1,3))

        bindefs = np.array([0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.014,
                            0.017, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07])-0.00001

        plotkey = 'Zgas'
        nm.plot_modelparametercollections_addhist(parametercollection_SF[oo][plotkey],parametercollection_AGN[oo][plotkey],
                                                  stat_SF[oo][plotkey],stat_AGN[oo][plotkey],
                                                  SFcol,AGNcol,LW,bindefs=bindefs,Nbins=None)

        plt.xscale('log')
        plt.xlabel('Z')
        plt.xlim([0.00001,0.1])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # logUs
        plt.subplot(Nrows, Ncols, (4,6))

        bindefs = np.arange(-4.75, -0.25, 0.5)

        plotkey = 'logUs'
        nm.plot_modelparametercollections_addhist(parametercollection_SF[oo][plotkey],parametercollection_AGN[oo][plotkey],
                                                  stat_SF[oo][plotkey],stat_AGN[oo][plotkey],
                                                  SFcol,AGNcol,LW,bindefs=bindefs,Nbins=None)

        plt.xlabel('log(U)')
        plt.xlim([-5,-0.5])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # xid
        plt.subplot(Nrows, Ncols, (7,9))

        bindefs = np.array([0.0, 0.2, 0.4, 0.6])

        plotkey = 'xid'
        nm.plot_modelparametercollections_addhist(parametercollection_SF[oo][plotkey],parametercollection_AGN[oo][plotkey],
                                                  stat_SF[oo][plotkey],stat_AGN[oo][plotkey],
                                                  SFcol,AGNcol,LW,bindefs=bindefs,Nbins=None)

        plt.xlabel('$\\xi_\\textrm{d}$')
        plt.xlim([-0.05,0.65])
        plt.ylabel(ylabel)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # nh
        plt.subplot(Nrows, Ncols, (10,12))

        bindefs = 10**np.array([1.5, 2.5, 3.5, 4.5])

        plotkey = 'nh'
        nm.plot_modelparametercollections_addhist(parametercollection_SF[oo][plotkey],parametercollection_AGN[oo][plotkey],
                                                  stat_SF[oo][plotkey],stat_AGN[oo][plotkey],
                                                  SFcol,AGNcol,LW,bindefs=bindefs,Nbins=None)
        plt.xscale('log')
        plt.xlabel('$n_\\textrm{H}$ [cm$^{-3}$]')
        plt.xlim([10,1e5])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # COCOsol
        plt.subplot(Nrows, Ncols, (13,14))

        #bindefs = np.array([0.10, 0.14, 0.20, 0.27, 0.38, 0.52, 0.72, 1.00, 1.40])
        bindefs = np.arange(0.05,1.5,0.06)

        plotkey = 'COCOsol'
        nm.plot_modelparametercollections_addhist(parametercollection_SF[oo][plotkey],None,
                                                  stat_SF[oo][plotkey],None,
                                                  SFcol,AGNcol,LW,bindefs=bindefs,Nbins=None)

        plt.xlabel('C/O [(C/O)$_\\textrm{sun}$]')
        plt.xlim([0.00,1.55])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # mup
        plt.subplot(Nrows, Ncols, (15,16))

        bindefs = np.array([0,200,400])

        plotkey = 'mup'
        nm.plot_modelparametercollections_addhist(parametercollection_SF[oo][plotkey],None,
                                                  stat_SF[oo][plotkey],None,
                                                  SFcol,AGNcol,LW,bindefs=bindefs,Nbins=None)

        plt.xlabel('m$_\\textrm{up}$ [M$_\\textrm{sun}$]')
        plt.xlim([-10,410])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # alpha
        plt.subplot(Nrows, Ncols, (17,18))

        bindefs = np.array([-2.15,-1.85,-1.55,-1.25,-0.95])

        plotkey = 'alpha'
        Nbins   = 10
        nm.plot_modelparametercollections_addhist(None,parametercollection_AGN[oo][plotkey],
                                                  None,stat_AGN[oo][plotkey],
                                                  SFcol,AGNcol,LW,bindefs=None,Nbins=Nbins)

        plt.xlabel('$\\alpha$')
        plt.xlim([-2.2,-0.9])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plt.savefig(plotname_obj)
        plt.clf()
        plt.close('all')
        if verbose: print(' - Successfully saved figure to file')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_modelparametercollections_addhist(paramcol_SF,paramcol_AGN,stat_SF,stat_AGN,
                                           SFcol,AGNcol,LW,bindefs=None,Nbins=10):
    """

    --- INPUT ---
    parametercollection_SF
    parametercollection_AGN
    stat_SF
    stat_AGN
    SFcol
    AGNcol
    LW
    bindefs
    Nbins

    --- EXAMPLE OF USE ---
    see nm.plot_modelparametercollections() above

    """
    if (paramcol_SF is None) & (paramcol_AGN is not None):
        if len(paramcol_AGN) > 0:
            if bindefs is None:
                xmin       = np.min(paramcol_AGN)*0.95
                xmax       = np.max(paramcol_AGN)*1.05
                binwidth_x = np.diff([xmin,xmax])/Nbins
                bindefs    = np.arange(xmin, xmax+binwidth_x, binwidth_x)

            plt.hist(paramcol_AGN, bins=bindefs,histtype='step',color=AGNcol)

            yminsys, ymaxsys = plt.ylim()

            plt.plot([stat_AGN[0]]*2, [yminsys, ymaxsys], '--',lw=LW,alpha=0.5,color=AGNcol)  # mean
            plt.plot([stat_AGN[2]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=AGNcol)  # median
            plt.plot([stat_AGN[4]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=AGNcol)  # 68% confidence interval lower
            plt.plot([stat_AGN[8]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=AGNcol)  # 68% confidence interval lower

    elif (paramcol_SF is not None) & (paramcol_AGN is None):
        if len(paramcol_SF) > 0:

            if bindefs is None:
                xmin       = np.min(paramcol_SF)*0.95
                xmax       = np.max(paramcol_SF)*1.05
                binwidth_x = np.diff([xmin,xmax])/Nbins
                bindefs    = np.arange(xmin, xmax+binwidth_x, binwidth_x)

            plt.hist(paramcol_SF,  bins=bindefs,histtype='step',color=SFcol)

            yminsys, ymaxsys = plt.ylim()

            plt.plot([stat_SF[0]]*2,  [yminsys, ymaxsys], '--',lw=LW,alpha=0.5,color=SFcol)  # mean
            plt.plot([stat_SF[2]]*2,  [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=SFcol)  # median
            plt.plot([stat_SF[4]]*2,  [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=SFcol)  # 68% confidence interval lower
            plt.plot([stat_SF[8]]*2,  [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=SFcol)  # 68% confidence interval lower

    else:

        if (len(paramcol_SF) == 0) & (len(paramcol_AGN) == 0):
            pass
        else:
            if bindefs is None:
                xmin       = np.min(np.append(paramcol_SF,paramcol_AGN))*0.95
                xmax       = np.max(np.append(paramcol_SF,paramcol_AGN))*1.05
                binwidth_x = np.diff([xmin,xmax])/Nbins
                bindefs    = np.arange(xmin, xmax+binwidth_x, binwidth_x)
                # if xlog:
                #     bindefs = np.logspace(np.log10(bindefs[0]),np.log10(bindefs[-1]),len(bindefs))
                #     axHistax.set_xscale('log')

            plt.hist(paramcol_SF,  bins=bindefs,histtype='step',color=SFcol)
            plt.hist(paramcol_AGN, bins=bindefs,histtype='step',color=AGNcol)

            yminsys, ymaxsys = plt.ylim()

            if len(paramcol_SF) > 0:
                plt.plot([stat_SF[0]]*2,  [yminsys, ymaxsys], '--',lw=LW,alpha=0.5,color=SFcol)  # mean
                plt.plot([stat_SF[2]]*2,  [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=SFcol)  # median
                plt.plot([stat_SF[4]]*2,  [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=SFcol)  # 68% confidence interval lower
                plt.plot([stat_SF[8]]*2,  [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=SFcol)  # 68% confidence interval lower

            if len(paramcol_AGN) > 0:
                plt.plot([stat_AGN[0]]*2, [yminsys, ymaxsys], '--',lw=LW,alpha=0.5,color=AGNcol)  # mean
                plt.plot([stat_AGN[2]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=AGNcol)  # median
                plt.plot([stat_AGN[4]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=AGNcol)  # 68% confidence interval lower
                plt.plot([stat_AGN[8]]*2, [yminsys, ymaxsys], ':', lw=LW,alpha=0.5,color=AGNcol)  # 68% confidence interval lower

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =