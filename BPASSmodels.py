# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Functions and scripts to analyze and deal with BPASS ionization models
# https://bpass.auckland.ac.nz/4.html (Xiao, Stanway and Eldridge 2018)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import numpy as np
import glob
import astropy.io.fits as afits
import astropy.io.ascii as aascii
from astropy.table import Table, vstack
import sys
import pdb
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import matplotlib
import BPASSmodels as bm
import itertools
import literaturecollection_emissionlinestrengths as lce
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def load_model(Zgas,filepath='/Users/kschmidt/work/catalogs/BPASSbasedNebularEmission/UV/',binaries=True,verbose=True):
    """
    Loat the content of a model output

    --- INPUT ---
    Zgas       The interstellar gas metallicity of the model to load. If set to 'combined'
               the combined model output [generated with bm.combine_modeloutputs()] will be loaded.
    filepath   The directory to load models from
    binaries   Load the binaries models (True) or the single star models (False)
    verbose    Toggle verbosity


    --- EXAMPLE OF USE ---
    import BPASSmodels as bm
    models_bin_uv, filename = bm.load_model(0.001,filepath='/Users/kschmidt/work/catalogs/BPASSbasedNebularEmission/UV/',binaries=True)
    models_bin_op, filename = bm.load_model(0.001,filepath='/Users/kschmidt/work/catalogs/BPASSbasedNebularEmission/Optical/',binaries=True)

    models_bin_uv, filename = bm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/BPASSbasedNebularEmission/UV/',binaries=True)
    models_bin_op, filename = bm.load_model('combined',filepath='/Users/kschmidt/work/catalogs/BPASSbasedNebularEmission/optical/',binaries=True)

    """
    if binaries:
        typestr = 'bin'
    else:
        typestr = 'sin'

    if Zgas == 'combined':
        Zgasstring = Zgas
    elif (Zgas == 0.0001):
        Zgasstring = 'em4'
    elif (Zgas == 0.00001):
        Zgasstring = 'em5'
    else:
        Zgasstring = str("%.3f" % Zgas).split('.')[-1]
    if 'UV' in filepath:
        colnames   = ['No','logU','lognH','logAge','HeII1640','EW_HeII1640','CIII1907','EW_CIII1907','CIII1910','EW_CIII1910',
                      'CIV1548','EW_CIV1548','CIV1551','EW_CIV1551','OI1357','EW_OI1357','OIII1661','EW_OIII1661',
                      'OIII1666','EW_OIII1666','SiII1263','EW_SiII1263','SiII1308','EW_SiII1308','SiII1531','EW_SiII1531']

        filename   = filepath+typestr+'/UV_data_'+typestr+'_z'+Zgasstring+'.dat'

    elif 'Optical' in filepath:
        colnames   =  ['No','logU','lognH','logAge','NII6548','EW_NII6548','NII6584','EW_NII6584','SII6716','EW_SII6716',
                       'SII6731','EW_SII6731','OI6300','EW_OI6300','OIII4959','EW_OIII4959','OIII5007','EW_OIII5007',
                       'Ha6563','EW_Ha6563','Hb4861','EW_Hb4861','HeII4686','EW_HeII4686']
        filename   = filepath+typestr+'/Optical_data_'+typestr+'_z'+Zgasstring+'.dat'

    else:
        sys.exit('Did not find "UV" or "Optical" in file path so cannot figure out what format to use')

    if verbose: print(' - Attempting to load data from:\n   '+filename)
    modeldata  = aascii.read(filename,names=colnames)
    Ncol  = str(len(colnames))
    Nrows = str(len(modeldata['No']))
    if verbose: print('   ... successfully read '+Ncol+' columns and '+Nrows+' rows')
    return modeldata, filename
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def linenames():
    """
    Return dictionary that can convert BPASS column names to flux catalog LINENAMEs

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
def combine_modeloutputs(outputname='xxRENAMExx_binary_Zcombined.txt',
                         datatype='bin',verbose=True):
    """
    Combine the model outputs to have a single 'master-model' with all variables included

    --- EXAMPLE OF USE ---
    import BPASSmodels as nm
    bm.combine_modeloutputs(outputname='nebular_emission_BPASS_binaries_Zcombined.txt',datatype='bin')
    bm.combine_modeloutputs(outputname='nebular_emission_BPASS_singles_Zcombined.txt',datatype='sin')

    """
    outputname   = '/Users/kschmidt/work/catalogs/BPASSbasedNebularEmission/'+outputname
    zvals        = [0.040,0.030,0.020,0.014,0.010,0.008,0.006,0.004,0.003,0.002,0.001,0.0001,0.00001]

    if datatype == 'bin':
        binaries = True
    elif datatype == 'sin':
        binaries = False

    outputtable = None
    filelist    = []
    for zval in zvals:
        models_bin_uv, fileloaded_uv = bm.load_model(zval,filepath='/Users/kschmidt/work/catalogs/BPASSbasedNebularEmission/UV/',
                                                     binaries=binaries)
        filelist.append(fileloaded_uv)
        models_bin_op, fileloaded_op = bm.load_model(zval,filepath='/Users/kschmidt/work/catalogs/BPASSbasedNebularEmission/Optical/',
                                                     binaries=binaries)
        filelist.append(fileloaded_op)

        col_Zval = Table.Column(name='Zgas', data=models_bin_uv['logU']*0.0+zval)
        models_bin_uv.add_column(col_Zval,index=0)

        for cn in models_bin_op.colnames:
            if cn in ['No','logU','lognH','logAge']:
                coldiff  = models_bin_op[cn]-models_bin_uv[cn]
                Ndiffobj = len(np.where(coldiff != 0.0)[0])
                if Ndiffobj > 0:
                    print('\n The column '+cn+' does not match for the UV and Optical fields for Z='+str(zval)+'; stopping for investigation.\n')
                    pdb.set_trace()
            else:
                coldat = Table.Column(name=cn, data=models_bin_op[cn])
                models_bin_uv.add_column(coldat)

        if outputtable is None:
            outputtable = models_bin_uv
        else:
            outputtable = vstack([outputtable, models_bin_uv])

    if verbose: print(' - Writing output to:\n   '+outputname)
    aascii.write(outputtable, output=outputname, format='commented_header')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
