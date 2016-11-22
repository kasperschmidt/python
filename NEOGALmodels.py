# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Functions and scripts to analyze and deal with NEOGAL ionization models
# http://www.iap.fr/neogal/models.html
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import numpy as np
import glob
import pyfits
import sys
import pdb
import crossmatch as cm
import matplotlib.pyplot as plt
import matplotlib
import NEOGALmodels as nm
import itertools
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

    """
    if Zgas == 'combined':
        Zgasstring = Zgas
    else:
        Zgasstring = str(Zgas).split('.')[-1]
    filename   = filepath+'nebular_emission_Z'+Zgasstring+'.txt'
    if verbose: print ' - Attempting to load data from:\n   '+filename
    modeldata  = np.genfromtxt(filename,names=True,dtype=None)
    if verbose: print '   ... successful'
    return modeldata
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_LvsL(modeldata,line1='CIV1548',line2='CIII1908',plotname='./TESTPLOT.pdf',
              Zgas=False,logU=False,xid=0.3,nh=100,COratio=0.38,Mcutoff=100,
              logx=False,logy=False,logp1=False,logp2=False,fixxrange=False,fixyrange=False,verbose=True):
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
    verbose   Toggle verbosity

    --- EXAMPLE OF USE ---
    import NEOGALmodels as nm
    modeldata = nm.load_model('combined',verbose=True)
    nm.plot_LvsL(modeldata,line1='CIV1548',line2='CIII1908',Zgas=False,logU=False,xid=0.3,nh=100,COratio=0.38,Mcutoff=100,logx=True,logy=True,logp1=True)

    nm.plot_LvsL(modeldata,line1='CIV1548',line2='CIII1908',Zgas=False,logU=False,xid=0.1,nh=10,COratio=0.1,Mcutoff=100,logx=True,logy=True,logp1=True)

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
        inforstr      = inforstr+' Mcutoff:vary, '
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
        if verbose: print ' - Getting data for '+str(Ngoodent)+' data points satisfying model selection '
        param1 = modeldata[freeparam[0]][goodent]
        if logp1:
            param1 = np.log10(param1)

        param2 = modeldata[freeparam[1]][goodent]
        if logp2:
            param2 = np.log10(param2)

        lum1   = modeldata[line1][goodent]
        lum2   = modeldata[line2][goodent]
    else:
        if verbose: print ' WARNING: Less than 2 model grid points to plot; no output generated'
        return

    # - - - - - - - - - - - PLOTTING - - - - - - - - - - -
    if verbose: print ' - Setting up and generating plot'
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

    plt.xlabel(r'L$_\textrm{'+line1+'}$ / [1e33 erg/s]/[Msun/yr]')
    plt.ylabel(r'L$_\textrm{'+line2+'}$ / [1e33 erg/s]/[Msun/yr]')

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

    if verbose: print '   Saving plot to',plotname
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

    if verbose: print ' - Putting together permutations of chosen setups for plotting'
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

    if verbose: print ' - With the restriction Nfalse=2 the setup will results in '+str(Nplots)+\
                      ' plots (if model data allows)'
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
            infostr = ' - Generating plot '+str("%.4d" % plotno)+'/'+str("%.4d" % Nplots)+': '+plotname+'  '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        if not Zval:
            logp1 = True
        else:
            logp1 = False

        nm.plot_LvsL(modeldata,line1=line1,line2=line2,logx=True,logy=True,logp1=logp1,logp2=False,verbose=False,
                     Zgas=Zval,logU=Uval,xid=Xival,nh=Nhval,COratio=COval,Mcutoff=Mval,
                     fixxrange=line1range,fixyrange=line2range,plotname=plotname)

    print '\n   ... done'
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def combine_modeloutputs(outputname='nebular_emission_Zcombined.txt',
                         filepath='/Users/kschmidt/work/catalogs/NEOGALlines/nebular_emission/',
                         verbose=True):
    """
    Combine the model outputs to have a single 'master-model' with all variables included

    --- EXAMPLE OF USE ---
    import NEOGALmodels as nm
    nm.combine_modeloutputs()

    """
    output     = filepath+outputname
    if verbose: print ' - Setting up output for:\n   '+output
    modelfiles = glob.glob(filepath+'nebular_emission_Z0*.txt')
    header     = open(modelfiles[0]).readline().rstrip()
    header     = header.replace('##','##  Zgas  ')

    fout = open(output, 'w')
    fout.write(header)
    if verbose: print ' - Writing the following files to ouput:'
    for mf in modelfiles:
        if verbose: print '   '+mf
        Zgasstring   = mf.split('/')[-1].split('emission_Z')[-1].split('.txt')[0]

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
def convert_Fline2Lbol(Fline,redshift,verbose=True):
    """
    Converting an observed integrated line flux [erg/s/cm2] to bolometric luminoisity [erg/s]
    """

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def convert_Lbol2Fline(Lbol,redshift,verbose=True):
    """
    Converting bolometric luminoisity [erg/s] to integrated line flux [erg/s/cm2]
    """

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =