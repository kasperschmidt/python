# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Scripts for manipulating template spectra, e.g., from the directory /Users/kschmidt/work/templatespecs/
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import numpy as np
import glob
import matplotlib.pyplot as plt
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_BC03templatespecs(xrange=[0,4e4],yrange=[0.01,100],metalonly=False,interactive=False,verbose=True):
    """
    Plot BC03 template specs together

    --- EXAMPLE OF USE ---
    import templatespecs as ts
    outname = ts.plot_BC03templatespecs(metalonly='02',interactive=True)

    """
    specdir  = '/Users/kschmidt/work/templatespecs/bc03/templates/'
    if verbose: print ' - Globbing for spectra in '+specdir
    speclist = glob.glob(specdir+'*spec')
    outname  = specdir+'BC03_templatespecs.pdf'
    if verbose: print '   Found '+str(len(speclist))+' to plot'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Setting up colors and linestyles'

    widthdic   = {'11Gyr':3.0,
                  '5Gyr':2.7,
                  '2.5Gyr':2.4,
                  '1.4Gyr':2.1,
                  '900Myr':1.8,
                  '640Myr':1.5,
                  '290Myr':1.2,
                  '100Myr':0.9,
                  '25Myr':0.6,
                  '5Myr':0.3,
                  '6gyr':1,
                  '12gyr':1,}

    coldic    = {'cst':'red','ssp':'black','t5e9':'green','t9e9':'blue'}
    linedic   = {'z008':'-','z02':'--','z05':'-.'}

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Plotting the BC03 template spectra to \n   '+outname

    fig = plt.figure(figsize=(10, 4.3))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.06, right=0.63, bottom=0.1, top=0.95)
    Fsize  = 10
    lthick = 1
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    if interactive:
        plt.ion()
    else:
        plt.ioff()
    #plt.title(plotname.split('/')[-1].replace('_','\_'),fontsize=Fsize)

    #--------------------------------------------------------------------------------------
    for ss, spec in enumerate(speclist):
        sname     = spec.split('/')[-1].split('.spe')[0]
        sdat      = np.genfromtxt(spec,comments='#')
        swave     = sdat[:,0]
        sflux     = sdat[:,1]

        colstr    = sname.split('_')[0]
        widthstr  = sname.split('_')[1]
        linestr   = sname.split('_')[2]
        scol      = coldic[colstr]
        sline     = linedic[linestr]
        swidth    = widthdic[widthstr]

        if metalonly == False:
            plt.plot(swave,sflux,color=scol,ls=sline,alpha=0.7,lw=swidth,label=sname.replace('_',', '))
        elif metalonly in sname:
            plt.plot(swave,sflux,color=scol,ls='-',alpha=0.7,lw=swidth,label=sname.replace('_',', '))
    #--------------------------------------------------------------------------------------

    plt.xlabel('Wavelength [\AA] ', fontsize=Fsize)
    plt.ylabel('log(Flux)', fontsize=Fsize)

    plt.xlim(xrange)
    plt.ylim(yrange)
    plt.yscale('log')

    leg = plt.legend(fancybox=True, loc='upper right',prop={'size':Fsize}, ncol=2, numpoints=1,
                     bbox_to_anchor=(1.65, 1.04))  # add the legend
    leg.get_frame().set_alpha(0.7)

    if not interactive:
        if verbose: print '   Saving plot to',outname
        plt.savefig(outname)
        plt.clf()
        plt.close('all')

    return outname
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =