#import astropy
#import scipy
from importlib import reload
import pdb
import matplotlib.pylab as plt
import numpy as np

def test(verbose=True):
    print('hello world')
    for ii, number in enumerate(np.array([10, 20, 30, 40])):
        print('This is just a test '+str(number)+'  and then the index '+str(ii))

    verbose = True

    xvalues = np.arange(1, 100, 1.0)
    yvalues = np.sin(xvalues)
    yerr    = np.ones(len(xvalues))*0.2 # np.random(len(xvalues)) # np.sqrt(xvalues)
    xerr    = np.asarray([None] * len(xvalues))

    # - - - - - - - - - - - - - - - - - - - - - - PLOTTING - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Setting up and generating plot')
    outdir   = '/Users/kaschm/plots/'
    plotname = outdir+'Testplot.pdf'
    fig = plt.figure(figsize=(7, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.2, right=0.97, bottom=0.10, top=0.9)
    Fsize = 10
    lthick = 2
    marksize = 4
    plt.rc('text', usetex=False)
    plt.rc('font', family='serif', size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    # plt.title(inforstr[:-2],fontsize=Fsize)

    plt.errorbar(xvalues, yvalues, xerr=None, yerr=yerr,
                 marker='o', lw=1, markersize=marksize, alpha=0.5,
                 markerfacecolor='gray', ecolor='k',
                 markeredgecolor='k', zorder=10)

    plt.xlabel(' these are the xvalues ')
    plt.ylabel(' and this is the sine curve values æåø')

    # --------- RANGES ---------
    goodent = np.where(np.isfinite(xvalues) & np.isfinite(yvalues))
    xmin = np.min(xvalues[goodent])
    xmax = np.max(xvalues[goodent])
    dx = xmax - xmin

    ymin = np.min(yvalues[goodent])
    ymax = np.max(yvalues[goodent])
    dy = ymax - ymin

    plt.xlim([xmin - dx * 0.05, xmax + dx * 0.05])
    plt.ylim([ymin - dy * 0.05, ymax + dy * 0.05])

    # if logx:
    #     plt.xscale('log')
    # if logy:
    #     plt.yscale('log')

    #leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.7}, ncol=3, numpoints=1,
    #                 bbox_to_anchor=(0.5, 1.1), )  # add the legend
    #leg.get_frame().set_alpha(0.7)
    # --------------------------

    saveplot = True
    if saveplot:
        if verbose: print('   Saving plot to', plotname)
        plt.savefig(plotname)
    else:
        plt.show()

    plt.clf()
    plt.close('all')


