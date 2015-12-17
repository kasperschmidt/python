# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import numpy as np
import matplotlib.pyplot as plt
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_zVSmag(objectcatalog,verbose=True):
    """
    Plotting high redshift sources in comparison plot

    --- EXAMPLE OF USE ---
    import plotRedshift10objects as pzo
    objcat = '/Users/kasperborelloschmidt/work/BoRG/redshiftGT8sources.txt'
    pzo.plot_zVSmag(objcat)

    """
    plotname = objectcatalog.replace('.txt','_zVSmag.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading redshift catalog \n   '+objectcatalog
    zdat = np.genfromtxt(objectcatalog,dtype=None,skip_header=2,names=True,comments='#')

    refsymbols = {}
    refsymbols['BE16'] = ['o','black'   ,'This work']
    refsymbols['CA15'] = ['^','green'   ,'Calvi et al. (2015)']
    refsymbols['BO15'] = ['s','red'     ,'Bouwens et al. (2015)']
    refsymbols['HO15'] = ['h','blue'    ,'Holwerda et al. (2015)']
    refsymbols['RB15'] = ['d','magenta' ,'Roberts-Borsani et al. (2015)']
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Plotting objects \n   '+plotname
    fig = plt.figure(figsize=(6, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.97, bottom=0.07, top=0.85)

    Fsize  = 12
    lthick = 2
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    labels_all = np.array([])
    for object in zdat:
        objmark  = refsymbols[object['ref']][0]
        objcol   = refsymbols[object['ref']][1]
        objlab   = refsymbols[object['ref']][2]
        if objlab in labels_all:
            objlab = None
        else:
            labels_all = np.append(labels_all, objlab)

        objz     = object['z_phot']
        objzerrp = object['z_phot_errp']
        objzerrm = object['z_phot_errm']
        if objzerrp == -99:
            zerr = None
        else:
            if objzerrm == objzerrp:
                zerr = objzerrm
            else:
                zerr = [[objzerrm,objzerrp]]

        mag      = object['F160W']
        magerrp  = object['F160Werr']
        magerrm  = object['F160Werr']
        if magerrp == -99:
            magerr = None
        else:
            if magerrm == magerrp:
                magerr = magerrm
            else:
                magerr = [[magerrm,magerrp]]

        size     = object['r_e']
        objs = size * 100
        # adjust marker sizes
        if objs < 0.0  : objmark = None
        if objs > 0.0  : objmark_size = 3
        if objs > 2.0  : objmark_size = 5
        if objs > 20.0 : objmark_size = 7
        if objs > 40.0 : objmark_size = 9
        if objs > 60.0 : objmark_size = 11
        if objs > 80.0 : objmark_size = 13
        if objs > 100.0: objmark_size = 15

        sizeerr  = object['r_eerr']

        try:
            plt.errorbar(objz,mag,xerr=zerr,yerr=magerr,
                         marker=objmark,ms=objmark_size,capsize=0,
                         color=objcol,ecolor=objcol,lw=lthick-1)
        except:
            pdb.set_trace()

    #create labels
    for key in refsymbols.keys():
        lab = refsymbols[key]
        plt.plot(-10,-10,'.',marker=lab[0],ms=8,color=lab[1],label=lab[2])

    #plt.title('Chi2 curve from Gaussian Convolution')
    plt.xlabel('$z$',  fontsize=Fsize)
    plt.ylabel('$H_{160}$', fontsize=Fsize)
    plt.xlim([8.5,11.5])
    plt.ylim([30.5,24])

    leg = plt.legend(fancybox=True,numpoints=1, loc='lower right',prop={'size':Fsize-2},
                     ncol=2,bbox_to_anchor=(0.94, 1.0))  # add the legend
    leg.get_frame().set_alpha(0.7)
    plt.savefig(plotname)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =