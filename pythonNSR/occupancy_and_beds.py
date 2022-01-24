#from importlib import reload
import pdb
import sys
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
import datetime

#sys.path.append('C:/Users/kaschm/GitHub/python/pythonNSR/')
import occupancy_and_beds as oab


# -----------------------------------------------------------------------------------------------------------------------
def generate_datastructure(verbose=True):
    """
    Function to generate the pivor data structure to use for plotting
    Example of use
    -------
    import occupancy_and_beds as oab
    dataframe = oab.generate_datastructure()

    """
    if verbose: print(' - Loading data from excel sheet Belægningshistorik_alle_afdelinger_2018til2022.xlsx')
    filepath  = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Belægningsprocenter/'
    filename  = 'Belægningshistorik_alle_afdelinger_2018til2022.xlsx'
    dataframe = pd.read_excel(filepath+filename, sheet_name='Belægningsoversigt')

    return dataframe
# -----------------------------------------------------------------------------------------------------------------------
def plot_beds(datemin='01-01-2018',datemax='01-01-2022',hour2show=23,SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'] ,verbose=True):
    """

    Parameters
    ----------
    datemin
    datemax
    hour2show
    verbose

    Example of use
    -------
    import occupancy_and_beds as oab
    oab.plot_beds(SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'])
    oab.plot_beds(SORsections=['SJ NAEOKI8, ORTOPÆDKIR. SENGEAFSNIT 8, NAE', 'SJ SLAOKIOT, O-KIR. SENGEAFS., TRAUME, SLA', 'SJ SLAANINT, INTENSIV SENGEAFS., SLA', 'SJ SLAANITM, INTERMEDIÆRT SENGEAFS., SLA'])
    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dataframe = oab.generate_datastructure(verbose=verbose)
    datemin = datetime.datetime.strptime(datemin, "%d-%m-%Y")
    datemax = datetime.datetime.strptime(datemax, "%d-%m-%Y")

    ptable = pd.pivot_table(dataframe, values='Belægning disponible senge', index=['Belægning tidspunkt'],
                            columns=['Ophold afsnit navn'], aggfunc=np.mean)

    if SORsections == 'all':
        SORsectionlist = ptable.columns.tolist()
    else:
        SORsectionlist = SORsections
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'disponiblesenge.pdf'
    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(15, 6))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.05, right=0.97, bottom=0.18,top=0.80)
    Fsize = 10
    lthick = 2
    marksize = 4

    plt.clf()
    plt.ioff()

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    #plt.gca().xaxis.set_major_locator(mdates.YearLocator())
    #plt.gcf().autofmt_xdate()  # set font and rotation for date tick labels
    plt.xticks(rotation=45)

    #ax = ptable.plot(kind='line', lw=2, colormap='viridis', zorder=10., alpha=0.5, figsize=(15, 6), fontsize=Fsize)
    #plt.subplots_adjust(wspace=0.1, hspace=0.1, left=0.05, right=0.97, bottom=0.02, top=0.80)

    plt.rc('text', usetex=False)
    plt.rc('font', family='serif', size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    # plt.title(inforstr[:-2],fontsize=Fsize)

    xerr = None
    yerr = None

    plt.grid(linestyle=':', linewidth=lthick/2.)
    # --------- RANGES ---------
    ymin = 0
    ymax = 50.0
    dy   = ymax - ymin
    plt.ylim([ymin,ymax])

    #plt.xscale('log')
    #plt.yscale('log')

    # --------- COLORMAP ---------
    cmap = plt.cm.get_cmap('viridis')
    cmin = 0
    cmax = 100

    colnorm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, num=50)
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(cmaparr)
    #cb = plt.colorbar(m)
    #cb.set_label('belægningsprocent')

    # --------- POINT AND CURVES ---------
    for bs, SORsection in enumerate(SORsectionlist):
        if len(SORsectionlist) == 1:
             pointcolor = cmap(colnorm(45))
        else:
            pointcolor = cmap(colnorm((cmax-cmin)/len(SORsectionlist)*(bs+1)))
        values2plot = ptable[SORsection][np.where((ptable[SORsection].index.hour == hour2show) &
                                                  (ptable[SORsection].index > datemin) &
                                                  (ptable[SORsection].index < datemax))[0]]

        plt.errorbar(values2plot.index, values2plot.values, xerr=xerr, yerr=yerr,
                     marker='.', lw=lthick, markersize=0, alpha=0.5, color=pointcolor,
                     markerfacecolor=pointcolor, ecolor=pointcolor,
                     markeredgecolor=pointcolor, zorder=18.,
                     label=SORsection)



    lineymin = ymax - dy*0.40
    lineymax = ymax - dy*0.28
    textymin = ymax - dy*0.25
    textymax = ymin - dy*0.10
    if 'SJ NAELUIN, LUNGEMED. SENGEAFS., NAE' in SORsectionlist:
        plt.plot([datetime.datetime.strptime("10-06-2021", "%d-%m-%Y"), datetime.datetime.strptime("10-06-2021", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
        plt.text(datetime.datetime.strptime("10-06-2021", "%d-%m-%Y"), textymin, '10-06-2021', fontsize=Fsize, rotation=90, color='gray',
                 horizontalalignment='center', verticalalignment='bottom')

        plt.plot([datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"), datetime.datetime.strptime("01-03-2021", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
        plt.text(datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"), textymin, '01-03-2021', fontsize=Fsize, rotation=90, color='gray',
                 horizontalalignment='center', verticalalignment='bottom')

    # --------- LABELS ---------
    plt.xlabel('Dato')
    plt.ylabel('Antal disponible senge')

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=3, numpoints=1,
                     bbox_to_anchor=(0.5, 1.32), )  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------
    plotdir   = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Belægningsprocenter/'
    outputfig = plotdir+plotname
    plt.savefig(outputfig)
    plt.clf()
    plt.close('all')
    if verbose: print(' - Saved figure to \n  '+outputfig)
#=======================================================================================================================

