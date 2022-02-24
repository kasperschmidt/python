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
def generate_datastructure(filename  = 'Belægningshistorik_alle_afdelinger_220124.xlsx',verbose=True):
    """
    Function to generate the pivor data structure to use for plotting
    Example of use
    -------
    import occupancy_and_beds as oab
    dataframe = oab.generate_datastructure()

    """
    if verbose: print(' - Loading data from excel sheet '+filename)
    filepath  = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Belægningsprocenter/'
    #filename  = 'Belægningshistorik_alle_afdelinger_2018til2022.xlsx'
    dataframe = pd.read_excel(filepath+filename, sheet_name='Belægningsoversigt')
    if verbose: print('   ... done')
    return dataframe
# -----------------------------------------------------------------------------------------------------------------------
def plot_beds(dataframe,datemin='01-01-2018',datemax='01-02-2022', plotname='disponiblesenge.pdf', hour2show=23,SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'] ,verbose=True):
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
    dataframe = oab.generate_datastructure(verbose=verbose)
    oab.plot_beds(dataframe,SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'])
    oab.plot_beds(dataframe,SORsections=['SJ NAEOKI8, ORTOPÆDKIR. SENGEAFSNIT 8, NAE', 'SJ SLAOKIOT, O-KIR. SENGEAFS., TRAUME, SLA', 'SJ SLAANINT, INTENSIV SENGEAFS., SLA', 'SJ SLAANITM, INTERMEDIÆRT SENGEAFS., SLA'])
    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datemin = datetime.datetime.strptime(datemin, "%d-%m-%Y")
    datemax = datetime.datetime.strptime(datemax, "%d-%m-%Y")

    ptable = pd.pivot_table(dataframe, values='Belægning disponible senge', index=['Belægning tidspunkt'],
                            columns=['Ophold afsnit navn'], aggfunc=np.mean, dropna=False)

    if SORsections == 'all':
        SORsectionlist = ptable.columns.tolist()
    else:
        SORsectionlist = SORsections
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
        elif len(SORsectionlist) == 2:
             pointcolor = [cmap(colnorm(30)), cmap(colnorm(80))][bs]
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


        if SORsections != 'all':
            legscale   = 1.0
            ncol       = 3
            diffvals   = np.diff(values2plot.values)
            changeent  = np.where(np.abs(diffvals) > 0)[0]
            if len(changeent) > 0:
                changeval  = diffvals[changeent]
                changedate = [datetime.datetime.strptime(str(dd).split(' ')[0], "%Y-%m-%d") for dd in values2plot.index[changeent+1]]

                lineymin = ymax - dy*0.60
                lineymax = ymax - dy*0.48
                textymin = ymax - dy*0.45
                for cc, cval in enumerate(changeval):
                    if cval < 0:
                        changecol = 'red'
                    else:
                        changecol = 'green'

                    plt.plot([changedate[cc], changedate[cc]],[lineymin, lineymax], '-', color=changecol, lw=lthick, zorder=5,alpha=0.5)
                    cvalstring = changedate[cc].strftime("%d-%m-%Y")+' ('+str(cval)+' senge)'
                    plt.text(changedate[cc], textymin, cvalstring, zorder=5,
                             fontsize=Fsize-2, rotation=90, color=changecol, horizontalalignment='center', verticalalignment='bottom')
        else:
            legscale = 1.5
            ncol     = 4

    # --------- LABELS ---------
    plt.xlabel('Dato')
    plt.ylabel('Antal disponible senge')

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / legscale}, ncol=ncol, numpoints=1,
                     bbox_to_anchor=(0.5, 1.32), )  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------
    plotdir   = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Belægningsprocenter/'
    outputfig = plotdir+plotname
    plt.savefig(outputfig)
    plt.clf()
    plt.close('all')
    if verbose: print(' - Saved figure to \n  '+outputfig)

# -----------------------------------------------------------------------------------------------------------------------
def plot_beds_insets(verbose=True):
    """
    Function to wrap around plotting to generate

    import occupancy_and_beds as oab
    oab.plot_beds_insets()

    """
    dataframe = oab.generate_datastructure(filename='Belægningshistorik_alle_afdelinger_220208.xlsx', verbose=verbose)
    oab.plot_beds(dataframe, SORsections='all', plotname = 'disponiblesenge_alle.pdf')
    # - - - - - - - - - - - - - -Order from Excel sheet - - - - - - - - - - - - - - - - - -
    oab.plot_beds(dataframe, SORsections=['SJ SLAAKI1, AKUT AFD.STUEN, SENGEAFS., SLA'],         plotname = 'disponiblesenge_SLAAKI1.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAAKI2, AKUT AFD. 1.SAL, SENGEAFS., SLA'],        plotname = 'disponiblesenge_SLAAKI2.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLALUI, MED. LUNGE SENGEAFS., SLA'],               plotname = 'disponiblesenge_SLALUI.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAMGIS, MED. GASTRO. SENGEAFS., SLA'],            plotname = 'disponiblesenge_SLAMGIS.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAKAI, KARDIOLOGISK SENGEAFS., SLA'],             plotname = 'disponiblesenge_SLAKAI.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAENIMS, HORMON-MULTISGD. SENGEAFS., SLA'],       plotname = 'disponiblesenge_SLAENIMS.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAMSID, MULTISYGDOM DAGAFSNIT, SLA'],             plotname = 'disponiblesenge_SLAMSID.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLANEI, NEUROLOGISK SENGEAFS., SLA'],              plotname = 'disponiblesenge_SLANEI.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAGEIG1, GERIATRISK SENGEAFS. G1, SLA'],          plotname = 'disponiblesenge_SLAGEIG1.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAGEIG2, GERIATRISK SENGEAFS. G2, SLA'],          plotname = 'disponiblesenge_SLAGEIG2.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAPÆI, PÆD. SENGEAFS., SLA'],                     plotname = 'disponiblesenge_SLAPÆI.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAPÆIN, PÆD. NEO SENGEAFS., SLA'],                plotname = 'disponiblesenge_SLAPÆIN.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAOKIOT, O-KIR. SENGEAFS., TRAUME, SLA'],         plotname = 'disponiblesenge_SLAOKIOT.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAGYI, GYN. SENGEAFSNIT, SLA'],                   plotname = 'disponiblesenge_SLAGYI.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAOBI, OBST. MOR-BARN SENGEAFS., SLA'],           plotname = 'disponiblesenge_SLAOBI.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAOBISV, OBST. GRAVIDITETSAFSNIT, SLA'],          plotname = 'disponiblesenge_SLAOBISV.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAOBIFØMO, OBST. FØDEMODT. SENGEAFS., SLA'],      plotname = 'disponiblesenge_SLAOBIFØMO.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAKGI, KIRURGISK SENGEAFSNIT, SLA'],              plotname = 'disponiblesenge_SLAKGI.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAANITM, INTERMEDIÆRT SENGEAFS., SLA'],           plotname = 'disponiblesenge_SLAANITM.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAANINT, INTENSIV SENGEAFS., SLA'],               plotname = 'disponiblesenge_SLAANINT.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'],               plotname = 'disponiblesenge_NAELUIN.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ NAENEICNN, CENTER NEUROREHAB. SENGEAFS., NAE'],    plotname = 'disponiblesenge_NAENEICNN.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ NAEOKI8, ORTOPÆDKIR. SENGEAFSNIT 8, NAE'],         plotname = 'disponiblesenge_NAEOKI8.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    oab.plot_beds(dataframe, SORsections=['SJ NAENEICN, CENTER NEUROREHAB.SENGEAFS.NAE'],        plotname ='disponiblesenge_NAENEICN.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ NAENEICN, CENTER NEUROREHAB.SENGEAFS.NAE', 'SJ NAENEICNN, CENTER NEUROREHAB. SENGEAFS., NAE'],
                  plotname='disponiblesenge_NAENEICNogNAENEICNN.pdf')

    oab.plot_beds(dataframe, SORsections=['SJ SLAMGIS, MED. GASTRO. SENGEAFS., SLA', 'SJ SLAMGI, MED. GASTRO. SENGEAFS., SLA'],            plotname = 'disponiblesenge_SLAMGISogSLAMGI.pdf')

    oab.plot_beds(dataframe, SORsections=['SJ SLAKAIHA, KARDIOLOGISK SENGEAFS., HA, SLA'],             plotname = 'disponiblesenge_SLAKAIHA.pdf')
    oab.plot_beds(dataframe, SORsections=['SJ SLAKAIHA, KARDIOLOGISK SENGEAFS., HA, SLA', 'SJ SLAKAI, KARDIOLOGISK SENGEAFS., SLA'],
                  plotname='disponiblesenge_SLAKAIogSLAKAIHA.pdf')
    #empty: oab.plot_beds(dataframe, SORsections=['SJ SLAKAIHA2, KARDIOLOGISK SENGEAFS., HA2, SLA'],             plotname = 'disponiblesenge_SLAKAIHA2.pdf')

    oab.plot_beds(dataframe, SORsections=['SJ SLAENI, ENDOKRINOL. SENGEAFS., SLA'], plotname='disponiblesenge_SLAENI.pdf')

# -----------------------------------------------------------------------------------------------------------------------
def plot_occupancy(dataframe,datemin='02-01-2018',datemax='01-01-2022', plotname='belægning.pdf',SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'] ,verbose=True):
    """

    Parameters
    ----------
    datemin
    datemax
    verbose

    Example of use
    -------
    import occupancy_and_beds as oab
    dataframe = oab.generate_datastructure(filename='Belægningshistorik_alle_afdelinger_220208.xlsx', verbose=verbose)
    oab.plot_occupancy(dataframe,SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'])

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datemin = datetime.datetime.strptime(datemin, "%d-%m-%Y")
    datemax = datetime.datetime.strptime(datemax, "%d-%m-%Y")

    df_time = dataframe.set_index('Belægning tidspunkt')
    occdic  = {}

    if SORsections == 'all':
        SORsectionlist = np.unique(dataframe['Ophold afsnit navn'].values.tolist())
    else:
        SORsectionlist = SORsections

    if verbose: print(' - Defined date range and starting to fill array with values')
    for ss, SORsection in enumerate(SORsectionlist):
        if verbose: print('\n   Looking at data for '+SORsection+'  ('+str(ss+1)+'/'+str(len(SORsectionlist))+')')
        badent   = np.where(dataframe['Ophold afsnit navn'] != SORsection)[0]
        df_time  = dataframe.drop(badent).set_index('Belægning tidspunkt')

        datemin_data = np.min(df_time.index)
        datemax_data = np.max(df_time.index)
        daterange = pd.date_range(start=np.max(np.asarray([datemin_data,datemin])), end=np.min(np.asarray([datemax_data,datemax])))
        Ndates = len(daterange)
        occarr   = np.zeros((24, Ndates))*np.nan

        for hh, hour2use in enumerate(np.arange(0, 24, 1)):
            for dd, entdate in enumerate(daterange):
                #if verbose: print('   datetime '+str(entdate)+' and hour '+str(hour2use))
                goodent = np.where((df_time.index.hour == hour2use) &
                                   (df_time.index.date == entdate))[0]
                Ngood = len(goodent)
                if Ngood == 1:
                    occarr[hh,dd] = df_time['Belægning i %'][goodent]
                else:
                    if verbose: print('  Found '+str(Ngood)+' for '+SORsection+' '+str(entdate)+' hour '+str(hour2use))

        occdic[SORsection] = occarr

    if verbose: print(' - Done generating data array; moving on to plotting')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(15, 6))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.05, right=0.97, bottom=0.18, top=0.80)
    Fsize = 10
    lthick = 2
    marksize = 4

    plt.clf()
    plt.ioff()

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    plt.xticks(rotation=45)

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
    ymax = 150.0
    dy   = ymax - ymin
    #plt.ylim([ymin,ymax])

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
    for ss, SORsection in enumerate(SORsectionlist):
        if verbose: print('   Painting curves for '+SORsection+'  ('+str(ss+1)+'/'+str(len(SORsectionlist))+')')
        if len(SORsectionlist) == 1:
             pointcolor = cmap(colnorm(45))
        elif len(SORsectionlist) == 2:
             pointcolor = [cmap(colnorm(30)), cmap(colnorm(80))][ss]
        else:
            pointcolor = cmap(colnorm((cmax-cmin)/len(SORsectionlist)*(ss+1)))

        occarr_hoursort = np.sort(occdic[SORsection], axis=0)
        medval          = np.nanmedian(occarr_hoursort, axis=0)
        CI68low         = occarr_hoursort[3, :]
        CI68high        = occarr_hoursort[19, :]
        CI95low         = occarr_hoursort[1, :]
        CI95high        = occarr_hoursort[22, :]

        duration = (np.max(daterange)-np.min(daterange)).days
        if duration <= 66:
            plt.fill_between(daterange, CI68high, CI68low, color=pointcolor, alpha=1.0, step='mid', zorder=15,
                             label='68% konfidensinterval')
            plt.fill_between(daterange, CI95high, CI95low, color=pointcolor, alpha=0.6, step='mid', zorder=10,
                             label='95% konfidensinterval')
            mediancol = 'black'
        else:
            mediancol = pointcolor

        plt.step(daterange, medval, where='mid',color=mediancol, alpha=1.0, zorder=20,
                 label=SORsection+'\nMedian af belægning over døgnets 24 timer', linestyle='-', lw=lthick)



        if SORsections != 'all':
            diffval_high = occdic[SORsection][23, :] - CI95high
            Nabove = len(np.where(diffval_high[~np.isnan(diffval_high)] > 0)[0])
            fracabove = Nabove / len(CI68high) * 100.

            diffval_low = CI95low - occdic[SORsection][23, :]
            Nbelow = len(np.where(diffval_low[~np.isnan(diffval_low)] > 0)[0])
            fracbelow = Nbelow/len(CI68low) * 100.

            plt.step(daterange, occdic[SORsection][23, :], where='mid', color='red', alpha=0.8, zorder=25,
                     label='Belægning kl. 23'+' ('+str('%.1f' % fracbelow)+'% < 95%KI; '+str('%.1f' % fracabove)+'% > 95%KI)', linestyle='-', lw=lthick)
            legscale   = 1.0
            ncol       = 4
        else:
            legscale = 1.5
            ncol     = 4

    # --------- LABELS ---------
    plt.xlabel('Dato')
    plt.ylabel('Belægningsprocent')

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / legscale}, ncol=ncol, numpoints=1,
                     bbox_to_anchor=(0.5, 1.32), )  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------
    plotdir   = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Belægningsprocenter/'
    outputfig = plotdir+plotname
    plt.savefig(outputfig)
    plt.clf()
    plt.close('all')
    if verbose: print(' - Saved figure to \n  '+outputfig)

# -----------------------------------------------------------------------------------------------------------------------
def plot_occupancy_insets(verbose=True):
    """
    Function to wrap around plotting to generate figures of bed occupancy percentages

    import occupancy_and_beds as oab
    oab.plot_occupancy_insets()

    """
    dataframe = oab.generate_datastructure(filename='Belægningshistorik_alle_afdelinger_220208.xlsx', verbose=verbose)
    # - - - - - - - - - - - - - -Order from Excel sheet - - - - - - - - - - - - - - - - - -
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAAKI1, AKUT AFD.STUEN, SENGEAFS., SLA'],         plotname = 'disponiblesenge_SLAAKI1.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAAKI2, AKUT AFD. 1.SAL, SENGEAFS., SLA'],        plotname = 'disponiblesenge_SLAAKI2.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLALUI, MED. LUNGE SENGEAFS., SLA'],               plotname = 'disponiblesenge_SLALUI.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAMGIS, MED. GASTRO. SENGEAFS., SLA'],            plotname = 'disponiblesenge_SLAMGIS.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAKAI, KARDIOLOGISK SENGEAFS., SLA'],             plotname = 'disponiblesenge_SLAKAI.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAENIMS, HORMON-MULTISGD. SENGEAFS., SLA'],       plotname = 'disponiblesenge_SLAENIMS.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAMSID, MULTISYGDOM DAGAFSNIT, SLA'],             plotname = 'disponiblesenge_SLAMSID.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLANEI, NEUROLOGISK SENGEAFS., SLA'],              plotname = 'disponiblesenge_SLANEI.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAGEIG1, GERIATRISK SENGEAFS. G1, SLA'],          plotname = 'disponiblesenge_SLAGEIG1.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAGEIG2, GERIATRISK SENGEAFS. G2, SLA'],          plotname = 'disponiblesenge_SLAGEIG2.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAPÆI, PÆD. SENGEAFS., SLA'],                     plotname = 'disponiblesenge_SLAPÆI.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAPÆIN, PÆD. NEO SENGEAFS., SLA'],                plotname = 'disponiblesenge_SLAPÆIN.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAOKIOT, O-KIR. SENGEAFS., TRAUME, SLA'],         plotname = 'disponiblesenge_SLAOKIOT.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAGYI, GYN. SENGEAFSNIT, SLA'],                   plotname = 'disponiblesenge_SLAGYI.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAOBI, OBST. MOR-BARN SENGEAFS., SLA'],           plotname = 'disponiblesenge_SLAOBI.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAOBISV, OBST. GRAVIDITETSAFSNIT, SLA'],          plotname = 'disponiblesenge_SLAOBISV.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAOBIFØMO, OBST. FØDEMODT. SENGEAFS., SLA'],      plotname = 'disponiblesenge_SLAOBIFØMO.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAKGI, KIRURGISK SENGEAFSNIT, SLA'],              plotname = 'disponiblesenge_SLAKGI.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAANITM, INTERMEDIÆRT SENGEAFS., SLA'],           plotname = 'disponiblesenge_SLAANITM.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ SLAANINT, INTENSIV SENGEAFS., SLA'],               plotname = 'disponiblesenge_SLAANINT.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'],               plotname = 'disponiblesenge_NAELUIN.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ NAENEICNN, CENTER NEUROREHAB. SENGEAFS., NAE'],    plotname = 'disponiblesenge_NAENEICNN.pdf')
    oab.plot_occupancy(dataframe, SORsections=['SJ NAEOKI8, ORTOPÆDKIR. SENGEAFSNIT 8, NAE'],         plotname = 'disponiblesenge_NAEOKI8.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#=======================================================================================================================

