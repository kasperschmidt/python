#from importlib import reload
import pdb
import sys
import os

import scipy.stats as st
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
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    #plt.gca().xaxis.set_major_locator(mdates.DateLocator(interval=1))
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
def plot_beds_wrapper(verbose=True):
    """
    Function to wrap around plotting to generate

    import occupancy_and_beds as oab
    oab.plot_beds_wrapper()

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

    oab.plot_occupancy(dataframe, datemin='01-12-2021', datemax='01-01-2022', SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'])

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    sys.exit('After building plot_occupancy_CIs it turns out that the manual grouping in this code is faulty; status 220309')
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
        'Belægning i %'
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

        plt.step(daterange, medval, where='mid', color=mediancol, alpha=1.0, zorder=20,
                 label=SORsection+'\nMedian af belægning over døgnets 24 timer', linestyle='-', lw=lthick)


        if SORsections != 'all':
            diffval_high = occdic[SORsection][23, :] - CI95high
            Nabove = len(np.where(diffval_high[~np.isnan(diffval_high)] > 0)[0])
            fracabove = Nabove / len(CI68high) * 100.

            diffval_low = CI95low - occdic[SORsection][23, :]
            Nbelow = len(np.where(diffval_low[~np.isnan(diffval_low)] > 0)[0])
            fracbelow = Nbelow/len(CI68low) * 100.

            plt.step(daterange, occdic[SORsection][23, :], where='mid', color='red', alpha=0.6, zorder=25,
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
def get_scipyCI(dataframe, col_group, col_measure, group_dt='year', location='mean', CI=95.0 ,verbose=True):
    """
    Confidence interval based on Student's T distribution
    For low-number normally distributed data.
    For high-number normally distributed data one can use stats.norm.interval(0.68, loc=mu, scale=sigma)


    import occupancy_and_beds as oab

    #### Occupancy at 23hours
    badentcomb    = np.where((dataframe['Ophold afsnit navn'] != 'SJ SLALUI, MED. LUNGE SENGEAFS., SLA') | (dataframe.set_index('Belægning tidspunkt').index.hour != 23))[0]
    df_time       = dataframe.drop(badentcomb)
    loc, lower, upper = oab.get_scipyCI(df_time, 'Belægning tidspunkt', 'Belægning i %', location='mean', CI=95.0, verbose=True)

    #### CI for a specific ward
    badentcomb    = np.where((dataframe['Ophold afsnit navn'] != 'SJ SLALUI, MED. LUNGE SENGEAFS., SLA'))[0]
    df_time       = dataframe.drop(badentcomb)
    loc, lower, upper = oab.get_scipyCI(df_time, 'Belægning tidspunkt', 'Belægning i %', location='mean', CI=95.0, verbose=True)

    """
    # options for grouping: dayofweek, date, month, year, quarter, week, day_og_week, etc:
    #
    if group_dt == 'year':
        # counts contains a pd.Series with sample size for each category
        counts = dataframe.groupby(dataframe[col_group].dt.year, as_index=False)[col_measure].count()
        # defining locations of central points for grouping
        loc = dataframe.groupby(dataframe[col_group].dt.year)[col_measure].agg(location.lower())
        # standard deviation and standard error
        std = dataframe.groupby(dataframe[col_group].dt.year)[col_measure].agg(np.std)
    elif group_dt == 'date':
        # counts contains a pd.Series with sample size for each category
        counts = dataframe.groupby(dataframe[col_group].dt.date, as_index=False)[col_measure].count()
        # defining locations of central points for grouping
        loc = dataframe.groupby(dataframe[col_group].dt.date)[col_measure].agg(location.lower())
        # standard deviation
        std = dataframe.groupby(dataframe[col_group].dt.date)[col_measure].agg(np.std)
    else:
        # counts contains a pd.Series with sample size for each category
        counts = dataframe.groupby(dataframe[col_group].dt.year, as_index=False)[col_measure].count()

        # cat has names of the categories, like 'category 1', 'category 2'
        #cat = list(dataframe.groupby(col_group, as_index=False)[col_measure].count()[col_group])
        #dataframe.groupby(dataframe[col_group].dt.year, as_index=False)[col_measure].count()
        #[np.unique(dataframe[col_group].dt.year)]

        # the average value of col2 across the categories
        #loc = dataframe.groupby(col_group)[col_measure].agg(location.lower())
        loc = dataframe.groupby(dataframe[col_group].dt.year)[col_measure].agg(location.lower())
        # standard deviation and standard error
        std = dataframe.groupby(dataframe[col_group].dt.year)[col_measure].agg(np.std)


    # standard error for grouping values
    #se  = std / np.sqrt(counts)
    se = std.values / np.sqrt(np.transpose(counts.values))

    lower, upper = st.t.interval(alpha=CI/100., df=np.transpose(counts.values) - 1, loc=loc, scale=se)

    return loc, lower, upper

# -----------------------------------------------------------------------------------------------------------------------
def plot_occupancy_CIs(dataframe, datemin='02-01-2018', datemax='01-01-2022', plotname='belægningCI.pdf',
                       SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'], verbose=True):
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
    oab.plot_occupancy_CIs(dataframe, datemin='01-12-2021', datemax='01-01-2022', SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'])
    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datemin = datetime.datetime.strptime(datemin, "%d-%m-%Y")
    datemax = datetime.datetime.strptime(datemax, "%d-%m-%Y")

    occdic = {}

    if SORsections == 'all':
        SORsectionlist = np.unique(dataframe['Ophold afsnit navn'].values.tolist())
    else:
        SORsectionlist = SORsections

    if verbose: print(' - Defined date range and starting to fill array with values')
    for ss, SORsection in enumerate(SORsectionlist):
        if verbose: print(
            '   Looking at data for ' + SORsection + '  (' + str(ss + 1) + '/' + str(len(SORsectionlist)) + ')')
        #badent  = np.where(dataframe['Ophold afsnit navn'] != SORsection)[0]
        #df_time = dataframe.drop(badent).set_index('Belægning tidspunkt')

        #### Occupancy at 23 o'clock
        badentcomb23    = np.where((dataframe['Ophold afsnit navn'] != SORsection) | (dataframe.set_index('Belægning tidspunkt').index.hour != 23))[0]
        df_drop23       = dataframe.drop(badentcomb23)
        datemin_data23  = np.min(df_drop23['Belægning tidspunkt'])
        datemax_data23  = np.max(df_drop23['Belægning tidspunkt'])
        df_range23      = df_drop23[(df_drop23.set_index('Belægning tidspunkt').index >= np.max(np.asarray([datemin_data23, datemin]))) &
                                (df_drop23.set_index('Belægning tidspunkt').index <= np.min(np.asarray([datemax_data23, datemax])))]

        #### CI for a specific ward on daily basis
        badentcomb    = np.where((dataframe['Ophold afsnit navn'] != SORsection))[0]
        df_drop       = dataframe.drop(badentcomb)
        datemin_data  = np.min(df_drop['Belægning tidspunkt'])
        datemax_data  = np.max(df_drop['Belægning tidspunkt'])

        df_range      = df_drop[(df_drop.set_index('Belægning tidspunkt').index >= np.max(np.asarray([datemin_data, datemin]))) &
                                (df_drop.set_index('Belægning tidspunkt').index <= np.min(np.asarray([datemax_data, datemax])))]

        dfgroups   = df_range.groupby(df_range['Belægning tidspunkt'].dt.date, sort=True)['Belægning i %']

        dfgroups_sort = df_range.sort_values(['Belægning i %'], ascending=True).groupby(df_range['Belægning tidspunkt'].dt.date)['Belægning i %']
        #dfgroups_sort.get_group(datetime.date(2021, 12, 25))

        #counts        = df_range.groupby(df_range['Belægning tidspunkt'].dt.date, as_index=False)['Belægning i %'].count()
        #dfgroups_core = df_range.groupby(df_range['Belægning tidspunkt'].dt.date, sort=True)['Belægning i %'].nth([0, -1])

        CI95index = [1, -2]
        CI68index = [3, -5]
        CI95low  = dfgroups_sort.nth([CI95index[0]])
        CI95high = dfgroups_sort.nth([CI95index[1]])
        CI68low  = dfgroups_sort.nth([CI68index[0]])
        CI68high = dfgroups_sort.nth([CI68index[1]])

        group_median = df_range.groupby(df_range['Belægning tidspunkt'].dt.date)['Belægning i %'].agg('median')
        groupnames   = [grkey for grkey in dict(list(dfgroups)).keys()]

        occdic[SORsection] = [df_range23, group_median, CI68low, CI68high, CI95low, CI95high]

    checkcalculations = False
    if checkcalculations:
        datevalues = df_range['Belægning i %'].values[np.where(df_range['Belægning tidspunkt'].dt.date.values == datetime.date(2021, 12, 25))]
        dv_mean    = np.mean(datevalues)
        dv_median  = np.median(datevalues)
        dv_std     = np.std(datevalues)
        print(dv_mean, dv_median, dv_std)

        counts     = df_range.groupby(df_range['Belægning tidspunkt'].dt.date, as_index=False)['Belægning i %'].count()
        loc_mean   = df_range.groupby(df_range['Belægning tidspunkt'].dt.date)['Belægning i %'].agg('mean')
        loc_median = df_range.groupby(df_range['Belægning tidspunkt'].dt.date)['Belægning i %'].agg('median')
        std        = df_range.groupby(df_range['Belægning tidspunkt'].dt.date)['Belægning i %'].agg(np.std)
        se = std.values / np.sqrt(np.transpose(counts.values))
        #st.sem() # standard error of men
        print(loc_mean, loc_median, std)

        lower, upper = st.t.interval(alpha=95 / 100., df=np.transpose(counts.values) - 1, loc=loc_mean, scale=se)

        print(datevalues)
        print(lower[0][0],upper[0][0])

        dateent  = np.where((dataframe['Belægning tidspunkt'].dt.date.values == datetime.date(2021, 12, 3)) & (dataframe['Ophold afsnit navn'] == SORsection))
        datedata = np.asarray(list(dataframe['Belægning i %'].values[dateent]))
        times    = np.asarray(list(dataframe['Belægning tidspunkt'].values[dateent]))
        #np.mean(datedata[~np.isnan(datedata)])
        val23 = datedata[-1]
        np.median(datedata[~np.isnan(datedata)])

    if verbose: print(' - Done generating data dictionary; moving on to plotting\n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - Initiating ' + plotname)

    fig = plt.figure(figsize=(15, 6))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.05, right=0.97, bottom=0.18, top=0.80)
    Fsize = 10
    lthick = 2
    marksize = 4

    plt.clf()
    plt.ioff()

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
    plt.xticks(rotation=45)

    plt.rc('text', usetex=False)
    plt.rc('font', family='serif', size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    # plt.title(inforstr[:-2],fontsize=Fsize)

    xerr = None
    yerr = None

    plt.grid(linestyle=':', linewidth=lthick / 2.)
    # --------- RANGES ---------
    ymin = 0
    ymax = 150.0
    dy = ymax - ymin
    # plt.ylim([ymin,ymax])

    # plt.xscale('log')
    # plt.yscale('log')

    # --------- COLORMAP ---------
    cmap = plt.cm.get_cmap('viridis')
    cmin = 0
    cmax = 100

    colnorm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, num=50)
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(cmaparr)
    # cb = plt.colorbar(m)
    # cb.set_label('belægningsprocent')

    # --------- POINT AND CURVES ---------
    for ss, SORsection in enumerate(SORsectionlist):
        if verbose: print(
            '   Painting curves for ' + SORsection + '  (' + str(ss + 1) + '/' + str(len(SORsectionlist)) + ')')
        if len(SORsectionlist) == 1:
            pointcolor = cmap(colnorm(45))
        elif len(SORsectionlist) == 2:
            pointcolor = [cmap(colnorm(30)), cmap(colnorm(80))][ss]
        else:
            pointcolor = cmap(colnorm((cmax - cmin) / len(SORsectionlist) * (ss + 1)))

        df_range23, group_median, CI68low, CI68high, CI95low, CI95high = occdic[SORsection]

        duration = (np.max(group_median.index) - np.min(group_median.index)).days
        if duration <= 100:
            plt.gcf().autofmt_xdate()  # set font and rotation for date tick labels
            plt.fill_between(list(CI68high.index), list(CI68high.values), list(CI68low.values), color=pointcolor, alpha=1.0, step='mid', zorder=15,
                             label='68% konfidensinterval')
            plt.fill_between(list(CI95high.index), list(CI95high.values), list(CI95low.values), color=pointcolor, alpha=0.6, step='mid', zorder=10,
                             label='95% konfidensinterval')

            centercol = 'black'
        else:
            plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))
            centercol = pointcolor

        plt.step(group_median.index, group_median.values, where='mid', color=centercol, alpha=1.0, zorder=20,
                 label=SORsection + '\nMedian af belægningen over døgnets 24 timer', linestyle='-', lw=lthick)

        if SORsections != 'all':
            diffval_high = np.asarray(list(df_range23['Belægning i %'].values - CI95high.values))
            Nabove = len(np.where(diffval_high[~np.isnan(diffval_high)] > 0)[0])
            fracabove = Nabove / len(CI68high) * 100.

            diffval_low = np.asarray(list(CI95low.values - df_range23['Belægning i %'].values))
            Nbelow = len(np.where(diffval_low[~np.isnan(diffval_low)] > 0)[0])
            fracbelow = Nbelow / len(CI68low) * 100.

            plt.step(df_range23['Belægning kalenderdato'].values, df_range23['Belægning i %'].values, where='mid', color='red', alpha=0.6, zorder=25,
                     label='Belægning kl. 23' + ' (' + str('%.1f' % fracbelow) + '% < 95%KI; ' + str(
                         '%.1f' % fracabove) + '% > 95%KI)', linestyle='-', lw=lthick)
            legscale = 1.0
            ncol = 4
        else:
            legscale = 1.5
            ncol = 4

    # --------- LABELS ---------
    plt.xlabel('Dato')
    plt.ylabel('Belægningsprocent')

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / legscale}, ncol=ncol, numpoints=1,
                     bbox_to_anchor=(0.5, 1.32), )  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------
    plotdir = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Belægningsprocenter/'
    outputfig = plotdir + plotname
    plt.savefig(outputfig)
    plt.clf()
    plt.close('all')
    if verbose: print(' - Saved figure to \n  ' + outputfig)

# -----------------------------------------------------------------------------------------------------------------------
def plot_occupancy_wrapper(verbose=True):
    """
    Function to wrap around plotting to generate figures of bed occupancy percentages

    import occupancy_and_beds as oab
    oab.plot_occupancy_wrapper()

    """
    dataframe = oab.generate_datastructure(filename='Belægningshistorik_alle_afdelinger_220208.xlsx', verbose=verbose)
    # - - - - - - - - - - - - - -Order from Excel sheet - - - - - - - - - - - - - - - - - -
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAAKI1, AKUT AFD.STUEN, SENGEAFS., SLA'],         plotname = 'belægningCIs_SLAAKI1.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAAKI2, AKUT AFD. 1.SAL, SENGEAFS., SLA'],        plotname = 'belægningCIs_SLAAKI2.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLALUI, MED. LUNGE SENGEAFS., SLA'],               plotname = 'belægningCIs_SLALUI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAMGIS, MED. GASTRO. SENGEAFS., SLA'],            plotname = 'belægningCIs_SLAMGIS.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAKAI, KARDIOLOGISK SENGEAFS., SLA'],             plotname = 'belægningCIs_SLAKAI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAENIMS, HORMON-MULTISGD. SENGEAFS., SLA'],       plotname = 'belægningCIs_SLAENIMS.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAMSID, MULTISYGDOM DAGAFSNIT, SLA'],             plotname = 'belægningCIs_SLAMSID.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLANEI, NEUROLOGISK SENGEAFS., SLA'],              plotname = 'belægningCIs_SLANEI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAGEIG1, GERIATRISK SENGEAFS. G1, SLA'],          plotname = 'belægningCIs_SLAGEIG1.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAGEIG2, GERIATRISK SENGEAFS. G2, SLA'],          plotname = 'belægningCIs_SLAGEIG2.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAPÆI, PÆD. SENGEAFS., SLA'],                     plotname = 'belægningCIs_SLAPÆI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAPÆIN, PÆD. NEO SENGEAFS., SLA'],                plotname = 'belægningCIs_SLAPÆIN.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAOKIOT, O-KIR. SENGEAFS., TRAUME, SLA'],         plotname = 'belægningCIs_SLAOKIOT.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAGYI, GYN. SENGEAFSNIT, SLA'],                   plotname = 'belægningCIs_SLAGYI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAOBI, OBST. MOR-BARN SENGEAFS., SLA'],           plotname = 'belægningCIs_SLAOBI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAOBISV, OBST. GRAVIDITETSAFSNIT, SLA'],          plotname = 'belægningCIs_SLAOBISV.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAOBIFØMO, OBST. FØDEMODT. SENGEAFS., SLA'],      plotname = 'belægningCIs_SLAOBIFØMO.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAKGI, KIRURGISK SENGEAFSNIT, SLA'],              plotname = 'belægningCIs_SLAKGI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAANITM, INTERMEDIÆRT SENGEAFS., SLA'],           plotname = 'belægningCIs_SLAANITM.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAANINT, INTENSIV SENGEAFS., SLA'],               plotname = 'belægningCIs_SLAANINT.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'],               plotname = 'belægningCIs_NAELUIN.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ NAENEICNN, CENTER NEUROREHAB. SENGEAFS., NAE'],    plotname = 'belægningCIs_NAENEICNN.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ NAEOKI8, ORTOPÆDKIR. SENGEAFSNIT 8, NAE'],         plotname = 'belægningCIs_NAEOKI8.pdf')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    datemin='01-11-2021'
    datemax='01-02-2022'

    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAAKI1, AKUT AFD.STUEN, SENGEAFS., SLA'],      datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAAKI1.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAAKI2, AKUT AFD. 1.SAL, SENGEAFS., SLA'],     datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAAKI2.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLALUI, MED. LUNGE SENGEAFS., SLA'],            datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLALUI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAMGIS, MED. GASTRO. SENGEAFS., SLA'],         datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAMGIS.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAKAI, KARDIOLOGISK SENGEAFS., SLA'],          datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAKAI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAENIMS, HORMON-MULTISGD. SENGEAFS., SLA'],    datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAENIMS.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAMSID, MULTISYGDOM DAGAFSNIT, SLA'],          datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAMSID.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLANEI, NEUROLOGISK SENGEAFS., SLA'],           datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLANEI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAGEIG1, GERIATRISK SENGEAFS. G1, SLA'],       datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAGEIG1.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAGEIG2, GERIATRISK SENGEAFS. G2, SLA'],       datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAGEIG2.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAPÆI, PÆD. SENGEAFS., SLA'],                  datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAPÆI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAPÆIN, PÆD. NEO SENGEAFS., SLA'],             datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAPÆIN.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAOKIOT, O-KIR. SENGEAFS., TRAUME, SLA'],      datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAOKIOT.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAGYI, GYN. SENGEAFSNIT, SLA'],                datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAGYI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAOBI, OBST. MOR-BARN SENGEAFS., SLA'],        datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAOBI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAOBISV, OBST. GRAVIDITETSAFSNIT, SLA'],       datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAOBISV.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAOBIFØMO, OBST. FØDEMODT. SENGEAFS., SLA'],   datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAOBIFØMO.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAKGI, KIRURGISK SENGEAFSNIT, SLA'],           datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAKGI.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAANITM, INTERMEDIÆRT SENGEAFS., SLA'],        datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAANITM.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ SLAANINT, INTENSIV SENGEAFS., SLA'],            datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_SLAANINT.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ NAELUIN, LUNGEMED. SENGEAFS., NAE'],            datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_NAELUIN.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ NAENEICNN, CENTER NEUROREHAB. SENGEAFS., NAE'], datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_NAENEICNN.pdf')
    oab.plot_occupancy_CIs(dataframe, SORsections=['SJ NAEOKI8, ORTOPÆDKIR. SENGEAFSNIT 8, NAE'],      datemin=datemin, datemax=datemax,   plotname = 'belægningCIs_Nov21Jan22_NAEOKI8.pdf')
#=======================================================================================================================

