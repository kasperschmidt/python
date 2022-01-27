#from importlib import reload
import pdb
import sys
import os

import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
import seaborn
import datetime

#sys.path.append('C:/Users/kaschm/GitHub/python/pythonNSR/')
import occupancy_and_beds as oab
import lungemed_contactdays_updates as lcu

# -----------------------------------------------------------------------------------------------------------------------
def load_dataframes_from_excel(verbose=True):
    """

    df_baseline, df_updates = lcu.load_dataframes_from_excel()

    """
    path_baseline = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/plots/'
    file_baseline = 'lungemedLPR3_SQLbeddays.xlsx'
    #file_baseline = 'lungemedLPR3dataframe.xlsx'
    if verbose: print(' - loading dataframe from '+file_baseline)
    df_baseline   = pd.read_excel(path_baseline + file_baseline, sheet_name='data output')

    path_updates = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/månedligetræk fra SP/'

    updatenames = 'lungemed_SPupdate*.xlsx'
    filelist_updates  = glob.glob(path_updates + updatenames)
    updates_framelist = []
    if verbose: print(' - loading and concatenating '+str(len(updates_framelist))+' dataframes with updates from ' + updatenames)
    for file_update in filelist_updates:
        if verbose: print('    - including '+file_update.split('/')[-1])
        df_update = pd.read_excel(file_update, sheet_name='data output')
        updates_framelist = updates_framelist + [df_update]

    df_updates = pd.concat(updates_framelist, ignore_index=True)

    if verbose: print(' - returning dataframes for baseline and updates')
    return df_baseline, df_updates

# -----------------------------------------------------------------------------------------------------------------------
def evaluate_beddays(datemin_mostrecent='01-01-2022',verbose=True):
    """

    lcu.evaluate_beddays()

    """
    plotdir = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/plots/'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    df_baseline, df_updates = lcu.load_dataframes_from_excel(verbose=verbose)
    df_baseline = df_baseline.sort_values('INDTIDSPUNKT_DRGKONTAKT')
    df_updates  = df_updates.sort_values('Kontakt startdato Dato-tid')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'beddays_stat_updates.pdf'
    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(9, 6))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.1, right=0.97, bottom=0.10, top=0.90)
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

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d-%Y'))
    plt.gca().xaxis.set_major_locator(mdates.YearLocator())

    xerr = None
    yerr = None

    plt.grid(linestyle=':', linewidth=lthick/2.)
    # --------- RANGES ---------
    ymin = -2
    ymax = 30
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
    pointcolor = cmap(colnorm(60))
    # df_baseline_month = df_baseline.set_index('INDTIDSPUNKT_DRGKONTAKT')
    # subsamp = df_baseline_month[df_baseline_month.index.month == 1]
    # plt.errorbar(subsamp['INDTIDSPUNKT_DRGKONTAKT'], subsamp['KONTAKTDAGE'], xerr=xerr, yerr=yerr,
    #              marker='o', lw=0, markersize=marksize, alpha=0.5, color=pointcolor,
    #              markerfacecolor=pointcolor, ecolor=pointcolor,
    #              markeredgecolor=pointcolor, zorder=10.,
    #              label='Kontaktdage Lungemed. NAE samplet på månedsbasis (Data fra SP)')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Plotting the baseline of the data 2019, 2020, 2021 (excluding december 2021) data
    showbaselinealone = True
    if showbaselinealone:
        pointcolor = cmap(colnorm(30))
        Ndaysavg = 30

        datemax = datetime.datetime.strptime('01-12-2021', "%d-%m-%Y")
        showval = np.where(df_baseline['INDTIDSPUNKT_DRGKONTAKT'] < datemax)[0]
        dropval = np.where(df_baseline['INDTIDSPUNKT_DRGKONTAKT'] >= datemax)[0]

        df_baseline['KONTAKTDAGE_movingavg'] = df_baseline['KONTAKTDAGE'].rolling(window=Ndaysavg).median()
        df_baseline['KONTAKTDAGE_movingavg_std'] = df_baseline['KONTAKTDAGE'].rolling(window=Ndaysavg).std()

        plt.errorbar(df_baseline['INDTIDSPUNKT_DRGKONTAKT'].values[showval], df_baseline['KONTAKTDAGE_movingavg'].values[showval], xerr=xerr, yerr=yerr,
                     marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                     markerfacecolor=pointcolor, ecolor=pointcolor,
                     markeredgecolor=pointcolor, zorder=25.,
                     label=str(Ndaysavg)+' dage glidende middel med \nstandard afvigelse (Data: LPR3)')

        ylow  = df_baseline['KONTAKTDAGE_movingavg']-df_baseline['KONTAKTDAGE_movingavg_std']
        yhigh = df_baseline['KONTAKTDAGE_movingavg']+df_baseline['KONTAKTDAGE_movingavg_std']
        plt.fill_between(df_baseline['INDTIDSPUNKT_DRGKONTAKT'].values[showval], ylow.values[showval], yhigh.values[showval], color=pointcolor, alpha=0.4, zorder=15)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Plotting data from SP since start 2022
    pointcolor = cmap(colnorm(70))
    Ndaysavg = 30
    df_baseline_trimmed = df_baseline.drop(df_baseline.index[dropval])
    dates_combined = df_baseline_trimmed['INDTIDSPUNKT_DRGKONTAKT'].append(df_updates['Kontakt startdato Dato-tid'])
    days_combined  = df_baseline_trimmed['KONTAKTDAGE'].append(df_updates['Forskel på kontakt start og slut (antal dage)'])

    days_combined_movingavg     = days_combined.rolling(window=Ndaysavg).median()
    days_combined_movingavg_std = days_combined.rolling(window=Ndaysavg).std()

    datemin = datetime.datetime.strptime('30-11-2021', "%d-%m-%Y")
    showval = np.where(dates_combined > datemin)[0]
    plt.errorbar(dates_combined.values[showval], days_combined_movingavg.values[showval], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label=str(Ndaysavg)+' dage glidende middel med \nstandard afvigelse (Data: SP/WebI)')

    ylow  = days_combined_movingavg-days_combined_movingavg_std
    yhigh = days_combined_movingavg+days_combined_movingavg_std
    plt.fill_between(dates_combined.values[showval], ylow.values[showval], yhigh.values[showval], color=pointcolor, alpha=0.3, zorder=15)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Plotting last months points
    datemin_mostrecent = datetime.datetime.strptime(datemin_mostrecent, "%d-%m-%Y")
    showval    = np.where(df_updates['Kontakt startdato Dato-tid'] > datemin_mostrecent)[0]
    pointcolor = cmap(colnorm(80))
    plt.errorbar(df_updates['Kontakt startdato Dato-tid'].values[showval], df_updates['Forskel på kontakt start og slut (antal dage)'].values[showval],
                 xerr=xerr, yerr=yerr, marker='o', lw=0, markersize=marksize, alpha=0.5, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=10.,
                 label='Kontaktdage Lungemed. \nNAE (Data: SP/WebI)')

    # --------- LABELS ---------
    plt.xlabel('Dato for måling')
    plt.ylabel('Kontaktdage Lungemd. NAE')

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=3, numpoints=1,
                     bbox_to_anchor=(0.5, 1.12), )  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    if verbose: print(' - saved plot to '+plotdir)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    #Locate the indices for the subsamples?? loc[pd.date_range(ts.index.min(), ts.index.max(), freq='M')]

# -----------------------------------------------------------------------------------------------------------------------
def evaluate_beddays_boxplot(datemin_mostrecent='01-01-2022',verbose=True):
    """

    lcu.evaluate_beddays()

    """
    plotdir = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/plots/'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    df_baseline, df_updates = lcu.load_dataframes_from_excel(verbose=verbose)
    df_baseline = df_baseline.sort_values('INDTIDSPUNKT_DRGKONTAKT')
    df_updates  = df_updates.sort_values('Kontakt startdato Dato-tid')

    datemax = datetime.datetime.strptime('01-12-2021', "%d-%m-%Y")
    dropval = np.where(df_baseline['INDTIDSPUNKT_DRGKONTAKT'] >= datemax)[0]

    df_baseline_trimmed = df_baseline.drop(df_baseline.index[dropval])
    dates_combined = df_baseline_trimmed['INDTIDSPUNKT_DRGKONTAKT'].append(df_updates['Kontakt startdato Dato-tid'])
    days_combined  = df_baseline_trimmed['KONTAKTDAGE'].append(df_updates['Forskel på kontakt start og slut (antal dage)'])

    xvals = dates_combined.dt.to_period('M').map(lambda s: s.strftime('%m-%Y'))
    yvals = days_combined
    # xvals = df_baseline['INDTIDSPUNKT_DRGKONTAKT'].dt.to_period('M').map(lambda s: s.strftime('%m-%Y'))
    # yvals = df_baseline['KONTAKTDAGE']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'beddays_stat_updates_boxplot.pdf'
    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(9, 6))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.1, right=0.97, bottom=0.20, top=0.98)
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

    xerr = None
    yerr = None

    plt.grid(linestyle=':', linewidth=lthick/2., zorder=100, alpha=1.0, color='gray')
    # --------- RANGES ---------
    ymin = -2
    ymax = 30
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
    pointcolor = cmap(colnorm(30))

    zorder_boxplot = 10
    flierprops   = dict(marker='o', markerfacecolor=pointcolor, markersize=marksize, markeredgecolor=pointcolor, alpha=0.3)
    boxprops     = dict(linewidth=lthick/2, color=pointcolor, facecolor=pointcolor, zorder=zorder_boxplot)
    whiskerprops = dict(linewidth=lthick/2, color=pointcolor, zorder=zorder_boxplot)
    medianprops  = dict(linewidth=lthick/2, color=pointcolor, zorder=zorder_boxplot)
    capprops     = dict(linewidth=lthick/2, color=pointcolor, zorder=zorder_boxplot)
    ax = seaborn.boxplot(x=xvals, y=yvals, linewidth=lthick, showfliers=True,
                         flierprops=flierprops, boxprops=boxprops, whiskerprops=whiskerprops,
                         medianprops=medianprops, capprops=capprops,
                         color=pointcolor)
    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0.5))

    showrawdata = False
    if showrawdata:
        pointcolor = 'red'
        plt.errorbar(xvals, yvals, xerr=xerr, yerr=yerr,
                     marker='.', lw=0, markersize=marksize, alpha=0.4, color=pointcolor,
                     markerfacecolor=pointcolor, ecolor=pointcolor,
                     markeredgecolor=pointcolor, zorder=25.,
                     label=' rå data')

    # plt.plot([datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"), datetime.datetime.strptime("01-03-2021", "%d-%m-%Y")],
    #          [ymin, ymax], '--', color='green', lw=lthick, zorder=5)
    plt.plot([33.5, 33.5], [ymin, ymax], '--', color='gray', lw=lthick, zorder=5)


    # --------- LABELS ---------
    #fig.autofmt_xdate()
    plt.xticks(rotation=90)
    plt.xlabel('Måned og år for måling')
    plt.ylabel('Kontaktdage Lungemd. NAE')

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=3, numpoints=1,
                     bbox_to_anchor=(0.5, 1.12), )  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    if verbose: print(' - saved plot to '+plotdir)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# -----------------------------------------------------------------------------------------------------------------------
def beddays_distributions(datemin_mostrecent = '01-01-2022', normalization="probability", verbose=True):
    """
    lcu.beddays_distributions(normalization="probability")
    lcu.beddays_distributions(normalization="count")
    """
    plotdir = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/plots/'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    df_baseline, df_updates = lcu.load_dataframes_from_excel()

    datemin     = datetime.datetime.strptime(datemin_mostrecent, "%d-%m-%Y")
    dropval     = np.where(df_updates['Kontakt startdato Dato-tid'] < datemin)[0]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'distributions_days_'+normalization+'.pdf'
    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(6, 6))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.12, right=0.98, bottom=0.10, top=0.98)
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

    # --------- RANGES ---------
    xmin = 0
    xmax = 40
    dy   = xmax - xmin
    plt.xlim([xmin,xmax])

    plt.grid(linestyle=':', linewidth=lthick/2.)

    # --------- COLORMAP ---------
    cmap = plt.cm.get_cmap('viridis')
    cmin = 0
    cmax = 100

    colnorm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, num=50)
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(cmaparr)

    # --------- POINT AND CURVES ---------
    distcolor = cmap(colnorm(30))
    seaborn.histplot(data=df_baseline, x='KONTAKTDAGE', color=distcolor, alpha=0.3,
                     kde=True, binwidth=1, discrete=True, stat=normalization, common_norm=False,
                     label='Baseline (2019-2021)')

    distcolor = cmap(colnorm(70))
    seaborn.histplot(data=df_updates, x='Forskel på kontakt start og slut (antal dage)', color=distcolor, alpha=0.3,
                     kde=True, binwidth=1, discrete=True, stat=normalization, common_norm=False,
                     label='Forløb siden 01-12-2021')

    distcolor = 'black'
    seaborn.histplot(data=df_updates.drop(df_baseline.index[dropval]), x='Forskel på kontakt start og slut (antal dage)', color=distcolor, alpha=0.3,
                     kde=True, binwidth=1, discrete=True, stat=normalization, common_norm=False,
                     label='Forløb siden '+datemin_mostrecent)

    # --------- LABELS ---------
    plt.xlabel('Kontaktdage')
    ylabelstr = 'Antal forløb '
    if normalization == 'probability':
        ylabelstr = ylabelstr+'(brøkdel af samlede antal forløb)'
    plt.ylabel(ylabelstr)

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=1, numpoints=1,
                     bbox_to_anchor=(0.7, 0.95))  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    if verbose: print(' - saved plot to '+plotdir)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# -----------------------------------------------------------------------------------------------------------------------
def diagnoses_distributions(datemin_mostrecent = '01-01-2022', normalization="probability",verbose=True):
    """
    lcu.diagnoses_distributions()
    """
    plotdir = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/plots/'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    df_baseline, df_updates = lcu.load_dataframes_from_excel()
    df_baseline['DIA01'] = pd.Categorical(df_baseline['DIA01'], np.sort(np.unique(np.asarray(df_baseline['DIA01']))))

    datemin     = datetime.datetime.strptime(datemin_mostrecent, "%d-%m-%Y")
    dropval     = np.where(df_updates['Kontakt startdato Dato-tid'] < datemin)[0]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'distributions_diagnoses_'+normalization+'.pdf'
    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(11, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.12, right=0.98, bottom=0.15, top=0.98)
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

    # --------- RANGES ---------
    xmin = 0
    xmax = 40
    dy   = xmax - xmin
    #plt.xlim([xmin,xmax])

    #plt.grid(linestyle=':', linewidth=lthick/2.)

    # --------- COLORMAP ---------
    cmap = plt.cm.get_cmap('viridis')
    cmin = 0
    cmax = 100

    colnorm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, num=50)
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(cmaparr)

    # --------- POINT AND CURVES ---------
    distcolor = cmap(colnorm(30))
    seaborn.histplot(data=df_baseline, x='DIA01', color=distcolor, alpha=0.3,
                     kde=True, binwidth=1, discrete=True, stat=normalization, common_norm=False,
                     label='Baseline (2019-2021)')

    distcolor = cmap(colnorm(70))
    seaborn.histplot(data=df_updates, x='Aktionsdiagnosekode', color=distcolor, alpha=0.3,
                     kde=True, binwidth=1, discrete=True, stat=normalization, common_norm=False,
                     label='Forløb siden 01-12-2021')

    distcolor = 'black'
    seaborn.histplot(data=df_updates.drop(df_baseline.index[dropval]), x='Aktionsdiagnosekode', color=distcolor, alpha=0.3,
                     kde=True, binwidth=1, discrete=True, stat=normalization, common_norm=False,
                     label='Forløb siden '+datemin_mostrecent)

    # --------- LABELS ---------
    plt.xticks(rotation=90, fontsize=1.4)
    plt.xlabel('Kontaktdage')
    ylabelstr = 'Antal forløb med diagnose '
    if normalization == 'probability':
        ylabelstr = ylabelstr+'(brøkdel af samlede antal forløb)'
    plt.ylabel(ylabelstr)

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=1, numpoints=1,
                     bbox_to_anchor=(0.7, 0.95))  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    if verbose: print(' - saved plot to '+plotdir)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



# -----------------------------------------------------------------------------------------------------------------------
def heatmap_diagnosisVSbeddays(verbose=True):
    """
    lcu.beddays_distributions()
    """


#=======================================================================================================================

