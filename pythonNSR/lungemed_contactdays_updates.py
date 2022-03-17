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
from matplotlib.ticker import PercentFormatter
import seaborn
import datetime
import loadMDCgroups
import scipy

#sys.path.append('C:/Users/kaschm/GitHub/python/pythonNSR/')
import occupancy_and_beds as oab
import lungemed_contactdays_updates as lcu
import lungemed_beddays_and_visitations as lbv

# -----------------------------------------------------------------------------------------------------------------------
def runall_for_updating(datemin_mostrecent='01-01-2022', verbose=True):
    """
    simple wrapper to re-generate plots after adding data updates to
    O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/månedligetræk fra SP/

    The datemin_mostrecent keyword indicates the first date of the chunk of data to treat as "last month".
    When defining these data, everything prior to that data is ignored.
    Format is "%d-%m-%Y"

    import lungemed_contactdays_updates as lcu
    lcu.runall_for_updating(datemin_mostrecent='01-02-2022', verbose=True)

    """

    if verbose: print('\n-------------------- Initiating full run with datemine_mostrecent = '+datemin_mostrecent+' --------------------\n')

    lbv.plot_perday_occupancy(measurehours=[23], include_updates=True, loaddatafile='lungemedLPR3dataframe.xlsx', datemin_mostrecent=datemin_mostrecent)
    lcu.evaluate_beddays_boxplot(verbose=verbose)
    lcu.beddays_distributions(datemin_mostrecent=datemin_mostrecent, normalization="probability", verbose=verbose)
    lcu.diagnoses_distributions(datemin_mostrecent=datemin_mostrecent, normalization="probability", verbose=verbose)
    lcu.heatmap_diagnosisVSbeddays(groupdiagnoses='mdc', datemin_mostrecent=datemin_mostrecent, fileformat='pdf', verbose=verbose)
    lcu.heatmap_diagnosisVSbeddays(groupdiagnoses=3, datemin_mostrecent=datemin_mostrecent, fileformat='pdf',verbose=verbose)
    lcu.heatmap_diagnosisVSbeddays(groupdiagnoses=None, datemin_mostrecent=datemin_mostrecent, fileformat='pdf',verbose=verbose)

    if verbose: print('\n-------------------- Full run with datemine_mostrecent = '+datemin_mostrecent+' completed --------------------')
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
def evaluate_beddays_boxplot(verbose=True):
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
    plt.plot([34.5, 34.5], [ymin, ymax], '--', color='gray', lw=lthick, zorder=5)


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
def get_kde_measures(kdedata,linenumber):
    """

    """
    xx, yy = kdedata.get_lines()[linenumber].get_data()
    cdf = scipy.integrate.cumtrapz(yy, xx, initial=0)
    nearest_50pct = np.abs(cdf - 0.5).argmin()
    x_median = xx[nearest_50pct]
    y_median = yy[nearest_50pct]
    y_max    = np.max(yy)
    x_max    = xx[np.where(yy == y_max)[0]]
    #print('        ->'+str(x_median)+'  '+str(x_max[0]))
    return x_median, x_max[0], y_max
# -----------------------------------------------------------------------------------------------------------------------
def beddays_distributions(datemin_mostrecent = '01-01-2022', normalization="probability", verbose=True):
    """
    lcu.beddays_distributions(normalization="probability",datemin_mostrecent = '01-02-2022')
    lcu.beddays_distributions(normalization="count")
    """
    plotdir = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/plots/'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    df_baseline, df_updates = lcu.load_dataframes_from_excel()

    datemin     = datetime.datetime.strptime(datemin_mostrecent, "%d-%m-%Y")
    dropval     = np.where(df_updates['Kontakt startdato Dato-tid'] < datemin)[0]

    for kdeonly in [True, False]:
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plotname = 'distributions_days_'+normalization+'.pdf'
        if kdeonly:
            plotname = plotname.replace('.pdf', '_kdeonly.pdf')
        if verbose: print(' - Initiating '+plotname)

        fig = plt.figure(figsize=(6, 6))
        fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.14, right=0.98, bottom=0.10, top=0.98)
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
        binwidth     = 1
        if kdeonly:
            distcolor = cmap(colnorm(30))
            kdedata_b = seaborn.kdeplot(data=df_baseline, x='KONTAKTDAGE', color=distcolor, alpha=1.0,
                                        common_norm=False,
                                        label='Baseline (2019-2021)' + '; ' + str(len(df_baseline)) + ' forløb')
            kde_median_baseline, kde_max_baseline, kde_ymax_baseline = lcu.get_kde_measures(kdedata_b,0)

            distcolor = cmap(colnorm(70))
            kdedata_u = seaborn.kdeplot(data=df_updates, x='Forskel på kontakt start og slut (antal dage)', color=distcolor, alpha=1.0,
                                        common_norm=False,
                                        label='Forløb siden 01-12-2021' + '; ' + str(len(df_updates)) + ' forløb')
            kde_median_updates, kde_max_updates, kde_ymax_updates = lcu.get_kde_measures(kdedata_u,1)

            distcolor = 'black'
            kdedata_l = seaborn.kdeplot(data=df_updates.drop(df_baseline.index[dropval]), x='Forskel på kontakt start og slut (antal dage)', color=distcolor, alpha=1.0,
                                        common_norm=False,
                                        label='Forløb siden ' + datemin_mostrecent + '; ' + str(len(df_updates.drop(df_baseline.index[dropval]))) + ' forløb')
            kde_median_lastmonth, kde_max_lastmonth, kde_ymax_lastmonth = lcu.get_kde_measures(kdedata_l,2)

            plt.ylim(0, plt.gca().get_ylim()[1] / binwidth)  # similar limits on the y-axis to align the plots
            plt.gca().yaxis.set_major_formatter(PercentFormatter(1 / binwidth))  # show axis such that 1/binwidth corresponds to 100%
        else:
            distcolor = cmap(colnorm(30))
            seaborn.histplot(data=df_baseline, x='KONTAKTDAGE', color=distcolor, alpha=0.3,
                             kde=True, binwidth=binwidth, discrete=True, stat=normalization, common_norm=False,
                             label='Baseline (2019-2021)'+'; '+str(len(df_baseline))+' forløb')

            distcolor = cmap(colnorm(70))
            seaborn.histplot(data=df_updates, x='Forskel på kontakt start og slut (antal dage)', color=distcolor, alpha=0.3,
                             kde=True, binwidth=binwidth, discrete=True, stat=normalization, common_norm=False,
                             label='Forløb siden 01-12-2021'+'; '+str(len(df_updates))+' forløb')

            distcolor = 'black'
            seaborn.histplot(data=df_updates.drop(df_baseline.index[dropval]), x='Forskel på kontakt start og slut (antal dage)', color=distcolor, alpha=0.3,
                             kde=True, binwidth=binwidth, discrete=True, stat=normalization, common_norm=False,
                             label='Forløb siden '+datemin_mostrecent+'; '+str(len(df_updates.drop(df_baseline.index[dropval])))+' forløb')

        # --------- ADD TEXT WITH VALUES ---------
        topleft_x = 0.4
        topleft_y = 0.75
        collist = [cmap(colnorm(30)), cmap(colnorm(70)), 'black']
        daystr  = ['KONTAKTDAGE', 'Forskel på kontakt start og slut (antal dage)', 'Forskel på kontakt start og slut (antal dage)',]
        plt.text(topleft_x, topleft_y, 'data: median,  genms. +/- stdafv.', fontsize=Fsize, rotation=0,
                 color='black', horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)

        for dd, datavals in enumerate([df_baseline, df_updates, df_updates.drop(df_baseline.index[dropval])]):
            datavals = datavals.assign(contactdays=datavals[daystr[dd]].values) # add column with name contactdays

            val_mean   = datavals.contactdays.mean()
            val_median = datavals.contactdays.median()
            val_std    = datavals.contactdays.std()

            outstr = str('%12.2f' % val_median)+' dage,  '+str('%4.1f' % val_mean)+' +/- '+str('%4.1f' % val_std)+' dage'
            plt.text(topleft_x, topleft_y-(dd+1)*0.05, outstr, fontsize=Fsize, rotation=0,
                     color=collist[dd], horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)

        if kdeonly:
            plt.text(topleft_x, topleft_y-4*0.05, 'kurve: median,  maks.', fontsize=Fsize, rotation=0,
                     color='black', horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)

            kde_max    = [kde_max_baseline, kde_max_updates, kde_max_lastmonth]
            kde_ymax   = [kde_ymax_baseline, kde_ymax_updates, kde_ymax_lastmonth]
            kde_median = [kde_median_baseline, kde_median_updates, kde_median_lastmonth]
            for dd, val_max in enumerate(kde_max):
                val_median = kde_median[dd]
                outstr = str('%12.1f' % val_median)+' dage,  ('+str('%4.1f' % val_max)+' dage, '+str('%4.1f' % (kde_ymax[dd]*100.))+'% )'
                plt.text(topleft_x, topleft_y-(dd+5)*0.05, outstr, fontsize=Fsize, rotation=0,
                         color=collist[dd], horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)

        # --------- LABELS ---------
        plt.xlabel('Kontaktdage')
        if kdeonly:
            ylabelstr = 'Procentvis fordeling af forløb'
        else:
            ylabelstr = 'Antal forløb '
            if normalization == 'probability':
                ylabelstr = ylabelstr+'(brøkdel af samlede antal forløb)'

        plt.ylabel(ylabelstr)

        # --------- LEGEND ---------
        leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=1, numpoints=1,
                         bbox_to_anchor=(0.63, 0.95))  # add the legend
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
    lcu.diagnoses_distributions(normalization="probability")
    lcu.diagnoses_distributions(normalization="count")
    """
    plotdir = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/plots/'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    df_baseline, df_updates = lcu.load_dataframes_from_excel()
    df_baseline['DIA01'] = pd.Categorical(df_baseline['DIA01'], np.sort(np.unique(np.asarray(df_baseline['DIA01']))))

    datemin     = datetime.datetime.strptime(datemin_mostrecent, "%d-%m-%Y")
    dropval     = np.where(df_updates['Kontakt startdato Dato-tid'] < datemin)[0]

    # -----------------------------------------------------------------------------------------
    if verbose: print('\n - Load MDC groups and assign them to baseline and updates data')
    mdcgroups                = loadMDCgroups.load_into_dataframe(verbose=True)
    groupnames, groupindices = loadMDCgroups.get_group_indices(mdcgroups, verbose=True)

    diacol = ['DIA01', 'Aktionsdiagnosekode']
    collist = []
    for ff, dframe in enumerate([df_baseline, df_updates]):
        mdcgroupcol  = np.zeros(len(dframe[diacol[ff]]))

        for groupno in np.arange(1, 27, 1):
            groupdia = np.asarray([dd.replace('\xa0', '') for dd in mdcgroups['diagnosekode'][groupindices['group' + str(groupno)]]])
            for dia, datadia in enumerate(dframe[diacol[ff]].values):
                if datadia in groupdia:
                    mdcgroupcol[dia] = groupno
                    #break # jump out of inner loop to advance to next diagnose

        collist.append(mdcgroupcol)

    df_baseline = df_baseline.assign(mdcgroup=collist[0]) # add column with MDC groups to dataframe
    df_updates  = df_updates.assign(mdcgroup=collist[1])  # add column with MDC groups to dataframe
    #-----------------------------------------------------------------------------------------

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'distributions_diagnoses_'+normalization+'.pdf'
    if verbose: print('\n - Initiating '+plotname)

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
    plt.xlabel('Aktionsdiagnosekoder (DIA01)')
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

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Build numpy array with diagnose stats')
    df_lastmonth    = df_updates.drop(df_updates.index[dropval])
    diaglist        = np.sort(np.unique(np.asarray(df_baseline['DIA01'].values.tolist()+ df_updates['Aktionsdiagnosekode'].values.tolist())))
    Ndiag           = len(diaglist)

    if verbose: print(' - Found '+str(Ndiag)+' diagnoses to group data by')
    Nbaseline   = len(df_baseline['DIA01'])
    Nupdates    = len(df_updates['Aktionsdiagnosekode'])
    Nlastmonth  = len(df_lastmonth['Aktionsdiagnosekode'])

    diag_baseline_count   = np.zeros(Ndiag)
    diag_update_count     = np.zeros(Ndiag)
    diag_lastmonth_count  = np.zeros(Ndiag)
    diag_baseline_frac    = np.zeros(Ndiag)
    diag_update_frac      = np.zeros(Ndiag)
    diag_lastmonth_frac   = np.zeros(Ndiag)

    for dd, diag in enumerate(diaglist):
        diag_baseline_count[dd]  = len(np.where(np.asarray(df_baseline['DIA01'].values) == diag)[0])
        diag_update_count[dd]    = len(np.where(np.asarray(df_updates['Aktionsdiagnosekode'].values) == diag)[0])
        diag_lastmonth_count[dd] = len(np.where(np.asarray(df_lastmonth['Aktionsdiagnosekode'].values) == diag)[0])
        diag_baseline_frac[dd]   = len(np.where(np.asarray(df_baseline['DIA01'].values) == diag)[0]) / Nbaseline
        diag_update_frac[dd]     = len(np.where(np.asarray(df_updates['Aktionsdiagnosekode'].values) == diag)[0]) / Nupdates
        diag_lastmonth_frac[dd]  = len(np.where(np.asarray(df_lastmonth['Aktionsdiagnosekode'].values) == diag)[0]) / Nlastmonth

    # df_diag = pd.DataFrame(list(zip(diag_baseline_count, diag_baseline_frac,
    #                                 diag_update_count, diag_update_frac,
    #                                 diag_lastmonth_count, diag_lastmonth_frac)),
    #                        columns=diaglist)

    entarray = np.arange(Ndiag)

    if normalization == 'probability':
        diag_baseline_val = diag_baseline_frac
        diag_update_val = diag_update_frac
        diag_lastmonth_val = diag_lastmonth_frac
    elif normalization == 'count':
        diag_baseline_val = diag_baseline_count
        diag_update_val = diag_update_count
        diag_lastmonth_val = diag_lastmonth_count

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'diagnoses_comparison_'+normalization+'.pdf'

    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(11, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.12, right=0.98, bottom=0.1, top=0.98)
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
    # plt.step(entarray, diag_baseline_val, where='mid', lw=lthick/2, markersize=0, alpha=0.5, color=distcolor,
    #          zorder=20., label='Baseline (2019-2021)')
    plt.fill_between(entarray, diag_baseline_val, step="mid", alpha=0.4, color=distcolor, linewidth=0,
                     zorder=20., label='Baseline (2019-2021)'+'; '+str(Nbaseline)+' forløb')

    distcolor = cmap(colnorm(70))
    # plt.step(entarray, diag_update_val, where='mid', lw=lthick/2, markersize=0, alpha=0.5, color=distcolor,
    #          zorder=20., label='Forløb siden 01-12-2021')
    plt.fill_between(entarray, diag_update_val, step="mid", alpha=0.4, color=distcolor, linewidth=0,
                     zorder=20., label='Forløb siden 01-12-2021'+'; '+str(Nupdates)+' forløb')

    distcolor = 'black'
    # plt.step(entarray, diag_lastmonth_val, where='mid', lw=lthick/2, markersize=0, alpha=0.5, color=distcolor,
    #          zorder=25., label='Forløb siden '+datemin_mostrecent)
    plt.fill_between(entarray, diag_lastmonth_val, step="mid", alpha=0.4, color=distcolor, linewidth=0,
                     zorder=25., label='Forløb siden '+datemin_mostrecent+'; '+str(Nlastmonth)+' forløb')

    spikecut = 100
    if normalization == 'probability':
        spikecut = 0.02
    diag_set  = [diag_baseline_frac, diag_update_frac, diag_lastmonth_frac]
    maxarray  = np.asarray([max(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(maxarray > spikecut)[0]

    for spike_ent in spike_ents:
        plt.text(spike_ent, maxarray[spike_ent], diaglist[spike_ent], fontsize=Fsize-2, rotation=30,
                 color='gray', horizontalalignment='left', verticalalignment='bottom')

    # --------- LABELS ---------
    plt.xticks(rotation=90, fontsize=Fsize)
    plt.xticks(entarray, diaglist, fontsize=1.4)
    #plt.xticks([])
    plt.xlabel('Aktionsdiagnosekoder (DIA01)')
    if normalization == 'probability':
        ylabelstr = 'Andelen af forløb med given diagnose'
        plt.ylim(0, plt.gca().get_ylim()[1])
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    else:
        ylabelstr = 'Antal forløb med given diagnose '

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

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'diagnoses_comparison_diff_probability.pdf'

    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(11, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.12, right=0.98, bottom=0.1, top=0.90)
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
    distcolor = 'red'
    plt.fill_between(entarray, diag_baseline_frac-diag_update_frac, step="mid", alpha=0.5, color=distcolor, linewidth=0,
                     zorder=20., label='Forskel: Baseline (2019-2021) og Forløb siden 01-12-2021')

    distcolor = 'black'
    plt.fill_between(entarray, diag_baseline_frac-diag_lastmonth_frac, step="mid", alpha=0.5, color=distcolor, linewidth=0,
                     zorder=20., label='Forskel: Baseline (2019-2021) og Forløb siden '+datemin_mostrecent)

    fraccut   = 0.02
    diag_set   = [diag_baseline_frac-diag_lastmonth_frac, diag_baseline_frac-diag_update_frac]

    spikearray = np.asarray([max(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(spikearray > fraccut)[0]
    for spike_ent in spike_ents:
        plt.text(spike_ent, spikearray[spike_ent], diaglist[spike_ent], fontsize=Fsize-2, rotation=30,
                 color='gray', horizontalalignment='left', verticalalignment='bottom')

    spikearray = np.asarray([min(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(spikearray < -1*fraccut)[0]
    for spike_ent in spike_ents:
        plt.text(spike_ent, spikearray[spike_ent], diaglist[spike_ent], fontsize=Fsize-2, rotation=30,
                 color='gray', horizontalalignment='left', verticalalignment='bottom')

    # --------- LABELS ---------
    plt.xticks(rotation=90, fontsize=Fsize)
    #plt.xticks(entarray, diaglist, fontsize=1.4)
    plt.xticks([])
    plt.xlabel('Aktionsdiagnosekoder (DIA01)')
    ylabelstr = 'Forskel i andelen af forløb med givne diagnoser '
    plt.ylabel(ylabelstr)

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=2, numpoints=1,
                     bbox_to_anchor=(0.48, 1.1))  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    if verbose: print(' - saved plot to '+plotdir)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'diagnoses_comparison_diff_probability_updates.pdf'

    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(11, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.12, right=0.98, bottom=0.1, top=0.90)
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
    distcolor = 'black'
    plt.fill_between(entarray, diag_baseline_frac-diag_update_frac, step="mid", alpha=0.5, color=distcolor, linewidth=0,
                     zorder=20., label='Forskel: Baseline (2019-2021) og Forløb siden 01-12-2021')

    fraccut   = 0.02
    diag_set   = [diag_baseline_frac-diag_update_frac]

    spikearray = np.asarray([max(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(spikearray > fraccut)[0]
    for spike_ent in spike_ents:
        plt.text(spike_ent, spikearray[spike_ent], diaglist[spike_ent], fontsize=Fsize-2, rotation=30,
                 color='gray', horizontalalignment='left', verticalalignment='bottom')

    spikearray = np.asarray([min(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(spikearray < -1*fraccut)[0]
    for spike_ent in spike_ents:
        plt.text(spike_ent, spikearray[spike_ent], diaglist[spike_ent], fontsize=Fsize-2, rotation=30,
                 color='gray', horizontalalignment='left', verticalalignment='bottom')

    # --------- LABELS ---------
    plt.xticks(rotation=90, fontsize=Fsize)
    #plt.xticks(entarray, diaglist, fontsize=1.4)
    plt.xticks([])
    plt.xlabel('Aktionsdiagnosekoder (DIA01)')
    ylabelstr = 'Forskel i andelen af forløb med givne diagnoser '
    plt.ylabel(ylabelstr)

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=2, numpoints=1,
                     bbox_to_anchor=(0.48, 1.1))  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    if verbose: print(' - saved plot to '+plotdir)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'diagnoses_comparison_diff_probability_lastmonth.pdf'

    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(11, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.12, right=0.98, bottom=0.1, top=0.90)
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
    distcolor = 'black'
    plt.fill_between(entarray, diag_baseline_frac-diag_lastmonth_frac, step="mid", alpha=0.5, color=distcolor, linewidth=0,
                     zorder=20., label='Forskel: Baseline (2019-2021) og Forløb siden '+datemin_mostrecent)

    fraccut   = 0.02
    diag_set   = [diag_baseline_frac-diag_lastmonth_frac]

    spikearray = np.asarray([max(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(spikearray > fraccut)[0]
    for spike_ent in spike_ents:
        plt.text(spike_ent, spikearray[spike_ent], diaglist[spike_ent], fontsize=Fsize-2, rotation=30,
                 color='gray', horizontalalignment='left', verticalalignment='bottom')

    spikearray = np.asarray([min(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(spikearray < -1*fraccut)[0]
    for spike_ent in spike_ents:
        plt.text(spike_ent, spikearray[spike_ent], diaglist[spike_ent], fontsize=Fsize-2, rotation=30,
                 color='gray', horizontalalignment='left', verticalalignment='bottom')

    # --------- LABELS ---------
    plt.xticks(rotation=90, fontsize=Fsize)
    #plt.xticks(entarray, diaglist, fontsize=1.4)
    plt.xticks([])
    plt.xlabel('Aktionsdiagnosekoder (DIA01)')
    ylabelstr = 'Forskel i andelen af forløb med givne diagnoser '
    plt.ylabel(ylabelstr)

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=2, numpoints=1,
                     bbox_to_anchor=(0.48, 1.1))  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    if verbose: print(' - saved plot to '+plotdir)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    #                                                   MDC grouping

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Build numpy array with diagnose MDC grouping stats')
    df_lastmonth    = df_updates.drop(df_updates.index[dropval])
    diaglist        = np.sort(np.unique(np.asarray(df_baseline['DIA01'].values.tolist()+ df_updates['Aktionsdiagnosekode'].values.tolist())))
    Ndiag           = len(diaglist)

    if verbose: print(' - Found '+str(Ndiag)+' diagnoses to group data by')
    Ngroups = 26

    diag_baseline_count   = np.zeros(Ngroups)
    diag_update_count     = np.zeros(Ngroups)
    diag_lastmonth_count  = np.zeros(Ngroups)
    diag_baseline_frac    = np.zeros(Ngroups)
    diag_update_frac      = np.zeros(Ngroups)
    diag_lastmonth_frac   = np.zeros(Ngroups)

    for dd, groupno in enumerate(np.arange(Ngroups)):
        diag_baseline_count[dd]  = len(np.where(np.asarray(df_baseline['mdcgroup'].values) == groupno+1.0)[0])
        diag_update_count[dd]    = len(np.where(np.asarray(df_updates['mdcgroup'].values) == groupno+1.0)[0])
        diag_lastmonth_count[dd] = len(np.where(np.asarray(df_lastmonth['mdcgroup'].values) == groupno+1.0)[0])
        diag_baseline_frac[dd]   = len(np.where(np.asarray(df_baseline['mdcgroup'].values) == groupno+1.0)[0]) / Nbaseline
        diag_update_frac[dd]     = len(np.where(np.asarray(df_updates['mdcgroup'].values) == groupno+1.0)[0]) / Nupdates
        diag_lastmonth_frac[dd]  = len(np.where(np.asarray(df_lastmonth['mdcgroup'].values) == groupno+1.0)[0]) / Nlastmonth

    entarray = np.arange(1,27,1)

    if normalization == 'probability':
        diag_baseline_val = diag_baseline_frac
        diag_update_val = diag_update_frac
        diag_lastmonth_val = diag_lastmonth_frac
    elif normalization == 'count':
        diag_baseline_val = diag_baseline_count
        diag_update_val = diag_update_count
        diag_lastmonth_val = diag_lastmonth_count

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'diagnoses_comparison_'+normalization+'_mdcgrouping.pdf'

    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(11, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.12, right=0.98, bottom=0.1, top=0.98)
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
    # plt.step(entarray, diag_baseline_val, where='mid', lw=lthick/2, markersize=0, alpha=0.5, color=distcolor,
    #          zorder=20., label='Baseline (2019-2021)')
    plt.fill_between(entarray, diag_baseline_val, step="pre", alpha=0.4, color=distcolor, linewidth=0,
                     zorder=20., label='Baseline (2019-2021)'+'; '+str(Nbaseline)+' forløb')

    distcolor = cmap(colnorm(70))
    # plt.step(entarray, diag_update_val, where='mid', lw=lthick/2, markersize=0, alpha=0.5, color=distcolor,
    #          zorder=20., label='Forløb siden 01-12-2021')
    plt.fill_between(entarray, diag_update_val, step="pre", alpha=0.4, color=distcolor, linewidth=0,
                     zorder=20., label='Forløb siden 01-12-2021'+'; '+str(Nupdates)+' forløb')

    distcolor = 'black'
    # plt.step(entarray, diag_lastmonth_val, where='mid', lw=lthick/2, markersize=0, alpha=0.5, color=distcolor,
    #          zorder=25., label='Forløb siden '+datemin_mostrecent)
    plt.fill_between(entarray, diag_lastmonth_val, step="pre", alpha=0.4, color=distcolor, linewidth=0,
                     zorder=25., label='Forløb siden '+datemin_mostrecent+'; '+str(Nlastmonth)+' forløb')

    spikecut = 100
    if normalization == 'probability':
        spikecut = 0.05
    diag_set  = [diag_baseline_frac, diag_update_frac, diag_lastmonth_frac]
    maxarray  = np.asarray([max(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(maxarray > spikecut)[0]

    for spike_ent in spike_ents:
        plt.text(spike_ent+0.5, maxarray[spike_ent], groupnames[spike_ent], fontsize=Fsize-2, rotation=0,
                 color='gray', horizontalalignment='center', verticalalignment='bottom')

    # --------- LABELS ---------
    #plt.xticks(rotation=90, fontsize=Fsize)
    plt.xticks(entarray-0.5, entarray.astype(str), fontsize=Fsize)
    #plt.xticks([])
    plt.xlabel('Aktionsdiagnosekoder (DIA01) fordelt på MDC grupper')
    if normalization == 'probability':
        ylabelstr = 'Andelen af forløb med given MDC gruppe'
        plt.ylim(0, plt.gca().get_ylim()[1])
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    else:
        ylabelstr = 'Antal forløb med given diagnose '

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

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'diagnoses_comparison_diff_probability_mdcgrouping.pdf'

    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(11, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.12, right=0.98, bottom=0.1, top=0.90)
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
    distcolor = 'red'
    plt.fill_between(entarray, diag_baseline_frac-diag_update_frac, step="pre", alpha=0.5, color=distcolor, linewidth=0,
                     zorder=20., label='Forskel: Baseline (2019-2021) og Forløb siden 01-12-2021')

    distcolor = 'black'
    plt.fill_between(entarray, diag_baseline_frac-diag_lastmonth_frac, step="pre", alpha=0.5, color=distcolor, linewidth=0,
                     zorder=20., label='Forskel: Baseline (2019-2021) og Forløb siden '+datemin_mostrecent)

    fraccut   = 0.02
    diag_set   = [diag_baseline_frac-diag_lastmonth_frac, diag_baseline_frac-diag_update_frac]

    spikearray = np.asarray([max(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(spikearray > fraccut)[0]
    for spike_ent in spike_ents:
        plt.text(spike_ent+0.5, spikearray[spike_ent], groupnames[spike_ent], fontsize=Fsize-2, rotation=0,
                 color='gray', horizontalalignment='center', verticalalignment='bottom')

    spikearray = np.asarray([min(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(spikearray < -1*fraccut)[0]
    for spike_ent in spike_ents:
        plt.text(spike_ent+0.5, spikearray[spike_ent], groupnames[spike_ent], fontsize=Fsize-2, rotation=0,
                 color='gray', horizontalalignment='center', verticalalignment='top')

    # --------- Indicating regions ---------
    xrange = plt.gca().get_xlim()
    plt.xlim(xrange)
    plt.plot(xrange , [0, 0], '--', color='black', zorder=30, linewidth=lthick/2.)
    plt.text(0.97, 0.96, 'Baseline > Forløb', fontsize=Fsize-2, rotation=0,
                 color='black', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes)

    plt.text(0.97, 0.04, 'Baseline < Forløb', fontsize=Fsize-2, rotation=0,
                 color='black', horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)

    # --------- LABELS ---------
    #plt.xticks(rotation=90, fontsize=Fsize)
    plt.xticks(entarray-0.5, entarray.astype(str), fontsize=Fsize)
    #plt.xticks([])
    plt.xlabel('Aktionsdiagnosekoder (DIA01) fordelt på MDC grupper')
    ylabelstr = 'Forskel i andelen af forløb med givne MDC grupper '
    plt.ylabel(ylabelstr)

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=2, numpoints=1,
                     bbox_to_anchor=(0.48, 1.1))  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    if verbose: print(' - saved plot to '+plotdir)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'diagnoses_comparison_diff_probability_updates_mdcgrouping.pdf'

    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(11, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.12, right=0.98, bottom=0.1, top=0.90)
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
    distcolor = 'black'
    plt.fill_between(entarray, diag_baseline_frac-diag_update_frac, step="pre", alpha=0.5, color=distcolor, linewidth=0,
                     zorder=20., label='Forskel: Baseline (2019-2021) og Forløb siden 01-12-2021')

    fraccut   = 0.02
    diag_set   = [diag_baseline_frac-diag_update_frac]

    spikearray = np.asarray([max(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(spikearray > fraccut)[0]
    for spike_ent in spike_ents:
        plt.text(spike_ent+0.5, spikearray[spike_ent], groupnames[spike_ent], fontsize=Fsize-2, rotation=0,
                 color='gray', horizontalalignment='center', verticalalignment='bottom')

    spikearray = np.asarray([min(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(spikearray < -1*fraccut)[0]
    for spike_ent in spike_ents:
        plt.text(spike_ent+0.5, spikearray[spike_ent], groupnames[spike_ent], fontsize=Fsize-2, rotation=0,
                 color='gray', horizontalalignment='center', verticalalignment='top')


    # --------- Indicating regions ---------
    xrange = plt.gca().get_xlim()
    plt.xlim(xrange)
    plt.plot(xrange , [0, 0], '--', color='black', zorder=30, linewidth=lthick/2.)
    plt.text(0.97, 0.96, 'Baseline > Forløb', fontsize=Fsize-2, rotation=0,
                 color='black', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes)

    plt.text(0.97, 0.04, 'Baseline < Forløb', fontsize=Fsize-2, rotation=0,
                 color='black', horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)

    # --------- LABELS ---------
    #plt.xticks(rotation=90, fontsize=Fsize)
    plt.xticks(entarray-0.5, entarray.astype(str), fontsize=Fsize)
    #plt.xticks([])
    plt.xlabel('Aktionsdiagnosekoder (DIA01) fordelt på MDC grupper')
    ylabelstr = 'Forskel i andelen af forløb med givne MDC grupper '
    plt.ylabel(ylabelstr)

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=2, numpoints=1,
                     bbox_to_anchor=(0.48, 1.1))  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    if verbose: print(' - saved plot to '+plotdir)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'diagnoses_comparison_diff_probability_lastmonth_mdcgrouping.pdf'

    if verbose: print(' - Initiating '+plotname)

    fig = plt.figure(figsize=(11, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.12, right=0.98, bottom=0.1, top=0.90)
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
    distcolor = 'black'
    plt.fill_between(entarray, diag_baseline_frac-diag_lastmonth_frac, step="pre", alpha=0.5, color=distcolor, linewidth=0,
                     zorder=20., label='Forskel: Baseline (2019-2021) og Forløb siden '+datemin_mostrecent)

    fraccut   = 0.02
    diag_set   = [diag_baseline_frac-diag_lastmonth_frac]

    spikearray = np.asarray([max(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(spikearray > fraccut)[0]
    for spike_ent in spike_ents:
        plt.text(spike_ent+0.5, spikearray[spike_ent], groupnames[spike_ent], fontsize=Fsize-2, rotation=0,
                 color='gray', horizontalalignment='center', verticalalignment='bottom')

    spikearray = np.asarray([min(idx) for idx in zip(*diag_set)])
    spike_ents = np.where(spikearray < -1*fraccut)[0]
    for spike_ent in spike_ents:
        plt.text(spike_ent+0.5, spikearray[spike_ent], groupnames[spike_ent], fontsize=Fsize-2, rotation=0,
                 color='gray', horizontalalignment='center', verticalalignment='top')

    # --------- Indicating regions ---------
    xrange = plt.gca().get_xlim()
    plt.xlim(xrange)
    plt.plot(xrange , [0, 0], '--', color='black', zorder=30, linewidth=lthick/2.)
    plt.text(0.97, 0.96, 'Baseline > Forløb', fontsize=Fsize-2, rotation=0,
                 color='black', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes)

    plt.text(0.97, 0.04, 'Baseline < Forløb', fontsize=Fsize-2, rotation=0,
                 color='black', horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)


    # --------- LABELS ---------
    #plt.xticks(rotation=90, fontsize=Fsize)
    plt.xticks(entarray-0.5, entarray.astype(str), fontsize=Fsize)
    #plt.xticks([])
    plt.xlabel('Aktionsdiagnosekoder (DIA01) fordelt på MDC grupper')
    ylabelstr = 'Forskel i andelen af forløb med givne MDC grupper '
    plt.ylabel(ylabelstr)

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=2, numpoints=1,
                     bbox_to_anchor=(0.48, 1.1))  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    if verbose: print(' - saved plot to '+plotdir)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# -----------------------------------------------------------------------------------------------------------------------
def heatmap_diagnosisVSbeddays(groupdiagnoses=4, datemin_mostrecent='01-01-2022', fileformat='pdf', verbose=True):
    """

    Parameters
    ----------
    groupdiagnoses        Possible choices are None, integer or 'mdc'. Integers mark the grouping based on diagnose level
                          whereas 'mdc' introduces a grouping according to MDC groups.
    datemin_mostrecent    minimum date of chunk of data to treat as "last month"
    fileformat            format of output figure - provide extesion fx. pdf or png
    verbose               toggle verbosity

    Returns
    -------
    Heatmap figures

    Example of use
    --------------

    lcu.heatmap_diagnosisVSbeddays(groupdiagnoses=2)


    """
    plotdir = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/plots/'
    df_baseline, df_updates = lcu.load_dataframes_from_excel(verbose=verbose)

    # -----------------------------------------------------------------------------------------
    if verbose: print('\n - Load MDC groups and assign them to baseline and updates data')
    mdcgroups                = loadMDCgroups.load_into_dataframe(verbose=True)
    groupnames, groupindices = loadMDCgroups.get_group_indices(mdcgroups, verbose=True)

    diacol = ['DIA01', 'Aktionsdiagnosekode']
    collist = []
    for ff, dframe in enumerate([df_baseline, df_updates]):
        mdcgroupcol  = np.zeros(len(dframe[diacol[ff]]))

        for groupno in np.arange(1, 27, 1):
            groupdia = np.asarray([dd.replace('\xa0', '') for dd in mdcgroups['diagnosekode'][groupindices['group' + str(groupno)]]])
            for dia, datadia in enumerate(dframe[diacol[ff]].values):
                if datadia in groupdia:
                    mdcgroupcol[dia] = groupno
                    #break # jump out of inner loop to advance to next diagnose

        collist.append(mdcgroupcol)

    df_baseline = df_baseline.assign(mdcgroup=collist[0]) # add column with MDC groups to dataframe
    df_updates  = df_updates.assign(mdcgroup=collist[1])  # add column with MDC groups to dataframe
    #-----------------------------------------------------------------------------------------

    datemin     = datetime.datetime.strptime(datemin_mostrecent, "%d-%m-%Y")
    dropval     = np.where(df_updates['Kontakt startdato Dato-tid'] < datemin)[0]

    df_lastmonth    = df_updates.drop(df_updates.index[dropval])

    diaglist        = np.sort(np.unique(np.asarray(df_baseline['DIA01'].values.tolist() + df_updates['Aktionsdiagnosekode'].values.tolist())))

    if groupdiagnoses is not None:
        if groupdiagnoses == 'mdc':
            Ncharacters = 2
            diaglist    = np.arange(1, 27, 1)
        else:
            Ncharacters = groupdiagnoses
            diaggroup   = np.sort(np.unique(np.asarray([diag[:Ncharacters] for diag in diaglist])))
            diaglist    = diaggroup
    else:
        Ncharacters = 10

    if groupdiagnoses == 'mdc':
        diag_datavalues_list = [df_baseline['mdcgroup'].values,
                                df_updates['mdcgroup'].values,
                                df_lastmonth['mdcgroup'].values]
    else:
        diag_datavalues_list = [df_baseline['DIA01'].values,
                                df_updates['Aktionsdiagnosekode'].values,
                                df_lastmonth['Aktionsdiagnosekode'].values]
    days_datavalues_list = [df_baseline['KONTAKTDAGE'].values,
                            df_updates['Forskel på kontakt start og slut (antal dage)'].values,
                            df_lastmonth['Forskel på kontakt start og slut (antal dage)'].values]
    colormap_list        = ['Blues', 'Greens', 'Greys']
    plotname_text        = ['baseline', 'updates', 'lastmonth']
    title_text           = ['Baseline (2019-2021)', 'Forløb siden 01-12-2021', 'Forøb siden '+datemin_mostrecent]

    for dent, diag_datavalues in enumerate(diag_datavalues_list):
        # diag_datavalues[diag_datavalues == np.nan] = 'tom'
        df_diaglist = np.asarray([str(diag)[:Ncharacters] for diag in diag_datavalues])
        Ncontacts   = len(df_diaglist)

        Ndiag        = len(diaglist)
        entarray     = np.arange(Ndiag)
        Nkontaktdage = 25
        kontaktdage_list = np.arange(Nkontaktdage)
        if verbose: print(' - Number of diagnoses to show: ' + str(Ndiag))

        if verbose: print(' - building array to use for heatmap')
        map2d = np.zeros([Ndiag, Nkontaktdage])
        for kk in np.arange(Nkontaktdage):
            for dd, diag in enumerate(diaglist):
                if groupdiagnoses == 'mdc':
                    df_diaglist = df_diaglist.astype(float).astype(int)
                Ninstances = len( np.where((np.asarray(days_datavalues_list[dent]) == kk) & (df_diaglist == diag))[0])

                map2d[dd, kk] = Ninstances / Ncontacts

        max_instances = np.max(map2d)
        #map2d[np.where(map2d == 0)] = np.nan

        if verbose: print('   Maximum fraction of instances (to normalize heatmap to): '+str(max_instances))
        df_map = pd.DataFrame(map2d, index=diaglist, columns=(np.arange(Nkontaktdage)))

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        plotname = 'heatmap_diagnosisVSbeddays_' + plotname_text[dent] + '.' + fileformat
        if groupdiagnoses is not None:
            if groupdiagnoses == 'mdc':
                plotname = plotname.replace('.' + fileformat, '_mdcgrouping.' + fileformat)
            elif groupdiagnoses > 0:
                plotname = plotname.replace('.' + fileformat, '_diagrouplevel'+str(groupdiagnoses)+'.' + fileformat)

        if verbose: print(' - Initiating '+plotname)

        if groupdiagnoses == 'mdc':
            ydim = 20
        else:
            ydim = 30
        fig = plt.figure(figsize=(7, ydim))
        fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.13, right=0.95, bottom=0.03, top=1.13)
        Fsize = 18
        lthick = 2
        marksize = 4
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif', size=Fsize)
        plt.rc('xtick', labelsize=Fsize)
        plt.rc('ytick', labelsize=Fsize)
        plt.clf()
        plt.ioff()
        plt.title(title_text[dent]+'; '+str(Ncontacts)+' forløb', fontsize=Fsize)

        # --------- COLOR MAP ---------
        #cmap = plt.cm.get_cmap('viridis')
        cmap = plt.cm.get_cmap(colormap_list[dent])
        cmin = 0.005
        if groupdiagnoses == 'mdc':
            cmax = 0.05
        else:
            cmax = 0.025

        colnorm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)
        cmaparr = np.linspace(cmin, cmax, num=50)
        m = plt.cm.ScalarMappable(cmap=cmap)
        m.set_array(cmaparr)

        # --------- COLOR BAR ---------
        clabel      = 'Andel af forløb'
        colshrink   = 1.0
        colaspect   = 80
        colanchor   = (0.5, 2.9)
        cextend     = 'neither'
        colbarscale = 2.1
        colaspect   = colaspect / colbarscale


        cb = plt.colorbar(m, extend=cextend, orientation='horizontal', location='top', # anchor=colanchor,
                          pad=0.02, aspect=colaspect, shrink=colshrink, use_gridspec=False, format=lambda x, _: f"{x:.2%}")
        cb.set_label(clabel, fontsize=Fsize)
        cb.ax.tick_params(labelsize=12)

        # --------- DRAW AND ADD AXIS ---------
        #seaborn.heatmap(ptable)
        plt.pcolormesh(df_map, cmap=cmap, norm=colnorm)

        plt.xticks(np.arange(0.5, len(df_map.columns), 1), kontaktdage_list, fontsize=12)
        plt.xlabel('Kontaktdage')

        if len(diaglist) > 250:
            yfontsize = 2
        elif len(diaglist) > 50:
            yfontsize = 8
        else:
            yfontsize = Fsize

        if groupdiagnoses is not None:
            if groupdiagnoses == 'mdc':
                plt.yticks(np.arange(0.5, len(diaglist), 1), [str(int(dstr)) for dstr in diaglist], fontsize=yfontsize)
                plt.ylabel('Aktionsdiagnosekoder (DIA01) fordelt på MDC grupper')
            elif groupdiagnoses > 0:
                plt.yticks(np.arange(0.5, len(diaglist), 1), [dstr+'*' for dstr in diaglist], fontsize=yfontsize)
                plt.ylabel('Aktionsdiagnosekodegrupper (DIA01)')
        else:
            plt.yticks(np.arange(0.5, len(diaglist), 1), diaglist, fontsize=yfontsize)
            plt.ylabel('Aktionsdiagnosekoder (DIA01)')



        plt.savefig(plotdir+plotname, dpi=800)
        plt.clf()
        plt.close('all')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - saved plots to '+plotdir)

#=======================================================================================================================

